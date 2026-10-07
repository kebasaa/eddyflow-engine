!***************************************************************************
! production_parallel.f90
! -----------------------
! Copyright © 2026, ETH Zurich, Jonathan Muller
!
! This file is part of EddyFlow®.
!
! EddyFlow (TM) is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version. You should have received a copy
! of the GNU General Public License along with EddyFlow (R). If not,
! see <http://www.gnu.org/licenses/>.
!
! EddyFlow® contains additional Open Source Components. The licenses
! and/or notices these Components can be found in the file LIBRARIES.txt.
!
! EddyFlow® is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
!***************************************************************************
!
! \brief       Splits the production pass - the main period loop that computes
!              the fluxes - across worker processes, with results identical to
!              a single pass.
!
!              The pre-passes split easily because their periods are
!              independent. This loop's are not quite: a few things carry from
!              one half-hour into the next, and each is rebuilt rather than
!              passed along.
!
!              - The output layout. Which columns every file has, the FLUXNET
!                slots, the headers - all of it is fixed by the first period
!                the run processes fully. The parent processes that "head"
!                first. A worker then processes the same head again, silently,
!                so every one of those latches fires in it exactly as it did
!                in the parent, without anyone having to list them.
!
!              - The storage terms need the half-hour before. A worker
!                processes the half-hour before its piece too, silently. Its
!                cuts are placed so that half-hour begins with a raw file of
!                its own, which makes the file search land on the same file a
!                single pass would - see ProdCutAllowed.
!
!              - The dynamic metadata merge keeps, for each setting, the last
!                valid value of every record applied so far. A worker applies
!                the records of the periods it skips, in order, and nothing
!                else of them.
!
!              Everything a worker writes before its piece begins lands in a
!              folder of its own and is dropped: the continuous files (full
!              output, statistics, FLUXNET...) are cut at the byte offset they
!              had when the piece began, and the parent appends what follows to
!              its own copy, piece by piece in time order. Files written per
!              period (spectra, raw-data levels) go straight into the parent's
!              folders once the piece has begun.
!
!              What the pre-passes worked out - the time-lag windows, the
!              planar-fit matrices - and the spectral frequency grid only the
!              parent surveys, the worker reads from a context file instead of
!              recomputing.
!
!              PWB time lags are classified against the last lag each gas
!              settled on, a memory with no limit. With PWB the split has three
!              phases: detection workers ('pd') read each piece only as far as
!              the PWB detection and hand back what each period's own data say
!              (m_pwb_stream); the parent feeds that through the classifier in
!              time order, continuing from where its own piece left the stream;
!              and production workers ('pr') then compute the fluxes with the
!              lags it settled, from their lead-in on.
!
!              Refused, so the loop simply runs on serially: the Billesbach
!              random error (one random stream for the whole run), raw data
!              from a shared link, embedded mode, and a run whose first fully
!              processed period left an output file still to be opened.
!
! \author      Jonathan Muller
! \sa          prepass_parallel.f90, eddyflow-rp_main.f90
!***************************************************************************
module m_production_parallel
    use m_rp_global_var
    use m_batch_pool, only: BatchDumpPath, PrepassChunkCount, ForcedPieceLength, &
        RemoveStaleWorkerRoots, WorkerRoot, RemoveWorkerRoots, NoSlashEnd, &
        AppendBytes, FileBytes
    use m_remote_source, only: RemoteSetWindow, RemotePositionOf
    use m_pwb_timelag, only: InitPwbTimelagCache, AppendPwbCacheRows
    use m_pwb_stream, only: PwbEvidenceType, PwbVerdictType, PwbReplayEvidence, &
        WritePwbEvidence, ReadPwbEvidence, WritePwbVerdict, ReadPwbVerdict
    implicit none
    private
    public :: ProdForcedPieceLength, PlanProdCuts
    public :: AdoptProdWorkerOutput, CaptureProdContext
    public :: WriteProdContext, ReadProdContext
    public :: ProdPieceBegins, FinishProdWorker, MergeProdPieces
    public :: KeepPwbEvidence, FinishEvidenceWorker
    public :: KeepPwbVerdict, ReplayProdEvidence
    public :: ReadProdVerdicts, FindProdVerdict, RemoveWorkerRoots
    public :: RemoveStaleWorkerRoots
    public :: SetProdRemoteWindow

    !> Both files are read back by the same binary in the same run, so the
    !> format never has to survive a version change; the magic guards against a
    !> file a crashed earlier run left behind.
    character(20), parameter :: CtxMagic  = 'EDDYFLOW_PRODCTX_01 '
    character(20), parameter :: DumpMagic = 'EDDYFLOW_PRODOUT_01 '
    character(20), parameter :: EvidMagic = 'EDDYFLOW_PRODEVID_01'
    character(20), parameter :: VerdMagic = 'EDDYFLOW_PRODVERD_01'

    !> The units of the continuous output files, uqc to uflxnt. Each is opened
    !> once, while the head is processed, and written a row per period.
    integer, parameter :: FirstOutUnit = uqc
    integer, parameter :: LastOutUnit = uflxnt

    !> Pieces per worker. Each piece costs its worker the head and one
    !> half-hour of lead-in on top of its own periods - twice with PWB - so
    !> fewer than the pre-passes' four; two still lets a fast core take more
    !> than a slow one. But only once a piece would hold LongPiece half-hours
    !> or more: below that the fixed cost dominates, and a second round of
    !> pieces pays it again. A Yatir day cut into 11 pieces for 8 workers ran
    !> its main pass 1.5 times faster than one process did.
    integer, parameter :: PiecesPerWorker = 2
    integer, parameter :: LongPiece = 16
    integer, parameter :: MaxPieces = 99
    real(kind = dbl), parameter :: EmptyPeriodWeight = 0.05d0

    !> What the main pass starts from, kept when it starts so the workers get
    !> exactly that - not what the head has made of it since.
    type(MethType) :: CtxMeth
    type(TOSetupType) :: CtxTOSetup
    type(PFSetupType) :: CtxPFSetup
    type(TimeLagType) :: CtxToPasGas(E2NumVar)
    type(TimeLagType) :: CtxToH2O(toMaxH2OClass, E2NumVar)
    logical :: CtxTimeLagOptSelected = .false.
    real(kind = dbl) :: CtxPFMat(3, 3, MaxNumWSect)
    real(kind = dbl) :: CtxPFb(3, MaxNumWSect)
    logical, allocatable :: CtxGoPlanarFit(:)
    real(kind = dbl) :: CtxBinGridAcFreq = -1d0
    real(kind = dbl), allocatable :: CtxBf(:)
    !> The PWB settings and time-lag table as the main pass starts: a table
    !> read from file, or the one the cache-generating pre-pass just settled,
    !> which a worker - skipping that pre-pass - would otherwise not have.
    type(PWBSetupType) :: CtxPWBSetup
    logical :: CtxPwbCacheLoaded = .false.
    integer :: CtxPwbN = 0
    type(PWBTimelagCacheEntryType), allocatable :: CtxPwbRows(:)

    !> PWB in a split run. A detection worker writes each period's evidence
    !> to its dump as it goes; the parent keeps every verdict from the head on,
    !> its own and those it replays; a production worker holds the verdicts of
    !> its lead-in and piece.
    !> newunit numbers are negative, so whether the dump is open is kept apart.
    integer :: EvidenceUnit = 0
    logical :: EvidenceOpen = .false.
    !> The parent keeps only its last verdict - the lead-in of the piece after
    !> it - and writes every other straight into the file of the piece it
    !> belongs to, so its memory does not grow with the run. A production
    !> worker holds the verdicts of its own piece.
    type(PwbVerdictType) :: LastVerdict
    logical :: HaveLastVerdict = .false.
    type(PwbVerdictType), allocatable :: Verdicts(:)
    integer :: nVerdicts = 0

    !> Worker: its parent's output folder, and where each continuous file stood
    !> when its piece began.
    character(PathLen) :: ParentMainOut = ''
    character(PathLen) :: WorkerMainOut = ''
    logical :: PieceOpen(FirstOutUnit:LastOutUnit) = .false.
    integer(8) :: PieceOffset(FirstOutUnit:LastOutUnit) = 0
    character(PathLen) :: PiecePath(FirstOutUnit:LastOutUnit) = ''

contains

    !***************************************************************************
    !> \brief Test hook: production pieces of this many periods. 0 if unset.
    !***************************************************************************
    integer function ProdForcedPieceLength()
        ProdForcedPieceLength = ForcedPieceLength('EDDYFLOW_PROD_PIECE_PERIODS')
    end function ProdForcedPieceLength

    !***************************************************************************
    !> \brief Whether a piece may begin at period p.
    !>
    !> A worker processes period p - 1 before its piece, silently, so that the
    !> storage terms of p see the half-hour before them. That half-hour has to
    !> come out exactly as in a single pass, and the one thing it inherits that
    !> a worker cannot rebuild is the raw-file cursor: a single pass searches
    !> for p - 1's first file from wherever the period before left off, a
    !> worker from wherever its head did, further back.
    !>
    !> Both searches find the same file when p - 1 begins cleanly with a file
    !> of its own: one that is relevant to p - 1, not to p - 2, and with no
    !> file just before it relevant to p - 1 as well. Then nothing a single
    !> pass imported earlier can have moved its cursor past that file, and
    !> nothing between the worker's cursor and that file is relevant. Files
    !> that overlap, duplicates, a file half-consumed by the period before -
    !> any of those makes the cut unsafe, and the cut moves on to a later
    !> period.
    !***************************************************************************
    logical function ProdCutAllowed(p, Series, nSeries, Files, nFiles, cursor)
        integer, intent(in) :: p
        integer, intent(in) :: nSeries
        integer, intent(in) :: nFiles
        type(DateType), intent(in) :: Series(nSeries)
        type(FileListType), intent(in) :: Files(nFiles)
        !> First file that may still be relevant to p - 1; advanced in place,
        !> so a forward sweep over p costs one pass over the files.
        integer, intent(inout) :: cursor
        integer :: q
        integer :: i
        integer :: j
        logical, external :: FileIsRelevantToCurrentPeriod

        ProdCutAllowed = .false.
        q = p - 1
        if (q < 2 .or. q + 1 > nSeries) return
        do while (cursor <= nFiles)
            if (Files(cursor)%timestamp + DatafileDateStep > Series(q - 1)) exit
            cursor = cursor + 1
        end do
        do i = cursor, nFiles
            if (Files(i)%timestamp >= Series(q + 1)) return
            if (.not. FileIsRelevantToCurrentPeriod(Files(i)%name, &
                Series(q), Series(q + 1))) cycle
            if (FileIsRelevantToCurrentPeriod(Files(i)%name, &
                Series(q - 1), Series(q))) return
            do j = max(1, i - 4), i - 1
                if (FileIsRelevantToCurrentPeriod(Files(j)%name, &
                    Series(q), Series(q + 1))) return
            end do
            ProdCutAllowed = .true.
            return
        end do
    end function ProdCutAllowed

    !***************************************************************************
    !> \brief Raw data from a shared link: which downloaded files this process
    !>        may delete, and which it may fetch ahead, in a split pass.
    !>
    !> Every process of a run shares one staging folder, and a single pass
    !> deletes each file once it is two behind the one being read. Split, that
    !> rule deletes files other processes still read: every worker re-reads the
    !> head, and the half-hour before each cut is read by the pieces on both
    !> sides of it. So each process deletes only the files its own piece alone
    !> reads, and fetches ahead only within its own piece, so it never fetches
    !> a file a neighbour has already finished with and deleted - which would
    !> download it twice, and leave it behind. The files no window covers - the
    !> head's, the boundaries' - stay until the parent's RemoteCleanup.
    !>
    !> role: 'head' - the parent before it splits, deleting nothing;
    !>       'p1'   - the parent in its own piece [pieceStart, pieceEnd);
    !>       'pd'   - a detection worker, deleting nothing: its production
    !>                twin reads the same piece afterwards;
    !>       'pr'   - a production worker.
    !***************************************************************************
    subroutine SetProdRemoteWindow(role, pieceStart, pieceEnd, lastPiece, headEnd, &
            Series, nSeries, Files, nFiles)
        character(*), intent(in) :: role
        integer, intent(in) :: pieceStart
        integer, intent(in) :: pieceEnd
        logical, intent(in) :: lastPiece
        integer, intent(in) :: headEnd
        integer, intent(in) :: nSeries
        integer, intent(in) :: nFiles
        type(DateType), intent(in) :: Series(nSeries)
        type(FileListType), intent(in) :: Files(nFiles)
        integer :: lo
        integer :: hi
        integer :: flo
        integer :: fhi

        if (nFiles < 1) return
        if (RemotePositionOf(Files(1)%path) == 0) return

        select case (role)
            case ('head')
                call RemoteSetWindow(0, 0, 1, huge(1))
                return
            case ('p1')
                lo = PosOf(LastFileOf(headEnd))
                flo = 1
            case default
                lo = PosOf(FirstFileOf(pieceStart - 1))
                flo = lo
        end select
        if (lastPiece) then
            hi = huge(1)
            fhi = huge(1)
        else
            hi = PosOf(FirstFileOf(pieceEnd - 1)) - 1
            fhi = PosOf(LastFileOf(pieceEnd - 1))
        end if
        !> An edge that could not be found deletes nothing rather than guess:
        !> lo at 0 would reach back into the head.
        if (lo <= 0 .or. hi <= lo) then
            lo = 0
            hi = 0
        end if
        if (role == 'pd') then
            call RemoteSetWindow(0, 0, flo, fhi)
        else
            call RemoteSetWindow(lo, hi, flo, fhi)
        end if

    contains

        !> The first and last file a period reads; 0 if none.
        integer function FirstFileOf(q)
            integer, intent(in) :: q
            integer :: i
            logical, external :: FileIsRelevantToCurrentPeriod

            FirstFileOf = 0
            if (q < 1 .or. q + 1 > nSeries) return
            do i = 1, nFiles
                if (FileIsRelevantToCurrentPeriod(Files(i)%name, Series(q), Series(q + 1))) then
                    FirstFileOf = i
                    return
                end if
            end do
        end function FirstFileOf

        integer function LastFileOf(q)
            integer, intent(in) :: q
            integer :: i
            logical, external :: FileIsRelevantToCurrentPeriod

            LastFileOf = 0
            if (q < 1 .or. q + 1 > nSeries) return
            do i = nFiles, 1, -1
                if (FileIsRelevantToCurrentPeriod(Files(i)%name, Series(q), Series(q + 1))) then
                    LastFileOf = i
                    return
                end if
            end do
        end function LastFileOf

        !> A list index's processing position. A period without a file has no
        !> window edge to give: nothing is deleted then, rather than guessing.
        integer function PosOf(i)
            integer, intent(in) :: i

            PosOf = 0
            if (i >= 1 .and. i <= nFiles) PosOf = RemotePositionOf(Files(i)%path)
        end function PosOf
    end subroutine SetProdRemoteWindow

    !***************************************************************************
    !> \brief Where to cut [iStart, iEnd) into pieces.
    !>
    !> Weighted by work as the pre-passes are (a period with a raw file weighs
    !> 1, one without almost nothing), and each cut moved forward to the next
    !> period ProdCutAllowed accepts. With forced, a cut every `forced`
    !> periods instead, moved forward the same way.
    !***************************************************************************
    subroutine PlanProdCuts(iStart, iEnd, nEff, forced, Series, nSeries, &
            Files, nFiles, cuts, nCuts)
        integer, intent(in) :: iStart
        integer, intent(in) :: iEnd
        integer, intent(in) :: nEff
        integer, intent(in) :: forced
        integer, intent(in) :: nSeries
        integer, intent(in) :: nFiles
        type(DateType), intent(in) :: Series(nSeries)
        type(FileListType), intent(in) :: Files(nFiles)
        integer, intent(out) :: cuts(MaxPieces)
        integer, intent(out) :: nCuts
        integer :: p
        integer :: j
        integer :: want
        integer :: cursor
        real(kind = dbl) :: total
        real(kind = dbl) :: cum
        real(kind = dbl) :: share
        real(kind = dbl), allocatable :: weight(:)
        logical, allocatable :: ok(:)

        nCuts = 0
        if (iEnd - iStart < 2) return

        allocate(weight(iStart:iEnd - 1), ok(iStart:iEnd - 1))
        j = 1
        do p = iStart, iEnd - 1
            do while (j <= nFiles)
                if (Files(j)%timestamp + DatafileDateStep > Series(p)) exit
                j = j + 1
            end do
            weight(p) = EmptyPeriodWeight
            if (j <= nFiles) then
                if (Files(j)%timestamp < Series(p + 1)) weight(p) = 1d0
            end if
        end do
        ok = .false.
        cursor = 1
        do p = iStart + 1, iEnd - 1
            ok(p) = ProdCutAllowed(p, Series, nSeries, Files, nFiles, cursor)
        end do

        if (forced > 0) then
            p = iStart + forced
            do while (p <= iEnd - 1 .and. nCuts < MaxPieces - 1)
                if (ok(p)) then
                    nCuts = nCuts + 1
                    cuts(nCuts) = p
                    p = p + forced
                else
                    p = p + 1
                end if
            end do
        else
            if ((iEnd - iStart) / (PiecesPerWorker * nEff) >= LongPiece) then
                want = min(PiecesPerWorker * nEff, MaxPieces)
            else
                want = min(nEff, MaxPieces, max(1, (iEnd - iStart) / 4))
            end if
            want = max(want, min(nEff, iEnd - iStart))
            total = sum(weight)
            cum = 0d0
            do p = iStart + 1, iEnd - 1
                cum = cum + weight(p - 1)
                if (nCuts + 1 >= want) exit
                share = total * dble(nCuts + 1) / dble(want)
                if (cum >= share .and. ok(p)) then
                    nCuts = nCuts + 1
                    cuts(nCuts) = p
                end if
            end do
        end if
        deallocate(weight, ok)
    end subroutine PlanProdCuts

    !***************************************************************************
    !> \brief Worker: write into a folder of its own, under its parent's names.
    !>
    !> Called before the output folder is created. Everything the worker
    !> writes there - the head, the lead-in, and its piece - is its own; the
    !> parent takes from it only what its piece wrote.
    !***************************************************************************
    subroutine AdoptProdWorkerOutput()
        integer :: u
        integer :: io_status
        character(20) :: magic

        open(newunit = u, file = trim(CtxPath()), form = 'unformatted', &
            access = 'stream', status = 'old', action = 'read', iostat = io_status)
        if (io_status /= 0) error stop 'Production worker could not open its context file.'
        read(u, iostat = io_status) magic
        if (io_status /= 0 .or. magic /= CtxMagic) &
            error stop 'Production worker context file is not one of this run.'
        read(u, iostat = io_status) Timestamp_FilePadding
        close(u)
        if (io_status /= 0) error stop 'Production worker context file is truncated.'

        ParentMainOut = Dir%main_out
        WorkerMainOut = WorkerRoot(BatchParentPid, BatchKind, BatchIndex)
        Dir%main_out = WorkerMainOut
        BatchOwnOutDir = WorkerMainOut
        !> Empty before it is used. A folder of the same name can only be a
        !> killed run's whose process id Windows has since handed to this
        !> parent - and a file opened over an old one keeps the old bytes past
        !> what is written, which the merge would then append.
        io_status = system(trim(comm_rmdir) // ' "' &
            // trim(NoSlashEnd(WorkerMainOut)) // '"' // comm_err_redirect)
    end subroutine AdoptProdWorkerOutput

    !***************************************************************************
    !> \brief Parent: keep what the main pass starts from.
    !***************************************************************************
    subroutine CaptureProdContext(GoPlanarFit, bf)
        logical, allocatable, intent(in) :: GoPlanarFit(:)
        real(kind = dbl), allocatable, intent(in) :: bf(:)

        CtxMeth = Meth
        CtxTOSetup = TOSetup
        CtxPFSetup = PFSetup
        CtxToPasGas = toPasGas
        CtxToH2O = toH2O
        CtxTimeLagOptSelected = TimeLagOptSelected
        CtxPFMat = PFMat
        CtxPFb = PFb
        if (allocated(CtxGoPlanarFit)) deallocate(CtxGoPlanarFit)
        if (allocated(GoPlanarFit)) then
            allocate(CtxGoPlanarFit(size(GoPlanarFit)))
            CtxGoPlanarFit = GoPlanarFit
        end if
        CtxBinGridAcFreq = BinGridAcFreq
        if (allocated(CtxBf)) deallocate(CtxBf)
        if (allocated(bf)) then
            allocate(CtxBf(size(bf)))
            CtxBf = bf
        end if
        CtxPWBSetup = PWBSetup
        CtxPwbCacheLoaded = PwbCacheLoaded
        CtxPwbN = PwbTimelagCacheN
        if (allocated(CtxPwbRows)) deallocate(CtxPwbRows)
        allocate(CtxPwbRows(max(1, CtxPwbN)))
        if (CtxPwbN > 0) CtxPwbRows(1:CtxPwbN) = PwbTimelagCache(1:CtxPwbN)
    end subroutine CaptureProdContext

    !***************************************************************************
    !> \brief Parent: write the context the workers start from.
    !***************************************************************************
    subroutine WriteProdContext(headEnd)
        integer, intent(in) :: headEnd
        integer :: u
        integer :: n
        integer :: io_status

        open(newunit = u, file = trim(CtxPath()), form = 'unformatted', &
            access = 'stream', status = 'replace', iostat = io_status)
        if (io_status /= 0) error stop 'Could not write the production context file.'
        write(u) CtxMagic
        write(u) Timestamp_FilePadding
        write(u) headEnd
        write(u) CtxMeth, CtxTOSetup, CtxPFSetup, CtxToPasGas, CtxToH2O
        write(u) CtxTimeLagOptSelected, CtxPFMat, CtxPFb
        n = -1
        if (allocated(CtxGoPlanarFit)) n = size(CtxGoPlanarFit)
        write(u) n
        if (n > 0) write(u) CtxGoPlanarFit
        write(u) CtxBinGridAcFreq
        n = -1
        if (allocated(CtxBf)) n = size(CtxBf)
        write(u) n
        if (n > 0) write(u) CtxBf
        write(u) CtxPWBSetup, CtxPwbCacheLoaded, CtxPwbN
        if (CtxPwbN > 0) write(u) CtxPwbRows(1:CtxPwbN)
        write(u) CtxMagic
        close(u)
    end subroutine WriteProdContext

    !***************************************************************************
    !> \brief Worker: take up the context, at the start of the main pass.
    !***************************************************************************
    subroutine ReadProdContext(headEnd, GoPlanarFit, bf)
        integer, intent(out) :: headEnd
        logical, allocatable, intent(inout) :: GoPlanarFit(:)
        real(kind = dbl), allocatable, intent(inout) :: bf(:)
        integer :: u
        integer :: n
        integer :: io_status
        character(20) :: magic
        character(len(Timestamp_FilePadding)) :: padding

        open(newunit = u, file = trim(CtxPath()), form = 'unformatted', &
            access = 'stream', status = 'old', action = 'read', iostat = io_status)
        if (io_status /= 0) error stop 'Production worker could not open its context file.'
        read(u, iostat = io_status) magic
        if (io_status /= 0 .or. magic /= CtxMagic) &
            error stop 'Production worker context file is not one of this run.'
        read(u) padding
        read(u) headEnd
        read(u) Meth, TOSetup, PFSetup, toPasGas, toH2O
        read(u) TimeLagOptSelected, PFMat, PFb
        read(u) n
        if (allocated(GoPlanarFit)) deallocate(GoPlanarFit)
        if (n >= 0) then
            allocate(GoPlanarFit(n))
            if (n > 0) read(u) GoPlanarFit
        end if
        read(u) BinGridAcFreq
        read(u) n
        if (allocated(bf)) deallocate(bf)
        if (n >= 0) then
            allocate(bf(n))
            if (n > 0) read(u) bf
        end if
        read(u) PWBSetup, CtxPwbCacheLoaded, CtxPwbN
        call InitPwbTimelagCache()
        if (CtxPwbN > 0) then
            if (allocated(CtxPwbRows)) deallocate(CtxPwbRows)
            allocate(CtxPwbRows(CtxPwbN))
            read(u) CtxPwbRows
            call AppendPwbCacheRows(CtxPwbRows, CtxPwbN)
            deallocate(CtxPwbRows)
        end if
        PwbCacheLoaded = CtxPwbCacheLoaded
        read(u, iostat = io_status) magic
        close(u)
        if (io_status /= 0 .or. magic /= CtxMagic) &
            error stop 'Production worker context file is truncated.'
    end subroutine ReadProdContext

    !***************************************************************************
    !> \brief Worker: its piece begins.
    !>
    !> Notes where each continuous file stands - what follows is the piece's -
    !> starts the accumulators afresh, and points the per-period folders at the
    !> parent's, so a period's spectra and raw-data levels land where a single
    !> pass would have put them.
    !***************************************************************************
    subroutine ProdPieceBegins(nOk)
        integer, intent(inout) :: nOk
        integer :: u
        integer :: k
        logical :: op
        logical :: any_full
        character(PathLen) :: nm

        nOk = 0
        StorCacheN = 0
        if (allocated(StorCache)) deallocate(StorCache)
        !> Which gas lent each PWB lag, tallied per period for the aggregate
        !> summary's choice of donor.
        PwbSummaryDonorCount = 0

        do u = FirstOutUnit, LastOutUnit
            inquire(unit = u, opened = op, name = nm)
            PieceOpen(u) = op
            PieceOffset(u) = 0
            PiecePath(u) = ''
            if (.not. op) cycle
            flush(u)
            PiecePath(u) = nm
            PieceOffset(u) = FileBytes(nm)
        end do

        do k = 1, 7
            if (RPsetup%out_raw(k)) call ToParentTree(RawSubDir(k))
        end do
        if (RPsetup%out_bin_sp) call ToParentTree(BinCospectraDir)
        if (RPsetup%out_bin_og) call ToParentTree(BinOgivesDir)
        any_full = .false.
        do k = 1, GHGNumVar
            if (RPsetup%out_full_sp(k) .or. RPsetup%out_full_cosp(k)) any_full = .true.
        end do
        if (any_full) call ToParentTree(CospectraDir)
    end subroutine ProdPieceBegins

    !***************************************************************************
    !> \brief Worker: its piece is done. Hand back where its rows begin in each
    !>        file, and what it added to the accumulators, then stop.
    !***************************************************************************
    subroutine FinishProdWorker(nOk, pwbRows, pwbN)
        integer, intent(in) :: nOk
        !> The rows of the PWB time-lag summary its piece added, in order.
        type(TimeLagOptType), intent(in) :: pwbRows(:)
        integer, intent(in) :: pwbN
        integer :: u
        integer :: n
        integer :: ud
        integer :: io_status
        logical :: op

        !> A file opened during the piece and not before has its header inside
        !> the piece's rows, or not, depending on where the cut fell. The
        !> parent refuses to split whenever that could happen; this makes sure.
        do u = FirstOutUnit, LastOutUnit
            inquire(unit = u, opened = op)
            if (op .and. .not. PieceOpen(u)) then
                call LogSay(' A production worker opened an output file during its')
                call LogSay(' piece that was not open when the piece began.')
                error stop 'Production worker output cannot be merged.'
            end if
            if (op) close(u)
        end do

        open(newunit = ud, file = trim(BatchOutPath), form = 'unformatted', &
            access = 'stream', status = 'replace', iostat = io_status)
        if (io_status /= 0) error stop 'Production worker could not write its dump.'
        write(ud) DumpMagic
        write(ud) nOk, NumUserVar, nbVars
        n = count(PieceOpen)
        write(ud) n
        do u = FirstOutUnit, LastOutUnit
            if (.not. PieceOpen(u)) cycle
            write(ud) u, PieceOffset(u), PiecePath(u)
        end do
        write(ud) StorCacheN
        if (StorCacheN > 0) write(ud) StorCache(1:StorCacheN)
        write(ud) pwbN
        if (pwbN > 0) write(ud) pwbRows(1:pwbN)
        write(ud) PwbSummaryDonorCount
        write(ud) DumpMagic
        close(ud)
    end subroutine FinishProdWorker

    !***************************************************************************
    !> \brief Parent: append every worker's rows to its own files, in piece
    !>        order, and add up what they counted.
    !>
    !> After this the files hold what a single pass would have written, and
    !> the end-of-run code runs on them unchanged. The two counts the end of
    !> the run reads from the last period - how many custom variables, how many
    !> biomet variables - are taken from the last piece, which is where a
    !> single pass would have left them.
    !***************************************************************************
    subroutine MergeProdPieces(nOk, pwbRows, pwbSize, pwbN)
        integer, intent(inout) :: nOk
        integer, intent(in) :: pwbSize
        type(TimeLagOptType), intent(inout) :: pwbRows(pwbSize)
        integer, intent(inout) :: pwbN
        integer :: wPwbN
        integer :: wDonors(E2NumVar, E2NumVar)
        integer :: k
        integer :: u
        integer :: i
        integer :: n
        integer :: ud
        integer :: io_status
        integer :: nPieces
        integer :: wOk
        integer :: wUser
        integer :: wBiomet
        integer :: wStorN
        integer(8) :: offset
        logical :: op
        logical :: seen(FirstOutUnit:LastOutUnit)
        logical :: parentOpen(FirstOutUnit:LastOutUnit)
        character(PathLen) :: parentPath(FirstOutUnit:LastOutUnit)
        character(PathLen) :: path
        character(20) :: magic
        type(StorCacheEntryType), allocatable :: rows(:)
        type(StorCacheEntryType), allocatable :: tmp(:)

        do u = FirstOutUnit, LastOutUnit
            inquire(unit = u, opened = op, name = path)
            parentOpen(u) = op
            parentPath(u) = ''
            if (.not. op) cycle
            parentPath(u) = path
            close(u)
        end do

        nPieces = PrepassChunkCount()
        do k = 2, nPieces
            open(newunit = ud, file = trim(BatchDumpPath('pr', k)), &
                form = 'unformatted', access = 'stream', status = 'old', &
                action = 'read', iostat = io_status)
            if (io_status /= 0) error stop 'Could not read a production worker dump.'
            read(ud, iostat = io_status) magic
            if (io_status /= 0 .or. magic /= DumpMagic) &
                error stop 'A production worker dump is not one of this run.'
            read(ud) wOk, wUser, wBiomet
            read(ud) n
            seen = .false.
            do i = 1, n
                read(ud) u, offset, path
                if (u < FirstOutUnit .or. u > LastOutUnit) &
                    error stop 'A production worker dump names an unknown file.'
                if (.not. parentOpen(u)) then
                    call LogSay(' A production worker wrote ' // trim(path))
                    call LogSay(' which this process never opened.')
                    error stop 'Production worker output cannot be merged.'
                end if
                seen(u) = .true.
                call AppendBytes(path, offset, parentPath(u))
            end do
            do u = FirstOutUnit, LastOutUnit
                if (parentOpen(u) .and. .not. seen(u)) then
                    call LogSay(' A production worker did not write ' // trim(parentPath(u)))
                    error stop 'Production worker output cannot be merged.'
                end if
            end do
            read(ud) wStorN
            if (wStorN > 0) then
                allocate(rows(wStorN))
                read(ud) rows
                allocate(tmp(StorCacheN + wStorN))
                if (StorCacheN > 0) tmp(1:StorCacheN) = StorCache(1:StorCacheN)
                tmp(StorCacheN + 1:StorCacheN + wStorN) = rows
                call move_alloc(tmp, StorCache)
                StorCacheN = StorCacheN + wStorN
                deallocate(rows)
            end if
            read(ud) wPwbN
            if (wPwbN > 0) then
                if (pwbN + wPwbN > pwbSize) &
                    error stop 'PWB time-lag summary rows from the workers do not fit.'
                read(ud) pwbRows(pwbN + 1:pwbN + wPwbN)
                pwbN = pwbN + wPwbN
            end if
            read(ud) wDonors
            PwbSummaryDonorCount = PwbSummaryDonorCount + wDonors
            read(ud, iostat = io_status) magic
            close(ud)
            if (io_status /= 0 .or. magic /= DumpMagic) &
                error stop 'A production worker dump is truncated.'
            nOk = nOk + wOk
            if (k == nPieces) then
                NumUserVar = wUser
                nbVars = wBiomet
            end if
        end do
        call RemoveWorkerRoots('pr')
    end subroutine MergeProdPieces

    !***************************************************************************
    !> \brief Detection worker: keep this period's PWB evidence.
    !***************************************************************************
    subroutine KeepPwbEvidence(ev)
        type(PwbEvidenceType), intent(in) :: ev

        call OpenEvidenceDump()
        write(EvidenceUnit) .true.
        call WritePwbEvidence(EvidenceUnit, ev)
    end subroutine KeepPwbEvidence

    subroutine OpenEvidenceDump()
        integer :: io_status

        if (EvidenceOpen) return
        open(newunit = EvidenceUnit, file = trim(BatchOutPath), form = 'unformatted', &
            access = 'stream', status = 'replace', iostat = io_status)
        if (io_status /= 0) error stop 'Detection worker could not write its dump.'
        EvidenceOpen = .true.
        write(EvidenceUnit) EvidMagic
    end subroutine OpenEvidenceDump

    !***************************************************************************
    !> \brief Detection worker: its piece is done; close the evidence off.
    !***************************************************************************
    subroutine FinishEvidenceWorker()
        call OpenEvidenceDump()
        write(EvidenceUnit) .false.
        write(EvidenceUnit) EvidMagic
        close(EvidenceUnit)
        EvidenceOpen = .false.
    end subroutine FinishEvidenceWorker

    !***************************************************************************
    !> \brief Parent: keep the verdict the detection call just made, for the
    !>        workers whose lead-in falls on this period.
    !***************************************************************************
    subroutine KeepPwbVerdict(v)
        type(PwbVerdictType), intent(in) :: v

        LastVerdict = v
        HaveLastVerdict = .true.
    end subroutine KeepPwbVerdict

    !***************************************************************************
    !> \brief Parent: classify the detection workers' evidence, piece by piece
    !>        in time order, continuing the stream this process left at the end
    !>        of its own piece, and give each production worker the verdicts
    !>        it needs - its lead-in's and its piece's - as they come.
    !***************************************************************************
    subroutine ReplayProdEvidence(cutStarts, nPieces)
        integer, intent(in) :: nPieces
        integer, intent(in) :: cutStarts(nPieces)
        integer :: k
        integer :: ud
        integer :: uv
        integer :: io_status
        logical :: more
        logical :: havePrev
        character(20) :: magic
        type(PwbEvidenceType) :: ev
        type(PwbVerdictType) :: v
        type(PwbVerdictType) :: prev

        do k = 2, nPieces
            open(newunit = ud, file = trim(BatchDumpPath('pd', k)), &
                form = 'unformatted', access = 'stream', status = 'old', &
                action = 'read', iostat = io_status)
            if (io_status /= 0) error stop 'Could not read a detection worker dump.'
            read(ud, iostat = io_status) magic
            if (io_status /= 0 .or. magic /= EvidMagic) &
                error stop 'A detection worker dump is not one of this run.'
            open(newunit = uv, file = trim(VerdictPath(k)), form = 'unformatted', &
                access = 'stream', status = 'replace', iostat = io_status)
            if (io_status /= 0) error stop 'Could not write a PWB verdict file.'
            write(uv) VerdMagic
            havePrev = .false.
            !> The lead-in's verdict, when the half-hour before the piece was
            !> classified at all
            if (HaveLastVerdict) then
                if (LastVerdict%pcount >= cutStarts(k) - 1) then
                    write(uv) .true.
                    call WritePwbVerdict(uv, LastVerdict)
                    prev = LastVerdict
                    havePrev = .true.
                end if
            end if
            do
                read(ud, iostat = io_status) more
                if (io_status /= 0) error stop 'A detection worker dump is truncated.'
                if (.not. more) exit
                call ReadPwbEvidence(ud, ev)
                call PwbReplayEvidence(ev, v)
                write(uv) .true.
                if (havePrev) then
                    call WritePwbVerdict(uv, v, prev)
                else
                    call WritePwbVerdict(uv, v)
                end if
                prev = v
                havePrev = .true.
                call KeepPwbVerdict(v)
            end do
            write(uv) .false.
            write(uv) VerdMagic
            close(uv)
            read(ud, iostat = io_status) magic
            close(ud)
            if (io_status /= 0 .or. magic /= EvidMagic) &
                error stop 'A detection worker dump is truncated.'
        end do
        call RemoveWorkerRoots('pd')
    end subroutine ReplayProdEvidence

    !***************************************************************************
    !> \brief Production worker: take up the verdicts its parent settled.
    !***************************************************************************
    subroutine ReadProdVerdicts()
        integer :: ud
        integer :: io_status
        logical :: more
        character(20) :: magic
        type(PwbVerdictType), allocatable :: grown(:)

        open(newunit = ud, file = trim(VerdictPath(BatchIndex)), &
            form = 'unformatted', access = 'stream', status = 'old', &
            action = 'read', iostat = io_status)
        if (io_status /= 0) error stop 'Production worker could not read its PWB verdicts.'
        read(ud, iostat = io_status) magic
        if (io_status /= 0 .or. magic /= VerdMagic) &
            error stop 'PWB verdict file is not one of this run.'
        if (allocated(Verdicts)) deallocate(Verdicts)
        allocate(Verdicts(64))
        nVerdicts = 0
        do
            read(ud, iostat = io_status) more
            if (io_status /= 0) error stop 'PWB verdict file is truncated.'
            if (.not. more) exit
            if (nVerdicts == size(Verdicts)) then
                allocate(grown(2 * size(Verdicts)))
                grown(1:nVerdicts) = Verdicts(1:nVerdicts)
                call move_alloc(grown, Verdicts)
            end if
            nVerdicts = nVerdicts + 1
            if (nVerdicts == 1) then
                call ReadPwbVerdict(ud, Verdicts(1))
            else
                call ReadPwbVerdict(ud, Verdicts(nVerdicts), Verdicts(nVerdicts - 1))
            end if
        end do
        read(ud, iostat = io_status) magic
        close(ud)
        if (io_status /= 0 .or. magic /= VerdMagic) &
            error stop 'PWB verdict file is truncated.'
    end subroutine ReadProdVerdicts

    !> Production worker: the verdict for period p. A period the parent never
    !> classified is one this worker reached and a single pass did not - so
    !> something differs, and nothing computed from here on could be trusted.
    subroutine FindProdVerdict(p, v)
        integer, intent(in) :: p
        type(PwbVerdictType), intent(out) :: v
        integer :: i
        character(16) :: LogString

        do i = 1, nVerdicts
            if (Verdicts(i)%pcount == p) then
                v = Verdicts(i)
                return
            end if
        end do
        write(LogString, '(i0)') p
        call LogSay(' No PWB verdict for period ' // trim(LogString) &
            // ', which this worker reached.')
        error stop 'Production worker is out of step with its parent.'
    end subroutine FindProdVerdict

    character(PathLen) function VerdictPath(k)
        integer, intent(in) :: k
        character(32) :: tag

        write(tag, '(a,i2.2,a)') 'batch_pv_b', k, '.bin'
        if (BatchIndex > 0) then
            VerdictPath = trim(BatchTmpDir) // slash // trim(tag)
        else
            VerdictPath = trim(TmpDir) // trim(tag)
        end if
    end function VerdictPath

    !> A per-period folder under this worker's own output folder, moved to the
    !> same place under the parent's.
    subroutine ToParentTree(dir)
        character(*), intent(inout) :: dir
        integer :: n

        n = len_trim(WorkerMainOut)
        if (n == 0) return
        if (dir(1:min(len(dir), n)) /= WorkerMainOut(1:n)) return
        dir = trim(ParentMainOut) // dir(n + 1:len_trim(dir))
    end subroutine ToParentTree

    !> The context file sits in the parent's temporary folder, which a worker
    !> knows as --batch-tmp.
    character(PathLen) function CtxPath()
        if (BatchIndex > 0) then
            CtxPath = trim(BatchTmpDir) // slash // 'batch_pr_ctx.bin'
        else
            CtxPath = trim(TmpDir) // 'batch_pr_ctx.bin'
        end if
    end function CtxPath

end module m_production_parallel
