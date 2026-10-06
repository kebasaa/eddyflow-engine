!***************************************************************************
! fcc_parallel.f90
! ----------------
! Copyright © 2026-    , ETH Zurich, Jonathan Muller
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
! \brief       Splits FCC's flux computation across worker processes.
!
!              FCC corrects one essentials record at a time, and a record's
!              fluxes depend on nothing but the record and the spectral
!              assessment done before the loop. So after the first record that
!              is written - the head, which opens the output files - the rest
!              of the file is cut into pieces and handed to copies of this
!              program (m_batch_pool), the parent taking the first itself.
!
!              A worker does not repeat the spectral assessment: it reads what
!              the parent's left - the assessment and the settings it changed -
!              from a context file. It processes the head again, into a folder
!              of its own, so its output files open exactly as the parent's
!              did; skips to its piece without parsing the records between;
!              and processes its piece. Then it says where in each output file
!              its piece's rows begin, and hands back the periods it kept for
!              post-flux despiking. The parent appends those rows to its own
!              files piece by piece in record order, so the files hold what a
!              single pass would have written, and the end of the run - the
!              despiking, the datasets - runs on them unchanged.
!
!              Refused, so the loop runs on serially: embedded mode, and inputs
!              from a shared link, which every worker would download again.
!
!              The binned (co)spectra the spectral assessment imports are read
!              by workers too ('sb'), but only read: a worker parses its files
!              and hands back their bins. The parent runs the import loop as it
!              always did, taking each file's bins from the workers in file
!              order instead of reading the file - so the lookup in the
!              essentials file, the quality screening and the sums, which
!              depend on the order the files come in, are its own and
!              unchanged. It reads the first piece itself while the workers
!              read the rest.
!
! \author      Jonathan Muller
! \sa          batch_pool.f90, eddyflow-fcc_main.f90,
!              production_parallel.f90 (the same pattern in RP)
!***************************************************************************
module m_fcc_parallel
    use m_fx_global_var
    use m_sa_rates
    use m_batch_pool, only: PlanPrepassBatches, StartPrepassBatches, &
        WaitPrepassBatches, PrepassChunkCount, BatchDumpPath, ForcedPieceLength, &
        RemoveStaleWorkerRoots, WorkerRoot, RemoveWorkerRoots, NoSlashEnd, &
        AppendBytes, FileBytes, MaxChunks, WaitForBatchPiece, CancelBatchPool, &
        PrepassChunk, FinishBatchWorker, StopIfParentGone
    use m_remote_source, only: RemoteFetchedAny
    implicit none
    private

    public :: CaptureFccContext, TryFccFluxSplit, FccWorkerStart
    public :: FccPieceBegins, FinishFccWorker, MergeFccPieces
    public :: StartBinnedSplit, GetBinnedFile, FinishBinnedSplit, RunBinnedReadWorker

    character(20), parameter :: CtxMagic = 'EDDYFLOW_FCCCTX_01  '
    character(20), parameter :: DumpMagic = 'EDDYFLOW_FCCOUT_01  '
    character(20), parameter :: BinMagic = 'EDDYFLOW_FCCBIN_01  '

    !> A binned-import piece holds this many files at least. One piece per
    !> worker: the files are much alike, and every further piece would pay
    !> for starting a process again.
    integer, parameter :: MinPieceFiles = 16

    !> A worker piece pays for starting a process and reading the essentials
    !> file once more, so it should hold this many records at least.
    integer, parameter :: MinPieceRecords = 64

    !> The continuous output files the flux loop writes. Named one by one: the
    !> essentials file being read, uex, lies between them in unit numbers.
    integer, parameter :: NOut = 3
    integer, parameter :: OutUnits(NOut) = [uflx, umd, uflxnt]

    !> What the flux loop starts from, kept when it starts.
    type(EddyFlowProjType) :: CtxProj
    type(FCCsetupType) :: CtxFCCsetup
    type(RegParType) :: CtxRegPar(GHGNumVar, MaxGasClasses)
    real(kind = dbl) :: CtxStPar(2)
    real(kind = dbl) :: CtxUnPar(2)
    character(32) :: CtxTFShape
    type(MassParType) :: CtxMassPar(GHGNumVar, 2)

    !> A worker's piece, as it began.
    logical :: PieceOpen(NOut) = .false.
    integer(8) :: PieceOffset(NOut) = 0
    character(PathLen) :: PiecePath(NOut) = ''

    !> The parent's side of a split binned import: whether one is under way,
    !> which piece it is reading, and that piece's dump.
    logical :: SbActive = .false.
    integer :: SbPiece = 0
    integer :: SbUnit = 0
    logical :: SbOpen = .false.

contains

    !***************************************************************************
    !> \brief Parent: keep what the flux loop starts from, once the spectral
    !>        assessment is done.
    !***************************************************************************
    subroutine CaptureFccContext()
        CtxProj = EddyFlowProj
        CtxFCCsetup = FCCsetup
        CtxRegPar = RegPar
        CtxStPar = StPar
        CtxUnPar = UnPar
        CtxTFShape = TFShape
        CtxMassPar = MassPar
    end subroutine CaptureFccContext

    !***************************************************************************
    !> \brief Parent, once the head record has opened the output files: cut
    !>        the records after it into pieces and start the workers.
    !>
    !> parentEnd is the record the parent's own piece ends before; 0 if the
    !> loop is not split. Records are cut by count: a record costs about the
    !> same as another, and those past the end of a selected period are left
    !> at once.
    !***************************************************************************
    subroutine TryFccFluxSplit(head, nRecords, parentEnd, nEff)
        integer, intent(in) :: head
        integer, intent(in) :: nRecords
        integer, intent(out) :: parentEnd
        integer, intent(out) :: nEff
        integer :: nLeft
        integer :: forced
        integer :: want
        integer :: k
        integer :: n
        integer :: cuts(MaxChunks)
        type(DateType) :: noSeries(1)
        type(FileListType) :: noFiles(1)

        parentEnd = 0
        nEff = 1
        if (EddyFlowProj%run_env == 'embedded') then
            call LogSay('  The flux computation is not split: it runs in embedded mode.')
            return
        end if
        if (RemoteFetchedAny()) then
            call LogSay('  The flux computation is not split: its inputs come from a shared link.')
            return
        end if

        !> [head + 1, nRecords + 1), half-open like the pre-passes' ranges
        nLeft = nRecords - head
        forced = ForcedPieceLength('EDDYFLOW_FCC_PIECE_PERIODS')
        if (forced > 0) then
            call PlanPrepassBatches(nLeft, .true., nEff, 'flux computation', 1)
        else
            call PlanPrepassBatches(nLeft, .true., nEff, 'flux computation', MinPieceRecords)
        end if
        if (nEff <= 1) return

        if (forced > 0) then
            want = min(MaxChunks, max(1, nLeft / forced))
        else
            want = nEff
        end if
        if (want <= 1) return
        n = 0
        do k = 1, want - 1
            n = n + 1
            cuts(n) = head + 1 + int(int(nLeft, 8) * k / want)
            if (n > 1) then
                if (cuts(n) <= cuts(n - 1)) n = n - 1
            end if
        end do
        if (n == 0) return

        call RemoveStaleWorkerRoots()
        call WriteFccContext(head)
        call StartPrepassBatches('fx', head + 1, nRecords + 1, nEff, noSeries, 1, &
            noFiles, 1, cuts(1:n))
        parentEnd = cuts(1)
    end subroutine TryFccFluxSplit

    !***************************************************************************
    !> \brief Parent: write the context the workers start from.
    !***************************************************************************
    subroutine WriteFccContext(head)
        integer, intent(in) :: head
        integer :: u
        integer :: io_status

        open(newunit = u, file = trim(CtxPath()), form = 'unformatted', &
            access = 'stream', status = 'replace', iostat = io_status)
        if (io_status /= 0) error stop 'Could not write the flux computation context file.'
        write(u) CtxMagic
        write(u) Timestamp_FilePadding
        write(u) head
        write(u) CtxProj, CtxFCCsetup, CtxRegPar, CtxStPar, CtxUnPar, CtxTFShape, CtxMassPar
        !> m_sa_rates: the assessment per acquisition rate
        write(u) nRateConfigs, ConfigFileRate, ConfigGasRate, ConfigPeriods
        write(u) nGasRates, GasRate, GasRatePeriods, nFileRates, FileRate
        write(u) MultiRateSA, SAPassConfig
        write(u) RegParR, SlotCnt, SlotFilled, GasRateUse
        write(u) UnParR, StParR, FileSlotFilled, FileRateUse
        write(u) CtxMagic
        close(u)
    end subroutine WriteFccContext

    !***************************************************************************
    !> \brief Worker, in place of the spectral assessment: take up the parent's
    !>        results, and write into a folder of its own.
    !***************************************************************************
    subroutine FccWorkerStart(head)
        integer, intent(out) :: head
        integer :: u
        integer :: io_status
        integer :: mkdir_status
        character(20) :: magic
        character(len(EddyFlowProj%caller)) :: ownCaller
        integer, external :: CreateDir

        open(newunit = u, file = trim(CtxPath()), form = 'unformatted', &
            access = 'stream', status = 'old', action = 'read', iostat = io_status)
        if (io_status /= 0) error stop 'Flux computation worker could not open its context file.'
        read(u, iostat = io_status) magic
        if (io_status /= 0 .or. magic /= CtxMagic) &
            error stop 'Flux computation context file is not one of this run.'
        ownCaller = EddyFlowProj%caller
        read(u) Timestamp_FilePadding
        read(u) head
        read(u) EddyFlowProj, FCCsetup, RegPar, StPar, UnPar, TFShape, MassPar
        read(u) nRateConfigs, ConfigFileRate, ConfigGasRate, ConfigPeriods
        read(u) nGasRates, GasRate, GasRatePeriods, nFileRates, FileRate
        read(u) MultiRateSA, SAPassConfig
        read(u) RegParR, SlotCnt, SlotFilled, GasRateUse
        read(u) UnParR, StParR, FileSlotFilled, FileRateUse
        read(u, iostat = io_status) magic
        close(u)
        if (io_status /= 0 .or. magic /= CtxMagic) &
            error stop 'Flux computation context file is truncated.'
        !> Its own: a worker reports as a console caller whatever started the run
        EddyFlowProj%caller = ownCaller

        Dir%main_out = WorkerRoot(BatchParentPid, 'fx', BatchIndex)
        BatchOwnOutDir = Dir%main_out
        !> Empty before it is used - see AdoptProdWorkerOutput.
        io_status = system(trim(comm_rmdir) // ' "' &
            // trim(NoSlashEnd(Dir%main_out)) // '"' // comm_err_redirect)
        mkdir_status = CreateDir('"' // trim(Dir%main_out) // '"')
    end subroutine FccWorkerStart

    !***************************************************************************
    !> \brief Worker: its piece begins. Note where each output file stands -
    !>        what follows is the piece's - and start the despiking cache
    !>        afresh.
    !***************************************************************************
    subroutine FccPieceBegins()
        integer :: j
        logical :: op
        character(PathLen) :: nm

        PfdCacheN = 0
        do j = 1, NOut
            inquire(unit = OutUnits(j), opened = op, name = nm)
            PieceOpen(j) = op
            PieceOffset(j) = 0
            PiecePath(j) = ''
            if (.not. op) cycle
            flush(OutUnits(j))
            PiecePath(j) = nm
            PieceOffset(j) = FileBytes(nm)
        end do
    end subroutine FccPieceBegins

    !***************************************************************************
    !> \brief Worker: its piece is done. Hand back where its rows begin in each
    !>        file and the periods it kept for despiking.
    !***************************************************************************
    subroutine FinishFccWorker()
        integer :: j
        integer :: ud
        integer :: io_status
        logical :: op

        do j = 1, NOut
            inquire(unit = OutUnits(j), opened = op)
            if (op .and. .not. PieceOpen(j)) then
                call LogSay(' A flux computation worker opened an output file during')
                call LogSay(' its piece that was not open when the piece began.')
                error stop 'Flux computation worker output cannot be merged.'
            end if
            if (op) close(OutUnits(j))
        end do

        open(newunit = ud, file = trim(BatchOutPath), form = 'unformatted', &
            access = 'stream', status = 'replace', iostat = io_status)
        if (io_status /= 0) error stop 'Flux computation worker could not write its dump.'
        write(ud) DumpMagic
        write(ud) count(PieceOpen)
        do j = 1, NOut
            if (.not. PieceOpen(j)) cycle
            write(ud) OutUnits(j), PieceOffset(j), PiecePath(j)
        end do
        write(ud) PfdCacheN
        if (PfdCacheN > 0) write(ud) PfdCache(1:PfdCacheN)
        write(ud) DumpMagic
        close(ud)
    end subroutine FinishFccWorker

    !***************************************************************************
    !> \brief Parent: append every worker's rows to its own files, in piece
    !>        order, and its despiking periods to the parent's.
    !***************************************************************************
    subroutine MergeFccPieces()
        integer :: k
        integer :: j
        integer :: i
        integer :: n
        integer :: u
        integer :: ud
        integer :: io_status
        integer :: wN
        integer(8) :: offset
        logical :: op
        logical :: seen(NOut)
        logical :: parentOpen(NOut)
        character(PathLen) :: parentPath(NOut)
        character(PathLen) :: path
        character(20) :: magic
        type(PfdCacheEntryType), allocatable :: rows(:)

        do j = 1, NOut
            inquire(unit = OutUnits(j), opened = op, name = path)
            parentOpen(j) = op
            parentPath(j) = ''
            if (.not. op) cycle
            parentPath(j) = path
            close(OutUnits(j))
        end do

        do k = 2, PrepassChunkCount()
            open(newunit = ud, file = trim(BatchDumpPath('fx', k)), &
                form = 'unformatted', access = 'stream', status = 'old', &
                action = 'read', iostat = io_status)
            if (io_status /= 0) error stop 'Could not read a flux computation worker dump.'
            read(ud, iostat = io_status) magic
            if (io_status /= 0 .or. magic /= DumpMagic) &
                error stop 'A flux computation worker dump is not one of this run.'
            read(ud) n
            seen = .false.
            do i = 1, n
                read(ud) u, offset, path
                j = findloc(OutUnits, u, dim = 1)
                if (j == 0) error stop 'A flux computation worker dump names an unknown file.'
                if (.not. parentOpen(j)) then
                    call LogSay(' A flux computation worker wrote ' // trim(path))
                    call LogSay(' which this process never opened.')
                    error stop 'Flux computation worker output cannot be merged.'
                end if
                seen(j) = .true.
                call AppendBytes(path, offset, parentPath(j))
            end do
            do j = 1, NOut
                if (parentOpen(j) .and. .not. seen(j)) then
                    call LogSay(' A flux computation worker did not write ' // trim(parentPath(j)))
                    error stop 'Flux computation worker output cannot be merged.'
                end if
            end do
            read(ud) wN
            if (wN > 0) then
                allocate(rows(wN))
                read(ud) rows
                do i = 1, wN
                    call StorePfdCache(rows(i)%date, rows(i)%time, rows(i)%nee, &
                        rows(i)%h, rows(i)%le)
                end do
                deallocate(rows)
            end if
            read(ud, iostat = io_status) magic
            close(ud)
            if (io_status /= 0 .or. magic /= DumpMagic) &
                error stop 'A flux computation worker dump is truncated.'
        end do
        call RemoveWorkerRoots('fx')
    end subroutine MergeFccPieces

    !***************************************************************************
    !> \brief Parent, before the import loop: cut the selected binned files
    !>        into pieces and start the workers reading them.
    !***************************************************************************
    subroutine StartBinnedSplit(first, last, Files, nFiles)
        integer, intent(in) :: first
        integer, intent(in) :: last
        integer, intent(in) :: nFiles
        type(FileListType), intent(in) :: Files(nFiles)
        integer :: nEff
        integer :: n
        integer :: c
        integer :: k
        integer :: want
        integer :: forced
        integer :: u
        integer :: io_status
        integer :: cuts(MaxChunks)
        type(DateType) :: noSeries(1)

        SbActive = .false.
        SbPiece = 0
        SbOpen = .false.
        if (BatchIndex > 0) return
        if (EddyFlowProj%run_env == 'embedded') return
        if (RemoteFetchedAny()) return

        n = last - first + 1
        forced = ForcedPieceLength('EDDYFLOW_FCC_PIECE_FILES')
        if (forced > 0) then
            call PlanPrepassBatches(n, .true., nEff, 'binned (co)spectra import', 1)
        else
            call PlanPrepassBatches(n, .true., nEff, 'binned (co)spectra import', MinPieceFiles)
        end if
        if (nEff <= 1) return
        if (forced > 0) then
            want = min(MaxChunks, max(1, n / forced))
        else
            want = min(MaxChunks, nEff, max(1, n / MinPieceFiles))
        end if
        if (want <= 1) return
        k = 0
        do c = 1, want - 1
            k = k + 1
            cuts(k) = first + int(int(n, 8) * c / want)
            if (k > 1) then
                if (cuts(k) <= cuts(k - 1)) k = k - 1
            end if
        end do
        if (k == 0) return

        !> The list the workers read from, so they need not list the folder
        open(newunit = u, file = trim(BinListPath()), form = 'unformatted', &
            access = 'stream', status = 'replace', iostat = io_status)
        if (io_status /= 0) error stop 'Could not write the binned file list for the workers.'
        write(u) BinMagic, first, last
        write(u) Files(first:last)
        close(u)

        call StartPrepassBatches('sb', first, last + 1, nEff, noSeries, 1, &
            Files, nFiles, cuts(1:k))
        SbActive = .true.
        SbPiece = 1
    end subroutine StartBinnedSplit

    !***************************************************************************
    !> \brief Parent, in place of ReadBinnedFile: file fcount's bins, read here
    !>        in the first piece and taken from a worker's dump after it.
    !>
    !> Called for every file of the loop in turn, as ReadBinnedFile was. Like
    !> it, leaves nbins alone for a file that could not be read, and says so
    !> with Error(62) - the worker's own message stays in its log.
    !***************************************************************************
    subroutine GetBinnedFile(fcount, InFile, BinSpec, BinCosp, nrow, nbins, skip)
        integer, intent(in) :: fcount
        type(FileListType), intent(in) :: InFile
        integer, intent(in) :: nrow
        type(SpectraSetType), intent(inout) :: BinSpec(nrow)
        type(SpectraSetType), intent(inout) :: BinCosp(nrow)
        integer, intent(inout) :: nbins
        logical, intent(out) :: skip
        integer :: s
        integer :: e
        integer :: fc
        integer :: nb
        logical :: more
        logical :: sk

        if (.not. SbActive) then
            call ReadBinnedFile(InFile, BinSpec, BinCosp, nrow, nbins, skip)
            return
        end if
        call PrepassChunk(1, s, e)
        if (fcount < e) then
            call ReadBinnedFile(InFile, BinSpec, BinCosp, nrow, nbins, skip)
            return
        end if

        !> The piece this file is in; the pieces before it are done with
        do
            call PrepassChunk(SbPiece, s, e)
            if (fcount < e .and. SbPiece > 1) exit
            if (SbOpen) close(SbUnit, status = 'delete')
            SbOpen = .false.
            SbPiece = SbPiece + 1
        end do
        if (.not. SbOpen) call OpenBinDump(SbPiece)

        read(SbUnit) more, fc, sk, nb
        if (.not. more .or. fc /= fcount) &
            error stop 'A binned-import worker dump is out of step with the import.'
        skip = sk
        if (skip) then
            call ExceptionHandler(62)
            return
        end if
        nbins = nb
        call ReadBins(SbUnit, BinSpec, nrow, nbins)
        call ReadBins(SbUnit, BinCosp, nrow, nbins)
    end subroutine GetBinnedFile

    !> Wait for piece k and open its dump, checking it was written for the
    !> gases this process expects.
    subroutine OpenBinDump(k)
        integer, intent(in) :: k
        integer :: io_status
        character(20) :: magic
        character(64) :: mine(GHGNumVar)
        character(64) :: theirs(GHGNumVar)
        include '../src_common/interfaces_1.inc'

        call WaitForBatchPiece(k)
        open(newunit = SbUnit, file = trim(BatchDumpPath('sb', k)), &
            form = 'unformatted', access = 'stream', status = 'old', &
            action = 'readwrite', iostat = io_status)
        if (io_status /= 0) error stop 'Could not read a binned-import worker dump.'
        read(SbUnit, iostat = io_status) magic, theirs
        if (io_status /= 0 .or. magic /= BinMagic) &
            error stop 'A binned-import worker dump is not one of this run.'
        call SpectralVarTags(mine)
        if (any(mine /= theirs)) &
            error stop 'A binned-import worker read the files for other gases.'
        SbOpen = .true.
    end subroutine OpenBinDump

    !***************************************************************************
    !> \brief Parent, after the import loop: let the workers still reading
    !>        finish - the loop may have stopped early, at the end of the
    !>        essentials file - and remove what is left of their dumps.
    !***************************************************************************
    subroutine FinishBinnedSplit()
        integer :: k
        integer :: u
        integer :: io_status
        logical :: ex

        if (.not. SbActive) return
        if (SbOpen) close(SbUnit, status = 'delete')
        SbOpen = .false.
        call CancelBatchPool()
        do k = 2, PrepassChunkCount()
            inquire(file = trim(BatchDumpPath('sb', k)), exist = ex)
            if (.not. ex) cycle
            open(newunit = u, file = trim(BatchDumpPath('sb', k)), status = 'old', &
                iostat = io_status)
            if (io_status == 0) close(u, status = 'delete')
        end do
        SbActive = .false.
    end subroutine FinishBinnedSplit

    !***************************************************************************
    !> \brief Binned-import worker: read its piece of the files and stop.
    !***************************************************************************
    subroutine RunBinnedReadWorker()
        integer :: u
        integer :: ud
        integer :: io_status
        integer :: first
        integer :: last
        integer :: fcount
        integer :: nbins
        logical :: skip
        character(20) :: magic
        character(64) :: tags(GHGNumVar)
        type(FileListType), allocatable :: Files(:)
        type(SpectraSetType) :: BinSpec(MaxNumBins)
        type(SpectraSetType) :: BinCosp(MaxNumBins)
        include '../src_common/interfaces_1.inc'

        open(newunit = u, file = trim(BinListPath()), form = 'unformatted', &
            access = 'stream', status = 'old', action = 'read', iostat = io_status)
        if (io_status /= 0) error stop 'Binned-import worker could not open its file list.'
        read(u, iostat = io_status) magic, first, last
        if (io_status /= 0 .or. magic /= BinMagic) &
            error stop 'Binned-import file list is not one of this run.'
        allocate(Files(first:last))
        read(u) Files
        close(u)
        if (BatchSliceStart < first .or. BatchSliceEnd - 1 > last) &
            error stop 'Binned-import worker piece lies outside the file list.'

        open(newunit = ud, file = trim(BatchOutPath), form = 'unformatted', &
            access = 'stream', status = 'replace', iostat = io_status)
        if (io_status /= 0) error stop 'Binned-import worker could not write its dump.'
        call SpectralVarTags(tags)
        write(ud) BinMagic, tags
        nbins = 0
        do fcount = BatchSliceStart, BatchSliceEnd - 1
            call StopIfParentGone()
            call ReadBinnedFile(Files(fcount), BinSpec, BinCosp, MaxNumBins, nbins, skip)
            write(ud) .true., fcount, skip, nbins
            if (skip) cycle
            call WriteBins(ud, BinSpec, MaxNumBins, nbins)
            call WriteBins(ud, BinCosp, MaxNumBins, nbins)
        end do
        write(ud) .false., 0, .true., 0
        close(ud)
        call FinishBatchWorker()
        stop
    end subroutine RunBinnedReadWorker

    !> One file's bins, rows 1:nbins and gas slots up to the last with a value:
    !> the rest of the array is ErrSpec, which is what ReadBinnedFile leaves.
    subroutine WriteBins(u, Bins, nrow, nbins)
        integer, intent(in) :: u
        integer, intent(in) :: nrow
        integer, intent(in) :: nbins
        type(SpectraSetType), intent(in) :: Bins(nrow)
        integer :: i
        integer :: j
        integer :: hi

        hi = 0
        do i = 1, nbins
            do j = size(Bins(i)%of), hi + 1, -1
                if (Bins(i)%of(j) /= error) then
                    hi = j
                    exit
                end if
            end do
        end do
        write(u) hi
        do i = 1, nbins
            write(u) Bins(i)%fnum, Bins(i)%fn, Bins(i)%fnorm
            if (hi > 0) write(u) Bins(i)%of(1:hi)
        end do
    end subroutine WriteBins

    subroutine ReadBins(u, Bins, nrow, nbins)
        integer, intent(in) :: u
        integer, intent(in) :: nrow
        integer, intent(in) :: nbins
        type(SpectraSetType), intent(inout) :: Bins(nrow)
        integer :: i
        integer :: hi

        Bins = ErrSpec
        read(u) hi
        do i = 1, nbins
            read(u) Bins(i)%fnum, Bins(i)%fn, Bins(i)%fnorm
            if (hi > 0) read(u) Bins(i)%of(1:hi)
        end do
    end subroutine ReadBins

    character(PathLen) function BinListPath()
        if (BatchIndex > 0) then
            BinListPath = trim(BatchTmpDir) // slash // 'batch_sb_list.bin'
        else
            BinListPath = trim(TmpDir) // 'batch_sb_list.bin'
        end if
    end function BinListPath

    !> The context file sits in the parent's temporary folder, which a worker
    !> knows as --batch-tmp.
    character(PathLen) function CtxPath()
        if (BatchIndex > 0) then
            CtxPath = trim(BatchTmpDir) // slash // 'batch_fx_ctx.bin'
        else
            CtxPath = trim(TmpDir) // 'batch_fx_ctx.bin'
        end if
    end function CtxPath

end module m_fcc_parallel
