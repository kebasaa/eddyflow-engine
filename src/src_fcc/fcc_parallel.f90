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
        AppendBytes, FileBytes, MaxChunks
    use m_remote_source, only: RemoteFetchedAny
    implicit none
    private

    public :: CaptureFccContext, TryFccFluxSplit, FccWorkerStart
    public :: FccPieceBegins, FinishFccWorker, MergeFccPieces

    character(20), parameter :: CtxMagic = 'EDDYFLOW_FCCCTX_01  '
    character(20), parameter :: DumpMagic = 'EDDYFLOW_FCCOUT_01  '

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
