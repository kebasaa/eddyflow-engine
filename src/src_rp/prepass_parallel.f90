!***************************************************************************
! prepass_parallel.f90
! --------------------
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
! \brief       Splits the assessment pre-passes across worker processes.
!
!              The planar fit and the time-lag optimiser both walk every flux
!              averaging period in their range, reduce the raw data, and append
!              one record to a flat array. The fit itself - the sector
!              regressions, OptimizeTimelags - runs
!              once at the end and is cheap. So the loop is the cost, and the
!              loop's periods are independent of one another.
!
!              This module cuts that loop into contiguous pieces of similar
!              work - several per worker - and runs each in a copy of this
!              program, which processes its piece and writes the records it
!              produced to a file instead of going on to compute fluxes. A
!              fixed number run at once, and each that finishes makes room for
!              the next, so a fast core takes more pieces than a slow one. The
!              parent reads the pieces back in order and hands the
!              concatenation to the unchanged finalisation code.
!
!              Processes rather than threads because the period loop reaches
!              most of the program's global state - Stats, E2Col, Essentials,
!              the metadata - and making that thread-safe would be a rewrite.
!              A process has its own copy of all of it by construction.
!
!              The concatenation is exact: the records are fixed-size derived
!              types with no allocatable components, written unformatted and
!              read back by the same binary, and appended in piece order -
!              whatever order the pieces finished in.
!
!              The PWB cache pre-pass may be split too, now - see the call
!              site in eddyflow-rp_main.f90 for why splitting it stayed unsafe
!              for as long as it did: the streaming classifier's verdict used
!              to depend on the last settled detection before a period, a
!              chain with no time limit that a slice starting cold could not
!              rebuild. That cross-period state (the terminal fallback's lag,
!              the aggregate dataset's membership, the donor tally) has since
!              moved into PostProcessPwbTimelagCache, which the parent runs
!              once, serially, over every slice's output after they finish -
!              so what a worker produces is evidence for that pass rather than
!              a verdict of its own, and does not depend on where its slice
!              began.
!
! \author      Jonathan Muller
! \note
! \sa          batch_pool.f90 (the worker pool itself),
!              eddyflow-rp_main.f90, pwb_timelag_handle.f90
!***************************************************************************
module m_prepass_parallel
    use m_rp_global_var
    use m_batch_pool
    use m_pwb_timelag, only: AppendPwbCacheRows
    implicit none
    private

    !> The pool is m_batch_pool's; its names are passed on so that the
    !> program's calls did not have to change.
    public :: PlanPrepassBatches, PrepassChunk, PrepassChunkCount
    public :: StartPrepassBatches, WaitPrepassBatches, TopUpPrepassBatches
    public :: BatchDumpPath
    public :: FinishBatchWorker, StopIfParentGone
    public :: WriteTlagBatchDump, MergeTlagBatchDumps
    public :: WritePwbBatchDump, MergePwbBatchDumps
    public :: WritePfBatchDump, MergePfBatchDumps

    !> Guards the unformatted dumps. Parent and workers are the same binary in
    !> the same run, so the format never has to survive a version change - but
    !> a file a crashed earlier run left behind would otherwise be read as data.
    !>
    !> 05 adds the pre-pass kind to the header: a planar-fit worker once wrote
    !> time-lag records under the planar-fit name, and the parent read them as
    !> wind. Both ends now say which pre-pass a file belongs to.
    !>
    !> 06: TimeLagOptType carries one humidity per slot, not one scalar, so
    !> the records it dumps changed size.
    character(20), parameter :: BatchMagic = 'EDDYFLOW_PREPASS_06 '

contains

    !***************************************************************************
    !> \brief A worker's time-lag records, on their way back to the parent.
    !***************************************************************************
    subroutine WriteTlagBatchDump(dataset, nmax, n)
        integer, intent(in) :: nmax
        integer, intent(in) :: n
        type(TimeLagOptType), intent(in) :: dataset(nmax)
        integer :: u
        integer :: io_status

        call RequireBatchKind('to')

        open(newunit = u, file = trim(BatchOutPath), form = 'unformatted', &
            access = 'stream', status = 'replace', iostat = io_status)
        if (io_status /= 0) &
            error stop 'A pre-pass worker could not write its records.'

        write(u) BatchMagic
        write(u) BatchIndex, BatchCount
        write(u) BatchKind
        write(u) n
        if (n > 0) write(u) dataset(1:n)
        close(u)
    end subroutine WriteTlagBatchDump

    !***************************************************************************
    !> \brief Read the workers' pieces back, in order, as if one loop ran.
    !>
    !> Piece 1 is already in dataset(1:n) - the parent ran it itself - so this
    !> appends pieces 2..nChunks after it, which is the order the serial loop
    !> would have produced them in.
    !***************************************************************************
    subroutine MergeTlagBatchDumps(kind, nChunks, dataset, nmax, n)
        character(*), intent(in) :: kind
        integer, intent(in) :: nChunks
        integer, intent(in) :: nmax
        type(TimeLagOptType), intent(inout) :: dataset(nmax)
        integer, intent(inout) :: n
        integer :: k
        integer :: i
        integer :: u
        integer :: io_status
        integer :: nrec
        integer :: idx
        integer :: idxCount
        character(20) :: magic
        character(2) :: dumpKind
        type(TimeLagOptType), allocatable :: slice(:)

        do k = 2, nChunks
            open(newunit = u, file = trim(BatchDumpPath(kind, k)), &
                form = 'unformatted', access = 'stream', status = 'old', &
                iostat = io_status)
            if (io_status /= 0) &
                error stop 'A pre-pass worker record file could not be read.'

            read(u) magic
            if (magic /= BatchMagic) &
                error stop 'A pre-pass worker record file is not one of ours.'
            read(u) idx, idxCount
            read(u) dumpKind
            if (dumpKind /= kind) &
                error stop 'A pre-pass worker record file belongs to another pre-pass.'
            read(u) nrec

            if (nrec > 0) then
                allocate(slice(nrec))
                read(u) slice
                do i = 1, nrec
                    if (n >= nmax) &
                        error stop 'Merged pre-pass dataset is larger than the run allowed for.'
                    n = n + 1
                    dataset(n) = slice(i)
                end do
                deallocate(slice)
            end if
            close(u)
        end do
    end subroutine MergeTlagBatchDumps

    !***************************************************************************
    !> \brief A worker's PWB slice, on its way back to the parent.
    !>
    !> Two things travel, and neither can be derived from the other. The cache
    !> rows are the evidence the post-pass settles from - one per gas per
    !> period, fixed-size and with no allocatable component, so they go over
    !> unformatted exactly as the time-lag records do. The aggregate dataset
    !> carries only what the table cannot say: the humidity, and which period
    !> each row belongs to.
    !>
    !> The worker does NOT settle anything. Classification depends on having
    !> read the whole run, which is the reason a slice could not be trusted to
    !> classify its own periods in the first place.
    !***************************************************************************
    subroutine WritePwbBatchDump(dataset, nmax, nOpt)
        integer, intent(in) :: nmax
        integer, intent(in) :: nOpt
        type(TimeLagOptType), intent(in) :: dataset(nmax)
        integer :: u
        integer :: io_status

        call RequireBatchKind('to')

        open(newunit = u, file = trim(BatchOutPath), form = 'unformatted', &
            access = 'stream', status = 'replace', iostat = io_status)
        if (io_status /= 0) &
            error stop 'A pre-pass worker could not write its PWB records.'

        write(u) BatchMagic
        write(u) BatchIndex, BatchCount
        write(u) BatchKind
        write(u) PwbTimelagCacheN
        if (PwbTimelagCacheN > 0) write(u) PwbTimelagCache(1:PwbTimelagCacheN)
        write(u) nOpt
        if (nOpt > 0) then
            write(u) dataset(1:nOpt)
            write(u) PwbOptDate(1:nOpt)
            write(u) PwbOptTime(1:nOpt)
        end if
        close(u)
    end subroutine WritePwbBatchDump

    !***************************************************************************
    !> \brief Read the workers' PWB pieces back, in order, as if one loop ran.
    !>
    !> Piece 1 is already here - the parent ran it - so this appends
    !> 2..nChunks after it. Order is the whole point: the post-pass sorts the
    !> table by timestamp with a stable insertion sort, so rows appended in
    !> piece order
    !> come out exactly as a single loop would have left them, and periods
    !> sharing a timestamp keep their gas order.
    !>
    !> The rows go through AppendPwbCacheRows, the one way rows are added,
    !> so the table's capacity and its order flag stay true to its contents.
    !***************************************************************************
    subroutine MergePwbBatchDumps(kind, nChunks, dataset, nmax, nOpt)
        character(*), intent(in) :: kind
        integer, intent(in) :: nChunks
        integer, intent(in) :: nmax
        type(TimeLagOptType), intent(inout) :: dataset(nmax)
        integer, intent(inout) :: nOpt
        integer :: k, i, u, io_status
        integer :: nrec, idx, idxCount
        character(20) :: magic
        character(2) :: dumpKind
        type(PWBTimelagCacheEntryType), allocatable :: rows(:)
        type(TimeLagOptType), allocatable :: slice(:)
        character(10), allocatable :: sdate(:)
        character(5), allocatable :: stime(:)

        do k = 2, nChunks
            open(newunit = u, file = trim(BatchDumpPath(kind, k)), &
                form = 'unformatted', access = 'stream', status = 'old', &
                iostat = io_status)
            if (io_status /= 0) &
                error stop 'A pre-pass worker PWB file could not be read.'

            read(u) magic
            if (magic /= BatchMagic) &
                error stop 'A pre-pass worker PWB file is not one of ours.'
            read(u) idx, idxCount
            read(u) dumpKind
            if (dumpKind /= kind) &
                error stop 'A pre-pass worker record file belongs to another pre-pass.'

            read(u) nrec
            if (nrec > 0) then
                allocate(rows(nrec))
                read(u) rows
                call AppendPwbCacheRows(rows, nrec)
                deallocate(rows)
            end if

            read(u) nrec
            if (nrec > 0) then
                allocate(slice(nrec), sdate(nrec), stime(nrec))
                read(u) slice
                read(u) sdate
                read(u) stime
                do i = 1, nrec
                    if (nOpt >= nmax) &
                        error stop 'Merged PWB dataset is larger than the run allowed for.'
                    nOpt = nOpt + 1
                    dataset(nOpt) = slice(i)
                    PwbOptDate(nOpt) = sdate(i)
                    PwbOptTime(nOpt) = stime(i)
                end do
                deallocate(slice, sdate, stime)
            end if
            close(u)
        end do
    end subroutine MergePwbBatchDumps

    !***************************************************************************
    !> \brief A worker's planar-fit wind means, on their way back.
    !***************************************************************************
    subroutine WritePfBatchDump(wind, nmax, n)
        integer, intent(in) :: nmax
        integer, intent(in) :: n
        real(kind = dbl), intent(in) :: wind(nmax, 3)
        integer :: u
        integer :: io_status

        call RequireBatchKind('pf')

        open(newunit = u, file = trim(BatchOutPath), form = 'unformatted', &
            access = 'stream', status = 'replace', iostat = io_status)
        if (io_status /= 0) &
            error stop 'A pre-pass worker could not write its records.'

        write(u) BatchMagic
        write(u) BatchIndex, BatchCount
        write(u) BatchKind
        write(u) n
        if (n > 0) write(u) wind(1:n, 1:3)
        close(u)
    end subroutine WritePfBatchDump

    !***************************************************************************
    !> \brief Read the planar-fit pieces back, in order, after the parent's own.
    !***************************************************************************
    subroutine MergePfBatchDumps(nChunks, wind, nmax, n)
        integer, intent(in) :: nChunks
        integer, intent(in) :: nmax
        real(kind = dbl), intent(inout) :: wind(nmax, 3)
        integer, intent(inout) :: n
        integer :: k
        integer :: u
        integer :: i
        integer :: io_status
        integer :: nrec
        integer :: idx
        integer :: idxCount
        character(20) :: magic
        character(2) :: dumpKind
        real(kind = dbl), allocatable :: slice(:, :)

        do k = 2, nChunks
            open(newunit = u, file = trim(BatchDumpPath('pf', k)), &
                form = 'unformatted', access = 'stream', status = 'old', &
                iostat = io_status)
            if (io_status /= 0) &
                error stop 'A pre-pass worker record file could not be read.'

            read(u) magic
            if (magic /= BatchMagic) &
                error stop 'A pre-pass worker record file is not one of ours.'
            read(u) idx, idxCount
            read(u) dumpKind
            if (dumpKind /= 'pf') &
                error stop 'A pre-pass worker record file belongs to another pre-pass.'
            read(u) nrec
            if (nrec > 0) then
                allocate(slice(nrec, 3))
                read(u) slice
                do i = 1, nrec
                    if (n >= nmax) &
                        error stop 'Merged pre-pass dataset is larger than the run allowed for.'
                    n = n + 1
                    wind(n, 1:3) = slice(i, 1:3)
                end do
                deallocate(slice)
            end if
            close(u)
        end do
    end subroutine MergePfBatchDumps

end module m_prepass_parallel
