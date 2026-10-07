!***************************************************************************
! batch_pool.f90
! --------------
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
! \brief       Runs a range of averaging periods - or files - as pieces, each in
!              a worker process, a fixed number at a time.
!
!              Shared by RP's pre-passes and production pass and by FCC. Cuts
!              a range into pieces, writes one script per piece that starts a
!              copy of the program on it (--batch), keeps a fixed number of
!              them running as they finish, and stops the run as soon as one
!              fails. A worker watches the process that started it and stops
!              when it has gone. What a piece produces, and how the parent
!              puts the pieces back together, is the caller's: see
!              prepass_parallel.f90, production_parallel.f90 and
!              fcc_parallel.f90.
!
! \author      Jonathan Muller
! \note
! \sa          prepass_parallel.f90, production_parallel.f90
!***************************************************************************
module m_batch_pool
    use m_common_global_var
    use m_ghg_prefetch, only: GhgPrefetchCleanup
    use m_remote_source, only: RemoteCleanup
    use m_process_os, only: ProcessSelfId, ParentGone, ProcessAlive
    implicit none
    private

    public :: PlanPrepassBatches, PlanPrepassChunks, PrepassChunk, PrepassChunkCount
    public :: StartPrepassBatches, WaitPrepassBatches, TopUpPrepassBatches
    public :: WaitForBatchPiece, CancelBatchPool
    public :: BatchDumpPath, RequireBatchKind, NoTrailingSlash
    public :: FinishBatchWorker, StopIfParentGone
    public :: ForcedPieceLength, RemoveStaleWorkerRoots, WorkerRoot
    public :: RemoveWorkerRoots, NoSlashEnd, AppendBytes, FileBytes
    public :: MaxChunks

    !> More workers than this is never a throughput win on a machine that also
    !> has to feed them raw data, and it multiplies the per-worker cost of
    !> listing the raw directory.
    integer, parameter :: MaxWorkers = 32

    !> Roughly one second per tick, so this is a full day of waiting. A
    !> pre-pass over a season legitimately takes hours; a worker that has hung
    !> should still not hold the parent for ever.
    integer, parameter :: MaxWaitTicks = 86400


    !> The range is cut into about this many pieces per worker, so a worker
    !> that finishes early - on a performance core, or over a stretch with no
    !> data - takes the next piece instead of idling until the slowest is done.
    integer, parameter :: ChunksPerWorker = 4

    !> Every per-piece file name carries the piece number in two digits.
    integer, parameter :: MaxChunks = 99

    !> What walking a period with no raw file costs, against 1 for a period
    !> that is read and reduced: it imports nothing, so almost nothing.
    real(kind = dbl), parameter :: EmptyPeriodWeight = 0.05d0

    !> The current pre-pass's pieces, half-open like the loops: piece k is
    !> [ChunkStart(k), ChunkEnd(k)). One pre-pass runs at a time, so one plan.
    integer :: NumChunks = 0
    integer, allocatable :: ChunkStart(:)
    integer, allocatable :: ChunkEnd(:)

    !> The next piece no worker has been given yet.
    integer :: NextChunk = 0

    !> The pool of the pass under way: which pass, how many processes may run
    !> at once (the parent included), which pieces a worker is running, and
    !> how many pieces are done - the parent's own counted once it waits.
    logical :: PoolActive = .false.
    character(2) :: PoolKind = ''
    integer :: PoolEff = 0
    integer :: PoolDone = 0
    logical, allocatable :: Running(:)

contains

    !***************************************************************************
    !> \brief How many workers to use.
    !>
    !> Returns nEff = 1 for "run it serially, as always". The caller does not
    !> have to special-case that: it simply takes the loop it already had.
    !>
    !> `allowed` is the caller's judgement that this particular pre-pass may
    !> be split at all. The PWB cache pre-pass may not: its classifier decides
    !> a period partly from the last settled detection before it, that chain
    !> has no time limit, and on real data the weak species never settle - COS
    !> reached no settled detection at all over two days of CH-LAE, so no
    !> lead-in of any length would rebuild its state. A split there cannot
    !> reproduce a single pass, so it is not offered.
    !***************************************************************************
    subroutine PlanPrepassBatches(nPeriods, allowed, nEff, what, minPerWorker)
        integer, intent(in) :: nPeriods
        logical, intent(in) :: allowed
        integer, intent(out) :: nEff
        !> What is being split, for the log line; the pre-pass by default.
        character(*), intent(in), optional :: what
        !> Fewest periods per worker that make a split worth it; 4 by default.
        integer, intent(in), optional :: minPerWorker
        integer :: requested
        integer :: perWorker
        character(64) :: LogString
        character(4096) :: cmdline

        nEff = 1

        !> A worker never splits again. Without this a worker would spawn its
        !> own workers, each of which would spawn more - and the growth is
        !> exponential, so it saturates the machine in seconds.
        if (BatchIndex > 0) return

        !> The same question asked of the raw command line rather than of the
        !> parsed result, because the parsed result is exactly what fails when
        !> something is wrong with the switch handling. A process carrying
        !> --batch that nevertheless believes it is a parent is one that
        !> failed to parse its own instructions, and it is about to start a
        !> fork bomb. Stop it here instead, where the cause is still legible.
        call get_command(cmdline)
        if (index(cmdline, '--batch') > 0) then
            call LogSay(' This process was launched as a pre-pass worker but did')
            call LogSay(' not recognise its own --batch argument.')
            call LogSay(' Command line: ' // trim(cmdline))
            error stop 'Pre-pass worker could not read its batch assignment.'
        end if

        !> Nothing that carries state between periods is split.
        if (.not. allowed) return

        requested = NumJobs
        if (requested <= 0) requested = DetectCoreCount()
        requested = min(requested, MaxWorkers)
        if (requested <= 1) return

        !> A slice has to be worth the cost of starting a process and listing
        !> the raw directory again. Below that the split loses to the serial
        !> loop, so do not make one.
        perWorker = 4
        if (present(minPerWorker)) perWorker = max(1, minPerWorker)
        nEff = max(1, min(requested, nPeriods / perWorker))
        if (nEff <= 1) then
            nEff = 1
            return
        end if

        write(LogString, '(i6)') nEff
        if (present(what)) then
            call LogSay('  Splitting the ' // trim(what) // ' across ' &
                // trim(adjustl(LogString)) // ' worker processes.')
        else
            call LogSay('  Splitting the pre-pass across ' &
                // trim(adjustl(LogString)) // ' worker processes.')
        end if
    end subroutine PlanPrepassBatches

    !***************************************************************************
    !> \brief Cut [iStart, iEnd) into pieces of about equal work.
    !>
    !> The range used to be cut into one slice per worker by period count. On
    !> the Yatir run that gave one slice 801 raw files and another none - the
    !> data thin out from June and stop for two weeks in July - and the run
    !> waited on the heaviest. So the cut is by work instead: each period
    !> weighs 1 if a raw file covers it and EmptyPeriodWeight if none does,
    !> and the pieces hold equal shares of the total. There are several per
    !> worker, handed out as workers come free, which takes care of what the
    !> weights cannot see: a core that is slower than another.
    !>
    !> The weight is only an estimate of cost - the file's time span from
    !> its name, not whether it holds usable records - and it decides nothing
    !> but where the cuts fall. Any cut is correct: the pieces tile the range
    !> and are merged back in order.
    !>
    !> The range is HALF-OPEN, [iStart, iEnd), because that is what the period
    !> loops themselves do: both exit on `pcount >= endIndex`, so the end index
    !> is one past the last period processed, and period p runs from
    !> Series(p) to Series(p + 1). Piece ends are exclusive for the same reason
    !> and go straight into the loop's end index.
    !>
    !> Piece 1 is the parent's, and the parent has to read at least one raw
    !> file in it: the code after the loop reads state that reading a file
    !> establishes (see StartPrepassBatches). So piece 1 always reaches the
    !> first period that has a file, however long the empty stretch before it.
    !***************************************************************************
    subroutine PlanPrepassChunks(iStart, iEnd, nEff, Series, nSeries, Files, nFiles)
        integer, intent(in) :: iStart
        integer, intent(in) :: iEnd
        integer, intent(in) :: nEff
        integer, intent(in) :: nSeries
        integer, intent(in) :: nFiles
        type(DateType), intent(in) :: Series(nSeries)
        type(FileListType), intent(in) :: Files(nFiles)
        integer :: p
        integer :: j
        integer :: c
        integer :: want
        integer :: firstData
        integer :: nCuts
        !> One fewer cut than pieces, at most.
        integer :: cuts(MaxChunks)
        real(kind = dbl) :: total
        real(kind = dbl) :: cum
        real(kind = dbl), allocatable :: weight(:)

        if (allocated(ChunkStart)) deallocate(ChunkStart, ChunkEnd)

        !> No more pieces than periods a few to each, as for the slices
        !> before: a piece still has to pay for starting a process.
        want = min(ChunksPerWorker * nEff, MaxChunks, (iEnd - iStart) / 4)
        want = max(want, min(nEff, iEnd - iStart))

        !> Which periods a raw file covers. Files are in time order, so one
        !> forward sweep: skip the files that end before the period starts,
        !> then the period has data if the next file starts before it ends.
        allocate(weight(iStart:iEnd - 1))
        firstData = 0
        j = 1
        do p = iStart, iEnd - 1
            do while (j <= nFiles)
                if (Files(j)%timestamp + DatafileDateStep > Series(p)) exit
                j = j + 1
            end do
            weight(p) = EmptyPeriodWeight
            if (j <= nFiles) then
                if (Files(j)%timestamp < Series(p + 1)) then
                    weight(p) = 1d0
                    if (firstData == 0) firstData = p
                end if
            end if
        end do
        total = sum(weight)

        !> Cut where the running total passes each equal share. A cut is the
        !> first period of the next piece, so it lies strictly inside the
        !> range, and there is at most one per period, so each lies strictly
        !> after the one before; a share too thin to hold a period of its own
        !> is absorbed into the next piece.
        nCuts = 0
        cum = 0d0
        do p = iStart + 1, iEnd - 1
            cum = cum + weight(p - 1)
            if (nCuts + 1 >= want) exit
            if (cum >= total * dble(nCuts + 1) / dble(want)) then
                nCuts = nCuts + 1
                cuts(nCuts) = p
            end if
        end do
        deallocate(weight)

        !> Piece 1 reaches the first period with a file.
        if (firstData > 0) then
            do while (nCuts > 0)
                if (cuts(1) > firstData) exit
                cuts(1:nCuts - 1) = cuts(2:nCuts)
                nCuts = nCuts - 1
            end do
        end if

        NumChunks = nCuts + 1
        allocate(ChunkStart(NumChunks), ChunkEnd(NumChunks))
        ChunkStart(1) = iStart
        do c = 1, nCuts
            ChunkEnd(c) = cuts(c)
            ChunkStart(c + 1) = cuts(c)
        end do
        ChunkEnd(NumChunks) = iEnd
    end subroutine PlanPrepassChunks

    !***************************************************************************
    !> \brief The index range of piece k of the current pre-pass, half-open.
    !***************************************************************************
    subroutine PrepassChunk(k, sliceStart, sliceEnd)
        integer, intent(in) :: k
        integer, intent(out) :: sliceStart
        integer, intent(out) :: sliceEnd

        sliceStart = ChunkStart(k)
        sliceEnd = ChunkEnd(k)
    end subroutine PrepassChunk

    !***************************************************************************
    !> \brief How many pieces the current pre-pass was cut into.
    !***************************************************************************
    integer function PrepassChunkCount()
        PrepassChunkCount = NumChunks
    end function PrepassChunkCount

    !***************************************************************************
    !> \brief Where piece k of pre-pass `kind` leaves its records.
    !***************************************************************************
    character(PathLen) function BatchDumpPath(kind, k)
        character(*), intent(in) :: kind
        integer, intent(in) :: k
        character(16) :: tag

        write(tag, '(a,a,i2.2)') trim(kind), '_b', k
        BatchDumpPath = trim(TmpDir) // 'batch_' // trim(tag) // '.bin'
    end function BatchDumpPath

    !***************************************************************************
    !> \brief Cut the range into pieces, start the first workers, and return at
    !>        once, leaving piece 1 to the caller.
    !>
    !> nEff is how many processes may run at once, parent included. While the
    !> parent works through piece 1, nEff - 1 workers take pieces 2, 3, ...;
    !> WaitPrepassBatches hands out the rest as they come free.
    !>
    !> The parent takes the first piece itself rather than waiting idle. That is
    !> not only one process fewer: the code after the period loop reads global
    !> state the loop itself established - SortWindBySector wants the north
    !> offset out of E2Col, which is filled when a raw file's metadata is read -
    !> and a parent that had skipped the loop would arrive there with that state
    !> unset and bin every period into the wrong wind sector. Nothing enumerates
    !> what the finalisation depends on, so the parent runs a piece and thereby
    !> has all of it, exactly as it always did.
    !***************************************************************************
    subroutine StartPrepassBatches(kind, iStart, iEnd, nEff, Series, nSeries, &
            Files, nFiles, cuts)
        character(*), intent(in) :: kind
        integer, intent(in) :: iStart
        integer, intent(in) :: iEnd
        integer, intent(in) :: nEff
        integer, intent(in) :: nSeries
        integer, intent(in) :: nFiles
        type(DateType), intent(in) :: Series(nSeries)
        type(FileListType), intent(in) :: Files(nFiles)
        !> Where to cut, if the caller has decided that itself: the first period
        !> of every piece but the first, ascending, strictly inside the range.
        integer, intent(in), optional :: cuts(:)
        integer :: k
        integer :: covered
        logical :: ex
        character(PathLen) :: childPath
        character(PathLen) :: rcPath
        character(PathLen) :: exePath
        character(PathLen) :: envPath
        character(64) :: LogString

        call get_command_argument(0, value = exePath)

        !> A trailing separator immediately before a closing quote is read as an
        !> escape by cmd.exe, which would swallow the quote and split the path at
        !> its first space. InitEnv puts the separator back.
        envPath = homedir
        do while (len_trim(envPath) > 1)
            if (envPath(len_trim(envPath):len_trim(envPath)) /= slash) exit
            envPath = envPath(1:len_trim(envPath) - 1)
        end do

        if (present(cuts)) then
            if (size(cuts) + 1 > MaxChunks) &
                error stop 'Too many pieces requested for a parallel pass.'
            if (allocated(ChunkStart)) deallocate(ChunkStart, ChunkEnd)
            NumChunks = size(cuts) + 1
            allocate(ChunkStart(NumChunks), ChunkEnd(NumChunks))
            ChunkStart(1) = iStart
            do k = 1, size(cuts)
                ChunkEnd(k) = cuts(k)
                ChunkStart(k + 1) = cuts(k)
            end do
            ChunkEnd(NumChunks) = iEnd
        else
            call PlanPrepassChunks(iStart, iEnd, nEff, Series, nSeries, Files, nFiles)
        end if

        !> The pieces have to tile [iStart, iEnd) exactly. A gap drops periods
        !> from the fit and an overlap counts them twice, and neither shows up as
        !> anything but a slightly different answer - so it is checked here
        !> rather than left to be discovered.
        covered = 0
        do k = 1, NumChunks
            if (ChunkEnd(k) <= ChunkStart(k)) covered = -huge(covered)
            if (k > 1) then
                if (ChunkStart(k) /= ChunkEnd(k - 1)) covered = -huge(covered)
            end if
            covered = covered + (ChunkEnd(k) - ChunkStart(k))
        end do
        if (ChunkStart(1) /= iStart .or. ChunkEnd(NumChunks) /= iEnd) &
            covered = -huge(covered)
        if (covered /= iEnd - iStart) &
            error stop 'Pre-pass slices do not tile the period range.'

        write(LogString, '(i6)') NumChunks
        call LogSay('  The range is cut into ' // trim(adjustl(LogString)) &
            // ' pieces of similar work; each worker takes the next as it')
        call LogSay('  comes free.')

        !> A stale return code from an earlier attempt would be read as a worker
        !> that had already finished.
        do k = 2, NumChunks
            rcPath = ReturnCodePath(kind, k)
            inquire(file = trim(rcPath), exist = ex)
            if (ex) call system(comm_del // '"' // trim(rcPath) // '"' &
                // comm_err_redirect)
        end do

        !> One script per piece, written now, started when its turn comes. The
        !> alternative - putting each command line inside a launcher's own
        !> start/cmd quoting - nests quotes three deep, and every raw data
        !> directory here has a space in its name.
        do k = 2, NumChunks
            call WriteChildScript(kind, k, NumChunks, ChunkStart(k), ChunkEnd(k), &
                exePath, envPath, childPath)
        end do

        PoolKind = kind
        PoolEff = nEff
        PoolDone = 0
        if (allocated(Running)) deallocate(Running)
        allocate(Running(NumChunks))
        Running = .false.
        NextChunk = 2
        call LaunchChunks(kind, min(nEff - 1, NumChunks - 1))
        Running(2:NextChunk - 1) = .true.
        PoolActive = .true.
    end subroutine StartPrepassBatches

    !***************************************************************************
    !> \brief Start the next n pieces, each in its own worker, and return.
    !>
    !> Through a launcher script, as the slices always were: it holds one
    !> start line per piece and exits at once, leaving the workers running.
    !> They are still this process's descendants for everything that matters -
    !> the interface's job object, and their own watch on this process.
    !***************************************************************************
    subroutine LaunchChunks(kind, n)
        character(*), intent(in) :: kind
        integer, intent(in) :: n
        integer :: k
        integer :: u
        integer :: io_status
        character(PathLen) :: masterPath
        character(PathLen) :: childPath
        character(16) :: tag
        character(2048) :: cmd

        if (n <= 0) return

        if (OS == 'win') then
            masterPath = trim(TmpDir) // 'batch_' // trim(kind) // '_run.bat'
        else
            masterPath = trim(TmpDir) // 'batch_' // trim(kind) // '_run.sh'
        end if

        open(newunit = u, file = trim(masterPath), status = 'replace', &
            iostat = io_status)
        if (io_status /= 0) then
            call LogSay(' Could not write the batch launcher to ' // trim(masterPath))
            error stop 'Parallel pre-pass could not be started.'
        end if
        if (OS == 'win') then
            write(u, '(a)') '@echo off'
        else
            write(u, '(a)') '#!/bin/sh'
        end if

        do k = NextChunk, NextChunk + n - 1
            write(tag, '(a,a,i2.2)') trim(kind), '_b', k
            if (OS == 'win') then
                childPath = trim(TmpDir) // 'batch_' // trim(tag) // '.bat'
                write(u, '(a)') 'start "" /B cmd /c "' // trim(childPath) // '"'
            else
                childPath = trim(TmpDir) // 'batch_' // trim(tag) // '.sh'
                write(u, '(a)') 'sh "' // trim(childPath) // '" &'
            end if
        end do
        close(u)
        NextChunk = NextChunk + n

        !> On Windows this returns as soon as the workers are started. On the
        !> others each line is backgrounded, so the launcher returns as soon
        !> as it has started them too.
        if (OS == 'win') then
            cmd = 'cmd /c "' // trim(masterPath) // '"'
        else
            cmd = 'sh "' // trim(masterPath) // '"'
        end if
        call system(trim(cmd))
    end subroutine LaunchChunks

    !***************************************************************************
    !> \brief Parent, between two of its own periods: hand a free core the next
    !>        piece, and fail loudly as soon as any worker does.
    !>
    !> A worker that finishes while the parent is still on piece 1 would
    !> otherwise leave its core idle until the parent gets to the wait - with
    !> PWB for a while, since a detection piece is cheaper than the parent's
    !> own, which computes fluxes too. Called once per period; it only looks
    !> for return-code files, so it costs nothing a period would notice. Quiet:
    !> the period loop may have a progress line open, and the pieces finished
    !> are counted when the parent starts waiting.
    !***************************************************************************
    subroutine TopUpPrepassBatches()
        logical :: progressed

        if (.not. PoolActive) return
        call PollPrepassPool(PoolEff - 1, .true., progressed)
    end subroutine TopUpPrepassBatches

    !***************************************************************************
    !> \brief Hand out the remaining pieces as workers come free, and fail
    !>        loudly as soon as any of them does.
    !>
    !> The parent has finished piece 1 when it gets here, so from now on it
    !> only dispatches, and nEff workers run at once. Each poll looks for the
    !> return codes of the running pieces; every one that has appeared frees a
    !> place for the next piece. A fast core therefore simply finishes more
    !> pieces than a slow one, and nobody idles while work is left.
    !>
    !> A worker that dies has almost always hit a data or configuration fault
    !> the serial run would have hit too. Carrying on with the pieces that did
    !> work would fit the planar fit, or the time-lag windows, to less data than
    !> was asked for and say so only in a line of log - so this stops instead,
    !> at once rather than after the rest have finished. The workers still
    !> running notice within one period that this process has gone, and stop.
    !***************************************************************************
    subroutine WaitPrepassBatches(kind, nEff)
        character(*), intent(in) :: kind
        integer, intent(in) :: nEff
        integer :: k
        integer :: ticks
        logical :: progressed

        if (.not. PoolActive .or. trim(kind) /= trim(PoolKind)) &
            error stop 'Waiting for workers of a pass that was not started.'
        PoolEff = nEff

        call LogSay('  Waiting for the workers:')
        !> Piece 1 is the parent's own, done now; and whatever finished while
        !> the parent was on it.
        PoolDone = PoolDone + 1
        if (PoolDone > 1) call SayPiecesDone()
        ticks = 0
        do
            if (PoolDone >= NumChunks) exit
            call PollPrepassPool(nEff, .false., progressed)
            if (progressed) ticks = 0
            if (PoolDone >= NumChunks) exit

            call system(comm_sleep)
            ticks = ticks + 1
            if (ticks > MaxWaitTicks) then
                call LogSay('')
                call LogSay(' No pre-pass worker has finished within a day.')
                error stop 'Parallel pre-pass timed out.'
            end if
        end do
        deallocate(Running)
        PoolActive = .false.

        !> Every worker's log, in piece order, so the run log reads as one walk
        !> through the range whatever order the pieces finished in.
        do k = 2, NumChunks
            call AppendWorkerLog(kind, k)
        end do
    end subroutine WaitPrepassBatches

    !***************************************************************************
    !> \brief One look at the pool: count the pieces whose workers have
    !>        returned, stop the run if one failed, and start pieces until
    !>        `slots` workers are running or none is left.
    !***************************************************************************
    subroutine PollPrepassPool(slots, quiet, progressed)
        integer, intent(in) :: slots
        logical, intent(in) :: quiet
        logical, intent(out) :: progressed
        integer :: k
        integer :: rc
        integer :: nRunning
        logical :: finished
        logical :: delivered
        character(64) :: LogString

        progressed = .false.
        do k = 2, NextChunk - 1
            if (.not. Running(k)) cycle
            call ChunkReturnCode(PoolKind, k, finished, rc)
            if (.not. finished) cycle
            Running(k) = .false.
            PoolDone = PoolDone + 1
            progressed = .true.

            delivered = rc == 0
            if (delivered) delivered = DumpExists(PoolKind, k)
            if (.not. delivered) then
                !> Ends the progress line the period loop may have left open.
                if (quiet) call LogSay('')
                call AppendWorkerLog(PoolKind, k)
                write(LogString, '(i6)') k
                if (rc /= 0) then
                    !> Its console output rather than its log: a worker that died
                    !> rather than returned never closed the log, so the last thing
                    !> it managed to say - which is the thing worth reading - is
                    !> only in what the launcher captured.
                    call LogSay(' Pre-pass worker ' // trim(adjustl(LogString)) &
                        // ' failed. What it printed before it stopped:')
                    call DumpWorkerStdout(PoolKind, k)
                    error stop 'A parallel pre-pass worker failed.'
                end if
                call LogSay(' Pre-pass worker ' // trim(adjustl(LogString)) &
                    // ' exited cleanly but wrote no records. What it printed:')
                call DumpWorkerStdout(PoolKind, k)
                error stop 'A parallel pre-pass worker produced no output.'
            end if

            if (.not. quiet) call SayPiecesDone()
        end do

        !> Keep `slots` workers busy while pieces are left.
        nRunning = count(Running)
        if (nRunning < slots .and. NextChunk <= NumChunks) then
            k = NextChunk
            call LaunchChunks(PoolKind, min(slots - nRunning, NumChunks - NextChunk + 1))
            Running(k:NextChunk - 1) = .true.
        end if
    end subroutine PollPrepassPool

    subroutine SayPiecesDone()
        character(64) :: LogString
        character(64) :: CountString

        write(LogString, '(i6)') PoolDone
        write(CountString, '(i6)') NumChunks
        call LogSay('   ' // trim(adjustl(LogString)) // ' of ' &
            // trim(adjustl(CountString)) // ' pieces done.')
    end subroutine SayPiecesDone

    !***************************************************************************
    !> rief Parent: block until piece k is done, keeping the other workers
    !>        busy meanwhile.
    !>
    !> For a parent that consumes the pieces in order while they are being
    !> produced, rather than waiting for all of them: it is busy with that
    !> itself, so nEff - 1 workers run. Quiet, like TopUpPrepassBatches.
    !***************************************************************************
    subroutine WaitForBatchPiece(k)
        integer, intent(in) :: k
        integer :: ticks
        logical :: progressed

        if (.not. PoolActive) return
        if (k < 2 .or. k > NumChunks) return
        ticks = 0
        do
            call PollPrepassPool(PoolEff - 1, .true., progressed)
            if (k < NextChunk .and. .not. Running(k)) exit
            if (progressed) ticks = 0
            call system(comm_sleep)
            ticks = ticks + 1
            if (ticks > MaxWaitTicks) then
                call LogSay('')
                call LogSay(' No worker has finished within a day.')
                error stop 'Parallel pass timed out.'
            end if
        end do
    end subroutine WaitForBatchPiece

    !***************************************************************************
    !> rief Parent: start no more pieces, wait for the running ones to end,
    !>        and close the pool, without folding their logs into the run's.
    !>
    !> For a parent that has taken what it needed from the pieces - all of
    !> them, or fewer, when its input ended early - and says itself whatever
    !> they said.
    !***************************************************************************
    subroutine CancelBatchPool()
        integer :: ticks
        logical :: progressed

        if (.not. PoolActive) return
        NextChunk = NumChunks + 1
        ticks = 0
        do
            if (count(Running) == 0) exit
            call PollPrepassPool(0, .true., progressed)
            if (count(Running) == 0) exit
            if (progressed) ticks = 0
            call system(comm_sleep)
            ticks = ticks + 1
            if (ticks > MaxWaitTicks) then
                call LogSay('')
                call LogSay(' No worker has finished within a day.')
                error stop 'Parallel pass timed out.'
            end if
        end do
        deallocate(Running)
        PoolActive = .false.
    end subroutine CancelBatchPool

    !***************************************************************************
    !> \brief Whether piece k's worker has finished, and with what code.
    !>
    !> The script creates the return-code file and then echoes the code into
    !> it, so a read that lands between the two finds an empty file. That reads
    !> as "not yet", and the next poll finds the code.
    !***************************************************************************
    subroutine ChunkReturnCode(kind, k, finished, rc)
        character(*), intent(in) :: kind
        integer, intent(in) :: k
        logical, intent(out) :: finished
        integer, intent(out) :: rc
        integer :: u
        integer :: io_status
        logical :: ex

        finished = .false.
        rc = -1
        inquire(file = trim(ReturnCodePath(kind, k)), exist = ex)
        if (.not. ex) return
        open(newunit = u, file = trim(ReturnCodePath(kind, k)), &
            status = 'old', iostat = io_status)
        if (io_status /= 0) return
        read(u, *, iostat = io_status) rc
        close(u)
        if (io_status /= 0) then
            rc = -1
            return
        end if
        finished = .true.
    end subroutine ChunkReturnCode

    logical function DumpExists(kind, k)
        character(*), intent(in) :: kind
        integer, intent(in) :: k

        inquire(file = trim(BatchDumpPath(kind, k)), exist = DumpExists)
    end function DumpExists

    !***************************************************************************
    !> \brief Write the command line for one worker into its own script.
    !***************************************************************************
    subroutine WriteChildScript(kind, k, nChunks, sliceStart, sliceEnd, &
            exePath, envPath, childPath)
        character(*), intent(in) :: kind
        integer, intent(in) :: k
        integer, intent(in) :: nChunks
        integer, intent(in) :: sliceStart
        integer, intent(in) :: sliceEnd
        character(*), intent(in) :: exePath
        character(*), intent(in) :: envPath
        character(PathLen), intent(out) :: childPath
        integer :: u
        integer :: io_status
        character(64) :: batchArg
        character(2048) :: cmd
        character(16) :: tag
        character(16) :: parentId

        write(tag, '(a,a,i2.2)') trim(kind), '_b', k
        if (OS == 'win') then
            childPath = trim(TmpDir) // 'batch_' // trim(tag) // '.bat'
        else
            childPath = trim(TmpDir) // 'batch_' // trim(tag) // '.sh'
        end if

        write(batchArg, '(a,a,i0,a,i0,a,i0,a,i0)') trim(kind), ':', k, ':', &
            nChunks, ':', sliceStart, ':', sliceEnd
        write(parentId, '(i0)') ProcessSelfId()

        !> PrjPath rather than the path this program was handed: an EddyPro
        !> project has already been imported into one of ours by now, and N
        !> workers all importing it again would race to write the same file.
        !>
        !> -c console rather than the caller's own mode, so a worker of a run
        !> started from the interface does not also write progress the
        !> interface would try to read as the parent's.
        !> Every switch comes BEFORE the project path. The argument loop reads
        !> the project path and the switches from the same list, and a switch
        !> that lands after the path is the one most likely to be mishandled;
        !> putting them first makes the worker's instructions independent of
        !> that. -j 1 is belt and braces on top of the --batch interlock.
        cmd = '"' // trim(exePath) // '"' &
            // ' -s ' // trim(OS) &
            // ' -e "' // trim(envPath) // '"' &
            // ' -m ' // trim(EddyFlowProj%run_env) &
            // ' -c console' &
            // ' -j 1' &
            // ' --batch ' // trim(batchArg) &
            // ' --batch-out "' // trim(BatchDumpPath(kind, k)) // '"' &
            // ' --batch-tmp "' // trim(NoTrailingSlash(TmpDir)) // '"' &
            // ' --batch-parent ' // trim(parentId) &
            // ' "' // trim(PrjPath) // '"'

        open(newunit = u, file = trim(childPath), status = 'replace', &
            iostat = io_status)
        if (io_status /= 0) then
            call LogSay(' Could not write the worker script ' // trim(childPath))
            error stop 'Parallel pre-pass could not be started.'
        end if
        if (OS == 'win') then
            write(u, '(a)') '@echo off'
            write(u, '(a)') trim(cmd) // ' > "' // trim(StdoutPath(kind, k)) &
                // '" 2>&1'
            write(u, '(a)') 'echo %errorlevel% > "' &
                // trim(ReturnCodePath(kind, k)) // '"'
        else
            write(u, '(a)') '#!/bin/sh'
            write(u, '(a)') trim(cmd) // ' > "' // trim(StdoutPath(kind, k)) &
                // '" 2>&1'
            write(u, '(a)') 'echo $? > "' // trim(ReturnCodePath(kind, k)) // '"'
        end if
        close(u)
    end subroutine WriteChildScript

    character(PathLen) function ReturnCodePath(kind, k)
        character(*), intent(in) :: kind
        integer, intent(in) :: k
        character(16) :: tag

        write(tag, '(a,a,i2.2)') trim(kind), '_b', k
        ReturnCodePath = trim(TmpDir) // 'batch_' // trim(tag) // '.rc'
    end function ReturnCodePath

    character(PathLen) function StdoutPath(kind, k)
        character(*), intent(in) :: kind
        integer, intent(in) :: k
        character(16) :: tag

        write(tag, '(a,a,i2.2)') trim(kind), '_b', k
        StdoutPath = trim(TmpDir) // 'batch_' // trim(tag) // '.out'
    end function StdoutPath

    !***************************************************************************
    !> \brief Fold a worker's run log into the parent's, so the run has one.
    !***************************************************************************
    subroutine AppendWorkerLog(kind, k)
        character(*), intent(in) :: kind
        integer, intent(in) :: k
        integer :: u
        integer :: io_status
        logical :: ex
        character(PathLen) :: logPath
        character(1024) :: dataline
        character(64) :: LogString

        logPath = trim(BatchDumpPath(kind, k)) // '.log'
        inquire(file = trim(logPath), exist = ex)
        if (.not. ex) return

        write(LogString, '(i6)') k
        call LogSay('')
        call LogSay(' ---- worker ' // trim(adjustl(LogString)) // ' ----')
        open(newunit = u, file = trim(logPath), status = 'old', iostat = io_status)
        if (io_status /= 0) return
        do
            read(u, '(a)', iostat = io_status) dataline
            if (io_status /= 0) exit
            write(ulog, '(a)') trim(dataline)
        end do
        close(u)
        call LogSay(' ---- end of worker ' // trim(adjustl(LogString)) // ' ----')
        call LogSay('')
    end subroutine AppendWorkerLog

    !***************************************************************************
    !> \brief Echo what a failed worker printed, which is where a crash lands.
    !***************************************************************************
    subroutine DumpWorkerStdout(kind, k)
        character(*), intent(in) :: kind
        integer, intent(in) :: k
        integer :: u
        integer :: io_status
        logical :: ex
        character(1024) :: dataline

        inquire(file = trim(StdoutPath(kind, k)), exist = ex)
        if (.not. ex) return
        open(newunit = u, file = trim(StdoutPath(kind, k)), status = 'old', &
            iostat = io_status)
        if (io_status /= 0) return
        do
            read(u, '(a)', iostat = io_status) dataline
            if (io_status /= 0) exit
            call LogSay('   ' // trim(dataline))
        end do
        close(u)
    end subroutine DumpWorkerStdout

    !***************************************************************************
    !> \brief How many cores this machine says it has.
    !***************************************************************************
    integer function DetectCoreCount()
        integer :: u
        integer :: io_status
        integer :: n
        character(32) :: nprocs

        n = 0
        if (OS == 'win') then
            call get_environment_variable('NUMBER_OF_PROCESSORS', nprocs, &
                status = io_status)
            if (io_status == 0) read(nprocs, *, iostat = io_status) n
            if (io_status /= 0) n = 0
        else
            call system('getconf _NPROCESSORS_ONLN > "' // trim(TmpDir) &
                // 'ncpu.tmp"' // comm_err_redirect)
            open(newunit = u, file = trim(TmpDir) // 'ncpu.tmp', &
                status = 'old', iostat = io_status)
            if (io_status == 0) then
                read(u, *, iostat = io_status) n
                close(u)
                if (io_status /= 0) n = 0
            end if
        end if
        if (n < 1) n = 1
        if (n > MaxWorkers) n = MaxWorkers
        DetectCoreCount = n
    end function DetectCoreCount

    !***************************************************************************
    !> \brief A directory without its trailing separator.
    !>
    !> For a quoted command-line argument on Windows: a backslash right before
    !> the closing quote escapes it, and the argument would run on into the
    !> next one.
    !***************************************************************************
    character(PathLen) function NoTrailingSlash(dir)
        character(*), intent(in) :: dir
        integer :: n

        NoTrailingSlash = adjustl(dir)
        n = len_trim(NoTrailingSlash)
        if (n > 1) then
            if (NoTrailingSlash(n:n) == slash) NoTrailingSlash(n:n) = ' '
        end if
    end function NoTrailingSlash

    !***************************************************************************
    !> \brief Stop this worker if the process that started it has gone.
    !>
    !> Called at the top of every period of a worker's slice. A worker is its
    !> own process, started through a shell script, and nothing ties its life
    !> to the parent's: killing the parent used to leave every worker running
    !> to the end of its slice - hours, on the Yatir run - writing records into
    !> a directory nobody would read again. Six were found still running two
    !> minutes after the interface's Stop.
    !>
    !> Nothing is written: a partial slice is not a record of anything. The
    !> worker tidies its temporary directory as it does on finishing, and exits
    !> with 3, distinct from both success and an ordinary failure, for whoever
    !> reads the .rc file later.
    !>
    !> A period is the granularity because it is the only safe point - between
    !> periods the worker holds nothing half-done - and it bounds the delay to
    !> one period's work, ~25 s at the worst measured so far.
    !***************************************************************************
    subroutine StopIfParentGone()
        integer :: u
        integer :: rmdir_status
        logical :: op
        character(PathLen) :: dir

        if (BatchIndex <= 0) return
        if (.not. ParentGone()) return
        !> Ends the progress line the period loop left open, so the reason
        !> does not read as part of a date.
        call LogSay('')
        call LogSay(' The process that started this worker has ended,')
        call LogSay(' so nothing will read its records. Stopping without them.')
        !> A production worker's output folder: its own files are closed first,
        !> since Windows will not delete a file that is open.
        if (len_trim(BatchOwnOutDir) > 0) then
            do u = uqc, uflxnt
                inquire(unit = u, opened = op)
                if (op) close(u)
            end do
            dir = NoTrailingSlash(BatchOwnOutDir)
            rmdir_status = system(trim(comm_rmdir) // ' "' // trim(dir) // '"' &
                // comm_err_redirect)
        end if
        call FinishBatchWorker()
        stop 3
    end subroutine StopIfParentGone

    !***************************************************************************
    !> \brief Stop unless this worker was launched for pre-pass `kind`.
    !>
    !> Every worker runs the program from the top, so each pre-pass it meets
    !> before its own has to stand aside for it - and the one time that did not
    !> happen, a planar-fit worker reported time-lag records as its result.
    !> Writing a dump for the wrong pre-pass is the moment that mistake becomes
    !> data, so this is where it is refused.
    !***************************************************************************
    subroutine RequireBatchKind(kind)
        character(*), intent(in) :: kind

        if (BatchKind == kind) return
        call LogSay(' This pre-pass worker was launched for the ' // trim(BatchKind) &
            // ' pre-pass but reached the ' // trim(kind) // ' one.')
        error stop 'A pre-pass worker ran the wrong pre-pass.'
    end subroutine RequireBatchKind

    !***************************************************************************
    !> \brief Tidy up after a worker's slice, before it stops.
    !>
    !> A worker stops as soon as its records are written, so it never reaches
    !> the cleanup at the end of a run and its temporary directory - extracted
    !> archives, file lists, shared-link scratch - was left behind, one per
    !> worker per pre-pass. Everything the parent needs from it, the records
    !> and the worker's log, is beside --batch-out in the parent's directory,
    !> not in this one.
    !>
    !> Only in desktop mode, like the parent's own: in embedded mode every
    !> process shares one temporary directory, and it is the parent's.
    !***************************************************************************
    subroutine FinishBatchWorker()
        integer :: rmdir_status

        call GhgPrefetchCleanup()
        call RemoteCleanup()
        if (EddyFlowProj%run_env == 'desktop') &
            rmdir_status = system(trim(comm_rmdir) // ' "' &
                // trim(adjustl(TmpDir)) // '"' // comm_err_redirect)
    end subroutine FinishBatchWorker

    !***************************************************************************
    !> \brief Test hook: cut a range into pieces of this many periods (or
    !>        files), whatever the core count, so that every fixture - most
    !>        span a day or less - is cut at nearly every half-hour. The value
    !>        of environment variable `name`; 0 if unset.
    !***************************************************************************
    integer function ForcedPieceLength(name)
        character(*), intent(in) :: name
        character(32) :: val
        integer :: stat
        integer :: io_status

        ForcedPieceLength = 0
        call get_environment_variable(name, val, status = stat)
        if (stat /= 0 .or. len_trim(val) == 0) return
        read(val, *, iostat = io_status) ForcedPieceLength
        if (io_status /= 0) ForcedPieceLength = 0
        ForcedPieceLength = max(0, ForcedPieceLength)
    end function ForcedPieceLength

    !***************************************************************************
    !> \brief Parent: remove the worker folders runs that were killed left.
    !>
    !> A run stopped from the interface is terminated whole, workers and all,
    !> and nothing gets the chance to tidy up - so each such run left one
    !> folder per worker under tmp, holding the rows its piece had written.
    !> A folder is a dead run's if the process id in its name is not running,
    !> or is this process's own: then the id has been reused, and the folder
    !> is about to be wanted again.
    !***************************************************************************
    subroutine RemoveStaleWorkerRoots()
        integer :: u
        integer :: io_status
        integer :: under
        integer :: pid
        integer :: rmdir_status
        character(PathLen) :: base
        character(PathLen) :: listing
        character(PathLen) :: entry
        character(2048) :: cmd

        base = trim(homedir) // 'tmp' // slash
        listing = trim(TmpDir) // 'pr_stale_roots.txt'
        if (OS == 'win') then
            cmd = 'dir /b /ad "' // trim(base) // 'p*_*" > "' // trim(listing) &
                // '"' // comm_err_redirect
        else
            cmd = 'cd "' // trim(base) // '" && ls -1d p*_* > "' // trim(listing) &
                // '"' // comm_err_redirect
        end if
        rmdir_status = system(trim(cmd))
        open(newunit = u, file = trim(listing), status = 'old', action = 'read', &
            iostat = io_status)
        if (io_status /= 0) return
        do
            read(u, '(a)', iostat = io_status) entry
            if (io_status /= 0) exit
            entry = adjustl(entry)
            !> p<pid>_<kind><k>: anything else under tmp is not ours.
            under = index(entry, '_')
            if (entry(1:1) /= 'p' .or. under < 3) cycle
            if (verify(entry(2:under - 1), '0123456789') /= 0) cycle
            if (index(' pd pr fx sb ', ' ' // entry(under + 1:under + 2) // ' ') == 0) cycle
            if (verify(trim(entry(under + 3:)), '0123456789') /= 0) cycle
            read(entry(2:under - 1), *, iostat = io_status) pid
            if (io_status /= 0) cycle
            if (pid /= ProcessSelfId()) then
                if (ProcessAlive(pid)) cycle
            end if
            rmdir_status = system(trim(comm_rmdir) // ' "' // trim(base) &
                // trim(entry) // '"' // comm_err_redirect)
        end do
        close(u, status = 'delete')
    end subroutine RemoveStaleWorkerRoots

    !***************************************************************************
    !> \brief Worker k's own output folder: tmp\p<parent's process id>_<kind><k>.
    !>
    !> Not inside the parent's temporary folder, which would be shorter to
    !> clean up: that one's name is 24 characters, and under it a worker's
    !> per-period file names - a timestamp, the run stamp, the run mode - went
    !> past Windows' 260-character path limit in a home folder where the
    !> parent's own output did not. The parent's process id keeps two runs
    !> sharing a home folder apart, and the kind the detection and production
    !> workers of one piece.
    !***************************************************************************
    character(PathLen) function WorkerRoot(parentPid, kind, k)
        integer, intent(in) :: parentPid
        character(*), intent(in) :: kind
        integer, intent(in) :: k
        character(32) :: tag

        write(tag, '(a,i0,a,a,i2.2)') 'p', parentPid, '_', trim(kind), k
        WorkerRoot = trim(homedir) // 'tmp' // slash // trim(tag) // slash
    end function WorkerRoot

    !> Parent: remove every worker folder of one kind.
    subroutine RemoveWorkerRoots(kind)
        character(*), intent(in) :: kind
        integer :: k
        integer :: rmdir_status
        character(PathLen) :: path

        do k = 2, PrepassChunkCount()
            path = WorkerRoot(ProcessSelfId(), kind, k)
            rmdir_status = system(trim(comm_rmdir) // ' "' &
                // trim(NoSlashEnd(path)) // '"' // comm_err_redirect)
        end do
    end subroutine RemoveWorkerRoots

    character(PathLen) function NoSlashEnd(dir)
        character(*), intent(in) :: dir
        integer :: n

        NoSlashEnd = dir
        n = len_trim(NoSlashEnd)
        if (n > 1) then
            if (NoSlashEnd(n:n) == slash) NoSlashEnd(n:n) = ' '
        end if
    end function NoSlashEnd

    !***************************************************************************
    !> \brief Copy bytes [offset, end) of `from` onto the end of `to`.
    !>
    !> Raw bytes, so nothing about a row - trailing blanks, line ends, UTF-8 -
    !> can change on the way.
    !***************************************************************************
    subroutine AppendBytes(from, offset, to)
        character(*), intent(in) :: from
        integer(8), intent(in) :: offset
        character(*), intent(in) :: to
        integer, parameter :: Chunk = 1048576
        integer :: ui
        integer :: uo
        integer :: io_status
        integer :: n
        integer(8) :: total
        integer(8) :: pos
        character(len = :), allocatable :: buf

        total = FileBytes(from)
        if (total <= offset) return
        open(newunit = ui, file = trim(from), access = 'stream', &
            form = 'unformatted', status = 'old', action = 'read', iostat = io_status)
        if (io_status /= 0) error stop 'Could not read a production worker file.'
        open(newunit = uo, file = trim(to), access = 'stream', &
            form = 'unformatted', status = 'old', position = 'append', &
            action = 'write', iostat = io_status)
        if (io_status /= 0) error stop 'Could not append to an output file.'
        allocate(character(len = Chunk) :: buf)
        pos = offset + 1
        do while (pos <= total)
            n = int(min(int(Chunk, 8), total - pos + 1))
            read(ui, pos = pos, iostat = io_status) buf(1:n)
            if (io_status /= 0) error stop 'Could not read a production worker file.'
            write(uo, iostat = io_status) buf(1:n)
            if (io_status /= 0) error stop 'Could not append to an output file.'
            pos = pos + n
        end do
        close(ui)
        close(uo)
    end subroutine AppendBytes

    integer(8) function FileBytes(path)
        character(*), intent(in) :: path
        logical :: ex

        FileBytes = 0
        inquire(file = trim(path), exist = ex, size = FileBytes)
        if (.not. ex .or. FileBytes < 0) FileBytes = 0
    end function FileBytes

end module m_batch_pool
