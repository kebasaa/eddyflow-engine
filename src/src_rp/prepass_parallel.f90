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
! \sa          eddyflow-rp_main.f90, pwb_timelag_handle.f90
!***************************************************************************
module m_prepass_parallel
    use m_rp_global_var
    use m_ghg_prefetch, only: GhgPrefetchCleanup
    use m_remote_source, only: RemoteCleanup
    use m_process_os, only: ProcessSelfId, ParentGone
    use m_pwb_timelag, only: AppendPwbCacheRows
    implicit none
    private

    public :: PlanPrepassBatches, PrepassChunk, PrepassChunkCount
    public :: StartPrepassBatches, WaitPrepassBatches
    public :: BatchDumpPath
    public :: WriteTlagBatchDump, MergeTlagBatchDumps
    public :: WritePwbBatchDump, MergePwbBatchDumps
    public :: WritePfBatchDump, MergePfBatchDumps
    public :: FinishBatchWorker, StopIfParentGone

    !> More workers than this is never a throughput win on a machine that also
    !> has to feed them raw data, and it multiplies the per-worker cost of
    !> listing the raw directory.
    integer, parameter :: MaxWorkers = 32

    !> Roughly one second per tick, so this is a full day of waiting. A
    !> pre-pass over a season legitimately takes hours; a worker that has hung
    !> should still not hold the parent for ever.
    integer, parameter :: MaxWaitTicks = 86400

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

        NextChunk = 2
        call LaunchChunks(kind, min(nEff - 1, NumChunks - 1))
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
        integer :: rc
        integer :: ticks
        integer :: nDone
        integer :: nRunning
        logical :: finished
        logical, allocatable :: running(:)
        character(64) :: LogString
        character(64) :: CountString

        allocate(running(NumChunks))
        running = .false.
        running(2:NextChunk - 1) = .true.

        write(CountString, '(i6)') NumChunks
        call LogSay('  Waiting for the workers:')
        nDone = 1
        ticks = 0
        do
            do k = 2, NextChunk - 1
                if (.not. running(k)) cycle
                call ChunkReturnCode(kind, k, finished, rc)
                if (.not. finished) cycle
                running(k) = .false.
                nDone = nDone + 1
                ticks = 0

                if (rc /= 0) then
                    call AppendWorkerLog(kind, k)
                    write(LogString, '(i6)') k
                    !> Its console output rather than its log: a worker that died
                    !> rather than returned never closed the log, so the last thing
                    !> it managed to say - which is the thing worth reading - is
                    !> only in what the launcher captured.
                    call LogSay(' Pre-pass worker ' // trim(adjustl(LogString)) &
                        // ' failed. What it printed before it stopped:')
                    call DumpWorkerStdout(kind, k)
                    error stop 'A parallel pre-pass worker failed.'
                end if

                if (.not. DumpExists(kind, k)) then
                    call AppendWorkerLog(kind, k)
                    write(LogString, '(i6)') k
                    call LogSay(' Pre-pass worker ' // trim(adjustl(LogString)) &
                        // ' exited cleanly but wrote no records. What it printed:')
                    call DumpWorkerStdout(kind, k)
                    error stop 'A parallel pre-pass worker produced no output.'
                end if

                write(LogString, '(i6)') nDone
                call LogSay('   ' // trim(adjustl(LogString)) // ' of ' &
                    // trim(adjustl(CountString)) // ' pieces done.')
            end do
            if (nDone >= NumChunks) exit

            !> Keep nEff workers busy while pieces are left.
            nRunning = count(running)
            if (nRunning < nEff .and. NextChunk <= NumChunks) then
                k = NextChunk
                call LaunchChunks(kind, min(nEff - nRunning, NumChunks - NextChunk + 1))
                running(k:NextChunk - 1) = .true.
            end if

            call system(comm_sleep)
            ticks = ticks + 1
            if (ticks > MaxWaitTicks) then
                call LogSay('')
                call LogSay(' No pre-pass worker has finished within a day.')
                error stop 'Parallel pre-pass timed out.'
            end if
        end do
        deallocate(running)

        !> Every worker's log, in piece order, so the run log reads as one walk
        !> through the range whatever order the pieces finished in.
        do k = 2, NumChunks
            call AppendWorkerLog(kind, k)
        end do
    end subroutine WaitPrepassBatches

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

end module m_prepass_parallel
