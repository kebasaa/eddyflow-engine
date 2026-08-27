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
!              This module splits that loop into contiguous slices and runs
!              each in a copy of this program, which processes its slice and
!              writes the records it produced to a file instead of going on to
!              compute fluxes. The parent reads the slices back in order and
!              hands the concatenation to the unchanged finalisation code.
!
!              Processes rather than threads because the period loop reaches
!              most of the program's global state - Stats, E2Col, Essentials,
!              the metadata - and making that thread-safe would be a rewrite.
!              A process has its own copy of all of it by construction.
!
!              The concatenation is exact: the records are fixed-size derived
!              types with no allocatable components, written unformatted and
!              read back by the same binary, and appended in slice order.
!
!              The PWB cache pre-pass is deliberately NOT split. Its classifier
!              decides a period partly from the last settled detection before
!              it, and that chain has no time limit - the streaming classifier
!              sets its flag once and never clears it. A lead-in could rebuild
!              the state only if it contained a settled detection for every
!              gas, and on real data the weak species do not oblige: over two
!              days of CH-LAE, COS reached no settled detection at all and
!              nitrous oxide reached two in 103 periods. So a split there could
!              not reproduce a single pass, and it is not offered.
!
! \author      Jonathan Muller
! \note
! \sa          eddyflow-rp_main.f90, pwb_timelag_handle.f90
!***************************************************************************
module m_prepass_parallel
    use m_rp_global_var
    implicit none
    private

    public :: PlanPrepassBatches, PrepassSlice
    public :: StartPrepassBatches, WaitPrepassBatches
    public :: BatchDumpPath
    public :: WriteTlagBatchDump, MergeTlagBatchDumps
    public :: WritePwbBatchDump, MergePwbBatchDumps
    public :: WritePfBatchDump, MergePfBatchDumps

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
    character(20), parameter :: BatchMagic = 'EDDYFLOW_PREPASS_04 '

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
    subroutine PlanPrepassBatches(nPeriods, allowed, nEff)
        integer, intent(in) :: nPeriods
        logical, intent(in) :: allowed
        integer, intent(out) :: nEff
        integer :: requested
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
        nEff = max(1, min(requested, nPeriods / 4))
        if (nEff <= 1) then
            nEff = 1
            return
        end if

        write(LogString, '(i6)') nEff
        call LogSay('  Splitting the pre-pass across ' &
            // trim(adjustl(LogString)) // ' worker processes.')
    end subroutine PlanPrepassBatches

    !***************************************************************************
    !> \brief The index range worker k owns.
    !>
    !> The range is HALF-OPEN, [iStart, iEnd), because that is what the period
    !> loops themselves do: both exit on `pcount >= endIndex`, so the end index
    !> is one past the last period processed. sliceEnd is exclusive for the
    !> same reason and can be assigned straight to the loop's end index.
    !>
    !> Slices tile the range exactly - no gaps, no overlaps - with the
    !> remainder spread over the first slices rather than dumped on the last.
    !***************************************************************************
    subroutine PrepassSlice(iStart, iEnd, nEff, k, sliceStart, sliceEnd)
        integer, intent(in) :: iStart
        integer, intent(in) :: iEnd
        integer, intent(in) :: nEff
        integer, intent(in) :: k
        integer, intent(out) :: sliceStart
        integer, intent(out) :: sliceEnd
        integer :: total
        integer :: base
        integer :: rem
        integer :: len

        total = iEnd - iStart
        base = total / nEff
        rem = total - base * nEff

        sliceStart = iStart + (k - 1) * base + min(k - 1, rem)
        len = base
        if (k <= rem) len = len + 1
        sliceEnd = sliceStart + len
    end subroutine PrepassSlice

    !***************************************************************************
    !> \brief Where worker k of pre-pass `kind` leaves its records.
    !***************************************************************************
    character(PathLen) function BatchDumpPath(kind, k)
        character(*), intent(in) :: kind
        integer, intent(in) :: k
        character(16) :: tag

        write(tag, '(a,a,i2.2)') trim(kind), '_b', k
        BatchDumpPath = trim(TmpDir) // 'batch_' // trim(tag) // '.bin'
    end function BatchDumpPath

    !***************************************************************************
    !> \brief Launch workers 2..nEff and return at once, leaving slice 1 to the
    !>        caller.
    !>
    !> The parent takes the first slice itself rather than waiting idle. That is
    !> not only one process fewer: the code after the period loop reads global
    !> state the loop itself established - SortWindBySector wants the north
    !> offset out of E2Col, which is filled when a raw file's metadata is read -
    !> and a parent that had skipped the loop would arrive there with that state
    !> unset and bin every period into the wrong wind sector. Nothing enumerates
    !> what the finalisation depends on, so the parent runs a slice and thereby
    !> has all of it, exactly as it always did.
    !***************************************************************************
    subroutine StartPrepassBatches(kind, iStart, iEnd, nEff)
        character(*), intent(in) :: kind
        integer, intent(in) :: iStart
        integer, intent(in) :: iEnd
        integer, intent(in) :: nEff
        integer :: k
        integer :: u
        integer :: io_status
        integer :: sliceStart
        integer :: sliceEnd
        integer :: covered
        logical :: ex
        character(PathLen) :: masterPath
        character(PathLen) :: childPath
        character(PathLen) :: rcPath
        character(PathLen) :: exePath
        character(PathLen) :: envPath
        character(2048) :: cmd

        call get_command_argument(0, value = exePath)

        !> A trailing separator immediately before a closing quote is read as an
        !> escape by cmd.exe, which would swallow the quote and split the path at
        !> its first space. InitEnv puts the separator back.
        envPath = homedir
        do while (len_trim(envPath) > 1)
            if (envPath(len_trim(envPath):len_trim(envPath)) /= slash) exit
            envPath = envPath(1:len_trim(envPath) - 1)
        end do

        !> The slices have to tile [iStart, iEnd) exactly. A gap drops periods
        !> from the fit and an overlap counts them twice, and neither shows up as
        !> anything but a slightly different answer - so it is checked here
        !> rather than left to be discovered.
        covered = 0
        do k = 1, nEff
            call PrepassSlice(iStart, iEnd, nEff, k, sliceStart, sliceEnd)
            covered = covered + (sliceEnd - sliceStart)
        end do
        if (covered /= iEnd - iStart) &
            error stop 'Pre-pass slices do not tile the period range.'

        !> A stale return code from an earlier attempt would be read as a worker
        !> that had already finished.
        do k = 2, nEff
            rcPath = ReturnCodePath(kind, k)
            inquire(file = trim(rcPath), exist = ex)
            if (ex) call system(comm_del // '"' // trim(rcPath) // '"' &
                // comm_err_redirect)
        end do

        !> One script per worker, and a master that only launches them. The
        !> alternative - putting each command line inside the master's own
        !> start/cmd quoting - nests quotes three deep, and every raw data
        !> directory here has a space in its name.
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

        do k = 2, nEff
            call PrepassSlice(iStart, iEnd, nEff, k, sliceStart, sliceEnd)
            call WriteChildScript(kind, k, nEff, sliceStart, sliceEnd, &
                exePath, envPath, childPath)
            if (OS == 'win') then
                write(u, '(a)') 'start "" /B cmd /c "' // trim(childPath) // '"'
            else
                write(u, '(a)') 'sh "' // trim(childPath) // '" &'
            end if
        end do
        close(u)

        !> On Windows this returns as soon as the workers are started, so the
        !> caller gets on with its own slice and the waiting happens later. On
        !> the others the launcher ends in `wait`, so this call blocks - the
        !> parent's slice then runs after the workers rather than beside them,
        !> which costs one slice of wall clock and nothing in correctness.
        if (OS == 'win') then
            cmd = 'cmd /c "' // trim(masterPath) // '"'
        else
            cmd = 'sh "' // trim(masterPath) // '"'
        end if
        call system(trim(cmd))
    end subroutine StartPrepassBatches

    !***************************************************************************
    !> \brief Wait for workers 2..nEff, and fail loudly if any of them did.
    !>
    !> A worker that dies has almost always hit a data or configuration fault
    !> the serial run would have hit too. Carrying on with the slices that did
    !> work would fit the planar fit, or the time-lag windows, to less data than
    !> was asked for and say so only in a line of log - so this stops instead.
    !***************************************************************************
    subroutine WaitPrepassBatches(kind, nEff)
        character(*), intent(in) :: kind
        integer, intent(in) :: nEff
        integer :: k
        integer :: u
        integer :: rc
        integer :: ticks
        integer :: io_status
        logical :: ex
        logical :: allDone
        character(64) :: LogString

        call LogSayNoAdv('  Waiting for the workers..')
        ticks = 0
        do
            allDone = .true.
            do k = 2, nEff
                inquire(file = trim(ReturnCodePath(kind, k)), exist = ex)
                if (.not. ex) allDone = .false.
            end do
            if (allDone) exit
            call system(comm_sleep)
            ticks = ticks + 1
            if (ticks > MaxWaitTicks) then
                call LogSay('')
                call LogSay(' A pre-pass worker has not finished within a day.')
                error stop 'Parallel pre-pass timed out.'
            end if
        end do
        !> The return code is echoed into the file after it is created, so a read
        !> arriving between the two sees an empty file.
        call system(comm_sleep)
        call LogSay(' Done.')

        do k = 2, nEff
            rc = -1
            open(newunit = u, file = trim(ReturnCodePath(kind, k)), &
                status = 'old', iostat = io_status)
            if (io_status == 0) then
                read(u, *, iostat = io_status) rc
                close(u)
                if (io_status /= 0) rc = -1
            end if

            call AppendWorkerLog(kind, k)

            if (rc /= 0) then
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

            inquire(file = trim(BatchDumpPath(kind, k)), exist = ex)
            if (.not. ex) then
                write(LogString, '(i6)') k
                call LogSay(' Pre-pass worker ' // trim(adjustl(LogString)) &
                    // ' exited cleanly but wrote no records. What it printed:')
                call DumpWorkerStdout(kind, k)
                error stop 'A parallel pre-pass worker produced no output.'
            end if
        end do
    end subroutine WaitPrepassBatches

    !***************************************************************************
    !> \brief Write the command line for one worker into its own script.
    !***************************************************************************
    subroutine WriteChildScript(kind, k, nEff, sliceStart, sliceEnd, &
            exePath, envPath, childPath)
        character(*), intent(in) :: kind
        integer, intent(in) :: k
        integer, intent(in) :: nEff
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

        write(tag, '(a,a,i2.2)') trim(kind), '_b', k
        if (OS == 'win') then
            childPath = trim(TmpDir) // 'batch_' // trim(tag) // '.bat'
        else
            childPath = trim(TmpDir) // 'batch_' // trim(tag) // '.sh'
        end if

        write(batchArg, '(a,a,i0,a,i0,a,i0,a,i0)') trim(kind), ':', k, ':', &
            nEff, ':', sliceStart, ':', sliceEnd

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

        open(newunit = u, file = trim(BatchOutPath), form = 'unformatted', &
            access = 'stream', status = 'replace', iostat = io_status)
        if (io_status /= 0) &
            error stop 'A pre-pass worker could not write its records.'

        write(u) BatchMagic
        write(u) BatchIndex, BatchCount
        write(u) n
        if (n > 0) write(u) dataset(1:n)
        close(u)
    end subroutine WriteTlagBatchDump

    !***************************************************************************
    !> rief Read the workers' slices back, in order, as if one loop ran.
    !>
    !> Slice 1 is already in dataset(1:n) - the parent ran it itself - so this
    !> appends slices 2..nEff after it, which is the order the serial loop would
    !> have produced them in.
    !***************************************************************************
    subroutine MergeTlagBatchDumps(kind, nEff, dataset, nmax, n)
        character(*), intent(in) :: kind
        integer, intent(in) :: nEff
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
        type(TimeLagOptType), allocatable :: slice(:)

        do k = 2, nEff
            open(newunit = u, file = trim(BatchDumpPath(kind, k)), &
                form = 'unformatted', access = 'stream', status = 'old', &
                iostat = io_status)
            if (io_status /= 0) &
                error stop 'A pre-pass worker record file could not be read.'

            read(u) magic
            if (magic /= BatchMagic) &
                error stop 'A pre-pass worker record file is not one of ours.'
            read(u) idx, idxCount
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

        open(newunit = u, file = trim(BatchOutPath), form = 'unformatted', &
            access = 'stream', status = 'replace', iostat = io_status)
        if (io_status /= 0) &
            error stop 'A pre-pass worker could not write its PWB records.'

        write(u) BatchMagic
        write(u) BatchIndex, BatchCount
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
    !> \brief Read the workers' PWB slices back, in order, as if one loop ran.
    !>
    !> Slice 1 is already here - the parent ran it - so this appends 2..nEff
    !> after it. Order is the whole point: the post-pass sorts the table by
    !> timestamp with a stable insertion sort, so rows appended in slice order
    !> come out exactly as a single loop would have left them, and periods
    !> sharing a timestamp keep their gas order.
    !>
    !> The cache grows once per slice rather than once per row. Storing a row
    !> at a time reallocates and copies the whole table each time, which is
    !> quadratic and is what StorePwbTimelagCacheAt does for the serial walk.
    !***************************************************************************
    subroutine MergePwbBatchDumps(kind, nEff, dataset, nmax, nOpt)
        character(*), intent(in) :: kind
        integer, intent(in) :: nEff
        integer, intent(in) :: nmax
        type(TimeLagOptType), intent(inout) :: dataset(nmax)
        integer, intent(inout) :: nOpt
        integer :: k, i, u, io_status
        integer :: nrec, idx, idxCount
        character(20) :: magic
        type(PWBTimelagCacheEntryType), allocatable :: rows(:), grown(:)
        type(TimeLagOptType), allocatable :: slice(:)
        character(10), allocatable :: sdate(:)
        character(5), allocatable :: stime(:)

        do k = 2, nEff
            open(newunit = u, file = trim(BatchDumpPath(kind, k)), &
                form = 'unformatted', access = 'stream', status = 'old', &
                iostat = io_status)
            if (io_status /= 0) &
                error stop 'A pre-pass worker PWB file could not be read.'

            read(u) magic
            if (magic /= BatchMagic) &
                error stop 'A pre-pass worker PWB file is not one of ours.'
            read(u) idx, idxCount

            read(u) nrec
            if (nrec > 0) then
                allocate(rows(nrec))
                read(u) rows
                allocate(grown(PwbTimelagCacheN + nrec))
                if (PwbTimelagCacheN > 0) &
                    grown(1:PwbTimelagCacheN) = PwbTimelagCache(1:PwbTimelagCacheN)
                grown(PwbTimelagCacheN + 1:PwbTimelagCacheN + nrec) = rows
                call move_alloc(grown, PwbTimelagCache)
                PwbTimelagCacheN = PwbTimelagCacheN + nrec
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

        open(newunit = u, file = trim(BatchOutPath), form = 'unformatted', &
            access = 'stream', status = 'replace', iostat = io_status)
        if (io_status /= 0) &
            error stop 'A pre-pass worker could not write its records.'

        write(u) BatchMagic
        write(u) BatchIndex, BatchCount
        write(u) n
        if (n > 0) write(u) wind(1:n, 1:3)
        close(u)
    end subroutine WritePfBatchDump

    !***************************************************************************
    !> \brief Read the planar-fit slices back, after the parent's own.
    !***************************************************************************
    subroutine MergePfBatchDumps(nEff, wind, nmax, n)
        integer, intent(in) :: nEff
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
        real(kind = dbl), allocatable :: slice(:, :)

        do k = 2, nEff
            open(newunit = u, file = trim(BatchDumpPath('pf', k)), &
                form = 'unformatted', access = 'stream', status = 'old', &
                iostat = io_status)
            if (io_status /= 0) &
                error stop 'A pre-pass worker record file could not be read.'

            read(u) magic
            if (magic /= BatchMagic) &
                error stop 'A pre-pass worker record file is not one of ours.'
            read(u) idx, idxCount
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
