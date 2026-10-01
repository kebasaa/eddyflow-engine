!***************************************************************************
! survey_ghg_ac_freq.f90
! ----------------------
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
! \brief       Find where the acquisition frequency of a project's GHG files
!              changes, and the highest one, without reading them all.
! \author      Jonathan Muller
! \note        Two uses. The binned-spectra grid is built up to the highest
!              rate's Nyquist. And the list of changes, each with the time the
!              new rate starts, is shown before processing begins - nothing is
!              shown when there is one rate only, so such a log is unchanged.
!
!              Neither is what processing relies on: every file's rate is
!              checked anyway when it is imported, so a change this survey
!              misses is still processed correctly. What a miss costs is its
!              line in the list, and, for a rate above the grid, the spectra
!              above the grid's Nyquist, which Warning(117) then reports.
!
!              Reading a file's rate means extracting its .metadata, about
!              0.1 s. Reading every one would cost 45 minutes on a year of
!              half-hourly files. What is read instead:
!               - the first and the last file, and NumEven spread evenly
!                 between them;
!               - where the compressed size steps up or down from one file to
!                 the next, the larger of the two. Size scales with the number
!                 of rows, so 10 to 20 Hz roughly doubles it, and asking the
!                 file system for a size costs nothing. A short file at a data
!                 gap trips this too, which only spends one read.
!              at most MaxPeeks of those, the largest size steps first. Then,
!              between any two files read that disagree, the one halfway,
!              until the two that disagree are neighbours: the later one is
!              where the new rate starts. That is about log2(files between)
!              reads per change - 15 for a year of half-hourly files.
!
!              A change that goes and comes back between two files read at
!              the same rate, with no size step to point at it, is not found.
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine SurveyGhgAcFreq(FileList, nfiles, fmax)
    use m_rp_global_var
    use m_remote_source, only: RemoteSizeOf
    implicit none
    !> in/out variables
    integer, intent(in) :: nfiles
    type(FileListType), intent(in) :: FileList(nfiles)
    !> On entry the rate already known (the preamble file's), on exit the
    !> highest rate seen, never lower than on entry.
    real(kind = dbl), intent(inout) :: fmax
    !> local variables
    integer, parameter :: NumEven = 20
    integer, parameter :: MaxPeeks = 60
    !> For the halving, on top of MaxPeeks: 15 changes in a year of files
    integer, parameter :: MaxRefine = 240
    real(kind = dbl), parameter :: StepUp = 1.4d0
    real(kind = dbl), parameter :: StepDown = 0.7d0
    integer :: i, k, best
    integer :: a, m
    integer :: nrefine
    integer :: nread
    integer :: nchanges
    integer(kind = 8), allocatable :: fsize(:)
    real(kind = dbl), allocatable :: step(:)
    real(kind = dbl), allocatable :: rate(:)
    logical, allocatable :: pick(:)
    logical, allocatable :: tried(:)
    real(kind = dbl) :: ratio
    real(kind = dbl) :: prev
    logical :: found
    character(16) :: RateString
    character(10) :: date
    character(5) :: time
    character(16) :: CountString
    type(DateType) :: tsFrom


    if (nfiles < 2) return
    allocate(fsize(nfiles), step(nfiles), pick(nfiles), tried(nfiles), &
        rate(nfiles))
    pick = .false.
    tried = .false.
    step = 0d0
    rate = -1d0

    !> File sizes, and the size steps between neighbours. Only the larger file
    !> of a step is worth a read, since only the highest rate is wanted.
    !> From a shared link, the size is the listed one (Dropbox lists it,
    !> Google does not, and then there are no size steps to follow)
    do i = 1, nfiles
        fsize(i) = RemoteSizeOf(FileList(i)%path)
        if (fsize(i) < 0) inquire(file = FileList(i)%path, size = fsize(i))
    end do
    do i = 2, nfiles
        if (fsize(i) <= 0 .or. fsize(i - 1) <= 0) cycle
        ratio = dble(fsize(i)) / dble(fsize(i - 1))
        if (ratio > StepUp) then
            step(i) = log(ratio)
        else if (ratio < StepDown) then
            step(i - 1) = max(step(i - 1), -log(ratio))
        end if
    end do

    !> First, last and evenly spaced
    pick(1) = .true.
    pick(nfiles) = .true.
    do k = 1, NumEven
        i = 1 + nint(dble(k) * dble(nfiles - 1) / dble(NumEven + 1))
        pick(min(max(i, 1), nfiles)) = .true.
    end do

    !> Then the size steps, largest first, until the budget is spent
    do while (count(pick) < MaxPeeks)
        best = 0
        do i = 1, nfiles
            if (pick(i) .or. step(i) <= 0d0) cycle
            if (best == 0) then
                best = i
            else if (step(i) > step(best)) then
                best = i
            end if
        end do
        if (best == 0) exit
        pick(best) = .true.
    end do

    !> Read the picked files' rates
    do i = 1, nfiles
        if (pick(i)) call Peek(i)
    end do

    !> Halve every interval whose two ends disagree until they are
    !> neighbours. An unreadable file is skipped over: the untried file
    !> nearest the middle is read instead, and an interval with none left
    !> is as narrow as it gets.
    nrefine = 0
    do while (nrefine < MaxRefine)
        found = .false.
        a = 0
        do i = 1, nfiles
            if (rate(i) <= 0d0) cycle
            if (a > 0) then
                if (i > a + 1 .and. Differ(rate(a), rate(i))) then
                    m = UntriedNear((a + i) / 2, a, i)
                    if (m > 0) then
                        call Peek(m)
                        nrefine = nrefine + 1
                        found = .true.
                        exit
                    end if
                end if
            end if
            a = i
        end do
        if (.not. found) exit
    end do

    !> The list of changes. Nothing at all for a project at one rate.
    nchanges = 0
    prev = -1d0
    do i = 1, nfiles
        if (rate(i) <= 0d0) cycle
        fmax = max(fmax, rate(i))
        if (prev > 0d0 .and. Differ(prev, rate(i))) nchanges = nchanges + 1
        prev = rate(i)
    end do

    if (nchanges > 0) then
        call LogSay(' The raw files are not all at one acquisition frequency:')
        prev = -1d0
        do i = 1, nfiles
            if (rate(i) <= 0d0) cycle
            if (prev < 0d0 .or. Differ(prev, rate(i))) then
                !> When the new rate starts: the start of that file's data
                tsFrom = FileList(i)%timestamp
                if (EddyFlowLog%tstamp_end) tsFrom = tsFrom - DatafileDateStep
                call DateTypeToDateTime(tsFrom, date, time)
                write(RateString, '(f8.3)') rate(i)
                call LogSay('  ' // RateString(1:8) // ' Hz from ' // date // ' ' // time)
            end if
            prev = rate(i)
        end do
        write(RateString, '(f0.3)') fmax
        nread = count(tried)
        write(CountString, '(i0, a, i0)') nread, ' of ', nfiles
        call LogSay(' (Frequencies read from ' // trim(CountString) // ' files.)')
        call LogSay(' Each averaging period is processed at its own files'' frequency.')
        call LogSay(' Binned spectra are gridded up to the Nyquist frequency of ' &
            // trim(RateString) // ' Hz.')
    end if

    deallocate(fsize, step, pick, tried, rate)

contains

    !> One file's rate into rate(j)
    subroutine Peek(j)
        integer, intent(in) :: j
        real(kind = dbl) :: freq
        logical :: ok

        tried(j) = .true.
        call PeekGhgAcFreq(FileList(j), freq, ok)
        if (ok) rate(j) = freq
    end subroutine Peek

    logical function Differ(f1, f2)
        real(kind = dbl), intent(in) :: f1, f2

        Differ = abs(f1 - f2) > 1d-6 * max(f1, f2)
    end function Differ

    !> The untried file nearest mid, strictly between lo and hi; 0 if none
    integer function UntriedNear(mid, lo, hi)
        integer, intent(in) :: mid, lo, hi
        integer :: d

        UntriedNear = 0
        do d = 0, hi - lo
            if (mid - d > lo) then
                if (.not. tried(mid - d)) then
                    UntriedNear = mid - d
                    return
                end if
            end if
            if (mid + d < hi) then
                if (.not. tried(mid + d)) then
                    UntriedNear = mid + d
                    return
                end if
            end if
        end do
    end function UntriedNear
end subroutine SurveyGhgAcFreq

!***************************************************************************
!
! \brief       The acquisition frequency a GHG archive's .metadata states,
!              extracting that file alone and reading that line alone.
! \author      Jonathan Muller
! \note        Leaves every global alone - Metadata, Col and the rest - which
!              ReadMetadataFile would not.
!
!              The member is streamed out (7-Zip's -so) into a file of a fixed
!              name, rather than extracted under its own: an archive renamed
!              after it was written - which people do - no longer tells what
!              its members are called, and the biomet sibling, which states no
!              acquisition frequency of the raw data, is excluded by pattern.
!***************************************************************************
subroutine PeekGhgAcFreq(GhgFile, freq, ok)
    use m_rp_global_var
    use m_remote_source, only: RemoteEnsure
    implicit none
    !> in/out variables
    type(FileListType), intent(in) :: GhgFile
    real(kind = dbl), intent(out) :: freq
    logical, intent(out) :: ok
    !> local variables
    integer :: umeta
    integer :: io_status
    integer :: unzip_status
    integer :: del_status
    integer :: eq
    character(PathLen) :: MetaFile
    character(CommLen) :: comm
    character(LongInstringLen) :: dataline


    ok = .false.
    freq = -1d0

    !> Not *.metadata: UnZipArchive lists TmpDir and would take it for an
    !> archive's own if one were ever left behind
    MetaFile = trim(adjustl(TmpDir)) // 'ac_freq_peek.tmp'
    !> Out of processing order, so nothing is fetched ahead of it
    call RemoteEnsure(GhgFile%path, .false.)
    comm = trim(comm_7zip) // ' e -so "' // trim(GhgFile%path) &
        // '" "*.metadata" "-x!*-biomet.metadata" > "' // trim(MetaFile) // '"' &
        // comm_err_redirect
    unzip_status = system(trim(comm))

    if (unzip_status == 0) then
        open(newunit = umeta, file = MetaFile, status = 'old', action = 'read', &
            iostat = io_status)
        if (io_status == 0) then
            do
                read(umeta, '(a)', iostat = io_status) dataline
                if (io_status /= 0) exit
                dataline = adjustl(dataline)
                if (index(dataline, 'acquisition_frequency') /= 1) cycle
                eq = index(dataline, '=')
                if (eq == 0) cycle
                read(dataline(eq + 1:), *, iostat = io_status) freq
                ok = io_status == 0 .and. freq > 0d0
                exit
            end do
            close(umeta)
        end if
    end if

    del_status = system(trim(comm_del) // ' "' // trim(MetaFile) // '"' &
        // comm_err_redirect)
end subroutine PeekGhgAcFreq
