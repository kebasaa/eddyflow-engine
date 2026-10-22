!***************************************************************************
! tlag_preassess.f90
! ------------------
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
! \brief       The preliminary time-lag assessment: a random sample of the
!              periods, and a table of every lag found in them.
!
!              The interface's "Timelag preliminary assessment" dialog runs
!              a throw-away copy of the project as an assessment-only time-lag
!              pre-pass over a wide search window, so the user can see where
!              the lags actually fall before choosing the nominal, minimum and
!              maximum lags in the metadata. Two things that pre-pass did not
!              do are here:
!
!              - walk a share of the periods rather than all of them
!                (to_sample_frac, at least to_sample_min_n, to_sample_seed).
!                The draw is from the periods a raw file covers, and is the
!                same in every process of the run, so worker pieces skip
!                exactly what the parent would have;
!              - write every period's detected lag, before the optimiser's
!                flux and window-edge gates (to_dump_lags), because a pile-up
!                at a window edge is precisely what the user needs to see.
!
!              None of these keys is ever written into a user's project.
!
! \author      Jonathan Muller
! \sa          batch_pool.f90 (SetPeriodSample), add_to_timelag_opt_dataset.f90
!***************************************************************************
module m_tlag_preassess
    use m_rp_global_var
    use m_batch_pool, only: PeriodsWithFile, SetPeriodSample
    use m_log, only: LogSay
    use m_pwb_timelag, only: GasLabel
    implicit none
    private

    public :: DrawTlagSample, WriteTlagPreassessTable, WritePwbPreassessTable

    !> Park and Miller's minimal standard generator. Not for its statistics -
    !> any draw will do - but because it is the same sequence on every
    !> compiler and every run, and its products fit a 64-bit integer.
    integer(8), parameter :: PmModulus = 2147483647_8
    integer(8), parameter :: PmMultiplier = 16807_8

    character(15), parameter :: Preassess_FilePadding = '_tlag_preassess'

contains

    !***************************************************************************
    !> \brief Draw the sample for the pre-pass over [iStart, iEnd) and set it.
    !>
    !> n = max(sample_min_n, round(sample_frac * periods with a file)), no more
    !> than there are. Does nothing when no sample was asked for, so the
    !> caller need not ask. nWalk is how many periods the pass will walk.
    !***************************************************************************
    subroutine DrawTlagSample(iStart, iEnd, Series, nSeries, Files, nFiles, nWalk)
        integer, intent(in) :: iStart
        integer, intent(in) :: iEnd
        integer, intent(in) :: nSeries
        integer, intent(in) :: nFiles
        type(DateType), intent(in) :: Series(nSeries)
        type(FileListType), intent(in) :: Files(nFiles)
        integer, intent(out) :: nWalk
        logical, allocatable :: has(:)
        logical, allocatable :: keep(:)
        integer, allocatable :: pool(:)
        integer :: nAvail
        integer :: n
        integer :: i
        integer :: j
        integer :: p
        integer :: tmp
        integer(8) :: state
        character(16) :: s1
        character(16) :: s2

        nWalk = iEnd - iStart
        if (TOSetup%sample_frac <= 0d0 .or. iEnd <= iStart) return

        allocate(has(iStart:iEnd - 1))
        call PeriodsWithFile(iStart, iEnd, Series, nSeries, Files, nFiles, has)
        nAvail = count(has)
        allocate(pool(max(1, nAvail)))
        n = 0
        do p = iStart, iEnd - 1
            if (.not. has(p)) cycle
            n = n + 1
            pool(n) = p
        end do

        n = max(TOSetup%sample_min_n, nint(TOSetup%sample_frac * dble(nAvail)))
        n = min(n, nAvail)

        !> A partial Fisher-Yates shuffle: the first n of the pool, drawn
        !> without replacement.
        state = mod(abs(int(TOSetup%sample_seed, 8)), PmModulus - 1_8) + 1_8
        do i = 1, n
            state = mod(PmMultiplier * state, PmModulus)
            j = i + int(mod(state, int(nAvail - i + 1, 8)))
            tmp = pool(i)
            pool(i) = pool(j)
            pool(j) = tmp
        end do

        allocate(keep(iEnd - iStart))
        keep = .false.
        do i = 1, n
            keep(pool(i) - iStart + 1) = .true.
        end do
        call SetPeriodSample(iStart, keep)
        nWalk = n

        write(s1, '(i0)') n
        write(s2, '(i0)') nAvail
        call LogSay('  Preliminary time-lag assessment: ' // trim(s1) &
            // ' of ' // trim(s2) // ' periods with raw data, drawn at random.')
        deallocate(has, keep, pool)
    end subroutine DrawTlagSample

    !***************************************************************************
    !> \brief The table of a covariance-maximisation pre-pass, one row per
    !>        period and gas that had a lag.
    !***************************************************************************
    subroutine WriteTlagPreassessTable(rows, nrow, n)
        integer, intent(in) :: nrow
        integer, intent(in) :: n
        type(TimeLagOptType), intent(in) :: rows(nrow)
        integer :: u
        integer :: i
        integer :: gas
        logical :: ok

        if (.not. TOSetup%dump_lags) return
        call OpenTable(u, ok)
        if (.not. ok) return
        do i = 1, n
            do gas = firstGas, lastGas
                if (rows(i)%RawLag(gas) == error) cycle
                call WriteRow(u, rows(i)%date, rows(i)%time, gas, 'maxcov', &
                    rows(i)%RawLag(gas), rows(i)%WinMin(gas), rows(i)%WinMax(gas), &
                    rows(i)%AtEdge(gas), rows(i)%Flux(gas), rows(i)%FluxOk(gas))
            end do
        end do
        call CloseTable(u)
    end subroutine WriteTlagPreassessTable

    !***************************************************************************
    !> \brief The same table from a PWB pre-pass, out of the half-hourly PWB
    !>        table: each period's own detection, before any gap-filling.
    !>
    !> PWB judges no flux, so every row passes that gate, and the flux column
    !> is left empty.
    !***************************************************************************
    subroutine WritePwbPreassessTable()
        integer :: u
        integer :: i
        logical :: ok

        if (.not. TOSetup%dump_lags) return
        call OpenTable(u, ok)
        if (.not. ok) return
        do i = 1, PwbTimelagCacheN
            if (PwbTimelagCache(i)%result%selected_lag == error) cycle
            call WriteRow(u, PwbTimelagCache(i)%date, PwbTimelagCache(i)%time, &
                PwbTimelagCache(i)%gas, 'pwb', &
                PwbTimelagCache(i)%result%selected_lag, &
                PwbTimelagCache(i)%result%effective_min_lag, &
                PwbTimelagCache(i)%result%effective_max_lag, &
                PwbTimelagCache(i)%result%edge_pinned, error, .true.)
        end do
        call CloseTable(u)
    end subroutine WritePwbPreassessTable

    subroutine OpenTable(u, ok)
        integer, intent(out) :: u
        logical, intent(out) :: ok
        integer :: ios
        character(PathLen) :: path

        ok = .false.
        u = -1
        if (Dir%main_out == 'none') return
        path = Dir%main_out(1:len_trim(Dir%main_out)) &
            // EddyFlowProj%id(1:len_trim(EddyFlowProj%id)) &
            // Preassess_FilePadding // Timestamp_FilePadding // CsvExt
        open(newunit = u, file = path, status = 'replace', iostat = ios, &
            encoding = 'utf-8')
        if (ios /= 0) return
        !> gas_record is the project's gas_<N>_* record of the gas, whose
        !> gas_<N>_col names its metadata column - which is how the interface
        !> finds the column to write the chosen lags to.
        write(u, '(a)') 'date,time,gas,gas_record,method,detected_lag_s,' &
            // 'window_min_s,window_max_s,at_edge,flux,flux_ok'
        call LogSay(' Preliminary time-lag table written to: ' // trim(path))
        ok = .true.
    end subroutine OpenTable

    subroutine CloseTable(u)
        integer, intent(in) :: u
        close(u)
    end subroutine CloseTable

    subroutine WriteRow(u, date, time, gas, method, lag, wmin, wmax, edge, flux, fluxOk)
        integer, intent(in) :: u
        character(*), intent(in) :: date
        character(*), intent(in) :: time
        integer, intent(in) :: gas
        character(*), intent(in) :: method
        real(kind = dbl), intent(in) :: lag
        real(kind = dbl), intent(in) :: wmin
        real(kind = dbl), intent(in) :: wmax
        logical, intent(in) :: edge
        real(kind = dbl), intent(in) :: flux
        logical, intent(in) :: fluxOk
        character(32) :: fluxString

        fluxString = ''
        if (flux /= error) write(fluxString, '(es14.6)') flux
        write(u, '(a,",",a,",",a,",",i0,",",a,",",f10.4,",",f10.4,",",f10.4,' &
            // '",",i1,",",a,",",i1)') &
            trim(date), trim(time), trim(GasLabel(gas)), gas - firstGas + 1, &
            trim(method), lag, wmin, wmax, merge(1, 0, edge), &
            trim(adjustl(fluxString)), merge(1, 0, fluxOk)
    end subroutine WriteRow
end module m_tlag_preassess
