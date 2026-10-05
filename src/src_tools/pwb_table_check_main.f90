!***************************************************************************
! pwb_table_check_main.f90
! ------------------------
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
! EddyFlow® is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
!***************************************************************************
!
! \brief       Proves the PWB table finds the row a scan from row 1 finds,
!              and times storing and looking up a year of it.
!
!              The table's store and lookup used to scan from row 1, and every
!              append reallocated it to one row more. Both now take a faster
!              path - bisection while the rows are in period order, a capacity
!              that doubles - and this checks the one thing that matters about
!              that: for every key, the row returned is the first row a scan
!              from row 1 would return. Every row is given a distinct
!              actual_lag, so "the same values" means "the same row".
!
!              Cases: a year of four gases stored period by period the way the
!              pre-pass stores them; every key looked up, plus keys that are
!              not there; a store to an existing key (it must update the first
!              match, not append); duplicated keys appended in order, which
!              the bisection must resolve to the first; and rows appended out
!              of order, which must turn bisection off.
!
! \note        Built by `mingw32-make pwbtablecheck`, not in `all`, linked
!              against every engine object but the main program, so it runs
!              the shipped routines. static_checks/test_pwb_table_static.py
!              runs it.
! \author      Jonathan Muller
!***************************************************************************
program pwb_table_check_main
    use, intrinsic :: iso_fortran_env, only: int64
    use m_rp_global_var, only: dbl, PWBTimelagCacheEntryType, PWBResultType, &
        PwbTimelagCache, PwbTimelagCacheN, PwbCacheInTimeOrder, PwbCacheLoaded, &
        PwbPeriodDate, PwbPeriodTime
    use m_pwb_timelag, only: InitPwbTimelagCache, StorePwbTimelagCache, &
        LookupPwbTimelagCache, SetPwbPeriodTimestamp, AppendPwbCacheRows, InitPwbResult
    implicit none
    integer, parameter :: ngas = 4
    integer, parameter :: nper = 17520
    integer :: iper, g, bad, checked, before
    integer(int64) :: c0, c1, rate
    real(kind = dbl) :: t_store, t_look
    character(10) :: date
    character(5) :: time
    type(PWBTimelagCacheEntryType) :: extra(6)

    bad = 0
    checked = 0
    call InitPwbTimelagCache()
    PwbCacheLoaded = .true.

    !> A year of half-hours, four gases each, stored as the pre-pass does.
    call system_clock(c0, rate)
    do iper = 1, nper
        call Stamp(iper, date, time)
        call SetPwbPeriodTimestamp(date, time)
        do g = 1, ngas
            call StoreOne(g, dble((iper - 1) * ngas + g))
        end do
    end do
    call system_clock(c1)
    t_store = dble(c1 - c0) / dble(rate)
    if (PwbTimelagCacheN /= nper * ngas) call Fail('row count after storing')
    if (.not. PwbCacheInTimeOrder) call Fail('a table stored in period order is not marked ordered')

    !> Every key, in period order as the production pass asks: timed alone,
    !> then checked against a scan from row 1 - which is itself quadratic, so
    !> it is kept out of the timing.
    call system_clock(c0, rate)
    do iper = 1, nper
        call Stamp(iper, date, time)
        call SetPwbPeriodTimestamp(date, time)
        do g = 1, ngas
            call LookupOnly(g)
        end do
    end do
    call system_clock(c1)
    t_look = dble(c1 - c0) / dble(rate)
    !> Every key again, and a gas that was never stored, against the scan.
    do iper = 1, nper
        call Stamp(iper, date, time)
        call SetPwbPeriodTimestamp(date, time)
        do g = 1, ngas + 1
            call CheckOne(g)
        end do
    end do
    !> A period after the last and one before the first.
    call SetPwbPeriodTimestamp('2027-01-01', '00:00')
    call CheckOne(1)
    call SetPwbPeriodTimestamp('2000-01-01', '00:00')
    call CheckOne(1)

    !> Storing to a key that exists updates that row and appends nothing.
    call Stamp(5000, date, time)
    call SetPwbPeriodTimestamp(date, time)
    before = PwbTimelagCacheN
    call StoreOne(2, -7d0)
    if (PwbTimelagCacheN /= before) call Fail('a store to an existing key appended a row')
    call CheckOne(2)

    !> Duplicates, appended in order at the end: the bisection must return
    !> the first of them, as a scan would.
    call Stamp(nper + 1, date, time)
    extra(1:3) = MakeRow(date, time, 1, 1d6)
    extra(2)%actual_lag = 1d6 + 1
    extra(3)%actual_lag = 1d6 + 2
    extra(4) = MakeRow(date, time, 2, 1d6 + 3)
    call AppendPwbCacheRows(extra(1:4), 4)
    if (.not. PwbCacheInTimeOrder) call Fail('in-order duplicates cleared the order flag')
    call SetPwbPeriodTimestamp(date, time)
    call CheckOne(1)
    call CheckOne(2)

    !> And the first period stored again, at the end, out of order: bisection
    !> has to stop, and the scan finds the original row.
    call Stamp(1, date, time)
    extra(5) = MakeRow(date, time, 1, 2d6)
    call AppendPwbCacheRows(extra(5:5), 1)
    if (PwbCacheInTimeOrder) call Fail('an out-of-order row left the order flag set')
    call SetPwbPeriodTimestamp(date, time)
    call CheckOne(1)
    call Stamp(nper + 1, date, time)
    call SetPwbPeriodTimestamp(date, time)
    call CheckOne(1)

    write(*, '(a,i0,a,i0)') ' keys checked: ', checked, '   wrong rows: ', bad
    write(*, '(a,i0,a,f8.3,a,f8.3,a)') ' a year of ', nper * ngas, ' rows: stored in ', &
        t_store, ' s, every key looked up in ', t_look, ' s'
    if (bad == 0) then
        write(*, '(a)') ' SAME ROW AS A SCAN FROM ROW 1'
    else
        write(*, '(a)') ' MISMATCH'
        error stop 1
    end if

contains

    subroutine Stamp(k, d, t)
        integer, intent(in) :: k
        character(10), intent(out) :: d
        character(5), intent(out) :: t
        integer :: day, half, month, mday
        integer, parameter :: mlen(12) = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

        day = (k - 1) / 48
        half = mod(k - 1, 48)
        month = 1
        mday = day
        do while (mday >= mlen(month))
            mday = mday - mlen(month)
            month = month + 1
            if (month > 12) then
                month = 1
            end if
        end do
        if (k > 365 * 48) then
            write(d, '(a,i2.2,a,i2.2)') '2020-', month, '-', mday + 1
        else
            write(d, '(a,i2.2,a,i2.2)') '2019-', month, '-', mday + 1
        end if
        write(t, '(i2.2,a,i2.2)') half / 2, ':', 30 * mod(half, 2)
    end subroutine Stamp

    type(PWBTimelagCacheEntryType) function MakeRow(d, t, g, lag)
        character(*), intent(in) :: d, t
        integer, intent(in) :: g
        real(kind = dbl), intent(in) :: lag

        MakeRow = PwbTimelagCache(1)
        MakeRow%date = d
        MakeRow%time = t
        MakeRow%gas = g
        MakeRow%actual_lag = lag
    end function MakeRow

    subroutine StoreOne(g, lag)
        integer, intent(in) :: g
        real(kind = dbl), intent(in) :: lag
        type(PWBResultType) :: res

        call InitPwbResult(res)
        call StorePwbTimelagCache(g, lag, lag, 0, .false., res)
    end subroutine StoreOne

    subroutine LookupOnly(g)
        integer, intent(in) :: g
        logical :: found, default_used
        real(kind = dbl) :: actual_lag, used_lag
        integer :: row_lag
        type(PWBResultType) :: res

        call LookupPwbTimelagCache(g, found, actual_lag, used_lag, row_lag, default_used, res)
        if (.not. found) call Miss('missed a key in the timed pass', g)
    end subroutine LookupOnly

    !> Lookup against a scan from row 1 over the same table.
    subroutine CheckOne(g)
        integer, intent(in) :: g
        logical :: found, default_used
        real(kind = dbl) :: actual_lag, used_lag
        integer :: row_lag, i, first
        type(PWBResultType) :: res

        first = 0
        do i = 1, PwbTimelagCacheN
            if (PwbTimelagCache(i)%date == PwbPeriodDate .and. PwbTimelagCache(i)%time == PwbPeriodTime &
                .and. PwbTimelagCache(i)%gas == g) then
                first = i
                exit
            end if
        end do
        call LookupPwbTimelagCache(g, found, actual_lag, used_lag, row_lag, default_used, res)
        checked = checked + 1
        if (first == 0) then
            if (found) call Miss('found a key that is not there', g)
        else
            if (.not. found) then
                call Miss('missed a key that is there', g)
            else if (transfer(actual_lag, 0_int64) /= &
                     transfer(PwbTimelagCache(first)%actual_lag, 0_int64)) then
                call Miss('returned a row other than the first match', g)
            end if
        end if
    end subroutine CheckOne

    subroutine Miss(what, g)
        character(*), intent(in) :: what
        integer, intent(in) :: g
        if (bad < 10) write(*, '(a,a,a,a,a,i0)') ' ', what, ': ', &
            PwbPeriodDate // ' ' // PwbPeriodTime, ' gas ', g
        bad = bad + 1
    end subroutine Miss

    subroutine Fail(what)
        character(*), intent(in) :: what
        write(*, '(a,a)') ' FAIL: ', what
        error stop 1
    end subroutine Fail

end program pwb_table_check_main
