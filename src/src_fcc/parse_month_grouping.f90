!***************************************************************************
! parse_month_grouping.f90
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
! \brief       Which months a gas pools before a transfer function is fitted.
! \author      Jonathan Muller
! \note        A gas states its grouping as `gas_<i>_sa_months`, a list of
!              month ranges: `1-12` is one group over the calendar, `1-6,7-12`
!              is two. A group's ordinal in the list is its class index, which
!              is what indexes RegPar and MeanBinSpec.
!
!              This replaces three flat tables of twelve start/stop pairs -
!              one each for CO2, CH4 and the fourth gas - which is why every
!              gas past the fourth had to inherit CO2's grouping.
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************

!***************************************************************************
!
! \brief       Read a month-grouping list into a per-month class map.
! \author      Jonathan Muller
! \note        All-or-nothing. A list the parser cannot read in full leaves
!              `cls` zeroed and `ok` false, rather than applying the groups it
!              managed to read: a half-applied grouping classes some months
!              and silently leaves the rest at 0, which reads downstream as
!              "this gas has no fitted cutoff for those months" and is the
!              worst of both answers.
!
!              The caller treats a refusal exactly as it treats an absent key.
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine ParseMonthGrouping(spec, cls, ncls, ok)
    use m_common_global_var
    implicit none
    !> in/out variables
    character(*), intent(in) :: spec
    !> Indexed by month; the value is the group's ordinal, 0 for a month in
    !> no group. Twelve because there are twelve months, not because of any
    !> capacity constant.
    integer, intent(out) :: cls(12)
    integer, intent(out) :: ncls
    logical, intent(out) :: ok
    !> local variables
    character(len(spec)) :: body
    character(len(spec)) :: token
    integer :: pos
    integer :: nxt
    integer :: dash
    integer :: first_month
    integer :: last_month
    integer :: month
    integer :: read_status

    cls = 0
    ncls = 0
    ok = .false.

    body = adjustl(spec)
    if (len_trim(body) == 0) return

    pos = 1
    do while (pos <= len_trim(body))
        nxt = index(body(pos:len_trim(body)), ',')
        if (nxt == 0) then
            token = body(pos:len_trim(body))
            pos = len_trim(body) + 1
        else
            token = body(pos:pos + nxt - 2)
            pos = pos + nxt
        end if
        token = adjustl(token)
        if (len_trim(token) == 0) return

        !> `3` is accepted as `3-3`. Readers widen, writers narrow: the
        !> interface always emits the range form, but the per-month grouping
        !> spelled out in full is a real configuration and the shorthand
        !> halves it.
        dash = index(trim(token), '-')
        if (dash == 0) then
            read(token, *, iostat = read_status) first_month
            if (read_status /= 0) return
            last_month = first_month
        else
            if (dash == 1 .or. dash == len_trim(token)) return
            read(token(1:dash - 1), *, iostat = read_status) first_month
            if (read_status /= 0) return
            read(token(dash + 1:len_trim(token)), *, iostat = read_status) last_month
            if (read_status /= 0) return
        end if

        if (first_month < 1 .or. last_month > 12) return
        if (first_month > last_month) return

        ncls = ncls + 1
        !> Twelve months admit at most twelve groups, so this bound is the
        !> calendar's, not a capacity choice that happens to coincide.
        if (ncls > MaxGasClasses) return

        do month = first_month, last_month
            !> A month in two groups has no answer to "which class is it".
            if (cls(month) /= 0) return
            cls(month) = ncls
        end do
    end do

    if (ncls == 0) return
    ok = .true.
end subroutine ParseMonthGrouping
