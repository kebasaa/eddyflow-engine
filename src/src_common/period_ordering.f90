!***************************************************************************
! period_ordering.f90
! --------------------
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
! \brief       Chronological ordering for a whole-run per-period cache:
!              sort by PeriodMinutes, then find the modal spacing between
!              consecutive periods.
!
! \details     Shared by every whole-run post-pass that accumulates one row
!              per averaging period and then needs it in time order with a
!              periods-per-day count - pfd_handle.f90's post-flux despiking
!              (FCC) and stor_clean_handle.f90's storage-flux cleaning (RP).
!              Originally written once inside pfd_handle.f90; moved here
!              when the RP side needed the same ordering rather than a
!              second copy that could drift from it.
!
! \author      Jonathan Muller, ETH Zurich
! \note
! \sa          period_minutes.f90, pfd_handle.f90, stor_clean_handle.f90
!***************************************************************************

!> Insertion-sort ord(1:n) by tmin, ascending - the caches this serves are
!> one row per averaging period, one run, so O(n^2) costs nothing here.
subroutine SortPeriodsByMinutes(ord, tmin, n)
    implicit none
    integer, intent(in) :: n
    integer, intent(inout) :: ord(n)
    integer(8), intent(in) :: tmin(n)
    integer :: i, j, key
    integer(8) :: keyval

    do i = 2, n
        key = ord(i)
        keyval = tmin(key)
        j = i - 1
        do while (j >= 1)
            if (tmin(ord(j)) <= keyval) exit
            ord(j+1) = ord(j)
            j = j - 1
        end do
        ord(j+1) = key
    end do
end subroutine SortPeriodsByMinutes


!> The most common step between consecutive sorted periods. A simple
!> "first minus second" would be wrong if the run's very first gap happens
!> at the start; this instead counts how often each observed step occurs
!> and keeps the winner, which is robust to any number of gaps as long as
!> most periods are still contiguous.
integer(8) function ModalPeriodStep(tmin, ord, n)
    implicit none
    integer, intent(in) :: n
    integer(8), intent(in) :: tmin(n)
    integer, intent(in) :: ord(n)
    integer(8), allocatable :: steps(:), uniq(:)
    integer, allocatable :: cnt(:)
    integer :: i, j, nu, best_i

    ModalPeriodStep = 0_8
    if (n < 2) return

    allocate(steps(n-1))
    do i = 1, n-1
        steps(i) = tmin(ord(i+1)) - tmin(ord(i))
    end do

    allocate(uniq(n-1), cnt(n-1))
    nu = 0
    do i = 1, n-1
        if (steps(i) <= 0_8) cycle
        do j = 1, nu
            if (uniq(j) == steps(i)) then
                cnt(j) = cnt(j) + 1
                exit
            end if
        end do
        if (j > nu) then
            nu = nu + 1
            uniq(nu) = steps(i)
            cnt(nu) = 1
        end if
    end do

    if (nu > 0) then
        best_i = 1
        do j = 2, nu
            if (cnt(j) > cnt(best_i)) best_i = j
        end do
        ModalPeriodStep = uniq(best_i)
    end if

    deallocate(steps, uniq, cnt)
end function ModalPeriodStep
