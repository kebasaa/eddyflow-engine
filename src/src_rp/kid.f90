!***************************************************************************
! kid.f90
! -------
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
! \brief       KID test: kurtosis on stochastically detrended data
! \author      Jonathan Muller, ETH Zurich
!
!***************************************************************************
! Reference: Vitale, D. et al. (2020). Biogeosciences Discussion.
!            KID = kurtosis of stochastic-residual; ZCD = zero-crossing density.
subroutine KID(Set, nrow, ncol)
    use m_rp_global_var
    implicit none
    integer, intent(in) :: nrow, ncol
    real(kind = dbl), intent(in) :: Set(nrow, ncol)
    integer :: icol
    integer, external :: CountZeroCrossings
    real(kind = dbl) :: residuals(nrow, ncol)

    do icol = u, ts
        call VariableStochasticDetrending(Set(:, icol), residuals(:, icol), nrow)
        call KurtosisNoError(residuals(:, icol), nrow, 1, Essentials%KID(icol), error)
        Essentials%ZCD(icol) = CountZeroCrossings(residuals(:, icol), nrow)
    end do
    do icol = co2, gas4
        if (E2Col(icol)%present) then
            call VariableStochasticDetrending(Set(:, icol), residuals(:, icol), nrow)
            call KurtosisNoError(residuals(:, icol), nrow, 1, Essentials%KID(icol), error)
            Essentials%ZCD(icol) = CountZeroCrossings(residuals(:, icol), nrow)
        else
            Essentials%KID(icol) = error
            Essentials%ZCD(icol) = ierror
        end if
    end do
end subroutine KID


integer function CountZeroCrossings(arr, nrow)
    use m_rp_global_var
    implicit none
    integer, intent(in) :: nrow
    real(kind = dbl), intent(in) :: arr(nrow)
    integer :: irow, current_sign, previous_sign

    CountZeroCrossings = 0
    previous_sign = 0
    do irow = 1, nrow
        if (arr(irow) == error .or. arr(irow) == 0d0) cycle

        if (arr(irow) > 0d0) then
            current_sign = 1
        else
            current_sign = -1
        end if

        if (previous_sign /= 0 .and. current_sign /= previous_sign) &
            CountZeroCrossings = CountZeroCrossings + 1
        previous_sign = current_sign
    end do
end function CountZeroCrossings
