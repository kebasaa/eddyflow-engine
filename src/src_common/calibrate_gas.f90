!***************************************************************************
! calibrate_gas.f90
! ------------------
! Copyright © 2007-2011, Eco2s team, Gerardo Fratini
! Copyright © 2011-2026, LI-COR Biosciences, Gerardo Fratini
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
! \brief       Calibrate 4th gas if a cal-ref column is available.
!              Note that so far the calibration procedure is fully customized
!              on the needs of a specific O3 analyzer
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine CalibrateGases(Set, nrow, ncol)
    use m_common_global_var
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    real(kind = dbl), intent(inout) :: Set(nrow, ncol)
    !> local variables
    integer :: j
    integer :: slot

    !> One calibration per reference column, against the gas that column
    !> names. This scaled Set(:, gas4) whichever gas the reference belonged
    !> to, because there could be only one reference and it was assumed to be
    !> the fourth slot's - so a site calibrating anything else had its fourth
    !> gas silently rescaled by an unrelated reference, and the gas it meant
    !> to calibrate left alone.
    do j = 1, NumUserVar
        if (UserCol(j)%var /= 'cal-ref') cycle
        slot = UserCalRefSlot(j)
        if (slot < firstGas .or. slot > lastGas) cycle
        if (slot > ncol) cycle
        if (Stats%Mean(slot) == error .or. Stats%Mean(slot) == 0d0) cycle
        if (UserStats%Mean(j) == error) cycle
        !> Converts mV in ppb and then to ppm (with "/ 1d3")
        Set(:, slot) = Set(:, slot) * UserStats%Mean(j) &
            / Stats%Mean(slot) / 1d3
    end do
end subroutine CalibrateGases
