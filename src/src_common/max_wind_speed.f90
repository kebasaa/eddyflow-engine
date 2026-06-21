!***************************************************************************
! max_wind_speed.f90
! ------------------
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
! \brief       Calculate maximal wind speed in a 3d wind array
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine MaxWindSpeed(Set, nrow, ncol, MaxSpeed)
    use m_common_global_var
    implicit none
    !> In/out variables
    integer, intent(in) :: nrow, ncol
    real(kind = dbl) , intent(in) :: Set(nrow, ncol)
    real(kind = dbl) , intent(out) :: MaxSpeed
    !> Local variables
    integer :: i
    real(kind = dbl)  :: CurrentSpeed

    MaxSpeed = 0d0
    do i = 1, nrow
        CurrentSpeed = dsqrt(Set(i, u)**2 + Set(i, v)**2 + Set(i, w)**2)
        if (Set(i, u) /= error .and. Set(i, v) /= error &
            .and. Set(i, w) /= error .and. &
         CurrentSpeed > MaxSpeed) MaxSpeed = CurrentSpeed
    end do

end subroutine MaxWindSpeed
