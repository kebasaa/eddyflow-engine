!***************************************************************************
! CRA.f90
! -------
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
! \brief       Applies centered running average smoothing.
! \author      Jonathan Muller
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine CRA(Set, nrow, ncol, fs, tconst, fcol)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow
    integer, intent(in) :: ncol
    integer, intent(in) :: tconst
    integer, intent(in) :: fcol
    real(kind = dbl), intent(in) :: fs
    real(kind = dbl), intent(inout) :: Set(nrow, ncol)
    !> local variables
    integer :: i
    integer :: np
    integer :: nnp
    real(kind = dbl) tmp(nrow)


    np = nint(fs * tconst)

    !> Centered average for first part of array
    do i = 1, np/2
        nnp = (i-1)*2
        call AverageNoError(Set(i-nnp/2:i+nnp/2, fcol), nnp+1, 1, tmp(i), error)
    end do
    !> Centered average for inner points
    do i = np/2+1, nrow-np/2
        call AverageNoError(Set(i-np/2:i+np/2, fcol), np+1, 1, tmp(i), error)
    end do
    !> Centered average for last part of array
    do i = nrow-np/2+1, nrow
        nnp = (nrow-i)*2
        call AverageNoError(Set(i-nnp/2:i+nnp/2, fcol), nnp+1, 1, tmp(i), error)
    end do

    Set(:, fcol) = tmp
end subroutine CRA
