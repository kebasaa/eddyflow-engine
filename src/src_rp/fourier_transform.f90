!***************************************************************************
! fourier_transform.f90
! ---------------------
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
! \brief       Provides a driver to rfftf. Fourier-transform the in/out array "xx".
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine FourierTransform(xx, N, M)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    integer, intent(in) :: M
    real(kind = dbl), intent(inout) :: xx(N, M)
    !> local variables
    integer :: i
    real :: xxx(N)
    real :: wsave(N*2 + 15)

    write(*, '(a)', advance = 'no') '   FFT-ing..'
    call rffti(N, wsave)
    do i = 1, M
        !> data in 1D vector
        xxx(:) = sngl(xx(:, i))
        !> fast fourier transform
        call rfftf(N, xxx, wsave)
        !> replace time data with spectral data
        xx(:, i) = dble(xxx(:))
    end do
    write(*,'(a)') ' Done.'
end subroutine FourierTransform
