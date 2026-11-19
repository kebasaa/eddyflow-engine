!***************************************************************************
! power_of_two.f90
! ----------------
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
! \brief       Calculates closest power of two to N
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine PowerOfTwo(N, n2)
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    integer, intent(out) :: n2
    !> local variables
    integer :: i = 0

    loop: do i = 1, 25
        if((2**i) > N) then
            n2 = 2**(i - 1)
            exit loop
        end if
    end do loop
end subroutine PowerOfTwo

!***************************************************************************
!
! \brief       The largest number of samples up to N with no prime factor
!              above 13
! \author      Jonathan Muller, ETH Zurich
! \note        The FFT (FFTPACK, src_rp/fft4) has fast passes for factors 2,
!              3, 4 and 5 and a generic one for any other, whose cost grows
!              with the factor. Measured with this FFT, per transform: 32768
!              rows 0.12 ms, 36000 0.14 ms, 35750 (2 5 5 5 11 13) 0.21 ms,
!              35760 (36000 less a 12 s lag at 20 Hz, with a factor of 149)
!              0.33 ms, 35761 (11 x 3251) 11 ms, 35759 (prime) 865 ms. So
!              only a large prime factor costs anything that matters.
!
!              A power of two avoids that by dropping up to half the period -
!              32768 of 36000 rows, 9 %, the last 2.7 minutes of a half-hour.
!              Allowing factors up to 13 keeps within 1.5 times its speed and
!              loses at most 1 % of the rows between 10000 and 72000 (0.15 %
!              on average): 36000 stays 36000, 35760 becomes 35750.
!***************************************************************************
subroutine FftFriendlyLength(N, n2)
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    integer, intent(out) :: n2
    !> local variables
    integer :: m
    integer :: k
    integer, parameter :: primes(6) = [2, 3, 5, 7, 11, 13]

    n2 = max(N, 1)
    do while (n2 > 1)
        m = n2
        do k = 1, size(primes)
            do while (mod(m, primes(k)) == 0)
                m = m / primes(k)
            end do
        end do
        if (m == 1) return
        n2 = n2 - 1
    end do
end subroutine FftFriendlyLength
