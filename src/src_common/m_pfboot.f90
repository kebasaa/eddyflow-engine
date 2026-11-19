!***************************************************************************
! m_pfboot.f90
! ------------
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
! \brief       Bootstrap uncertainty of the planar-fit tilt angles, as innFLUX
!              reports it for each wind direction (Striednig et al. 2020,
!              Fig. 1; innFLUX_step0.m, 1000 resamples), on the planar fit of
!              Wilczak et al. (2001).
!
!              The periods of a sector are resampled with replacement, the
!              plane w = b0 + b1 u + b2 v (or, without velocity bias, w =
!              b1 u + b2 v) is fitted again to each resample, and the standard
!              deviations of the pitch alpha = asin(p31) and roll beta =
!              asin(-p32 / sqrt(p32^2 + p33^2)) over the resamples are the
!              uncertainties. p3j is the last row of the rotation matrix,
!              (-b1, -b2, 1) / sqrt(b1^2 + b2^2 + 1) (Wilczak et al. 2001, Eq. 42;
!              innFLUX calc_tilt_correction.m).
!
!              The generator is xorshift64 (Marsaglia 2003), seeded per sector
!              from the project's seed, so the result does not depend on the
!              order sectors are fitted in, nor on the platform.
!
! \note        Pure, and a function of its arguments only.
! \author      Jonathan Muller
!***************************************************************************
module m_pfboot
    use m_numeric_kinds
    implicit none
    private

    public :: TiltAnglesDeg
    public :: FitPlane
    public :: BootstrapTilt

contains

    !> Pitch and roll in degrees from the plane slopes b1, b2.
    pure subroutine TiltAnglesDeg(b1, b2, alpha, beta)
        real(kind = dbl), intent(in) :: b1
        real(kind = dbl), intent(in) :: b2
        real(kind = dbl), intent(out) :: alpha
        real(kind = dbl), intent(out) :: beta
        real(kind = dbl), parameter :: rad2deg = 57.295779513082320876d0
        real(kind = dbl) :: aux
        real(kind = dbl) :: p31
        real(kind = dbl) :: p32
        real(kind = dbl) :: p33

        aux = sqrt(b1 ** 2 + b2 ** 2 + 1d0)
        p31 = -b1 / aux
        p32 = -b2 / aux
        p33 = 1d0 / aux
        alpha = asin(p31) * rad2deg
        beta = asin(-p32 / sqrt(p32 ** 2 + p33 ** 2)) * rad2deg
    end subroutine TiltAnglesDeg

    !***************************************************************************
    !> \brief Least-squares plane through the points idx(1:m) of (u, v, w):
    !>        b = (b0, b1, b2); b0 = 0 when no_bias. ok false when singular.
    !***************************************************************************
    pure subroutine FitPlane(u, v, w, n, idx, m, no_bias, b, ok)
        integer, intent(in) :: n
        integer, intent(in) :: m
        real(kind = dbl), intent(in) :: u(n)
        real(kind = dbl), intent(in) :: v(n)
        real(kind = dbl), intent(in) :: w(n)
        integer, intent(in) :: idx(m)
        logical, intent(in) :: no_bias
        real(kind = dbl), intent(out) :: b(3)
        logical, intent(out) :: ok
        real(kind = dbl) :: a(3, 3)
        real(kind = dbl) :: r(3)
        real(kind = dbl) :: x(3)
        integer :: k
        integer :: i
        integer :: j
        integer :: p
        integer :: nd
        real(kind = dbl) :: f
        real(kind = dbl) :: t(3)
        real(kind = dbl) :: tr

        a = 0d0
        r = 0d0
        do k = 1, m
            x = [1d0, u(idx(k)), v(idx(k))]
            do i = 1, 3
                do j = 1, 3
                    a(i, j) = a(i, j) + x(i) * x(j)
                end do
                r(i) = r(i) + x(i) * w(idx(k))
            end do
        end do
        b = 0d0
        ok = .false.
        !> Without velocity bias the intercept row and column drop out.
        if (no_bias) then
            p = 2
        else
            p = 1
        end if
        nd = 3
        !> Gaussian elimination with partial pivoting on rows p..3.
        do i = p, nd
            j = i - 1 + maxloc(abs(a(i:nd, i)), 1)
            if (abs(a(j, i)) <= 1d-12 * max(1d0, maxval(abs(a(p:nd, p:nd))))) return
            if (j /= i) then
                t = a(i, :)
                a(i, :) = a(j, :)
                a(j, :) = t
                tr = r(i)
                r(i) = r(j)
                r(j) = tr
            end if
            do j = i + 1, nd
                f = a(j, i) / a(i, i)
                a(j, i:nd) = a(j, i:nd) - f * a(i, i:nd)
                r(j) = r(j) - f * r(i)
            end do
        end do
        do i = nd, p, -1
            b(i) = (r(i) - sum(a(i, i + 1:nd) * b(i + 1:nd))) / a(i, i)
        end do
        ok = .true.
    end subroutine FitPlane

    !> xorshift64 (Marsaglia 2003): the next state.
    pure function XorShift(s) result(t)
        integer(kind = 8), intent(in) :: s
        integer(kind = 8) :: t

        t = s
        t = ieor(t, ishft(t, 13))
        t = ieor(t, ishft(t, -7))
        t = ieor(t, ishft(t, 17))
    end function XorShift

    !***************************************************************************
    !> \brief Standard deviations, in degrees, of the pitch and roll fitted to
    !>        nboot resamples of the n points, with replacement. err when
    !>        fewer than two resamples could be fitted.
    !***************************************************************************
    pure subroutine BootstrapTilt(u, v, w, n, nboot, seed, no_bias, err, &
        sd_alpha, sd_beta)
        integer, intent(in) :: n
        integer, intent(in) :: nboot
        integer, intent(in) :: seed
        real(kind = dbl), intent(in) :: u(n)
        real(kind = dbl), intent(in) :: v(n)
        real(kind = dbl), intent(in) :: w(n)
        logical, intent(in) :: no_bias
        real(kind = dbl), intent(in) :: err
        real(kind = dbl), intent(out) :: sd_alpha
        real(kind = dbl), intent(out) :: sd_beta
        integer :: idx(n)
        integer(kind = 8) :: s
        integer :: k
        integer :: i
        integer :: nfit
        real(kind = dbl) :: b(3)
        real(kind = dbl) :: al
        real(kind = dbl) :: be
        real(kind = dbl) :: sa
        real(kind = dbl) :: saa
        real(kind = dbl) :: sb
        real(kind = dbl) :: sbb
        logical :: ok

        sd_alpha = err
        sd_beta = err
        if (n < 3 .or. nboot < 2) return
        !> Never zero: xorshift would stay there.
        s = 88172645463325252_8 + int(seed, 8) * 1000003_8
        if (s == 0) s = 88172645463325252_8
        sa = 0d0
        saa = 0d0
        sb = 0d0
        sbb = 0d0
        nfit = 0
        do k = 1, nboot
            do i = 1, n
                s = XorShift(s)
                idx(i) = 1 + int(modulo(s, int(n, 8)))
            end do
            call FitPlane(u, v, w, n, idx, n, no_bias, b, ok)
            if (.not. ok) cycle
            call TiltAnglesDeg(b(2), b(3), al, be)
            nfit = nfit + 1
            sa = sa + al
            saa = saa + al ** 2
            sb = sb + be
            sbb = sbb + be ** 2
        end do
        if (nfit < 2) return
        !> Sample standard deviation, MATLAB std's n - 1.
        sd_alpha = sqrt(max(0d0, (saa - sa ** 2 / dble(nfit)) / dble(nfit - 1)))
        sd_beta = sqrt(max(0d0, (sbb - sb ** 2 / dble(nfit)) / dble(nfit - 1)))
    end subroutine BootstrapTilt

end module m_pfboot
