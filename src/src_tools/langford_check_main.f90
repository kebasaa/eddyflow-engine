!***************************************************************************
! langford_check_main.f90
! -----------------------
! Copyright (C) 2026, ETH Zurich, Jonathan Muller
!
! This file is part of EddyFlow.
!
! EddyFlow (TM) is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! EddyFlow (TM) is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with EddyFlow (TM).  If not, see <http://www.gnu.org/licenses/>.
!
!***************************************************************************
!
! \brief       Known-answer checks of the Langford et al. (2015) far-lag
!              noise in m_covmax_core: FarLagCovarianceNoise,
!              IntegralScaleFromCorr and LimitOfDetection.
!
!              1. RE_RMSE^2 = RE_SD^2 + mean^2, Langford Eq. 9, exactly.
!              2. Two independent white-noise series: the covariance at any
!                 lag scatters with sigma_w sigma_c / sqrt(N), so that is
!                 what both statistics must return, to sampling error.
!              3. A slow oscillation common to both series keeps the far-lag
!                 covariance off zero: RMSE must then exceed SD clearly -
!                 the offset Langford's Fig. 2b shows the SD missing.
!              4. Windows outside the series, or with too few lags, give
!                 the error code.
!              5. An exponential correlation exp(-t/T) integrates to T.
!              6. LoD = alpha RE, and the error code passes through.
!
!              Prints one line per check and stops with code 1 on the first
!              failure. Built by mingw32-make langfordcheck.
! \author      Jonathan Muller
!***************************************************************************
program langford_check
    use m_numeric_kinds
    use m_covmax_core, only: FarLagCovarianceNoise, IntegralScaleFromCorr, &
        LimitOfDetection
    implicit none
    integer, parameter :: n = 36000
    real(kind = dbl), parameter :: err = -9999d0
    real(kind = dbl), parameter :: pi = 3.14159265358979d0
    real(kind = dbl) :: w(n)
    real(kind = dbl) :: c(n)
    real(kind = dbl) :: rmse
    real(kind = dbl) :: sd
    real(kind = dbl) :: expected
    real(kind = dbl) :: corr(0:2000)
    real(kind = dbl) :: its
    integer :: nl_rmse
    integer :: nl_sd
    integer :: k
    integer, allocatable :: seed(:)

    call random_seed(size = k)
    allocate(seed(k))
    seed = 20151215
    call random_seed(put = seed)

    !> 1 and 2: independent white noise, unit variance in w, 3 in c.
    call Gaussian(w)
    call Gaussian(c)
    c = 3d0 * c
    call FarLagCovarianceNoise(w, c, n, 0, 3000, 3600, err, 'rmse', rmse, nl_rmse)
    call FarLagCovarianceNoise(w, c, n, 0, 3000, 3600, err, 'sd', sd, nl_sd)
    call Check(nl_rmse == 2 * 601 .and. nl_sd == nl_rmse, 'both windows, every lag counted')
    !> The mean of the pooled covariances, recovered from the two statistics.
    call Check(abs(rmse**2 - sd**2 - MeanSquared()) <= 1d-12 * rmse**2, &
        'RMSE^2 = SD^2 + mean^2 (Langford Eq. 9)')
    expected = 1d0 * 3d0 / sqrt(dble(n))
    call Check(abs(rmse / expected - 1d0) < 0.25d0, &
        'white noise: RE_RMSE ~ sigma_w sigma_c / sqrt(N)')
    call Check(abs(sd / expected - 1d0) < 0.25d0, &
        'white noise: RE_SD ~ sigma_w sigma_c / sqrt(N)')

    !> 3: a common oscillation of period 4 N rows keeps the far-lag
    !> covariance positive and nearly constant over the windows.
    do k = 1, n
        w(k) = w(k) + sin(2d0 * pi * dble(k) / dble(4 * n))
        c(k) = c(k) + sin(2d0 * pi * dble(k) / dble(4 * n))
    end do
    call FarLagCovarianceNoise(w, c, n, 0, 3000, 3600, err, 'rmse', rmse, nl_rmse)
    call FarLagCovarianceNoise(w, c, n, 0, 3000, 3600, err, 'sd', sd, nl_sd)
    call Check(rmse > 5d0 * sd, 'offset: RE_RMSE counts it, RE_SD does not')

    !> The windows follow the centre: c delayed by 500 rows, measured about
    !> lag 500, sees the same noise as the undelayed pair about lag 0.
    call Gaussian(w)
    call Gaussian(c)
    call FarLagCovarianceNoise(w, c, n, 0, 3000, 3600, err, 'rmse', sd, nl_sd)
    c(501:n) = c(1:n - 500)
    c(1:500) = err
    call FarLagCovarianceNoise(w, c, n, 500, 3000, 3600, err, 'rmse', rmse, nl_rmse)
    call Check(abs(rmse / sd - 1d0) < 0.05d0, 'windows are placed about the centre')

    !> 4: outside the series, and too few lags.
    call FarLagCovarianceNoise(w, c, n, 0, 3000, n, err, 'rmse', rmse, nl_rmse)
    call Check(rmse == err, 'a window past the series end gives the error code')
    call FarLagCovarianceNoise(w, c, n, 0, 3000, 3000, err, 'rmse', rmse, nl_rmse)
    call Check(rmse == err .and. nl_rmse == 2, 'fewer than three lags gives the error code')
    c = err
    call FarLagCovarianceNoise(w, c, n, 0, 3000, 3600, err, 'rmse', rmse, nl_rmse)
    call Check(rmse == err .and. nl_rmse == 0, 'an all-error series gives the error code')

    !> 5: exp(-t/T) with T = 5 s at 20 Hz integrates to T, to the
    !> rectangle rule's dt/2 and the 100 s cut-off.
    do k = 0, 2000
        corr(k) = exp(-dble(k) / 20d0 / 5d0)
    end do
    its = IntegralScaleFromCorr(corr, 2000, 1d0 / 20d0, 'full_integral')
    call Check(abs(its - 5d0) < 0.05d0, 'integral scale of exp(-t/T) is T')
    its = IntegralScaleFromCorr(corr, 2000, 1d0 / 20d0, 'cross_e')
    call Check(its > 0d0 .and. its < 5d0, 'cross_e stops at 1/e')

    !> 6
    call Check(LimitOfDetection(2d0, 1.96d0, err) == 3.92d0, 'LoD = alpha RE')
    call Check(LimitOfDetection(err, 1.96d0, err) == err, 'LoD keeps the error code')

    write(*, '(a)') 'langford_check: all passed'

contains

    subroutine Gaussian(x)
        real(kind = dbl), intent(out) :: x(:)
        real(kind = dbl) :: u1(size(x))
        real(kind = dbl) :: u2(size(x))

        call random_number(u1)
        call random_number(u2)
        u1 = max(u1, 1d-300)
        x = sqrt(-2d0 * log(u1)) * cos(2d0 * pi * u2)
    end subroutine Gaussian

    !> The squared mean of the pooled covariances, from RMSE and SD with the
    !> 'sd' branch's own arithmetic undone - an independent route to it.
    real(kind = dbl) function MeanSquared()
        real(kind = dbl) :: cl(-3600:-3000)
        real(kind = dbl) :: ch(3000:3600)
        real(kind = dbl) :: s
        integer :: j

        call SeriesOf(-3600, -3000, cl)
        call SeriesOf(3000, 3600, ch)
        s = 0d0
        do j = -3600, -3000
            s = s + cl(j)
        end do
        do j = 3000, 3600
            s = s + ch(j)
        end do
        MeanSquared = (s / dble(2 * 601)) ** 2
    end function MeanSquared

    subroutine SeriesOf(lmin, lmax, cv)
        use m_covmax_core, only: CrossCovarianceSeries
        integer, intent(in) :: lmin
        integer, intent(in) :: lmax
        real(kind = dbl), intent(out) :: cv(lmin:lmax)

        call CrossCovarianceSeries(w, c, n, lmin, lmax, err, cv)
    end subroutine SeriesOf

    subroutine Check(ok, what)
        logical, intent(in) :: ok
        character(*), intent(in) :: what

        if (ok) then
            write(*, '(a, a)') 'ok    ', what
        else
            write(*, '(a, a)') 'FAIL  ', what
            write(*, '(a, 2es16.8)') '      rmse, sd = ', rmse, sd
            stop 1
        end if
    end subroutine Check
end program langford_check
