!***************************************************************************
! m_covmax_core.f90
! -----------------
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
! \brief       The lagged cross-covariance covariance maximisation searches,
!              computed without building each shifted pair.
!
!              CovMax is what time-lag optimisation runs for every gas of every
!              period, and what PWB runs for its terminal fallback. For each lag
!              it used to allocate two arrays, copy the shifted pair into one
!              and that into the other, and hand the pair to
!              CovarianceMatrixNoError - which computes the full 2x2 matrix,
!              four passes over the data, of which only Cov(1,2) was read.
!
!              This computes Cov(1,2) alone, reading the shifted elements where
!              they already are. It performs exactly the operations the matrix
!              routine performs for that entry: the same pairwise test against
!              the error code, the same products, sums and count accumulated in
!              the same order, and the same C/N - (si/N)(sj/N) at the end. So
!              every value is bitwise the one it replaces - which
!              src_tools/covmax_bitexact_main.f90 checks against a frozen copy
!              of the old path, sparse columns included.
!
! \note        Pure, and a function of its arguments only, so it can be checked
!              outside the engine like m_pwb_core.
! \author      Jonathan Muller
! \sa          timelag_handle.f90, stats_operator_no_error.f90
!***************************************************************************
module m_covmax_core
    use m_numeric_kinds
    implicit none
    private

    public :: CrossCovarianceSeries
    public :: FarLagCovarianceNoise
    public :: IntegralScaleFromCorr
    public :: LimitOfDetection

contains

    !***************************************************************************
    !> \brief Cov(1,2) of the pair (col1, col2) shifted by every lag in
    !>        [lagmin, lagmax], counting only pairs where neither is err.
    !>
    !> A negative lag pairs col1(k - lag) with col2(k); a non-negative one
    !> pairs col1(k) with col2(k + lag); k runs from 1 to n - |lag|. A lag with
    !> no valid pair - including one as wide as the series - gives err.
    !***************************************************************************
    pure subroutine CrossCovarianceSeries(col1, col2, n, lagmin, lagmax, err, cov)
        integer, intent(in) :: n
        integer, intent(in) :: lagmin
        integer, intent(in) :: lagmax
        real(kind = dbl), intent(in) :: col1(n)
        real(kind = dbl), intent(in) :: col2(n)
        real(kind = dbl), intent(in) :: err
        real(kind = dbl), intent(out) :: cov(lagmin:lagmax)
        integer :: lag
        integer :: k
        integer :: n2
        integer :: nact
        real(kind = dbl) :: a
        real(kind = dbl) :: b
        real(kind = dbl) :: c
        real(kind = dbl) :: si
        real(kind = dbl) :: sj

        do lag = lagmin, lagmax
            n2 = n - abs(lag)
            c = 0d0
            si = 0d0
            sj = 0d0
            nact = 0
            if (lag < 0) then
                do k = 1, n2
                    a = col1(k - lag)
                    b = col2(k)
                    if (a /= err .and. b /= err) then
                        nact = nact + 1
                        c = c + a * b
                        si = si + a
                        sj = sj + b
                    end if
                end do
            else
                do k = 1, n2
                    a = col1(k)
                    b = col2(k + lag)
                    if (a /= err .and. b /= err) then
                        nact = nact + 1
                        c = c + a * b
                        si = si + a
                        sj = sj + b
                    end if
                end do
            end if
            !> The matrix routine's normalisation, step for step.
            if (nact /= 0) then
                si = si / dble(nact)
                sj = sj / dble(nact)
                c = c / dble(nact)
                c = c - si * sj
            else
                c = err
            end if
            cov(lag) = c
        end do
    end subroutine CrossCovarianceSeries

    !***************************************************************************
    !> rief Random flux error from the cross-covariance far from the lag,
    !>        Langford et al. (2015), Sect. 2.2.
    !>
    !> The cross-covariance of (col1, col2) is taken over the two windows
    !> [centre - hi, centre - lo] and [centre + lo, centre + hi], in rows and
    !> in CrossCovarianceSeries' sign convention, and the two are pooled.
    !> Out there w and the scalar should no longer be correlated, so what the
    !> function does is noise. With stat 'rmse' this returns
    !> sqrt(mean f^2) - Langford's RE_RMSE, Eq. 9, which also counts an
    !> offset of the function from zero - and with 'sd' the population
    !> standard deviation, so that RE_RMSE^2 = RE_SD^2 + mean^2 exactly.
    !>
    !> Lags with no valid pair are skipped. Fewer than three valid lags, or a
    !> window that does not fit inside the series, gives err; nlags reports
    !> how many were used.
    !***************************************************************************
    pure subroutine FarLagCovarianceNoise(col1, col2, n, centre, lo, hi, err, &
        stat, re, nlags)
        integer, intent(in) :: n
        integer, intent(in) :: centre
        integer, intent(in) :: lo
        integer, intent(in) :: hi
        real(kind = dbl), intent(in) :: col1(n)
        real(kind = dbl), intent(in) :: col2(n)
        real(kind = dbl), intent(in) :: err
        character(*), intent(in) :: stat
        real(kind = dbl), intent(out) :: re
        integer, intent(out) :: nlags
        real(kind = dbl) :: cov_lo(centre - hi:centre - lo)
        real(kind = dbl) :: cov_hi(centre + lo:centre + hi)
        real(kind = dbl) :: s1
        real(kind = dbl) :: s2
        real(kind = dbl) :: mu
        integer :: lag

        re = err
        nlags = 0
        if (lo < 0 .or. hi < lo) return
        if (abs(centre - hi) >= n .or. abs(centre + hi) >= n) return

        call CrossCovarianceSeries(col1, col2, n, centre - hi, centre - lo, err, cov_lo)
        call CrossCovarianceSeries(col1, col2, n, centre + lo, centre + hi, err, cov_hi)

        s1 = 0d0
        s2 = 0d0
        do lag = centre - hi, centre - lo
            if (cov_lo(lag) == err) cycle
            nlags = nlags + 1
            s1 = s1 + cov_lo(lag)
            s2 = s2 + cov_lo(lag) ** 2
        end do
        do lag = centre + lo, centre + hi
            if (cov_hi(lag) == err) cycle
            nlags = nlags + 1
            s1 = s1 + cov_hi(lag)
            s2 = s2 + cov_hi(lag) ** 2
        end do
        if (nlags < 3) return

        if (stat == 'sd') then
            mu = s1 / dble(nlags)
            re = sqrt(max(0d0, s2 / dble(nlags) - mu ** 2))
        else
            re = sqrt(s2 / dble(nlags))
        end if
    end subroutine FarLagCovarianceNoise

    !***************************************************************************
    !> rief Integral scale of a cross-correlation already normalised to one
    !>        at lag zero, by the rule its_meth names.
    !>
    !> The three rules IntegralTurbulenceScale has always applied, moved here
    !> so the automatic Langford window measures the same scale: integrate
    !> |corr| up to its first crossing of zero ('cross_0'), of 1/e
    !> ('cross_e'), or over every lag given ('full_integral'). Same
    !> operations in the same order as before, so the scale it returns there
    !> is bitwise unchanged.
    !***************************************************************************
    pure function IntegralScaleFromCorr(corr, lagmax, dt, its_meth) result(its)
        integer, intent(in) :: lagmax
        real(kind = dbl), intent(in) :: corr(0:lagmax)
        real(kind = dbl), intent(in) :: dt
        character(*), intent(in) :: its_meth
        real(kind = dbl) :: its
        integer :: lag

        its = 0d0
        select case(its_meth)
            case('cross_0')
                do lag = 0, lagmax
                    if (corr(lag) < 0d0) exit
                    its = its + dt * dabs(corr(lag))
                end do
            case('cross_e')
                do lag = 0, lagmax
                    if (corr(lag) <  1d0/exp(1d0)) exit
                    its = its + dt * dabs(corr(lag))
                end do
            case('full_integral')
                do lag = 0, lagmax
                    its = its + dt * dabs(corr(lag))
                end do
        end select
    end function IntegralScaleFromCorr

    !***************************************************************************
    !> rief Limit of detection LoD = alpha * RE, Langford et al. (2015),
    !>        keeping the error code where there is no random error.
    !***************************************************************************
    pure elemental function LimitOfDetection(re, alpha, err) result(lod)
        real(kind = dbl), intent(in) :: re
        real(kind = dbl), intent(in) :: alpha
        real(kind = dbl), intent(in) :: err
        real(kind = dbl) :: lod

        if (re == err) then
            lod = err
        else
            lod = alpha * re
        end if
    end function LimitOfDetection

end module m_covmax_core
