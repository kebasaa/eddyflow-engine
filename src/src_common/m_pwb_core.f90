!***************************************************************************
! m_pwb_core.f90
! --------------
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
! \brief       The arithmetic of pre-whitening block-bootstrap time-lag
!              detection, with no dependency on engine state.
!
! \details     Everything here is a function of its arguments: no project
!              settings, no column table, no logging. That is the point. The
!              deterministic half of this chain is pinned to RFlux v3.2.0 by
!              static_checks/test_pwb_reference_static.py, which drives it
!              through pwb_reference_main.f90 - a program that could not link
!              at all if these routines reached for PWBSetup or E2Col.
!
!              The reference is dyco (dyco/pwb.py), whose own test suite pins
!              the same chain to the same R output at twelve significant
!              digits and whose module docstring catalogues every deliberate
!              departure from it. Where a comment here says "R does X", dyco
!              is where that was read.
!
! \author      Jonathan Muller
! \note
! \sa
! \bug
! \deprecated
! \test        static_checks/test_pwb_reference_static.py
! \todo
!***************************************************************************
module m_pwb_core
    use m_numeric_kinds
    implicit none
    private
    public :: LinearDetrend, IsStationary, DifferenceSeries
    public :: FitArAic, ApplyArFilter
    public :: ComputeCcfWindow, ComputeCcovWindow, SmoothAndFill
    public :: ArgmaxAbs, Hdi95, MapLagEstimate, MedianOf
    public :: CopyBlock, MixIn, RandBelow
    public :: PwbPreWhiten, PwbPreWhitenType, BartlettCv99

    !> Everything the deterministic half of the chain produces, which is also
    !> everything the reference test compares.
    type :: PwbPreWhitenType
        logical :: differenced = .false.
        !> AR orders and first coefficients, in the order scalar, w, tsonic.
        integer :: p_scalar = 0, p_w = 0, p_t = 0
        real(kind = dbl) :: phi1_scalar = 0d0, phi1_w = 0d0, phi1_t = 0d0
        !> Peak of the unsmoothed pre-whitened CCF (scalar AR, scalar x w),
        !> in records, and the CCF value there. R: tl_pww, cor_pww.
        integer :: tlag_pw_rl = 0
        real(kind = dbl) :: corr_pw = 0d0
        !> Peak of the raw cross-covariance and its value. R: mcw, cov_mcw.
        integer :: ccov_rl = 0
        real(kind = dbl) :: cov_mcw = 0d0
        !> Records actually used, after the differencing branch.
        integer :: n_eff = 0
    end type PwbPreWhitenType

contains

!***************************************************************************
!> Remove a least-squares straight line, in place.
!>
!> R detrends linearly before the raw cross-covariance. The engine's own
!> CalculateTrend/Detrend pair carries a missing-value sentinel through, which
!> these arrays no longer need - they are gap-filled by the time they arrive -
!> and dividing that work in two just to reuse it would cost more than the
!> eight lines it saves.
!***************************************************************************
subroutine LinearDetrend(x, n)
    integer, intent(in) :: n
    real(kind = dbl), intent(inout) :: x(n)
    integer :: i
    real(kind = dbl) :: sx, st, stt, sxt, tbar, xbar, slope, tt

    if (n < 2) return
    sx = 0d0; st = 0d0; stt = 0d0; sxt = 0d0
    do i = 1, n
        tt = dble(i - 1)
        sx = sx + x(i)
        st = st + tt
        stt = stt + tt * tt
        sxt = sxt + x(i) * tt
    end do
    xbar = sx / dble(n)
    tbar = st / dble(n)
    if (stt - st * tbar == 0d0) then
        x = x - xbar
        return
    end if
    slope = (sxt - st * xbar) / (stt - st * tbar)
    do i = 1, n
        x(i) = x(i) - (xbar + slope * (dble(i - 1) - tbar))
    end do
end subroutine LinearDetrend

!***************************************************************************
!> Breitung (2002) variance-ratio unit-root test: stationary or not.
!>
!> R: egcm::bvr.test at alpha = 0.01. Because 0.01 is itself a tabulated
!> quantile, comparing the p-value to it is the same as comparing the
!> statistic to the 1% critical value; for any EC averaging period egcm
!> clamps the sample size to its n=1250 column, whose entry is the constant
!> below. Turbulent series give rho of order 1e-5 and pass; a drifting one
!> gives 0.05 to 0.1 and fails.
!>
!> Failing on any one of the three series differences all three. That is not
!> the edge case it looks like - sonic temperature following the diurnal
!> cycle over half an hour is enough, and dyco measures p = 0.057 on its
!> bundled half hour.
!***************************************************************************
logical function IsStationary(x, n)
    integer, intent(in) :: n
    real(kind = dbl), intent(in) :: x(n)
    integer :: i
    real(kind = dbl) :: meanx, sse, cum, rho
    real(kind = dbl), parameter :: cv_1pct = 0.00537748023783321d0

    meanx = sum(x) / dble(n)
    sse = 0d0
    cum = 0d0
    rho = 0d0
    do i = 1, n
        sse = sse + (x(i) - meanx)**2
    end do
    if (sse <= 0d0) then
        IsStationary = .true.
        return
    end if
    do i = 1, n
        cum = cum + x(i) - meanx
        rho = rho + cum**2
    end do
    rho = rho / (dble(n)**2 * sse)
    IsStationary = rho < cv_1pct
end function IsStationary

!> First differences, compacted to the front; ne is n-1 on return.
!>
!> R and numpy both return n-1 values here. Writing them back in place over
!> 2..n and setting x(1) = 0 keeps the length at n and feeds the AR fit one
!> fabricated sample: at n = 6000 that moves the first AR coefficient in the
!> sixth significant digit against RFlux, which the reference test sees.
subroutine DifferenceSeries(x, n, ne)
    integer, intent(in) :: n
    integer, intent(out) :: ne
    real(kind = dbl), intent(inout) :: x(n)
    integer :: i
    do i = 1, n - 1
        x(i) = x(i+1) - x(i)
    end do
    ne = n - 1
    x(n) = 0d0
end subroutine DifferenceSeries

!***************************************************************************
!> AR(p) by AIC, order searched up to floor(100*log10(N)).
!>
!> That ceiling is R's, and real chunks use a good deal of it: dyco reports
!> orders of 133, 87 and 312 on a 20 Hz half hour. An order below the optimum
!> leaves autocorrelation in the residuals, which is exactly what broadens
!> the CCF peak that pre-whitening exists to sharpen - which is why the cap
!> this once offered as a speed option is gone.
!***************************************************************************
subroutine FitArAic(x, n, phi_best, p_best)
    integer, intent(in) :: n
    real(kind = dbl), intent(in) :: x(n)
    real(kind = dbl), allocatable, intent(out) :: phi_best(:)
    integer, intent(out) :: p_best
    integer :: max_lag, p, i
    real(kind = dbl) :: meanx, sigma2, kappa, best_aic, aic
    real(kind = dbl), allocatable :: acf(:), phi(:), phi_old(:)

    max_lag = min(int(floor(100d0 * log10(dble(max(n, 2))))), n - 1)
    if (max_lag < 1) then
        allocate(phi_best(0))
        p_best = 0
        return
    end if
    allocate(acf(0:max_lag))
    meanx = sum(x) / dble(n)
    do p = 0, max_lag
        acf(p) = 0d0
        do i = 1, n - p
            acf(p) = acf(p) + (x(i) - meanx) * (x(i+p) - meanx)
        end do
        acf(p) = acf(p) / dble(n)
    end do
    if (acf(0) <= 0d0) then
        allocate(phi_best(0))
        p_best = 0
        deallocate(acf)
        return
    end if

    allocate(phi(1:max_lag), phi_old(1:max_lag), phi_best(0))
    p_best = 0
    best_aic = dble(n) * log(acf(0))
    sigma2 = acf(0)
    phi = 0d0
    do p = 1, max_lag
        if (p == 1) then
            kappa = acf(1) / sigma2
        else
            kappa = (acf(p) - dot_product(phi(1:p-1), acf(p-1:1:-1))) / sigma2
        end if
        phi_old = phi
        if (p > 1) phi(1:p-1) = phi_old(1:p-1) - kappa * phi_old(p-1:1:-1)
        phi(p) = kappa
        sigma2 = sigma2 * (1d0 - kappa**2)
        if (sigma2 <= 0d0) exit
        aic = dble(n) * log(sigma2) + 2d0 * dble(p)
        if (aic < best_aic) then
            best_aic = aic
            p_best = p
            if (allocated(phi_best)) deallocate(phi_best)
            allocate(phi_best(p_best))
            phi_best = phi(1:p_best)
        end if
    end do
    deallocate(acf, phi, phi_old)
end subroutine FitArAic

!> R: stats::filter(x, c(1, -phi), method="convolution", sides=1). The first
!> p outputs have no past to filter against; R marks them NA and dyco zeroes
!> them before resampling, which is what this does.
subroutine ApplyArFilter(x, n, phi, p, y)
    integer, intent(in) :: n, p
    real(kind = dbl), intent(in) :: x(n)
    real(kind = dbl), intent(in) :: phi(:)
    real(kind = dbl), intent(out) :: y(n)
    integer :: i, j
    real(kind = dbl) :: meanx

    meanx = sum(x) / dble(n)
    y = x - meanx
    if (p <= 0) return
    do i = n, 1, -1
        y(i) = x(i) - meanx
        do j = 1, min(p, i - 1)
            y(i) = y(i) - phi(j) * (x(i-j) - meanx)
        end do
        if (i <= p) y(i) = 0d0
    end do
end subroutine ApplyArFilter

!***************************************************************************
!> Normalised cross-correlation over [min_rl, max_rl].
!>
!> The divisor is constant across lags, so within one pairing it cannot move
!> an argmax. It is not optional all the same: the four pre-whitening
!> combinations are compared against each other at their modes, and they pair
!> the scalar against vertical wind and against sonic temperature, whose
!> covariances carry different physical units. Unnormalised, the winner is
!> decided by that unit scale.
!>
!> xc and yc are caller-owned scratch of length n, so a bootstrap loop does
!> not allocate once per replicate.
!***************************************************************************
subroutine ComputeCcfWindow(x, y, n, min_rl, max_rl, ccf, xc, yc)
    integer, intent(in) :: n, min_rl, max_rl
    real(kind = dbl), intent(in) :: x(n), y(n)
    real(kind = dbl), intent(out) :: ccf(min_rl:max_rl)
    real(kind = dbl), intent(inout) :: xc(n), yc(n)
    integer :: lag, i, nn
    real(kind = dbl) :: mx, my, vx, vy, denom, cov

    mx = sum(x) / dble(n)
    my = sum(y) / dble(n)
    xc = x - mx
    yc = y - my
    vx = sum(xc * xc)
    vy = sum(yc * yc)
    denom = sqrt(vx * vy)
    if (denom <= 0d0) then
        ccf = 0d0
        return
    end if

    do lag = min_rl, max_rl
        nn = n - abs(lag)
        if (nn <= 1) then
            ccf(lag) = 0d0
            cycle
        end if
        cov = 0d0
        if (lag >= 0) then
            do i = 1, nn
                cov = cov + xc(i) * yc(i + lag)
            end do
        else
            do i = 1, nn
                cov = cov + xc(i - lag) * yc(i)
            end do
        end if
        ccf(lag) = cov / denom
    end do
end subroutine ComputeCcfWindow

!***************************************************************************
!> Un-normalised cross-covariance, the biased estimator R's
!> ccf(type="covariance") returns.
!>
!> Biased means two specific things, and getting either wrong shows up
!> immediately against the reference. The series are centred ONCE on their
!> own full-length means, not per lag on the means of whatever overlaps at
!> that lag. And every lag divides by n, not by its overlap count: dividing
!> by the overlap inflates the covariance by n/(n-|lag|), which at lag 169 of
!> 6000 is 2.9% - the exact error this had.
!>
!> xc and yc are caller-owned scratch of length n. `missing` is written where
!> a lag has too little overlap to mean anything.
!***************************************************************************
subroutine ComputeCcovWindow(x, y, n, min_rl, max_rl, ccov, missing, xc, yc)
    integer, intent(in) :: n, min_rl, max_rl
    real(kind = dbl), intent(in) :: x(n), y(n)
    real(kind = dbl), intent(in) :: missing
    real(kind = dbl), intent(out) :: ccov(min_rl:max_rl)
    real(kind = dbl), intent(inout) :: xc(n), yc(n)
    integer :: lag, i, nn
    real(kind = dbl) :: cov

    xc(1:n) = x - sum(x) / dble(n)
    yc(1:n) = y - sum(y) / dble(n)
    do lag = min_rl, max_rl
        nn = n - abs(lag)
        if (nn <= 1) then
            ccov(lag) = missing
            cycle
        end if
        cov = 0d0
        if (lag >= 0) then
            do i = 1, nn
                cov = cov + xc(i) * yc(i + lag)
            end do
        else
            do i = 1, nn
                cov = cov + xc(i - lag) * yc(i)
            end do
        end if
        ccov(lag) = cov / dble(n)
    end do
end subroutine ComputeCcovWindow

!***************************************************************************
!> Centred rolling mean, edges carried in from the nearest computed value.
!>
!> The fill is R's two-pass zoo::na.locf, and it exists so that a peak at the
!> boundary of the search window is still a candidate for the argmax rather
!> than being skipped as missing.
!>
!> The divisor is the number of terms actually summed, which for an even
!> width is width+1 - a centred window has no midpoint otherwise.
!***************************************************************************
subroutine SmoothAndFill(x, min_rl, max_rl, width, y)
    integer, intent(in) :: min_rl, max_rl, width
    real(kind = dbl), intent(in) :: x(min_rl:max_rl)
    real(kind = dbl), intent(out) :: y(min_rl:max_rl)
    integer :: i, j, half, first_valid, last_valid
    half = width / 2
    first_valid = min_rl + half
    last_valid = max_rl - half
    do i = first_valid, last_valid
        y(i) = 0d0
        do j = i - half, i + half
            y(i) = y(i) + x(j)
        end do
        y(i) = y(i) / dble(2 * half + 1)
    end do
    if (first_valid <= last_valid) then
        y(min_rl:first_valid - 1) = y(first_valid)
        y(last_valid + 1:max_rl) = y(last_valid)
    else
        y = x
    end if
end subroutine SmoothAndFill

integer function ArgmaxAbs(x, min_rl, max_rl)
    integer, intent(in) :: min_rl, max_rl
    real(kind = dbl), intent(in) :: x(min_rl:max_rl)
    integer :: i
    ArgmaxAbs = min_rl
    do i = min_rl, max_rl
        if (abs(x(i)) > abs(x(ArgmaxAbs))) ArgmaxAbs = i
    end do
end function ArgmaxAbs

!> Shortest interval containing 95% of the samples. Sorts in place.
subroutine Hdi95(x, n, lo, hi)
    integer, intent(in) :: n
    real(kind = dbl), intent(inout) :: x(n)
    real(kind = dbl), intent(out) :: lo, hi
    integer :: i, j, m, best
    real(kind = dbl) :: tmp, width
    do i = 2, n
        tmp = x(i)
        j = i - 1
        do while (j >= 1)
            if (x(j) <= tmp) exit
            x(j+1) = x(j)
            j = j - 1
        end do
        x(j+1) = tmp
    end do
    m = max(1, int(floor(0.95d0 * dble(n))))
    if (m >= n) then
        lo = x(1); hi = x(n); return
    end if
    best = 1
    width = x(1+m) - x(1)
    do i = 2, n - m
        if (x(i+m) - x(i) < width) then
            best = i
            width = x(i+m) - x(i)
        end if
    end do
    lo = x(best)
    hi = x(best + m)
end subroutine Hdi95

!> Mode of the bootstrap lag distribution, by Gaussian KDE on the integer
!> grid. R uses bayestestR::map_estimate on jittered samples; dyco uses
!> scipy's gaussian_kde. All three differ in bandwidth and grid and agree
!> well inside bootstrap noise.
integer function MapLagEstimate(samples, n)
    integer, intent(in) :: n
    integer, intent(in) :: samples(n)
    integer :: i, grid, lo, hi, best
    real(kind = dbl) :: mean_s, var_s, sd_s, bw, dens, best_dens, z

    lo = minval(samples)
    hi = maxval(samples)
    if (lo == hi) then
        MapLagEstimate = lo
        return
    end if

    mean_s = sum(dble(samples)) / dble(n)
    var_s = 0d0
    do i = 1, n
        var_s = var_s + (dble(samples(i)) - mean_s)**2
    end do
    sd_s = sqrt(max(0d0, var_s / max(1d0, dble(n - 1))))
    bw = max(1d0, 1.06d0 * sd_s * dble(n)**(-0.2d0))

    best = lo
    best_dens = -1d0
    do grid = lo, hi
        dens = 0d0
        do i = 1, n
            z = (dble(grid) - dble(samples(i))) / bw
            dens = dens + exp(-0.5d0 * z * z)
        end do
        if (dens > best_dens) then
            best_dens = dens
            best = grid
        end if
    end do
    MapLagEstimate = best
end function MapLagEstimate

real(kind = dbl) function MedianOf(x, n)
    integer, intent(in) :: n
    real(kind = dbl), intent(inout) :: x(n)
    integer :: i, j
    real(kind = dbl) :: tmp

    do i = 2, n
        tmp = x(i)
        j = i - 1
        do while (j >= 1)
            if (x(j) <= tmp) exit
            x(j+1) = x(j)
            j = j - 1
        end do
        x(j+1) = tmp
    end do
    if (mod(n, 2) == 1) then
        MedianOf = x((n + 1) / 2)
    else
        MedianOf = 0.5d0 * (x(n/2) + x(n/2 + 1))
    end if
end function MedianOf

subroutine CopyBlock(x, y, n, start, block_len, xb, yb, pos)
    integer, intent(in) :: n, start, block_len
    integer, intent(inout) :: pos
    real(kind = dbl), intent(in) :: x(n), y(n)
    real(kind = dbl), intent(inout) :: xb(n), yb(n)
    integer :: j, src
    do j = 0, block_len - 1
        if (pos > n) exit
        src = min(n, start + j)
        xb(pos) = x(src)
        yb(pos) = y(src)
        pos = pos + 1
    end do
end subroutine CopyBlock

!> One xorshift round over an absorbed value. Shifts and exclusive-or only,
!> so there is no multiplication to overflow.
subroutine MixIn(state, val)
    integer(8), intent(inout) :: state
    integer(8), intent(in) :: val

    state = ieor(state, val)
    state = ieor(state, ishft(state, 13))
    state = ieor(state, ishft(state, -7))
    state = ieor(state, ishft(state, 17))
end subroutine MixIn

!***************************************************************************
!> A uniform integer in [0, upper), from the top bits of xorshift64.
!>
!> The low bits of a linear congruential generator are its worst, and a
!> bootstrap replicate draws about ninety block starts in a row: correlated
!> starts make replicates resemble one another, which narrows the HDI and
!> hands out S1 to detections that have not earned it. ISHFT is a logical
!> shift, so the shifted result is non-negative whatever the sign of state.
!***************************************************************************
integer function RandBelow(state, upper)
    integer(8), intent(inout) :: state
    integer, intent(in) :: upper

    state = ieor(state, ishft(state, 13))
    state = ieor(state, ishft(state, -7))
    state = ieor(state, ishft(state, 17))
    RandBelow = int(mod(ishft(state, -33), int(max(1, upper), 8)), 4)
end function RandBelow

!> R's Bartlett 99% significance band for the pre-whitened CCF, 3.291 over
!> the root of n*13. The 13 is R's smoothing factor; ported as R states it.
real(kind = dbl) function BartlettCv99(n_eff)
    integer, intent(in) :: n_eff
    BartlettCv99 = 3.291d0 / sqrt(dble(max(1, n_eff)) * 13d0)
end function BartlettCv99

!***************************************************************************
!> The deterministic half of the PWB chain: stationarity, AR fits, filtered
!> series, the raw cross-covariance and the full-data pre-whitened CCF.
!>
!> Everything the block bootstrap needs comes out of here, and so does
!> everything the reference test compares. The three inputs arrive gap-filled
!> and are consumed, not preserved.
!>
!> The differencing branch feeds the AR filters ONLY. The raw cross-covariance
!> is R's ccf(detrend(scalar), detrend(w), type="covariance") over the
!> UNDIFFERENCED series: read off the differenced arrays instead it becomes a
!> covariance of increments, which on a drifting record is smaller by two
!> orders of magnitude and free to carry the opposite sign. dyco had that bug
!> and fixed it; so did this.
!***************************************************************************
subroutine PwbPreWhiten(ss, ww, tt, n, min_rl, max_rl, missing, res, &
    s_fs, w_fs, t_fs, s_fw, w_fw, s_ft, t_ft, ccov, xc, yc)
    integer, intent(in) :: n, min_rl, max_rl
    real(kind = dbl), intent(inout) :: ss(n), ww(n), tt(n)
    real(kind = dbl), intent(in) :: missing
    type(PwbPreWhitenType), intent(out) :: res
    real(kind = dbl), intent(out) :: s_fs(n), w_fs(n), t_fs(n)
    real(kind = dbl), intent(out) :: s_fw(n), w_fw(n), s_ft(n), t_ft(n)
    real(kind = dbl), intent(out) :: ccov(min_rl:max_rl)
    real(kind = dbl), intent(inout) :: xc(n), yc(n)
    real(kind = dbl), allocatable :: phi_s(:), phi_w(:), phi_t(:)
    real(kind = dbl), allocatable :: ss_undiff(:), ww_undiff(:), pw_ccf(:)
    integer :: p_s, p_w, p_t, ne, n1, n2, n3, ntrim

    !> Kept before anything differences them.
    allocate(ss_undiff(n), ww_undiff(n))
    ss_undiff = ss
    ww_undiff = ww

    ne = n
    res%differenced = .not. (IsStationary(ss, n) .and. IsStationary(ww, n) &
        .and. IsStationary(tt, n))
    if (res%differenced) then
        call DifferenceSeries(ss, n, n1)
        call DifferenceSeries(ww, n, n2)
        call DifferenceSeries(tt, n, n3)
        ne = min(n1, n2, n3)
    end if
    res%n_eff = ne

    call FitArAic(ss(1:ne), ne, phi_s, p_s)
    call FitArAic(ww(1:ne), ne, phi_w, p_w)
    call FitArAic(tt(1:ne), ne, phi_t, p_t)
    res%p_scalar = p_s
    res%p_w = p_w
    res%p_t = p_t
    if (p_s > 0) res%phi1_scalar = phi_s(1)
    if (p_w > 0) res%phi1_w = phi_w(1)
    if (p_t > 0) res%phi1_t = phi_t(1)

    call ApplyArFilter(ss(1:ne), ne, phi_s, p_s, s_fs(1:ne))
    call ApplyArFilter(ww(1:ne), ne, phi_s, p_s, w_fs(1:ne))
    call ApplyArFilter(tt(1:ne), ne, phi_s, p_s, t_fs(1:ne))
    call ApplyArFilter(ss(1:ne), ne, phi_w, p_w, s_fw(1:ne))
    call ApplyArFilter(ww(1:ne), ne, phi_w, p_w, w_fw(1:ne))
    call ApplyArFilter(ss(1:ne), ne, phi_t, p_t, s_ft(1:ne))
    call ApplyArFilter(tt(1:ne), ne, phi_t, p_t, t_ft(1:ne))

    !> The raw cross-covariance, on the UNDIFFERENCED series at full length.
    !> R: ccf(detrend(scalar), detrend(w), type="covariance").
    call LinearDetrend(ww_undiff, n)
    call LinearDetrend(ss_undiff, n)
    call ComputeCcovWindow(ww_undiff, ss_undiff, n, min_rl, max_rl, ccov, missing, xc, yc)
    res%ccov_rl = ArgmaxAbs(ccov, min_rl, max_rl)
    res%cov_mcw = ccov(res%ccov_rl)

    !> The full-data pre-whitened CCF, scalar AR, scalar against w.
    !> R: ccf_pww, and tl_pww is the argmax of the UNSMOOTHED one.
    !>
    !> The first p_s filtered values have no past to filter against - R marks
    !> them NA and drops them here with na.action = na.omit. Keeping them as
    !> zeros instead shifts the means, the variances and every sum that
    !> follows, which at an AR order of 67 moves cor_pww in the fourth
    !> significant digit. The BOOTSTRAP keeps the zeros, as R does; only this
    !> full-data diagnostic trims.
    ntrim = max(1, ne - p_s)
    allocate(pw_ccf(min_rl:max_rl))
    call ComputeCcfWindow(w_fs(p_s+1:ne), s_fs(p_s+1:ne), ntrim, &
        min_rl, max_rl, pw_ccf, xc(1:ntrim), yc(1:ntrim))
    res%tlag_pw_rl = ArgmaxAbs(pw_ccf, min_rl, max_rl)
    res%corr_pw = pw_ccf(res%tlag_pw_rl)

    deallocate(ss_undiff, ww_undiff, pw_ccf)
    if (allocated(phi_s)) deallocate(phi_s)
    if (allocated(phi_w)) deallocate(phi_w)
    if (allocated(phi_t)) deallocate(phi_t)
end subroutine PwbPreWhiten

end module m_pwb_core
