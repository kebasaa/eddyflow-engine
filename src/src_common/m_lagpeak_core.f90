!***************************************************************************
! m_lagpeak_core.f90
! ------------------
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
! \brief       Peak search on a cross-covariance function after smoothing,
!              as innFLUX does it (Striednig et al. 2020, Sect. 2.3.2; code
!              subroutines/find_covariance_peak.m, innFLUX v1.2.0, University
!              of Innsbruck). Smoothing the function before the extremum is
!              taken follows Taipale et al. (2010).
!
!              innFLUX smooths with MATLAB filtfilt and a normalised Hamming
!              window. filtfilt of an FIR filter b is, away from the two ends
!              of the series, exactly a convolution with the zero-phase kernel
!              h(m) = sum_{j-k=m} b(j) b(k): the forward pass convolves with b,
!              the backward pass with b reversed. The ends are where filtfilt's
!              reflection padding and initial conditions act. Rather than
!              reproduce them, every caller computes the covariance over a
!              window wider than the search by PeakPadding lags on each side,
!              so the search window only ever sees interior values - which
!              are the same as innFLUX's whatever the padding.
!
!              Lags are EddyFlow's: a positive lag means the second series
!              (the scalar) is delayed against the first (w), as in
!              m_covmax_core's CrossCovarianceSeries. MATLAB's xcov(w, c)
!              runs the other way, m = -lag; innFLUX's
!              lag_samples = -(idx - range - 1) undoes that, so its lag in
!              samples is this lag.
!
!              Everything here is off unless asked for: with smooth_len = 0
!              and sign_mode = 0 FindCovPeak is not called at all and CovMax
!              keeps its own search.
!
! \note        Pure, and a function of its arguments only, so it can be checked
!              outside the engine (src_tools/innflux_compare_main.f90).
! \author      Jonathan Muller
! \sa          m_covmax_core.f90, timelag_handle.f90
!***************************************************************************
module m_lagpeak_core
    use m_numeric_kinds
    implicit none
    private

    public :: HammingKernel
    public :: ZeroPhaseKernel
    public :: ZeroPhaseSmooth
    public :: PeakPadding
    public :: FindCovPeak
    public :: FokenRelDev

    !> Length of the moving average innFLUX runs over the derivative of the
    !> covariance function to decide whether its peak is a maximum or a
    !> minimum. Fixed in find_covariance_peak.m, so fixed here.
    integer, parameter, public :: DerivFilterLen = 100

    !> Sign rules. 0: the largest |cov| (EddyFlow's covariance maximisation).
    !> 1: the largest cov - innFLUX's positive_peak_only, used on summed
    !> absolute functions. 2: innFLUX's derivative rule.
    integer, parameter, public :: PeakAbs = 0
    integer, parameter, public :: PeakPositive = 1
    integer, parameter, public :: PeakDerivative = 2

contains

    !***************************************************************************
    !> \brief MATLAB hamming(L), the symmetric window, normalised to sum one,
    !>        which is how innFLUX applies it: hamming(L)/sum(hamming(L)).
    !***************************************************************************
    pure subroutine HammingKernel(l, b)
        integer, intent(in) :: l
        real(kind = dbl), intent(out) :: b(l)
        real(kind = dbl), parameter :: pi = 3.14159265358979323846d0
        integer :: k

        if (l == 1) then
            b(1) = 1d0
            return
        end if
        do k = 0, l - 1
            b(k + 1) = 0.54d0 - 0.46d0 * cos(2d0 * pi * dble(k) / dble(l - 1))
        end do
        b = b / sum(b)
    end subroutine HammingKernel

    !***************************************************************************
    !> \brief The kernel filtfilt(b, 1, .) amounts to away from the ends:
    !>        h(m) = sum over j - k = m of b(j) b(k), m = -(l-1) .. l-1.
    !***************************************************************************
    pure subroutine ZeroPhaseKernel(b, l, h)
        integer, intent(in) :: l
        real(kind = dbl), intent(in) :: b(l)
        real(kind = dbl), intent(out) :: h(-(l - 1):(l - 1))
        integer :: j
        integer :: k

        h = 0d0
        do j = 1, l
            do k = 1, l
                h(j - k) = h(j - k) + b(j) * b(k)
            end do
        end do
    end subroutine ZeroPhaseKernel

    !***************************************************************************
    !> \brief y(i) = sum_m h(m) x(i + m) for i in [olo, ohi].
    !>
    !> x is given on [lo, hi]. A point whose support reaches outside it, or
    !> holds an error-coded value, is err: a covariance with no valid pair
    !> cannot be smoothed into its neighbours.
    !***************************************************************************
    pure subroutine ZeroPhaseSmooth(x, lo, hi, h, hw, err, olo, ohi, y)
        integer, intent(in) :: lo
        integer, intent(in) :: hi
        integer, intent(in) :: hw
        integer, intent(in) :: olo
        integer, intent(in) :: ohi
        real(kind = dbl), intent(in) :: x(lo:hi)
        real(kind = dbl), intent(in) :: h(-hw:hw)
        real(kind = dbl), intent(in) :: err
        real(kind = dbl), intent(out) :: y(olo:ohi)
        integer :: i
        integer :: m
        real(kind = dbl) :: s
        logical :: bad

        do i = olo, ohi
            if (i - hw < lo .or. i + hw > hi) then
                y(i) = err
                cycle
            end if
            s = 0d0
            bad = .false.
            do m = -hw, hw
                if (x(i + m) == err) then
                    bad = .true.
                    exit
                end if
                s = s + h(m) * x(i + m)
            end do
            if (bad) then
                y(i) = err
            else
                y(i) = s
            end if
        end do
    end subroutine ZeroPhaseSmooth

    !***************************************************************************
    !> \brief How many lags the covariance must extend past each end of the
    !>        search window for FindCovPeak to see interior values only.
    !***************************************************************************
    pure function PeakPadding(smooth_len, sign_mode) result(pad)
        integer, intent(in) :: smooth_len
        integer, intent(in) :: sign_mode
        integer :: pad

        pad = 0
        if (smooth_len > 1) pad = smooth_len - 1
        !> The derivative is one lag shorter, and its own moving average
        !> reaches DerivFilterLen - 1 further.
        if (sign_mode == PeakDerivative) pad = max(pad, DerivFilterLen)
    end function PeakPadding

    !***************************************************************************
    !> \brief The lag of the covariance peak inside [slo, shi].
    !>
    !> cov is the raw cross-covariance on [lo, hi], which must reach at least
    !> PeakPadding(smooth_len, sign_mode) past [slo, shi] on both sides for
    !> the result to equal innFLUX's; nearer an end the smoothed values are
    !> error-coded and simply not chosen.
    !>
    !> smooth_len > 1 smooths with a normalised Hamming window of that many
    !> lags, innFLUX's COVPEAK_FILTER_LENGTH. sign_mode picks the extremum:
    !>   PeakAbs        the largest |cov|, the earliest lag winning a tie, as
    !>                  CovMax always has;
    !>   PeakPositive   the largest cov;
    !>   PeakDerivative innFLUX's rule: the raw function's first difference,
    !>                  smoothed by a DerivFilterLen moving average, rises to
    !>                  its maximum before it falls to its minimum on a
    !>                  positive peak, so the maximum is taken if the
    !>                  derivative's maximum comes first, else the minimum.
    !>                  The rule reads the same on EddyFlow's lag axis as on
    !>                  MATLAB's reversed one: both orderings flip together.
    !> With PeakPositive and PeakDerivative a tie goes to the latest lag,
    !> which is the one MATLAB's max/min - first index on its reversed axis -
    !> returns.
    !>
    !> debaseline subtracts the chord joining the smoothed function's values
    !> at slo and shi before the extremum is taken, CovMax's baseline
    !> criterion. found is false when no lag in the window has a value.
    !***************************************************************************
    pure subroutine FindCovPeak(cov, lo, hi, slo, shi, smooth_len, sign_mode, &
        debaseline, err, lag, found)
        integer, intent(in) :: lo
        integer, intent(in) :: hi
        integer, intent(in) :: slo
        integer, intent(in) :: shi
        integer, intent(in) :: smooth_len
        integer, intent(in) :: sign_mode
        logical, intent(in) :: debaseline
        real(kind = dbl), intent(in) :: cov(lo:hi)
        real(kind = dbl), intent(in) :: err
        integer, intent(out) :: lag
        logical, intent(out) :: found
        real(kind = dbl) :: sm(slo:shi)
        real(kind = dbl) :: score(slo:shi)
        real(kind = dbl), allocatable :: b(:)
        real(kind = dbl), allocatable :: h(:)
        real(kind = dbl), allocatable :: d(:)
        real(kind = dbl), allocatable :: ds(:)
        real(kind = dbl) :: best
        real(kind = dbl) :: dmax
        real(kind = dbl) :: dmin
        integer :: imax
        integer :: imin
        integer :: k
        logical :: positive
        logical :: any_d

        lag = slo
        found = .false.
        if (shi < slo) return

        !> The function the extremum is taken on.
        if (smooth_len > 1) then
            allocate(b(smooth_len), h(-(smooth_len - 1):(smooth_len - 1)))
            call HammingKernel(smooth_len, b)
            call ZeroPhaseKernel(b, smooth_len, h)
            call ZeroPhaseSmooth(cov, lo, hi, h, smooth_len - 1, err, slo, shi, sm)
        else
            do k = slo, shi
                if (k < lo .or. k > hi) then
                    sm(k) = err
                else
                    sm(k) = cov(k)
                end if
            end do
        end if

        score = sm
        if (debaseline .and. shi - slo >= 2) then
            if (sm(slo) /= err .and. sm(shi) /= err) then
                do k = slo, shi
                    if (sm(k) == err) cycle
                    score(k) = sm(k) - (sm(slo) + (sm(shi) - sm(slo)) &
                        * dble(k - slo) / dble(shi - slo))
                end do
            end if
        end if

        !> Which sign the peak has.
        positive = .true.
        if (sign_mode == PeakDerivative) then
            !> d(j) = cov(j + 1) - cov(j), the first difference, on [lo, hi - 1].
            allocate(d(lo:hi - 1), ds(slo:shi))
            do k = lo, hi - 1
                if (cov(k) == err .or. cov(k + 1) == err) then
                    d(k) = err
                else
                    d(k) = cov(k + 1) - cov(k)
                end if
            end do
            if (allocated(b)) deallocate(b)
            if (allocated(h)) deallocate(h)
            allocate(b(DerivFilterLen))
            allocate(h(-(DerivFilterLen - 1):(DerivFilterLen - 1)))
            b = 1d0 / dble(DerivFilterLen)
            call ZeroPhaseKernel(b, DerivFilterLen, h)
            call ZeroPhaseSmooth(d, lo, hi - 1, h, DerivFilterLen - 1, err, slo, shi, ds)
            any_d = .false.
            imax = slo
            imin = slo
            dmax = 0d0
            dmin = 0d0
            do k = slo, shi
                if (ds(k) == err) cycle
                if (.not. any_d) then
                    any_d = .true.
                    dmax = ds(k)
                    dmin = ds(k)
                    imax = k
                    imin = k
                    cycle
                end if
                !> Latest lag on a tie: MATLAB's first index on its axis.
                if (ds(k) >= dmax) then
                    dmax = ds(k)
                    imax = k
                end if
                if (ds(k) <= dmin) then
                    dmin = ds(k)
                    imin = k
                end if
            end do
            if (any_d) positive = imax < imin
        end if

        best = 0d0
        do k = slo, shi
            if (score(k) == err .or. sm(k) == err) cycle
            select case (sign_mode)
                case (PeakAbs)
                    !> Strictly greater: the earliest lag wins, as in CovMax.
                    if (.not. found .or. dabs(score(k)) > best) then
                        best = dabs(score(k))
                        lag = k
                        found = .true.
                    end if
                case default
                    if (positive) then
                        if (.not. found .or. score(k) >= best) then
                            best = score(k)
                            lag = k
                            found = .true.
                        end if
                    else
                        if (.not. found .or. score(k) <= best) then
                            best = score(k)
                            lag = k
                            found = .true.
                        end if
                    end if
            end select
        end do
    end subroutine FindCovPeak

    !***************************************************************************
    !> \brief Steady-state test as innFLUX computes it (steady_state_test.m;
    !>        Foken and Wichura 1996): |mean of the nsub sub-period covariances
    !>        - the whole period's| / |the whole period's|, on w and c with c
    !>        shifted by lag rows. innFLUX uses five sub-periods.
    !>
    !> Pairs where either value is err are skipped, and each covariance takes
    !> the mean of its own valid pairs. err when the whole-period covariance is
    !> zero or a sub-period has fewer than two pairs.
    !***************************************************************************
    pure function FokenRelDev(w, c, n, lag, nsub, err) result(sst)
        integer, intent(in) :: n
        integer, intent(in) :: lag
        integer, intent(in) :: nsub
        real(kind = dbl), intent(in) :: w(n)
        real(kind = dbl), intent(in) :: c(n)
        real(kind = dbl), intent(in) :: err
        real(kind = dbl) :: sst
        real(kind = dbl) :: total
        real(kind = dbl) :: sub
        real(kind = dbl) :: msub
        integer :: k0
        integer :: k1
        integer :: m
        integer :: j
        integer :: a
        integer :: b
        logical :: ok

        sst = err
        if (nsub < 1 .or. abs(lag) >= n) return
        !> The pairs w(k), c(k + lag) that exist.
        k0 = max(1, 1 - lag)
        k1 = min(n, n - lag)
        m = (k1 - k0 + 1) / nsub
        if (m < 2) return
        call PairCov(k0, k1, total, ok)
        if (.not. ok .or. total == 0d0) return
        msub = 0d0
        do j = 1, nsub
            a = k0 + (j - 1) * m
            b = a + m - 1
            call PairCov(a, b, sub, ok)
            if (.not. ok) return
            msub = msub + sub
        end do
        msub = msub / dble(nsub)
        sst = abs((msub - total) / total)

    contains

        pure subroutine PairCov(ka, kb, cv, okk)
            integer, intent(in) :: ka
            integer, intent(in) :: kb
            real(kind = dbl), intent(out) :: cv
            logical, intent(out) :: okk
            real(kind = dbl) :: sw
            real(kind = dbl) :: sc
            real(kind = dbl) :: swc
            integer :: k
            integer :: np

            sw = 0d0
            sc = 0d0
            swc = 0d0
            np = 0
            do k = ka, kb
                if (w(k) == err .or. c(k + lag) == err) cycle
                np = np + 1
                sw = sw + w(k)
                sc = sc + c(k + lag)
                swc = swc + w(k) * c(k + lag)
            end do
            okk = np >= 2
            cv = 0d0
            if (okk) cv = swc / dble(np) - (sw / dble(np)) * (sc / dble(np))
        end subroutine PairCov
    end function FokenRelDev

end module m_lagpeak_core
