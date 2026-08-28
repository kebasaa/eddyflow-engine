!***************************************************************************
! test_rflux_diagnostics.f90
! ---------------------------
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
! \brief       Extra raw-signal instrument-malfunction diagnostics, ported
!              from RFlux's inst_prob_test() (Vitale et al. 2020,
!              Biogeosciences): lag-1 autocorrelation (AL1), the
!              discrete/dominant-value test (DDI), and robust-scaled spike
!              counts on the raw series (HF5/HF10) and on its first
!              differences with flat-lined transitions excluded (HD5/HD10).
! \author      Jonathan Muller, ETH Zurich
! \note        RFlux scales HF5/HF10/HD5/HD10 by the Qn estimator
!              (Rousseeuw & Croux 1993). This engine already carries a Qn
!              implementation (qn_estimator.f90), but its all-pairs
!              construction is O(n^2) in memory and time - at a 20 Hz,
!              30-minute period (36000 samples) that is ~6.5e8 pairs, well
!              past what a per-period raw-level test can spend. This test
!              uses the median/MAD scale (median of |x - median(x)|,
!              normalised by 0.6745 as in TestSpikeDetectionMauder13) instead
!              - the same O(n log n) robust-scale idiom already used
!              elsewhere in this file's neighbourhood, at the same
!              statistical role Qn plays in RFlux. A true O(n log n) Qn is a
!              separate follow-up.
!***************************************************************************
subroutine TestRFluxDiagnostics(Set, N, printout)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    logical, intent(in) :: printout
    real(kind = dbl), intent(in) :: Set(N, E2NumVar)
    !> local variables
    integer :: icol

    if (printout) write(*, '(a)', advance = 'no') '   RFlux raw-signal diagnostics test..'
    if (printout) write(ulog, '(a)', advance = 'no') '   RFlux raw-signal diagnostics test..'

    do icol = u, ts
        call RFluxDiagColumn(Set(:, icol), N, icol)
    end do
    do icol = firstGas, lastGas
        if (E2Col(icol)%present) then
            call RFluxDiagColumn(Set(:, icol), N, icol)
        else
            Essentials%AL1(icol)  = error
            Essentials%DDI(icol)  = error
            Essentials%HF5(icol)  = error
            Essentials%HF10(icol) = error
            Essentials%HD5(icol)  = error
            Essentials%HD10(icol) = error
        end if
    end do

    if (printout) write(*,'(a)') ' Done.'
    if (printout) write(ulog,'(a)') ' Done.'
end subroutine TestRFluxDiagnostics


!> One raw column's worth of RFlux-derived diagnostics, written into
!> Essentials(icol). Mirrors inst_prob_test()'s per-variable computations.
subroutine RFluxDiagColumn(xcol, N, icol)
    use m_rp_global_var
    implicit none
    integer, intent(in) :: N
    integer, intent(in) :: icol
    real(kind = dbl), intent(in) :: xcol(N)
    !> local variables
    integer :: i, ii, lN
    real(kind = dbl), allocatable :: valid(:), sorted_valid(:)
    real(kind = dbl) :: col_mean, dip_stat
    real(kind = dbl), external :: LaggedCovarianceNoError
    real(kind = dbl) :: variance(1)
    real(kind = dbl) :: mean1(1)

    !> Count of usable (non error-coded) samples
    lN = count(xcol /= error)
    if (lN < 3) then
        Essentials%AL1(icol)  = error
        Essentials%DDI(icol)  = error
        Essentials%HF5(icol)  = error
        Essentials%HF10(icol) = error
        Essentials%HD5(icol)  = error
        Essentials%HD10(icol) = error
        Essentials%DIP(icol)  = error
        return
    end if

    !> AL1: lag-1 autocorrelation = lag-1 autocovariance / variance
    call AverageNoError(reshape(xcol, [N, 1]), N, 1, mean1, error)
    call StDevNoError(reshape(xcol, [N, 1]), N, 1, variance, error)
    col_mean = mean1(1)
    if (variance(1) /= error .and. variance(1) > 0d0) then
        Essentials%AL1(icol) = LaggedCovarianceNoError(xcol, xcol, N, 1, error) &
            / (variance(1) ** 2)
    else
        Essentials%AL1(icol) = error
    end if

    !> Squeeze out the error-coded samples once; every remaining diagnostic
    !> works on this compact, gap-free copy.
    allocate(valid(lN))
    ii = 0
    do i = 1, N
        if (xcol(i) == error) cycle
        ii = ii + 1
        valid(ii) = xcol(i)
    end do

    call RFluxDDI(valid, lN, Essentials%DDI(icol))
    call RFluxSpikeCounts(valid, lN, col_mean, &
        Essentials%HF5(icol), Essentials%HF10(icol), &
        Essentials%HD5(icol), Essentials%HD10(icol))

    !> DIP: Hartigan's dip test p-value on the fluctuations (RFlux: on
    !> flucts <- x - mean(x)). The dip statistic is shift-invariant - it
    !> depends only on the relative spacing of sorted values - so sorting
    !> `valid` directly gives the identical result without subtracting the
    !> mean first.
    allocate(sorted_valid(lN))
    sorted_valid = valid
    call HPSORT(lN, sorted_valid)
    call DipStatistic(sorted_valid, lN, dip_stat)
    call DipPvalue(lN, dip_stat, Essentials%DIP(icol))
    deallocate(sorted_valid)

    deallocate(valid)
end subroutine RFluxDiagColumn


!> Discrete/Dominant-value test: the largest histogram bin count under a
!> Freedman-Diaconis bin width. A high count means many samples collapsing
!> onto few values - flat-lining or a coarsened instrument resolution.
!> RFlux: DDI <- max(hist(x, breaks="FD")$counts), or n if x is constant.
subroutine RFluxDDI(x, n, ddi)
    use m_rp_global_var
    implicit none
    integer, intent(in) :: n
    real(kind = dbl), intent(in) :: x(n)
    real(kind = dbl), intent(out) :: ddi
    real(kind = dbl) :: xcopy(n)
    real(kind = dbl) :: q1(1), q3(1), iqr, binwidth, xmin, xmax
    integer :: nbins, bin, i
    integer, allocatable :: counts(:)
    real(kind = dbl), external :: quantile_sas5

    xmin = minval(x)
    xmax = maxval(x)
    if (xmax <= xmin) then
        !> A constant series is maximally discrete - every sample is the
        !> same "bin".
        ddi = dble(n)
        return
    end if

    xcopy = x
    q1(1) = quantile_sas5(xcopy, n, 0.25d0)
    xcopy = x
    q3(1) = quantile_sas5(xcopy, n, 0.75d0)
    iqr = q3(1) - q1(1)

    if (iqr <= 0d0) then
        binwidth = (xmax - xmin) / dble(n)
    else
        binwidth = 2d0 * iqr / (dble(n) ** (1d0 / 3d0))
    end if
    if (binwidth <= 0d0) then
        ddi = dble(n)
        return
    end if

    nbins = max(1, ceiling((xmax - xmin) / binwidth))
    allocate(counts(nbins))
    counts = 0
    do i = 1, n
        bin = min(nbins, 1 + int((x(i) - xmin) / binwidth))
        counts(bin) = counts(bin) + 1
    end do
    ddi = dble(maxval(counts))
    deallocate(counts)
end subroutine RFluxDDI


!> HF5/HF10: count of raw-fluctuation samples beyond 5/10 robust sigmas.
!> HD5/HD10: count of first-difference samples beyond 5/10 robust sigmas of
!> the differences, with transitions into a repeated (flat-lined) value
!> excluded from the scale estimate - RFlux replaces the sample right after
!> a near-zero difference with NA before differencing again, so a flat-lined
!> stretch cannot deflate its own detection threshold.
!> RFlux requires more than 1000 usable difference samples to trust the
!> scale; below that HD5/HD10 report as unavailable, same as inst_prob_test.
subroutine RFluxSpikeCounts(x, n, x_mean, hf5, hf10, hd5, hd10)
    use m_rp_global_var
    implicit none
    integer, intent(in) :: n
    real(kind = dbl), intent(in) :: x(n)
    real(kind = dbl), intent(in) :: x_mean
    real(kind = dbl), intent(out) :: hf5, hf10, hd5, hd10
    real(kind = dbl) :: fluct(n)
    real(kind = dbl) :: sigma_f, sigma_d, medx, mad
    real(kind = dbl) :: d0(n - 1)
    real(kind = dbl), allocatable :: d1(:)
    logical :: repeat_next(n)
    integer :: i, nd1

    !> HF5 / HF10, scaled by the robust sigma of the fluctuations
    fluct = x - x_mean
    call RobustSigma(fluct, n, sigma_f)
    hf5  = dble(count(dabs(fluct) > 5d0 * sigma_f))
    hf10 = dble(count(dabs(fluct) > 10d0 * sigma_f))

    if (n < 2) then
        hd5 = error
        hd10 = error
        return
    end if

    do i = 1, n - 1
        d0(i) = x(i + 1) - x(i)
    end do

    !> Mark the sample after each near-zero step as part of a repeated run,
    !> so the scale estimate for HD5/HD10 is not deflated by flat-lining.
    repeat_next = .false.
    do i = 1, n - 1
        if (dabs(d0(i)) < 1d-3) repeat_next(i + 1) = .true.
    end do

    allocate(d1(n - 1))
    nd1 = 0
    do i = 1, n - 1
        if (repeat_next(i) .or. repeat_next(i + 1)) cycle
        nd1 = nd1 + 1
        d1(nd1) = x(i + 1) - x(i)
    end do

    if (nd1 > 1000) then
        call RobustSigma(d1(1:nd1), nd1, sigma_d)
        hd5  = dble(count(dabs(d0) > 5d0 * sigma_d))
        hd10 = dble(count(dabs(d0) > 10d0 * sigma_d))
    else
        hd5 = error
        hd10 = error
    end if
    deallocate(d1)
end subroutine RFluxSpikeCounts


!> Median/MAD robust scale, normalised to be consistent at the normal
!> distribution (the same 0.6745 convention TestSpikeDetectionMauder13
!> uses), floored at 0.01 as RFlux floors its Qn-based sigma.
subroutine RobustSigma(x, n, sigma)
    use m_rp_global_var
    implicit none
    integer, intent(in) :: n
    real(kind = dbl), intent(in) :: x(n)
    real(kind = dbl), intent(out) :: sigma
    real(kind = dbl) :: xcopy(n), medx, mad

    xcopy = x
    call median(xcopy, n, medx)
    xcopy = dabs(x - medx)
    call median(xcopy, n, mad)
    sigma = max(0.01d0, mad / 0.6745d0)
end subroutine RobustSigma
