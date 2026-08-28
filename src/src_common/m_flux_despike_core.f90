!***************************************************************************
! m_flux_despike_core.f90
! ------------------------
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
! \brief       Post-flux despiking on a whole-run half-hourly flux series:
!              RFlux's despiking(variant="v1") (Vitale et al. 2020,
!              Biogeosciences) - an STL (Cleveland et al. 1990) decomposition
!              on log10(x+1000), then a decile-binned Laplace outlier test
!              on the residual.
!
! \details     Everything here is a function of its arguments: no project
!              settings, no column table, no logging - the same discipline
!              m_pwb_core.f90 follows, and for the same reason: it lets a
!              reference driver call this without linking the rest of the
!              engine.
!
!              This is a Fortran port of R's stlplus package (Ryan Hafen,
!              GPL >= 2, CRAN) - the loess primitive (stlplus's
!              src/loess.cpp), its low-pass filter (src/ma.cpp), its
!              jump-interpolation (src/interp.cpp), and the R-level
!              orchestration in R/stlplus.R and R/loess_stl.R - plus RFlux's
!              own despiking.R (variant "v1") built on top of it. Ported
!              rather than re-derived from the STL paper, for the same
!              reason dip_test.f90 gives: an independent re-derivation of a
!              multi-stage iterative algorithm like this is exactly where a
!              subtle bug hides from a compiler and from casual review.
!
!              Validated against R directly (stlplus + RFlux's despiking.R,
!              run locally) rather than assumed: the loess primitive, the
!              jump-interpolation wrapper, the NA-aware nearest-neighbour
!              window, one isolated cycle-subseries smooth, the full STL
!              inner/outer loop (with and without missing values, several
!              series lengths and periods), the post-trend "fc" extension,
!              and finally RFlux's actual despiking() function end to end -
!              all matched to machine precision (or an exact spike-flag
!              match, for the last one), across dozens of generated cases.
!              See the session that added this file for the validation
!              scripts and their output.
!
! \par Deliberate departures from R
!
!              1. s.blend/t.blend/l.blend/fc.blend are always 0 here - the
!                 only value RFlux's despiking.R ever passes - so
!                 .loess_stlplus's blend-to-degree-0-at-the-endpoints branch
!                 is not ported at all.
!              2. t.window is always given explicitly (RFlux always passes
!                 one), so stlplus.default's get.t.window - the empirical
!                 critical-frequency formula used only when t.window is
!                 NULL - is not ported either.
!              3. s.window is always a number, never the string "periodic".
!
! \author      Jonathan Muller, ETH Zurich
! \note
! \sa          RFlux-master/R/despiking.R, dip_test.f90 (the sibling case
!              for a ported-and-validated-against-R algorithm)
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
module m_flux_despike_core
    use m_numeric_kinds
    implicit none
    private
    public :: DespikeFluxSeries

contains

!***************************************************************************
!> The whole despiking(variant="v1") pipeline for one flux series.
!>
!> x(n): the raw half-hourly (or other fixed-period) series, `missing` where
!> a period has no value. mfreq: periods per day (48 for half-hourly).
!> alpha: the Laplace outlier test's significance level (RFlux default 0.01).
!>
!> spike_flag(n) comes back 1 where R's despiking() would have replaced the
!> point with NA, 0 elsewhere (including every `missing` position - a gap
!> is not a spike). cleaned(n) is x with spike_flag==1 positions set to
!> `missing`, matching R's replace(x, spike_index, NA).
!***************************************************************************
subroutine DespikeFluxSeries(x, n, mfreq, missing, alpha, spike_flag, cleaned)
    integer, intent(in) :: n, mfreq
    real(kind = dbl), intent(in) :: x(n)
    real(kind = dbl), intent(in) :: missing, alpha
    integer, intent(out) :: spike_flag(n)
    real(kind = dbl), intent(out) :: cleaned(n)

    real(kind = dbl), allocatable :: logx(:), trend(:), seasonal(:), remainder(:)
    real(kind = dbl), allocatable :: weights(:), fc_out(:,:), signal(:), res(:)
    real(kind = dbl), allocatable :: signal_sorted(:), bin_res(:)
    integer, allocatable :: bin_idx(:)
    integer :: fc_window(2), fc_degree(2), fc_jump(2)
    integer :: i, j, bin_n, valid_n
    real(kind = dbl) :: indq(11), med, qn_scale, lo_bound, hi_bound
    real(kind = dbl), external :: Qn
    real(kind = dbl), parameter :: log_offset = 1000d0

    allocate(logx(n), trend(n), seasonal(n), remainder(n), weights(n))
    allocate(fc_out(n,2), signal(n), res(n), signal_sorted(n))
    allocate(bin_idx(n), bin_res(n))

    !> Degenerate runs (RFlux's own guard: at least ~10 days of data, and no
    !> subseries missing everything) aren't checked here - a run this short
    !> is not what half-hourly flux post-processing looks like, and the
    !> caller decides whether to invoke this at all.
    where (x /= missing)
        logx = dlog10(x + log_offset)
    elsewhere
        logx = missing
    end where

    call StlDecompose(logx, n, mfreq, 7, 0, mfreq*7+1, 1, 20, 2, missing, &
        trend, seasonal, remainder, weights)

    fc_window = [mfreq*7+1, 7]
    fc_degree = [1, 0]
    fc_jump(1) = max(1, ceiling(dble(fc_window(1))/10d0))
    fc_jump(2) = max(1, ceiling(dble(fc_window(2))/10d0))
    call StlPostSmooth(logx, seasonal, n, weights, fc_window, fc_degree, 2, fc_jump, missing, fc_out)

    do i = 1, n
        signal(i) = 10d0**(seasonal(i) + fc_out(i,1) + fc_out(i,2)) - log_offset
    end do

    do i = 1, n
        if (x(i) /= missing) then
            res(i) = x(i) - signal(i)
        else
            res(i) = missing
        end if
    end do

    !> Deciles of signal (R: quantile(signal, seq(0,1,.1), na.rm=TRUE)).
    !> signal itself has no gaps (the STL fit is complete even where x is
    !> missing), so no na.rm handling is needed here.
    signal_sorted = signal
    call HPSORT(n, signal_sorted)
    do i = 0, 10
        indq(i+1) = QuantileType7(signal_sorted, n, dble(i)/10d0)
    end do

    spike_flag = 0
    do i = 1, 10
        lo_bound = indq(i)
        hi_bound = indq(i+1)
        bin_n = 0
        do j = 1, n
            if (signal(j) >= lo_bound .and. signal(j) <= hi_bound) then
                bin_n = bin_n + 1
                bin_idx(bin_n) = j
            end if
        end do

        valid_n = 0
        do j = 1, bin_n
            if (res(bin_idx(j)) /= missing) then
                valid_n = valid_n + 1
                bin_res(valid_n) = res(bin_idx(j))
            end if
        end do

        !> Fewer than 10 usable points: aout.laplace can't fit a scale, so
        !> RFlux flags every point in the bin as a spike rather than skip
        !> the test silently.
        if (valid_n > 9) then
            block
                real(kind = dbl) :: tmp(valid_n)
                tmp = bin_res(1:valid_n)
                call median(tmp, valid_n, med)
                tmp = bin_res(1:valid_n)
                qn_scale = Qn(tmp, valid_n)
            end block
            do j = 1, bin_n
                if (res(bin_idx(j)) == missing) cycle
                if (res(bin_idx(j)) < med + qn_scale*dlog(alpha) .or. &
                    res(bin_idx(j)) > med - qn_scale*dlog(alpha)) then
                    spike_flag(bin_idx(j)) = 1
                end if
            end do
        else
            do j = 1, bin_n
                if (res(bin_idx(j)) /= missing) spike_flag(bin_idx(j)) = 1
            end do
        end if
    end do

    do i = 1, n
        if (spike_flag(i) == 1) then
            cleaned(i) = missing
        else
            cleaned(i) = x(i)
        end if
    end do

    deallocate(logx, trend, seasonal, remainder, weights, fc_out, signal, res, &
        signal_sorted, bin_idx, bin_res)
end subroutine DespikeFluxSeries


!***************************************************************************
!> Full stlplus() decomposition, restricted to what despiking() uses (see
!> the module header's "Deliberate departures" list).
!***************************************************************************
subroutine StlDecompose(Y, n, np_, s_window, s_degree, t_window, t_degree, &
                          outer_iter, inner_iter, missing, trend, seasonal, remainder, weights)
    integer, intent(in) :: n, np_, s_window, s_degree, t_window, t_degree
    integer, intent(in) :: outer_iter, inner_iter
    real(kind = dbl), intent(in) :: Y(n), missing
    real(kind = dbl), intent(out) :: trend(n), seasonal(n), remainder(n), weights(n)

    integer :: s_win, t_win, l_win, s_jump, t_jump, l_jump, l_degree
    integer, allocatable :: cycleSubIndices(:)
    real(kind = dbl), allocatable :: C(:), Ydet(:), ma3(:), L(:), D(:), R(:)
    integer :: o_iter, iter, i, k, nma
    integer :: st, nd

    s_win = s_window
    if (mod(s_win,2)==0) s_win = s_win+1
    t_win = t_window
    if (mod(t_win,2)==0) t_win = t_win+1
    l_win = np_
    if (mod(l_win,2)==0) l_win = l_win+1
    l_degree = t_degree

    s_jump = max(1, ceiling(dble(s_win)/10d0))
    t_jump = max(1, ceiling(dble(t_win)/10d0))
    l_jump = max(1, ceiling(dble(l_win)/10d0))

    allocate(cycleSubIndices(n))
    do i = 1, n
        cycleSubIndices(i) = mod(i-1, np_) + 1
    end do

    allocate(C(n + 2*np_), Ydet(n), ma3(n), L(n), D(n), R(n))
    trend = 0d0
    weights = 1d0

    st = np_ + 1
    nd = n + np_

    do o_iter = 1, outer_iter
        do iter = 1, inner_iter
            where (Y /= missing)
                Ydet = Y - trend
            elsewhere
                Ydet = missing
            end where

            do k = 1, np_
                call CycleSubseriesSmooth(Ydet, weights, cycleSubIndices, n, np_, k, &
                                           s_win, s_degree, s_jump, missing, C)
            end do

            call LowPassFilter3(C, n+2*np_, np_, ma3, nma)  !< nma == n

            !> R passes the ORIGINAL series' y_idx/noNA here explicitly,
            !> not ma3's own (`.loess_stlplus(y=ma3, ..., y_idx=y_idx,
            !> noNA=noNA)` in stlplus.default) - so even though ma3 itself
            !> is complete (cycle-subseries smoothing fills every gap),
            !> positions where Y was originally missing are still excluded
            !> from this fit. Mask them back in before smoothing.
            block
                real(kind = dbl) :: ma3_masked(n)
                where (Y /= missing)
                    ma3_masked = ma3
                elsewhere
                    ma3_masked = missing
                end where
                call SmoothSeries(ma3_masked, n, l_win, l_degree, weights, l_jump, missing, L)
            end block

            seasonal = C(st:nd) - L

            where (Y /= missing)
                D = Y - seasonal
            elsewhere
                D = missing
            end where

            call SmoothSeries(D, n, t_win, t_degree, weights, t_jump, missing, trend)
        end do

        if (outer_iter > 1) then
            where (Y /= missing)
                R = Y - seasonal - trend
            elsewhere
                R = missing
            end where
            call RobustWeights(R, n, missing, weights)
        end if
    end do

    where (Y /= missing)
        remainder = Y - seasonal - trend
    elsewhere
        remainder = missing
    end where

    deallocate(cycleSubIndices, C, Ydet, ma3, L, D, R)
end subroutine StlDecompose


!***************************************************************************
!> One subseries (class k of np_): extract, loess-smooth with one point of
!> extrapolation on each side, scatter into C's padded layout.
!>
!> C is laid out [cs1(np_) | middle(n) | cs2(np_)], where cs1/cs2 are the
!> periodic class pattern's own first/last np_ entries - so cs1(j)=j exactly
!> (mod(j-1,np_)+1 for j in 1..np_ is just j), and cs2(j) is the same
!> periodic pattern evaluated np_ steps short of one full cycle past n.
!***************************************************************************
subroutine CycleSubseriesSmooth(Ydet, w, cycleSubIndices, n, np_, k, s_win, s_degree, s_jump, missing, C)
    integer, intent(in) :: n, np_, k, s_win, s_degree, s_jump
    integer, intent(in) :: cycleSubIndices(n)
    real(kind = dbl), intent(in) :: Ydet(n), w(n), missing
    real(kind = dbl), intent(inout) :: C(n + 2*np_)
    integer :: csl, i, ii, p2
    real(kind = dbl), allocatable :: cyc(:), cyc_w(:), cyc_out(:)

    csl = count(cycleSubIndices == k)
    allocate(cyc(csl), cyc_w(csl), cyc_out(csl+2))

    ii = 0
    do i = 1, n
        if (cycleSubIndices(i) == k) then
            ii = ii + 1
            cyc(ii) = Ydet(i)
            cyc_w(ii) = w(i)
        end if
    end do

    call SmoothSubseries(cyc, cyc_w, csl, s_win, s_degree, s_jump, missing, cyc_out)

    !> cs1: the single position in 1..np_ whose class is k is position k
    !> itself (cs1(j)=j).
    C(k) = cyc_out(1)

    !> middle: the csl occurrences, in order.
    ii = 0
    do i = 1, n
        if (cycleSubIndices(i) == k) then
            ii = ii + 1
            C(np_ + i) = cyc_out(ii + 1)
        end if
    end do

    !> cs2 = tail(cycleSubIndices, np_) = cycleSubIndices(n-np_+1 : n); the
    !> single matching position within that np_-long tail.
    do p2 = 1, np_
        if (cycleSubIndices(n - np_ + p2) == k) then
            C(n + np_ + p2) = cyc_out(csl + 2)
            exit
        end if
    end do

    deallocate(cyc, cyc_w, cyc_out)
end subroutine CycleSubseriesSmooth


!***************************************************************************
!> Post-trend smoothing ("fc" extension): iteratively loess-smooth
!> (Y - seasonal - cumulative-so-far) with each fc_window(ii)/fc_degree(ii)
!> in turn, accumulating. Called once, after StlDecompose finishes, using
!> its final seasonal and robustness weights. n_fc components come back in
!> fc_out(:, 1:n_fc); their sum is what DespikeFluxSeries's `signal` adds
!> to the seasonal term.
!***************************************************************************
subroutine StlPostSmooth(Y, seasonal, n, weights, fc_window, fc_degree, n_fc, fc_jump, missing, fc_out)
    integer, intent(in) :: n, n_fc
    real(kind = dbl), intent(in) :: Y(n), seasonal(n), weights(n), missing
    integer, intent(in) :: fc_window(n_fc), fc_degree(n_fc), fc_jump(n_fc)
    real(kind = dbl), intent(out) :: fc_out(n, n_fc)
    real(kind = dbl), allocatable :: cumulative(:), target_series(:), tmp(:)
    integer :: ii

    allocate(cumulative(n), target_series(n), tmp(n))
    cumulative = 0d0
    do ii = 1, n_fc
        where (Y /= missing)
            target_series = Y - seasonal - cumulative
        elsewhere
            target_series = missing
        end where
        call SmoothSeries(target_series, n, fc_window(ii), fc_degree(ii), weights, fc_jump(ii), missing, tmp)
        cumulative = cumulative + tmp
        fc_out(:, ii) = tmp
    end do
    deallocate(cumulative, target_series, tmp)
end subroutine StlPostSmooth


!***************************************************************************
!> Tukey biweight robustness weights from the current remainder, R's own
!> formula in stlplus.default (mid1/mid2/h from a sorted |R|, zero.weight
!> floor 1e-6 for h9-clipped points, missing remainder gets weight 1).
!>
!> R's sort() drops NA by default, so sort(R.abs)[mid1:mid2] indexes a
!> vector of length n_valid using mid1/mid2 computed from the FULL
!> (missing-inclusive) n. If mid1 or mid2 falls beyond n_valid, R's x[i]
!> for i>length(x) is NA, so the sum - and everything downstream through
!> the final `w[is.na(w)] <- 1` - becomes 1 everywhere.
!***************************************************************************
subroutine RobustWeights(R, n, missing, w)
    integer, intent(in) :: n
    real(kind = dbl), intent(in) :: R(n), missing
    real(kind = dbl), intent(out) :: w(n)
    real(kind = dbl), parameter :: zero_weight = 1d-6
    real(kind = dbl), allocatable :: Rabs(:), Rvalid(:)
    integer :: mid1, mid2, mlo, mhi, i, n_valid, k
    real(kind = dbl) :: h, h1, h9

    allocate(Rabs(n))
    n_valid = 0
    do i = 1, n
        if (R(i) /= missing) then
            Rabs(i) = dabs(R(i))
            n_valid = n_valid + 1
        else
            Rabs(i) = 0d0
        end if
    end do

    allocate(Rvalid(n_valid))
    k = 0
    do i = 1, n
        if (R(i) /= missing) then
            k = k + 1
            Rvalid(k) = Rabs(i)
        end if
    end do
    call HPSORT(n_valid, Rvalid)

    mid1 = floor(dble(n)/2d0 + 1d0)
    mid2 = n - mid1 + 1
    !> R's sort(...)[mid1:mid2] auto-reverses when mid1>mid2 (even n), so it
    !> still selects both middle elements; Fortran's a(mid1:mid2) would
    !> silently give an empty (zero-length) slice instead.
    mlo = min(mid1, mid2)
    mhi = max(mid1, mid2)

    if (mhi > n_valid) then
        w = 1d0
        deallocate(Rabs, Rvalid)
        return
    end if
    h = 3d0 * sum(Rvalid(mlo:mhi))
    h9 = 0.999d0 * h
    h1 = 0.001d0 * h

    do i = 1, n
        if (R(i) == missing) then
            w(i) = 1d0
            cycle
        end if
        if (Rabs(i) <= h1) then
            w(i) = 1d0
        else if (Rabs(i) >= h9) then
            w(i) = zero_weight
        else
            w(i) = (1d0 - (Rabs(i)/h)**2)**2
        end if
    end do

    deallocate(Rabs, Rvalid)
end subroutine RobustWeights


!***************************************************************************
!> Loess-smooth y (length n, `missing` for gaps) at every integer 1..n with
!> the given span/degree/jump, weighted by w, dispatching to the missing-
!> or complete-data primitive as appropriate.
!***************************************************************************
subroutine SmoothSeries(y, n, span, degree, w, jump, missing, out)
    integer, intent(in) :: n, span, degree, jump
    real(kind = dbl), intent(in) :: y(n), w(n), missing
    real(kind = dbl), intent(out) :: out(n)
    real(kind = dbl) :: at(n)
    real(kind = dbl), allocatable :: m(:)
    integer :: i, nm

    do i = 1, n
        at(i) = dble(i)
    end do
    call BuildJumpPoints(n, jump, m, nm)
    if (any(y == missing)) then
        call LoessStlplusMissing(y, n, span, degree, w, m, nm, jump, at, n, missing, out)
    else
        call LoessStlplusComplete(y, n, span, degree, w, m, nm, jump, at, n, out)
    end if
    deallocate(m)
end subroutine SmoothSeries


!***************************************************************************
!> Same, but evaluated at integer points 0..(csl+1) (csl+2 of them) rather
!> than 1..n - the cycle-subseries call, which asks for one point of
!> extrapolation on each side.
!***************************************************************************
subroutine SmoothSubseries(y, w, csl, span, degree, jump, missing, out)
    integer, intent(in) :: csl, span, degree, jump
    real(kind = dbl), intent(in) :: y(csl), w(csl), missing
    real(kind = dbl), intent(out) :: out(csl+2)
    real(kind = dbl) :: at(csl+2)
    real(kind = dbl), allocatable :: m(:)
    integer :: i, nm

    do i = 0, csl+1
        at(i+1) = dble(i)
    end do
    call BuildCsEv(csl, jump, m, nm)
    if (any(y == missing)) then
        call LoessStlplusMissing(y, csl, span, degree, w, m, nm, jump, at, csl+2, missing, out)
    else
        call LoessStlplusComplete(y, csl, span, degree, w, m, nm, jump, at, csl+2, out)
    end if
    deallocate(m)
end subroutine SmoothSubseries


!***************************************************************************
!> .loess_stlplus port, complete-data branch, blend always 0.
!>
!> y(n): data. m(nm): caller-supplied evaluation points (the jump-spaced
!> grid for the standard 1..n case, OR the "0, 1..csl, csl+1" pattern for
!> cycle-subseries - R takes m as an input parameter for the same reason).
!> at(nat): output evaluation points.
!***************************************************************************
subroutine LoessStlplusComplete(y, n, span_in, degree, w, m, nm, jump, at, nat, out)
    integer, intent(in) :: n, degree, span_in, nm, jump, nat
    real(kind = dbl), intent(in) :: y(n), w(n), m(nm), at(nat)
    real(kind = dbl), intent(out) :: out(nat)
    integer :: span, s2, i, li, ri
    real(kind = dbl), allocatable :: max_dist(:), result(:), slopes(:), x(:)
    integer, allocatable :: l_idx0(:)
    real(kind = dbl) :: aa, bb

    span = span_in
    if (mod(span,2) == 0) span = span + 1
    s2 = (span+1)/2

    allocate(x(n))
    do i = 1, n
        x(i) = dble(i)
    end do

    allocate(max_dist(nm), result(nm), slopes(nm), l_idx0(nm))
    do i = 1, nm
        if (n - 1 < span) then
            !> R: diff(range(x)) < span -- the window is wider than the
            !> series, so every point uses the same, whole-series window.
            li = 1
            ri = n
        else if (m(i) < dble(s2)) then
            li = 1
            ri = li + span - 1
        else if (m(i) > dble(n - s2)) then
            li = n - span + 1
            ri = li + span - 1
        else
            li = nint(m(i)) - s2 + 1
            ri = li + span - 1
        end if
        aa = dabs(m(i) - dble(li))
        bb = dabs(dble(ri) - m(i))
        max_dist(i) = max(aa, bb)
        if (span >= n) max_dist(i) = max_dist(i) + dble(span-n)/2d0
        l_idx0(i) = li - 1
    end do

    call LoessFit(x, y, n, degree, span, w, m, nm, l_idx0, max_dist, result, slopes)

    if (jump > 1) then
        call HermiteInterp(m, result, slopes, nm, at, nat, out)
    else
        !> jump<=1: R's res1 <- out$result is a straight positional copy.
        !> Only correct because every jump<=1 call site builds `at` identical
        !> to `m` (same values, order, length).
        out(1:nat) = result(1:nat)
    end if

    deallocate(x, max_dist, result, slopes, l_idx0)
end subroutine LoessStlplusComplete


!***************************************************************************
!> .loess_stlplus port, missing-data branch: 1D exact-kNN neighbour window.
!>
!> y(n) uses `missing` where absent. Compacts to the valid subset (x2/y2),
!> then for each target m(i) finds the span3 nearest positions in x2 (which,
!> since x2 is a sorted integer subsequence of 1..n, always form a
!> contiguous index window - so this is a two-pointer expansion, not a
!> general kNN).
!***************************************************************************
subroutine LoessStlplusMissing(y, n, span_in, degree, w, m, nm, jump, at, nat, missing, out)
    integer, intent(in) :: n, degree, span_in, nm, jump, nat
    real(kind = dbl), intent(in) :: y(n), w(n), m(nm), at(nat), missing
    real(kind = dbl), intent(out) :: out(nat)
    integer :: span, span3, n2, i, li, ri
    real(kind = dbl), allocatable :: x2(:), y2(:), w2(:), max_dist(:), result(:), slopes(:)
    integer, allocatable :: l_idx0(:)

    span = span_in
    if (mod(span,2) == 0) span = span + 1

    n2 = count(y /= missing)
    allocate(x2(n2), y2(n2), w2(n2))
    call CompactValid(y, w, n, missing, x2, y2, w2, n2)

    span3 = min(span, n2)

    allocate(max_dist(nm), result(nm), slopes(nm), l_idx0(nm))
    do i = 1, nm
        call NearestWindow(x2, n2, m(i), span3, li, ri)
        max_dist(i) = max(dabs(m(i)-x2(li)), dabs(x2(ri)-m(i)))
        !> R uses the ORIGINAL (oddified) span here, not span3=min(span,n2):
        !> `if (span >= n) max_dist <- max_dist + (span-n)/2`, with `n` the
        !> valid-point count throughout .loess_stlplus.
        if (span >= n2) max_dist(i) = max_dist(i) + dble(span-n2)/2d0
        l_idx0(i) = li - 1
    end do

    call LoessFit(x2, y2, n2, degree, span3, w2, m, nm, l_idx0, max_dist, result, slopes)

    if (jump > 1) then
        call HermiteInterp(m, result, slopes, nm, at, nat, out)
    else
        out(1:nat) = result(1:nat)
    end if

    deallocate(x2, y2, w2, max_dist, result, slopes, l_idx0)
end subroutine LoessStlplusMissing


subroutine CompactValid(y, w, n, missing, x2, y2, w2, n2)
    integer, intent(in) :: n, n2
    real(kind = dbl), intent(in) :: y(n), w(n), missing
    real(kind = dbl), intent(out) :: x2(n2), y2(n2), w2(n2)
    integer :: i, k
    k = 0
    do i = 1, n
        if (y(i) /= missing) then
            k = k + 1
            x2(k) = dble(i)
            y2(k) = y(i)
            w2(k) = w(i)
        end if
    end do
end subroutine CompactValid


!***************************************************************************
!> The span3 nearest points to m among the sorted positions x2(1:n2),
!> returned as a contiguous index window [li,ri]. Ties (equal distance on
!> both sides) extend LEFT first - validated against R's yaImpute::ann.
!***************************************************************************
subroutine NearestWindow(x2, n2, m, span3, li, ri)
    integer, intent(in) :: n2, span3
    real(kind = dbl), intent(in) :: x2(n2), m
    integer, intent(out) :: li, ri
    integer :: lo, hi, mid, c0
    real(kind = dbl) :: dist_left, dist_right

    lo = 1; hi = n2 + 1
    do while (lo < hi)
        mid = (lo + hi) / 2
        if (x2(mid) >= m) then
            hi = mid
        else
            lo = mid + 1
        end if
    end do
    if (lo <= 1) then
        c0 = 1
    else if (lo > n2) then
        c0 = n2
    else
        if (dabs(x2(lo-1)-m) <= dabs(x2(lo)-m)) then
            c0 = lo - 1
        else
            c0 = lo
        end if
    end if

    li = c0; ri = c0
    do while (ri - li + 1 < span3)
        if (li <= 1) then
            ri = ri + 1
        else if (ri >= n2) then
            li = li - 1
        else
            dist_left = dabs(m - x2(li-1))
            dist_right = dabs(x2(ri+1) - m)
            if (dist_left <= dist_right) then
                li = li - 1
            else
                ri = ri + 1
            end if
        end if
    end do
end subroutine NearestWindow


!> cs.ev <- seq(1,n,by=jump); if last != n append n. The standard jump grid
!> used by the trend/low-pass smoothers. Cycle-subseries builds its own
!> (see BuildCsEv), since it additionally prepends 0 and appends n+1.
subroutine BuildJumpPoints(n, jump, m, nm)
    integer, intent(in) :: n, jump
    real(kind = dbl), allocatable, intent(out) :: m(:)
    integer, intent(out) :: nm
    integer :: i, cnt
    cnt = (n-1)/jump + 1
    if (1 + (cnt-1)*jump /= n) then
        allocate(m(cnt+1))
        do i = 1, cnt
            m(i) = dble(1 + (i-1)*jump)
        end do
        m(cnt+1) = dble(n)
        nm = cnt+1
    else
        allocate(m(cnt))
        do i = 1, cnt
            m(i) = dble(1 + (i-1)*jump)
        end do
        nm = cnt
    end if
end subroutine BuildJumpPoints


!> cs.ev <- seq(1,csl,by=jump); if last != csl append csl; then prepend 0
!> and append csl+1 - the cycle-subseries evaluation grid (one point of
!> extrapolation on each side of the subseries).
subroutine BuildCsEv(csl, jump, m, nm)
    integer, intent(in) :: csl, jump
    real(kind = dbl), allocatable, intent(out) :: m(:)
    integer, intent(out) :: nm
    real(kind = dbl), allocatable :: inner(:)
    integer :: nm_inner, i

    call BuildJumpPoints(csl, jump, inner, nm_inner)
    nm = nm_inner + 2
    allocate(m(nm))
    m(1) = 0d0
    do i = 1, nm_inner
        m(i+1) = inner(i)
    end do
    m(nm) = dble(csl+1)
    deallocate(inner)
end subroutine BuildCsEv


!***************************************************************************
!> Direct port of stlplus's c_loess (src/loess.cpp): local polynomial
!> (degree 0/1/2) regression with a tricube weight kernel.
!>
!> x,y: data restricted to valid points (length n). w: weights (length n).
!> m(nm): integer target points (as reals). l_idx0(nm): 0-based left window
!> start index into x/y. max_dist(nm): per-point max neighbour distance.
!***************************************************************************
subroutine LoessFit(x, y, n, degree, span_in, w, m, nm, l_idx0, max_dist, result, slopes)
    integer, intent(in) :: n, degree, span_in, nm
    real(kind = dbl), intent(in) :: x(n), y(n), w(n), m(nm), max_dist(nm)
    integer, intent(in) :: l_idx0(nm)
    real(kind = dbl), intent(out) :: result(nm), slopes(nm)
    integer :: span, i, j, base
    real(kind = dbl) :: xj(span_in), wj(span_in), xwj(span_in), x2wj(span_in), x3wj(span_in)
    real(kind = dbl) :: r, tmp1, tmp2
    real(kind = dbl) :: a, b, c, d, e, a1, b1, c1, a2, b2, c2, det

    span = span_in
    if (span > n) span = n

    result = 0d0
    slopes = 0d0

    do i = 1, nm
        base = l_idx0(i)  !< 0-based
        a = 0d0
        do j = 1, span
            xj(j) = x(base + j) - m(i)
            r = dabs(xj(j))
            tmp1 = r / max_dist(i)
            tmp2 = 1d0 - tmp1*tmp1*tmp1
            wj(j) = tmp2*tmp2*tmp2
            wj(j) = wj(j) * w(base + j)
            a = a + wj(j)
        end do

        if (degree == 0) then
            a1 = 1d0 / a
            do j = 1, span
                result(i) = result(i) + wj(j)*a1*y(base+j)
            end do
        else
            b = 0d0; c = 0d0
            do j = 1, span
                xwj(j) = xj(j)*wj(j)
                x2wj(j) = xj(j)*xwj(j)
                b = b + xwj(j)
                c = c + x2wj(j)
            end do
            if (degree == 1) then
                det = 1d0 / (a*c - b*b)
                a1 = c*det; b1 = -b*det; c1 = a*det
                do j = 1, span
                    result(i) = result(i) + (wj(j)*a1 + xwj(j)*b1)*y(base+j)
                    slopes(i) = slopes(i) + (wj(j)*b1 + xwj(j)*c1)*y(base+j)
                end do
            else
                d = 0d0; e = 0d0
                do j = 1, span
                    x3wj(j) = xj(j)*x2wj(j)
                    d = d + x3wj(j)
                    e = e + x3wj(j)*xj(j)
                end do
                a1 = e*c - d*d
                b1 = c*d - e*b
                c1 = b*d - c*c
                a2 = c*d - e*b
                b2 = e*a - c*c
                c2 = b*c - d*a
                det = 1d0 / (a*a1 + b*b1 + c*c1)
                a1=a1*det; b1=b1*det; c1=c1*det
                a2=a2*det; b2=b2*det; c2=c2*det
                do j = 1, span
                    result(i) = result(i) + (wj(j)*a1 + xwj(j)*b1 + x2wj(j)*c1)*y(base+j)
                    slopes(i) = slopes(i) + (wj(j)*a2 + xwj(j)*b2 + x2wj(j)*c2)*y(base+j)
                end do
            end if
        end if
    end do
end subroutine LoessFit


!> Direct port of stlplus's c_ma (src/ma.cpp): two width-np_ moving
!> averages, then one width-3 moving average - the STL low-pass filter.
subroutine LowPassFilter3(x, n, np_, ans, nans)
    integer, intent(in) :: n, np_
    real(kind = dbl), intent(in) :: x(n)
    integer, intent(out) :: nans
    real(kind = dbl), intent(out) :: ans(n - 2*np_)
    integer :: nn, i
    real(kind = dbl) :: ma_tmp
    real(kind = dbl), allocatable :: ma(:), ma2(:)

    nn = n - np_*2
    nans = nn
    allocate(ma(nn+np_+1), ma2(nn+2))

    ma_tmp = 0d0
    do i = 1, np_
        ma_tmp = ma_tmp + x(i)
    end do
    ma(1) = ma_tmp/dble(np_)
    do i = np_+1, nn+2*np_
        ma_tmp = ma_tmp - x(i-np_) + x(i)
        ma(i-np_+1) = ma_tmp/dble(np_)
    end do

    ma_tmp = 0d0
    do i = 1, np_
        ma_tmp = ma_tmp + ma(i)
    end do
    ma2(1) = ma_tmp/dble(np_)
    do i = np_+1, nn+np_+1
        ma_tmp = ma_tmp - ma(i-np_) + ma(i)
        ma2(i-np_+1) = ma_tmp/dble(np_)
    end do

    ma_tmp = 0d0
    do i = 1, 3
        ma_tmp = ma_tmp + ma2(i)
    end do
    ans(1) = ma_tmp/3d0
    do i = 4, nn+2
        ma_tmp = ma_tmp - ma2(i-3) + ma2(i)
        ans(i-2) = ma_tmp/3d0
    end do

    deallocate(ma, ma2)
end subroutine LowPassFilter3


!> Direct port of stlplus's c_interp (src/interp.cpp): cubic Hermite
!> interpolation between the sparse loess evaluation points, using both
!> the fitted value and slope at each.
subroutine HermiteInterp(m, fits, slopes, nm, at, nat, ans)
    integer, intent(in) :: nm, nat
    real(kind = dbl), intent(in) :: m(nm), fits(nm), slopes(nm), at(nat)
    real(kind = dbl), intent(out) :: ans(nat)
    integer :: i, j
    real(kind = dbl) :: u, h, u2, u3

    j = 1  !< 1-based "leftmost vertex" index; R's j starts 0-based at 0
    do i = 1, nat
        if (at(i) > m(j+1)) j = j + 1
        h = m(j+1) - m(j)
        u = (at(i) - m(j)) / h
        u2 = u*u
        u3 = u2*u
        ans(i) = (2d0*u3 - 3d0*u2 + 1d0)*fits(j) + (3d0*u2 - 2d0*u3)*fits(j+1) + &
                 (u3 - 2d0*u2 + u)*slopes(j)*h + (u3 - u2)*slopes(j+1)*h
    end do
end subroutine HermiteInterp


!> R's default quantile() algorithm (Hyndman & Fan 1996, type 7).
double precision function QuantileType7(x_sorted, n, p) result(q)
    integer, intent(in) :: n
    real(kind = dbl), intent(in) :: x_sorted(n), p
    real(kind = dbl) :: h
    integer :: lo, hi

    h = dble(n-1)*p + 1d0
    lo = int(h)
    if (lo < 1) lo = 1
    if (lo >= n) then
        q = x_sorted(n)
        return
    end if
    hi = lo + 1
    q = x_sorted(lo) + (h - dble(lo)) * (x_sorted(hi) - x_sorted(lo))
end function QuantileType7

end module m_flux_despike_core
