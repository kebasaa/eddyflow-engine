!***************************************************************************
! m_driftfix_core.f90
! -------------------
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
! \brief       Clock-drift estimate from a time series of time lags, after
!              innFLUX's driftfix tool (innFLUX v1.2.0, 2023, University of
!              Innsbruck, unpublished; tools/driftfix/innFLUX_driftfix.m,
!              based on a routine by M. Graus, adapted by M. Striednig). It is
!              not described in Striednig et al. (2020).
!
!              A data logger's clock that drifts against the analyser's
!              makes the time lag drift; a clock that is reset makes it jump.
!              driftfix smooths the per-period lags of a few well-resolved
!              tracers with a robust LOWESS (Cleveland 1979; MATLAB
!              smooth(..., 'rlowess')) over one day, finds the jumps where
!              that one-day curve rises faster than a threshold, and smooths
!              again, over several days, inside each stretch between jumps.
!              The result is the drift: the lag, as a function of time, that
!              the clock adds.
!
!              RobustLowess follows MATLAB's definition: for each point a
!              straight line is fitted to the span nearest points in index
!              order, weighted by the tricube of the distance in x scaled by
!              the largest distance in that window, then refitted five times
!              with bisquare robustness weights on the residuals divided by
!              six times their median absolute value.
!              DetectJumps stands in for MATLAB's risetime, which innFLUX
!              uses to place each jump: here the jump is at the midpoint
!              of the first 10 % and 90 % crossings of the smoothed lag
!              within two days of where the rise began.
!
! \note        Pure, and a function of its arguments only, so it can be checked
!              outside the engine (src_tools/innflux_compare_main.f90).
! \author      Jonathan Muller
! \sa          m_cumcov.f90
!***************************************************************************
module m_driftfix_core
    use m_numeric_kinds
    implicit none
    private

    public :: RobustLowess
    public :: InterpLinear
    public :: DetectJumps
    public :: SegmentedDrift
    public :: DaysFromCivil

contains

    !***************************************************************************
    !> \brief Days since 1970-01-01 of a proleptic Gregorian date, plus the
    !>        fraction of the day. Howard Hinnant's days_from_civil.
    !***************************************************************************
    pure function DaysFromCivil(y, m, d, hh, mm) result(days)
        integer, intent(in) :: y
        integer, intent(in) :: m
        integer, intent(in) :: d
        integer, intent(in) :: hh
        integer, intent(in) :: mm
        real(kind = dbl) :: days
        integer :: yy
        integer :: era
        integer :: yoe
        integer :: doy
        integer :: doe
        integer :: mp

        yy = y
        if (m <= 2) yy = y - 1
        if (yy >= 0) then
            era = yy / 400
        else
            era = (yy - 399) / 400
        end if
        yoe = yy - era * 400
        mp = mod(m + 9, 12)
        doy = (153 * mp + 2) / 5 + d - 1
        doe = yoe * 365 + yoe / 4 - yoe / 100 + doy
        days = dble(era * 146097 + doe - 719468) + (dble(hh) * 60d0 + dble(mm)) / 1440d0
    end function DaysFromCivil

    !***************************************************************************
    !> \brief Robust locally weighted linear regression, MATLAB smooth(x, y,
    !>        span, 'rlowess'). x ascending. span points per window, odd
    !>        or even as given, clipped to [2, n].
    !***************************************************************************
    pure subroutine RobustLowess(x, y, n, span, niter, ys)
        integer, intent(in) :: n
        integer, intent(in) :: span
        integer, intent(in) :: niter
        real(kind = dbl), intent(in) :: x(n)
        real(kind = dbl), intent(in) :: y(n)
        real(kind = dbl), intent(out) :: ys(n)
        real(kind = dbl) :: rw(n)
        real(kind = dbl) :: res(n)
        real(kind = dbl) :: absres(n)
        real(kind = dbl) :: mad
        integer :: it
        integer :: i
        integer :: k

        ys = y
        if (n < 3) return
        k = max(2, min(span, n))
        rw = 1d0
        do it = 0, niter
            do i = 1, n
                ys(i) = LocalFit(i)
            end do
            if (it == niter) exit
            res = y - ys
            absres = abs(res)
            mad = Median(absres, n)
            if (mad <= 0d0) exit
            do i = 1, n
                rw(i) = abs(res(i)) / (6d0 * mad)
                if (rw(i) < 1d0) then
                    rw(i) = (1d0 - rw(i) ** 2) ** 2
                else
                    rw(i) = 0d0
                end if
            end do
        end do

    contains

        pure function LocalFit(i) result(fit)
            integer, intent(in) :: i
            real(kind = dbl) :: fit
            integer :: lo
            integer :: hi
            integer :: j
            real(kind = dbl) :: dmax
            real(kind = dbl) :: wj
            real(kind = dbl) :: sw
            real(kind = dbl) :: sx
            real(kind = dbl) :: sy
            real(kind = dbl) :: sxx
            real(kind = dbl) :: sxy
            real(kind = dbl) :: den

            !> The span nearest points in index order, centred on i where the
            !> ends allow and pushed inwards where they do not.
            lo = i - (k - 1) / 2
            hi = lo + k - 1
            if (lo < 1) then
                lo = 1
                hi = k
            end if
            if (hi > n) then
                hi = n
                lo = n - k + 1
            end if
            dmax = max(abs(x(i) - x(lo)), abs(x(hi) - x(i)))
            sw = 0d0
            sx = 0d0
            sy = 0d0
            sxx = 0d0
            sxy = 0d0
            do j = lo, hi
                if (dmax > 0d0) then
                    wj = abs(x(j) - x(i)) / (dmax * 1.0000001d0)
                    wj = (1d0 - wj ** 3) ** 3
                else
                    wj = 1d0
                end if
                wj = wj * rw(j)
                sw = sw + wj
                sx = sx + wj * x(j)
                sy = sy + wj * y(j)
                sxx = sxx + wj * x(j) ** 2
                sxy = sxy + wj * x(j) * y(j)
            end do
            if (sw <= 0d0) then
                fit = y(i)
                return
            end if
            den = sw * sxx - sx ** 2
            if (abs(den) <= 1d-12 * max(1d0, abs(sw * sxx))) then
                fit = sy / sw
            else
                fit = (sy * sxx - sx * sxy) / den + (sw * sxy - sx * sy) / den * x(i)
            end if
        end function LocalFit
    end subroutine RobustLowess

    !> Median of a(1:n), by a copy and insertion sort: n is a few thousand.
    pure function Median(a, n) result(med)
        integer, intent(in) :: n
        real(kind = dbl), intent(in) :: a(n)
        real(kind = dbl) :: med
        real(kind = dbl) :: b(n)
        real(kind = dbl) :: t
        integer :: i
        integer :: j

        b = a
        do i = 2, n
            t = b(i)
            j = i - 1
            do while (j >= 1)
                if (b(j) <= t) exit
                b(j + 1) = b(j)
                j = j - 1
            end do
            b(j + 1) = t
        end do
        if (mod(n, 2) == 1) then
            med = b((n + 1) / 2)
        else
            med = 0.5d0 * (b(n / 2) + b(n / 2 + 1))
        end if
    end function Median

    !***************************************************************************
    !> \brief y at xq by linear interpolation in (x, y), x ascending; constant
    !>        beyond the ends.
    !***************************************************************************
    pure function InterpLinear(x, y, n, xq) result(yq)
        integer, intent(in) :: n
        real(kind = dbl), intent(in) :: x(n)
        real(kind = dbl), intent(in) :: y(n)
        real(kind = dbl), intent(in) :: xq
        real(kind = dbl) :: yq
        integer :: lo
        integer :: hi
        integer :: mid

        if (n == 1 .or. xq <= x(1)) then
            yq = y(1)
            return
        end if
        if (xq >= x(n)) then
            yq = y(n)
            return
        end if
        lo = 1
        hi = n
        do while (hi - lo > 1)
            mid = (lo + hi) / 2
            if (x(mid) <= xq) then
                lo = mid
            else
                hi = mid
            end if
        end do
        if (x(hi) == x(lo)) then
            yq = y(lo)
        else
            yq = y(lo) + (y(hi) - y(lo)) * (xq - x(lo)) / (x(hi) - x(lo))
        end if
    end function InterpLinear

    !***************************************************************************
    !> \brief Where the one-day smoothed lag rises faster than thr allows.
    !>
    !> t (days, ascending) and y1d (the one-day robust LOWESS) on n periods.
    !> A period is in a rise when y1d(t) - y1d(t + 1 day) < thr, innFLUX's
    !> test (thr negative, in s per day). Each run of such periods is a jump,
    !> placed at the midpoint of the first 10 % and 90 % crossings of y1d
    !> within two days of the run's start. ok is false when the runs do not
    !> close inside the series, which innFLUX treats as "no jumps".
    !***************************************************************************
    pure subroutine DetectJumps(t, y1d, n, thr, njump, jump_at, ok)
        integer, intent(in) :: n
        real(kind = dbl), intent(in) :: t(n)
        real(kind = dbl), intent(in) :: y1d(n)
        real(kind = dbl), intent(in) :: thr
        integer, intent(out) :: njump
        integer, intent(out) :: jump_at(n)
        logical, intent(out) :: ok
        logical :: rising(n)
        integer :: i
        integer :: s
        integer :: e
        integer :: lt
        integer :: ut
        real(kind = dbl) :: ymin
        real(kind = dbl) :: ymax
        real(kind = dbl) :: lvl10
        real(kind = dbl) :: lvl90

        njump = 0
        jump_at = 0
        ok = .true.
        do i = 1, n
            if (t(i) + 1d0 > t(n)) then
                rising(i) = .false.
            else
                rising(i) = y1d(i) - InterpLinear(t, y1d, n, t(i) + 1d0) < thr
            end if
        end do
        if (n < 2) return
        !> A run open at either end has no start or no stop, innFLUX's
        !> unequal counts.
        if (rising(1) .or. rising(n)) then
            ok = .false.
            return
        end if
        i = 1
        do while (i < n)
            if (.not. rising(i) .and. rising(i + 1)) then
                s = i + 1
                e = s
                do while (e < n)
                    if (.not. rising(e + 1)) exit
                    e = e + 1
                end do
                !> The rise, within two days of its start.
                e = s
                do while (e < n)
                    if (t(e + 1) > t(s) + 2d0) exit
                    e = e + 1
                end do
                ymin = minval(y1d(s:e))
                ymax = maxval(y1d(s:e))
                lvl10 = ymin + 0.1d0 * (ymax - ymin)
                lvl90 = ymin + 0.9d0 * (ymax - ymin)
                lt = s + minloc(y1d(s:e), 1) - 1
                do while (lt < e)
                    if (y1d(lt) >= lvl10) exit
                    lt = lt + 1
                end do
                ut = lt
                do while (ut < e)
                    if (y1d(ut) >= lvl90) exit
                    ut = ut + 1
                end do
                njump = njump + 1
                jump_at(njump) = (lt + ut + 1) / 2
                !> On past this run.
                do while (i < n)
                    if (.not. rising(i + 1)) exit
                    i = i + 1
                end do
            end if
            i = i + 1
        end do
    end subroutine DetectJumps

    !***************************************************************************
    !> \brief The drift at each of the n periods t (days, ascending).
    !>
    !> Points (px, py), npt of them, px ascending, are every reference gas's
    !> lag in every good period. They are smoothed over one day
    !> (span1d points) to find jumps, then over span points within each
    !> stretch between jumps, or over all of them when the jumps cannot be
    !> placed. drift is error where no point exists at all.
    !***************************************************************************
    pure subroutine SegmentedDrift(px, py, npt, t, n, span1d, span, thr, err, &
        drift, njump, jumps_ok)
        integer, intent(in) :: npt
        integer, intent(in) :: n
        integer, intent(in) :: span1d
        integer, intent(in) :: span
        real(kind = dbl), intent(in) :: px(npt)
        real(kind = dbl), intent(in) :: py(npt)
        real(kind = dbl), intent(in) :: t(n)
        real(kind = dbl), intent(in) :: thr
        real(kind = dbl), intent(in) :: err
        real(kind = dbl), intent(out) :: drift(n)
        integer, intent(out) :: njump
        logical, intent(out) :: jumps_ok
        real(kind = dbl) :: s1d(npt)
        real(kind = dbl) :: y1d(n)
        real(kind = dbl) :: sw(npt)
        integer :: jump_at(n)
        integer :: seg_lo
        integer :: seg_hi
        integer :: k
        integer :: i
        integer :: p_lo
        integer :: p_hi
        real(kind = dbl) :: t_lo
        real(kind = dbl) :: t_hi

        drift = err
        njump = 0
        jumps_ok = .false.
        if (npt < 1 .or. n < 1) return
        if (npt < 3) then
            drift = sum(py) / dble(npt)
            return
        end if

        call RobustLowess(px, py, npt, span1d, 5, s1d)
        do i = 1, n
            y1d(i) = InterpLinear(px, s1d, npt, t(i))
        end do
        call DetectJumps(t, y1d, n, thr, njump, jump_at, jumps_ok)
        if (.not. jumps_ok) njump = 0

        !> Smooth each stretch between jumps on its own.
        seg_lo = 1
        do k = 1, njump + 1
            if (k <= njump) then
                seg_hi = jump_at(k) - 1
            else
                seg_hi = n
            end if
            if (seg_hi < seg_lo) cycle
            t_lo = t(seg_lo)
            if (k <= njump) then
                t_hi = t(jump_at(k))
            else
                t_hi = huge(1d0)
            end if
            if (k == 1) t_lo = -huge(1d0)
            !> The points inside this stretch.
            p_lo = 0
            p_hi = -1
            do i = 1, npt
                if (px(i) >= t_lo .and. px(i) < t_hi) then
                    if (p_lo == 0) p_lo = i
                    p_hi = i
                end if
            end do
            if (p_lo == 0) then
                seg_lo = seg_hi + 1
                cycle
            end if
            if (p_hi - p_lo + 1 >= 3) then
                call RobustLowess(px(p_lo:p_hi), py(p_lo:p_hi), p_hi - p_lo + 1, &
                    span, 5, sw(p_lo:p_hi))
            else
                sw(p_lo:p_hi) = py(p_lo:p_hi)
            end if
            do i = seg_lo, seg_hi
                drift(i) = InterpLinear(px(p_lo:p_hi), sw(p_lo:p_hi), &
                    p_hi - p_lo + 1, t(i))
            end do
            seg_lo = seg_hi + 1
        end do
    end subroutine SegmentedDrift

end module m_driftfix_core
