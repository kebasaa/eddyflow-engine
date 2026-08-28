!***************************************************************************
! m_storage_clean_core.f90
! -------------------------
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
! \brief       Whole-run storage-flux cleaning: RFlux's cleanFlux() storage
!              branch (Vitale et al. 2020, Biogeosciences) - a Tukey boxplot
!              (range = 3, "far out" fence) outlier test on the storage term,
!              binned by time-of-day across the whole run, followed by
!              linear interpolation over interior gaps.
!
! \details     Everything here is a function of its arguments: no project
!              settings, no column table, no logging - the same discipline
!              m_pwb_core.f90 and m_flux_despike_core.f90 follow, for the
!              same reason: it lets a reference driver call this without
!              linking the rest of the engine.
!
!              RFlux's own steps (cleanFlux.R):
!              1. Clip physically implausible values (its own hard-coded
!                 CO2 range, > 70 or < -100 umol m-2 s-1) to NA before
!                 anything else sees them. Not ported: EddyFlow supports
!                 gases RFlux does not, and a hard-coded CO2 range would be
!                 silently wrong for any other one. If a project wants a
!                 range check on its storage term, EddyFlow's own absolute
!                 limits machinery is the place for it, not this pass.
!              2. hod <- gl(48,1,N): the time-of-day class is a FIXED 48
!                 half-hourly slots per day, an assumption that only holds
!                 for exactly-30-minute averaging. Ported generalised
!                 instead, over mfreq classes derived from the run's own
!                 modal period spacing (the caller resolves this the same
!                 way PostProcessFluxDespiking does), so the same code
!                 handles any regular averaging period.
!              3. Boxplot(..., range=3): Tukey's "far out" fence, 3 times
!                 the interquartile range beyond each quartile. car::Boxplot
!                 delegates to grDevices::boxplot.stats, which takes its
!                 hinges from stats::fivenum() - Tukey's own hinge
!                 definition, NOT quantile()'s default type-7 interpolation
!                 (confirmed against R: type 7 matched 479/480 generated
!                 test values and missed exactly the one case where the two
!                 definitions diverge). FivenumHinge below is fivenum()'s
!                 own d <- c(1, n4, (n+1)/2, n+1-n4, n); 0.5*(x[floor(d)] +
!                 x[ceiling(d)]) formula, ported rather than approximated.
!              4. na.approx(..., na.rm = FALSE): linear interpolation over
!                 gaps that have a valid point on both sides; a leading or
!                 trailing gap is left as `missing`, matching zoo's default
!                 rule (no extrapolation past the data).
!              5. A period `missing` in the input is never fabricated by
!                 interpolation into a value - it comes back `missing`
!                 regardless of what surrounds it, matching cleanFlux.R's
!                 own final `replace(..., which(is.na(original)), NA)`.
!
! \author      Jonathan Muller, ETH Zurich
! \note
! \sa          RFlux-master/R/cleanFlux.R (the storage==TRUE branch),
!              m_flux_despike_core.f90 (the sibling whole-run post-pass)
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
module m_storage_clean_core
    use m_numeric_kinds
    implicit none
    private
    public :: CleanStorageSeries

contains

!***************************************************************************
!> x(n): the storage series, one value per period, in chronological order.
!> `missing` where a period has no value. hod(n): each period's time-of-day
!> class, 1..mfreq (the caller resolves this from the actual timestamps, not
!> from array position, so a gap in the run cannot shift the phase).
!>
!> spike_flag(n) comes back 1 where the value was outside its class's Tukey
!> fence, 0 elsewhere (including every `missing` position). cleaned(n) is
!> x with spike positions linearly interpolated where an interior neighbour
!> pair allows it, and `missing` everywhere x itself was `missing` -
!> regardless of whether interpolation could otherwise have filled it.
!***************************************************************************
subroutine CleanStorageSeries(x, n, hod, mfreq, missing, spike_flag, cleaned)
    integer, intent(in) :: n, mfreq
    real(kind = dbl), intent(in) :: x(n)
    integer, intent(in) :: hod(n)
    real(kind = dbl), intent(in) :: missing
    integer, intent(out) :: spike_flag(n)
    real(kind = dbl), intent(out) :: cleaned(n)

    real(kind = dbl), allocatable :: class_vals(:)
    integer :: c, i, nv, left, right
    real(kind = dbl) :: q1, q3, iqr, lo_fence, hi_fence
    logical :: gap(n)

    spike_flag = 0
    cleaned = x

    !> One Tukey fence per time-of-day class, from every non-missing value
    !> the whole run has in that class.
    allocate(class_vals(n))
    do c = 1, mfreq
        nv = 0
        do i = 1, n
            if (hod(i) == c .and. x(i) /= missing) then
                nv = nv + 1
                class_vals(nv) = x(i)
            end if
        end do
        if (nv < 4) cycle  !> Too few points for quartiles to mean anything.

        call HPSORT(nv, class_vals(1:nv))
        q1 = FivenumHinge(class_vals(1:nv), nv, .false.)
        q3 = FivenumHinge(class_vals(1:nv), nv, .true.)
        iqr = q3 - q1
        lo_fence = q1 - 3d0 * iqr
        hi_fence = q3 + 3d0 * iqr

        do i = 1, n
            if (hod(i) /= c .or. x(i) == missing) cycle
            if (x(i) < lo_fence .or. x(i) > hi_fence) spike_flag(i) = 1
        end do
    end do
    deallocate(class_vals)

    !> Linear interpolation over every spike or already-missing position,
    !> restricted to interior gaps (a valid point on both sides) - matching
    !> na.approx(..., na.rm = FALSE)'s default of not extrapolating past
    !> the data.
    gap = (spike_flag == 1) .or. (x == missing)
    i = 1
    do while (i <= n)
        if (.not. gap(i)) then
            i = i + 1
            cycle
        end if
        left = i - 1
        right = i
        do while (right <= n .and. gap(right))
            right = right + 1
        end do
        !> [left+1, right-1] is one contiguous gap run. Fill it only if it
        !> has a valid point on both sides.
        if (left >= 1 .and. right <= n) then
            do i = left + 1, right - 1
                cleaned(i) = x(left) + (x(right) - x(left)) &
                    * dble(i - left) / dble(right - left)
            end do
        end if
        i = right
    end do

    !> A genuinely missing period is never fabricated into a value, no
    !> matter what interpolation would have produced for it.
    where (x == missing) cleaned = missing

contains

    !> Tukey's upper (upper=.true.) or lower hinge, stats::fivenum()'s own
    !> d <- c(1, n4, (n+1)/2, n+1-n4, n); 0.5*(x[floor(d)] + x[ceiling(d)]).
    !> x_sorted must already be ascending.
    double precision function FivenumHinge(x_sorted, nn, upper) result(q)
        integer, intent(in) :: nn
        real(kind = dbl), intent(in) :: x_sorted(nn)
        logical, intent(in) :: upper
        real(kind = dbl) :: n4, d
        integer :: dlo, dhi

        n4 = dble(floor((dble(nn) + 3d0) / 2d0)) / 2d0
        if (upper) then
            d = dble(nn) + 1d0 - n4
        else
            d = n4
        end if
        dlo = floor(d)
        dhi = ceiling(d)
        q = 0.5d0 * (x_sorted(dlo) + x_sorted(dhi))
    end function FivenumHinge

end subroutine CleanStorageSeries

end module m_storage_clean_core
