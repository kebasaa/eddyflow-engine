!***************************************************************************
! langford_noise.f90
! ------------------
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
! \brief       Random flux error and limit of detection after Langford et
!              al. (2015), Sect. 2.2.
! \author      Jonathan Muller
!
! \note
! Langford, B., Acton, W., Ammann, C., Valach, A., and Nemitz, E.: Eddy-
! covariance data with low signal-to-noise ratio: time-lag determination,
! uncertainties and limit of detection, Atmos. Meas. Tech., 8, 4197-4213,
! 2015.
!
! Far from the transport lag, w and the scalar are no longer correlated, so
! whatever the cross-covariance function f(D) does out there is noise:
! instrument noise and the chance correlation of a finite record. Its
! root mean square over that region, RE_RMSE = sqrt(sigma^2 + mu^2) (their
! Eq. 9), is a random error of the covariance; its standard deviation alone
! is RE_SD. LoD = alpha * RE, with alpha 1.96 for 95 %.
!
! One routine with two callers, so that both measure the same thing:
!
!  - the random error method langford_15 (RandomUncertaintyHandle), on the
!    detrended series the flux is computed from, centred on lag zero since
!    those are already aligned;
!  - the detection limit langford_15 (FluxDetectionLimit), on the raw series
!    before the lag is applied, centred on each gas's lag, because time-lag
!    borrowing must test against it before the shift.
!
! The windows, RUsetup%lod_window:
!   0  +/-150..180 s, Langford et al. (2015), Sect. 2.2
!   1  +/-160..180 s, Maier et al. (2022), after Spirig et al. (2005)
!   2  automatic. Spirig et al. (2005) chose their window because it lay
!      "several times the integral timescale" away; neither they nor
!      Langford give a rule. Ours: start at lod_its_mult times the integral
!      turbulence scale, or 10 s beyond the variable's time-lag search
!      half-width, whichever is later, and run 30 s. Falls back to window 0
!      when the scale cannot be measured.
!
! \sa          m_covmax_core.f90 (FarLagCovarianceNoise, IntegralScaleFromCorr)
!***************************************************************************
subroutine LangfordNoise(Set, nrow, ncol, centre, want, re)
    use m_rp_global_var
    use m_covmax_core, only: FarLagCovarianceNoise, CrossCovarianceSeries, &
        IntegralScaleFromCorr
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow
    integer, intent(in) :: ncol
    real(kind = dbl), intent(in) :: Set(nrow, ncol)
    !> Row lag each variable's windows are placed about.
    integer, intent(in) :: centre(ncol)
    !> Which variables to measure; the rest are left at error.
    logical, intent(in) :: want(ncol)
    real(kind = dbl), intent(out) :: re(ncol)
    !> local variables
    integer :: var
    integer :: lo
    integer :: hi
    integer :: nlags
    integer :: lagmax
    integer :: half_rows
    real(kind = dbl) :: f
    real(kind = dbl) :: itscale
    real(kind = dbl), allocatable :: corr(:)
    character(256) :: msg

    re = error
    f = Metadata%ac_freq
    if (f <= 0d0) return

    do var = u, lastGas
        if (var == v .or. var == w) cycle
        if (var > ncol) exit
        if (.not. want(var)) cycle
        if (.not. E2Col(var)%present) cycle

        select case (RUsetup%lod_window)
            case (1)
                lo = nint(160d0 * f)
                hi = nint(180d0 * f)
            case (2)
                !> The integral turbulence itscale, by the run's own rule, from
                !> the cross-covariance on [centre, centre + ru_tlag_max]
                !> normalised to one at the centre - the function
                !> IntegralTurbulenceScale integrates, placed at this lag.
                itscale = error
                lagmax = nint(RUsetup%tlag_max * f)
                if (lagmax > 0 .and. abs(centre(var)) + lagmax < nrow) then
                    allocate(corr(0:lagmax))
                    call CrossCovarianceSeries(Set(:, w), Set(:, var), nrow, &
                        centre(var), centre(var) + lagmax, error, corr)
                    if (corr(0) /= error .and. corr(0) /= 0d0 &
                        .and. all(corr /= error)) then
                        corr = corr / corr(0)
                        itscale = IntegralScaleFromCorr(corr, lagmax, 1d0 / f, &
                            RUsetup%its_meth)
                    end if
                    deallocate(corr)
                end if
                !> The same sanity bound IntegralTurbulenceScale applies.
                if (itscale /= error .and. itscale > 0d0 &
                    .and. itscale <= 2d0 * RUsetup%tlag_max) then
                    if (var >= ts) then
                        half_rows = nint(0.5d0 * abs(E2Col(var)%max_tl &
                            - E2Col(var)%min_tl) * f)
                    else
                        half_rows = 0
                    end if
                    lo = max(nint(RUsetup%lod_its_mult * itscale * f), &
                        half_rows + nint(10d0 * f))
                    hi = lo + nint(30d0 * f)
                    write(msg, '(a, a, a, f7.1, a, f7.1, a)') &
                        '   Langford window (automatic) for ', &
                        trim(E2Col(var)%label), ': ', dble(lo) / f, &
                        ' to ', dble(hi) / f, ' s'
                    call LogSay(trim(msg))
                else
                    lo = nint(150d0 * f)
                    hi = nint(180d0 * f)
                end if
            case default
                lo = nint(150d0 * f)
                hi = nint(180d0 * f)
        end select

        !> Keep at least half the record behind every lag, so the far end of
        !> the window is not a covariance of a handful of samples. A window
        !> that no longer fits is no estimate at all.
        hi = min(hi, nrow / 2 - abs(centre(var)))
        if (hi - lo < 2) cycle

        call FarLagCovarianceNoise(Set(:, w), Set(:, var), nrow, centre(var), &
            lo, hi, error, RUsetup%lod_stat, re(var), nlags)
    end do
end subroutine LangfordNoise
