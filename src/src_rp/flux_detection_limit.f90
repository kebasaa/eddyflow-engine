!***************************************************************************
! flux_detection_limit.f90
! ------------------------
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
! \brief       Flux detection limit after Wienhold et al. (1994).
! \author      Jonathan Muller
!
! \note
! The cross-covariance function of w with a scalar carries the flux in a
! peak near the transport lag and nothing but noise far away from it. The
! scatter of the function out there is therefore a noise floor on the
! covariance: a flux smaller than it cannot be distinguished from zero.
!
! Wienhold et al. (1994), "Measurements of N2O fluxes from fertilized
! grassland using a fast response tunable diode laser spectrometer",
! measure that scatter in two windows placed symmetrically either side of
! the peak and average them. This routine does the same, taking the peak to
! be the lag the run actually used for that gas.
!
! The value is a standard deviation of a covariance, and is reported in
! covariance units - the same units as the covariance it qualifies, before
! any conversion to a flux. Comparing a flux against it therefore means
! comparing like with like only at level 0; a corrected flux has been scaled
! and the limit has not.
!
! Two windows, not one: an offset large enough to clear the peak on the late
! side can still be inside it on the early side if the cross-covariance is
! asymmetric, and averaging the two is what Wienhold does. Where only one
! side yields a value - the early window can fall off the start of a short
! record - that side is used alone rather than the period losing its limit.
!
! \sa          EC_Software_Preproc/EddyUH_detlim_Preproc.m in EddyUH 1.7b,
!              which is the same method with a 50 s window at +/-100 s.
!***************************************************************************
subroutine FluxDetectionLimit(Set, nrow, ncol)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow
    integer, intent(in) :: ncol
    real(kind = dbl), intent(in) :: Set(nrow, ncol)
    !> local variables
    integer :: gas
    integer :: side
    integer :: nside
    integer :: offset_rows
    integer :: half_rows
    real(kind = dbl) :: ColW(nrow)
    real(kind = dbl) :: ColGas(nrow)
    real(kind = dbl) :: side_sd
    real(kind = dbl) :: acc

    Essentials%detlim = error
    if (RPSetup%detlim_meth /= 'wienhold_94') return
    if (Metadata%ac_freq <= 0d0) return

    !> Window geometry on the row grid. Both are at least one row, so a
    !> setting finer than the sample interval degrades to the smallest
    !> window that exists rather than to an empty one.
    offset_rows = max(1, nint(RPSetup%detlim_offset_s * Metadata%ac_freq))
    half_rows   = max(1, nint(RPSetup%detlim_window_s * Metadata%ac_freq / 2d0))

    ColW(1:nrow) = Set(1:nrow, w)

    do gas = firstGas, lastGas
        if (.not. E2Col(gas)%present) cycle
        !> RowLags is the lag this period settled on for this gas, whatever
        !> method chose it, and it may be zero or negative - zero when the
        !> run compensates no lag at all, negative for an open path. The
        !> windows are placed relative to it either way, so there is nothing
        !> to reject here.
        ColGas(1:nrow) = Set(1:nrow, gas)

        acc = 0d0
        nside = 0
        do side = -1, 1, 2
            call CrossCovarianceScatter(ColW, ColGas, nrow, &
                RowLags(gas) + side * offset_rows, half_rows, side_sd)
            if (side_sd /= error) then
                acc = acc + side_sd
                nside = nside + 1
            end if
        end do

        if (nside > 0) Essentials%detlim(gas) = acc / dble(nside)
    end do
end subroutine FluxDetectionLimit

!***************************************************************************
!
! \brief       Standard deviation of the cross-covariance function over a
!              window of lags centred on centre_lag.
! \author      Jonathan Muller
!
! \note
! The window is evaluated lag by lag with CovarianceW, which is the same
! estimator the rest of the time-lag code uses, so the scatter measured here
! is the scatter of the quantity the flux is actually taken from.
!
! Lags for which the covariance cannot be formed are skipped rather than
! counted as zero; a window that keeps fewer than three of them is refused,
! because the standard deviation of one or two points says nothing about the
! spread of the function.
!***************************************************************************
subroutine CrossCovarianceScatter(ColW, ColGas, nrow, centre_lag, half_rows, sd)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow
    integer, intent(in) :: centre_lag
    integer, intent(in) :: half_rows
    real(kind = dbl), intent(in) :: ColW(nrow)
    real(kind = dbl), intent(in) :: ColGas(nrow)
    real(kind = dbl), intent(out) :: sd
    !> local variables
    integer, parameter :: min_lags = 3
    integer :: lag
    integer :: n
    real(kind = dbl) :: cov
    real(kind = dbl) :: s1
    real(kind = dbl) :: s2
    real(kind = dbl) :: mean

    sd = error
    n = 0
    s1 = 0d0
    s2 = 0d0

    do lag = centre_lag - half_rows, centre_lag + half_rows
        !> A window running off either end of the record leaves nothing to
        !> correlate. Asking for it is not an error - a short period simply
        !> has one usable side - so the lag is skipped and the caller counts
        !> how many survived.
        if (abs(lag) >= nrow) cycle
        call CovarianceW(ColW, ColGas, nrow, lag, cov)
        if (cov == error) cycle
        n = n + 1
        s1 = s1 + cov
        s2 = s2 + cov * cov
    end do

    if (n < min_lags) return

    mean = s1 / dble(n)
    !> Sample standard deviation, n-1, matching StDevNoError elsewhere in the
    !> tree. Clamped at zero because the algebraic form can go very slightly
    !> negative when every covariance in the window is the same number.
    sd = dsqrt(max(0d0, (s2 - dble(n) * mean * mean) / dble(n - 1)))
end subroutine CrossCovarianceScatter
