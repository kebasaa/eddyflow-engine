!***************************************************************************
! kid.f90
! -------
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
! \brief       KID test: kurtosis on stochastically detrended data
! \author      Jonathan Muller, ETH Zurich
!
!***************************************************************************
! Reference: Vitale, D. et al. (2020). Biogeosciences Discussion.
!            KID = kurtosis of stochastic-residual; ZCD = zero-crossing density.
!
! \note        RFlux's inst_prob_test() computes two kurtosis-of-increment
!              variants: KID0 = 3+kurtosis(diff(x)), and KID =
!              3+kurtosis(diff(na.omit(x_r))), where x_r is x with the
!              sample right after every near-zero step replaced by NA -
!              i.e. flat-lined transitions are excluded before
!              differencing again, so a repeated-value run does not
!              contribute a string of exact-zero (or spuriously large,
!              depending on where the flatline sits) differences to the
!              kurtosis. Essentials%KID(:) is this second, flatline-
!              corrected variant: it is what cleanFlux.R's own 30/50
!              severe/moderate thresholds are calibrated against (index 4
!              of its threshold vectors; KID0 at index 3 carries no
!              threshold at all), and the FLUXNET output column is
!              labelled plain "KID", not "KID0". ZCD (zero-crossing
!              density) still counts crossings on the plain, uncorrected
!              residual - RFlux has no flatline-corrected counterpart for
!              it, so there is nothing to port there.
subroutine KID(Set, nrow, ncol)
    use m_rp_global_var
    implicit none
    integer, intent(in) :: nrow, ncol
    real(kind = dbl), intent(in) :: Set(nrow, ncol)
    integer :: icol, nout
    integer, external :: CountZeroCrossings
    real(kind = dbl) :: residuals(nrow, ncol)
    real(kind = dbl) :: kid_diffs(nrow)

    do icol = u, ts
        call VariableStochasticDetrending(Set(:, icol), residuals(:, icol), nrow)
        Essentials%ZCD(icol) = CountZeroCrossings(residuals(:, icol), nrow)
        call VariableFlatlineCorrectedDiff(Set(:, icol), kid_diffs, nrow, nout)
        call KurtosisNoError(kid_diffs, nrow, 1, Essentials%KID(icol), error)
    end do
    !> Every configured gas. The FLUXNET writer emits a KID and a ZCD per
    !> configured gas, so a producer that stopped at the fourth left the rest
    !> reporting whatever Essentials held - a kurtosis index of exactly zero,
    !> which reads as a measurement. %present already carries the
    !> was-configured guard, and its else arm is the sentinel.
    do icol = firstGas, lastGas
        if (E2Col(icol)%present) then
            call VariableStochasticDetrending(Set(:, icol), residuals(:, icol), nrow)
            Essentials%ZCD(icol) = CountZeroCrossings(residuals(:, icol), nrow)
            call VariableFlatlineCorrectedDiff(Set(:, icol), kid_diffs, nrow, nout)
            call KurtosisNoError(kid_diffs, nrow, 1, Essentials%KID(icol), error)
        else
            Essentials%KID(icol) = error
            Essentials%ZCD(icol) = ierror
        end if
    end do
end subroutine KID


!***************************************************************************
!
! \brief       Flatline-corrected first difference, for KID only.
! \author      Jonathan Muller, ETH Zurich
! \note        RFlux (inst_prob_test.R): d0 <- diff(x); ind <- which(abs(d0)
!              < 1e-3); x_r <- replace(x, ind+1, NA); diff(na.omit(x_r)).
!              This compacts x FIRST (dropping every marked-NA and every
!              already-error sample) and only then differences the
!              shortened series - so a removed run is bridged by one
!              difference spanning the gap, not left as a separate NA at
!              each dropped position the way VariableStochasticDetrending's
!              plain diff (used elsewhere for time-lag pre-whitening, and
!              deliberately left untouched by this) would leave it.
!***************************************************************************
subroutine VariableFlatlineCorrectedDiff(Var, Primes, N, Nout)
    use m_common_global_var
    implicit none
    integer, intent(in) :: N
    real(kind = dbl), intent(in) :: Var(N)
    real(kind = dbl), intent(out) :: Primes(N)
    integer, intent(out) :: Nout
    logical :: keep(N)
    real(kind = dbl) :: compact(N)
    integer :: i, m

    keep = .true.
    do i = 1, N
        if (Var(i) == error) keep(i) = .false.
    end do
    do i = 2, N
        if (Var(i) == error .or. Var(i-1) == error) cycle
        if (dabs(Var(i) - Var(i-1)) < 1d-3) keep(i) = .false.
    end do

    m = 0
    do i = 1, N
        if (.not. keep(i)) cycle
        m = m + 1
        compact(m) = Var(i)
    end do

    Primes = error
    Nout = max(0, m - 1)
    do i = 1, Nout
        Primes(i) = compact(i + 1) - compact(i)
    end do
end subroutine VariableFlatlineCorrectedDiff


integer function CountZeroCrossings(arr, nrow)
    use m_rp_global_var
    implicit none
    integer, intent(in) :: nrow
    real(kind = dbl), intent(in) :: arr(nrow)
    integer :: irow, current_sign, previous_sign

    CountZeroCrossings = 0
    previous_sign = 0
    do irow = 1, nrow
        if (arr(irow) == error .or. arr(irow) == 0d0) cycle

        if (arr(irow) > 0d0) then
            current_sign = 1
        else
            current_sign = -1
        end if

        if (previous_sign /= 0 .and. current_sign /= previous_sign) &
            CountZeroCrossings = CountZeroCrossings + 1
        previous_sign = current_sign
    end do
end function CountZeroCrossings
