!***************************************************************************
! bpcf_eugster_senn_95.f90
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
! \brief       Calculates spectral correction factors for gas fluxes
!              according to Eugster and Senn (1995, BLM 74, 321-340)
! \author      Jonathan Muller
! \note        The damping of the gas measurement is modelled as a first-order
!              low-pass filter, 1 / (1 + 4 pi^2 f^2 L^2), applied to a
!              Kaimal-type cospectrum. L ("inductance" in the paper) is a
!              system time constant [s] lumping tube, chemistry and analyser.
!              The measurable fraction xi of the flux has a closed form, so
!              no transfer function is integrated; the factor is 1/xi.
!              The paper builds it on NO2, but it holds for any scalar whose
!              sensor shows frequency-dependent damping, so it is applied to
!              every gas slot. Anemometric fluxes are left as
!              BPCF_AnemometricFluxes set them.
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine BPCF_EugsterSenn95(measuring_height, displ_height, loc_var_present, &
        LocInstr, wind_speed, zL, ac_frequency, avrg_length, &
        detrending_time_constant, detrending_method, printout)
    use m_common_global_var
    implicit none
    !> in/out variables
    real(kind = dbl), intent(in) :: measuring_height
    real(kind = dbl), intent(in) :: displ_height
    logical, intent(in) :: loc_var_present(GHGNumVar)
    type(InstrumentType), intent(in) :: LocInstr(GHGNumVar)
    real(kind = dbl), intent(in) :: wind_speed
    real(kind = dbl), intent(in) :: zL
    real(kind = dbl), intent(in) :: ac_frequency
    integer, intent(in) :: avrg_length
    integer, intent(in) :: detrending_time_constant
    character(2), intent(in) :: detrending_method
    logical, intent(in) :: printout
    !> local variables
    integer :: gas
    real(kind = dbl) :: z_over_u
    real(kind = dbl) :: tau_L
    real(kind = dbl) :: xi
    !> Once per gas and run, not once per period: the setting does not change
    !> between periods, and a line per half hour would bury the log.
    logical, save :: warned_no_tau(GHGNumVar) = .false.
    character(32), external :: GasOutputLabel
    real(kind = dbl), external :: EugsterSenn95Fraction


    if (printout) write(*,'(a)') '   Band-pass correction for gas fluxes. Method: Eugster and Senn (1995)..'
    if (printout) write(ulog,'(a)') '   Band-pass correction for gas fluxes. Method: Eugster and Senn (1995)..'

    !> The paper corrects the high-frequency damping only. Low-frequency
    !> losses are added the way the 'none' method adds them, so switching from
    !> Moncrieff changes the high-frequency part and nothing else.
    if (EddyFlowProj%lf_meth == 'analytic') &
        call BPCF_OnlyLowFrequencyCorrection(measuring_height, &
            displ_height, loc_var_present, wind_speed, zL, ac_frequency, &
            avrg_length, detrending_time_constant, detrending_method)

    !> z/u and stability are the model's only atmospheric inputs. An error
    !> code would pass for a very unstable or very fast period, so refuse.
    if (wind_speed == error .or. wind_speed <= 0d0 .or. zL == error &
        .or. measuring_height == error .or. displ_height == error &
        .or. measuring_height - displ_height <= 0d0) then
        BPCF%of(firstGas:lastGas) = error
        if (printout) write(*,'(a)') '   Done.'
        if (printout) write(ulog,'(a)') '   Done.'
        return
    end if
    z_over_u = (measuring_height - displ_height) / wind_speed

    do gas = firstGas, lastGas
        if (.not. loc_var_present(gas)) cycle
        if (BPCF%of(gas) == error) cycle

        !> The gas's own L when the project states one, else the response
        !> time of its analyser. The latter usually leaves out tube and
        !> chemistry damping, which is what the paper's L is for.
        tau_L = EddyFlowProj%gas(gas - firstGas + 1)%es_tau
        if (tau_L < 0d0) tau_L = LocInstr(gas)%tau
        if (tau_L < 0d0) then
            if (.not. warned_no_tau(gas)) then
                write(*,'(a)') '   Eugster and Senn (1995): no time constant for ' &
                    // trim(GasOutputLabel(gas)) // ', its flux is not corrected for high-frequency damping.'
                write(ulog,'(a)') '   Eugster and Senn (1995): no time constant for ' &
                    // trim(GasOutputLabel(gas)) // ', its flux is not corrected for high-frequency damping.'
                warned_no_tau(gas) = .true.
            end if
            tau_L = 0d0
        end if

        xi = EugsterSenn95Fraction(z_over_u, zL, tau_L)
        BPCF%of(gas) = BPCF%of(gas) / xi
        !> Eugster and Senn discard every period in which less than 40% of
        !> the flux is measurable. The factor stays in BPCF so the output
        !> shows why the flux is missing.
        if (xi < 0.4d0) BPCF_discard(gas) = .true.
    end do

    if (printout) write(*,'(a)') '   Done.'
    if (printout) write(ulog,'(a)') '   Done.'
end subroutine BPCF_EugsterSenn95

!***************************************************************************
!
! \brief       Measurable fraction xi of a flux damped with time constant
!              tau_L, Eugster and Senn (1995, Eqs. 32 and 33)
! \author      Jonathan Muller
! \note        z_over_u is (z - d) / u [s]; zL is the stability parameter.
!              Both equations were fitted for -2 <= zL <= 2. Below -2 the
!              unstable fit does not depend on zL anyway; above 2, zL is
!              held at 2.
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
function EugsterSenn95Fraction(z_over_u, zL, tau_L) result(xi)
    use m_common_global_var
    implicit none
    !> in/out variables
    real(kind = dbl), intent(in) :: z_over_u
    real(kind = dbl), intent(in) :: zL
    real(kind = dbl), intent(in) :: tau_L
    real(kind = dbl) :: xi
    !> local variables
    real(kind = dbl) :: zeta
    real(kind = dbl) :: n0


    !> An undamped sensor measures the whole flux.
    if (tau_L <= 0d0) then
        xi = 1d0
        return
    end if

    if (zL <= 0d0) then
        !> Unstable and neutral, Eq. 32: a least-squares fit to the numerical
        !> integral of Eq. 29. The fit tends to 0.725 pi/2 = 1.14 as L goes to
        !> zero, past the 1 the integral cannot exceed, so it is capped there.
        xi = 0.725d0 * datan(1.24d0 * dlog(z_over_u / tau_L + 1d0) + 0.21d0)
        xi = min(xi, 1d0)
    else
        !> Stable, Eqs. 28 and 33. As printed, Eq. 33 is
        !>   3/sqrt(6) (z/u) (sqrt(6) z/u - 4 pi L n0) / (3 (z/u)^2 - 8 pi^2 L^2 n0^2)
        !> and numerator and denominator share the root z/u = 4 pi L n0 / sqrt(6),
        !> where it reads 0/0. Cancelling that common factor leaves the form
        !> below, which is the same function everywhere else and defined there.
        zeta = min(zL, 2d0)
        n0 = 0.23d0 * (1d0 + 6.4d0 * zeta)**0.75d0
        xi = z_over_u / (z_over_u + 4d0 * p / dsqrt(6d0) * tau_L * n0)
    end if
end function EugsterSenn95Fraction
