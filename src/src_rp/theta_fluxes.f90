!***************************************************************************
! theta_fluxes.f90
! ----------------
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
! \brief       The potential temperature flux w'theta' and the virtual
!              potential temperature (buoyancy) flux w'theta_v', kinematic,
!              as innFLUX reports them (Striednig et al. 2020, Atmos. Meas.
!              Tech. 13, 1447-1465; MET.wtheta and MET.wtheta_v of
!              innFLUX_step1.m). Opt-in: RPSetup%out_theta_flux.
!
!              Potential temperature is theta = T (p0 / p)^kappa with
!              p0 = 1000 hPa and kappa = R/cp = 0.286, innFLUX's constants
!              (Stull 1988, Eq. 1.5.1c). Over one averaging period p is a
!              constant, so w'theta' = w'T' (p0/p)^kappa.
!
!              w'T' is the humidity-corrected kinematic heat flux, made the
!              way Fluxes23_rp corrects H (van Dijk et al. 2004, Eq. 3.53, the
!              revision of Schotanus et al. 1983) but from the level 0
!              covariances: like innFLUX's, these fluxes carry no spectral
!              correction. Without a water flux the sonic temperature cannot
!              be separated from humidity and both fluxes are the error code;
!              from a thermocouple, w'T' is the covariance itself.
!
!              theta_v = theta (1 + 0.61 q), q the specific humidity, so to
!              first order w'theta_v' = (p0/p)^kappa [w'T' (1 + 0.61 q)
!              + 0.61 T w'q'], with w'q' = E / rho the kinematic water vapour
!              flux. innFLUX writes the same factor as 0.38 on the water mole
!              fraction (0.61 x 0.622).
! \author      Jonathan Muller
! \sa          fluxes23_rp.f90, write_out_fluxnet.f90
!***************************************************************************
subroutine ThetaFluxes()
    use m_rp_global_var
    implicit none
    !> local variables
    real(kind = dbl), parameter :: p0_hpa = 1000d0
    real(kind = dbl), parameter :: kappa = 0.286d0
    real(kind = dbl) :: exner
    real(kind = dbl) :: wT
    real(kind = dbl) :: wq

    Essentials%wtheta = error
    Essentials%wthetav = error
    if (.not. RPSetup%out_theta_flux) return
    !> Stats%Pr is in Pa.
    if (.not. Stats%Pr > 0d0) return
    if (Stats%Cov(w, ts) == error) return
    exner = (p0_hpa / (Stats%Pr * 1d-2)) ** kappa

    wT = error
    wq = error
    if (Flux0%E /= error .and. RHO%a > 0d0) wq = Flux0%E / RHO%a
    if (E2Col(ts)%instr%category == 'sonic') then
        if (wq /= error .and. Ambient%Q >= 0d0 .and. Ambient%alpha /= error) &
            wT = Stats%Cov(w, ts) &
                - Ambient%alpha * Stats%Mean(ts) * wq &
                - Ambient%alpha * Ambient%Q * Stats%Cov(w, ts)
    else
        wT = Stats%Cov(w, ts)
    end if
    if (wT == error) return

    Essentials%wtheta = wT * exner
    if (wq /= error .and. Ambient%Q >= 0d0 .and. Ambient%Ta > 0d0) &
        Essentials%wthetav = (wT * (1d0 + 0.61d0 * Ambient%Q) &
            + 0.61d0 * Ambient%Ta * wq) * exner
end subroutine ThetaFluxes
