!***************************************************************************
! fluxes0_rp.f90
! --------------
! Copyright © 2011-2026, LI-COR Biosciences, Gerardo Fratini
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
! \brief       Calculates uncorrected fluxes (refer to pseudo-code)
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine Fluxes0_rp(printout)
    use m_rp_global_var
    implicit none
    !> in/out variables
    logical, intent(in) :: printout
    !> local variables
    real(kind = dbl) :: Tp
    real(kind = dbl) :: dens_gain
    integer :: msl
    integer :: gas
    integer :: wsl
    integer :: wsx
    include '../src_common/interfaces_1.inc'

    !> The site's latent heat flux and evapotranspiration come from the
    !> first water record, not from the h2o slot.
    wsl = PrimaryWaterSlot()
    !> An always-in-bounds stand-in for wsl inside guard expressions.
    !>
    !> Fortran does not mandate short-circuit `.and.`, so
    !> `wsl >= firstGas .and. X(wsl)` still evaluates X(wsl) - and with no
    !> hygrometer wsl is 0, which is out of bounds. The wsl >= firstGas test
    !> still decides the outcome; wsx only keeps the subscript legal while it
    !> is being decided.
    wsx = max(wsl, firstGas)

    if (printout) write(*,'(a)', advance = 'no') &
        '  Calculating fluxes Level 0..'

    Flux0 = errFlux

    !> Sensible heat flux, H in [W m-2], Cp in [J Kg-1K-1]
    if (Ambient%RhoCp > 0d0 .and. Stats%Cov(w, ts) /= error) then
        Flux0%H = Ambient%RhoCp * Stats%Cov(w, ts)
    else
        Flux0%H = error
    end if

    !> Random error on sensible heat flux
    if (RUsetup%meth /= 'none') then
        if (Ambient%RhoCp > 0d0 .and. Essentials%rand_uncer(ts) /= error) &
            Essentials%rand_uncer(ts) = &
                Essentials%rand_uncer(ts) * Ambient%RhoCp
    end if

    !> Internal sensible heat flux, Hint in [W m-2], Cp in [J Kg-1K-1]
    !>
    !> One pass over the gas slots. The four arms this replaces read four
    !> named scalars, and the water arm already had to route its result to the
    !> resolved slot while reading the field named for the historical one -
    !> the covariances are now indexed by slot, so both ends agree.
    Flux0%Hi_gas(firstGas:lastGas) = error
    if (Ambient%RhoCp > 0d0) then
        do msl = firstGas, lastGas
            if (Stats%tc_cov_tl(msl) /= error) &
                Flux0%Hi_gas(msl) = Ambient%RhoCp * Stats%tc_cov_tl(msl)
        end do
    end if

    !> Uncorrected flux of each gas.
    !>
    !> One loop over the configured gases, replacing four near-identical
    !> blocks. They differed only in the molar-density scale factor: H2O is
    !> reported in mmol m-2 s-1 and the other gases in umol m-2 s-1. That is a
    !> property of the species, not of the slot, so it is decided from the
    !> species - a slot past the fourth holds whichever gas the project put
    !> there.
    do msl = firstGas, lastGas
        if (.not. E2Col(msl)%present) then
            Flux0%gas(msl) = error
            cycle
        end if
        if (GasSlotIsWater(msl)) then
            dens_gain = 1d0
        else
            dens_gain = 1d3
        end if

        select case (E2Col(msl)%measure_type)
            case ('molar_density')
                if (Stats%Cov(w, msl) /= error) then
                    Flux0%gas(msl) = Stats%Cov(w, msl) * dens_gain
                else
                    Flux0%gas(msl) = error
                end if
                if (RUsetup%meth /= 'none') then
                    if (Essentials%rand_uncer(msl) /= error) &
                        Essentials%rand_uncer(msl) = &
                            Essentials%rand_uncer(msl) * dens_gain
                end if

            case ('mixing_ratio')
                if (Ambient%Vd > 0d0 .and. Stats%Cov(w, msl) /= error) then
                    Flux0%gas(msl) = Stats%Cov(w, msl) / Ambient%Vd
                else
                    Flux0%gas(msl) = error
                end if
                if (RUsetup%meth /= 'none') then
                    if (Essentials%rand_uncer(msl) /= error &
                        .and. Ambient%Vd > 0d0) then
                        Essentials%rand_uncer(msl) = &
                            Essentials%rand_uncer(msl) / Ambient%Vd
                    else
                        Essentials%rand_uncer(msl) = error
                    end if
                end if

            case ('mole_fraction')
                if (Ambient%Va > 0d0 .and. Stats%Cov(w, msl) /= error) then
                    Flux0%gas(msl) = Stats%Cov(w, msl) / Ambient%Va
                else
                    Flux0%gas(msl) = error
                end if
                if (RUsetup%meth /= 'none') then
                    if (Essentials%rand_uncer(msl) /= error &
                        .and. Ambient%Va > 0d0) then
                        Essentials%rand_uncer(msl) = &
                            Essentials%rand_uncer(msl) / Ambient%Va
                    else
                        Essentials%rand_uncer(msl) = error
                    end if
                end if
        end select
    end do

    !> Latent heat flux and evapotranspiration. These are one-per-site
    !> quantities and come from the primary water record; a project that
    !> describes no water performs none of them, rather than computing them
    !> from whatever gas happens to sit in the h2o slot.
    if (wsl >= firstGas .and. Flux0%gas(wsx) /= error &
        .and. Ambient%lambda > 0d0) then
        Flux0%LE = Flux0%gas(wsl) * Ambient%lambda * MW_H2O * 1d-3
        Flux0%E  = Flux0%gas(wsl) * MW_H2O * 1d-3
        Flux0%ET = Flux0%gas(wsl) * h2o_to_ET
    else
        Flux0%LE = error
        Flux0%E  = error
        Flux0%ET = error
    end if

    !> Random uncertainty on latent heat flux, lambda in [J+1kg-1]
    if (RUsetup%meth /= 'none') then
        if (wsl >= firstGas .and. Essentials%rand_uncer(wsx) /= error &
            .and. Ambient%lambda > 0d0) then
            Essentials%rand_uncer_LE = &
                Essentials%rand_uncer(wsl) * Ambient%lambda * MW_H2O * 1d-3
            Essentials%rand_uncer_ET = &
                Essentials%rand_uncer(wsl) * h2o_to_ET
        else
            Essentials%rand_uncer_LE = error
            Essentials%rand_uncer_ET = error
        end if
    end if

    !> Level 0 evapotranspiration flux [kg m-2 -1]
    !> with H2O covariances at timelags of other scalars
    !>
    !> Nine arms before - three measure types by three named gas slots - which
    !> differed only in the divisor. The divisor is a property of the water
    !> record's measure type, so it is chosen once and the gases loop inside.
    !> The water slot itself has no entry in h2ocov_tl, so the loop skips it
    !> the same way the unrolled arms did by omission.
    !> The divisor belongs to the hygrometer whose covariance this is, and
    !> that is the one the gas names - not the site's. Chosen inside the loop
    !> for that reason: with two analysers the two hygrometers can differ in
    !> measure type, and the gate on the primary's path type meant a gas on a
    !> closed-path analyser lost this term entirely whenever the *primary*
    !> hygrometer happened to be open-path.
    !>
    !> TimeLagHandle has already declined to compute a covariance where the
    !> reference is unusable, so an `error` here is simply that answer
    !> arriving.
    Flux0%E_gas(firstGas:lastGas) = error
    do gas = firstGas, lastGas
        if (Stats%h2ocov_tl(gas) == error) cycle
        msl = E2Col(gas)%moist_ref
        if (msl < firstGas .or. msl > lastGas) cycle

        dens_gain = error
        select case (E2Col(msl)%measure_type)
            case ('molar_density')
                dens_gain = 1d0
            case ('mole_fraction')
                if (Ambient%Va > 0d0) dens_gain = 1d0 / Ambient%Va
            case ('mixing_ratio')
                if (Ambient%Vd > 0d0) dens_gain = 1d0 / Ambient%Vd
        end select
        if (dens_gain == error) cycle

        Flux0%E_gas(gas) = Stats%h2ocov_tl(gas) * MW_H2O * 1d-3 * dens_gain
    end do

    !> Friction velocity [m s-1]
    if (Stats%Cov(u, w) /= error .and. Stats%Cov(v, w) /= error) then
        Ambient%us = (Stats%Cov(u, w)**2 + Stats%Cov(v, w)**2)**(0.25d0)
    else
        Ambient%us = error
    end if
    Flux0%ustar = Ambient%us
    Essentials%ustar = Ambient%us

    !> Momentum flux [kg m-1 s-2], after Van Dijk et al. 2004 Eq. 2.44
    if (RHO%a > 0d0 .and. Ambient%us >= 0d0) then
        Flux0%tau = sign(RHO%a * Ambient%us ** 2d0, Stats%Cov(u, w))
    else
        Flux0%tau = error
    end if

    !> Random error on momentum flux
    if (RUsetup%meth /= 'none') then
        if (Essentials%rand_uncer(u) /= error &
            .and. RHO%a > 0d0 .and. Ambient%us >= 0d0) then
            Essentials%rand_uncer(u) = Essentials%rand_uncer(u) * RHO%a
        else
            Essentials%rand_uncer(u) = error
        end if
    end if

    !> Potential temperature
    if (Stats%Pr > 0d0 .and. Stats%T > 0d0) then
        Tp = Stats%T * (1d5 / Stats%Pr)**(.286d0)
    else
        Tp = error
    end if

    !> Monin-Obukhov length (L = - (Tp^ /(k*g))*(ustar**3/(w'Tp')^ in m)
    if (Stats%Cov(w, ts) /= 0d0 .and. Stats%Cov(w, ts) /= error .and. &
        Ambient%us > 0d0 .and. Tp > 0d0) then
        Ambient%L = -Tp * (Ambient%us**3) / (vk * g * Stats%Cov(w, ts))
    else
        Ambient%L = error
    end if
    Essentials%L = Ambient%L

    !> Monin-Obukhov stability parameter (zL = z/L)
    if (Ambient%L /= 0d0 .and. Ambient%L /= error) then
        Ambient%zL = (E2Col(u)%Instr%height - Metadata%d) / Ambient%L
    else
        Ambient%zL = error
    end if
    Essentials%zL = Ambient%zL

    !> SL dynamic temperature(T*), see e.g. Foken and Wichura (1996)
    if (Ambient%us > 0d0 .and. Stats%Cov(w, ts) /= error) then
        Ambient%Ts = - Stats%Cov(w, ts) / Ambient%us
    else
        Ambient%Ts = error
    end if

    !> Bowen ration (Bowen, 1926, Phyis Rev)
    if (Flux0%LE /= 0d0 .and. Flux0%LE /= error .and. Flux0%H /= error) then
        Ambient%Bowen = Flux0%H / Flux0%LE
    else
        Ambient%Bowen = error
    end if
    if (printout) write(*,'(a)')   ' Done.'
end subroutine Fluxes0_rp
