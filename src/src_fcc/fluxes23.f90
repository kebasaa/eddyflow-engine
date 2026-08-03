!***************************************************************************
! fluxes23.f90
! ------------
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
! \brief       Completes flux correction (Level 2/3). Applies WPL and \n
!              spectral corrections in the appropriate order
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine Fluxes23(lEx)
    use m_fx_global_var
    implicit none
    !> In/out variables
    type(ExType), intent(inout) :: lEx
    !> local variables
    real(kind = dbl) :: Tp
    real(kind = dbl) :: E_nowpl
    !> Water vapour terms of the H2O that corrects each gas. Equal to the
    !> global lEx%sigma / lEx%RHO%w unless the project gives that gas its own
    !> moisture reference, so a single-analyser site is unaffected.
    real(kind = dbl) :: sigma_g, rhow_g
    integer :: msl
    integer :: wsl
    include '../src_common/interfaces_1.inc'
    real(kind = dbl), parameter :: alpha = 0.51d0


    !> Water's own slot, from the records. FCC is what publishes the FLUXNET
    !> row under fcc_follows, so leaving this at the h2o constant undid the
    !> RP-side work entirely: on a project with water anywhere but record
    !> two, FCC wrote water's flux into slot 6 - overwriting whatever gas is
    !> actually there - and left the real water on the generic trace-gas path.
    wsl = PrimaryWaterSlot()

    Flux2 = errFlux
    Flux3 = errFlux

    !> Level 2 end 3 internal sensible heat, do nothing
    !>
    !> A pass-through, so it carries the whole gas block. Spelled out for
    !> co2, the water slot, ch4 and gas4, it dropped every gas past the
    !> fourth record: Fluxes0_rp computes Hi_gas for each configured gas and
    !> the FLUXNET file writes an H_CELL_* column for each, but only four
    !> survived the level 1 -> 2 -> 3 chain, so the rest were error codes
    !> whatever the analyser reported. It also relocated them - which four
    !> slots those names pick out depends on where water sits, so two
    !> projects differing only in record order disagreed about H_CELL.
    Flux2%Hi_gas(firstGas:lastGas) = Flux1%Hi_gas(firstGas:lastGas)
    Flux3%Hi_gas(firstGas:lastGas) = Flux2%Hi_gas(firstGas:lastGas)

    !> Level 2 evapotranspiration WPL corrected, including Burba if the case
    if (EddyFlowProj%wpl) then
        if (wsl >= firstGas .and. lEx%gas_instr(wsl)%path_type == 'open') then
            if (lEx%RhoCp > 0d0 .and. lEx%Ta > 0d0 &
                .and. Flux1%E /= error .and. Flux1%H /= error &
                .and. lEx%sigma /= error) then
                    !> Open-path uses Webb et al. (1980)
                    !> Note that Burba terms are forced to zero
                    !> if analyzer is /= LI-7500
                    Flux2%E = (1d0 + mu * lEx%sigma) * Flux1%E &
                            + (1d0 + mu * lEx%sigma) &
                            * (Flux1%H + lEx%Burba%h_top + lEx%Burba%h_bot + lEx%Burba%h_spar) &
                            * lEx%RHO%w / (lEx%RhoCp * lEx%Ta)
            else
                Flux2%E = error
            end if
        else
            !> Closed-path uses Ibrom et al. (2007) if conversion to mixing
            !> ratio did not already occur (which implies that some variables
            !> were missing)
            select case(lEx%measure_type(wsl))
                case ('molar_density', 'mole_fraction')
                    if (Flux1%E /= error .and. lEx%sigma /= error &
                        .and. lEx%Vcell(wsl) > 0d0 .and. lEx%Va > 0d0) then

                        if (Flux1%Hi_gas(wsl) /= error &
                            .and. lEx%cov_w_pcell(wsl) /= error) then
                            !> Complete formulation, should actually never be
                            !> used cause conversion to mixing ratio should have
                            !> already happened if everything is available
                            Flux2%E = (1d0 + mu * lEx%sigma) * Flux1%E &
                                * lEx%Vcell(wsl) / lEx%Va &
                                + (1d0 + mu * lEx%sigma) * Flux1%Hi_gas(wsl) &
                                * lEx%RHO%w / (lEx%RhoCp * lEx%Tcell_at(wsl)) &
                                - (1d0 + mu * lEx%sigma) * lEx%cov_w_pcell(wsl) &
                                * lEx%RHO%w / (lEx%Pcell_at(wsl))

                        elseif (Flux1%Hi_gas(wsl) /= error) then
                            !> Correct only for effect of T
                            Flux2%E = (1d0 + mu * lEx%sigma) * Flux1%E &
                                * lEx%Vcell(wsl) / lEx%Va &
                                + (1d0 + mu * lEx%sigma) * Flux1%Hi_gas(wsl) &
                                * lEx%RHO%w / (lEx%RhoCp * lEx%Tcell_at(wsl))

                        elseif (lEx%cov_w_pcell(wsl)  /= error) then
                            !> Correct only for effect of P
                            Flux2%E = (1d0 + mu * lEx%sigma) * Flux1%E &
                                * lEx%Vcell(wsl) / lEx%Va &
                                - (1d0 + mu * lEx%sigma) * lEx%cov_w_pcell(wsl) &
                                * lEx%RHO%w / (lEx%Pcell_at(wsl))
                        else
                            !> Can't correct for T and P
                            Flux2%E = Flux1%E * lEx%Vcell(wsl) / lEx%Va
                        end if
                    else
                        Flux2%E = error
                    end if
                case ('mixing_ratio')
                    Flux2%E = Flux1%E
            end select
        end if
    else
        !> If WPL should not be applied
        Flux2%E = Flux1%E
    end if
    !> If Flux1 was error, then set also 2 to error
    if (Flux1%E == error) Flux2%E = error

    !> Level 2 h2o and latent heat flux
    if (Flux2%E /= error) then
        Flux2%gas(wsl) = Flux2%E * 1d3 / MW_H2O
        Flux2%ET = Flux2%gas(wsl) * h2o_to_ET
        if (lEx%lambda /= error) then
            Flux2%LE = Flux2%E * lEx%lambda
        else
            Flux2%LE = error
        end if
    else
        Flux2%gas(wsl) = error
        Flux2%LE  = error
        Flux2%ET  = error
    end if

    !> Level 2 evapotranspiration fluxes with H2O covariances
    !> at time-lags of other scalars. Do nothing, WPL is deleterious here.
    !> Another pass-through, so it carries the whole gas block; the water
    !> slot has no entry by construction and copying it costs nothing.
    Flux2%E_gas(firstGas:lastGas) = Flux1%E_gas(firstGas:lastGas)

    !> Level 2 Sensible heat
    if (lEx%instr(sonic)%category == 'sonic') then
        !> Corrected for humidity, after Van Dyjk et al. (2004) eq. 3.53
        !> revising Schotanus et al. (1983)
        if (Flux1%H /= error) then
            if(lEx%Flux0%E /= error .and. lEx%cov_w(ts) /= error &
                .and. lEx%RHO%a > 0d0 .and. lEx%Q >= 0d0 &
                .and. lEx%RhoCp > 0d0 .and. alpha /= error) then
                Flux2%H = Flux1%H &
                    - lEx%RhoCp * alpha * lEx%Ts * lEx%Flux0%E / lEx%RHO%a &
                    - lEx%RhoCp * alpha * lEx%Q * lEx%cov_w(ts)
                    !> alternative
                    !- lEx%RhoCp * alpha * lEx%Ta * lEx%Flux0%E / lEx%RHO%a
            else
                Flux2%H = Flux1%H
            end if
        else
            Flux2%H = error
        end if
    else
        !> Equal Level 1 if Ts is not coming from a sonic, rather from
        !> a fast temperature sensor (such as a thermocouple)
        Flux2%H = Flux1%H
    end if

    !> Map for temperature factor (Van Dijk et al. 2004, eq.3.1)
    !> Currently not applied (to investigate better)
    !if(lEx%Tmap /= error .and. Flux2%H /= error) then
        !Flux2%H  = Flux2%H * lEx%Tmap
    !end if

    !> Level 3 sensible heat, spectral corrected
    if(Flux2%H /= error .and. BPCF%of(w_ts) /= error) then
        Flux3%H = Flux2%H * BPCF%of(w_ts)
    else
        Flux3%H = error
    end if

    !> Level 3 for evapotranspiration: for open path, WPL again with corrected H
    !> Starts again from Level 1 of E, Level 2 was only used to calculate H Level 3.
    if(EddyFlowProj%wpl .and. wsl >= firstGas &
        .and. lEx%gas_instr(wsl)%path_type == 'open') then
        if (lEx%RhoCp > 0d0 .and. lEx%Ta > 0d0 .and. Flux1%E /= error &
            .and. Flux1%H /= error .and. lEx%sigma /= error) then
            Flux3%E = (1d0 + mu * lEx%sigma) * Flux1%E &
                + (1d0 + mu * lEx%sigma) &
                * (Flux3%H + lEx%Burba%h_top + lEx%Burba%h_bot + lEx%Burba%h_spar)&
                * lEx%RHO%w / (lEx%RhoCp * lEx%Ta)
        else
            Flux3%E = Flux2%E
        end if
    else
        Flux3%E = Flux2%E
    end if

    !> Level 3 latent heat fluxes with H2O covariances at
    !> timelags of other scalars
    !> Do nothing
    Flux3%E_gas(firstGas:lastGas) = Flux2%E_gas(firstGas:lastGas)

    !> Level 3 h2o flux and latent heat flux
    if (Flux3%E /= error) then
        Flux3%gas(wsl) = Flux3%E * 1d3 / MW_H2O
        Flux3%ET = Flux3%gas(wsl) * h2o_to_ET
        if (lEx%lambda /= error) then
            Flux3%LE = Flux3%E * lEx%lambda
        else
            Flux3%LE = error
        end if
    else
        Flux3%gas(wsl) = error
        Flux3%LE  = error
        Flux3%ET  = error
    end if

    !> Calculate E_nowpl for closed and open path systems
    if (wsl >= firstGas .and. lEx%gas_instr(wsl)%path_type == 'closed') then
        if (Flux1%E /= error .and. BPCF%of(wsl) /= error) then
            E_nowpl = Flux1%E * BPCF%of(wsl)
        elseif(Flux1%E /= error) then
            E_nowpl = Flux1%E
        else
            E_nowpl = error
        end if
    else
        E_nowpl = Flux1%E
    end if

    !> Apply spectral correction to h2o and E/LE Level 3 fluxes for closed path
    if (wsl >= firstGas .and. Flux3%E /= error) then
        if (lEx%gas_instr(wsl)%path_type == 'closed') then
            !> Level 3, spectral correction
            Flux3%gas(wsl) = Flux3%gas(wsl) * BPCF%of(wsl)
            Flux3%E   = Flux3%E   * BPCF%of(wsl)
            Flux3%LE  = Flux3%LE  * BPCF%of(wsl)
            Flux3%ET  = Flux3%ET  * BPCF%of(wsl)
        end if
    end if

    if (wsl < firstGas) then
        Flux3%E = error; Flux3%LE = error; Flux3%ET = error
    elseif (.not. lEx%var_present(wsl)) then
        Flux3%gas(wsl) = error
        Flux3%E   = error
        Flux3%LE  = error
        Flux3%ET  = error
    end if


    !> Level 2 other gases.
    !>
    !> One loop over the gases, replacing three near-duplicate blocks that had
    !> drifted apart; see Level2GasFlux for the two differences resolved.
    !>
    !> Runs the full gas range: lEx%gas_instr is indexed by gas slot and the
    !> ex file now carries an analyser for every configured gas.
    do msl = firstGas, lastGas
        !> H2O's own flux is the evapotranspiration handled above.
        if (msl == wsl) cycle
        call MoistTerms(msl, sigma_g, rhow_g)
        call Level2GasFlux(msl, sigma_g, rhow_g)
    end do

    !> Level 3 other gases. For closed path apply the spectral correction now
    !> (e.g. Ibrom et al. 2007); for open path it is already included.
    do msl = firstGas, lastGas
        if (msl == wsl) cycle
        if (lEx%gas_instr(msl)%path_type == 'closed' &
            .and. Flux2%gas(msl) /= error) then
            !> No correction factor means the corrected flux is
            !> unavailable, not that it equals the uncorrected one - and
            !> certainly not Flux2 times the error sentinel, which is what the
            !> unguarded multiply produced for any gas the spectral chain
            !> could not reach.
            if (BPCF%of(msl) /= error) then
                Flux3%gas(msl) = Flux2%gas(msl) * BPCF%of(msl)
            else
                Flux3%gas(msl) = error
            end if
        else
            Flux3%gas(msl) = Flux2%gas(msl)
        end if
    end do

    !> Potential temperature
    !> If condition fails, previous value (from Fluxes0) holds for z/Ambient%L
    if (lEx%Pa > 0d0) then
        Tp = lEx%Ta * (1d5 / lEx%Pa)**(.286d0)
    else
        Tp = error
    end if

    !> Momentum flux and friction velocity
    Flux2%tau = Flux1%tau
    Flux3%tau = Flux1%tau
    Flux2%ustar = Flux1%ustar
    Flux3%ustar = Flux1%ustar

    !> Monin-Obukhov length (L = - (Tp^ /(k*g))*(ustar**3/(w'Tp')^ in m)
    if (Flux3%H /= 0d0 .and. Flux3%H /= error .and. &
        lEx%RhoCp > 0d0 .and. Flux3%ustar >= 0d0 .and. Tp > 0d0) then
        lEx%L = -Tp * (Flux3%ustar**3) / (vk * g * Flux3%H / lEx%RhoCp)
    else
        lEx%L = error
    end if

    !> Monin-Obukhov stability parameter (zL = z/L)
    !> If condition fails, previous value (from Fluxes0) holds
    if (lEx%L /= 0d0 .and. lEx%L /= error) &
        lEx%zL = (lEx%instr(sonic)%height - lEx%disp_height) / lEx%L

    !> scale temperature(T*)
    !> If condition fails, previous value (from Fluxes0) holds
    if (Flux3%ustar > 0d0 .and. Flux3%H /= error .and. lEx%RhoCp > 0d0) &
        lEx%Tstar = Flux3%H / (lEx%RhoCp * Flux3%ustar)

    !> Bowen ration (Bowen, 1926, Phyis Rev)
    if (Flux3%LE /= 0d0 .and. Flux3%LE /= error .and. Flux3%H /= error) then
        lEx%Bowen = Flux3%H / Flux3%LE
    else
        lEx%Bowen = error
    end if
contains

    !***********************************************************************
    !> Level 2 flux of one gas: the WPL / density correction.
    !>
    !> Mirrors Level2GasFlux in src_rp/fluxes23_rp.f90, reading from the ex
    !> record rather than the in-memory statistics. It replaces three
    !> near-identical blocks which had drifted; both differences are resolved
    !> in favour of the safe form, as agreed for the RP twin:
    !>
    !>  - the cell-volume conversion (lEx%Vcell / lEx%Va) is applied in every
    !>    closed-path branch, not all but one;
    !>  - the cell pressure is guarded with `> 0d0` rather than `/= error`,
    !>    which admitted zero and then divided by it.
    !>
    !> The original cascades enumerated all eight subsets of the (E, T, P)
    !> terms; accumulating each onto the base when its own inputs are present
    !> is equivalent, and keeps the floating-point association identical.
    subroutine Level2GasFlux(gas, sigma_gas, rhow_gas)
        implicit none
        integer, intent(in) :: gas
        real(kind = dbl), intent(in) :: sigma_gas
        real(kind = dbl), intent(in) :: rhow_gas
        real(kind = dbl) :: wpl

        if (Flux1%gas(gas) == error) then
            Flux2%gas(gas) = error
            return
        end if

        if (lEx%gas_instr(gas)%path_type == 'closed') then
            !> Closed path, after Ibrom et al. (2007) Tellus eq. 3a, with the
            !> H contribution from WPL24.
            select case (lEx%measure_type(gas))
                case ('mixing_ratio')
                    Flux2%gas(gas) = Flux1%gas(gas)
                    return
                case ('molar_density')
                    if (lEx%Va <= 0d0) then
                        Flux2%gas(gas) = error
                        return
                    end if
                    wpl = Flux1%gas(gas) * lEx%Vcell(gas) / lEx%Va
                case ('mole_fraction')
                    wpl = Flux1%gas(gas)
                case default
                    Flux2%gas(gas) = error
                    return
            end select

            if (sigma_gas >= 0d0 .and. lEx%Va > 0d0 .and. lEx%chi(gas) > 0d0) then
                !> Effect of the water vapour flux
                if (Flux3%E_gas(gas) /= error .and. rhow_gas > 0d0) &
                    wpl = wpl + Flux3%E_gas(gas) * mu * sigma_gas / rhow_gas &
                        * lEx%chi(gas) / lEx%Va
                !> Effect of cell temperature.
                !>
                !> From this gas's own analyser, and in SI. lEx%Tcell is
                !> instrument 1's *and* passes through the writer's degC gain
                !> without the reader inverting it, so this term used to divide
                !> by a temperature of about 27 instead of about 300.
                if (Flux3%Hi_gas(gas) /= error .and. lEx%RhoCp > 0d0 &
                    .and. lEx%Tcell_at(gas) > 0d0) &
                    wpl = wpl + (1d0 + mu * sigma_gas) * Flux3%Hi_gas(gas) &
                        / (lEx%RhoCp * lEx%Tcell_at(gas)) * lEx%chi(gas) / lEx%Va
                !> Effect of cell pressure, from this gas's own analyser. Same
                !> unit trap: lEx%Pcell carries the writer's kPa gain.
                if (lEx%cov_w_pcell(gas) /= error .and. lEx%Pcell_at(gas) > 0d0) &
                    wpl = wpl - (1d0 + mu * sigma_gas) * lEx%cov_w_pcell(gas) &
                        / (lEx%Pcell_at(gas)) * lEx%chi(gas) / lEx%Va
            end if
            Flux2%gas(gas) = wpl
        else
            !> Open path, after e.g. Burba et al. (2008, GCB, eq. 1)
            wpl = Flux1%gas(gas)
            if (Flux3%E /= error .and. lEx%RHO%d > 0d0 .and. sigma_gas /= error) &
                wpl = wpl + mu * Flux3%E * lEx%d(gas) * 1d3 &
                    / ((1d0 + mu * sigma_gas) * lEx%RHO%d)
            if (Flux3%H /= error .and. lEx%RhoCp > 0d0 .and. lEx%Ta > 0d0) &
                wpl = wpl + (Flux3%H + lEx%Burba%h_top + lEx%Burba%h_bot &
                    + lEx%Burba%h_spar) * lEx%d(gas) * 1d3 / (lEx%RhoCp * lEx%Ta)
            Flux2%gas(gas) = wpl
        end if

        if (.not. lEx%var_present(gas)) Flux2%gas(gas) = error
    end subroutine Level2GasFlux

    !> Water vapour terms to use when correcting `gas`.
    !>
    !> RP resolves each gas's moisture reference and writes the resulting
    !> terms into the ex file; they are read back per gas slot. Where a gas has
    !> no resolved reference the values are `error` and the single global
    !> sigma / RHO%w apply, which is the single-analyser case and leaves those
    !> projects unchanged.
    subroutine MoistTerms(gas, sigma_out, rhow_out)
        implicit none
        integer, intent(in) :: gas
        real(kind = dbl), intent(out) :: sigma_out
        real(kind = dbl), intent(out) :: rhow_out

        sigma_out = lEx%sigma
        rhow_out  = lEx%RHO%w

        if (gas < firstGas .or. gas > lastGas) return

        !> Only override where the referenced H2O actually yielded values; a
        !> partial record must not silently zero the correction.
        if (lEx%sigma_at(gas) /= error) sigma_out = lEx%sigma_at(gas)
        if (lEx%rhow_at(gas)  /= error) rhow_out  = lEx%rhow_at(gas)
    end subroutine MoistTerms

end subroutine Fluxes23
