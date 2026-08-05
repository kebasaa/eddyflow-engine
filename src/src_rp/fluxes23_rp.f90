!***************************************************************************
! fluxes23_rp.f90
! ---------------
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
! \todo        Merge with corresponding FCC sub into a common one
!***************************************************************************
subroutine Fluxes23_rp()
    use m_rp_global_var
    implicit none
    integer, external :: cellPressureSlot
    !> local variables
    real(kind = dbl) :: Tp
    real(kind = dbl) :: E_nowpl
    !> Humidity terms for the WPL correction of each gas, taken from the H2O
    !> measurement that gas names (E2Col(gas)%moist_ref). With one H2O these
    !> equal the global Ambient%sigma / RHO%w, so a single-analyser site is
    !> unaffected; with two analysers each gas is corrected with its own.
    real(kind = dbl) :: sigma_g, rhow_g
    integer :: msl
    integer :: wsl
    integer :: wsx
    integer :: gas
    include '../src_common/interfaces_1.inc'

    !> Water's own slot, resolved from the records. The terms below are about
    !> water vapour itself - latent heat, evapotranspiration, and the flux of
    !> water the other gases' dilution correction is built on - so they follow
    !> the project's primary water record rather than the h2o slot constant,
    !> which is record two and holds water only by convention.
    wsl = PrimaryWaterSlot()
    !> An always-in-bounds stand-in for wsl inside guard expressions.
    !>
    !> Fortran does not mandate short-circuit `.and.`, so
    !> `wsl >= firstGas .and. X(wsl)` still evaluates X(wsl) - and with no
    !> hygrometer wsl is 0, which is out of bounds. The wsl >= firstGas test
    !> still decides the outcome; wsx only keeps the subscript legal while it
    !> is being decided.
    wsx = max(wsl, firstGas)

    write(*,'(a)', advance = 'no') '  Calculating fluxes Level 2 and 3..'

    Flux2 = errFlux
    Flux3 = errFlux

    !> Level 2 end 3 internal sensible heat, do nothing
    !> A pass-through, so it carries the whole gas block. Named for four
    !> slots it dropped every gas past the fourth record - Fluxes0_rp computes
    !> Hi_gas per configured gas and the FLUXNET file writes an H_CELL_*
    !> column for each, but only four survived the chain.
    Flux2%Hi_gas(firstGas:lastGas) = Flux1%Hi_gas(firstGas:lastGas)
    Flux3%Hi_gas(firstGas:lastGas) = Flux2%Hi_gas(firstGas:lastGas)

    !> Level 2 evapotranspiration WPL corrected ,including Burba if the case
    if (EddyFlowProj%wpl) then
        if (wsl < firstGas) then
            !> No hygrometer, so there is no evapotranspiration to correct.
            !>
            !> Tested on its own rather than folded into the open-path test
            !> below: Fortran does not mandate short-circuit evaluation of
            !> `.and.`, so `wsl >= firstGas .and. E2Col(wsl)%...` is not a
            !> guard, and the closed-path arm reached E2Col(0) unconditionally.
            Flux2%E = error
        elseif (E2Col(wsl)%Instr%path_type == 'open') then
            !> Open-path uses Webb et al. (1980)
            !> Note that Burba terms are forced to zero
            !> if analyzer is /= LI-7500
            if (Ambient%RhoCp > 0d0 .and. Ambient%Ta > 0d0 &
                .and. Flux1%E /= error .and. Flux1%H /= error &
                .and. Ambient%sigma /= error) then
                Flux2%E = (1d0 + mu * Ambient%sigma) * Flux1%E &
                    + (1d0 + mu * Ambient%sigma) &
                        * (Flux1%H + Burba%h_top + Burba%h_bot + Burba%h_spar)&
                        * RHO%w / (Ambient%RhoCp * Ambient%Ta)
            else
                Flux2%E = error
            end if
        else
            !> Closed-path uses Ibrom et al. (2007) if conversion to mixing
            !> ratio did not already occur (which implies that some variables
            !> were missing)
            select case(E2Col(wsl)%measure_type)
                case ('molar_density', 'mole_fraction')
                    if (Flux1%E /= error .and. Ambient%sigma /= error &
                        .and. E2Col(wsl)%Va > 0d0 .and. Ambient%Va > 0d0) then

                        if (Flux1%Hi_gas(wsl) /= error &
                            .and. Stats%cov(w, pi) /= error) then
                            !> Complete formulation, should actually never be
                            !> used cause conversion to mixing ratio should have
                            !> already happened if everything is available
                            Flux2%E = (1d0 + mu * Ambient%sigma) * Flux1%E &
                                * E2Col(wsl)%Va / Ambient%Va &
                                + (1d0 + mu * Ambient%sigma) * Flux1%Hi_gas(wsl) &
                                * RHO%w / (Ambient%RhoCp * Ambient%Tcell) &
                                - (1d0 + mu * Ambient%sigma) * Stats%cov(w, pi) &
                                * RHO%w / (Ambient%Pcell)

                        elseif (Flux1%Hi_gas(wsl) /= error) then
                            !> Correct only for effect of T
                            Flux2%E = (1d0 + mu * Ambient%sigma) * Flux1%E &
                                * E2Col(wsl)%Va / Ambient%Va &
                                + (1d0 + mu * Ambient%sigma) * Flux1%Hi_gas(wsl) &
                                * RHO%w / (Ambient%RhoCp * Ambient%Tcell)

                        elseif (Stats%cov(w, pi)  /= error) then
                            !> Correct only for effect of P
                            Flux2%E = (1d0 + mu * Ambient%sigma) * Flux1%E &
                                * E2Col(wsl)%Va / Ambient%Va &
                                - (1d0 + mu * Ambient%sigma) * Stats%cov(w, pi) &
                                * RHO%w / (Ambient%Pcell)
                        else
                            !> Can't correct for T and P
                            Flux2%E = Flux1%E * E2Col(wsl)%Va / Ambient%Va
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
    !>
    !> `wsl` is a slot only when a hygrometer exists; with none it is 0 and
    !> Flux2%gas(0) writes past the low end of the array.
    if (wsl < firstGas) then
        Flux2%LE  = error
        Flux2%ET  = error
    elseif (Flux2%E /= error) then
        Flux2%gas(wsl) = Flux2%E * 1d3 / MW_H2O
        Flux2%ET = Flux2%gas(wsl) * h2o_to_ET
        if (Ambient%lambda /= error) then
            Flux2%LE = Flux2%E * Ambient%lambda
        else
            Flux2%LE = error
        end if
    else
        Flux2%gas(wsl) = error
        Flux2%LE  = error
        Flux2%ET  = error
    end if

    !> Level 2 evapotranspiration fluxes with H2O covariances
    !> at time-lags of other scalars. Do nothing, WPL is deleterious here
    Flux2%E_gas(firstGas:lastGas) = Flux1%E_gas(firstGas:lastGas)

    !> Level 2 Sensible heat
    if (E2Col(ts)%instr%category == 'sonic') then
        !> Corrected for humidity, after Van Dyjk et al. (2004) eq. 3.53
        !> revising Schotanus et al. (1983)
        if (Flux1%H /= error) then
            if(Flux0%E /= error .and. Stats%Cov(w, ts) /= error &
                .and. RHO%a > 0d0 .and. Ambient%Q >= 0d0 &
                .and. Ambient%RhoCp > 0d0 .and. Ambient%alpha /= error) then
                Flux2%H = Flux1%H &
                    - Ambient%RhoCp * Ambient%alpha * Stats%Mean(ts) * Flux0%E / RHO%a &
                    - Ambient%RhoCp * Ambient%alpha * Ambient%Q * Stats%Cov(w, ts)
                    !> alternative
                    !- Ambient%RhoCp * Ambient%alpha * Ambient%Ta * Flux0%E / RHO%a
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
    !if(Ambient%Tmap /= error .and. Flux2%H /= error) then
        !Flux2%H  = Flux2%H * Ambient%Tmap
    !end if

    !> Level 3 sensible heat, spectral corrected
    if(Flux2%H /= error .and. BPCF%of(w_ts) /= error) then
        Flux3%H = Flux2%H * BPCF%of(w_ts)
    else
        Flux3%H = error
    end if

    !> Level 3 for evapotranspiration: for open path, WPL again with corrected H
    !> Starts again from Level 1 of E, Level 2 was only used to
    !> calculate H Level 3.
    if(EddyFlowProj%wpl .and. wsl >= firstGas &
        .and. E2Col(wsx)%Instr%path_type == 'open') then
        if (Ambient%RhoCp > 0d0 .and. Ambient%Ta > 0d0 .and. Flux1%E /= error &
            .and. Flux1%H /= error .and. Ambient%sigma /= error) then
            Flux3%E = (1d0 + mu * Ambient%sigma) * Flux1%E &
                + (1d0 + mu * Ambient%sigma) &
                * (Flux3%H + Burba%h_top + Burba%h_bot + Burba%h_spar)&
                * RHO%w / (Ambient%RhoCp * Ambient%Ta)
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
    if (wsl < firstGas) then
        Flux3%LE  = error
        Flux3%ET  = error
    elseif (Flux3%E /= error) then
        Flux3%gas(wsl) = Flux3%E * 1d3 / MW_H2O
        Flux3%ET = Flux3%gas(wsl) * h2o_to_ET
        if (Ambient%lambda /= error) then
            Flux3%LE = Flux3%E * Ambient%lambda
        else
            Flux3%LE = error
        end if
    else
        Flux3%gas(wsl) = error
        Flux3%LE  = error
        Flux3%ET  = error
    end if

    !> Calculate E_nowpl for closed and open path systems
    if (wsl >= firstGas .and. E2Col(wsx)%Instr%path_type == 'closed') then
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
        if (E2Col(wsl)%Instr%path_type == 'closed') then
            !> Level 3, spectral correction
            Flux3%gas(wsl) = Flux3%gas(wsl) * BPCF%of(wsl)
            Flux3%E   = Flux3%E   * BPCF%of(wsl)
            Flux3%LE  = Flux3%LE  * BPCF%of(wsl)
            Flux3%ET  = Flux3%ET  * BPCF%of(wsl)
        end if
    end if

    if (wsl < firstGas) then
        !> No water record at all: nothing water-derived is performed.
        Flux2%E = error; Flux2%LE = error; Flux2%ET = error
        Flux3%E = error; Flux3%LE = error; Flux3%ET = error
    elseif (.not. E2Col(wsl)%present) then
        Flux3%gas(wsl) = error
        Flux3%E   = error
        Flux3%LE  = error
        Flux3%ET  = error
    end if

    !> Level 2 other gases.
    !>
    !> One loop over the configured gases, replacing three near-duplicate
    !> blocks. Those blocks enumerated every subset of the (E, T, P) WPL terms
    !> as an if/elseif cascade; accumulating each term when its inputs are
    !> available is equivalent, and does not have to be written once per gas.
    do msl = firstGas, lastGas
        !> The primary water's own flux is the evapotranspiration handled
        !> above, not a WPL-corrected trace gas flux. Keyed on the resolved
        !> slot, not the h2o constant: on a project that declares its water
        !> anywhere but record two, that constant skipped a trace gas from
        !> this loop - leaving its Level 2 and 3 fluxes unset - and sent the
        !> real water through the trace-gas path instead.
        !>
        !> A *second* hygrometer is deliberately still treated as a trace gas
        !> here: only one evapotranspiration is computed above, so skipping
        !> it would leave its flux unreported.
        if (msl == wsl) cycle
        if (.not. E2Col(msl)%present) then
            Flux2%gas(msl) = error
            cycle
        end if
        call MoistTerms(msl, sigma_g, rhow_g)
        call Level2GasFlux(msl, sigma_g, rhow_g)
    end do

    !> Level 3 other gases. For closed path apply the spectral correction
    !> now (e.g. Ibrom et al. 2007); for open path it is already included.
    !> BPCF%of is indexed by the w_* covariance labels, which carry the same
    !> numbering as the gas slots, so the slot indexes it directly.
    do msl = firstGas, lastGas
        if (msl == wsl) cycle
        if (E2Col(msl)%Instr%path_type == 'closed' .and. Flux2%gas(msl) /= error) then
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
    if (Stats%Pr > 0d0) then
        Tp = Ambient%Ta * (1d5 / Stats%Pr)**(.286d0)
    else
        Tp = error
    end if

    !> Monin-Obukhov length (L = - (Tp^ /(k*g))*(ustar**3/(w'Tp') in m)
    if (Flux3%H /= 0d0 .and. Flux3%H /= error .and. &
        Ambient%RhoCp > 0d0 .and. Ambient%us >= 0d0 .and. Tp > error) then
        Ambient%L = -Tp * (Ambient%us**3) / (vk * g * Flux3%H / Ambient%RhoCp)
    else
        Ambient%L = error
    end if

    !> Monin-Obukhov stability parameter (zL = z/L)
    !> If condition fails, previous value (from Fluxes0) holds
    if (Ambient%L /= 0d0 .and. Ambient%L /= error) &
        Ambient%zL = (E2Col(u)%Instr%height - Metadata%d) / Ambient%L

    !> scale temperature(T*)
    !> If condition fails, previous value (from Fluxes0) holds
    if (Ambient%us > 0d0 .and. Flux3%H /= error .and. Ambient%RhoCp > 0d0) &
        Ambient%Ts = - Flux3%H / (Ambient%RhoCp * Ambient%us)

    !> Bowen ration (Bowen, 1926, Phyis Rev)
    if (Flux3%LE /= 0d0 .and. Flux3%LE /= error .and. Flux3%H /= error) then
        Ambient%Bowen = Flux3%H / Flux3%LE
    else
        Ambient%Bowen = error
    end if

    !> Momentum flux
    Flux2%tau = Flux1%tau
    Flux3%tau = Flux1%tau
    Flux2%ustar = Flux1%ustar
    Flux3%ustar = Flux1%ustar

    !> If fluxes are error, set also time-lags to error, just for clarity
    !>
    !> Every configured gas. Four named slots left a fifth gas reporting a
    !> plausible time-lag beside an errored flux, and the column is written per
    !> configured gas at both ends. The water entry also read Flux2%gas(wsl)
    !> unguarded, so a project with no water at all indexed element zero.
    do gas = firstGas, lastGas
        if (gas - firstGas + 1 > min(EddyFlowProj%gas_num, MaxNumGases)) exit
        if (Flux2%gas(gas) == error) Essentials%used_timelag(gas) = error
    end do

    write(*,'(a)')   ' Done.'

contains

!***************************************************************************
!> Level 2 flux of one gas: the WPL / density correction.
!>
!> Replaces three near-identical blocks (co2, ch4, gas4) that between them had
!> drifted apart. Two differences were resolved in favour of the safe form:
!>
!>  - the ch4 block was missing the cell-volume conversion
!>    (`E2Col(gas)%Va / Ambient%Va`) in one closed-path branch, where co2 and
!>    gas4 applied it. It is applied uniformly here.
!>  - co2 and ch4 guarded the cell pressure with `Ambient%Pcell /= error`,
!>    which admits zero and then divides by it; gas4 used `> 0d0`. The
!>    positive test is used throughout.
!>
!> The original cascades enumerated all eight subsets of the (E, T, P) terms.
!> Adding each term when its own inputs are available is equivalent, and is
!> what lets this be written once rather than once per gas.
subroutine Level2GasFlux(gas, sigma_gas, rhow_gas)
    implicit none
    integer, intent(in) :: gas
    real(kind = dbl), intent(in) :: sigma_gas
    real(kind = dbl), intent(in) :: rhow_gas
    real(kind = dbl) :: base
    real(kind = dbl) :: wpl
    !> Cell conditions of this gas's own analyser, falling back to the global
    !> pair when it has none - which is the single-analyser case.
    real(kind = dbl) :: Tcell_g
    real(kind = dbl) :: Pcell_g

    if (Flux1%gas(gas) == error) then
        Flux2%gas(gas) = error
        return
    end if

    Tcell_g = Ambient%Tcell
    Pcell_g = Ambient%Pcell
    if (gas >= firstGas .and. gas <= lastGas) then
        if (Ambient%Tcell_at(gas) /= error) Tcell_g = Ambient%Tcell_at(gas)
        if (Ambient%Pcell_at(gas) /= error) Pcell_g = Ambient%Pcell_at(gas)
    end if

    if (E2Col(gas)%Instr%path_type == 'closed') then
        !> Closed path, after Ibrom et al. (2007) Tellus eq. 3a, with the H
        !> contribution from WPL24.
        select case (E2Col(gas)%measure_type)
            case ('mixing_ratio')
                !> Already a mixing ratio: no density correction applies.
                Flux2%gas(gas) = Flux1%gas(gas)
                return
            case ('molar_density')
                if (Ambient%Va <= 0d0) then
                    Flux2%gas(gas) = error
                    return
                end if
                base = Flux1%gas(gas) * E2Col(gas)%Va / Ambient%Va
            case ('mole_fraction')
                base = Flux1%gas(gas)
            case default
                Flux2%gas(gas) = error
                return
        end select

        !> Accumulated onto the base in the order the original cascade wrote
        !> the terms, so the floating-point association matches and the result
        !> is bit-identical to the code this replaces.
        wpl = base
        if (sigma_gas >= 0d0 .and. Ambient%Va > 0d0 .and. Stats%chi(gas) > 0d0) then
            !> Effect of the water vapour flux
            if (Flux3%E_gas(gas) /= error .and. rhow_gas > 0d0) &
                wpl = wpl + Flux3%E_gas(gas) * mu * sigma_gas / rhow_gas &
                    * Stats%chi(gas) / Ambient%Va
            !> Effect of cell temperature
            if (Flux3%Hi_gas(gas) /= error .and. Ambient%RhoCp > 0d0 &
                .and. Tcell_g > 0d0) &
                wpl = wpl + (1d0 + mu * sigma_gas) * Flux3%Hi_gas(gas) &
                    / (Ambient%RhoCp * Tcell_g) &
                    * Stats%chi(gas) / Ambient%Va
            !> Effect of cell pressure
            !> The pressure covariance is read from this gas's own cell
            !> pressure slot too, not from the first instrument's.
            if (Stats%cov(w, cellPressureSlot(gas)) /= error .and. Pcell_g > 0d0) &
                wpl = wpl - (1d0 + mu * sigma_gas) &
                    * Stats%cov(w, cellPressureSlot(gas)) &
                    / (Pcell_g) &
                    * Stats%chi(gas) / Ambient%Va
        end if
        Flux2%gas(gas) = wpl
    else
        !> Open path, after e.g. Burba et al. (2008, GCB, eq. 1)
        wpl = Flux1%gas(gas)
        if (Flux3%E /= error .and. RHO%d > 0d0 .and. sigma_gas /= error) &
            wpl = wpl + mu * Flux3%E * Stats%d(gas) * 1d3 &
                / ((1d0 + mu * sigma_gas) * RHO%d)
        if (Flux3%H /= error .and. Ambient%RhoCp > 0d0 .and. Ambient%Ta > 0d0) &
            wpl = wpl + (Flux3%H + Burba%h_top + Burba%h_bot + Burba%h_spar) &
                * Stats%d(gas) * 1d3 / (Ambient%RhoCp * Ambient%Ta)
        Flux2%gas(gas) = wpl
    end if

    if (.not. E2Col(gas)%present) Flux2%gas(gas) = error
end subroutine Level2GasFlux


!***************************************************************************
!> Humidity terms used to WPL-correct one gas.
!>
!> A gas names the H2O measurement it should be corrected with
!> (E2Col(gas)%moist_ref, resolved in DefineE2Set). When that reference is
!> unset - a legacy project, or a gas with no H2O available - fall back to the
!> single global pair, which is exactly the historical behaviour.
!>
!> With one H2O record the reference resolves to that same slot and these
!> return the global values, so single-analyser sites are bit-identical.
!***************************************************************************
subroutine MoistTerms(gas, sigma_out, rhow_out)
    implicit none
    integer, intent(in) :: gas
    real(kind = dbl), intent(out) :: sigma_out
    real(kind = dbl), intent(out) :: rhow_out
    integer :: msl

    sigma_out = Ambient%sigma
    rhow_out  = RHO%w

    msl = E2Col(gas)%moist_ref
    if (msl < firstGas .or. msl > lastGas) return
    if (.not. E2Col(msl)%present) return

    !> Only override when the referenced H2O actually yielded values; a
    !> partial record must not silently zero the correction.
    if (Ambient%sigma_at(msl) /= error) sigma_out = Ambient%sigma_at(msl)
    if (RHO%w_at(msl) /= error)         rhow_out  = RHO%w_at(msl)
end subroutine MoistTerms

end subroutine Fluxes23_rp
