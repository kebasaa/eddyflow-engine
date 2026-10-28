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
    !> Every hygrometer, and each one's fully corrected water flux and its
    !> flux without WPL - what a trace gas's correction reads from the
    !> hygrometer its record names. See the RP twin.
    integer :: w_slots(GHGNumVar)
    character(8) :: w_tags(GHGNumVar)
    integer :: n_w_slots
    integer :: iw
    logical :: isHygro(GHGNumVar)
    real(kind = dbl) :: E3w(GHGNumVar)
    real(kind = dbl) :: Enowplw(GHGNumVar)
    real(kind = dbl) :: E3_ref, Enowpl_ref
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

    !> Sensible heat comes first. It needs only the level 0 water flux, and
    !> every hygrometer's open-path WPL at level 3 needs it.

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

    !> Every hygrometer, and each one's sensible heat on its own air.
    call WaterOutSlots(w_slots, w_tags, n_w_slots)
    isHygro = .false.
    do iw = 1, n_w_slots
        isHygro(w_slots(iw)) = .true.
    end do
    call PerHygrometerHeat()

    !> Level 2 and 3 evapotranspiration of the designated hygrometer, WPL
    !> corrected after its own path type, including Burba if the case.
    if (wsl < firstGas) then
        !> No hygrometer, so there is no evapotranspiration to correct.
        Flux2%E = error
        Flux3%E = error
        E_nowpl = Flux1%E
    else
        call HygrometerWater(wsl, Flux1%E, lEx%sigma, lEx%RHO%w, lEx%RhoCp, &
            Flux3%H, Flux2%E, Flux3%E, E_nowpl)
    end if

    !> Level 2 h2o and latent heat flux
    if (wsl < firstGas) then
        Flux2%LE  = error
        Flux2%ET  = error
    elseif (Flux2%E /= error) then
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

    !> Level 2 and 3 evapotranspiration fluxes with H2O covariances at
    !> time-lags of other scalars. Do nothing, WPL is deleterious here.
    !> Another pass-through, so it carries the whole gas block; the water
    !> slot has no entry by construction and copying it costs nothing.
    Flux2%E_gas(firstGas:lastGas) = Flux1%E_gas(firstGas:lastGas)
    Flux3%E_gas(firstGas:lastGas) = Flux2%E_gas(firstGas:lastGas)

    !> Level 3 h2o flux and latent heat flux
    if (wsl < firstGas) then
        Flux3%LE  = error
        Flux3%ET  = error
    elseif (Flux3%E /= error) then
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

    !> Apply spectral correction to h2o and E/LE Level 3 fluxes for closed path
    if (wsl >= firstGas .and. Flux3%E /= error) then
        if (lEx%gas_instr(wsl)%path_type == 'closed') then
            !> No correction factor means the corrected flux is unavailable -
            !> not the flux times the error sentinel, which is what the
            !> unguarded multiply produced.
            if (BPCF%of(wsl) /= error) then
                !> Level 3, spectral correction
                Flux3%gas(wsl) = Flux3%gas(wsl) * BPCF%of(wsl)
                Flux3%E   = Flux3%E   * BPCF%of(wsl)
                Flux3%LE  = Flux3%LE  * BPCF%of(wsl)
                Flux3%ET  = Flux3%ET  * BPCF%of(wsl)
            else
                Flux3%gas(wsl) = error
                Flux3%E   = error
                Flux3%LE  = error
                Flux3%ET  = error
            end if
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

    !> Every other hygrometer: the same correction, after its own path type and
    !> on its own analyser and air. See the RP twin for what they had before.
    E3w = error
    Enowplw = error
    if (wsl >= firstGas .and. wsl <= lastGas) then
        E3w(wsl) = Flux3%E
        Enowplw(wsl) = E_nowpl
    end if
    do iw = 1, n_w_slots
        if (w_slots(iw) == wsl) cycle
        call OtherHygrometerWater(w_slots(iw))
    end do

    !> Level 2 other gases.
    !>
    !> One loop over the gases, replacing three near-duplicate blocks that had
    !> drifted apart; see Level2GasFlux for the two differences resolved.
    !>
    !> Runs the full gas range: lEx%gas_instr is indexed by gas slot and the
    !> ex file now carries an analyser for every configured gas.
    do msl = firstGas, lastGas
        !> Every hygrometer's own flux is its evapotranspiration, above.
        if (msl == wsl .or. isHygro(msl)) cycle
        call MoistTerms(msl, sigma_g, rhow_g)
        call WaterFluxFor(msl, E3_ref, Enowpl_ref)
        call Level2GasFlux(msl, sigma_g, rhow_g, E3_ref, Enowpl_ref)
    end do

    !> Level 3 other gases. For closed path apply the spectral correction now
    !> (e.g. Ibrom et al. 2007); for open path it is already included.
    do msl = firstGas, lastGas
        if (msl == wsl .or. isHygro(msl)) cycle
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

    !> Eugster and Senn (1995) discard a gas when less than 40% of its flux
    !> was measurable. Done last, so the other gases' WPL terms have already
    !> used this water flux; only the discarded gas's own fluxes go.
    do msl = firstGas, lastGas
        if (.not. BPCF_discard(msl)) cycle
        Flux3%gas(msl) = error
        if (msl == wsl) then
            Flux3%E  = error
            Flux3%LE = error
            Flux3%ET = error
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

    !> scale temperature(T*), see e.g. Foken and Wichura (1996)
    !>
    !> The minus is the definition - T* = -(w'T')/u* - not a convention. With
    !> a downward heat flux w'T' is negative and T* positive; without it this
    !> line returned the negative of the quantity its own header names, on
    !> every run that reached FCC, and overwrote the correctly signed value
    !> RP had written into the essentials record. Both RP sites (Fluxes0_rp,
    !> Fluxes23_rp) always carried it.
    !>
    !> If the condition fails, RP's value from the essentials record holds:
    !> FCC has no Fluxes0 of its own to fall back on.
    if (Flux3%ustar > 0d0 .and. Flux3%H /= error .and. lEx%RhoCp > 0d0) &
        lEx%Tstar = - Flux3%H / (lEx%RhoCp * Flux3%ustar)

    !> Bowen ration (Bowen, 1926, Phyis Rev)
    if (Flux3%LE /= 0d0 .and. Flux3%LE /= error .and. Flux3%H /= error) then
        lEx%Bowen = Flux3%H / Flux3%LE
    else
        lEx%Bowen = error
    end if

    !> The same quantities, once per hygrometer.
    !>
    !> Everything above computes one latent heat flux, one evapotranspiration,
    !> one sensible heat flux and one stability, from the designated
    !> hygrometer. A site fielding two produced a second water flux and nothing
    !> derived from it - so the two hygrometers could not be compared on any
    !> quantity a reader actually uses, and the disagreement that motivates
    !> fielding two was invisible.
    !>
    !> Each entry is the block above replayed on that hygrometer's own air, as
    !> RP resolved it into the moisture terms this record carries. The
    !> designated hygrometer's entry is then *assigned to* the scalars rather
    !> than recomputed, which is what keeps the bare columns and the numbered
    !> column describing the same instrument from ever disagreeing.
    call PerHygrometerFluxes()
contains

    !***********************************************************************
    !> Sensible heat, once per hygrometer.
    !>
    !> Split from the rest of the per-hygrometer set because it is needed
    !> earlier: each hygrometer's open-path WPL redoes its level 3
    !> evapotranspiration with its own level 3 H, as the designated hygrometer
    !> has always done with Flux3%H. Spectrally corrected, as Flux3%H is -
    !> Flux1%H carries no correction factor, so the numbered H of a second
    !> hygrometer was the uncorrected flux beside a corrected one.
    subroutine PerHygrometerHeat()
        implicit none
        !> `wslot`, not `w`: w is the module's wind-component index and a loop
        !> variable of that name would shadow it inside this routine, silently
        !> retargeting lEx%cov_w(ts) and every other use below.
        integer :: iw, wslot
        real(kind = dbl) :: rhocp_w, rhoa_w, q_w, e0_w

        Flux3%H_at = error

        do iw = 1, n_w_slots
            wslot = w_slots(iw)
            if (.not. lEx%var_present(wslot)) cycle

            !> A hygrometer's own moisture reference is itself, so these are
            !> its own air rather than some other instrument's.
            rhocp_w = lEx%rhocp_at(wslot)
            rhoa_w  = lEx%rhoa_at(wslot)
            q_w     = lEx%q_at(wslot)
            if (rhocp_w == error) rhocp_w = lEx%RhoCp
            if (rhoa_w  == error) rhoa_w  = lEx%RHO%a
            if (q_w     == error) q_w     = lEx%Q

            !> The humidity correction of sensible heat takes the *level 0*
            !> water flux, which is what the scalar above uses - lEx%Flux0%E is
            !> exactly Flux0%gas(designated) in these units. Reaching for the
            !> corrected flux here instead would make the numbered columns
            !> disagree with the bare ones by more than the hygrometers do.
            e0_w = error
            if (lEx%Flux0%gas(wslot) /= error) &
                e0_w = lEx%Flux0%gas(wslot) * MW_H2O * 1d-3

            !> Sensible heat corrected for humidity, after Van Dijk et al.
            !> (2004) eq. 3.53 revising Schotanus et al. (1983) - the block at
            !> the top of this routine, on this hygrometer's air.
            if (Flux1%H /= error) then
                if (e0_w /= error .and. lEx%cov_w(ts) /= error &
                    .and. rhoa_w > 0d0 .and. q_w >= 0d0 &
                    .and. rhocp_w > 0d0) then
                    Flux3%H_at(wslot) = Flux1%H &
                        - rhocp_w * alpha * lEx%Ts * e0_w / rhoa_w &
                        - rhocp_w * alpha * q_w * lEx%cov_w(ts)
                else
                    Flux3%H_at(wslot) = Flux1%H
                end if
                !> Level 3, spectral corrected, as the bare H.
                if (BPCF%of(w_ts) /= error) then
                    Flux3%H_at(wslot) = Flux3%H_at(wslot) * BPCF%of(w_ts)
                else
                    Flux3%H_at(wslot) = error
                end if
            end if
        end do
    end subroutine PerHygrometerHeat

    !***********************************************************************
    !> LE, ET, tau and the stability, one set per hygrometer.
    !>
    !> Called after the scalars are settled, and deliberately reads them: the
    !> water flux Flux3%gas(w) is already corrected after that hygrometer's
    !> own path type and spectrally adjusted, H_at is in place from
    !> PerHygrometerHeat, and Tp and ustar carry no humidity, so nothing here
    !> repeats work done above.
    subroutine PerHygrometerFluxes()
        implicit none
        !> `wslot`, not `w`: see PerHygrometerHeat.
        integer :: iw, wslot
        real(kind = dbl) :: rhocp_w, rhoa_w

        Flux3%E_at = error
        Flux3%LE_at = error
        Flux3%ET_at = error
        Flux3%tau_at = error
        Flux3%L_at = error
        Flux3%zL_at = error

        do iw = 1, n_w_slots
            wslot = w_slots(iw)
            if (.not. lEx%var_present(wslot)) cycle

            rhocp_w = lEx%rhocp_at(wslot)
            rhoa_w  = lEx%rhoa_at(wslot)
            if (rhocp_w == error) rhocp_w = lEx%RhoCp
            if (rhoa_w  == error) rhoa_w  = lEx%RHO%a

            !> Reported water flux: fully corrected, as LE and ET are.
            if (Flux3%gas(wslot) /= error) then
                Flux3%E_at(wslot) = Flux3%gas(wslot) * MW_H2O * 1d-3
                Flux3%ET_at(wslot) = Flux3%gas(wslot) * h2o_to_ET
                if (lEx%lambda /= error) &
                    Flux3%LE_at(wslot) = Flux3%E_at(wslot) * lEx%lambda
            end if

            !> Momentum, which follows the humidity only through air density.
            !> u* is not here: it comes from the wind covariances alone and is
            !> the same number for every hygrometer.
            if (rhoa_w > 0d0 .and. Flux3%ustar /= error) &
                Flux3%tau_at(wslot) = &
                    sign(rhoa_w * Flux3%ustar**2d0, Flux3%tau)

            !> Monin-Obukhov length and stability parameter
            if (Flux3%H_at(wslot) /= 0d0 .and. Flux3%H_at(wslot) /= error &
                .and. rhocp_w > 0d0 .and. Flux3%ustar >= 0d0 .and. Tp > 0d0) then
                Flux3%L_at(wslot) = -Tp * (Flux3%ustar**3) &
                    / (vk * g * Flux3%H_at(wslot) / rhocp_w)
                if (Flux3%L_at(wslot) /= 0d0) &
                    Flux3%zL_at(wslot) = &
                        (lEx%instr(sonic)%height - lEx%disp_height) &
                        / Flux3%L_at(wslot)
            end if
        end do

        !> The designated hygrometer's entry is *overwritten by* the scalars,
        !> not the other way round.
        !>
        !> Both directions would give the same numbers if this loop were a
        !> perfect replay, and it very nearly is. But the scalars have been
        !> through the storage chain, and this loop reproduces none of that -
        !> it does not need to, because every other hygrometer reaches it by
        !> the same route. Writing the loop's answer into the bare columns
        !> would move output that has been verified against v7.2.5; writing
        !> the scalar into the entry cannot.
        !>
        !> It also gives the row writer a single rule - walk WaterOutSlots and
        !> emit the `_at` entry - with no special case for the designated slot,
        !> which is one fewer place for header and row to disagree.
        if (wsl >= firstGas .and. wsl <= lastGas) then
            Flux3%E_at(wsl)   = Flux3%E
            Flux3%LE_at(wsl)  = Flux3%LE
            Flux3%ET_at(wsl)  = Flux3%ET
            Flux3%H_at(wsl)   = Flux3%H
            Flux3%tau_at(wsl) = Flux3%tau
            Flux3%L_at(wsl)   = lEx%L
            Flux3%zL_at(wsl)  = lEx%zL
        end if
    end subroutine PerHygrometerFluxes

    !***********************************************************************
    !> Level 2 and 3 evapotranspiration of one hygrometer, WPL corrected after
    !> its own path type. The RP twin, HygrometerWater in fluxes23_rp.f90,
    !> says why this is one routine.
    !>
    !> The caller supplies the air; the cell is the hygrometer's own record in
    !> the ex file, which is what the designated hygrometer always read. E3 is
    !> returned before the closed-path spectral correction, which the caller
    !> applies. With WPL off, E2 = E1, as EddyPro does.
    subroutine HygrometerWater(slot, E1, sigma, rhow, rhocp, H3, E2, E3, Enw)
        implicit none
        integer, intent(in) :: slot
        real(kind = dbl), intent(in) :: E1
        real(kind = dbl), intent(in) :: sigma, rhow, rhocp
        real(kind = dbl), intent(in) :: H3
        real(kind = dbl), intent(out) :: E2, E3, Enw

        E2 = error
        if (EddyFlowProj%wpl) then
            if (lEx%gas_instr(slot)%path_type == 'open') then
                if (rhocp > 0d0 .and. lEx%Ta > 0d0 &
                    .and. E1 /= error .and. Flux1%H /= error &
                    .and. sigma /= error) then
                        !> Open-path uses Webb et al. (1980)
                        !> Note that Burba terms are forced to zero
                        !> if analyzer is /= LI-7500
                        E2 = (1d0 + mu * sigma) * E1 &
                                + (1d0 + mu * sigma) &
                                * (Flux1%H + BurbaHeatFor(slot)) &
                                * rhow / (rhocp * lEx%Ta)
                end if
            else
                !> Closed-path uses Ibrom et al. (2007) if conversion to mixing
                !> ratio did not already occur (which implies that some
                !> variables were missing)
                select case(lEx%measure_type(slot))
                    case ('molar_density', 'mole_fraction')
                        if (E1 /= error .and. sigma /= error &
                            .and. lEx%Vcell(slot) > 0d0 .and. lEx%Va > 0d0) then

                            if (Flux1%Hi_gas(slot) /= error &
                                .and. lEx%cov_w_pcell(slot) /= error) then
                                !> Complete formulation, should actually never
                                !> be used cause conversion to mixing ratio
                                !> should have already happened if everything
                                !> is available
                                E2 = (1d0 + mu * sigma) * E1 &
                                    * lEx%Vcell(slot) / lEx%Va &
                                    + (1d0 + mu * sigma) * Flux1%Hi_gas(slot) &
                                    * rhow / (rhocp * lEx%Tcell_at(slot)) &
                                    - (1d0 + mu * sigma) * lEx%cov_w_pcell(slot) &
                                    * rhow / (lEx%Pcell_at(slot))

                            elseif (Flux1%Hi_gas(slot) /= error) then
                                !> Correct only for effect of T
                                E2 = (1d0 + mu * sigma) * E1 &
                                    * lEx%Vcell(slot) / lEx%Va &
                                    + (1d0 + mu * sigma) * Flux1%Hi_gas(slot) &
                                    * rhow / (rhocp * lEx%Tcell_at(slot))

                            elseif (lEx%cov_w_pcell(slot)  /= error) then
                                !> Correct only for effect of P
                                E2 = (1d0 + mu * sigma) * E1 &
                                    * lEx%Vcell(slot) / lEx%Va &
                                    - (1d0 + mu * sigma) * lEx%cov_w_pcell(slot) &
                                    * rhow / (lEx%Pcell_at(slot))
                            else
                                !> Can't correct for T and P
                                E2 = E1 * lEx%Vcell(slot) / lEx%Va
                            end if
                        end if
                    case ('mixing_ratio')
                        E2 = E1
                end select
            end if
        else
            !> If WPL should not be applied
            E2 = E1
        end if
        !> If Flux1 was error, then set also 2 to error
        if (E1 == error) E2 = error

        !> Level 3 for evapotranspiration: for open path, WPL again with
        !> corrected H. Starts again from Level 1 of E, Level 2 was only used
        !> to calculate H Level 3. A level 3 H that is error leaves level 2
        !> standing rather than entering the formula as -9999.
        if (EddyFlowProj%wpl .and. lEx%gas_instr(slot)%path_type == 'open') then
            if (rhocp > 0d0 .and. lEx%Ta > 0d0 .and. E1 /= error &
                .and. Flux1%H /= error .and. H3 /= error .and. sigma /= error) then
                E3 = (1d0 + mu * sigma) * E1 &
                    + (1d0 + mu * sigma) &
                    * (H3 + BurbaHeatFor(slot))&
                    * rhow / (rhocp * lEx%Ta)
            else
                E3 = E2
            end if
        else
            E3 = E2
        end if

        !> E without WPL, spectrally corrected, for the LI-7700 formulation
        if (lEx%gas_instr(slot)%path_type == 'closed') then
            if (E1 /= error .and. BPCF%of(slot) /= error) then
                Enw = E1 * BPCF%of(slot)
            elseif(E1 /= error) then
                Enw = E1
            else
                Enw = error
            end if
        else
            Enw = E1
        end if
    end subroutine HygrometerWater

    !***********************************************************************
    !> Level 2 and 3 water flux of a hygrometer other than the designated one:
    !> HygrometerWater on its own air and analyser, then the closed-path
    !> spectral correction, guarded. See the RP twin.
    subroutine OtherHygrometerWater(ws)
        implicit none
        integer, intent(in) :: ws
        real(kind = dbl) :: E1, E2, E3, Enw
        real(kind = dbl) :: sig, rw, rcp

        Flux2%gas(ws) = error
        Flux3%gas(ws) = error
        if (.not. lEx%var_present(ws)) return

        E1 = error
        if (Flux1%gas(ws) /= error) E1 = Flux1%gas(ws) * MW_H2O * 1d-3

        sig = lEx%sigma_at(ws)
        rw  = lEx%rhow_at(ws)
        rcp = lEx%rhocp_at(ws)
        if (sig == error) sig = lEx%sigma
        if (rw  == error) rw  = lEx%RHO%w
        if (rcp == error) rcp = lEx%RhoCp

        call HygrometerWater(ws, E1, sig, rw, rcp, Flux3%H_at(ws), E2, E3, Enw)

        if (E2 /= error) Flux2%gas(ws) = E2 * 1d3 / MW_H2O
        if (E3 /= error .and. lEx%gas_instr(ws)%path_type == 'closed') then
            if (BPCF%of(ws) /= error) then
                E3 = E3 * BPCF%of(ws)
            else
                E3 = error
            end if
        end if
        if (E3 /= error) Flux3%gas(ws) = E3 * 1d3 / MW_H2O

        E3w(ws) = E3
        Enowplw(ws) = Enw
    end subroutine OtherHygrometerWater

    !***********************************************************************
    !> The water flux a trace gas's WPL correction reads: that of the
    !> hygrometer its record names, as RP resolved it and wrote it into the ex
    !> file. Falls back to the designated hygrometer when that is not a fast
    !> hygrometer, or has no flux this period. See the RP twin.
    subroutine WaterFluxFor(gas, e3, enw)
        implicit none
        integer, intent(in) :: gas
        real(kind = dbl), intent(out) :: e3, enw
        integer :: ref

        e3 = Flux3%E
        enw = E_nowpl
        ref = lEx%moist_ref_at(gas)
        if (ref == wsl) return
        if (ref < firstGas .or. ref > lastGas) return
        if (.not. isHygro(ref)) return
        if (E3w(ref) == error) return
        e3 = E3w(ref)
        enw = Enowplw(ref)
    end subroutine WaterFluxFor

    !***********************************************************************
    !> The instrument-body heating terms, for the gas being corrected.
    !>
    !> Burba et al. (2008) is the LI-7500's own body warming the air in its
    !> path, so it belongs only to a gas that analyser measures.
    !> OverrideSettings switches the correction off for a site with no LI-7500,
    !> but that is site-wide: with an LI-7500 and an LI-7700 side by side it
    !> stays on, and the generic per-gas WPL was adding the LI-7500's heating
    !> to the methane flux. See the RP twin.
    real(kind = dbl) function BurbaHeatFor(gas)
        implicit none
        integer, intent(in) :: gas

        if (index(lEx%gas_instr(gas)%model, 'li7500') /= 0) then
            BurbaHeatFor = lEx%Burba%h_top + lEx%Burba%h_bot + lEx%Burba%h_spar
        else
            BurbaHeatFor = 0d0
        end if
    end function BurbaHeatFor

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
    !>
    !> `e3_ref` and `enw_ref` are the water flux of the hygrometer the gas's
    !> record names, from WaterFluxFor; and the WPL switch is honoured, as in
    !> EddyPro 7. Both are explained in the RP twin.
    subroutine Level2GasFlux(gas, sigma_gas, rhow_gas, e3_ref, enw_ref)
        implicit none
        integer, intent(in) :: gas
        real(kind = dbl), intent(in) :: sigma_gas
        real(kind = dbl), intent(in) :: rhow_gas
        real(kind = dbl), intent(in) :: e3_ref
        real(kind = dbl), intent(in) :: enw_ref
        real(kind = dbl) :: wpl
        real(kind = dbl) :: dens_to_chi

        if (Flux1%gas(gas) == error) then
            Flux2%gas(gas) = error
            return
        end if

        !> If WPL should not be applied: every path and every analyser alike.
        if (.not. EddyFlowProj%wpl) then
            Flux2%gas(gas) = Flux1%gas(gas)
            if (.not. lEx%var_present(gas)) Flux2%gas(gas) = error
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
            !> Both terms want chi/Va and reach it from the molar density, and
            !> the factor that recovers it differs by species: a trace gas has
            !> d = chi/Va * 1d-3 and water d = chi/Va. A bare 1d3 here is the
            !> trace-gas case, and a SECOND hygrometer used to come through
            !> this arm as a trace gas - so its open-path WPL term was a
            !> thousandfold too large. Hygrometers now have their own
            !> HygrometerWater and never arrive here; the water case stays
            !> correct regardless.
            if (GasSlotIsWater(gas)) then
                dens_to_chi = 1d0
            else
                dens_to_chi = 1d3
            end if

            if (IsLi7700(lEx%gas_instr(gas)%model)) then
                !> Webb et al. (1980) scaled by the LI-7700's own
                !> spectroscopic multipliers - see the RP twin for why B and C
                !> went missing and what it cost.
                !> No multipliers means no corrected flux - see the RP twin.
                if (lEx%Mul7700(gas)%A == error .or. lEx%Mul7700(gas)%B == error &
                    .or. lEx%Mul7700(gas)%C == error) then
                    Flux2%gas(gas) = error
                    return
                end if
                wpl = Flux1%gas(gas)
                if (enw_ref /= error .and. lEx%RHO%d > 0d0) &
                    wpl = wpl + lEx%Mul7700(gas)%B * mu * lEx%d(gas) * dens_to_chi &
                        * enw_ref / lEx%RHO%d
                if (Flux3%H /= error .and. lEx%RhoCp > 0d0 .and. lEx%Ta > 0d0 &
                    .and. sigma_gas /= error) &
                    wpl = wpl + lEx%Mul7700(gas)%C * (1d0 + mu * sigma_gas) * Flux3%H &
                        * lEx%d(gas) * dens_to_chi / (lEx%RhoCp * lEx%Ta)
                Flux2%gas(gas) = lEx%Mul7700(gas)%A * wpl
            else
                !> Open path, after e.g. Burba et al. (2008, GCB, eq. 1)
                wpl = Flux1%gas(gas)
                if (e3_ref /= error .and. lEx%RHO%d > 0d0 .and. sigma_gas /= error) &
                    wpl = wpl + mu * e3_ref * lEx%d(gas) * dens_to_chi &
                        / ((1d0 + mu * sigma_gas) * lEx%RHO%d)
                if (Flux3%H /= error .and. lEx%RhoCp > 0d0 .and. lEx%Ta > 0d0) &
                    wpl = wpl + (Flux3%H + BurbaHeatFor(gas)) &
                        * lEx%d(gas) * dens_to_chi / (lEx%RhoCp * lEx%Ta)
                Flux2%gas(gas) = wpl
            end if
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
