!***************************************************************************
! flux_params.f90
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
! \brief       Calculate micromet and auxilary params useful for \n
!              flux computation and correction, and for user analysis
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine FluxParams(printout)
    use m_rp_global_var
    implicit none
    !> in/out variables
    logical, intent(in) :: printout
    !> local variables
    ! real(kind = dbl) :: Ma
    real(kind = dbl) :: Cpd
    real(kind = dbl) :: Cpv
    !> Per-hygrometer counterparts of Cpv and the dry air partial pressure,
    !> both of which follow that hygrometer's own humidity.
    real(kind = dbl) :: Cpv_at
    real(kind = dbl) :: p_d_at
    integer :: msl
    integer :: wsl
    include '../src_common/interfaces_1.inc'

    if (printout) write(*,'(a)', advance = 'no') &
        '  Calculating auxiliary variables..'

    !> The scalar humidity, and everything derived from it, comes from the
    !> designated water record - not from the h2o slot, which is record two and
    !> holds water only by convention. Every hygrometer's own regime is
    !> computed into the per-slot arrays at the end of this routine; the
    !> scalars are the designated one's, so they are what a single-hygrometer
    !> site has always had.
    wsl = PrimaryWaterSlot()

    Ambient%alpha = 0.51d0

    !> Water vapour partial pressure at saturation [Pa]
    !> (this formula gives same results as that in Buck (1986)
    !> Buck (1996), Buck Research CR-1A User's Manual, Appendix 1
    !> Ambient%es = 611.21 * np.exp( (18.678 - T / 234.5) * (T / (257.14 + T)) )
    !> Where T is in Celsius!
    if (Stats%T > 0d0) then
        Ambient%es = (dexp(77.345d0 + 0.0057d0 * Stats%T &
                      - 7235.d0 / Stats%T)) / Stats%T**(8.2d0)
    else
        Ambient%es = error
    end if

    !> Humidity, from whichever source the site has.
    !>
    !> The order of these tests matters. It asks "is any humidity available?"
    !> before "is there a hygrometer?", because a biomet RH sensor is enough
    !> for the moist-air correction and a site may have one without an IRGA
    !> water channel. Testing for the hygrometer first discarded that: the
    !> guard that stopped the engine reading a non-water slot ran instead of
    !> the biomet branch, so such a site got no correction at all.
    if (biomet%val(bRH) > 0d0 .and. biomet%val(bRH) < RHmax) then
        !> If meteo RH is available, uses it for all slow parameters,
        !> including redefining chi, r and d of H2O
        Stats%RH = biomet%val(bRH)
        !> Water vapour partial pressure [Pa]
        if (Ambient%es /= error) then
            Ambient%e = Stats%RH * 1d-2 * Ambient%es
        else
            Ambient%e = error
        end if
        !> vapor pressure deficit [Pa]
        if (Ambient%e /= error) then
            Ambient%VPD = Ambient%es - Ambient%e
        else
            Ambient%VPD = error
        end if
        !> Water vapour mass density [kg_w m-3]
        if (Ambient%e /= error .and. Stats%T /= error .and. Stats%T /= 0d0) then
            RHO%w = Ambient%e / (Rw * Stats%T)
        else
            RHO%w = error
        end if
        !> Water vapour concentrations and densities, for EVERY hygrometer.
        !>
        !> This wrote the primary's slot alone, so on a site with two
        !> hygrometers one of them reported biomet and the other reported what
        !> it measured - and which was which followed the primary designation,
        !> a naming choice that has no business deciding whose numbers are
        !> real. On CH-LAE the tell was a mixing ratio of 19.9081 that followed
        !> the primary slot between two runs and matched neither instrument
        !> (the LI-7200 read 17.1089, the MIRO 16.354).
        !>
        !> There is one humidity in the air above the tower. If the biomet
        !> sensor is the better source - over a long deployment it usually is,
        !> which is why this override exists at all - it is the better source
        !> for every hygrometer.
        !>
        !> Per slot, not per site, for the cell molar volume and the path type:
        !> a closed-path hygrometer's molar density goes through its own cell,
        !> and two analysers do not share one.
        !>
        !> Everything downstream follows without further help. RHO%w_at is
        !> built from Stats%chi further down this routine, so sigma_at, Q_at,
        !> rho_a_at, RhoCp_at and RH_at all become the same regime - which is
        !> the point: a gas is WPL-corrected with the humidity of the air, not
        !> with whichever instrument happened to be listed first.
        do msl = firstGas, lastGas
            if (.not. GasSlotIsWater(msl)) cycle
            if (.not. E2Col(msl)%present) cycle
            if (RHO%w /= error .and. Ambient%Va /= error) then
                Stats%chi(msl) = RHO%w * Ambient%Va / MW_H2O * 1d3
                !> Water vapour mixing ratio
                Stats%r(msl)   = Stats%chi(msl) / (1.d0 - Stats%chi(msl) * 1d-3)
                !> Water vapour molar density
                if (E2Col(msl)%instr%path_type == 'closed') then
                    if (E2Col(msl)%Va > 0d0) then
                        Stats%d(msl) = Stats%chi(msl) / E2Col(msl)%Va
                    else
                        Stats%d(msl) = error
                        Stats%r(msl) = error
                        Stats%chi(msl) = error
                    end if
                else
                    Stats%d(msl) = Stats%chi(msl) / Ambient%Va
                end if
            else
                Stats%chi(msl) = error
                Stats%r(msl) = error
                Stats%d(msl) = error
            end if
        end do
    elseif (wsl >= firstGas) then
        !> If meteo RH is not available or out of range, uses H2O from raw data
        !> Molecular weight of wet air:
        !> Ma = chi(h2o) * MW_H2O + chi(dry_air) * Md
        !> if chi(dry_air) = 1 - chi(h2o) (assumes chi(h2o) in mmol mol_a-1)
        ! if (Stats%chi(h2o) > 0d0) then
        !     Ma = (Stats%chi(h2o) * 1d-3) * MW_H2O &
        !        + (1d0 - Stats%chi(h2o) * 1d-3) * Md
        ! else
        !     Ma = error
        ! end if

        !> Water vapour mass density [kg_w m-3]
        !> from mole fraction [mmol_w / mol_a]
        !> (good also when native is molar density)
        if (Stats%chi(wsl) > 0d0 .and. Ambient%Va > 0d0) then
            RHO%w = (Stats%chi(wsl) / Ambient%Va) * MW_H2O * 1d-3
        else
            RHO%w = error
        end if

        if (Stats%T > 0d0 .and. RHO%w >= 0d0) then
            !> Water vapour partial pressure [Pa]
            Ambient%e  =  RHO%w * Rw * Stats%T
            if (Ambient%es > 0d0) then
                !> Relative huimidity [%]
                Stats%RH = Ambient%e * 1d2 / Ambient%es
                !> vapor pressure deficit [hPa]
                Ambient%VPD = Ambient%es - Ambient%e
            else
                Stats%RH    = 0d0
                Ambient%VPD = 0d0
            end if
        else
            Ambient%e = error
            Ambient%es = error
            Stats%RH = error
            Ambient%VPD = error
        end if
        if (Stats%RH < 0d0 .or. Stats%RH > RHmax) then
            Stats%RH = error
            Ambient%VPD = error
        end if
        if (Stats%RH > 100d0 .and. Stats%RH < RHmax) then
            Stats%RH = 100d0 !< RH slightly higher than 100% is set to 100%
            Ambient%VPD = 0d0 !< RH slightly higher than 100%, VPD is set to 0
        end if
    else
        !> No humidity from any source - no hygrometer and no usable biomet
        !> RH. Every humidity-derived quantity is *not performed* rather than
        !> computed from a slot that holds something else, which is what
        !> reading the h2o slot did. WarnIfNoHumidity below says what that
        !> costs.
        RHO%w = error
        Ambient%e = error
        Ambient%VPD = error
        Stats%RH = error
    end if

    !> Water vapour mass density for every H2O measurement the project
    !> describes, so a gas can be corrected with the humidity from its own
    !> analyser. With a single H2O this reduces to RHO%w above, which is why
    !> the existing configuration comes out unchanged.
    !>
    !> Outside the branches, because it belongs to all of them. It used to sit
    !> inside the raw-data arm alone, so a site with biomet RH kept whatever
    !> the *previous* averaging period left here - RHO is a module global with
    !> no per-period reset - and every gas was then WPL-corrected with last
    !> half-hour's humidity. Where humidity comes from biomet and there is no
    !> hygrometer at all there is nothing to fill, and these stay `error`.
    RHO%w_at = error
    do msl = firstGas, lastGas
        if (.not. E2Col(msl)%present) cycle
        if (.not. GasSlotIsWater(msl)) cycle
        if (Stats%chi(msl) > 0d0 .and. Ambient%Va > 0d0) &
            RHO%w_at(msl) = (Stats%chi(msl) / Ambient%Va) * MW_H2O * 1d-3
    end do

    !> Dew-point temperature [K], after Campbell and Norman (1998)
    !> Environmental Biophysics. Here e is in Pa, thus it must be divided
    !> by 10^3 to get kPa as in the formula.
    if (Ambient%e > 0d0) then
        Ambient%Td = (240.97d0 * dlog(Ambient%e * 1d-3 /0.611d0) &
            / (17.502d0 - dlog(Ambient%e * 1d-3 / 0.611d0))) + 273.15d0
    else
        Ambient%Td = error
    end if

    !> Dry air partial pressure [Pa], as:
    !> ambient P minus water vapor partial pressure
    if (Stats%Pr > 0d0) then
        if (Ambient%e > 0d0) then
            Ambient%p_d =  Stats%Pr - Ambient%e
        else
            Ambient%p_d = Stats%Pr
        end if
    else
        Ambient%p_d = error
    end if

    !> Molar volume of dry air [m+3 mol_d-1] after Ibrom et al. (2007, Tellus B)
    if (Ambient%p_d > 0d0) then
        Ambient%Vd = (Stats%Pr * Ambient%Va) / Ambient%p_d
    else
        Ambient%Vd = error
    end if

    !> Density of dry air [kg_d m-3]
    if (Stats%T > 0d0) then
        RHO%d = Ambient%p_d / (Rd * Stats%T)
    else
        RHO%d = error
    end if

    !> Dry air heat capacity at costant pressure [J+1kg-1K-1],
    !> as a function of temperature
    Cpd = 1005d0 + (Stats%T - 273.15d0 + 23.12d0)**2 / 3364d0

    !> Density of wet air [kg_a m-3]
    if (RHO%d > 0d0) then
        if (RHO%w >= 0d0) then
            RHO%a = RHO%d + RHO%w
            !> alternative: analytically derived from RHO%a = Pa * Ma / (Ru * T)
            !> Gives identical result.
            !RHO%a = (Stats%Pr - (1d0-MW_H2O/Md) * Ambient%e) / (Rd*Stats%T)
        else
            RHO%a = RHO%d
        end if
    else
        RHO%a = error
    end if

    !> Specific humidity [kg_w kg_a-1]
    if (RHO%a > 0d0 .and. RHO%w >= 0d0) then
        Ambient%Q = RHO%w / RHO%a
    else
        Ambient%Q = error
    end if

    !> Air temperature = sonic temperature (corrected for side-wind) \n
    !> corrected for humidity (T in K), or = biomet T
    !> Condition is posed on either biomet T or raw air T (Mean(te))
    if (Stats%Mean(te) > 0d0 .or. biomet%val(bTa) > 0d0) then
        Ambient%Ta = Stats%T
        if (Stats%Mean(ts) > 0d0) then
            !> temperature mapping factor (Van Dijk et al. 2004, eq.3.1)
            Ambient%Tmap = Ambient%Ta / Stats%Mean(ts)
        else
            Ambient%Tmap = error
        end if
    elseif (E2Col(ts)%instr%category == 'fast_t_sensor') then
        !> If Ts was actually from a fast temperature sensor,
        !> do not apply Q correction
        Ambient%Ta = Stats%Mean(ts)
        Ambient%Tmap = 1d0
    else
        if (Ambient%Q > 0d0 .and. Ambient%alpha /= error &
            .and. Stats%Mean(ts) > 0d0) then
            Ambient%Ta = Stats%Mean(ts) / (1.d0 + Ambient%alpha * Ambient%Q)
            Ambient%Tmap = Ambient%Ta / Stats%Mean(ts)
        else
            Ambient%Ta = Stats%Mean(ts)
            Ambient%Tmap = 1d0
        end if

        !> Iterate the calculation of main quantities,
        !> after having better estimated air T
        if (Ambient%Ta > 0d0) then
            Ambient%es = (dexp(77.345d0 + 0.0057d0 * Ambient%Ta &
                          - 7235.d0 / Ambient%Ta)) / Ambient%Ta**(8.2d0)

            if (RHO%w >= 0d0) then
                Ambient%e  =  RHO%w * Rw * Ambient%Ta
                Stats%RH = Ambient%e * 1d2 / Ambient%es
                Ambient%VPD = Ambient%es - Ambient%e
                if (Stats%RH < 0d0 .or. Stats%RH > RHmax) then
                    Stats%RH = error
                    Ambient%VPD = error
                end if
                if (Stats%RH > 100d0 .and. Stats%RH < RHmax) then
                    Stats%RH = 100d0 !< RH slightly higher than 100% is set to 100%
                    Ambient%VPD = 0d0 !< RH slightly higher than 100% VPD is set to 0
                end if
                Ambient%Td = (240.97d0 * dlog(Ambient%e * 1d-3 /0.611d0) &
                    / (17.502d0 - dlog(Ambient%e * 1d-3 / 0.611d0))) + 273.15d0
                Ambient%p_d =  Stats%Pr - Ambient%e
                if (Ambient%p_d > 0d0) then
                    Ambient%Vd = (Stats%Pr * Ambient%Va) / Ambient%p_d
                else
                    Ambient%Vd = error
                end if
                RHO%d = Ambient%p_d / (Rd * Ambient%Ta)
                RHO%a = RHO%d + RHO%w
                Cpd = 1005d0 + (Ambient%Ta - 273.15d0 + 23.12d0)**2 / 3364d0
                Ambient%Q = RHO%w / RHO%a
                if (E2Col(ts)%instr%category == 'fast_t_sensor') then
                    !> If Ts was actually from a fast temperature sensor,
                    !> do not apply Q correction
                    Ambient%Ta = Stats%Mean(ts)
                    Ambient%Tmap = 1d0
                else
                    Ambient%Ta = Stats%Mean(ts) &
                        / (1.d0 + Ambient%alpha * Ambient%Q)
                    Ambient%Tmap = Ambient%Ta / Stats%Mean(ts)
                end if
            else
                Ambient%e = error
                Stats%RH = error
                Ambient%Q = error
                Ambient%Td = error
                Ambient%p_d = Stats%Pr
                Ambient%Vd = (Stats%Pr * Ambient%Va) / Ambient%p_d
                RHO%d = Ambient%p_d / (Rd * Ambient%Ta)
                RHO%a = RHO%d
                Cpd = 1005d0 + (Ambient%Ta - 273.15d0 + 23.12d0)**2 &
                    / 3364d0
                Ambient%Ta  = Stats%Mean(ts)
                Ambient%Tmap = 1
            end if
        end if
    end if

    !> Cell Temperature, if applicable
    if (Stats%Mean(tc) > 0d0) then
        Ambient%Tcell = Stats%Mean(tc)
    else
        Ambient%Tcell = Ambient%Ta
    end if

    !> Water vapour heat capacity at costant pressure [J+1kg-1K-1],
    !> as a function of temperature and RH
    Cpv = 1859d0 + 0.13d0 * Stats%RH &
        + (0.193d0 + 5.6d-3 * Stats%RH) * (Ambient%Ta - 273.15d0) &
        + (1d-3 + 5d-5 * Stats%RH) * (Ambient%Ta - 273.15d0)**2
    !> RhoAir by Cp (this is wet air Cp), in [J+1K-1m-3]
    if (RHO%d > 0d0 .and. RHO%w >= 0d0) then
            Ambient%RhoCp = Cpv * RHO%w + Cpd * Rho%d
            !> Alternative formulation (gives identical result)
            !Ambient%RhoCp = RHO%a * (Cpd * (1d0 - Ambient%Q) + Cpv * Ambient%Q)
        elseif (RHO%d > 0d0) then
            !> If RHO%d exists but RHO%w not, RhoCp is
            !> calculated as referred to dry air
            Ambient%RhoCp = Cpd * Rho%d
        else
            Ambient%RhoCp = error
    end if

    !> Specific heat of evaporation [J kg-1 K-1]
    if (Ambient%Ta > 0d0) then
        !> Gives same result of:
        !> lambda = − 0.0000614342*T^3 + 0.00158927*T^2 − 2.36418*T
        !> + 2500.79 in a large range (-30 to 50 °C)
        Ambient%lambda = (3147.5d0 - 2.37d0 * Ambient%Ta) * 1d3
    else
        Ambient%lambda = error
    end if

    !> water to dry air density ratio [adim.]
    if (RHO%d >0d0 .and. RHO%w > 0d0) then
        Ambient%sigma = RHO%w / RHO%d
    else
        Ambient%sigma = error
    end if

    !> The whole moisture regime, once per hygrometer.
    !>
    !> Each of these is the corresponding scalar above, recomputed from one
    !> hygrometer's own vapour density. Not merely sigma: a second hygrometer
    !> reads a different humidity, which is a different water vapour partial
    !> pressure, which leaves a different dry-air partial pressure and so its
    !> own dry- and wet-air density, specific humidity and heat capacity. Take
    !> sigma from one hygrometer and RhoCp from another and the two halves of
    !> the same correction describe different air.
    !>
    !> All on the one settled Ambient%Ta, which is why this sits at the end of
    !> the routine rather than inside the branches that establish it. Cpd is a
    !> function of Ta alone and so is shared; Cpv depends on RH and is not.
    !>
    !> This loop very nearly reproduces the scalars for the designated
    !> hygrometer, but not to the last digit - see the assignment below it,
    !> which is what actually makes them equal.
    Ambient%e_at     = error
    Ambient%RH_at    = error
    Ambient%rho_d_at = error
    Ambient%rho_a_at = error
    Ambient%Q_at     = error
    Ambient%RhoCp_at = error
    Ambient%sigma_at = error

    do msl = firstGas, lastGas
        if (RHO%w_at(msl) < 0d0) cycle
        if (Ambient%Ta <= 0d0 .or. Stats%Pr <= 0d0) cycle

        !> Water vapour partial pressure [Pa] and relative humidity [%]
        Ambient%e_at(msl) = RHO%w_at(msl) * Rw * Ambient%Ta
        if (Ambient%es > 0d0) then
            Ambient%RH_at(msl) = Ambient%e_at(msl) * 1d2 / Ambient%es
            !> Clamped as Stats%RH is: slightly over saturation is measurement
            !> noise, far over is a fault and disqualifies the value.
            if (Ambient%RH_at(msl) < 0d0 &
                .or. Ambient%RH_at(msl) > RHmax) then
                Ambient%RH_at(msl) = error
            elseif (Ambient%RH_at(msl) > 1d2) then
                Ambient%RH_at(msl) = 1d2
            end if
        end if

        !> Dry and wet air density [kg m-3]
        p_d_at = Stats%Pr - Ambient%e_at(msl)
        if (p_d_at <= 0d0) cycle
        Ambient%rho_d_at(msl) = p_d_at / (Rd * Ambient%Ta)
        if (Ambient%rho_d_at(msl) <= 0d0) cycle
        Ambient%rho_a_at(msl) = Ambient%rho_d_at(msl) + RHO%w_at(msl)

        !> Specific humidity [kg_w kg_a-1] and density ratio [adim.]
        if (Ambient%rho_a_at(msl) > 0d0) &
            Ambient%Q_at(msl) = RHO%w_at(msl) / Ambient%rho_a_at(msl)
        if (RHO%w_at(msl) > 0d0) &
            Ambient%sigma_at(msl) = RHO%w_at(msl) / Ambient%rho_d_at(msl)

        !> Wet air heat capacity [J K-1 m-3]. Cpv follows this hygrometer's RH;
        !> where that is unavailable the dry-air term stands alone, as it does
        !> in the scalar path when there is no vapour density at all.
        if (Ambient%RH_at(msl) /= error) then
            Cpv_at = 1859d0 + 0.13d0 * Ambient%RH_at(msl) &
                + (0.193d0 + 5.6d-3 * Ambient%RH_at(msl)) &
                  * (Ambient%Ta - 273.15d0) &
                + (1d-3 + 5d-5 * Ambient%RH_at(msl)) &
                  * (Ambient%Ta - 273.15d0)**2
            Ambient%RhoCp_at(msl) = Cpv_at * RHO%w_at(msl) &
                + Cpd * Ambient%rho_d_at(msl)
        else
            Ambient%RhoCp_at(msl) = Cpd * Ambient%rho_d_at(msl)
        end if
    end do

    !> The designated hygrometer's entries are the scalars themselves.
    !>
    !> The loop above is meant to reproduce them and very nearly does, but not
    !> to the last digit: it takes RHO%w_at(wsl), built from Stats%chi, where
    !> the scalar path reaches RHO%w. On CH-LAE the two RH values differ in the
    !> fourth decimal - nothing, until it decides which side of an RH-class
    !> boundary a period falls on, and then it moves a cutoff frequency and
    !> every flux behind it.
    !>
    !> Assigning removes the question. A project with one hygrometer gets
    !> exactly the numbers it got before this existed, because they *are* those
    !> numbers rather than a recomputation that agrees to four places.
    !>
    !> sigma_at is deliberately not in this list: it predates the per-
    !> hygrometer work, the WPL dilution already reads it, and its designated
    !> entry has been carrying the loop's value all along. Assigning it here
    !> would move output that has been verified against v7.2.5.
    if (wsl >= firstGas .and. wsl <= lastGas) then
        Ambient%e_at(wsl)     = Ambient%e
        Ambient%RH_at(wsl)    = Stats%RH
        Ambient%rho_d_at(wsl) = RHO%d
        Ambient%rho_a_at(wsl) = RHO%a
        Ambient%Q_at(wsl)     = Ambient%Q
        Ambient%RhoCp_at(wsl) = Ambient%RhoCp
    end if

    if (printout) write(*,'(a)') ' Done.'
end subroutine FluxParams

