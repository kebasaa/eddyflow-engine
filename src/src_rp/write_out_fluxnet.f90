!***************************************************************************
! write_outfiles.f90
! ------------------
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
! \brief       Write all results on (temporary) output files
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine WriteOutFluxnet(StDiff, DtDiff, STFlg, DTFlg)
    use m_rp_global_var
    implicit none
    !> in/out variables
    type(QCType), intent(in) :: StDiff
    type(QCType), intent(in) :: DtDiff
    integer, intent(in) :: STFlg(GHGNumVar)
    integer, intent(in) :: DTFlg(GHGNumVar)
    !> local variables
    integer :: var
    integer :: gas
    integer :: gas1
    integer :: gas2
    integer :: j
    integer :: i
    integer :: indx
    integer :: int_doy
    real(kind = dbl) :: float_doy
    real(kind = dbl), allocatable :: bAggrOut(:)
    real(kind = dbl) :: lrad
    character(16000) :: csv_row
    character(32) :: char_doy
    character(14) :: tsIso
    character(9) :: vm97flags(GHGNumVar)
    !> FluxnetGasScale / FluxnetGasAdvScale arrive with interfaces.inc, whose
    !> first line pulls in interfaces_1.inc.
    include '../src_common/interfaces.inc'


    !> write FLUXNET output file (csv) 
    call clearstr(csv_row)

    !> Start/end imestamps
    tsIso = Stats%start_date(1:4) // Stats%start_date(6:7) // Stats%start_date(9:10) &
                // Stats%start_time(1:2) // Stats%start_time(4:5)
    call AddDatum(csv_row, trim(adjustl(tsIso)), separator)
    tsIso = Stats%date(1:4) // Stats%date(6:7) // Stats%date(9:10) &
                // Stats%time(1:2) // Stats%time(4:5)
    call AddDatum(csv_row, trim(adjustl(tsIso)), separator)

    !> DOYs
    !>  Start
    call DateTimeToDOY(Stats%start_date, Stats%start_time, int_doy, float_doy)
    write(char_doy, *) float_doy
    call AddDatum(csv_row, trim(adjustl(char_doy(1: index(char_doy, '.')+ 4))), separator)
    !>  End
    call DateTimeToDOY(Stats%date, Stats%time, int_doy, float_doy)
    write(char_doy, *) float_doy
    call AddDatum(csv_row, trim(adjustl(char_doy(1: index(char_doy, '.')+ 4))), separator)

    !> Filename
    call AddDatum(csv_row, trim(adjustl(Essentials%fname)), separator)

    !> Potential Radiations
    indx = DateTimeToHalfHourNumber(Stats%date, Stats%time) - 1
    indx = max(indx, 2)
    lrad = (PotRad(indx) + PotRad(indx - 1)) / 2
    call AddFloatDatumToDataline(lrad, csv_row, EddyFlowProj%err_label)

    !> Daytime
    if (Stats%daytime) then
        call AddDatum(csv_row, '0', separator)
    else
        call AddDatum(csv_row, '1', separator)
    endif

    !> Number of records
    !> Number of records teoretically available for current Averaging Interval
    call AddIntDatumToDataline(MaxPeriodNumRecords, csv_row, EddyFlowProj%err_label)
    !> Number of records actually available for current Averaging Interval given length of actual files
    call AddIntDatumToDataline(Essentials%n_in, csv_row, EddyFlowProj%err_label)
    !> Number of records actually available after custom flags filtering
    call AddIntDatumToDataline(Essentials%n_after_custom_flags, csv_row, EddyFlowProj%err_label)
    !> Number of records actually available after wind direction filtering
    call AddIntDatumToDataline(Essentials%n_after_wdf, csv_row, EddyFlowProj%err_label)
    !> Number of valid records for anemometric data
    call AddIntDatumToDataline(Essentials%n(w), csv_row, EddyFlowProj%err_label)
    !> Number of valid records for IRGA data  (N_in – M_diag_IRGA)
    do var = ts, gas4
        call AddIntDatumToDataline(Essentials%n(var), csv_row, EddyFlowProj%err_label)
        end do
    !> Number of valid records available for each main covariance (w/u, w/ts, w/co2, w/h2o, w/ch4, w/gas4)
    call AddIntDatumToDataline(Essentials%n_wcov(u), csv_row, EddyFlowProj%err_label)
    do var = ts, gas4
        call AddIntDatumToDataline(Essentials%n_wcov(var), csv_row, EddyFlowProj%err_label)
        end do

    !> Fluxes
    !> Fluxes level 3 (final fluxes) 
    call AddFloatDatumToDataline(Flux3%tau, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Flux3%H, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Flux3%LE, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Flux3%ET, csv_row, EddyFlowProj%err_label)
    !> Scaled to the FLUXNET basis of each gas's species, not of its slot.
    !> FluxnetGasScale returns 1 for CO2 and H2O, so passing it for every gas
    !> reproduces what the four fixed calls emitted.
    do gas = co2, gas4
        call AddFloatDatumToDataline(Flux3%gas(gas), csv_row, EddyFlowProj%err_label, &
            gain=FluxnetGasScale(gas), offset=0d0)
    end do

    !> Flux random uncertainties
    call AddFloatDatumToDataline(Essentials%rand_uncer(u), csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Essentials%rand_uncer(ts), csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Essentials%rand_uncer_LE, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Essentials%rand_uncer_ET, csv_row, EddyFlowProj%err_label)
    do gas = co2, gas4
        call AddFloatDatumToDataline(Essentials%rand_uncer(gas), csv_row, &
            EddyFlowProj%err_label, gain=FluxnetGasScale(gas), offset=0d0)
    end do

    !> Additional flux terms (single-point calculation)
    !> Storage fluxes
    call AddFloatDatumToDataline(Stor%H, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Stor%LE, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Stor%ET, csv_row, EddyFlowProj%err_label)
    do gas = co2, gas4
        call AddFloatDatumToDataline(Stor%of(gas), csv_row, EddyFlowProj%err_label, &
            gain=FluxnetGasScale(gas), offset=0d0)
    end do
    !> Advection fluxes.
    !>
    !> The three-way slot switch this replaces had no final else: a slot past
    !> gas4 matched none of its arms and emitted no datum at all, which would
    !> shift every later column once the loop widened.
    do gas = co2, gas4
        if (Stats5%Mean(w) /= error .and. Stats%d(gas) >= 0d0) then
            if (Stats5%Mean(w) /= error .and. Stats%d(gas) /= error) then
                call AddFloatDatumToDataline(Stats5%Mean(w) * Stats%d(gas), &
                    csv_row, EddyFlowProj%err_label, &
                    gain=FluxnetGasAdvScale(gas), offset=0d0)
            else
                call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
            end if
        else
            call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        end if
    end do

    !> Turbulence and micromet
    !> Unrotated and rotated wind components
    call AddFloatDatumToDataline(Stats4%Mean(u), csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Stats4%Mean(v), csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Stats4%Mean(w), csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Stats5%Mean(u), csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Stats5%Mean(v), csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Stats5%Mean(w), csv_row, EddyFlowProj%err_label)
    !> wind speed, wind direction, u*, stability, bowen ratio
    call AddFloatDatumToDataline(Ambient%WS, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Ambient%MWS, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Stats4%wind_dir, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Stats4%wind_dir_stdev, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Ambient%us, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Stats%TKE, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Ambient%L, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Ambient%zL, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Ambient%bowen, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Ambient%Ts, csv_row, EddyFlowProj%err_label)

    !> Termodynamics 
    !> Temperature, pressure, RH, VPD, e, es, etc.
    call AddFloatDatumToDataline(Stats7%Mean(ts), csv_row, EddyFlowProj%err_label, gain=1d0, offset=-273.15d0)
    call AddFloatDatumToDataline(Ambient%Ta, csv_row, EddyFlowProj%err_label, gain=1d0, offset=-273.15d0)
    call AddFloatDatumToDataline(Stats%Pr, csv_row, EddyFlowProj%err_label, gain=1d-3, offset=0d0)
    call AddFloatDatumToDataline(Stats%RH, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Ambient%Va, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(RHO%a, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Ambient%RhoCp, csv_row, EddyFlowProj%err_label)
    if (RHO%a > 0) then
        call AddFloatDatumToDataline(Ambient%RhoCp / RHO%a, csv_row, EddyFlowProj%err_label)
    else
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
    end if

    !> Water
    call AddFloatDatumToDataline(RHO%w, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Ambient%e, csv_row, EddyFlowProj%err_label, gain=1d-2, offset=0d0)
    call AddFloatDatumToDataline(Ambient%es, csv_row, EddyFlowProj%err_label, gain=1d-2, offset=0d0)
    call AddFloatDatumToDataline(Ambient%Q, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Ambient%VPD, csv_row, EddyFlowProj%err_label, gain=1d-2, offset=0d0)
    call AddFloatDatumToDataline(Ambient%Td, csv_row, EddyFlowProj%err_label, gain=1d0, offset=-273.15d0)
    !> Dry air
    call AddFloatDatumToDataline(Ambient%p_d, csv_row, EddyFlowProj%err_label, gain=1d-3, offset=0d0)
    call AddFloatDatumToDataline(RHO%d, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Ambient%Vd, csv_row, EddyFlowProj%err_label)
    !> Specific heat of evaporation
    call AddFloatDatumToDataline(Ambient%lambda, csv_row, EddyFlowProj%err_label)
    !> Wet to dry air density ratio
    call AddFloatDatumToDataline(Ambient%sigma, csv_row, EddyFlowProj%err_label)
    !> Water Use Efficiency
    !>!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*

    !> Gases
    !> Concentrations, densities and "nature" of the raw data 
    !> (mixing ratio, mole fraction, molar density)
    !> Gas concentrations, densities and timelags
    do gas = co2, gas4
        if (E2Col(gas)%present) then
            select case (E2Col(gas)%measure_type)
                case('mixing_ratio')
                    call AddDatum(csv_row, '0', separator)
                case('mole_fraction')
                    call AddDatum(csv_row, '1', separator)
                case('molar_density')
                    call AddDatum(csv_row, '2', separator)
            end select
            !> Molar density stays on the internal mmol m-3 basis for every
            !> gas, as it always has in this file; only the mole-basis pair
            !> follows the species.
            call AddFloatDatumToDataline(Stats%d(gas), csv_row, EddyFlowProj%err_label)
            call AddFloatDatumToDataline(Stats%r(gas), csv_row, EddyFlowProj%err_label, &
                gain=FluxnetGasScale(gas), offset=0d0)
            call AddFloatDatumToDataline(Stats%chi(gas), csv_row, EddyFlowProj%err_label, &
                gain=FluxnetGasScale(gas), offset=0d0)
        else
            do i = 1, 4
                call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
            end do
        end if
    end do
    !> Timelags (calculated, used, min/max/nominal) for all gases
    !> Gas timelags
    do gas = co2, gas4
        if (E2Col(gas)%present) then
            call AddFloatDatumToDataline(Essentials%actual_timelag(gas), csv_row, EddyFlowProj%err_label)
            call AddFloatDatumToDataline(Essentials%used_timelag(gas), csv_row, EddyFlowProj%err_label)
            call AddFloatDatumToDataline(E2Col(gas)%def_tl, csv_row, EddyFlowProj%err_label)
            call AddFloatDatumToDataline(E2Col(gas)%min_tl, csv_row, EddyFlowProj%err_label)
            call AddFloatDatumToDataline(E2Col(gas)%max_tl, csv_row, EddyFlowProj%err_label)
        else
            call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
            call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
            call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
            call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
            call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        end if
    end do
    !> PWB lag source (0=native, 1=S3 carry-forward, 2=instrument_shared, 3=nominal/default, 4=maxcov/default)
    do gas = co2, gas4
        if (Meth%tlag == 'pwb' .and. E2Col(gas)%present) then
            select case(trim(PWBResult(gas)%fallback_source))
            case('native')
                call AddDatum(csv_row, '0', separator)
            case('S3_carryforward')
                call AddDatum(csv_row, '1', separator)
            case('instrument_shared')
                call AddDatum(csv_row, '2', separator)
            case('nominal/default')
                call AddDatum(csv_row, '3', separator)
            case('maxcov_default')
                call AddDatum(csv_row, '4', separator)
            case default
                call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
            end select
        else
            call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        end if
    end do

!> Basic stats
    !> 25-50-75%
    do var = u, gas4
        if (var == ts) then
            call AddFloatDatumToDataline(Stats6%Median(var), csv_row, &
                EddyFlowProj%err_label, gain=1d0, offset=-273.15d0)
        else
            call AddFloatDatumToDataline(Stats6%Median(var), csv_row, EddyFlowProj%err_label)
                end if
    end do
    do var = u, gas4
        if (var == ts) then
            call AddFloatDatumToDataline(Stats6%Q1(var), csv_row, &
                EddyFlowProj%err_label, gain=1d0, offset=-273.15d0)
        else
            call AddFloatDatumToDataline(Stats6%Q1(var), csv_row, EddyFlowProj%err_label)
                end if
    end do
    do var = u, gas4
        if (var == ts) then
            call AddFloatDatumToDataline(Stats6%Q3(var), csv_row, &
                EddyFlowProj%err_label, gain=1d0, offset=-273.15d0)
        else
            call AddFloatDatumToDataline(Stats6%Q3(var), csv_row, EddyFlowProj%err_label)
                end if
    end do
    !> Standard deviation
    do var = u, gas4
        call AddFloatDatumToDataline(sqrt(Stats7%Cov(var, var)), csv_row, EddyFlowProj%err_label)
        end do
    !> Skwenesses
    do var = u, gas4
        call AddFloatDatumToDataline(Stats7%Skw(var), csv_row, EddyFlowProj%err_label)
        end do
    !> Kurtosis
    do var = u, gas4
        call AddFloatDatumToDataline(Stats7%Kur(var), csv_row, EddyFlowProj%err_label)
        end do
    !> w-covariances 
    call AddFloatDatumToDataline(Stats7%Cov(u, w), csv_row, EddyFlowProj%err_label)
    do var = ts, gas4
        call AddFloatDatumToDataline(Stats7%Cov(w, var), csv_row, EddyFlowProj%err_label)
        end do
    !> Gases covariance matrix
    do gas1 = co2, ch4
        do gas2 = gas1 + 1, gas4 
            call AddFloatDatumToDataline(Stats7%Cov(gas1, gas2), csv_row, EddyFlowProj%err_label)
        end do
    end do

!> Footprint
    call AddFloatDatumToDataline(Foot%peak, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Foot%offset, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Foot%x10, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Foot%x30, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Foot%x50, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Foot%x70, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Foot%x80, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Foot%x90, csv_row, EddyFlowProj%err_label)

!> Intermediate results
    !> Fluxes level 0 (uncorrected fluxes)
    call AddFloatDatumToDataline(Essentials%ustar, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Essentials%L, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Essentials%zL, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Flux0%tau, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Flux0%H, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Flux0%LE, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Flux0%ET, csv_row, EddyFlowProj%err_label)
    do gas = co2, gas4
        call AddFloatDatumToDataline(Flux0%gas(gas), csv_row, EddyFlowProj%err_label, &
            gain=FluxnetGasScale(gas), offset=0d0)
    end do
    !> Fluxes level 1
    call AddFloatDatumToDataline(Flux1%tau, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Flux1%H, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Flux1%LE, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Flux1%ET, csv_row, EddyFlowProj%err_label)
    do gas = co2, gas4
        call AddFloatDatumToDataline(Flux1%gas(gas), csv_row, EddyFlowProj%err_label, &
            gain=FluxnetGasScale(gas), offset=0d0)
    end do
    !> Fluxes level 2
    call AddFloatDatumToDataline(Flux2%tau, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Flux2%H, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Flux2%LE, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Flux2%ET, csv_row, EddyFlowProj%err_label)
    do gas = co2, gas4
        call AddFloatDatumToDataline(Flux2%gas(gas), csv_row, EddyFlowProj%err_label, &
            gain=FluxnetGasScale(gas), offset=0d0)
    end do

    !> Tin and Tout                 ******************************************** Add

    !> Temperature, pressure and molar volume 
    !> in the cell of closed-paths, for all gases
    !> Cell parameters              ******************************************** Mke it gas specific like molar volume
    call AddFloatDatumToDataline(Ambient%Tcell, csv_row, EddyFlowProj%err_label, gain=1d0, offset=-273.15d0)
    call AddFloatDatumToDataline(Ambient%Pcell, csv_row, EddyFlowProj%err_label, gain=1d-3, offset=0d0)

    !> Molar volume
    do gas = co2, gas4
        call AddFloatDatumToDataline(E2Col(gas)%Va, csv_row, EddyFlowProj%err_label)
        end do
    !> Evapotranspiration and sensible heat fluxes in the cell of 
    !> closed-paths (for WPL), with timelags of other gases
    call AddFloatDatumToDataline(Flux3%E_gas(co2), csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Flux3%E_gas(ch4), csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Flux3%E_gas(gas4), csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Flux3%Hi_gas(co2), csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Flux3%Hi_gas(h2o), csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Flux3%Hi_gas(ch4), csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Flux3%Hi_gas(gas4), csv_row, EddyFlowProj%err_label)
    !> Burba Terms 
    if (RPsetup%bu_corr /= 'none') then 
        call AddFloatDatumToDataline(Burba%h_bot, csv_row, EddyFlowProj%err_label)
            call AddFloatDatumToDataline(Burba%h_top, csv_row, EddyFlowProj%err_label)
            call AddFloatDatumToDataline(Burba%h_spar, csv_row, EddyFlowProj%err_label)
        else
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
    end if
    !> LI-7700 multipliers
    if (E2Col(ch4)%Instr%model(1:len_trim(E2Col(ch4)%Instr%model) - 2) &
        == 'li7700') then
        call AddFloatDatumToDataline(Mul7700%A, csv_row, EddyFlowProj%err_label)
            call AddFloatDatumToDataline(Mul7700%B, csv_row, EddyFlowProj%err_label)
            call AddFloatDatumToDataline(Mul7700%C, csv_row, EddyFlowProj%err_label)
        else
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
    end if
    !> WPL Terms                    ********************************************(Individual: H, LE, Pressure)
    !>!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*!*
    !> Spectral correction factors
    if (E2Col(u)%present) then
        call AddFloatDatumToDataline(BPCF%of(w_u), csv_row, EddyFlowProj%err_label)
    else
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
    end if
    if (E2Col(ts)%present) then
        call AddFloatDatumToDataline(BPCF%of(w_ts), csv_row, EddyFlowProj%err_label)
    else
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
    end if
    if (E2Col(h2o)%present) then
        call AddFloatDatumToDataline(BPCF%of(w_h2o), csv_row, EddyFlowProj%err_label)
    else
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
    end if
    if (E2Col(h2o)%present) then
        call AddFloatDatumToDataline(BPCF%of(w_h2o), csv_row, EddyFlowProj%err_label)
    else
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
    end if
    if (E2Col(co2)%present) then
        call AddFloatDatumToDataline(BPCF%of(w_co2), csv_row, EddyFlowProj%err_label)
    else
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
    end if
    if (E2Col(h2o)%present) then
        call AddFloatDatumToDataline(BPCF%of(w_h2o), csv_row, EddyFlowProj%err_label)
    else
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
    end if
    if (E2Col(ch4)%present) then
        call AddFloatDatumToDataline(BPCF%of(w_ch4), csv_row, EddyFlowProj%err_label)
    else
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
    end if
    if (E2Col(gas4)%present) then
        call AddFloatDatumToDataline(BPCF%of(w_gas4), csv_row, EddyFlowProj%err_label)
    else
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
    end if

    !> Increasingly filtered w/T covariances (for spectral assessment)
    call AddFloatDatumToDataline(Essentials%degH(NumDegH + 1), csv_row, EddyFlowProj%err_label)
    do j = 1, NumDegH
        call AddFloatDatumToDataline(Essentials%degH(j), csv_row, EddyFlowProj%err_label)
    end do
    !>> Number of spikes per variable
    do var = u, gas4
        call AddIntDatumToDataline(Essentials%e2spikes(var), csv_row, EddyFlowProj%err_label)
    end do
!> QC
    !>> Number or records eliminated based on custom flags
    call AddIntDatumToDataline(Essentials%m_custom_flags, csv_row, EddyFlowProj%err_label)
    !>> Number or records eliminated based on wind direction filter
    call AddIntDatumToDataline(Essentials%m_wdf, csv_row, EddyFlowProj%err_label)
    !> Summary of data values eliminated based on diagnostics
    !>> Number or records whose anemometric data was eliminated based on Anemometer diagnostics
    call AddIntDatumToDataline(Essentials%m_diag_anem, csv_row, EddyFlowProj%err_label)
    !>> Number or records whose IRGA data was eliminated based on IRGA diagnostics
    do j = 1, nFluxnetLayoutSlots
        call AddIntDatumToDataline(Essentials%m_diag_irga(FluxnetLayoutSlots(j)), &
                                   csv_row, EddyFlowProj%err_label)
    end do
    !>> Number of values eliminated by the Spike test. The wind components and
    !>> sonic temperature first, then one per configured gas.
    do var = u, ts
        call AddIntDatumToDataline(Essentials%m_despiking(var), csv_row, EddyFlowProj%err_label)
    end do
    do j = 1, nFluxnetLayoutSlots
        call AddIntDatumToDataline(Essentials%m_despiking(FluxnetLayoutSlots(j)), &
                                   csv_row, EddyFlowProj%err_label)
    end do
    !>> Number of values eliminated by the Absolute Limits test
    do var = u, ts
        call AddIntDatumToDataline(Essentials%al_s(var), csv_row, EddyFlowProj%err_label)
    end do
    do j = 1, nFluxnetLayoutSlots
        call AddIntDatumToDataline(Essentials%al_s(FluxnetLayoutSlots(j)), &
                                   csv_row, EddyFlowProj%err_label)
    end do

    !> Uncomment to reintroduce VM details
    ! !> VM97 Stats used to calculate flags
    ! !>> Amplitude resolution
    ! do j = u, gas4
    !     call AddFloatDatumToDataline(Essentials%ar_s(j), csv_row, EddyFlowProj%err_label)
    !     end do
    ! !>> Dropouts central
    ! do j = u, gas4
    !     call AddFloatDatumToDataline(Essentials%do_s_ctr(j), csv_row, EddyFlowProj%err_label)
    !     end do
    ! !>> Dropouts extremes
    ! do j = u, gas4
    !     call AddFloatDatumToDataline(Essentials%do_s_ext(j), csv_row, EddyFlowProj%err_label)
    !     end do
    ! !>> Higher moments Skewness
    ! do j = u, gas4
    !     call AddFloatDatumToDataline(Essentials%sk_s_skw(j), csv_row, EddyFlowProj%err_label)
    !     end do                          
    ! !>> Higher moments Kurtosis     
    ! do j = u, gas4
    !     call AddFloatDatumToDataline(Essentials%sk_s_kur(j), csv_row, EddyFlowProj%err_label)
    !     end do
    ! !>> AoA
    ! call AddFloatDatumToDataline(Essentials%aa_s, csv_row, EddyFlowProj%err_label)
    ! !>> Non-steady wind
    ! call AddFloatDatumToDataline(Essentials%ns_s_rnv(1), csv_row, EddyFlowProj%err_label)
    ! call AddFloatDatumToDataline(Essentials%ns_s_rnv(2), csv_row, EddyFlowProj%err_label)
    ! call AddFloatDatumToDataline(Essentials%ns_s_rns, csv_row, EddyFlowProj%err_label)


    !> VM97 flags, here organized per variable instead of per test.
    !>
    !> The wind components and sonic temperature first, then one field per
    !> configured gas, matching the header. Each field stays 9 characters - a
    !> filler digit plus the 8 test outcomes - because that width counts tests,
    !> not gases. The CharHF/CharSF strings are indexed by variable and are
    !> already FlagStrLen wide, so a gas slot past the fourth reads out of them
    !> the same way as any other.
    !> One contiguous range: firstGas is ts + 1, and SelectFluxnetGasSlots
    !> assigns the layout its slots as firstGas + j - 1, so "u,v,w,ts then the
    !> configured gases" is exactly u .. ts + nFluxnetLayoutSlots. This is the
    !> same span the retired nExVar named as gas4 - u + 1, now following the
    !> project instead of the fixed four.
    do var = u, ts + nFluxnetLayoutSlots
        vm97flags(var)(1 : 1) = '8'
        vm97flags(var)(2 : 2) = CharHF%sr(var + 1 : var + 1)
        vm97flags(var)(3 : 3) = CharHF%ar(var + 1 : var + 1)
        vm97flags(var)(4 : 4) = CharHF%do(var + 1 : var + 1)
        vm97flags(var)(5 : 5) = CharHF%al(var + 1 : var + 1)
        vm97flags(var)(6 : 6) = CharHF%sk(var + 1 : var + 1)
        vm97flags(var)(7 : 7) = CharSF%sk(var + 1 : var + 1)
        vm97flags(var)(8 : 8) = CharHF%ds(var + 1 : var + 1)
        vm97flags(var)(9 : 9) = CharSF%ds(var + 1 : var + 1)
        call AddCharDatumToDataline(vm97flags(var), csv_row, EddyFlowProj%err_label)
    end do

    !> Uncomment to reintroduce flags for last 3 tests
    call AddDatum(csv_row, '8'//CharHF%tl(FlagStrLen-3:FlagStrLen), separator)
    call AddDatum(csv_row, '8'//CharSF%tl(FlagStrLen-3:FlagStrLen), separator)
    call AddDatum(csv_row, '8'//CharHF%aa(FlagStrLen:FlagStrLen), separator)
    call AddDatum(csv_row, '8'//CharHF%ns(FlagStrLen:FlagStrLen), separator)

    !> Quality test results. Three variable-shaped families - the wind
    !> components and sonic temperature, then one field per configured gas,
    !> matching AddVariableFamily in the header.
    !> Longest Gap Duration (LGDs)
    do var = u, ts + nFluxnetLayoutSlots
        call AddFloatDatumToDataline(Essentials%LGD(var), csv_row, EddyFlowProj%err_label)
    end do
    !> Kurtosis Index on Differenced variables (KIDs)
    do var = u, ts + nFluxnetLayoutSlots
        call AddFloatDatumToDataline(Essentials%KID(var), csv_row, EddyFlowProj%err_label)
    end do
    !> Zero-Counts on Differences variables
    do var = u, ts + nFluxnetLayoutSlots
        call AddIntDatumToDataline(Essentials%ZCD(var), csv_row, EddyFlowProj%err_label)
    end do
    !> Correlation differences with and without repeated values.
    !>
    !> Flux-shaped, so the first four are the momentum, heat and the two water
    !> fluxes - and LE, ET and FH2O are all the same correlation, because all
    !> three are water fluxes. Only the tail from FC onwards is per gas.
    call AddFloatDatumToDataline(Essentials%CorrDiff(u, w), csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Essentials%CorrDiff(u, ts), csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Essentials%CorrDiff(u, h2o), csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Essentials%CorrDiff(u, h2o), csv_row, EddyFlowProj%err_label)
    do j = 1, nFluxnetLayoutSlots
        call AddFloatDatumToDataline(Essentials%CorrDiff(u, FluxnetLayoutSlots(j)), &
                                     csv_row, EddyFlowProj%err_label)
    end do
    !> Mahrt 1998 Nonstationarity Ratios. Flux-shaped without the water
    !> fluxes: momentum and sensible heat, then one per configured gas.
    !> Indexed by the w_* covariance enumeration, which is numerically the
    !> variable enumeration - w_co2 and co2 are both 5 - so a gas slot indexes
    !> this array directly.
    call AddFloatDatumToDataline(Essentials%mahrt98_NR(w_u), csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Essentials%mahrt98_NR(w_ts), csv_row, EddyFlowProj%err_label)
    do j = 1, nFluxnetLayoutSlots
        call AddFloatDatumToDataline(Essentials%mahrt98_NR(FluxnetLayoutSlots(j)), &
                                     csv_row, EddyFlowProj%err_label)
    end do
    !> Foken stats used to calculate flags: the steady-state statistic per
    !> flux, then the integral turbulence characteristics on u/w/ts.
    call AddIntDatumToDataline(STDiff%w_u, csv_row, EddyFlowProj%err_label)
    call AddIntDatumToDataline(STDiff%w_ts, csv_row, EddyFlowProj%err_label)
    do j = 1, nFluxnetLayoutSlots
        call AddIntDatumToDataline(StDiff%w_gas(FluxnetLayoutSlots(j)), &
                                   csv_row, EddyFlowProj%err_label)
    end do
    call AddIntDatumToDataline(DtDiff%u, csv_row, EddyFlowProj%err_label)
    call AddIntDatumToDataline(DtDiff%w, csv_row, EddyFlowProj%err_label)
    call AddIntDatumToDataline(DtDiff%ts, csv_row, EddyFlowProj%err_label)
    !> Partial Foken flags. STFlg is indexed by the w_* covariance enumeration
    !> and is already GHGNumVar wide, so a gas slot indexes it directly.
    call AddIntDatumToDataline(STFlg(w_u), csv_row, EddyFlowProj%err_label)
    call AddIntDatumToDataline(STFlg(w_ts), csv_row, EddyFlowProj%err_label)
    do j = 1, nFluxnetLayoutSlots
        call AddIntDatumToDataline(STFlg(FluxnetLayoutSlots(j)), csv_row, EddyFlowProj%err_label)
    end do
    call AddIntDatumToDataline(DTFlg(u), csv_row, EddyFlowProj%err_label)
    call AddIntDatumToDataline(DTFlg(w), csv_row, EddyFlowProj%err_label)
    call AddIntDatumToDataline(DTFlg(ts), csv_row, EddyFlowProj%err_label)
    !> Final Foken flags. Flux-shaped: momentum, sensible heat and the two
    !> water fluxes carry the water flag, then one per configured gas.
    call AddIntDatumToDataline(QCFlag%tau, csv_row, EddyFlowProj%err_label)
    call AddIntDatumToDataline(QCFlag%H, csv_row, EddyFlowProj%err_label)
    call AddIntDatumToDataline(QCFlag%gas(h2o), csv_row, EddyFlowProj%err_label)     !< This is for LE
    call AddIntDatumToDataline(QCFlag%gas(h2o), csv_row, EddyFlowProj%err_label)     !< This is for ET
    do j = 1, nFluxnetLayoutSlots
        call AddIntDatumToDataline(QCFlag%gas(FluxnetLayoutSlots(j)), &
                                   csv_row, EddyFlowProj%err_label)
    end do

    !> LI-7x00 diagnostics breakdown
    if (Diag7200%present) then
        call AddIntDatumToDataline(Diag7200%head_detect, csv_row, EddyFlowProj%err_label)
            call AddIntDatumToDataline(Diag7200%t_out, csv_row, EddyFlowProj%err_label)
            call AddIntDatumToDataline(Diag7200%t_in, csv_row, EddyFlowProj%err_label)
            call AddIntDatumToDataline(Diag7200%aux_in, csv_row, EddyFlowProj%err_label)
            call AddIntDatumToDataline(Diag7200%delta_p, csv_row, EddyFlowProj%err_label)
            call AddIntDatumToDataline(Diag7200%chopper, csv_row, EddyFlowProj%err_label)
            call AddIntDatumToDataline(Diag7200%detector, csv_row, EddyFlowProj%err_label)
            call AddIntDatumToDataline(Diag7200%pll, csv_row, EddyFlowProj%err_label)
            call AddIntDatumToDataline(Diag7200%sync, csv_row, EddyFlowProj%err_label)
        else
        do i = 1, 9
            call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        end do
    end if
    if (Diag7500%present) then
        call AddIntDatumToDataline(Diag7500%chopper, csv_row, EddyFlowProj%err_label)
            call AddIntDatumToDataline(Diag7500%detector, csv_row, EddyFlowProj%err_label)
            call AddIntDatumToDataline(Diag7500%pll, csv_row, EddyFlowProj%err_label)
            call AddIntDatumToDataline(Diag7500%sync, csv_row, EddyFlowProj%err_label)
        else
        do i = 1, 4
            call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        end do
    end if
    if (Diag7700%present) then
        call AddIntDatumToDataline(Diag7700%not_ready, csv_row, EddyFlowProj%err_label)
            call AddIntDatumToDataline(Diag7700%no_signal, csv_row, EddyFlowProj%err_label)
            call AddIntDatumToDataline(Diag7700%re_unlocked, csv_row, EddyFlowProj%err_label)
            call AddIntDatumToDataline(Diag7700%bad_temp, csv_row, EddyFlowProj%err_label)
            call AddIntDatumToDataline(Diag7700%laser_temp_unregulated, csv_row, EddyFlowProj%err_label)
            call AddIntDatumToDataline(Diag7700%block_temp_unregulated, csv_row, EddyFlowProj%err_label)
            call AddIntDatumToDataline(Diag7700%motor_spinning, csv_row, EddyFlowProj%err_label)
            call AddIntDatumToDataline(Diag7700%pump_on, csv_row, EddyFlowProj%err_label)
            call AddIntDatumToDataline(Diag7700%top_heater_on, csv_row, EddyFlowProj%err_label)
            call AddIntDatumToDataline(Diag7700%bottom_heater_on, csv_row, EddyFlowProj%err_label)
            call AddIntDatumToDataline(Diag7700%calibrating, csv_row, EddyFlowProj%err_label)
            call AddIntDatumToDataline(Diag7700%motor_failure, csv_row, EddyFlowProj%err_label)
            call AddIntDatumToDataline(Diag7700%bad_aux_tc1, csv_row, EddyFlowProj%err_label)
            call AddIntDatumToDataline(Diag7700%bad_aux_tc2, csv_row, EddyFlowProj%err_label)
            call AddIntDatumToDataline(Diag7700%bad_aux_tc3, csv_row, EddyFlowProj%err_label)
            call AddIntDatumToDataline(Diag7700%box_connected, csv_row, EddyFlowProj%err_label)
        else
        do i = 1, 16
            call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        end do
    end if
    !> AGC/RSSI. Header may need to adapt to AGC/RSSI for 7200/7500.
    if(CompareSwVer(E2Col(co2)%instr%sw_ver, SwVerFromString('6.0.0'))) then
        call AddIntDatumToDataline(nint(Essentials%AGC72), csv_row, EddyFlowProj%err_label)
    else
        call AddIntDatumToDataline(-nint(Essentials%AGC72), csv_row, EddyFlowProj%err_label)
    end if
    !> LI-7500
    if(CompareSwVer(E2Col(co2)%instr%sw_ver, SwVerFromString('6.0.0'))) then
        call AddIntDatumToDataline(nint(Essentials%AGC75), csv_row, EddyFlowProj%err_label)
    else
        call AddIntDatumToDataline(-nint(Essentials%AGC75), csv_row, EddyFlowProj%err_label)
    end if
    !> LI-7700
    call AddIntDatumToDataline(nint(Essentials%RSSI77), csv_row, EddyFlowProj%err_label)

!> Processing settings
    !> Whether w-boost calibration was applied
    if (RPsetup%calib_wboost) then
        call AddIntDatumToDataline(1, csv_row, EddyFlowProj%err_label)
    else
        call AddIntDatumToDataline(0, csv_row, EddyFlowProj%err_label)
    end if
    !> Whether AoA calibration was applied
    select case(trim(adjustl(RPsetup%calib_aoa)))
        case('automatic')
            call AddIntDatumToDataline(-1, csv_row, EddyFlowProj%err_label)
        case('none')
            call AddIntDatumToDataline(0, csv_row, EddyFlowProj%err_label)
        case('nakai_06')
            call AddIntDatumToDataline(1, csv_row, EddyFlowProj%err_label)
        case('nakai_12')
            call AddIntDatumToDataline(2, csv_row, EddyFlowProj%err_label)
    end select
    !> Tilt compensation method
    select case(trim(adjustl(Meth%rot)))
        case('none')
            call AddIntDatumToDataline(0, csv_row, EddyFlowProj%err_label)
        case('double_rotation')
            call AddIntDatumToDataline(1, csv_row, EddyFlowProj%err_label)
        case('triple_rotation')
            call AddIntDatumToDataline(2, csv_row, EddyFlowProj%err_label)
        case('planar_fit')
            call AddIntDatumToDataline(3, csv_row, EddyFlowProj%err_label)
        case('planar_fit_no_bias')
            call AddIntDatumToDataline(4, csv_row, EddyFlowProj%err_label)
    end select
    !> Rotation angles
    call AddFloatDatumToDataline(Essentials%yaw, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Essentials%pitch, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Essentials%roll, csv_row, EddyFlowProj%err_label)
    !> Detrending method and time constant
    select case(trim(adjustl(Meth%det)))
        case('ba')
            call AddIntDatumToDataline(0, csv_row, EddyFlowProj%err_label)
        case('ld')
            call AddIntDatumToDataline(1, csv_row, EddyFlowProj%err_label)
        case('rm')
            call AddIntDatumToDataline(2, csv_row, EddyFlowProj%err_label)
        case('ew')
            call AddIntDatumToDataline(3, csv_row, EddyFlowProj%err_label)
    end select
    call AddIntDatumToDataline(RPsetup%Tconst, csv_row, EddyFlowProj%err_label)
    !> Time lag detection method
    select case(trim(adjustl(Meth%tlag)))
        case('none')
            call AddIntDatumToDataline(0, csv_row, EddyFlowProj%err_label)
        case('constant')
            call AddIntDatumToDataline(1, csv_row, EddyFlowProj%err_label)
        case('maxcov&default')
            call AddIntDatumToDataline(2, csv_row, EddyFlowProj%err_label)
        case('maxcov')
            call AddIntDatumToDataline(3, csv_row, EddyFlowProj%err_label)
        case('tlag_opt')
            call AddIntDatumToDataline(4, csv_row, EddyFlowProj%err_label)
        case('pwb')
            call AddIntDatumToDataline(5, csv_row, EddyFlowProj%err_label)
    end select
    !> WPL terms
    if (EddyFlowProj%wpl) then
        call AddIntDatumToDataline(1, csv_row, EddyFlowProj%err_label)
    else
        call AddIntDatumToDataline(0, csv_row, EddyFlowProj%err_label)
    end if
    !> Burba terms
    if (trim(adjustl(RPSetup%bu_corr)) == 'yes') then
        if (RPSetup%bu_multi) then
            call AddIntDatumToDataline(2, csv_row, EddyFlowProj%err_label)
        else
            call AddIntDatumToDataline(1, csv_row, EddyFlowProj%err_label)
        end if
    else
        call AddIntDatumToDataline(0, csv_row, EddyFlowProj%err_label)
    end if
    !> Spectral correction method
    select case(trim(adjustl(EddyFlowProj%hf_meth)))
        case('none')
            call AddIntDatumToDataline(0, csv_row, EddyFlowProj%err_label)
        case('moncrieff_97')
            call AddIntDatumToDataline(1, csv_row, EddyFlowProj%err_label)
        case('horst_97')
            call AddIntDatumToDataline(2, csv_row, EddyFlowProj%err_label)
        case('ibrom_07')
            call AddIntDatumToDataline(3, csv_row, EddyFlowProj%err_label)
        case('fratini_12')
            call AddIntDatumToDataline(4, csv_row, EddyFlowProj%err_label)
        case('massman_00')
            call AddIntDatumToDataline(5, csv_row, EddyFlowProj%err_label)
    end select
    !> Footprint model
    select case(trim(adjustl(foot_model_used)))
        case('none')
            call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        case('kljun_04')
            call AddIntDatumToDataline(0, csv_row, EddyFlowProj%err_label)
        case('kormann_meixner_01')
            call AddIntDatumToDataline(1, csv_row, EddyFlowProj%err_label)
        case('hsieh_00')
            call AddIntDatumToDataline(2, csv_row, EddyFlowProj%err_label)
    end select

!> Metadata
    !> Data logger software version
    call AddIntDatumToDataline(Metadata%logger_swver%major, csv_row, EddyFlowProj%err_label)
    call AddIntDatumToDataline(Metadata%logger_swver%minor, csv_row, EddyFlowProj%err_label)
    call AddIntDatumToDataline(Metadata%logger_swver%revision, csv_row, EddyFlowProj%err_label)
    !> Site location and features
    call AddFloatDatumToDataline(Metadata%lat, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Metadata%lon, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Metadata%alt, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Metadata%canopy_height, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Metadata%d, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Metadata%z0, csv_row, EddyFlowProj%err_label)
    !> Data acquisition settings
    call AddIntDatumToDataline(nint(Metadata%file_length), csv_row, EddyFlowProj%err_label)
    call AddIntDatumToDataline(nint(Metadata%ac_freq), csv_row, EddyFlowProj%err_label)
    !> Flux averaging interval
    call AddIntDatumToDataline(RPsetup%avrg_len, csv_row, EddyFlowProj%err_label)
    !> master anemometer
    call AddCharDatumToDataline(E2Col(u)%instr%firm, csv_row, EddyFlowProj%err_label)
    call AddCharDatumToDataline(E2Col(u)%Instr%model, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(E2Col(u)%Instr%height, csv_row, EddyFlowProj%err_label)
    call AddCharDatumToDataline(E2Col(u)%Instr%wformat, csv_row, EddyFlowProj%err_label)
    call AddCharDatumToDataline(E2Col(u)%Instr%wref, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(E2Col(u)%Instr%north_offset, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(E2Col(u)%Instr%hpath_length, csv_row, EddyFlowProj%err_label, gain=1d2, offset=0d0)
    call AddFloatDatumToDataline(E2Col(u)%Instr%vpath_length, csv_row, EddyFlowProj%err_label, gain=1d2, offset=0d0)
    call AddFloatDatumToDataline(E2Col(u)%Instr%tau, csv_row, EddyFlowProj%err_label)
    !> Gas analyser details, one block per configured gas.
    !>
    !> Iterates the layout list rather than co2:gas4, so a fifth gas gets its
    !> analyser columns here like any other. That list carries every gas the
    !> project names - including one selected without a column - which is what
    !> keeps this block exactly as wide as the header describing it.
    do j = 1, nFluxnetLayoutSlots
        gas = FluxnetLayoutSlots(j)
        call AddCharDatumToDataline(E2Col(gas)%Instr%firm, csv_row, EddyFlowProj%err_label)
        call AddCharDatumToDataline(E2Col(gas)%Instr%model, csv_row, EddyFlowProj%err_label)
        call AddFloatDatumToDataline(E2Col(gas)%Instr%nsep, csv_row, EddyFlowProj%err_label, gain=1d2, offset=0d0)
        call AddFloatDatumToDataline(E2Col(gas)%Instr%esep, csv_row, EddyFlowProj%err_label, gain=1d2, offset=0d0)
        call AddFloatDatumToDataline(E2Col(gas)%Instr%vsep, csv_row, EddyFlowProj%err_label, gain=1d2, offset=0d0)
        call AddFloatDatumToDataline(E2Col(gas)%Instr%tube_l, csv_row, EddyFlowProj%err_label, gain=1d2, offset=0d0)
        call AddFloatDatumToDataline(E2Col(gas)%Instr%tube_d, csv_row, EddyFlowProj%err_label, gain=1d3, offset=0d0)
        call AddFloatDatumToDataline(E2Col(gas)%Instr%tube_f, csv_row, EddyFlowProj%err_label, gain=6d4, offset=0d0)
        if (gas == h2o) then
            call AddFloatDatumToDataline(E2Col(gas)%Instr%kw, csv_row, EddyFlowProj%err_label)
            call AddFloatDatumToDataline(E2Col(gas)%Instr%ko, csv_row, EddyFlowProj%err_label)
        end if
        call AddFloatDatumToDataline(E2Col(gas)%Instr%hpath_length, csv_row, EddyFlowProj%err_label, gain=1d2, offset=0d0)
        call AddFloatDatumToDataline(E2Col(gas)%Instr%vpath_length, csv_row, EddyFlowProj%err_label, gain=1d2, offset=0d0)
        call AddFloatDatumToDataline(E2Col(gas)%Instr%tau, csv_row, EddyFlowProj%err_label)
    end do

    !> Number and mean values of custom variables
    call AddIntDatumToDataline(NumUserVar, csv_row, EddyFlowProj%err_label)
    if (NumUserVar > 0) then
        do var = 1, NumUserVar
            call AddFloatDatumToDataline(UserStats%Mean(var), csv_row, EddyFlowProj%err_label)
        end do
    end if

    !> Per-gas water vapour terms, one pair per configured gas, in the same
    !> order in which InitFluxnetFile_rp writes the matching headers. Keep this
    !> block, that routine and ReadExRecord in step: the file is positional.
    call AddIntDatumToDataline(nFluxnetGasSlots, csv_row, EddyFlowProj%err_label)
    do j = 1, nFluxnetGasSlots
        var = FluxnetGasSlots(j)
        !> The slot is written explicitly: configured gases are not
        !> necessarily contiguous, so position alone would not identify them.
        call AddIntDatumToDataline(var, csv_row, EddyFlowProj%err_label)
        !> The terms of the H2O that corrects this gas, not of the gas itself -
        !> RHO%w_at is only populated for H2O slots, so writing it directly
        !> would send nothing but error for every non-H2O gas.
        indx = E2Col(var)%moist_ref
        if (indx >= firstGas .and. indx <= lastGas) then
            call AddFloatDatumToDataline(RHO%w_at(indx), csv_row, EddyFlowProj%err_label)
            call AddFloatDatumToDataline(Ambient%sigma_at(indx), csv_row, EddyFlowProj%err_label)
        else
            call AddFloatDatumToDataline(error, csv_row, EddyFlowProj%err_label)
            call AddFloatDatumToDataline(error, csv_row, EddyFlowProj%err_label)
        end if
    end do

    !> Analyser of each gas past the four historical slots, in SI units and in
    !> the same order as InitFluxnetFile_rp writes the matching headers.
    call AddIntDatumToDataline(nFluxnetInstrSlots, csv_row, EddyFlowProj%err_label)
    do j = 1, nFluxnetInstrSlots
        var = FluxnetInstrSlots(j)
        call AddIntDatumToDataline(var, csv_row, EddyFlowProj%err_label)
        call AddCharDatumToDataline(E2Col(var)%Instr%firm, csv_row, EddyFlowProj%err_label)
        call AddCharDatumToDataline(E2Col(var)%Instr%model, csv_row, EddyFlowProj%err_label)
        call AddFloatDatumToDataline(E2Col(var)%Instr%nsep, csv_row, EddyFlowProj%err_label)
        call AddFloatDatumToDataline(E2Col(var)%Instr%esep, csv_row, EddyFlowProj%err_label)
        call AddFloatDatumToDataline(E2Col(var)%Instr%vsep, csv_row, EddyFlowProj%err_label)
        call AddFloatDatumToDataline(E2Col(var)%Instr%tube_l, csv_row, EddyFlowProj%err_label)
        call AddFloatDatumToDataline(E2Col(var)%Instr%tube_d, csv_row, EddyFlowProj%err_label)
        call AddFloatDatumToDataline(E2Col(var)%Instr%tube_f, csv_row, EddyFlowProj%err_label)
        call AddFloatDatumToDataline(E2Col(var)%Instr%hpath_length, csv_row, EddyFlowProj%err_label)
        call AddFloatDatumToDataline(E2Col(var)%Instr%vpath_length, csv_row, EddyFlowProj%err_label)
        call AddFloatDatumToDataline(E2Col(var)%Instr%tau, csv_row, EddyFlowProj%err_label)
        call AddFloatDatumToDataline(E2Col(var)%Instr%kw, csv_row, EddyFlowProj%err_label)
        call AddFloatDatumToDataline(E2Col(var)%Instr%ko, csv_row, EddyFlowProj%err_label)
    end do

    !> Per-gas families for the gases past the four historical slots, in the
    !> same order as InitFluxnetFile_rp writes the matching headers.
    call AddIntDatumToDataline(nFluxnetInstrSlots, csv_row, EddyFlowProj%err_label)
    do j = 1, nFluxnetInstrSlots
        var = FluxnetInstrSlots(j)
        call AddIntDatumToDataline(var, csv_row, EddyFlowProj%err_label)
        call AddIntDatumToDataline(Essentials%n(var), csv_row, EddyFlowProj%err_label)
        call AddIntDatumToDataline(Essentials%n_wcov(var), csv_row, EddyFlowProj%err_label)
        select case (E2Col(var)%measure_type)
            case('mixing_ratio');  call AddDatum(csv_row, '0', separator)
            case('mole_fraction'); call AddDatum(csv_row, '1', separator)
            case('molar_density'); call AddDatum(csv_row, '2', separator)
            case default
                call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        end select
        call AddFloatDatumToDataline(Stats%d(var), csv_row, EddyFlowProj%err_label)
        call AddFloatDatumToDataline(Stats%r(var), csv_row, EddyFlowProj%err_label)
        call AddFloatDatumToDataline(Stats%chi(var), csv_row, EddyFlowProj%err_label)
        !> Level 0 too: with fcc_follows, RP leaves the corrected fluxes at
        !> error and FCC recomputes them, so the uncorrected flux is the only
        !> one that can cross the file.
        call AddFloatDatumToDataline(Flux0%gas(var), csv_row, EddyFlowProj%err_label)
        call AddFloatDatumToDataline(Flux3%gas(var), csv_row, EddyFlowProj%err_label)
        call AddFloatDatumToDataline(Flux1%gas(var), csv_row, EddyFlowProj%err_label)
        call AddFloatDatumToDataline(Flux2%gas(var), csv_row, EddyFlowProj%err_label)
        call AddFloatDatumToDataline(BPCF%of(var), csv_row, EddyFlowProj%err_label)
        call AddFloatDatumToDataline(Essentials%rand_uncer(var), csv_row, EddyFlowProj%err_label)
        call AddFloatDatumToDataline(Stor%of(var), csv_row, EddyFlowProj%err_label)
        call AddFloatDatumToDataline(Essentials%actual_timelag(var), csv_row, EddyFlowProj%err_label)
        call AddFloatDatumToDataline(Essentials%used_timelag(var), csv_row, EddyFlowProj%err_label)
        call AddFloatDatumToDataline(E2Col(var)%def_tl, csv_row, EddyFlowProj%err_label)
        call AddFloatDatumToDataline(E2Col(var)%min_tl, csv_row, EddyFlowProj%err_label)
        call AddFloatDatumToDataline(E2Col(var)%max_tl, csv_row, EddyFlowProj%err_label)
        if (Meth%tlag == 'pwb' .and. E2Col(var)%present) then
            select case(trim(PWBResult(var)%fallback_source))
            case('native');            call AddDatum(csv_row, '0', separator)
            case('S3_carryforward');   call AddDatum(csv_row, '1', separator)
            case('instrument_shared'); call AddDatum(csv_row, '2', separator)
            case('nominal/default');   call AddDatum(csv_row, '3', separator)
            case('maxcov_default');    call AddDatum(csv_row, '4', separator)
            case default
                call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
            end select
        else
            call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        end if
        call AddFloatDatumToDataline(Stats6%Median(var), csv_row, EddyFlowProj%err_label)
        call AddFloatDatumToDataline(Stats6%Q1(var), csv_row, EddyFlowProj%err_label)
        call AddFloatDatumToDataline(Stats6%Q3(var), csv_row, EddyFlowProj%err_label)
        call AddFloatDatumToDataline(sqrt(Stats7%Cov(var, var)), csv_row, EddyFlowProj%err_label)
        call AddFloatDatumToDataline(Stats7%Skw(var), csv_row, EddyFlowProj%err_label)
        call AddFloatDatumToDataline(Stats7%Kur(var), csv_row, EddyFlowProj%err_label)
        call AddFloatDatumToDataline(Stats7%Cov(w, var), csv_row, EddyFlowProj%err_label)
        call AddFloatDatumToDataline(E2Col(var)%Va, csv_row, EddyFlowProj%err_label)
        call AddIntDatumToDataline(Essentials%e2spikes(var), csv_row, EddyFlowProj%err_label)
    end do

    !> All aggregated biomet values in FLUXNET units
    call AddIntDatumToDataline(nbVars, csv_row, EddyFlowProj%err_label)
    if (nbVars > 0) then
        if (.not. allocated(bAggrOut)) allocate(bAggrOut(size(bAggr)))
        if (EddyFlowProj%fluxnet_standardize_biomet) then
            bAggrOut = bAggrFluxnet
        else
            bAggrOut = bAggr
        end if

        do i = 1, nbVars
            call AddFloatDatumToDataline(bAggrOut(i), csv_row, EddyFlowProj%err_label)
        end do

        if (allocated(bAggrOut)) deallocate(bAggrOut)
    end if

    !> CEC partitioning ratios (always written; error when do_cec=0)
    call AddFloatDatumToDataline(CECDescriptor%r_ET, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(CECDescriptor%r_Fc, csv_row, EddyFlowProj%err_label)
    call AddIntDatumToDataline(CECDescriptor%n_valid, csv_row, EddyFlowProj%err_label)
    call AddIntDatumToDataline(CECDescriptor%n_O1, csv_row, EddyFlowProj%err_label)
    call AddIntDatumToDataline(CECDescriptor%n_O2, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(CECDescriptor%frac_O1, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(CECDescriptor%frac_O2, csv_row, EddyFlowProj%err_label)
    call AddIntDatumToDataline(merge(1, 0, CECDescriptor%h2o_valid), csv_row, EddyFlowProj%err_label)
    call AddIntDatumToDataline(merge(1, 0, CECDescriptor%co2_valid), csv_row, EddyFlowProj%err_label)
    call AddIntDatumToDataline(CECDescriptor%h2o_status, csv_row, EddyFlowProj%err_label)
    call AddIntDatumToDataline(CECDescriptor%co2_status, csv_row, EddyFlowProj%err_label)

    !> Replace error codes with user-defined error code
    csv_row = replace2(csv_row, ',-9999,', ',' // trim(EddyFlowProj%err_label) // ',')
    csv_row = replace2(csv_row, ',NaN,',   ',' // trim(EddyFlowProj%err_label) // ',')
    csv_row = replace2(csv_row, ',+Inf,', ',' // trim(EddyFlowProj%err_label) // ',')
    csv_row = replace2(csv_row, ',-Inf,', ',' // trim(EddyFlowProj%err_label) // ',')
    csv_row = replace2(csv_row, ',Inf,', ',' // trim(EddyFlowProj%err_label) // ',')
    csv_row = replace2(csv_row, ',+Infinity,', ',' // trim(EddyFlowProj%err_label) // ',')
    csv_row = replace2(csv_row, ',-Infinity,', ',' // trim(EddyFlowProj%err_label) // ',')
    csv_row = replace2(csv_row, ',Infinity,', ',' // trim(EddyFlowProj%err_label) // ',')

    write(uflxnt, '(a)') csv_row(1:len_trim(csv_row) - 1)
end subroutine WriteOutFluxnet
