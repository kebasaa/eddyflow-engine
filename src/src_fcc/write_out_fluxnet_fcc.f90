!***************************************************************************
! write_out_fluxnet_fcc.f90
! -------------------------
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
! \brief       Write results on output files
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine WriteOutFluxnetFcc(lEx)
    use m_fx_global_var
    implicit none
    !> in/out variables
    Type(ExType), intent(in) :: lEx
    character(16000) :: csv_row

    !> local variables
    integer :: var
    integer :: i
    integer :: gas
    integer :: igas
    integer :: vi
    character(9) :: vm97flags(GHGNumVar)
    include '../src_common/interfaces_1.inc'


    call clearstr(csv_row)
    !> Timestamp
    !> Start/end imestamps
    call AddDatum(csv_row, trim(adjustl(lEx%start_timestamp)), separator)
    call AddDatum(csv_row, trim(adjustl(lEx%end_timestamp)), separator)
    call AddFloatDatumToDataline(lEx%DOY_start, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%DOY_end, csv_row, EddyFlowProj%err_label)

    !> Filename
    call AddCharDatumToDataline(lEx%fname, csv_row, EddyFlowProj%err_label)

    !> Potential radiation and daytime
    call AddFloatDatumToDataline(lEx%RP, csv_row, EddyFlowProj%err_label)
    call AddIntDatumToDataline(lEx%nighttime_int, csv_row, EddyFlowProj%err_label)

    !> Number of records
    call AddIntDatumToDataline(lEx%nr_theor, csv_row, EddyFlowProj%err_label)        
    call AddIntDatumToDataline(lEx%nr_files, csv_row, EddyFlowProj%err_label)
    call AddIntDatumToDataline(lEx%nr_after_custom_flags, csv_row, EddyFlowProj%err_label)
    call AddIntDatumToDataline(lEx%nr_after_wdf, csv_row, EddyFlowProj%err_label)
    call AddIntDatumToDataline(lEx%nr(u), csv_row, EddyFlowProj%err_label)
    do var = ts, gas4
        call AddIntDatumToDataline(lEx%nr(var), csv_row, EddyFlowProj%err_label)
    end do
    call AddIntDatumToDataline(lEx%nr_w(u), csv_row, EddyFlowProj%err_label)
    do var = ts, gas4
        call AddIntDatumToDataline(lEx%nr_w(var), csv_row, EddyFlowProj%err_label)
    end do

    !> Final fluxes
    call AddFloatDatumToDataline(Flux3%Tau, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Flux3%H, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Flux3%LE, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Flux3%ET, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Flux3%co2, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Flux3%h2o, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Flux3%ch4, csv_row, EddyFlowProj%err_label, gain=1d3, offset=0d0)
    call AddFloatDatumToDataline(Flux3%gas4, csv_row, EddyFlowProj%err_label, gain=1d3, offset=0d0)

    !> Random uncertainties
    call AddFloatDatumToDataline(lEx%rand_uncer(u), csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%rand_uncer(ts), csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%rand_uncer_LE, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%rand_uncer_ET, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%rand_uncer(co2), csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%rand_uncer(h2o), csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%rand_uncer(ch4), csv_row, EddyFlowProj%err_label, gain=1d3, offset=0d0)
    call AddFloatDatumToDataline(lEx%rand_uncer(gas4), csv_row, EddyFlowProj%err_label, gain=1d3, offset=0d0)

    !> Storage fluxes
    call AddFloatDatumToDataline(lEx%Stor%H, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%Stor%LE, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%Stor%ET, csv_row, EddyFlowProj%err_label)
    do gas = co2, h2o
        call AddFloatDatumToDataline(lEx%Stor%of(gas), csv_row, EddyFlowProj%err_label)
        end do
    do gas = ch4, gas4
        call AddFloatDatumToDataline(lEx%Stor%of(gas), csv_row, EddyFlowProj%err_label)
    end do

    !> Advection fluxes
    do gas = co2, gas4
        if (lEx%rot_w /= error .and. lEx%d(gas) >= 0d0) then
            if (lEx%rot_w /= error .and. lEx%d(gas) /= error) then
                if (gas == co2) then
                    call AddFloatDatumToDataline(lEx%rot_w * lEx%d(gas), &
                        csv_row, EddyFlowProj%err_label, gain=1d3, offset=0d0)
                else if (gas == h2o) then
                    call AddFloatDatumToDataline(lEx%rot_w * lEx%d(gas), csv_row, EddyFlowProj%err_label)
                else if (gas == ch4 .or. gas == gas4) then
                    call AddFloatDatumToDataline(lEx%rot_w * lEx%d(gas), &
                        csv_row, EddyFlowProj%err_label, gain=1d6, offset=0d0)
                end if
            else
                call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
            end if
        else
            call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        end if
    end do

    !> Turbulence and micromet
    !> Unrotated and rotated wind components
    call AddFloatDatumToDataline(lEx%unrot_u, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%unrot_v, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%unrot_w, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%rot_u, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%rot_v, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%rot_w, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%WS, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%MWS, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%WD, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%WD_SIGMA, csv_row, EddyFlowProj%err_label)

    !> Turbulence
    call AddFloatDatumToDataline(Flux3%ustar, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%TKE, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%L, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%zL, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%bowen, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%Tstar, csv_row, EddyFlowProj%err_label)

    !> Thermodynamics
    !> Temperature, pressure, RH, VPD, e, es, etc.
    call AddFloatDatumToDataline(lEx%Ts, csv_row, EddyFlowProj%err_label, gain=1d0, offset=-273.15d0)
    call AddFloatDatumToDataline(lEx%Ta, csv_row, EddyFlowProj%err_label, gain=1d0, offset=-273.15d0)
    call AddFloatDatumToDataline(lEx%Pa, csv_row, EddyFlowProj%err_label, gain=1d-3, offset=0d0)
    call AddFloatDatumToDataline(lEx%RH, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%Va, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%RHO%a, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%RhoCp, csv_row, EddyFlowProj%err_label)
    if (lEx%RHO%a > 0) then
        call AddFloatDatumToDataline(lEx%RhoCp / lEx%RHO%a, csv_row, EddyFlowProj%err_label)
    else
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
    end if

    !> Water
    call AddFloatDatumToDataline(lEx%RHO%w, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%e, csv_row, EddyFlowProj%err_label, gain=1d-2, offset=0d0)
    call AddFloatDatumToDataline(lEx%es, csv_row, EddyFlowProj%err_label, gain=1d-2, offset=0d0)
    call AddFloatDatumToDataline(lEx%Q, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%VPD, csv_row, EddyFlowProj%err_label, gain=1d-2, offset=0d0)
    call AddFloatDatumToDataline(lEx%Tdew, csv_row, EddyFlowProj%err_label, gain=1d0, offset=-273.15d0)
    !> Dry air
    call AddFloatDatumToDataline(lEx%Pd, csv_row, EddyFlowProj%err_label, gain=1d-3, offset=0d0)
    call AddFloatDatumToDataline(lEx%RHO%d, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%Vd, csv_row, EddyFlowProj%err_label)
    !> Specific heat of evaporation
    call AddFloatDatumToDataline(lEx%lambda, csv_row, EddyFlowProj%err_label)
    !> Wet to dry air density ratio
    call AddFloatDatumToDataline(lEx%sigma, csv_row, EddyFlowProj%err_label)

    !> Gas concentrations/densities
    do gas = co2, gas4
        call AddIntDatumToDataline(lEx%measure_type_int(gas), csv_row, EddyFlowProj%err_label)
        call AddFloatDatumToDataline(lEx%d(gas), csv_row, EddyFlowProj%err_label)
        if (gas == ch4 .or. gas == gas4) then
            call AddFloatDatumToDataline(lEx%r(gas), csv_row, EddyFlowProj%err_label, gain=1d3, offset=0d0)
            call AddFloatDatumToDataline(lEx%chi(gas), csv_row, EddyFlowProj%err_label, gain=1d3, offset=0d0)
        else
            call AddFloatDatumToDataline(lEx%r(gas), csv_row, EddyFlowProj%err_label)
            call AddFloatDatumToDataline(lEx%chi(gas), csv_row, EddyFlowProj%err_label)
        end if
    end do

    !> Time lags
    do gas = co2, gas4
        call AddFloatDatumToDataline(lEx%act_tlag(gas), csv_row, EddyFlowProj%err_label)
        call AddFloatDatumToDataline(lEx%used_tlag(gas), csv_row, EddyFlowProj%err_label)
        call AddFloatDatumToDataline(lEx%nom_tlag(gas), csv_row, EddyFlowProj%err_label)
        call AddFloatDatumToDataline(lEx%min_tlag(gas), csv_row, EddyFlowProj%err_label)
        call AddFloatDatumToDataline(lEx%max_tlag(gas), csv_row, EddyFlowProj%err_label)
    end do

    !> Stats
    do var = u, gas4
        if (var == ts) then
            call AddFloatDatumToDataline(lEx%stats%median(var), csv_row, &
                EddyFlowProj%err_label, gain=1d0, offset=-273.15d0)
        else
            call AddFloatDatumToDataline(lEx%stats%median(var), csv_row, EddyFlowProj%err_label)
        end if
    end do
    do var = u, gas4
        if (var == ts) then
            call AddFloatDatumToDataline(lEx%stats%Q1(var), csv_row, &
                EddyFlowProj%err_label, gain=1d0, offset=-273.15d0)
        else
            call AddFloatDatumToDataline(lEx%stats%Q1(var), csv_row, EddyFlowProj%err_label)
        end if
    end do
    do var = u, gas4
        if (var == ts) then
            call AddFloatDatumToDataline(lEx%stats%Q3(var), csv_row, &
                EddyFlowProj%err_label, gain=1d0, offset=-273.15d0)
        else
            call AddFloatDatumToDataline(lEx%stats%Q3(var), csv_row, EddyFlowProj%err_label)
        end if
    end do
    do var = u, gas4
        call AddFloatDatumToDataline(sqrt(lEx%stats%Cov(var, var)), csv_row, EddyFlowProj%err_label)
    end do
    do var = u, gas4
        call AddFloatDatumToDataline(lEx%stats%Skw(var), csv_row, EddyFlowProj%err_label)
    end do
    do var = u, gas4
        call AddFloatDatumToDataline(lEx%stats%Kur(var), csv_row, EddyFlowProj%err_label)
    end do
    call AddFloatDatumToDataline(lEx%stats%Cov(w, u), csv_row, EddyFlowProj%err_label)
    do var = ts, gas4
        call AddFloatDatumToDataline(lEx%stats%Cov(w, var), csv_row, EddyFlowProj%err_label)
    end do
    do var = h2o, gas4
        call AddFloatDatumToDataline(lEx%stats%Cov(co2, var), csv_row, EddyFlowProj%err_label)
    end do
    do var = ch4, gas4
        call AddFloatDatumToDataline(lEx%stats%Cov(h2o, var), csv_row, EddyFlowProj%err_label)
    end do
    call AddFloatDatumToDataline(lEx%stats%Cov(ch4, gas4), csv_row, EddyFlowProj%err_label)

    !> Footprint
    call AddFloatDatumToDataline(Foot%peak, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Foot%offset, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Foot%x10, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Foot%x30, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Foot%x50, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Foot%x70, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Foot%x80, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Foot%x90, csv_row, EddyFlowProj%err_label)

    !> Fluxes Level 0 (uncorrected)
    call AddFloatDatumToDataline(lEx%Flux0%ustar, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%Flux0%L, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%Flux0%zL, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%Flux0%Tau, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%Flux0%H, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%Flux0%LE, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%Flux0%ET, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%Flux0%co2, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%Flux0%h2o, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%Flux0%ch4, csv_row, EddyFlowProj%err_label, gain=1d3, offset=0d0)
    call AddFloatDatumToDataline(lEx%Flux0%gas4, csv_row, EddyFlowProj%err_label, gain=1d3, offset=0d0)
    !> Fluxes Level 1 
    call AddFloatDatumToDataline(Flux1%Tau, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Flux1%H, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Flux1%LE, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Flux1%ET, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Flux1%co2, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Flux1%h2o, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Flux1%ch4, csv_row, EddyFlowProj%err_label, gain=1d3, offset=0d0)
    call AddFloatDatumToDataline(Flux1%gas4, csv_row, EddyFlowProj%err_label, gain=1d3, offset=0d0)
    !> Fluxes Level 2
    call AddFloatDatumToDataline(Flux2%Tau, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Flux2%H, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Flux2%LE, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Flux2%ET, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Flux2%co2, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Flux2%h2o, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Flux2%ch4, csv_row, EddyFlowProj%err_label, gain=1d3, offset=0d0)
    call AddFloatDatumToDataline(Flux2%gas4, csv_row, EddyFlowProj%err_label, gain=1d3, offset=0d0)

    !> Cell values
    call AddFloatDatumToDataline(lEx%Tcell, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%Pcell, csv_row, EddyFlowProj%err_label)
    do gas = co2, gas4
        call AddFloatDatumToDataline(lEx%Vcell(gas), csv_row, EddyFlowProj%err_label)
    end do
    call AddFloatDatumToDataline(lEx%Flux0%E_co2, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%Flux0%E_ch4, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%Flux0%E_gas4, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%Flux0%Hi_co2, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%Flux0%Hi_h2o, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%Flux0%Hi_ch4, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%Flux0%Hi_gas4, csv_row, EddyFlowProj%err_label)

    !> Burba terms
    if (lEx%Burba%h_bot + lEx%Burba%h_top + lEx%Burba%h_spar /= 0.0) then
        call AddFloatDatumToDataline(lEx%Burba%h_bot, csv_row, EddyFlowProj%err_label)
        call AddFloatDatumToDataline(lEx%Burba%h_top, csv_row, EddyFlowProj%err_label)
        call AddFloatDatumToDataline(lEx%Burba%h_spar, csv_row, EddyFlowProj%err_label)
    else
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
    end if
    !> LI-7700 multipliers
    call AddFloatDatumToDataline(lEx%Mul7700%A, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%Mul7700%B, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%Mul7700%C, csv_row, EddyFlowProj%err_label)

    !> Spectral correction factors
    call AddFloatDatumToDataline(BPCF%of(w_u), csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(BPCF%of(w_ts), csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(BPCF%of(w_h2o), csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(BPCF%of(w_h2o), csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(BPCF%of(w_co2), csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(BPCF%of(w_h2o), csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(BPCF%of(w_ch4), csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(BPCF%of(w_gas4), csv_row, EddyFlowProj%err_label)

    !> Degraded covariances
    call AddFloatDatumToDataline(lEx%degT%cov, csv_row, EddyFlowProj%err_label)
    do i = 1, 9
        call AddFloatDatumToDataline(lEx%degT%dcov(i), csv_row, EddyFlowProj%err_label)
    end do
    do var = u, gas4
        call AddIntDatumToDataline(lEx%spikes(var), csv_row, EddyFlowProj%err_label)
    end do

    !> Write first string from Chunks
    !> M_CUSTOM_FLAGS thru VM97_NSW_RNS
    call AddDatum(csv_row, trim(fluxnetChunks%s(1)), separator)

    !> VM97 flags, here organized per variable instead of per test
    if (lEx%vm_flags(1) == '-9999') then
        do var = u, gas4
            call AddCharDatumToDataline(EddyFlowProj%err_label, csv_row, EddyFlowProj%err_label)
        end do
    else
        do var = u, gas4
            vi = var + 1
            vm97flags(var)(1 : 1) = '8'
            vm97flags(var)(2 : 2) = lEx%vm_flags(1)(vi:vi)
            vm97flags(var)(3 : 3) = lEx%vm_flags(2)(vi:vi)
            vm97flags(var)(4 : 4) = lEx%vm_flags(3)(vi:vi)
            vm97flags(var)(5 : 5) = lEx%vm_flags(4)(vi:vi)
            vm97flags(var)(6 : 6) = lEx%vm_flags(5)(vi:vi)
            vm97flags(var)(7 : 7) = lEx%vm_flags(6)(vi:vi)
            vm97flags(var)(8 : 8) = lEx%vm_flags(7)(vi:vi)
            vm97flags(var)(9 : 9) = lEx%vm_flags(8)(vi:vi)
            call AddCharDatumToDataline(trim(vm97flags(var)), csv_row, EddyFlowProj%err_label)
        end do
    end if

    !> Uncomment to reintroduce flags for last 3 tests
    call AddCharDatumToDataline(lEx%vm_tlag_hf, csv_row, separator)
    call AddCharDatumToDataline(lEx%vm_tlag_sf, csv_row, separator)
    call AddCharDatumToDataline(lEx%vm_aoa_hf, csv_row, separator)
    call AddCharDatumToDataline(lEx%vm_nshw_hf, csv_row, separator)

    !> Write second string from Chunks
    call AddDatum(csv_row, fluxnetChunks%s(2), separator)

    !> Foken's QC details
    call AddFloatDatumToDataline(lEx%TAU_SS, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%H_SS, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%FC_SS, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%FH2O_SS, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%FCH4_SS, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%FGS4_SS, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%U_ITC, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%W_ITC, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%TS_ITC, csv_row, EddyFlowProj%err_label)

    !> Write second string from Chunks
    call AddDatum(csv_row, fluxnetChunks%s(3), separator)

    !> Foken's final flags
    call AddIntDatumToDataline(QCFlag%tau, csv_row, EddyFlowProj%err_label)
    call AddIntDatumToDataline(QCFlag%H, csv_row, EddyFlowProj%err_label)
    call AddIntDatumToDataline(QCFlag%h2o, csv_row, EddyFlowProj%err_label)
    call AddIntDatumToDataline(QCFlag%h2o, csv_row, EddyFlowProj%err_label)
    call AddIntDatumToDataline(QCFlag%co2, csv_row, EddyFlowProj%err_label)
    call AddIntDatumToDataline(QCFlag%h2o, csv_row, EddyFlowProj%err_label)
    call AddIntDatumToDataline(QCFlag%ch4, csv_row, EddyFlowProj%err_label)
    call AddIntDatumToDataline(QCFlag%gas4, csv_row, EddyFlowProj%err_label)

    !> LI-COR's IRGAs diagnostics breakdown
    do i = 1, 29
        call AddFloatDatumToDataline(lEx%licor_flags(i), csv_row, EddyFlowProj%err_label)
    end do

    !> AGC/RSSI
    call AddFloatDatumToDataline(lEx%agc72, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%agc75, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%rssi77, csv_row, EddyFlowProj%err_label)

    !> Write third string from Chunks
    !> WBOOST_APPLIED thru AXES_ROTATION_METHOD
    call AddDatum(csv_row, fluxnetChunks%s(4), separator)

    !> Rotation angles
    call AddFloatDatumToDataline(lEx%yaw, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%pitch, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%roll, csv_row, EddyFlowProj%err_label)

    !> Detrending method and time constant
    call AddIntDatumToDataline(lEx%det_meth_int, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%det_timec, csv_row, EddyFlowProj%err_label)

    !> Write forth string from Chunks
    !> TIMELAG_DETECTION_METHOD thru FOOTPRINT_MODEL
    call AddDatum(csv_row, fluxnetChunks%s(5), separator)

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
    call AddIntDatumToDataline(lEx%logger_swver%major, csv_row, EddyFlowProj%err_label)
    call AddIntDatumToDataline(lEx%logger_swver%minor, csv_row, EddyFlowProj%err_label)
    call AddIntDatumToDataline(lEx%logger_swver%revision, csv_row, EddyFlowProj%err_label)
    !>> Site info
    call AddFloatDatumToDataline(lEx%lat, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%lon, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%alt, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%canopy_height, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%disp_height, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%rough_length, csv_row, EddyFlowProj%err_label)
    !>> Acquisition setup
    call AddFloatDatumToDataline(lEx%file_length, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%ac_freq, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%avrg_length, csv_row, EddyFlowProj%err_label)
    !>> Master sonic height and north offset
    call AddDatum(csv_row, trim(lEx%instr(sonic)%firm), separator)
    call AddDatum(csv_row, trim(lEx%instr(sonic)%model), separator)
    call AddFloatDatumToDataline(lEx%instr(sonic)%height, csv_row, EddyFlowProj%err_label)
    call AddDatum(csv_row, lEx%instr(sonic)%wformat, separator)
    call AddDatum(csv_row, lEx%instr(sonic)%wref, separator)
    call AddFloatDatumToDataline(lEx%instr(sonic)%north_offset, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%instr(sonic)%hpath_length, csv_row, EddyFlowProj%err_label, gain=1d2, offset=0d0)
    call AddFloatDatumToDataline(lEx%instr(sonic)%vpath_length, csv_row, EddyFlowProj%err_label, gain=1d2, offset=0d0)
    call AddFloatDatumToDataline(lEx%instr(sonic)%tau, csv_row, EddyFlowProj%err_label)

    !>> irgas
    do igas = ico2, igas4
        call AddDatum(csv_row, trim(lEx%instr(igas)%firm), separator)
        call AddDatum(csv_row, trim(lEx%instr(igas)%model), separator)
        call AddFloatDatumToDataline(lEx%instr(igas)%nsep, csv_row, EddyFlowProj%err_label, gain=1d2, offset=0d0)
        call AddFloatDatumToDataline(lEx%instr(igas)%esep, csv_row, EddyFlowProj%err_label, gain=1d2, offset=0d0)
        call AddFloatDatumToDataline(lEx%instr(igas)%vsep, csv_row, EddyFlowProj%err_label, gain=1d2, offset=0d0)
        call AddFloatDatumToDataline(lEx%instr(igas)%tube_l, csv_row, EddyFlowProj%err_label, gain=1d2, offset=0d0)
        call AddFloatDatumToDataline(lEx%instr(igas)%tube_d, csv_row, EddyFlowProj%err_label, gain=1d3, offset=0d0)
        call AddFloatDatumToDataline(lEx%instr(igas)%tube_f, csv_row, EddyFlowProj%err_label, gain=6d4, offset=0d0)
        if (igas == ih2o) then
            call AddFloatDatumToDataline(lEx%instr(igas)%kw, csv_row, EddyFlowProj%err_label)
            call AddFloatDatumToDataline(lEx%instr(igas)%ko, csv_row, EddyFlowProj%err_label)
        end if
        call AddFloatDatumToDataline(lEx%instr(igas)%hpath_length, csv_row, EddyFlowProj%err_label, gain=1d2, offset=0d0)
        call AddFloatDatumToDataline(lEx%instr(igas)%vpath_length, csv_row, EddyFlowProj%err_label, gain=1d2, offset=0d0)
        call AddFloatDatumToDataline(lEx%instr(igas)%tau, csv_row, EddyFlowProj%err_label)
    end do

    !> Custom variables

    call AddIntDatumToDataline(lEx%ncustom, csv_row, EddyFlowProj%err_label)
    if (lEx%ncustom > 0) then
        do i = 1, lEx%ncustom
            call AddFloatDatumToDataline(lEx%user_var(i), csv_row, EddyFlowProj%err_label)
        end do
    end if

    !> Write sisxth string from Chunks
    !> Biomet data
    call AddDatum(csv_row, fluxnetChunks%s(6), separator)

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

end subroutine WriteOutFluxnetFcc
