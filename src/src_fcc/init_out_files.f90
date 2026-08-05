!***************************************************************************
! init_out_files.f90
! ------------------
! Copyright © 2007-2011, Eco2s team, Gerardo Fratini
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
! \brief       Initializes output files
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine InitOutFiles(lEx)
    use m_fx_global_var
    implicit none
    !> in/out variables
    Type(ExType), intent(in) :: lEx
    !> local variables
    integer :: mkdir_status
    integer :: open_status
    integer :: dot
    integer :: gas
    integer :: i
    integer :: k
    character(PathLen) :: Test_Path
    character(64) :: e2sg(E2NumVar)
    character(64) :: gas_tag(GHGNumVar)
    !> Full-output layout: the gas slots this file carries a block for.
    !> WriteOutFullFcc walks the same list.
    integer :: fo_slots(GHGNumVar)
    integer :: n_fo_slots
    !> Statistical-flag legends, shared with WriteOutFullFcc.
    character(LongOutstringLen) :: flag_legend
    character(LongOutstringLen) :: tl_legend
    integer :: n_flag_vars
    integer :: n_tl_vars
    character(LongOutstringLen) :: header1
    character(LongOutstringLen) :: header2
    character(LongOutstringLen) :: header3
    character(64) :: custom_label
    character(32) :: custom_unit
    character(2) :: utf8_mu
    integer, external :: CreateDir
    include '../src_common/interfaces_1.inc'

    utf8_mu = char(194) // char(181)

    e2sg(u)    = 'u_'
    e2sg(v)    = 'v_'
    e2sg(w)    = 'w_'
    e2sg(ts)   = 'ts_'
    !> Every configured gas gets a name, from the same helper RP uses, so the
    !> two halves of the full output cannot name the same column differently.
    call FullOutputGasTags(gas_tag)
    do gas = firstGas, lastGas
        e2sg(gas) = gas_tag(gas)
    end do

    !> The gas slots the full output carries a block for, from the helper the
    !> row writer also calls. The fixed format names four blocks and the row
    !> loop emitted sixty-four; RP's twin had the same split.
    call FullOutputGasSlots(fo_slots, n_fo_slots)

    call StatisticalFlagVars(n_flag_vars, flag_legend)
    call TimelagFlagLegend(n_tl_vars, tl_legend)

    !> Full output file
    if (EddyFlowProj%out_full) then
        !> Create output directory if it does not exist
        mkdir_status = CreateDir('"' // Dir%main_out(1:len_trim(Dir%main_out)) // '"')

        !> Open full output file and writes header
        Test_Path = Dir%main_out(1:len_trim(Dir%main_out)) &
                  // EddyFlowProj%id(1:len_trim(EddyFlowProj%id)) &
                  // FullOut_FilePadding // Timestamp_FilePadding // CsvExt
        dot = index(Test_Path, CsvExt, .true.) - 1
        FullOut_Path = Test_Path(1:dot) // CsvTmpExt
        open(uflx, file = FullOut_Path, iostat = open_status, encoding = 'utf-8')

        !> Initialize header strings to void
        call clearstr(header1)
        call clearstr(header2)
        call clearstr(header3)
        !> Initial file and timestamp info
        call AddDatum(header1,'file_info,,,,,,', separator)
        call AddDatum(header2,'filename,date,time,DOY,daytime,file_records,used_records', separator)
        call AddDatum(header3,',[yyyy-mm-dd],[HH:MM],[ddd.ddd],[1=daytime],[#],[#]', separator)

        !> Corrected fluxes (Level 3) and quality flags
        !> Tau
        call AddDatum(header1, 'corrected_fluxes_and_quality_flags,', separator)
        call AddDatum(header2,'Tau,qc_Tau', separator)
        call AddDatum(header3,'[kg+1m-1s-2],[#]', separator)
        if (RUsetup%meth /= 'none') then
            call AddDatum(header1, '', separator)
            call AddDatum(header2,'rand_err_Tau', separator)
            call AddDatum(header3,'[kg+1m-1s-2]', separator)
        end if

        !> H
        call AddDatum(header1, ',', separator)
        call AddDatum(header2, 'H,qc_H', separator)
        call AddDatum(header3, '[W+1m-2],[#]', separator)
        if (RUsetup%meth /= 'none') then
            call AddDatum(header1, '', separator)
            call AddDatum(header2, 'rand_err_H', separator)
            call AddDatum(header3, '[W+1m-2]', separator)
        end if

        !> LE
        if(fcc_var_present(PrimaryWaterOutSlot())) then
            call AddDatum(header1, ',', separator)
            call AddDatum(header2, 'LE,qc_LE', separator)
            call AddDatum(header3, '[W+1m-2],[#]', separator)
            if (RUsetup%meth /= 'none') then
                call AddDatum(header1, '', separator)
                call AddDatum(header2, 'rand_err_LE', separator)
                call AddDatum(header3, '[W+1m-2]', separator)
            end if
        end if

        !> Corrected co2 fluxes
        !> Corrected gas fluxes, one block per configured gas.
        do k = 1, n_fo_slots
            gas = fo_slots(k)
            if(.not. fcc_var_present(gas)) cycle
            call AddDatum(header1, ',', separator)
            call AddDatum(header2, e2sg(gas)(1:len_trim(e2sg(gas))) &
                // 'flux,qc_' // e2sg(gas)(1:len_trim(e2sg(gas))) // 'flux', separator)
            call AddDatum(header3, trim(gas_full_flux_label(gas)) // ',[#]', separator)
            if (RUsetup%meth /= 'none') then
                call AddDatum(header1, '', separator)
                call AddDatum(header2, 'rand_err_' // e2sg(gas)(1:len_trim(e2sg(gas))) // 'flux', separator)
                call AddDatum(header3, gas_full_flux_label(gas), separator)
            end if
        end do

        !> Storage
        call AddDatum(header1, 'storage_fluxes', separator)
        call AddDatum(header2,'H_strg', separator)
        call AddDatum(header3,'[W+1m-2]', separator)
        if(fcc_var_present(PrimaryWaterOutSlot())) call AddDatum(header1, '', separator)
        if(fcc_var_present(PrimaryWaterOutSlot())) call AddDatum(header2,'LE_strg', separator)
        if(fcc_var_present(PrimaryWaterOutSlot())) call AddDatum(header3,'[W+1m-2]', separator)
        do k = 1, n_fo_slots
            gas = fo_slots(k)
            if(.not. fcc_var_present(gas)) cycle
            call AddDatum(header1, '', separator)
            call AddDatum(header2, e2sg(gas)(1:len_trim(e2sg(gas))) // 'strg', separator)
            call AddDatum(header3, gas_full_flux_label(gas), separator)
        end do

        !> Advection fluxes
        header1 = header1(1:len_trim(header1)) // 'vertical_advection_fluxes'
        do k = 1, n_fo_slots
            gas = fo_slots(k)
            if(.not. fcc_var_present(gas)) cycle
            call AddDatum(header1, '', separator)
            call AddDatum(header2, e2sg(gas)(1:len_trim(e2sg(gas))) // 'v-adv', separator)
            call AddDatum(header3, gas_full_flux_label(gas), separator)
        end do

        !> Average gas concentrations
        call AddDatum(header1,'gas_densities_concentrations_and_timelags', separator)
        do k = 1, n_fo_slots
            gas = fo_slots(k)
            if(.not. fcc_var_present(gas)) cycle
            call AddDatum(header1, ',,,,', separator)
            call AddDatum(header2, e2sg(gas)(1:len_trim(e2sg(gas))) // 'molar_density,' &
                // e2sg(gas)(1:len_trim(e2sg(gas))) // 'mole_fraction,' &
                // e2sg(gas)(1:len_trim(e2sg(gas))) // 'mixing_ratio,' &
                // e2sg(gas)(1:len_trim(e2sg(gas))) // 'time_lag,' &
                // e2sg(gas)(1:len_trim(e2sg(gas))) // 'def_timelag', separator)
            call AddDatum(header3, &
                trim(gas_full_dens_label(gas)) // ',' // trim(gas_full_conc_label(gas)) &
                // ',' // trim(gas_full_mixr_label(gas)) // ',[s],[1=default]', separator)
        end do
        !> In Header 1 there is one comma too much, take it away
        header1 = header1(1:len_trim(header1) - 1)

        !> Air properties, wind components and rotation angles
        call AddDatum(header1, 'air_properties,,,,,,,,,,,,,,unrotated_wind,,,rotated_wind&
                      &,,,,,,rotation_angles_for_tilt_correction,,', separator)
        call AddDatum(header2,'sonic_temperature,air_temperature,air_pressure,air_density,air_heat_capacity,air_molar_volume,&
                      &ET,water_vapor_density,e,es,specific_humidity,RH,VPD,Tdew&
                      &,u_unrot,v_unrot,w_unrot,u_rot,v_rot,w_rot,wind_speed,max_wind_speed,wind_dir,yaw,pitch,roll', separator)
        call AddDatum(header3,'[K],[K],[Pa],[kg+1m-3],[J+1kg-1K-1],[m+3mol-1],&
                      &[mm+1hour-1],[kg+1m-3],[Pa],[Pa],[kg+1kg-1],[%],[Pa],[K],&
                      &[m+1s-1],[m+1s-1],[m+1s-1],[m+1s-1],[m+1s-1],&
                      &[m+1s-1],[m+1s-1],[m+1s-1],[deg_from_north],[deg],[deg],[deg]', separator)

        !> Turbulence
        call AddDatum(header1, 'turbulence,,,,,', separator)
        call AddDatum(header2,'u*,TKE,L,(z-d)/L,bowen_ratio,T*', separator)
        call AddDatum(header3,'[m+1s-1],[m+2s-2],[m],[#],[#],[K]', separator)

        !> Footprint, if requested
        if (Meth%foot /= 'none') then
            call AddDatum(header1, 'footprint,,,,,,,', separator)
            call AddDatum(header2,'model,x_peak,x_offset,x_10%,x_30%,x_50%,x_70%,x_90%', separator)
            call AddDatum(header3,'[0=KJ/1=KM/2=HS],[m],[m],[m],[m],[m],[m],[m]', separator)
        end if

        !> uncorrected fluxes
        !> Tau and H
        call AddDatum(header1, 'uncorrected_fluxes,,,', separator)
        call AddDatum(header2,'un_Tau,Tau_scf,un_H,H_scf', separator)
        call AddDatum(header3,'[kg+1m-1s-2],[#],[W+1m-2],[#]', separator)
        !> LE
        if(fcc_var_present(PrimaryWaterOutSlot())) call AddDatum(header1, ',', separator)
        if(fcc_var_present(PrimaryWaterOutSlot())) call AddDatum(header2,'un_LE,LE_scf', separator)
        if(fcc_var_present(PrimaryWaterOutSlot())) call AddDatum(header3,'[W+1m-2],[#]', separator)
        !> Uncorrected gas fluxes (Level 0)
        do k = 1, n_fo_slots
            gas = fo_slots(k)
            if(.not. fcc_var_present(gas)) cycle
            call AddDatum(header1, ',', separator)
            call AddDatum(header2, 'un_' // e2sg(gas)(1:len_trim(e2sg(gas))) &
                // 'flux,' // e2sg(gas)(1:len_trim(e2sg(gas))) // 'scf', separator)
            call AddDatum(header3, trim(gas_full_flux_label(gas)) // ',[#]', separator)
        end do

        !> Vickers and Mahrt 97 hard and soft flags
        call AddDatum(header1,'statistical_flags,,,,,,,,,,,', separator)
        call AddDatum(header2,'spikes_hf,amplitude_resolution_hf,drop_out_hf,absolute_limits_hf,&
            &skewness_kurtosis_hf,skewness_kurtosis_sf,discontinuities_hf,discontinuities_sf,timelag_hf,&
            &timelag_sf,attack_angle_hf,non_steady_wind_hf', separator)
        call AddDatum(header3, &
            '8' // trim(flag_legend) // ',8' // trim(flag_legend) &
            // ',8' // trim(flag_legend) // ',8' // trim(flag_legend) &
            // ',8' // trim(flag_legend) // ',8' // trim(flag_legend) &
            // ',8' // trim(flag_legend) // ',8' // trim(flag_legend) &
            // ',8' // trim(tl_legend) // ',8' // trim(tl_legend) &
            // ',8aa,8U', separator)

        !> Add spikes for EddyFlow variables
        call AddDatum(header1,'spikes,,,', separator)
        call AddDatum(header2,'u_spikes,v_spikes,w_spikes,ts_spikes', separator)
        call AddDatum(header3,'[#],[#],[#],[#]', separator)
        do k = 1, n_fo_slots
            gas = fo_slots(k)
            if(fcc_var_present(gas)) then
                call AddDatum(header1, '', separator)
                call AddDatum(header2, e2sg(gas)(1:len_trim(e2sg(gas))) // 'spikes' , separator)
                call AddDatum(header3, '[#]', separator)
            end if
        end do

       !> LI-COR's diagnostic flags
        if (Diag7200%present) then
            call AddDatum(header1,'diagnostic_flags_LI-7200,,,,,,,,', separator)
            call AddDatum(header2,'head_detect_LI-7200,t_out_LI-7200,t_in_LI-7200,aux_in_LI-7200,delta_p_LI-7200,&
                &chopper_LI-7200,detector_LI-7200,pll_LI-7200,sync_LI-7200', separator)
            call AddDatum(header3,'[#_flagged_recs],[#_flagged_recs],[#_flagged_recs],[#_flagged_recs],[#_flagged_recs],&
                &[#_flagged_recs],[#_flagged_recs],[#_flagged_recs],[#_flagged_recs]', separator)
        end if
        if (Diag7500%present) then
            call AddDatum(header1,'diagnostic_flags_LI-7500,,,', separator)
            call AddDatum(header2,'chopper_LI-7500,detector_LI-7500,pll_LI-7500,sync_LI-7500', separator)
            call AddDatum(header3,'[#_flagged_recs],[#_flagged_recs],[#_flagged_recs],[#_flagged_recs]', separator)
        end if
        if (Diag7700%present) then
            call AddDatum(header1,'diagnostic_flags_LI-7700,,,,,,,,,,,,,,,', separator)
            call AddDatum(header2,'not_ready_LI-7700,no_signal_LI-7700,re_unlocked_LI-7700,bad_temp_LI-7700,&
                &laser_temp_unregulated_LI-7700,block_temp_unregulated_LI-7700,motor_spinning_LI-7700,&
                &pump_on_LI-7700,top_heater_on_LI-7700,bottom_heater_on_LI-7700,calibrating_LI-7700,&
                &motor_failure_LI-7700,bad_aux_tc1_LI-7700,bad_aux_tc2_LI-7700,&
                &bad_aux_tc3_LI-7700,box_connected_LI-7700', separator)
            call AddDatum(header3,'[#_flagged_recs],[#_flagged_recs],[#_flagged_recs],[#_flagged_recs],[#_flagged_recs],&
                &[#_flagged_recs],[#_flagged_recs],[#_flagged_recs],[#_flagged_recs],[#_flagged_recs],[#_flagged_recs],&
                &[#_flagged_recs],[#_flagged_recs],[#_flagged_recs],[#_flagged_recs],[#_flagged_recs]', separator)
        end if

        !> AGCs and RSSIs for LI-7200 and LI-7500
        if (Diag7200%present) then
            if(co2_new_sw_ver) then
                call AddDatum(header1,'RSSI_LI-7200', separator)
                call AddDatum(header2,'mean_value_RSSI_LI-7200', separator)
                call AddDatum(header3,'[#]', separator)
            else
                call AddDatum(header1,'AGC_LI-7200', separator)
                call AddDatum(header2,'mean_value_AGC_LI-7200', separator)
                call AddDatum(header3,'[#]', separator)
            end if
        end if
        if (Diag7500%present) then
            if(co2_new_sw_ver) then
                call AddDatum(header1,'RSSI_LI-7500', separator)
                call AddDatum(header2,'mean_value_RSSI_LI-7500', separator)
                call AddDatum(header3,'[#]', separator)
            else
                call AddDatum(header1,'AGC_LI-7500', separator)
                call AddDatum(header2,'mean_value_AGC_LI-7500', separator)
                call AddDatum(header3,'[#]', separator)
            end if
        end if

        !> Variances
        call AddDatum(header1, 'variances,,,', separator)
        call AddDatum(header2, 'u_var,v_var,w_var,ts_var', separator)
        call AddDatum(header3, '[m+2s-2],[m+2s-2],[m+2s-2],[K+2]', separator)
        do k = 1, n_fo_slots
            gas = fo_slots(k)
            if(fcc_var_present(gas)) call AddDatum(header1, '', separator)
            if(fcc_var_present(gas)) call AddDatum(header2, e2sg(gas)(1:len_trim(e2sg(gas))) // 'var', separator)
            if(fcc_var_present(gas)) call AddDatum(header3, '--', separator)
        end do
        !> w/ts covariance
        call AddDatum(header1, 'covariances', separator)
        call AddDatum(header2,'w/ts_cov', separator)
        call AddDatum(header3,'[m+1K+1s-1]', separator)
        !> w-gases covariances
        do k = 1, n_fo_slots
            gas = fo_slots(k)
            if(fcc_var_present(gas)) call AddDatum(header1, '', separator)
            if(fcc_var_present(gas)) call AddDatum(header2, 'w/' // e2sg(gas)(1:len_trim(e2sg(gas))) // 'cov', separator)
            if(fcc_var_present(gas)) call AddDatum(header3, '--', separator)
        end do

        !> Mean values of user variables
        if (lEx%ncustom > 0) then
            call AddDatum(header1, 'custom_variables', separator)
            do i = 1, lEx%ncustom
                if (i > 1) call AddDatum(header1, '', separator)
                call clearstr(custom_label)
                call clearstr(custom_unit)
                if (i <= MaxUserVar) custom_label = UserVarHeader(i)
                if (len_trim(custom_label) == 0) write(custom_label, '("custom_", i0, "_mean")') i
                custom_unit = CustomUnitFromLabel(custom_label)
                call AddDatum(header2, custom_label(1:len_trim(custom_label)), separator)
                call AddDatum(header3, custom_unit(1:len_trim(custom_unit)), separator)
            end do
        end if

        !> Conditional Eddy Covariance outputs (Zahn et al. 2022)
        if (EddyFlowProj%do_cec == 1 .or. EddyFlowProj%do_cec == 2) then
            call AddDatum(header1, 'conditional_eddy_covariance_(H2O),,,,,', separator)
            call AddDatum(header2, 'E_cec,Tr_cec,E_cec_ET,Tr_cec_ET,r_ET_cec,qc_cec_h2o', separator)
            call AddDatum(header3, &
                '[mmol+1m-2s-1],[mmol+1m-2s-1],[mm+1hour-1],[mm+1hour-1],[#],[#]', separator)
        end if
        if (EddyFlowProj%do_cec == 1 .or. EddyFlowProj%do_cec == 3) then
            call AddDatum(header1, 'conditional_eddy_covariance_(CO2),,,,', separator)
            call AddDatum(header2, 'Reco_cec,P_cec,NEE_cec,r_Fc_cec,qc_cec_co2', separator)
            call AddDatum(header3, '[umol+1m-2s-1],[umol+1m-2s-1],[umol+1m-2s-1],[#],[#]', separator)
        end if

        !> Write on output file
        write(uflx, '(a)') header1(1:len_trim(header1) - 1)
        write(uflx, '(a)') header2(1:len_trim(header2) - 1)
        write(uflx, '(a)') header3(1:len_trim(header3) - 1)
    end if

    !************************************************************************************************************************************
    !************************************************************************************************************************************
    !> METADATA file
    if (EddyFlowProj%out_md) then
        !> Create metadata output file name
        !> Open dynamic matadata file
        Test_Path = Dir%main_out(1:len_trim(Dir%main_out)) &
                  // EddyFlowProj%id(1:len_trim(EddyFlowProj%id)) &
                  // MetaData_FilePadding // Timestamp_FilePadding // CsvExt
        dot = index(Test_Path, CsvExt, .true.) - 1
        Metadata_Path = Test_Path(1:dot) // CsvTmpExt
        open(umd, file = Metadata_Path, iostat = open_status, encoding = 'utf-8')
        call Clearstr(header1)
        call AddDatum(header1,'filename,date,time,latitude,longitude,altitude,canopy_height,displacement_height,&
            &roughness_length,file_length,acquisition_frequency,&
            &master_sonic_manufacturer,master_sonic_model,master_sonic_height,&
            &master_sonic_wformat,master_sonic_wref,master_sonic_north_offset,&
            &master_sonic_hpath_length,master_sonic_vpath_length,master_sonic_tau', separator)
        !> One block per configured gas, named from the same e2sg tags the
        !> full output uses. These were three literals and a fourth built by
        !> concatenation, so a gas past the fourth had no metadata columns -
        !> while WriteOutMetadataFcc, once widened, emits fourteen fields for
        !> it. Header and row move together or the file shifts.
        do gas = firstGas, lastGas
            if (.not. fcc_var_present(gas)) cycle
            call AddDatum(header1, &
                  trim(e2sg(gas)) // 'irga_manufacturer,' &
               // trim(e2sg(gas)) // 'irga_model,' &
               // trim(e2sg(gas)) // 'measure_type,' &
               // trim(e2sg(gas)) // 'irga_northward_separation,' &
               // trim(e2sg(gas)) // 'irga_eastward_separation,' &
               // trim(e2sg(gas)) // 'irga_vertical_separation,' &
               // trim(e2sg(gas)) // 'irga_tube_length,' &
               // trim(e2sg(gas)) // 'irga_tube_diameter,' &
               // trim(e2sg(gas)) // 'irga_tube_flowrate,' &
               // trim(e2sg(gas)) // 'irga_kw,' &
               // trim(e2sg(gas)) // 'irga_ko,' &
               // trim(e2sg(gas)) // 'irga_hpath_length,' &
               // trim(e2sg(gas)) // 'irga_vpath_length,' &
               // trim(e2sg(gas)) // 'irga_tau', separator)
        end do
        write(umd, '(a)') header1(1:len_trim(header1) - 1)
    end if

    !***************************************************************************
    !***************************************************************************

    if (EddyFlowProj%out_fluxnet) then
        !> Create output directory if it does not exist
        mkdir_status = CreateDir('"' // Dir%main_out(1:len_trim(Dir%main_out)) // '"')

        Test_Path = Dir%main_out(1:len_trim(Dir%main_out)) &
                  // EddyFlowProj%id(1:len_trim(EddyFlowProj%id)) &
                  // FLUXNET_FilePadding // Timestamp_FilePadding // CsvExt
        dot = index(Test_Path, CsvExt, .true.) - 1
        FLUXNET_Path = Test_Path(1:dot) // CsvTmpExt
        open(uflxnt, file = FLUXNET_Path, iostat = open_status, encoding = 'utf-8')

        !> Write on output file
        write(uflxnt, '(a)') trim(fluxnet_header)
    end if

contains

function CustomUnitFromLabel(label) result(unit_label)
    character(*), intent(in) :: label
    character(32) :: unit_label
    character(64) :: clean_label

    unit_label = '--'
    clean_label = label
    call lowercase(clean_label)

    if (index(clean_label, 'flowrate_') == 1) then
        unit_label = '[m+3s-1]'
    elseif (index(clean_label, 'co2_') == 1 &
        .or. index(clean_label, 'n2o_') == 1 &
        .or. index(clean_label, 'ch4_') == 1) then
        unit_label = '[' // utf8_mu // 'mol+1mol_a-1]'
    elseif (index(clean_label, 'h2o_') == 1) then
        unit_label = '[mmol+1mol_a-1]'
    elseif (index(clean_label, 'int_t_') == 1 &
        .or. index(clean_label, 'cell_t_') == 1) then
        unit_label = '[K]'
    elseif (index(clean_label, 'int_p_') == 1) then
        unit_label = '[Pa]'
    end if
end function CustomUnitFromLabel

end subroutine InitOutFiles
