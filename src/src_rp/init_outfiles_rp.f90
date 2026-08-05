!***************************************************************************
! init_outfiles_rp.f90
! --------------------
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
! \brief       Initializes EddyFlow output files
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine InitOutFiles_rp()
    use m_rp_global_var
    use iso_fortran_env
    implicit none
    !> in/out variables
    integer, external :: CreateDir
    !> local variables
    integer :: mkdir_status = 1      ! initializing to false
    integer :: open_status = 1      ! initializing to false
    integer :: dot
    integer :: var
    integer :: gas
    integer :: i
    integer :: j
    integer :: k
    character(PathLen) :: Test_Path
    character(64) :: e2sg(E2NumVar)
    character(64) :: gas_tag(GHGNumVar)
    !> Full-output layout: the gas slots this file carries a block for. The
    !> row writer walks the same list, so the two cannot disagree.
    integer :: fo_slots(GHGNumVar)
    integer :: n_fo_slots
    !> Statistical-flag legends, shared with WriteOutFull so the units row
    !> names exactly the digits the data row carries.
    character(LongOutstringLen) :: flag_legend
    character(LongOutstringLen) :: tl_legend
    integer :: n_flag_vars
    integer :: n_tl_vars
    !> Statistics-file layout: the slots st1..st7 carry, and the name each
    !> one goes by in their header.
    integer :: st_slots(E2NumVar)
    integer :: n_st_slots
    integer :: st_gas(GHGNumVar)
    integer :: n_st_gas
    character(64) :: st_name(E2NumVar)
    character(LongOutstringLen) :: stats_header
    character(64) :: usg(NumUserVar)
    character(64) :: user_header(NumUserVar)
    character(32) :: user_unit(NumUserVar)
    character(32) :: gas_flux_label(GHGNumVar), gas_conc_label(GHGNumVar)
    character(32) :: gas_mixr_label(GHGNumVar), gas_dens_label(GHGNumVar)
    real(kind = dbl) :: gas_flux_sc(GHGNumVar), gas_dens_sc(GHGNumVar)
    character(2) :: utf8_mu
    character(LongOutstringLen) :: header1
    character(LongOutstringLen) :: header2
    character(LongOutstringLen) :: header3
    character(LongOutstringLen) :: dataline
    logical :: proceed
    include '../src_common/interfaces.inc'


    utf8_mu = char(194) // char(181)

    !> Convenient strings
    e2sg(u)   = 'u_'
    e2sg(v)   = 'v_'
    e2sg(w)   = 'w_'
    e2sg(ts)  = 'ts_'
    !> Every configured gas gets a name, not just the historical four. The
    !> stems reproduce 'co2_'/'h2o_'/'ch4_' exactly, because those slots'
    !> records name those species and the label resolves to the metadata
    !> column's variable, which is what the hard-coded strings spelled.
    call FullOutputGasTags(gas_tag)
    do gas = firstGas, lastGas
        e2sg(gas) = gas_tag(gas)
    end do
    e2sg(tc)  = 'cell_t_'
    e2sg(ti1) = 'inlet_t_'
    e2sg(ti2) = 'outlet_t_'
    e2sg(pi)  = 'cell_p_'
    e2sg(te)  = 'air_t_'
    e2sg(pe)  = 'air_p_'

    !> Unit basis per gas, from the unit the project declares for it. Water
    !> keeps its own mmol labels; every other species asks GasFullOutputUnits,
    !> whose default arm is the umol basis the non-water slots used to have
    !> hard-coded, so a four-gas project is unchanged.
    call GasFullOutputUnitsAll(gas_flux_sc, gas_dens_sc, &
        gas_flux_label, gas_conc_label, gas_mixr_label, gas_dens_label)

    !> The gas slots the full output carries a block for. WriteOutFull walks
    !> the same list; when the two were spelled out separately the row ran
    !> sixty gas blocks past its own header.
    call StatisticalFlagVars(n_flag_vars, flag_legend)
    call TimelagFlagLegend(n_tl_vars, tl_legend)

    call FullOutputGasSlots(fo_slots, n_fo_slots)

    do j = 1, NumUserVar
        user_header(j) = FullOutputCustomLabel(j)
        usg(j) = user_header(j)
        if (index(usg(j), '_mean') > 0) usg(j) = usg(j)(1:index(usg(j), '_mean') - 1)
        user_unit(j) = '--'
        select case (trim(UserCol(j)%var))
            case ('flowrate')
                user_unit(j) = '[m+3s-1]'
            case ('cell_t', 'int_t_1', 'int_t_2')
                user_unit(j) = '[K]'
            case ('int_p')
                user_unit(j) = '[Pa]'
            case ('co2', 'n2o', 'ch4')
                user_unit(j) = '[' // utf8_mu // 'mol+1mol_a-1]'
            case ('h2o')
                user_unit(j) = '[mmol+1mol_a-1]'
        end select
    end do

    !> Create sub-directory
    !> Stats dir
    proceed = .false.
    do i = 1, 7
        if (RPsetup%out_st(i)) then
            proceed = .true.
            exit
        end if
    end do
    if (proceed) then
        StatsDir = Dir%main_out(1:len_trim(Dir%main_out)) // SubDirStats // slash
        mkdir_status = CreateDir('"' // StatsDir(1:len_trim(StatsDir)) // '"')
    end if
    !> Raw dataset dir
    proceed = .false.
    do i = 1, 7
        if (RPsetup%out_raw(i)) then
            proceed = .true.
            exit
        end if
    end do
    if (proceed) then
        RawDir = Dir%main_out(1:len_trim(Dir%main_out)) // SubDirRaw // slash
        mkdir_status = CreateDir('"' // RawDir(1:len_trim(RawDir)) // '"')
        !> Create subfolders for selected outputs
        if (RPsetup%out_raw(1)) then
            RawSubDir(1) = RawDir(1:len_trim(RawDir)) // 'level_1' // slash
            mkdir_status = CreateDir('"' // RawSubDir(1)(1:len_trim(RawSubDir(1))) // '"')
        end if
        if (RPsetup%out_raw(2)) then
            RawSubDir(2) = RawDir(1:len_trim(RawDir)) // 'level_2' // slash
            mkdir_status = CreateDir('"' // RawSubDir(2)(1:len_trim(RawSubDir(2))) // '"')
        end if
        if (RPsetup%out_raw(3)) then
            RawSubDir(3) = RawDir(1:len_trim(RawDir)) // 'level_3' // slash
            mkdir_status = CreateDir('"' // RawSubDir(3)(1:len_trim(RawSubDir(3))) // '"')
        end if
        if (RPsetup%out_raw(4)) then
            RawSubDir(4) = RawDir(1:len_trim(RawDir)) // 'level_4' // slash
            mkdir_status = CreateDir('"' // RawSubDir(4)(1:len_trim(RawSubDir(4))) // '"')
        end if
        if (RPsetup%out_raw(5)) then
            RawSubDir(5) = RawDir(1:len_trim(RawDir)) // 'level_5' // slash
            mkdir_status = CreateDir('"' // RawSubDir(5)(1:len_trim(RawSubDir(5))) // '"')
        end if
        if (RPsetup%out_raw(6)) then
            RawSubDir(6) = RawDir(1:len_trim(RawDir)) // 'level_6' // slash
            mkdir_status = CreateDir('"' // RawSubDir(6)(1:len_trim(RawSubDir(6))) // '"')
        end if
        if (RPsetup%out_raw(7)) then
            RawSubDir(7) = RawDir(1:len_trim(RawDir)) // 'level_7' // slash
            mkdir_status = CreateDir('"' // RawSubDir(7)(1:len_trim(RawSubDir(7))) // '"')
        end if
    end if

    !> Binned cospectral dir
    if (RPsetup%out_bin_sp) then
        BinCospectraDir = Dir%main_out(1:len_trim(Dir%main_out)) // SubDirBinCospectra // slash
        mkdir_status = CreateDir('"' // BinCospectraDir(1:len_trim(BinCospectraDir)) // '"')
    end if
    !> Binned ogive dir
    if (RPsetup%out_bin_og) then
        BinOgivesDir = Dir%main_out(1:len_trim(Dir%main_out)) // SubDirBinOgives // slash
        mkdir_status = CreateDir('"' // BinOgivesDir(1:len_trim(BinOgivesDir)) // '"')
    end if
    !> Full cospectra dir
    !> (First determine if at least one full (co)spectrum has to be written on output)
    proceed = .false.
    do var = 1, GHGNumVar
        if (RPsetup%out_full_sp(var) .or. RPsetup%out_full_cosp(var)) then
            proceed = .true.
            exit
        end if
    end do
    if (proceed) then
        CospectraDir = Dir%main_out(1:len_trim(Dir%main_out)) // SubDirCospectra // slash
        mkdir_status = CreateDir('"' // CospectraDir(1:len_trim(CospectraDir)) // '"')
    end if

    !> Open full output file and writes header
    if (EddyFlowProj%out_full) then
        Test_Path = Dir%main_out(1:len_trim(Dir%main_out)) &
                  // EddyFlowProj%id(1:len_trim(EddyFlowProj%id)) &
                  // FullOut_FilePadding // Timestamp_FilePadding // CsvExt
        dot = index(Test_Path, CsvExt, .true.) - 1
        FullOut_Path = Test_Path(1:dot) // CsvTmpExt
        open(uflx, file = FullOut_Path, iostat = open_status, encoding = 'utf-8')

        !> Initialize header strings to void
        call Clearstr(header1)
        call Clearstr(header2)
        call Clearstr(header3)

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
        if(OutVarPresent(PrimaryWaterOutSlot())) then
            call AddDatum(header1, ',', separator)
            call AddDatum(header2, 'LE,qc_LE', separator)
            call AddDatum(header3, '[W+1m-2],[#]', separator)
            if (RUsetup%meth /= 'none') then
                call AddDatum(header1, '', separator)
                call AddDatum(header2, 'rand_err_LE', separator)
                call AddDatum(header3, '[W+1m-2]', separator)
            end if
        end if

        !> Corrected gas fluxes, one block per configured gas. The four
        !> historical blocks were unrolled and only the fourth took its
        !> units from the project, which is why a fifth gas had no columns
        !> here at all rather than merely mislabelled ones.
        do k = 1, n_fo_slots
            gas = fo_slots(k)
            if(.not. OutVarPresent(gas)) cycle
            call AddDatum(header1, ',', separator)
            call AddDatum(header2, e2sg(gas)(1:len_trim(e2sg(gas))) &
                // 'flux,qc_' // e2sg(gas)(1:len_trim(e2sg(gas))) // 'flux', separator)
            call AddDatum(header3, gas_flux_label(gas)(1:len_trim(gas_flux_label(gas))) &
                // ',[#]', separator)
            if (RUsetup%meth /= 'none') then
                call AddDatum(header1, '', separator)
                call AddDatum(header2, 'rand_err_' // e2sg(gas)(1:len_trim(e2sg(gas))) // 'flux', separator)
                call AddDatum(header3, gas_flux_label(gas), separator)
            end if
        end do

        !> Storage
        call AddDatum(header1, 'storage_fluxes', separator)
        call AddDatum(header2,'H_strg', separator)
        call AddDatum(header3,'[W+1m-2]', separator)
        if(OutVarPresent(PrimaryWaterOutSlot())) call AddDatum(header1, '', separator)
        if(OutVarPresent(PrimaryWaterOutSlot())) call AddDatum(header2,'LE_strg', separator)
        if(OutVarPresent(PrimaryWaterOutSlot())) call AddDatum(header3,'[W+1m-2]', separator)
        do k = 1, n_fo_slots
            gas = fo_slots(k)
            if(.not. OutVarPresent(gas)) cycle
            call AddDatum(header1, '', separator)
            call AddDatum(header2, e2sg(gas)(1:len_trim(e2sg(gas))) // 'strg', separator)
            call AddDatum(header3, gas_flux_label(gas), separator)
        end do

        !> Advection fluxes
        header1 = header1(1:len_trim(header1)) // 'vertical_advection_fluxes'
        do k = 1, n_fo_slots
            gas = fo_slots(k)
            if(.not. OutVarPresent(gas)) cycle
            call AddDatum(header1, '', separator)
            call AddDatum(header2, e2sg(gas)(1:len_trim(e2sg(gas))) // 'v-adv', separator)
            call AddDatum(header3, gas_flux_label(gas), separator)
        end do

        !> Average gas concentrations
        call AddDatum(header1,'gas_densities_concentrations_and_timelags', separator)
        do k = 1, n_fo_slots
            gas = fo_slots(k)
            if(OutVarPresent(gas)) call AddDatum(header1, ',,,,', separator)
            if(OutVarPresent(gas)) call AddDatum(header2, e2sg(gas)(1:len_trim(e2sg(gas))) // 'molar_density,' &
                // e2sg(gas)(1:len_trim(e2sg(gas))) // 'mole_fraction,' &
                // e2sg(gas)(1:len_trim(e2sg(gas))) // 'mixing_ratio,' &
                // e2sg(gas)(1:len_trim(e2sg(gas))) // 'time_lag,' &
                // e2sg(gas)(1:len_trim(e2sg(gas))) // 'def_timelag', separator)
            !> Trimmed, unlike the fourth slot's arm this replaces: those
            !> labels are character(32) and were concatenated unpadded,
            !> so the units row carried trailing blanks inside three of
            !> its fields. Every gas now spells them the same way.
            if(OutVarPresent(gas)) call AddDatum(header3, &
                trim(gas_dens_label(gas)) // ',' // trim(gas_conc_label(gas)) &
                // ',' // trim(gas_mixr_label(gas)) // ',[s],[1=default]', separator)
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
                      &[m+1s-1],[m+1s-1],[m+1s-1],[m+1s-1],[m+1s-1],[m+1s-1],[m+1s-1],&
                      &[m+1s-1],[deg_from_north],[deg],[deg],[deg]', separator)

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
        if(OutVarPresent(PrimaryWaterOutSlot())) call AddDatum(header1, ',', separator)
        if(OutVarPresent(PrimaryWaterOutSlot())) call AddDatum(header2,'un_LE,LE_scf', separator)
        if(OutVarPresent(PrimaryWaterOutSlot())) call AddDatum(header3,'[W+1m-2],[#]', separator)
        !> Uncorrected gas fluxes (Level 0) and spectral correction factors
        do k = 1, n_fo_slots
            gas = fo_slots(k)
            if(.not. OutVarPresent(gas)) cycle
            call AddDatum(header1, ',', separator)
            call AddDatum(header2, 'un_' // e2sg(gas)(1:len_trim(e2sg(gas))) &
                // 'flux,' // e2sg(gas)(1:len_trim(e2sg(gas))) // 'scf', separator)
            call AddDatum(header3, trim(gas_flux_label(gas)) // ',[#]', separator)
        end do

        !> Vickers and Mahrt 97 hard and soft flags
        call AddDatum(header1,'statistical_flags,,,,,,,,,,,', separator)
        call AddDatum(header2,'spikes_hf,amplitude_resolution_hf,drop_out_hf,absolute_limits_hf,&
            &skewness_kurtosis_hf,skewness_kurtosis_sf,discontinuities_hf,discontinuities_sf,timelag_hf,&
            &timelag_sf,attack_angle_hf,non_steady_wind_hf', separator)
        !> Eight per-variable flags, then two gas-only time-lag flags,
        !> then two single-variable ones. Built from the shared legend so
        !> the units row names exactly the digits the data row carries -
        !> it was ten copies of a literal naming eight variables while
        !> the row emitted sixty-nine characters.
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
            var = fo_slots(k)
            if(OutVarPresent(var)) then
                call AddDatum(header1, '', separator)
                call AddDatum(header2, e2sg(var)(1:len_trim(e2sg(var))) // 'spikes' , separator)
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
        !>
        !> The word each column is headed with follows the firmware of the
        !> analyser that column names. Both branches read E2Col(co2), so a
        !> site with a 7200 and a 7500 of different vintage had both headed
        !> from whichever held slot five - and the row writer, now asking
        !> each analyser separately, would have disagreed with this header.
        if (Diag7200%present) then
            if(CompareSwVer(InstrSwVerFor('li7200'), SwVerFromString('6.0.0'))) then
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
            if(CompareSwVer(InstrSwVerFor('li7500'), SwVerFromString('6.0.0'))) then
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
            if(OutVarPresent(gas)) call AddDatum(header1, '', separator)
            if(OutVarPresent(gas)) call AddDatum(header2, e2sg(gas)(1:len_trim(e2sg(gas))) // 'var', separator)
            if(OutVarPresent(gas)) call AddDatum(header3, '--', separator)
        end do
        !> w/ts covariance
        call AddDatum(header1, 'covariances', separator)
        call AddDatum(header2,'w/ts_cov', separator)
        call AddDatum(header3,'[m+1K+1s-1]', separator)
        !> w-gases covariances
        do k = 1, n_fo_slots
            gas = fo_slots(k)
            if(OutVarPresent(gas)) call AddDatum(header1, '', separator)
            if(OutVarPresent(gas)) call AddDatum(header2, 'w/' // e2sg(gas)(1:len_trim(e2sg(gas))) // 'cov', separator)
            if(OutVarPresent(gas)) call AddDatum(header3, '--', separator)
        end do

        !> Mean values of user variables
        if (NumUserVar > 0) then
            call AddDatum(header1, 'custom_variables', separator)
            do var = 1, NumUserVar
                call AddDatum(header1, '', separator)
                call AddDatum(header2, user_header(var)(1:len_trim(user_header(var))), separator)
                call AddDatum(header3, user_unit(var)(1:len_trim(user_unit(var))), separator)
            end do
        end if

        !> Conditional Eddy Covariance outputs (Zahn et al. 2022)
        if (EddyFlowProj%do_cec == 1) then
            call AddDatum(header1, 'conditional_eddy_covariance,,,,,,,,,,', separator)
        else if (EddyFlowProj%do_cec == 2) then
            call AddDatum(header1, 'conditional_eddy_covariance,,,,,', separator)
        else if (EddyFlowProj%do_cec == 3) then
            call AddDatum(header1, 'conditional_eddy_covariance,,,,', separator)
        end if
        if (EddyFlowProj%do_cec == 1 .or. EddyFlowProj%do_cec == 2) then
            call AddDatum(header2, 'E_cec,Tr_cec,E_cec_ET,Tr_cec_ET,r_ET_cec,qc_cec_h2o', separator)
            call AddDatum(header3, &
                '[mmol+1m-2s-1],[mmol+1m-2s-1],[mm+1hour-1],[mm+1hour-1],[#],[#]', separator)
        end if
        if (EddyFlowProj%do_cec == 1 .or. EddyFlowProj%do_cec == 3) then
            call AddDatum(header2, 'Reco_cec,P_cec,NEE_cec,r_Fc_cec,qc_cec_co2', separator)
            call AddDatum(header3, '[umol+1m-2s-1],[umol+1m-2s-1],[umol+1m-2s-1],[#],[#]', separator)
        end if

        !> Write on output file
        write(uflx, '(a)') header1(1:len_trim(header1) - 1)
        write(uflx, '(a)') header2(1:len_trim(header2) - 1)
        write(uflx, '(a)') header3(1:len_trim(header3) - 1)

    end if

    !>==========================================================================
    !>==========================================================================
    !> Metadata output
    if (EddyFlowProj%out_md) then
        Test_Path = Dir%main_out(1:len_trim(Dir%main_out)) &
                  // EddyFlowProj%id(1:len_trim(EddyFlowProj%id)) &
                  // MetaData_FilePadding // Timestamp_FilePadding // CsvExt
        dot = index(Test_Path, CsvExt, .true.) - 1
        Metadata_Path = Test_Path(1:dot) // CsvTmpExt
        open(umd, file = Metadata_Path, iostat = open_status, encoding = 'utf-8')

        call Clearstr(dataline)
        call AddDatum(dataline,'filename,date,time,DOY,latitude,longitude,altitude,canopy_height,displacement_height,&
            &roughness_length,file_length,acquisition_frequency,&
            &master_sonic_manufacturer,master_sonic_model,master_sonic_height,&
            &master_sonic_wformat,master_sonic_wref,master_sonic_north_offset,&
            &master_sonic_hpath_length,master_sonic_vpath_length,master_sonic_tau', separator)
        !> One block per configured gas, named from e2sg - the same tags the
        !> full output uses, already filled from FullOutputGasTags above.
        !>
        !> This was four hand-written blocks: three literals and a fourth
        !> built by concatenation. WriteOutMetadata emits fourteen fields per
        !> gas, and the fourth block declared only nine names - `tube_flowrate`
        !> and the five after it were concatenated with no comma between them,
        !> so they arrived as one run-together name. Header and row have
        !> disagreed by five fields for any project with a fourth gas, and
        !> every column after that point was read under the wrong name.
        do gas = firstGas, lastGas
            if (.not. OutVarPresent(gas)) cycle
            call AddDatum(dataline, &
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
        write(umd, '(a)') dataline(1:len_trim(dataline) - 1)
    end if

    !>*********************************************************************************************
    !>*********************************************************************************************

    !> Details of stationarity and integral turbulence tests
    if(RPsetup%out_qc_details .and. Meth%qcflag /= 'none') then
        !> Open file
        Test_Path = Dir%main_out(1:len_trim(Dir%main_out)) &
                  // EddyFlowProj%id(1:len_trim(EddyFlowProj%id)) &
                  // QCdetails_FilePadding // Timestamp_FilePadding // CsvExt
        dot = index(Test_Path, CsvExt, .true.) - 1
        QCdetails_Path = Test_Path(1:dot) // CsvTmpExt
        open(uqc, file = QCdetails_Path, iostat = open_status, encoding = 'utf-8')

        call Clearstr(header1)
        call Clearstr(header2)
        call Clearstr(header3)
        call AddDatum(header1,'file_info,,,,stationarity test,,', separator)
        call AddDatum(header2,'filename,date,time,DOY,dev(u),dev(w),dev(ts)', separator)
        call AddDatum(header3,',[yyyy-mm-dd],[HH:MM],[ddd.ddd],[%],[%],[%]', separator)
        !> One column per configured gas. These were four blocks
        !> naming co2/h2o/ch4 and whatever sat in the fourth slot,
        !> so a gas past the fourth had no dev or flag column at
        !> all. WriteOutQCDetails had the same four, which is what
        !> kept the two in step; they move together or the file
        !> shifts.
        do gas = firstGas, lastGas
            if (.not. OutVarPresent(gas)) cycle
            call AddDatum(header1,'', separator)
            call AddDatum(header2,'dev(' // e2sg(gas)(1:len_trim(e2sg(gas)) - 1) // ')', separator)
            call AddDatum(header3,'[%]', separator)
        end do

        call AddDatum(header1,',', separator)
        call AddDatum(header2,'dev(w/u),dev(w/ts)', separator)
        call AddDatum(header3,'[%],[%]', separator)

        do gas = firstGas, lastGas
            if (.not. OutVarPresent(gas)) cycle
            call AddDatum(header1,'', separator)
            call AddDatum(header2,'dev(w/' // e2sg(gas)(1:len_trim(e2sg(gas)) - 1) // ')', separator)
            call AddDatum(header3,'[%]', separator)
        end do

        call AddDatum(header1,',', separator)
        call AddDatum(header2,'flag(w/u),flag(w/ts)', separator)
        call AddDatum(header3,'[#],[#]', separator)

        do gas = firstGas, lastGas
            if (.not. OutVarPresent(gas)) cycle
            call AddDatum(header1,'', separator)
            call AddDatum(header2,'flag(w/' // e2sg(gas)(1:len_trim(e2sg(gas)) - 1) // ')', separator)
            call AddDatum(header3,'[#]', separator)
        end do

        call AddDatum(header1,'well-developed_turbulence_test,,,,,', separator)
        call AddDatum(header2,'dev(u),dev(w),dev(ts),flag(u),flag(w),flag(ts)', separator)
        call AddDatum(header3,'[%],[%],[%],[#],[#],[#]', separator)

        !> Write on output file
        write(uqc, '(a)') header1(1:len_trim(header1) - 1)
        write(uqc, '(a)') header2(1:len_trim(header2) - 1)
        write(uqc, '(a)') header3(1:len_trim(header3) - 1)
    end if

    !>*********************************************************************************************
    !>*********************************************************************************************

    !> Statistics files Level 1
    !> Statistics files st1..st7 all carry the same columns, so the header is
    !> built once here and written to whichever units are opened below.
    !>
    !> It used to be seven byte-identical literals naming twelve variables,
    !> while WriteOutStats looped `u, pe`. Those agreed while E2NumVar was 14;
    !> at 102 the rows carried 528 fields against a header of 88 and the files
    !> could not be read at all. Both sides now walk StatsLayoutSlots.
    call StatsLayoutSlots(st_slots, n_st_slots)
    n_st_gas = 0
    do i = 1, n_st_slots
        j = st_slots(i)
        if (j >= firstGas .and. j <= lastGas) then
            n_st_gas = n_st_gas + 1
            st_gas(n_st_gas) = j
            st_name(j) = e2sg(j)(1:len_trim(e2sg(j)) - 1)
        else if (j == u) then
            st_name(j) = 'u'
        else if (j == v) then
            st_name(j) = 'v'
        else if (j == w) then
            st_name(j) = 'w'
        else if (j == ts) then
            st_name(j) = 'ts'
        else if (j == tc) then
            st_name(j) = 'tc'
        else if (j == pi) then
            st_name(j) = 'pc'
        else if (j == te) then
            st_name(j) = 'te'
        else
            st_name(j) = 'pe'
        end if
    end do

    call clearstr(stats_header)
    call AddDatum(stats_header, 'filename,date,time,DOY,used_records', separator)
    do i = 1, n_st_slots
        call AddDatum(stats_header, 'mean(' // trim(st_name(st_slots(i))) // ')', separator)
    end do
    call AddDatum(stats_header, 'WindDirection', separator)
    do i = 1, n_st_slots
        call AddDatum(stats_header, 'var(' // trim(st_name(st_slots(i))) // ')', separator)
    end do
    call AddDatum(stats_header, 'cov(u/v),cov(u/w),cov(u/ts)', separator)
    do i = 1, n_st_gas
        call AddDatum(stats_header, 'cov(u/' // trim(st_name(st_gas(i))) // ')', separator)
    end do
    call AddDatum(stats_header, 'cov(v/w),cov(v/ts)', separator)
    do i = 1, n_st_gas
        call AddDatum(stats_header, 'cov(v/' // trim(st_name(st_gas(i))) // ')', separator)
    end do
    call AddDatum(stats_header, 'cov(w/ts)', separator)
    do i = 1, n_st_gas
        call AddDatum(stats_header, 'cov(w/' // trim(st_name(st_gas(i))) // ')', separator)
    end do
    call AddDatum(stats_header, 'cov(w/tc),cov(w/pc),cov(w/te),cov(w/pe)', separator)
    do i = 1, n_st_slots
        call AddDatum(stats_header, 'st_dev(' // trim(st_name(st_slots(i))) // ')', separator)
    end do
    do i = 1, n_st_slots
        call AddDatum(stats_header, 'skw(' // trim(st_name(st_slots(i))) // ')', separator)
    end do
    do i = 1, n_st_slots
        call AddDatum(stats_header, 'kur(' // trim(st_name(st_slots(i))) // ')', separator)
    end do

    if (RPsetup%out_st(1)) then
        Test_Path = StatsDir(1:len_trim(StatsDir)) &
                  // EddyFlowProj%id(1:len_trim(EddyFlowProj%id)) &
                  // Stats1_FilePadding // Timestamp_FilePadding // CsvExt
        dot = index(Test_Path, CsvExt, .true.) - 1
        St1_Path = Test_Path(1:dot) // CsvTmpExt
        open(ust1, file = St1_Path, iostat = open_status, encoding = 'utf-8')

        write(ust1, '(a)') 'first_statistics:_on_raw_data'
        write(ust1, '(a)') stats_header(1:len_trim(stats_header) - 1)
    end if

    !> Statistics files Level 2
    if (RPsetup%out_st(2)) then
        Test_Path = StatsDir(1:len_trim(StatsDir)) &
                  // EddyFlowProj%id(1:len_trim(EddyFlowProj%id)) &
                  // Stats2_FilePadding // Timestamp_FilePadding // CsvExt
        dot = index(Test_Path, CsvExt, .true.) - 1
        St2_Path = Test_Path(1:dot) // CsvTmpExt
        open(ust2, file = St2_Path, iostat = open_status, encoding = 'utf-8')

        write(ust2, '(a)') 'second_statistics:_on_raw_data_after_after_despiking'
        write(ust2, '(a)') stats_header(1:len_trim(stats_header) - 1)
    end if

    !> Statistics files Level 3
    if (RPsetup%out_st(3)) then
        Test_Path = StatsDir(1:len_trim(StatsDir)) &
                  // EddyFlowProj%id(1:len_trim(EddyFlowProj%id)) &
                  // Stats3_FilePadding // Timestamp_FilePadding // CsvExt
        dot = index(Test_Path, CsvExt, .true.) - 1
        St3_Path = Test_Path(1:dot) // CsvTmpExt
        open(ust3, file = St3_Path, iostat = open_status, encoding = 'utf-8')

        write(ust3, '(a)') 'third_statistics:_on_raw_data_after_after_despiking_and_cross-wind_correction'
        write(ust3, '(a)') stats_header(1:len_trim(stats_header) - 1)
    end if

    !> Statistics files Level 4
    if (RPsetup%out_st(4)) then
        Test_Path = StatsDir(1:len_trim(StatsDir)) &
                  // EddyFlowProj%id(1:len_trim(EddyFlowProj%id)) &
                  // Stats4_FilePadding // Timestamp_FilePadding // CsvExt
        dot = index(Test_Path, CsvExt, .true.) - 1
        St4_Path = Test_Path(1:dot) // CsvTmpExt
        open(ust4, file = St4_Path, iostat = open_status, encoding = 'utf-8')

        write(ust4, '(a)') 'forth statistics:_on_raw_data_after_despiking_cross_wind_correction&
                            &_and_angle-of-attack_correction'
        write(ust4, '(a)') stats_header(1:len_trim(stats_header) - 1)
    end if

    !> Statistics files Level 5
    if (RPsetup%out_st(5)) then
        Test_Path = StatsDir(1:len_trim(StatsDir)) &
                  // EddyFlowProj%id(1:len_trim(EddyFlowProj%id)) &
                  // Stats5_FilePadding // Timestamp_FilePadding // CsvExt
        dot = index(Test_Path, CsvExt, .true.) - 1
        St5_Path = Test_Path(1:dot) // CsvTmpExt
        open(ust5, file = St5_Path, iostat = open_status, encoding = 'utf-8')

        write(ust5, '(a)') 'fifth_statistics:_on_raw_data_after_despiking_cross_wind_correction&
                            &_angle-of-attack_correction_and_tilt_correction'
        write(ust5, '(a)') stats_header(1:len_trim(stats_header) - 1)
    end if

    !> Statistics files Level 6
    if (RPsetup%out_st(6)) then
        Test_Path = StatsDir(1:len_trim(StatsDir)) &
                  // EddyFlowProj%id(1:len_trim(EddyFlowProj%id)) &
                  // Stats6_FilePadding // Timestamp_FilePadding // CsvExt
        dot = index(Test_Path, CsvExt, .true.) - 1
        St6_Path = Test_Path(1:dot) // CsvTmpExt
        open(ust6, file = St6_Path, iostat = open_status, encoding = 'utf-8')

        write(ust6, '(a)') 'sixth statistics:_on_raw_data_after_despiking_cross_wind_correction&
            &_angle-of-attack_correction_tilt_correction_and_time-lag_compensation'
        write(ust6, '(a)') stats_header(1:len_trim(stats_header) - 1)
    end if

    !> Statistics files Level 7
    if (RPsetup%out_st(7)) then
        Test_Path = StatsDir(1:len_trim(StatsDir)) &
                  // EddyFlowProj%id(1:len_trim(EddyFlowProj%id)) &
                  // Stats7_FilePadding // Timestamp_FilePadding // CsvExt
        dot = index(Test_Path, CsvExt, .true.) - 1
        St7_Path = Test_Path(1:dot) // CsvTmpExt
        open(ust7, file = St7_Path, iostat = open_status, encoding = 'utf-8')

        write(ust7, '(a)') 'seventh_statistics:seventh_statistics:_on_raw_data_after_despiking_cross_wind_correction&
            &_angle-of-attack_correction_tilt_correction_time-lag_compensation_and_detrending'
        write(ust7, '(a)') stats_header(1:len_trim(stats_header) - 1)
    end if

contains

function FullOutputCustomLabel(ordinal) result(clean_label)
    integer, intent(in) :: ordinal
    character(64) :: clean_label

    character(32) :: model_token
    character(32) :: var_token
    character(64) :: label_token
    character(16) :: ordinal_label

    call clearstr(clean_label)
    var_token = SanitizeOutputToken(UserCol(ordinal)%var)
    model_token = CustomModelToken(UserCol(ordinal)%instr%model, var_token)
    var_token = CustomVarToken(var_token)

    select case (trim(var_token))
        case ('flowrate', 'co2', 'h2o', 'ch4', 'n2o', 'int_t', 'int_p')
            if (LabelHasAlpha(model_token)) then
                clean_label = trim(var_token) // '_' // trim(model_token)
            else
                clean_label = trim(var_token)
            end if
        case default
            label_token = UserCol(ordinal)%label
            call lowercase(label_token)
            label_token = replace2(label_token, 'custom_', '')
            label_token = replace2(label_token, '_mean', '')
            clean_label = SanitizeOutputToken(label_token)
            if (.not. LabelHasAlpha(clean_label)) then
                write(ordinal_label, '(i0)') ordinal
                clean_label = 'custom_' // trim(adjustl(ordinal_label))
            end if
    end select

    if (len_trim(clean_label) <= len(clean_label) - 5) &
        clean_label = trim(clean_label) // '_mean'
end function FullOutputCustomLabel

function CustomVarToken(raw_var_token) result(var_token)
    character(*), intent(in) :: raw_var_token
    character(32) :: var_token

    var_token = raw_var_token
    select case (trim(var_token))
        case ('cell_t', 'int_t_1', 'int_t_2')
            var_token = 'int_t'
    end select
end function CustomVarToken

function CustomModelToken(raw_model, var_token) result(model_token)
    character(*), intent(in) :: raw_model
    character(*), intent(in) :: var_token
    character(32) :: model_token

    model_token = SanitizeOutputToken(raw_model)
    select case (trim(var_token))
        case ('co2', 'h2o', 'ch4', 'n2o')
            if (index(model_token, 'miro_') == 1 .and. len_trim(model_token) > 5) &
                model_token = model_token(6:len_trim(model_token))
    end select
end function CustomModelToken

function SanitizeOutputToken(raw_token) result(clean_token)
    character(*), intent(in) :: raw_token
    character(32) :: clean_token

    integer :: i
    integer :: out_pos
    character(32) :: tmp

    call clearstr(clean_token)
    tmp = raw_token
    call lowercase(tmp)

    out_pos = 0
    do i = 1, len_trim(tmp)
        select case (tmp(i:i))
            case ('a':'z', '0':'9', '_', '-')
                if (out_pos < len(clean_token)) then
                    out_pos = out_pos + 1
                    clean_token(out_pos:out_pos) = tmp(i:i)
                end if
            case default
                if (out_pos < len(clean_token)) then
                    out_pos = out_pos + 1
                    clean_token(out_pos:out_pos) = '_'
                end if
        end select
    end do

    do while (index(clean_token, '__') > 0)
        clean_token = replace2(clean_token, '__', '_')
    end do
    do while (len_trim(clean_token) > 0 .and. clean_token(1:1) == '_')
        clean_token = clean_token(2:len_trim(clean_token))
    end do
    do while (len_trim(clean_token) > 0 &
        .and. clean_token(len_trim(clean_token):len_trim(clean_token)) == '_')
        clean_token(len_trim(clean_token):len_trim(clean_token)) = ' '
    end do
end function SanitizeOutputToken

logical function LabelHasAlpha(label)
    character(*), intent(in) :: label

    LabelHasAlpha = index(label, 'a') > 0 .or. index(label, 'b') > 0 &
        .or. index(label, 'c') > 0 .or. index(label, 'd') > 0 &
        .or. index(label, 'e') > 0 .or. index(label, 'f') > 0 &
        .or. index(label, 'g') > 0 .or. index(label, 'h') > 0 &
        .or. index(label, 'i') > 0 .or. index(label, 'j') > 0 &
        .or. index(label, 'k') > 0 .or. index(label, 'l') > 0 &
        .or. index(label, 'm') > 0 .or. index(label, 'n') > 0 &
        .or. index(label, 'o') > 0 .or. index(label, 'p') > 0 &
        .or. index(label, 'q') > 0 .or. index(label, 'r') > 0 &
        .or. index(label, 's') > 0 .or. index(label, 't') > 0 &
        .or. index(label, 'u') > 0 .or. index(label, 'v') > 0 &
        .or. index(label, 'w') > 0 .or. index(label, 'x') > 0 &
        .or. index(label, 'y') > 0 .or. index(label, 'z') > 0
end function LabelHasAlpha
end subroutine InitOutFiles_rp
