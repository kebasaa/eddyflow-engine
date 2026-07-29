!***************************************************************************
! init_outfiles_rp.f90
! --------------------
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
! \brief       Initializes Fluxnet output file
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine InitFluxnetFile_rp()
    use m_rp_global_var
    use iso_fortran_env
    implicit none
    !> in/out variables
    integer, external :: CreateDir
    !> local variables
    integer :: open_status = 1      ! initializing to false
    integer :: dot
    integer :: i
    integer :: j
    character(PathLen) :: Test_Path
    character(32) :: g4label
    character(64) :: e2sg(E2NumVar)
    character(64) :: usg(NumUserVar)
    character(LongOutstringLen) :: csv_row
    include '../src_common/interfaces.inc'


    !> Convenient strings
    e2sg(u)   = 'u_'
    e2sg(v)   = 'v_'
    e2sg(w)   = 'w_'
    e2sg(ts)  = 'ts_'
    e2sg(co2) = 'co2_'
    e2sg(h2o) = 'h2o_'
    e2sg(ch4) = 'ch4_'
    g4label = FourthGasLabel()
    e2sg(gas4) = g4label(1:len_trim(g4label)) // '_'
    e2sg(tc)  = 'cell_t_'
    e2sg(ti1) = 'inlet_t_'
    e2sg(ti2) = 'outlet_t_'
    e2sg(pi)  = 'cell_p_'
    e2sg(te)  = 'air_t_'
    e2sg(pe)  = 'air_p_'

    call lowercase(e2sg(gas4))
    
    do j = 1, NumUserVar
        usg(j) = SafeFluxnetCustomLabel(j)
        call lowercase(usg(j))
    end do

    Test_Path = Dir%main_out(1:len_trim(Dir%main_out)) &
                // EddyFlowProj%id(1:len_trim(EddyFlowProj%id)) &
                // FLUXNET_FilePadding // Timestamp_FilePadding // CsvExt
    dot = index(Test_Path, CsvExt, .true.) - 1
    FLUXNET_Path = Test_Path(1:dot) // CsvTmpExt
    open(uflxnt, file = FLUXNET_Path, iostat = open_status, encoding = 'utf-8')

    !> The gas set comes from the project, not from E2Col: this runs before
    !> the first DefineE2Set of the run, so E2Col is not final yet, and a set
    !> derived from it would not match what the row writers later emit.
    !> Records take slot firstGas+i-1 whether or not they name a column, which
    !> is how ApplyGasRecords assigns them; a record without a column occupies
    !> its slot but is absent, so it gets no column here.
    call SelectFluxnetGasSlots()

    call clearstr(csv_row)
    csv_row = 'TIMESTAMP_START,TIMESTAMP_END,DOY_START,DOY_END,FILENAME_HF,SW_IN_POT,NIGHT,EXPECT_NR,&
                &FILE_NR,CUSTOM_FILTER_NR,WD_FILTER_NR,SONIC_NR,T_SONIC_NR,CO2_NR,H2O_NR,CH4_NR,GS4_NR,&
                &TAU_NR,H_NR,FC_NR,LE_NR,FCH4_NR,FGS4_NR,&
                &TAU,H,LE,ET,FC,FH2O,FCH4,FGS4,TAU_RANDUNC_HF,H_RANDUNC_HF,LE_RANDUNC_HF,ET_RANDUNC_HF,&
                &FC_RANDUNC_HF,FH2O_RANDUNC_HF,FCH4_RANDUNC_HF,FGS4_RANDUNC_HF,&
                &SH_SINGLE,SLE_SINGLE,SET_SINGLE,SC_SINGLE,SH2O_SINGLE,SCH4_SINGLE,SGS4_SINGLE,&
                &FC_VADV,FH2O_VADV,FCH4_VADV,FGS4_VADV,&
                &U_UNROT,V_UNROT,W_UNROT,U,V,W,&
                &WS,WS_MAX,WD,WD_SIGMA,USTAR,TKE,MO_LENGTH,ZL,BOWEN,TSTAR,&
                &T_SONIC,TA_EP,PA_EP,RH_EP,AIR_MV,AIR_DENSITY,AIR_RHO_CP,AIR_CP,&
                &VAPOR_DENSITY,VAPOR_PARTIAL_PRESSURE,VAPOR_PARTIAL_PRESSURE_SAT,SPECIFIC_HUMIDITY,VPD_EP,TDEW,&
                &DRYAIR_PARTIAL_PRESSURE,DRYAIR_DENSITY,DRYAIR_MV,SPECIFIC_HEAT_EVAP,VAPOR_DRYAIR_RATIO,&
                &CO2_MEAS_TYPE,CO2_MOLAR_DENSITY,CO2_MIXING_RATIO,CO2,&
                &H2O_MEAS_TYPE,H2O_MOLAR_DENSITY,H2O_MIXING_RATIO,H2O,&
                &CH4_MEAS_TYPE,CH4_MOLAR_DENSITY,CH4_MIXING_RATIO,CH4,&
                &GS4_MEAS_TYPE,GS4_MOLAR_DENSITY,GS4_MIXING_RATIO,GS4,&
                &CO2_TLAG_ACTUAL,CO2_TLAG_USED,CO2_TLAG_NOMINAL,CO2_TLAG_MIN,CO2_TLAG_MAX,&
                &H2O_TLAG_ACTUAL,H2O_TLAG_USED,H2O_TLAG_NOMINAL,H2O_TLAG_MIN,H2O_TLAG_MAX,&
                &CH4_TLAG_ACTUAL,CH4_TLAG_USED,CH4_TLAG_NOMINAL,CH4_TLAG_MIN,CH4_TLAG_MAX,&
                &GS4_TLAG_ACTUAL,GS4_TLAG_USED,GS4_TLAG_NOMINAL,GS4_TLAG_MIN,GS4_TLAG_MAX,&
                &CO2_TLAG_PWB_SOURCE,H2O_TLAG_PWB_SOURCE,CH4_TLAG_PWB_SOURCE,GS4_TLAG_PWB_SOURCE,&
                &U_MEDIAN,V_MEDIAN,W_MEDIAN,T_SONIC_MEDIAN,&
                &CO2_MEAS_MEDIAN,H2O_MEAS_MEDIAN,CH4_MEAS_MEDIAN,GS4_MEAS_MEDIAN,&
                &U_P25,V_P25,W_P25,T_SONIC_P25,CO2_MEAS_P25,H2O_MEAS_P25,CH4_MEAS_P25,GS4_MEAS_P25,&
                &U_P75,V_P75,W_P75,T_SONIC_P75,CO2_MEAS_P75,H2O_MEAS_P75,CH4_MEAS_P75,GS4_MEAS_P75,&
                &U_SIGMA,V_SIGMA,W_SIGMA,T_SONIC_SIGMA,CO2_MEAS_SIGMA,H2O_MEAS_SIGMA,CH4_MEAS_SIGMA,GS4_MEAS_SIGMA,&
                &U_SKW,V_SKW,W_SKW,T_SONIC_SKW,CO2_MEAS_SKW,H2O_MEAS_SKW,CH4_MEAS_SKW,GS4_MEAS_SKW,&
                &U_KUR,V_KUR,W_KUR,T_SONIC_KUR,CO2_MEAS_KUR,H2O_MEAS_KUR,CH4_MEAS_KUR,GS4_MEAS_KUR,&
                &W_U_COV,W_T_SONIC_COV,W_CO2_MEAS_COV,W_H2O_MEAS_COV,W_CH4_MEAS_COV,W_GS4_MEAS_COV,&
                &CO2_MEAS_H2O_MEAS_COV,CO2_MEAS_CH4_MEAS_COV,CO2_MEAS_GS4_MEAS_COV,&
                &H2O_MEAS_CH4_MEAS_COV,H2O_MEAS_GS4_MEAS_COV,CH4_MEAS_GS4_MEAS_COV,&
                &FETCH_MAX,FETCH_OFFSET,FETCH_10,FETCH_30,FETCH_50,FETCH_70,FETCH_80,FETCH_90,&
                &USTAR_UNCORR,MO_LENGTH_UNCORR,ZL_UNCORR,&
                &TAU_UNCORR,H_UNCORR,LE_UNCORR,ET_UNCORR,FC_UNCORR,FH2O_UNCORR,FCH4_UNCORR,FGS4_UNCORR,&
                &TAU_STAGE1,H_STAGE1,LE_STAGE1,ET_STAGE1,FC_STAGE1,FH2O_STAGE1,FCH4_STAGE1,FGS4_STAGE1,&
                &TAU_STAGE2,H_STAGE2,LE_STAGE2,ET_STAGE2,FC_STAGE2,FH2O_STAGE2,FCH4_STAGE2,FGS4_STAGE2,&
                &T_CELL,PA_CELL,MV_AIR_CELL_CO2,MV_AIR_CELL_H2O,MV_AIR_CELL_CH4,MV_AIR_CELL_GS4,&
                &FH2O_CELL_CO2,FH2O_CELL_CH4,FH2O_CELL_GS4,H_CELL_CO2,H_CELL_H2O,H_CELL_CH4,H_CELL_GS4,&
                &H_BU_BOT,H_BU_TOP,H_BU_SPAR,SPEC_CORR_LI7700_A,SPEC_CORR_LI7700_B,SPEC_CORR_LI7700_C,&
                &TAU_SCF,H_SCF,LE_SCF,ET_SCF,FC_SCF,FH2O_SCF,FCH4_SCF,FGS4_SCF,&
                &W_T_SONIC_COV_IBROM,W_T_SONIC_COV_IBROM_N1626,W_T_SONIC_COV_IBROM_N0614,&
                &W_T_SONIC_COV_IBROM_N0277,W_T_SONIC_COV_IBROM_N0133,&
                &W_T_SONIC_COV_IBROM_N0065,W_T_SONIC_COV_IBROM_N0032,&
                &W_T_SONIC_COV_IBROM_N0016,W_T_SONIC_COV_IBROM_N0008,W_T_SONIC_COV_IBROM_N0004,&
                &U_NUM_SPIKES,V_NUM_SPIKES,W_NUM_SPIKES,T_SONIC_NUM_SPIKES,&
                &CO2_NUM_SPIKES,H2O_NUM_SPIKES,CH4_NUM_SPIKES,GS4_NUM_SPIKES,&
                &'
    !> Records excluded by each screening test.
    !>
    !> Three per-gas runs - diagnostics, the spike test and the absolute
    !> limits test - each of which used to be a fixed CO2/H2O/CH4/GS4
    !> quadruple. FCC copies this whole chunk verbatim and echoes it back
    !> without parsing it, so only its width has to follow the gas count.
    call AddDatum(csv_row, 'CUSTOM_FILTER_NREX', separator)
    call AddDatum(csv_row, 'WD_FILTER_NREX', separator)
    call AddDatum(csv_row, 'SONIC_DIAG_NREX', separator)
    do j = 1, nFluxnetLayoutSlots
        call AddDatum(csv_row, trim(FluxnetLayoutTags(j)) // '_DIAG_NREX', separator)
    end do
    call AddDatum(csv_row, 'U_SPIKE_NREX', separator)
    call AddDatum(csv_row, 'V_SPIKE_NREX', separator)
    call AddDatum(csv_row, 'W_SPIKE_NREX', separator)
    call AddDatum(csv_row, 'T_SONIC_SPIKE_NREX', separator)
    do j = 1, nFluxnetLayoutSlots
        call AddDatum(csv_row, trim(FluxnetLayoutTags(j)) // '_SPIKE_NREX', separator)
    end do
    call AddDatum(csv_row, 'U_ABSLIM_NREX', separator)
    call AddDatum(csv_row, 'V_ABSLIM_NREX', separator)
    call AddDatum(csv_row, 'W_ABSLIM_NREX', separator)
    call AddDatum(csv_row, 'T_SONIC_ABSLIM_NREX', separator)
    do j = 1, nFluxnetLayoutSlots
        call AddDatum(csv_row, trim(FluxnetLayoutTags(j)) // '_ABSLIM_NREX', separator)
    end do

    csv_row = trim(csv_row) // 'U_VM97_TEST,V_VM97_TEST,W_VM97_TEST,T_SONIC_VM97_TEST,&
                &CO2_VM97_TEST,H2O_VM97_TEST,CH4_VM97_TEST,GS4_VM97_TEST,&
                &VM97_TLAG_HF,VM97_TLAG_SF,VM97_AOA_HF,VM97_NSHW_HF,&
                &U_LGD,V_LGD,W_LGD,T_SONIC_LGD,CO2_LGD,H2O_LGD,CH4_LGD,GS4_LGD,&
                &U_KID,V_KID,W_KID,T_SONIC_KID,CO2_KID,H2O_KID,CH4_KID,GS4_KID,&
                &U_ZCD,V_ZCD,W_ZCD,T_SONIC_ZCD,CO2_ZCD,H2O_ZCD,CH4_ZCD,GS4_ZCD,&
                &TAU_CORRDIFF,H_CORRDIFF,LE_CORRDIFF,ET_CORRDIFF,FC_CORRDIFF,&
                &FH2O_CORRDIFF,FCH4_CORRDIFF,FGS4_CORRDIFF,&
                &TAU_NSR,H_NSR,FC_NSR,FH2O_NSR,FCH4_NSR,FGS4_NSR,&
                &TAU_SS,H_SS,FC_SS,FH2O_SS,FCH4_SS,FGS4_SS,&
                &U_ITC,W_ITC,T_SONIC_ITC,TAU_SS_TEST,H_SS_TEST,FC_SS_TEST,&
                &FH2O_SS_TEST,FCH4_SS_TEST,FGS4_SS_TEST,&
                &U_ITC_TEST,W_ITC_TEST,T_SONIC_ITC_TEST,&
                &TAU_SSITC_TEST,H_SSITC_TEST,LE_SSITC_TEST,ET_SSITC_TEST,FC_SSITC_TEST,&
                &FH2O_SSITC_TEST,FCH4_SSITC_TEST,FGS4_SSITC_TEST,&
                &INST_LI7200_HEAD_DETECT,INST_LI7200_T_OUT,INST_LI7200_T_IN,INST_LI7200_AUX_IN,&
                &INST_LI7200_DELTA_P,INST_LI7200_CHOPPER,INST_LI7200_DETECTOR,INST_LI7200_PLL,INST_LI7200_SYNC,&
                &INST_LI7500_CHOPPER,INST_LI7500_DETECTOR,INST_LI7500_PLL,INST_LI7500_SYNC,&
                &INST_LI7700_NOT_READY,INST_LI7700_NO_SIGNAL,INST_LI7700_RE_UNLOCKED,&
                &INST_LI7700_BAD_TEMP,INST_LI7700_LASER_T_UNREG,INST_LI7700_BLOCK_T_UNREG,&
                &INST_LI7700_MOTOR_SPINNING,INST_LI7700_PUMP_ON,INST_LI7700_TOP_HEATER_ON,&
                &INST_LI7700_BOTTOM_HEATER_ON,INST_LI7700_CALIBRATING,INST_LI7700_MOTOR_FAILURE,&
                &INST_LI7700_BAD_AUX_TC1,INST_LI7700_BAD_AUX_TC2,INST_LI7700_BAD_AUX_TC3,INST_LI7700_BOX_CONNECTED,&
                &INST_LI7200_AGC_OR_RSSI,INST_LI7500_AGC_OR_RSSI,INST_LI7700_RSSI,&
                &WBOOST_APPLIED,AOA_METHOD,&
                &AXES_ROTATION_METHOD,ROT_YAW,ROT_PITCH,ROT_ROLL,&
                &DETRENDING_METHOD,DENTRENDING_TIME_CONSTANT,&
                &TIMELAG_DETECTION_METHOD,WPL_APPLIED,BURBA_METHOD,&
                &SPECTRAL_CORRECTION_METHOD,FOOTPRINT_MODEL,&
                &LOGGER_SWVER_MAJOR,LOGGER_SWVER_MINOR,LOGGER_SWVER_REVISION,&
                &BADM_LOCATION_LAT,BADM_LOCATION_LONG,BADM_LOCATION_ELEV,BADM_HEIGHTC,&
                &DISPLACEMENT_HEIGHT,ROUGHNESS_LENGTH,&
                &FILE_TIME_DURATION,BADM_INST_SAMPLING_INT,BADM_INST_AVERAGING_INT,&
                &MANUFACTURER_SA,BADM_INST_MODEL_SA,BADM_INST_HEIGHT_SA,&
                &BADM_INST_SA_WIND_FORMAT,BADM_INST_SA_GILL_ALIGN,BADM_SA_OFFSET_NORTH,&
                &HPATH_SA,VPATH_SA,RESPONSE_TIME_SA,&
                &'
    !> Analyser of every configured gas, one block each, in slot order.
    !>
    !> These used to be four hard-coded blocks named CO2/H2O/CH4/GS4, so a
    !> fifth gas had no analyser columns here at all and had to carry them in
    !> a separate self-describing block further down. Generated per gas, the
    !> two are the same thing and the separate block is redundant.
    do j = 1, nFluxnetLayoutSlots
        g4label = FluxnetLayoutTags(j)
        call AddDatum(csv_row, 'MANUFACTURER_GA_' // trim(g4label), separator)
        call AddDatum(csv_row, 'BADM_INST_MODEL_GA_' // trim(g4label), separator)
        call AddDatum(csv_row, 'BADM_INSTPAIR_NORTHWARD_SEP_GA_' // trim(g4label), separator)
        call AddDatum(csv_row, 'BADM_INSTPAIR_EASTWARD_SEP_GA_' // trim(g4label), separator)
        call AddDatum(csv_row, 'BADM_INSTPAIR_HEIGHT_SEP_GA_' // trim(g4label), separator)
        call AddDatum(csv_row, 'BADM_INST_GA_CP_TUBE_LENGTH_GA_' // trim(g4label), separator)
        call AddDatum(csv_row, 'BADM_INST_GA_CP_TUBE_IN_DIAM_GA_' // trim(g4label), separator)
        call AddDatum(csv_row, 'BADM_INST_GA_CP_TUBE_FLOW_RATE_GA_' // trim(g4label), separator)
        !> The krypton coefficients belong to a hygrometer, so only the water
        !> slot carries them - as it always has.
        if (FluxnetLayoutSlots(j) == h2o) then
            call AddDatum(csv_row, 'KRYPTON_HYDRO_KH2O_GA_' // trim(g4label), separator)
            call AddDatum(csv_row, 'KRYPTON_HYDRO_KO2_GA_' // trim(g4label), separator)
        end if
        call AddDatum(csv_row, 'HPATH_GA_' // trim(g4label), separator)
        call AddDatum(csv_row, 'VPATH_GA_' // trim(g4label), separator)
        call AddDatum(csv_row, 'RESPONSE_TIME_GA_' // trim(g4label), separator)
    end do

            !> If need to reitroduce details of VM, paste this after line:
            !>  "&CO2_ABSLIM_NREX,H2O_ABSLIM_NREX,CH4_ABSLIM_NREX,GS4_ABSLIM_NREX,&"

            !   &U_SPIKE_VM97_NR,V_SPIKE_VM97_NR,W_SPIKE_VM97_NR,T_SONIC_SPIKE_VM97_NR,&
            !   &CO2_SPIKE_VM97_NR,H2O_SPIKE_VM97_NR,CH4_SPIKE_VM97_NR,GS4_SPIKE_VM97_NR,&
            !   &U_AMPRES_VM97,V_AMPRES_VM97,W_AMPRES_VM97,T_SONIC_AMPRES_VM97,&
            !   &CO2_AMPRES_VM97,H2O_AMPRES_VM97,CH4_AMPRES_VM97,GS4_AMPRES_VM97,&
            !   &U_DRPOUT_C_VM97_NR,V_DRPOUT_C_VM97_NR,W_DRPOUT_C_VM97_NR,T_SONIC_DRPOUT_C_VM97_NR,&
            !   &CO2_DRPOUT_C_VM97_NR,H2O_DRPOUT_C_VM97_NR,CH4_DRPOUT_C_VM97_NR,GS4_DRPOUT_C_VM97_NR,&
            !   &U_DRPOUT_X_VM97_NR,V_DRPOUT_X_VM97_NR,W_DRPOUT_X_VM97_NR,T_SONIC_DRPOUT_X_VM97_NR,&
            !   &CO2_DRPOUT_X_VM97_NR,H2O_DRPOUT_X_VM97_NR,CH4_DRPOUT_X_VM97_NR,GS4_DRPOUT_X_VM97_NR,&
            !   &U_SKW_VM97_NR,V_SKW_VM97_NR,W_SKW_VM97_NR,T_SONIC_SKW_VM97_NR,&
            !   &CO2_SKW_VM97_NR,H2O_SKW_VM97_NR,CH4_SKW_VM97_NR,GS4_SKW_VM97_NR,&
            !   &U_KUR_VM97_NR,V_KUR_VM97_NR,W_KUR_VM97_NR,T_SONIC_KUR_VM97_NR,&
            !   &CO2_KUR_VM97_NR,H2O_KUR_VM97_NR,CH4_KUR_VM97_NR,GS4_KUR_VM97_NR,&
            !   &AOA_VM97_NR,WS_SS_ALONG_VM97,WS_SS_CROSS_VM97,WS_SS_VM97,&

            !> If need to reitroduce VM flags for last 3 tests, paste this after:
            !> "&VM97_DISCON_HFLAG,VM97_DISCON_SFLAG,&"

            !   &VM97_TIMELAG_HFLAG,VM97_TIMELAG_SFLAG,VM97_AOA_HFLAG,VM97_NSW_HFLAG,&


    !> Width of the fixed part of the row, counted from what was just built.
    !> Every datum above was appended with a trailing separator, so the field
    !> count is the number of separators. WriteOutFluxnetOnlyBiomet pads to
    !> this rather than to a literal of its own.
    nFluxnetFixedCols = 0
    do i = 1, len_trim(csv_row)
        if (csv_row(i:i) == separator) nFluxnetFixedCols = nFluxnetFixedCols + 1
    end do

    !> Add custom variables
    call AddDatum(csv_row, 'NUM_CUSTOM_VARS', separator)
    if (NumUserVar > 0) then
        do i = 1, NumUserVar
            call uppercase(usg(i)) 
            call AddDatum(csv_row, 'CUSTOM_' // usg(i)(1:len_trim(usg(i))) &
                // '_MEAN', separator)
        end do
    end if

    !> Per-gas water vapour terms.
    !>
    !> The fixed part of the row carries four gas slots. The H2O used to
    !> correct each gas is now per-gas though - with two analysers, CO2 from
    !> one is corrected with H2O from that same one - so the resolved terms are
    !> written here, one pair per configured gas. FCC recomputes the fluxes
    !> from this file and reads them back rather than repeating the resolution.
    !>
    !> The count is always emitted so the block is self-describing and the
    !> reader needs no knowledge of the project configuration.
    call AddDatum(csv_row, 'NUM_GAS_MOIST', separator)
    do j = 1, nFluxnetGasSlots
        i = FluxnetGasSlots(j)
        call AddDatum(csv_row, trim(FluxnetGasTags(j)) // '_MOIST_SLOT', separator)
        call AddDatum(csv_row, trim(FluxnetGasTags(j)) // '_MOIST_RHOW', separator)
        call AddDatum(csv_row, trim(FluxnetGasTags(j)) // '_MOIST_SIGMA', separator)
    end do

    !> Analyser describing each gas beyond the four historical slots.
    !>
    !> The MANUFACTURER_GA_* blocks above cover CO2/H2O/CH4/GS4 only, and FCC
    !> reaches them by instrument role - past the fourth gas that role index
    !> addresses an unrelated instrument. These columns carry the same
    !> metadata for the remaining gases so FCC can correct them too.
    !>
    !> Unlike the GA_* blocks, these are written in SI units and read back
    !> unchanged: the legacy columns are in metadata units and are converted
    !> on read, which is a double-conversion waiting to happen.
    call AddDatum(csv_row, 'NUM_GAS_INSTR', separator)
    do j = 1, nFluxnetInstrSlots
        i = FluxnetInstrSlots(j)
        call AddDatum(csv_row, trim(FluxnetInstrTags(j)) // '_INSTR_SLOT', separator)
        call AddDatum(csv_row, trim(FluxnetInstrTags(j)) // '_INSTR_MANUFACTURER', separator)
        call AddDatum(csv_row, trim(FluxnetInstrTags(j)) // '_INSTR_MODEL', separator)
        call AddDatum(csv_row, trim(FluxnetInstrTags(j)) // '_INSTR_NSEP', separator)
        call AddDatum(csv_row, trim(FluxnetInstrTags(j)) // '_INSTR_ESEP', separator)
        call AddDatum(csv_row, trim(FluxnetInstrTags(j)) // '_INSTR_VSEP', separator)
        call AddDatum(csv_row, trim(FluxnetInstrTags(j)) // '_INSTR_TUBE_LENGTH', separator)
        call AddDatum(csv_row, trim(FluxnetInstrTags(j)) // '_INSTR_TUBE_IN_DIAM', separator)
        call AddDatum(csv_row, trim(FluxnetInstrTags(j)) // '_INSTR_TUBE_FLOW_RATE', separator)
        call AddDatum(csv_row, trim(FluxnetInstrTags(j)) // '_INSTR_HPATH', separator)
        call AddDatum(csv_row, trim(FluxnetInstrTags(j)) // '_INSTR_VPATH', separator)
        call AddDatum(csv_row, trim(FluxnetInstrTags(j)) // '_INSTR_RESPONSE_TIME', separator)
        call AddDatum(csv_row, trim(FluxnetInstrTags(j)) // '_INSTR_KH2O', separator)
        call AddDatum(csv_row, trim(FluxnetInstrTags(j)) // '_INSTR_KO2', separator)
    end do

    !> Per-gas families for each gas past the four historical slots.
    !>
    !> Omits the NREX counts and the VM97 flag string: those reach FCC only
    !> inside the raw chunks it echoes verbatim, not as per-slot values, so it
    !> could not reproduce them here. They stay four-gas until those chunks are
    !> parsed per slot.
    !>
    !> The fixed part of the row carries these families for CO2/H2O/CH4/GS4
    !> only, interleaved and in fixed positions. Rather than renumber ~300
    !> columns, the same families are emitted here for the remaining gases,
    !> grouped per gas. Consumers read this file by column name, so grouping
    !> rather than interleaving costs nothing.
    call AddDatum(csv_row, 'NUM_GAS_EXTRA', separator)
    do j = 1, nFluxnetInstrSlots
        call AddDatum(csv_row, trim(FluxnetInstrTags(j)) // '_SLOT', separator)
        call AddDatum(csv_row, trim(FluxnetInstrTags(j)) // '_NR', separator)
        call AddDatum(csv_row, trim(FluxnetInstrTags(j)) // '_NR_W', separator)
        call AddDatum(csv_row, trim(FluxnetInstrTags(j)) // '_MEAS_TYPE', separator)
        call AddDatum(csv_row, trim(FluxnetInstrTags(j)) // '_MOLAR_DENSITY', separator)
        call AddDatum(csv_row, trim(FluxnetInstrTags(j)) // '_MIXING_RATIO', separator)
        call AddDatum(csv_row, trim(FluxnetInstrTags(j)) // '_MEAS', separator)
        call AddDatum(csv_row, trim(FluxnetInstrTags(j)) // '_FLUX_LEVEL0', separator)
        call AddDatum(csv_row, trim(FluxnetInstrTags(j)) // '_FLUX', separator)
        call AddDatum(csv_row, trim(FluxnetInstrTags(j)) // '_FLUX_STAGE1', separator)
        call AddDatum(csv_row, trim(FluxnetInstrTags(j)) // '_FLUX_STAGE2', separator)
        call AddDatum(csv_row, trim(FluxnetInstrTags(j)) // '_SCF', separator)
        call AddDatum(csv_row, trim(FluxnetInstrTags(j)) // '_RANDUNC_HF', separator)
        call AddDatum(csv_row, trim(FluxnetInstrTags(j)) // '_STORAGE', separator)
        call AddDatum(csv_row, trim(FluxnetInstrTags(j)) // '_TLAG_ACTUAL', separator)
        call AddDatum(csv_row, trim(FluxnetInstrTags(j)) // '_TLAG_USED', separator)
        call AddDatum(csv_row, trim(FluxnetInstrTags(j)) // '_TLAG_NOMINAL', separator)
        call AddDatum(csv_row, trim(FluxnetInstrTags(j)) // '_TLAG_MIN', separator)
        call AddDatum(csv_row, trim(FluxnetInstrTags(j)) // '_TLAG_MAX', separator)
        call AddDatum(csv_row, trim(FluxnetInstrTags(j)) // '_TLAG_PWB_SOURCE', separator)
        call AddDatum(csv_row, trim(FluxnetInstrTags(j)) // '_MEAS_MEDIAN', separator)
        call AddDatum(csv_row, trim(FluxnetInstrTags(j)) // '_MEAS_P25', separator)
        call AddDatum(csv_row, trim(FluxnetInstrTags(j)) // '_MEAS_P75', separator)
        call AddDatum(csv_row, trim(FluxnetInstrTags(j)) // '_MEAS_SIGMA', separator)
        call AddDatum(csv_row, trim(FluxnetInstrTags(j)) // '_MEAS_SKW', separator)
        call AddDatum(csv_row, trim(FluxnetInstrTags(j)) // '_MEAS_KUR', separator)
        call AddDatum(csv_row, trim(FluxnetInstrTags(j)) // '_W_MEAS_COV', separator)
        call AddDatum(csv_row, trim(FluxnetInstrTags(j)) // '_MV_AIR_CELL', separator)
        call AddDatum(csv_row, trim(FluxnetInstrTags(j)) // '_NUM_SPIKES', separator)
    end do

    !> Add biomet variables
    call AddDatum(csv_row, 'NUM_BIOMET_VARS', separator)

    if (nbVars > 0) then
        do i = 1, nbVars
            if (EddyFlowProj%fluxnet_standardize_biomet) then
                call AddDatum(csv_row, trim(bVars(i)%fluxnet_label), separator)
            else
                call AddDatum(csv_row, trim(bVars(i)%label), separator)
            end if
        end do
    end if

    call uppercase(e2sg(gas4))
    csv_row = replace2(csv_row, 'GS4', e2sg(gas4)(1:len_trim(e2sg(gas4)) - 1))

    !> CEC partitioning ratios (always present; error when do_cec=0)
    call AddDatum(csv_row, 'r_ET_cec', separator)
    call AddDatum(csv_row, 'r_Fc_cec', separator)
    call AddDatum(csv_row, 'CEC_N_VALID', separator)
    call AddDatum(csv_row, 'CEC_N_O1', separator)
    call AddDatum(csv_row, 'CEC_N_O2', separator)
    call AddDatum(csv_row, 'CEC_FRAC_O1', separator)
    call AddDatum(csv_row, 'CEC_FRAC_O2', separator)
    call AddDatum(csv_row, 'CEC_H2O_VALID', separator)
    call AddDatum(csv_row, 'CEC_CO2_VALID', separator)
    call AddDatum(csv_row, 'CEC_H2O_STATUS', separator)
    call AddDatum(csv_row, 'CEC_CO2_STATUS', separator)

    write(uflxnt, '(a)') csv_row(1:len_trim(csv_row) - 1)

contains

function SafeFluxnetCustomLabel(ordinal) result(clean_label)
    integer, intent(in) :: ordinal
    character(64) :: clean_label

    character(64) :: tmp
    character(32) :: model_token
    character(32) :: var_token
    character(16) :: ordinal_label

    call clearstr(clean_label)
    var_token = SanitizeFluxnetToken(UserCol(ordinal)%var)
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
            tmp = UserCol(ordinal)%label
            call lowercase(tmp)
            tmp = replace2(tmp, 'custom_', '')
            tmp = replace2(tmp, '_mean', '')
            clean_label = SanitizeFluxnetToken(tmp)
            if (.not. LabelHasAlpha(clean_label)) then
                write(ordinal_label, '(i0)') ordinal
                clean_label = 'custom_' // trim(adjustl(ordinal_label))
            end if
    end select

    if (.not. LabelHasAlpha(clean_label)) then
        write(ordinal_label, '(i0)') ordinal
        if (len_trim(var_token) > 0) then
            clean_label = trim(var_token) // '_' // trim(adjustl(ordinal_label))
        else
            clean_label = 'custom_' // trim(adjustl(ordinal_label))
        end if
    end if
end function SafeFluxnetCustomLabel

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

    model_token = SanitizeFluxnetToken(raw_model)
    select case (trim(var_token))
        case ('co2', 'h2o', 'ch4', 'n2o')
            if (index(model_token, 'miro_') == 1 .and. len_trim(model_token) > 5) &
                model_token = model_token(6:len_trim(model_token))
    end select
end function CustomModelToken

!> Column-name tag for a gas slot.
!>
!> Slots 5-8 keep the names the FLUXNET row has always used, so existing files
!> and the consumers that read them are unaffected; GS4 in particular is the
!> historical name of the fourth slot whatever gas occupies it. Gases beyond
!> those four are named after the gas itself.
!> Gas slots that get columns, and which of them carry an analyser block.
!>
!> Derived from the project configuration so it is knowable before any data is
!> read. Mirrors the slot assignment ApplyGasRecords performs.
subroutine SelectFluxnetGasSlots()
    integer :: k
    integer :: slot
    integer :: dup
    integer :: n
    character(32) :: tag
    character(8) :: ord

    nFluxnetGasSlots = 0
    nFluxnetInstrSlots = 0
    nFluxnetLayoutSlots = 0

    !> The layout list first: every configured gas, whether or not it has a
    !> column. The fixed part of the row is sized from this, so a gas selected
    !> without data still gets its column set and the fields after it stay put.
    !> Names are assigned here, once, and the present-gas list borrows them -
    !> deriving them twice is how two lists of the same gases end up disagreeing.
    if (EddyFlowProj%gas_num > 0) then
        do k = 1, min(EddyFlowProj%gas_num, MaxNumGases)
            slot = firstGas + k - 1
            if (slot > lastGas) exit
            tag = EddyFlowProj%gas(k)%var
            call uppercase(tag)
            if (len_trim(tag) == 0) tag = 'GAS'
            dup = 0
            do n = 1, nFluxnetLayoutSlots
                if (trim(FluxnetLayoutTags(n)) == trim(tag) .or. &
                    index(trim(FluxnetLayoutTags(n)), trim(tag) // '_') == 1) &
                    dup = dup + 1
            end do
            if (dup > 0) then
                write(ord, '(i0)') dup + 1
                tag = trim(tag) // '_' // trim(ord)
            end if
            nFluxnetLayoutSlots = nFluxnetLayoutSlots + 1
            FluxnetLayoutSlots(nFluxnetLayoutSlots) = slot
            FluxnetLayoutTags(nFluxnetLayoutSlots) = tag
        end do
    end if

    if (EddyFlowProj%gas_num > 0) then
        do k = 1, min(EddyFlowProj%gas_num, MaxNumGases)
            slot = firstGas + k - 1
            if (slot > lastGas) exit
            if (EddyFlowProj%gas(k)%col <= 0) cycle
            nFluxnetGasSlots = nFluxnetGasSlots + 1
            FluxnetGasSlots(nFluxnetGasSlots) = slot
            !> Every slot is named for its species, the fourth included.
            !>
            !> The fourth used to be tagged 'GS4' here and have its real label
            !> substituted into the finished header much later. That left this
            !> loop blind to the name the slot would actually carry, so the
            !> duplicate check below compared later gases against the literal
            !> 'GS4' instead of against a species. A project measuring the same
            !> gas in slot four and slot five therefore emitted two identical
            !> sets of columns - 20 duplicate names in one row - with no _2
            !> suffix to tell them apart.
            !>
            !> Slots one to three are unaffected: their records say co2, h2o
            !> and ch4, which is what HistoricGasTag returned for them.
            !> Borrow the name the layout pass already assigned to this slot,
            !> including any _2 suffix it needed. Re-deriving it here would
            !> disambiguate against a different set of gases.
            call clearstr(tag)
            do n = 1, nFluxnetLayoutSlots
                if (FluxnetLayoutSlots(n) /= slot) cycle
                tag = FluxnetLayoutTags(n)
                exit
            end do
            if (len_trim(tag) == 0) tag = 'GAS'
            FluxnetGasTags(nFluxnetGasSlots) = tag
            !> Only gases past the four historical slots need an analyser
            !> block; CO2/H2O/CH4/GS4 already have their GA_* columns.
            if (slot > gas4) then
                nFluxnetInstrSlots = nFluxnetInstrSlots + 1
                FluxnetInstrSlots(nFluxnetInstrSlots) = slot
                FluxnetInstrTags(nFluxnetInstrSlots) = &
                    FluxnetGasTags(nFluxnetGasSlots)
            end if
        end do
    else
        !> Legacy projects name their gases by fixed slot.
        do slot = co2, gas4
            if (EddyFlowProj%Col(slot) <= 0) cycle
            nFluxnetGasSlots = nFluxnetGasSlots + 1
            FluxnetGasSlots(nFluxnetGasSlots) = slot
            FluxnetGasTags(nFluxnetGasSlots) = HistoricGasTag(slot)
        end do
    end if
end subroutine SelectFluxnetGasSlots

function HistoricGasTag(gas_slot) result(tag)
    integer, intent(in) :: gas_slot
    character(32) :: tag

    call clearstr(tag)
    select case (gas_slot)
        case (co2);  tag = 'CO2'
        case (h2o);  tag = 'H2O'
        case (ch4);  tag = 'CH4'
        case default; tag = 'GS4'
    end select
end function HistoricGasTag

function FluxnetGasTag(gas_slot) result(tag)
    integer, intent(in) :: gas_slot
    character(32) :: tag

    call clearstr(tag)
    select case (gas_slot)
        case (co2);  tag = 'CO2'
        case (h2o);  tag = 'H2O'
        case (ch4);  tag = 'CH4'
        case (gas4); tag = 'GS4'
        case default
            tag = SanitizeFluxnetToken(E2Col(gas_slot)%var)
            call uppercase(tag)
    end select
end function FluxnetGasTag

function SanitizeFluxnetToken(raw_token) result(clean_token)
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
end function SanitizeFluxnetToken

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

end subroutine InitFluxnetFile_rp
