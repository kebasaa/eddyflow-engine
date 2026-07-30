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
    integer :: k
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
                &FILE_NR,CUSTOM_FILTER_NR,WD_FILTER_NR,'

    !> Record counts: the sonic, the sonic temperature, then one per gas.
    call AddDatum(csv_row, 'SONIC_NR', separator)
    call AddDatum(csv_row, 'T_SONIC_NR', separator)
    do j = 1, nFluxnetLayoutSlots
        call AddDatum(csv_row, trim(FluxnetLayoutTags(j)) // '_NR', separator)
    end do
    !> Records behind each covariance with w. Water is LE in this family.
    call AddDatum(csv_row, 'TAU_NR', separator)
    call AddDatum(csv_row, 'H_NR', separator)
    do j = 1, nFluxnetLayoutSlots
        call AddDatum(csv_row, trim(FluxnetNrwTag(j)) // '_NR', separator)
    end do

    call AddFluxFamily('')
    call AddFluxFamily('_RANDUNC_HF')

    !> Single-point storage. CO2 is SC, not SCO2.
    call AddDatum(csv_row, 'SH_SINGLE', separator)
    call AddDatum(csv_row, 'SLE_SINGLE', separator)
    call AddDatum(csv_row, 'SET_SINGLE', separator)
    do j = 1, nFluxnetLayoutSlots
        call AddDatum(csv_row, trim(FluxnetStorTag(j)) // '_SINGLE', separator)
    end do
    !> Vertical advection, per gas only.
    do j = 1, nFluxnetLayoutSlots
        call AddDatum(csv_row, trim(FluxnetFluxTag(j)) // '_VADV', separator)
    end do

    csv_row = trim(csv_row) // 'U_UNROT,V_UNROT,W_UNROT,U,V,W,&
                &WS,WS_MAX,WD,WD_SIGMA,USTAR,TKE,MO_LENGTH,ZL,BOWEN,TSTAR,&
                &T_SONIC,TA_EP,PA_EP,RH_EP,AIR_MV,AIR_DENSITY,AIR_RHO_CP,AIR_CP,&
                &VAPOR_DENSITY,VAPOR_PARTIAL_PRESSURE,VAPOR_PARTIAL_PRESSURE_SAT,SPECIFIC_HUMIDITY,VPD_EP,TDEW,&
                &DRYAIR_PARTIAL_PRESSURE,DRYAIR_DENSITY,DRYAIR_MV,SPECIFIC_HEAT_EVAP,VAPOR_DRYAIR_RATIO,&
                &'

    !> Concentration quadruple, then the timelag quintuplet, then the PWB
    !> lag source - each one block per configured gas.
    do j = 1, nFluxnetLayoutSlots
        g4label = FluxnetLayoutTags(j)
        call AddDatum(csv_row, trim(g4label) // '_MEAS_TYPE', separator)
        call AddDatum(csv_row, trim(g4label) // '_MOLAR_DENSITY', separator)
        call AddDatum(csv_row, trim(g4label) // '_MIXING_RATIO', separator)
        call AddDatum(csv_row, trim(g4label), separator)
    end do
    do j = 1, nFluxnetLayoutSlots
        g4label = FluxnetLayoutTags(j)
        call AddDatum(csv_row, trim(g4label) // '_TLAG_ACTUAL', separator)
        call AddDatum(csv_row, trim(g4label) // '_TLAG_USED', separator)
        call AddDatum(csv_row, trim(g4label) // '_TLAG_NOMINAL', separator)
        call AddDatum(csv_row, trim(g4label) // '_TLAG_MIN', separator)
        call AddDatum(csv_row, trim(g4label) // '_TLAG_MAX', separator)
    end do
    do j = 1, nFluxnetLayoutSlots
        call AddDatum(csv_row, trim(FluxnetLayoutTags(j)) // '_TLAG_PWB_SOURCE', separator)
    end do

    !> Per-variable statistics. The wind half and the gas half take different
    !> suffixes: U_MEDIAN beside CO2_MEAS_MEDIAN.
    call AddStatFamily('_MEDIAN', '_MEAS_MEDIAN')
    call AddStatFamily('_P25', '_MEAS_P25')
    call AddStatFamily('_P75', '_MEAS_P75')
    call AddStatFamily('_SIGMA', '_MEAS_SIGMA')
    call AddStatFamily('_SKW', '_MEAS_SKW')
    call AddStatFamily('_KUR', '_MEAS_KUR')

    !> Covariances with w, then the gas-gas triangle.
    call AddDatum(csv_row, 'W_U_COV', separator)
    call AddDatum(csv_row, 'W_T_SONIC_COV', separator)
    do j = 1, nFluxnetLayoutSlots
        call AddDatum(csv_row, 'W_' // trim(FluxnetLayoutTags(j)) // '_MEAS_COV', separator)
    end do
    do j = 1, nFluxnetLayoutSlots - 1
        do k = j + 1, nFluxnetLayoutSlots
            call AddDatum(csv_row, trim(FluxnetLayoutTags(j)) // '_MEAS_' &
                // trim(FluxnetLayoutTags(k)) // '_MEAS_COV', separator)
        end do
    end do

    csv_row = trim(csv_row) // 'FETCH_MAX,FETCH_OFFSET,FETCH_10,FETCH_30,FETCH_50,FETCH_70,FETCH_80,FETCH_90,&
                &USTAR_UNCORR,MO_LENGTH_UNCORR,ZL_UNCORR,&
                &'

    call AddFluxFamily('_UNCORR')
    call AddFluxFamily('_STAGE1')
    call AddFluxFamily('_STAGE2')

    !> Cell quantities: the shared cell T and P, the per-gas cell molar volume,
    !> the per-gas water flux in the cell - which the water slot itself does not
    !> have - and the per-gas cell sensible heat.
    call AddDatum(csv_row, 'T_CELL', separator)
    call AddDatum(csv_row, 'PA_CELL', separator)
    do j = 1, nFluxnetLayoutSlots
        call AddDatum(csv_row, 'MV_AIR_CELL_' // trim(FluxnetLayoutTags(j)), separator)
    end do
    do j = 1, nFluxnetLayoutSlots
        if (FluxnetLayoutSlots(j) == h2o) cycle
        call AddDatum(csv_row, 'FH2O_CELL_' // trim(FluxnetLayoutTags(j)), separator)
    end do
    do j = 1, nFluxnetLayoutSlots
        call AddDatum(csv_row, 'H_CELL_' // trim(FluxnetLayoutTags(j)), separator)
    end do

    csv_row = trim(csv_row) // 'H_BU_BOT,H_BU_TOP,H_BU_SPAR,&
                &SPEC_CORR_LI7700_A,SPEC_CORR_LI7700_B,SPEC_CORR_LI7700_C,&
                &'

    call AddFluxFamily('_SCF')

    csv_row = trim(csv_row) // 'W_T_SONIC_COV_IBROM,W_T_SONIC_COV_IBROM_N1626,W_T_SONIC_COV_IBROM_N0614,&
                &W_T_SONIC_COV_IBROM_N0277,W_T_SONIC_COV_IBROM_N0133,&
                &W_T_SONIC_COV_IBROM_N0065,W_T_SONIC_COV_IBROM_N0032,&
                &W_T_SONIC_COV_IBROM_N0016,W_T_SONIC_COV_IBROM_N0008,W_T_SONIC_COV_IBROM_N0004,&
                &'

    call AddStatFamily('_NUM_SPIKES', '_NUM_SPIKES')
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

    !> Vickers and Mahrt (1997) test outcomes, one field per variable: the
    !> wind components and sonic temperature, then one per configured gas.
    !>
    !> Unlike the NREX chunk above, FCC does not echo this one - it parses the
    !> fields into lEx%vm_flags and writes them out again - so the reader and
    !> FCC's re-emit had to move with this loop, not just its width.
    call AddDatum(csv_row, 'U_VM97_TEST', separator)
    call AddDatum(csv_row, 'V_VM97_TEST', separator)
    call AddDatum(csv_row, 'W_VM97_TEST', separator)
    call AddDatum(csv_row, 'T_SONIC_VM97_TEST', separator)
    do j = 1, nFluxnetLayoutSlots
        call AddDatum(csv_row, trim(FluxnetLayoutTags(j)) // '_VM97_TEST', separator)
    end do

    csv_row = trim(csv_row) // 'VM97_TLAG_HF,VM97_TLAG_SF,VM97_AOA_HF,VM97_NSHW_HF,'

    !> Longest gap duration, kurtosis index on differenced variables, and zero
    !> crossing distance. Variable-shaped, like the NREX counts: the wind
    !> components and sonic temperature, then one field per configured gas.
    !>
    !> FCC copies this whole chunk verbatim - it is fluxnetChunks%s(2) - so only
    !> the header, RP's writer and the chunk's width had to follow the gas count.
    call AddVariableFamily('_LGD')
    call AddVariableFamily('_KID')
    call AddVariableFamily('_ZCD')

    !> Correlation differences with and without repeated values. Flux-shaped
    !> rather than variable-shaped: the momentum and heat fluxes first, then one
    !> per configured gas.
    call AddDatum(csv_row, 'TAU_CORRDIFF', separator)
    call AddDatum(csv_row, 'H_CORRDIFF', separator)
    call AddDatum(csv_row, 'LE_CORRDIFF', separator)
    call AddDatum(csv_row, 'ET_CORRDIFF', separator)
    do j = 1, nFluxnetLayoutSlots
        call AddDatum(csv_row, trim(FluxnetFluxTag(j)) // '_CORRDIFF', separator)
    end do

    !> Mahrt 1998 nonstationarity ratios. Flux-shaped, without the two water
    !> fluxes: momentum and sensible heat, then one per configured gas.
    call AddDatum(csv_row, 'TAU_NSR', separator)
    call AddDatum(csv_row, 'H_NSR', separator)
    do j = 1, nFluxnetLayoutSlots
        call AddDatum(csv_row, trim(FluxnetFluxTag(j)) // '_NSR', separator)
    end do

    !> Foken statistics: the steady-state measure per flux, then the integral
    !> turbulence characteristics on the three anemometric variables.
    call AddDatum(csv_row, 'TAU_SS', separator)
    call AddDatum(csv_row, 'H_SS', separator)
    do j = 1, nFluxnetLayoutSlots
        call AddDatum(csv_row, trim(FluxnetFluxTag(j)) // '_SS', separator)
    end do
    call AddDatum(csv_row, 'U_ITC', separator)
    call AddDatum(csv_row, 'W_ITC', separator)
    call AddDatum(csv_row, 'T_SONIC_ITC', separator)

    !> Partial Foken flags: the steady-state test per flux, then the integral
    !> turbulence characteristics test on the three anemometric variables.
    !> Copied verbatim by FCC as fluxnetChunks%s(3).
    call AddDatum(csv_row, 'TAU_SS_TEST', separator)
    call AddDatum(csv_row, 'H_SS_TEST', separator)
    do j = 1, nFluxnetLayoutSlots
        call AddDatum(csv_row, trim(FluxnetFluxTag(j)) // '_SS_TEST', separator)
    end do
    call AddDatum(csv_row, 'U_ITC_TEST', separator)
    call AddDatum(csv_row, 'W_ITC_TEST', separator)
    call AddDatum(csv_row, 'T_SONIC_ITC_TEST', separator)

    !> Final Foken flags, flux-shaped: momentum, sensible heat and the two
    !> water fluxes, then one per configured gas.
    call AddDatum(csv_row, 'TAU_SSITC_TEST', separator)
    call AddDatum(csv_row, 'H_SSITC_TEST', separator)
    call AddDatum(csv_row, 'LE_SSITC_TEST', separator)
    call AddDatum(csv_row, 'ET_SSITC_TEST', separator)
    do j = 1, nFluxnetLayoutSlots
        call AddDatum(csv_row, trim(FluxnetFluxTag(j)) // '_SSITC_TEST', separator)
    end do

    csv_row = trim(csv_row) // 'INST_LI7200_HEAD_DETECT,INST_LI7200_T_OUT,INST_LI7200_T_IN,INST_LI7200_AUX_IN,&
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

!> Emit one column per variable for a family: the wind components and sonic
!> temperature, then one per configured gas.
!>
!> Three families in this row share that shape exactly, and writing the loop
!> out three times is how the four-gas quadruples came to be duplicated in the
!> first place.
subroutine AddVariableFamily(suffix)
    character(*), intent(in) :: suffix
    integer :: n

    call AddDatum(csv_row, 'U' // suffix, separator)
    call AddDatum(csv_row, 'V' // suffix, separator)
    call AddDatum(csv_row, 'W' // suffix, separator)
    call AddDatum(csv_row, 'T_SONIC' // suffix, separator)
    do n = 1, nFluxnetLayoutSlots
        call AddDatum(csv_row, trim(FluxnetLayoutTags(n)) // suffix, separator)
    end do
end subroutine AddVariableFamily

!> Per-gas prefix of the record-count family: FC, LE, FCH4, F<TAG>.
!>
!> Water is LE here and FH2O everywhere else - the column counts the records
!> behind the latent-heat flux. One of three per-gas naming conventions in this
!> row; they cannot be collapsed without renaming shipped columns.
function FluxnetNrwTag(layout_index) result(tag)
    integer, intent(in) :: layout_index
    character(32) :: tag

    if (FluxnetLayoutSlots(layout_index) == h2o) then
        tag = 'LE'
    else
        tag = FluxnetFluxTag(layout_index)
    end if
end function FluxnetNrwTag

!> Per-gas prefix of the single-point storage family: SC, SH2O, SCH4, S<TAG>.
!>
!> Carbon dioxide is SC, not SCO2 - the S is prefixed to the flux name with its
!> F dropped, which only shows up as a difference for CO2.
function FluxnetStorTag(layout_index) result(tag)
    integer, intent(in) :: layout_index
    character(32) :: tag

    if (FluxnetLayoutSlots(layout_index) == co2) then
        tag = 'SC'
    else
        tag = 'S' // trim(FluxnetLayoutTags(layout_index))
    end if
end function FluxnetStorTag

!> One column per flux: momentum, sensible heat, the two water fluxes, then one
!> per configured gas under the flux naming.
subroutine AddFluxFamily(suffix)
    character(*), intent(in) :: suffix
    integer :: n

    call AddDatum(csv_row, 'TAU' // suffix, separator)
    call AddDatum(csv_row, 'H' // suffix, separator)
    call AddDatum(csv_row, 'LE' // suffix, separator)
    call AddDatum(csv_row, 'ET' // suffix, separator)
    do n = 1, nFluxnetLayoutSlots
        call AddDatum(csv_row, trim(FluxnetFluxTag(n)) // suffix, separator)
    end do
end subroutine AddFluxFamily

!> A per-variable statistic whose wind and gas halves take different suffixes:
!> U_MEDIAN beside CO2_MEAS_MEDIAN.
subroutine AddStatFamily(wind_suffix, gas_suffix)
    character(*), intent(in) :: wind_suffix
    character(*), intent(in) :: gas_suffix
    integer :: n

    call AddDatum(csv_row, 'U' // wind_suffix, separator)
    call AddDatum(csv_row, 'V' // wind_suffix, separator)
    call AddDatum(csv_row, 'W' // wind_suffix, separator)
    call AddDatum(csv_row, 'T_SONIC' // wind_suffix, separator)
    do n = 1, nFluxnetLayoutSlots
        call AddDatum(csv_row, trim(FluxnetLayoutTags(n)) // gas_suffix, separator)
    end do
end subroutine AddStatFamily

!> Flux-family prefix of a layout slot: FC, FH2O, FCH4, FCOS, ...
!>
!> These columns are named for the flux rather than for the species, and carbon
!> dioxide's flux is FC, not FCO2 - so the tag alone does not answer this. Every
!> other gas takes F followed by its tag, which is what the third and fourth
!> slots already spelled out as FCH4 and FGS4.
function FluxnetFluxTag(layout_index) result(tag)
    integer, intent(in) :: layout_index
    character(32) :: tag

    if (FluxnetLayoutSlots(layout_index) == co2) then
        tag = 'FC'
    else
        tag = 'F' // trim(FluxnetLayoutTags(layout_index))
    end if
end function FluxnetFluxTag

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
