!***************************************************************************
! m_common_global_var.f90
! -----------------------
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
! \brief       Contain declaration of all variables common to EddyFlow projects.
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
module m_common_global_var
    use m_typedef
    use m_dates
    use m_methane_tables
    use m_index_parameters
    use libdate
    use m_fp2_to_float
    implicit none
    save

    !> arrays and files defaults
    !> Max key lines in one parsed INI section. Section selection is a prefix
    !> match, so ReadIniRP('RawProcess') sweeps every [RawProcess_*] group into
    !> one buffer; the per-gas settings make that buffer grow with the gas count.
    integer, parameter :: MaxNLinesIni = 8000
    integer, parameter :: MaxNumAnem = 7
    integer, parameter :: MaxNumAdc = 10
    integer, parameter :: MaxNumVar = MaxNumAnem + MaxNumAdc
    integer, parameter :: NumSpecVar = 3
    integer, parameter :: MaxSpecRow = 72000
    integer, parameter :: NumPar = 3
    integer, parameter :: MaxNumDP = 100
    integer, parameter :: NTlagH2Oclasses = 20
    integer, parameter :: Nt = 9

    !> Validated variables
    integer :: NumInstruments
    integer :: NumRawFlags
    integer :: NumCol
    integer :: NumBiometCol
    integer :: nbVars
    integer :: nbItems
    integer :: NumAllVar
    integer :: NumVar
    integer :: NumDiag
    integer :: NumUserVar
    logical :: FileWithFlags
    integer :: n_cstm_biomet
    integer :: Gas4CalRefCol

    !> Platform management
    character(8) :: OS
    character(8) :: root
    character(1) :: slash
    character(1) :: escape
    character(16) :: comm_del
    character(16) :: comm_rmdir
    character(16) :: comm_err_redirect
    character(16) :: comm_out_redirect
    character(16) :: comm_7zip
    character(16) :: comm_7zip_x_opt
    character(16) :: comm_copy
    character(16) :: comm_move
    character(16) :: comm_force_opt
    character(15) :: comm_dir
    character(PathLen) :: homedir
    character(PathLen) :: IniDir
    character(PathLen) :: TmpDir
    character(PathLen) :: PrjPath

    character(19), parameter :: PrjFile   = 'processing.eddyflow'
    character(6), parameter :: licor_appdata = '.licor'
    character(22)  :: Timestamp_FilePadding
    character(8), parameter  :: EDDYFLOW_FilePadding    = 'eddyflow'
    character(12), parameter :: RP_FilePadding = 'eddyflow-rp_'
    character(13), parameter :: FX_FilePadding = 'eddyflow-fcc_'
    character(8),  parameter :: SC_FilePadding = '_spectra'
    character(7),  parameter :: EC_FilePadding = '_fluxes'
    character(8),  parameter :: PF_FilePadding = '_tilting'
    character(8),  parameter :: TO_FilePadding = '_timelag'
    character(25), parameter :: SubDirBinCospectra      = 'eddyflow_binned_cospectra'
    character(23), parameter :: SubDirCospectra         = 'eddyflow_full_cospectra'
    character(22), parameter :: RS_flags_FilePadding    = '_statistical_screening'
    character(13), parameter :: RS_spike_FilePadding    = '_spike_counts'
    character(23), parameter :: Rot2D_FilePadding       = '_double_rotation_angles'
    character(9),  parameter :: Metadata_FilePadding    = '_metadata'
    character(11), parameter :: QCdetails_FilePadding  = '_qc_details'
    character(16), parameter :: H2OCov_FilePadding      = '_h2o_covariances'
    character(9),  parameter :: Tlag_FilePadding        = '_timelags'
    character(14), parameter :: RH_FilePadding          = '_timelag_vs_rh'
    character(17), parameter :: Flux_FilePadding        = '_tentative_fluxes'
    character(17), parameter :: BinCospec_FilePadding   = '_binned_cospectra'
    character(14), parameter :: BinOgives_FilePadding   = '_binned_ogives'
    character(15), parameter :: Cospec_FilePadding      = '_full_cospectra'
    character(29), parameter :: DegT_FilePadding        = '_degraded_wt_covariances_time'
    character(24), parameter :: vDegT_FilePadding       = '_degraded_wt_covariances'
    character(18), parameter :: QC_FilePadding          = '_stationarity_test'
    character(12), parameter :: FullOut_FilePadding     = '_full_output'
    character(11), parameter :: PlanarFit_FilePadding   = '_planar_fit'
    character(12), parameter :: TimelagOpt_FilePadding  = '_timelag_opt'
    character(16), parameter :: PwbTimelagDiag_FilePadding = '_pwb_diagnostics'
    character(12), parameter :: PwbSummary_FilePadding     = '_pwb_summary'
    character(8),  parameter :: FLUXNET_FilePadding     = '_fluxnet'
    character(7),  parameter  :: Biomet_FilePadding     = '_biomet'
    character(14), parameter :: Quality_FilePadding     = '_quality_check'
    character(18), parameter :: WPL_FilePadding         = '_wpl_contributions'
    character(20), parameter :: BPCF_FilePadding        = '_bandpass_correction'
    character(21), parameter :: H2OAvrg_FilePadding     = '_h2o_ensemble_spectra'
    character(27), parameter :: Cosp_FilePadding        = '_ensemble_cospectra_by_time'
    character(29), parameter :: Stability_FilePadding   = '_ensemble_and_model_cospectra'
    character(31), parameter :: PASGAS_Avrg_FilePadding = '_passive_gases_ensemble_spectra'
    character(21), parameter :: CH4Avrg_FilePadding     = '_ch4_ensemble_spectra'
    character(21), parameter :: N2OAvrg_FilePadding     = '_n2o_ensemble_spectra'
    character(23), parameter :: LPCF_FilePadding        = '_spec_corr_model_params'
    character(22), parameter :: RHFCO_FilePadding       = '_h2o_cutoff_freq_vs_rh'
    character(20), parameter :: SA_FilePadding          = '_spectral_assessment'
    character(12),  parameter  :: Raw_FilePadding        = '_raw_dataset'
    character(4),  parameter  :: Stats1_FilePadding     = '_st1'
    character(4),  parameter  :: Stats2_FilePadding     = '_st2'
    character(4),  parameter  :: Stats3_FilePadding     = '_st3'
    character(4),  parameter  :: Stats4_FilePadding     = '_st4'
    character(4),  parameter  :: Stats5_FilePadding     = '_st5'
    character(4),  parameter  :: Stats6_FilePadding     = '_st6'
    character(4),  parameter  :: Stats7_FilePadding     = '_st7'
    character(9),  parameter  :: UserStats1_FilePadding = '_user_st1'
    character(9),  parameter  :: UserStats2_FilePadding = '_user_st2'
    character(9),  parameter  :: UserStats3_FilePadding = '_user_st3'
    character(9),  parameter  :: UserStats4_FilePadding = '_user_st4'
    character(9),  parameter  :: UserStats5_FilePadding = '_user_st5'
    character(9),  parameter  :: UserStats6_FilePadding = '_user_st6'
    character(9),  parameter  :: UserStats7_FilePadding = '_user_st7'
    character(17), parameter  :: PF1FilePadding         = '_pf_fitting_plane'
    character(21), parameter  :: PF2FilePadding         = '_pf_rotation_matrices'
    character(17), parameter  :: TlagOpt_FilePadding    = '_optimal_timelags'
    character(56), parameter  :: BinnedFilePrototype    = 'yyyymmdd-HHMM_xxxxxx_xxxxxxxxx_xxxx-xx-xxTxxxxxx_xxx.csv'
    character(54), parameter  :: FullFilePrototype      = 'yyyymmdd-HHMM_xxxx_xxxxxxxxx_xxxx-xx-xxTxxxxxx_xxx.csv'

    character(PathLen) :: FLUXNET_Path
    character(PathLen) :: FullOut_Path
    character(PathLen) :: Metadata_Path

    !> physical params and other useful numbers
    integer :: mmm
    real(kind = dbl) :: Dc(E2NumVar) !< Diffus. coeff. of gases in air [m+2 s-1]
    data (Dc(mmm), mmm = co2, gas4) / 0.00001381d0, 0.00002178d0, 0.00001952d0, 0.00001436d0/ !--> Massman (1998, Atm Env, Table 2)
    real(kind = sgl) :: MW(E2NumVar) !< Molecular weights
    data (MW(mmm), mmm = co2, gas4) / 44.01e-3, 18.02e-3, 16.04e-3, 44.01e-3/
    real(kind = dbl), parameter :: h2o_to_ET =  0.0648d0  !< To convert between H2O flux [mmol m-2 s-1] and ET flux (mm  hour-1)
    real(kind = dbl), parameter :: p = 3.14159265358979323846d0 !< Greek pi
    real(kind = dbl), parameter :: StdVair = 0.02245d0  !< gas molar volume at 25 �C and 101.325 kPa
    real(kind = dbl), parameter :: vk = 0.41d0 !< Von Karman constant
    real(kind = dbl), parameter :: Md = 0.02897d0 !< molecular weight of dry air [kg_d/mol_d]
    real(kind = dbl), parameter :: mu = Md / 18.02d-3
    real(kind = dbl), parameter :: kg_gamma = 0.95d0 !< for H correction after Kaimal and Gaynor (1991).
    real(kind = dbl), parameter :: g  = 9.81d0 !< gravity
    real(kind = dbl), parameter :: Ru = 8.314d0 !< universal gas constant J/[mol K]
    real(kind = dbl), parameter :: Rd = 287.04d0 !< gas constant for dry air [J/kg K]
    real(kind = dbl), parameter :: Rw = 461.5d0 !< gas constant for water vapour [J/kg K]
    real(kind = dbl), parameter :: RHmax = 130d0 !< max acceptable RH for keep doing calculations
    real(kind = dbl), parameter :: kj_us_min = 0.2d0 !< minimum ustar for Kljun model
    real(kind = dbl), parameter :: kj_zL_min = -200d0 !< minimum zL for Kljun model
    real(kind = dbl), parameter :: kj_zL_max = 1d0 !< minimum zL for Kljun model
    real(kind = dbl), parameter :: error = -9999.d0 !< main error label float
    integer, parameter :: ierror = -9999 !< main error label int
    real(kind = dbl), parameter :: aflx_error = -6999.d0 !< ameriflux error label
    real(kind = dbl), parameter :: MaxNormSpecValue = 1d4 !< maximum plausible value for a normalized spectral value
    real(kind = dbl), parameter :: MaxSpecValue = 1d4 !< maximum plausible value for an un-normalized spectral value
    real(kind = dbl), parameter :: MaxWindIntensity = 5d2 !< maximum plausible value for wind speed
    real(kind = dbl), parameter :: MaxWTCov = 100d0 !< maximum plausible value for wind speed
    !> Co-spectral model parameters (Runkle et al. 2012, Eq. 3)
    real(kind = dbl) , parameter :: beta1 = 1.05d0
    real(kind = dbl) , parameter :: beta2 = 1.33d0
    real(kind = dbl) , parameter :: beta3 = 0.387d0
    real(kind = dbl) , parameter :: beta4 = 0.38d0
    real(kind = dbl) , allocatable :: xFit(:)
    real(kind = dbl) , allocatable :: yFit(:)
    real(kind = dbl) , allocatable :: zFit(:)
    real(kind = dbl) , allocatable :: zzFit(:)
    real(kind = dbl) , allocatable :: ddum(:)

    type(RUsetupType) :: RUsetup
    type(FootType) :: Foot
    type(EddyFlowLogType)   :: EddyFlowLog
    type(EddyFlowProjType) :: EddyFlowProj
    type(SpectralType) :: BPCF
    type(SpectralType), parameter :: &
        errBPCF = spectraltype(error)
    type(SpectralType) :: ADDCF
    type(fluxnetChunksType) :: fluxnetChunks

    !> Variables to be validate
    real(kind = dbl) :: PFMat(3, 3, MaxNumWSect) = 0.d0
    real(kind = dbl) :: PFb(3, MaxNumWSect) = 0.d0
    real(kind = dbl) :: ITS(E2NumVar)

    !> filename tools
    character(4), parameter  :: CsvExt                  = '.csv'
    character(4), parameter  :: TmpExt                  = '.tmp'
    character(8), parameter  :: CsvTmpExt               = '.csv.tmp'
    character(4), parameter  :: TxtExt                  = '.txt'

    !> logging variables and parameters
    character(10) :: LogInteger
    logical :: LogAll = .false. !< working variable, for debug only
    logical :: co2_new_sw_ver = .false.

    integer, parameter :: ErrLab1 = 2

    !> labels for standard global set
    integer, parameter :: gU   = 1
    integer, parameter :: gV   = 2
    integer, parameter :: gW   = 3
    integer, parameter :: gTs  = 4
    integer, parameter :: gSoS = 5
    integer, parameter :: gIntC= 6
    integer, parameter :: gTa  = 7

    integer, parameter :: gCO2 = 8
    integer, parameter :: gH2O = 9
    integer, parameter :: gCH4 = 10
    integer, parameter :: gN2O = 11
    integer, parameter :: gTc  = 12
    integer, parameter :: gTi1 = 13
    integer, parameter :: gTi2 = 14
    integer, parameter :: gPi  = 15
    integer, parameter :: gTe  = 16
    integer, parameter :: gPe  = 17

    type(DateType), parameter :: &
        nullTimestamp = DateType(0, 0, 0, 0, 0)
    type(SpectraSetType), parameter :: &
        ErrSpec = SpectraSetType(0, error, error, error)
    type(SpectraSetType), parameter :: &
        NullSpec = SpectraSetType(0, 0d0, 0d0, 0d0)
    type(MeanSpectraType), parameter :: &
        NullMeanSpec = MeanSpectraType(0, 0, 0d0, 0d0, 0d0)
    type(FitSpectraType), parameter :: &
        NullFitCosp = FitSpectraType(0d0, 0d0)

    real(kind = dbl) :: StdFco(9)
    data (StdFco(mmm), mmm = 1, 9) / 0.004d0, 0.008d0, 0.016d0, 0.032d0, 0.065d0, 0.133d0, &
                                0.277d0, 0.614d0, 1.626d0 /

    integer, parameter :: NumStdDynMDVars = 75
    character(64) :: StdDynMDVars(NumStdDynMDVars)
    data (StdDynMDVars(mmm), mmm = 1, NumStdDynMDVars) /'date', 'time', 'latitude', 'longitude', 'altitude',&
                'file_length', 'acquisition_frequency', &
                'canopy_height', 'displacement_height', 'roughness_length', &
                'master_sonic_manufacturer', 'master_sonic_model', 'master_sonic_height', &
                'master_sonic_wformat', 'master_sonic_wref', 'master_sonic_north_offset', &
                'master_sonic_hpath_length', 'master_sonic_vpath_length', 'master_sonic_tau', &
                'co2_irga_manufacturer', 'co2_irga_model', 'co2_measure_type', &
                'co2_irga_northward_separation', 'co2_irga_eastward_separation', 'co2_irga_vertical_separation', &
                'co2_irga_tube_length', 'co2_irga_tube_diameter', 'co2_irga_tube_flowrate',  &
                'co2_irga_kw', 'co2_irga_ko', 'co2_irga_hpath_length', 'co2_irga_vpath_length', 'co2_irga_tau',  &
                'h2o_irga_manufacturer', 'h2o_irga_model', 'h2o_measure_type', &
                'h2o_irga_northward_separation', 'h2o_irga_eastward_separation', 'h2o_irga_vertical_separation', &
                'h2o_irga_tube_length', 'h2o_irga_tube_diameter', 'h2o_irga_tube_flowrate',  &
                'h2o_irga_kw', 'h2o_irga_ko', 'h2o_irga_hpath_length', 'h2o_irga_vpath_length', 'h2o_irga_tau',  &
                'ch4_irga_manufacturer', 'ch4_irga_model', 'ch4_measure_type', &
                'ch4_irga_northward_separation', 'ch4_irga_eastward_separation', 'ch4_irga_vertical_separation', &
                'ch4_irga_tube_length', 'ch4_irga_tube_diameter', 'ch4_irga_tube_flowrate',  &
                'ch4_irga_kw', 'ch4_irga_ko', 'ch4_irga_hpath_length', 'ch4_irga_vpath_length', 'ch4_irga_tau',  &
                'gas4_irga_manufacturer', 'gas4_irga_model', 'gas4_measure_type', &
                'gas4_irga_northward_separation', 'gas4_irga_eastward_separation', 'gas4_irga_vertical_separation', &
                'gas4_irga_tube_length', 'gas4_irga_tube_diameter', 'gas4_irga_tube_flowrate',  &
                'gas4_irga_kw', 'gas4_irga_ko', 'gas4_irga_hpath_length', 'gas4_irga_vpath_length', 'gas4_irga_tau' /


!    integer, parameter :: NumStdUnits = 107
!    character(32) :: StdUnits(NumStdUnits)
!    data (StdUnits(mmm), mmm = 1, NumStdUnits) &
!    /'K','PA','%','W+1M-2','UMOL+1M-2S-1','M','M+1S-1','PPM', 'PPT','PPB','DEG','NONE', &   !> standard units
!    'N+1M-2','NM^-2','N/M2','N/M^2','M/S','MS^-1','MS-1',& !> units that don't need to be changed
!    'WM^-2','W/M2','WM-2','W/M^2','WATTM^-2','WATT/M2','WATT/M^2',& !> units that don't need to be changed
!    'J/M2S','JM-2S-1','JM^-2S^-1','J/(M^2*S)', & !> units that don't need to be changed
!    'UMOLM-2S-1','UMOL/(M^2*S)','UMOLM^-2*S^-1','UMOL/M^2/S^1',  & !> units that don't need to be changed
!    'UE+1M-2S-1','UE/(M^2*S)','UEM^-2*S^-1',  & !> units that don't need to be changed
!    'MICROEINSTEIN+1M-2S-1','MICROEINSTEINM-2S-1','MICROEINSTEIN/(M^2*S)','MICROEINSTEINM^-2*S^-1',  &
!    'UEINSTEIN+1M-2S-1','UEINSTEINM-2S-1','UEINSTEIN/(M^2*S)','UEINSTEINM^-2*S^-1',  &
!    'PPMD','UMOLMOL-1','UMOL/MOL','UMOLMOL^-1', & !> units that don't need to be changed
!    '�','�DEG','DEGREES','DEGREESFROMNORTH', & !> units that don't need to be changed
!    '#','PERCENT','%VOL', & !> units that don't need to be changed
!    'M+3M-3','M3/M3','M^3M^-3','M^3/M^3','M3M-3',& !> units that don't need to be changed
!    'M+2M-2','M2/M2','M^2M^-2','M2M-2',& !> units that don't need to be changed
!    'NUMBER','#','DIMENSIONLESS','OTHER' ,'OTHERS', & !> units that don't need to be changed
!    'C','�C','F','�F','CK','CC','C�C','CF','C�F', & !> units that need change
!    'HPA','KPA','MMHG','PSI', 'BAR', 'ATM', 'TORR', & !> units that need change
!    'NM','UM','MM','CM','KM', & !> units that need change
!    'CM+1S-1','CM/S','CMS^-1','CMS-1','MM+1S-1','MM/S','MMS^-1','MMS-1', & !> units that need change
!    'PPTD','MMOLMOL-1','MMOL/MOL','MMOLMOL^-1', & !> units that need change
!    'PPBD','NMOLMOL-1','NMOL/MOL','NMOLMOL^-1'/ !> units that need change

    !> indentations and other strings
    character(0) :: indent0 = ''
    character(1) :: indent1 = ' '
    character(2) :: indent2 = '  '
    character(3) :: indent3 = '   '
    character(4) :: indent4 = '    '
    character(1) :: separator = ','

    type(SwVerType), parameter :: &
        errSwVer = SwVerType(ierror, ierror, ierror)

    !> Newest project-file ('ini_version') format this engine can read. The GUI
    !> stamps the format version into every project file it writes; a file
    !> newer than this is refused rather than half-parsed. Bump this in step
    !> with Defs::PROJECT_FILE_VERSION_STR on the GUI side.
    character(5), parameter :: MaxSupportedIniVer = '4.5.0'
    type(DateType), parameter :: &
        tsNull= DateType(0, 0, 0, 0, 0)
    type(InstrumentType), parameter :: &
        NullInstrument = InstrumentType('none', 'none', 'none', 'none', 'none', &
        error, error, error, error, &
        error, error, error, error, error, error, error, error, error, error, 'none', 'none', &
        'none', 'none', .false., .false., errSwVer)
    type(RawFlagType), parameter :: NullRawFlag = RawFlagType(0, error, .false.)
    type(ColType), parameter :: &
        NullCol = Coltype('none', 'none', 'none', 'none', '', '', '', &
        0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, NullInstrument, NullRawFlag, .false., .false., error, &
        0, 0, 0)
    type(BiometColType), parameter :: &
        NullBiometCol = BiometColType('none', 'none', 'none', 'none', 'none', error, error)

    !> Keyword form rather than positional: the gas members are arrays now, so
    !> a positional list would bind the wrong components without complaint if
    !> the type ever gains or reorders a field.
    type(fluxtype), parameter :: &
        errFlux = fluxtype(date = '', time = '', &
            gas = error, E = error, ET = error, E_gas = error, &
            LE = error, H = error, Hi_gas = error, &
            tau = error, ustar = error, L = error, zL = error)

    integer :: RowLags(E2NumVar)
    integer :: UserRowLags(MaxUserVar)
    type(MethType) :: Meth
    type(DirType) :: Dir
    type(FileType) :: AuxFile
    type(BinaryType) :: Binary
    type(MetadataType) :: Metadata

    integer :: fnbRecs, nbRecs
    type(BiometFileMetadataType) :: bFileMetadata
    type(BiometVarsType), allocatable :: bVars(:)
    real(kind = dbl), allocatable :: fbSet(:, :)
    real(kind = dbl), allocatable :: bSet(:, :)
    real(kind = dbl), allocatable :: auxbSet(:, :)
    real(kind = dbl), allocatable :: bAggr(:)
    real(kind = dbl), allocatable :: bAggrFluxnet(:)
    real(kind = dbl), allocatable :: bAggrEddyFlow(:)
    type(DateType), allocatable :: fbTs(:), auxbTs(:)
    type(DateType), allocatable :: bTs(:)
    type(BiometVarsType), parameter :: &
        nullbVar = BiometVarsType(nint(error), nint(error), nint(error), &
            error, error, 'none', 'none', 'none', 'none', 'none', 'none', 'none', 'none', &
            'none', 'none', 'none', 'none','none')
    type(FootType), parameter :: &
        errFootprint = FootType(error, error, error, error, error, &
            error, error, error)
    !> variables from metadata file
    type(RawFlagType)   :: RawFlag(MaxNumRawFlags)
    type(InstrumentType)   :: Instr(MaxNumInstruments)
    type(FileInterpreterType)   :: FileInterpreter
    type(ColType) :: Col(MaxNumCol)
    type(ColType) :: E2Col(E2NumVar)
    type(ColType) :: SpecCol(E2NumVar)
    type(ColType), allocatable :: UserCol(:)
    type(EssentialsType)  :: Essentials
    logical :: ArchiveIsCreated = .false.

    !> other shared variables
    real(kind = dbl)  :: PotRad(17568)
    real(kind = dbl)  :: magnetic_declination
    real(kind = dbl) BuMultiPar(2, 3, 4)
    character(12)  :: app
    character(32)  :: TFShape
    character(32) :: foot_model_used
    Type(StatsType)  :: Stats
    Type(UserStatsType)  :: UserStats
    type(AmbientStateType)    :: Ambient
    type(RegParType) :: RegPar(GHGNumVar, MaxGasClasses)
    type(DateType)   :: DateStep
    type(DateType)   :: DatafileDateStep
    type(GenericE2Var) :: OptTLagDef
    type(GenericE2Var) :: OptTLagStDev
    type(Diag7200Type) :: Diag7200
    type(Diag7500Type) :: Diag7500
    type(Diag7700Type) :: Diag7700
    type(DiagAnemType) :: DiagAnemometer
    type (QCType) :: QCFlag
    real(kind = dbl) :: f_c(GHGNumVar)
    real(kind = dbl) :: f_2(GHGNumVar)
    real(kind = dbl) :: StPar(2) = error
    real(kind = dbl) :: UnPar(2) = error

    !> tags of the [Project] group of processing.eddypro file
    integer, parameter :: Npn = 403
    integer, parameter :: Npc = 274
    !> BEGIN GENERATED ProjectRecordOrigins - edit gen_project_tags.py, not this block
    !> Slot origins for the appended gas/cell/diag records. The value
    !> is the index of the FIRST field of record 1, so record i field f
    !> is <origin> + (i-1)*<leap> + f, with f zero-based.
    integer, parameter :: gasNumTag = 33
    integer, parameter :: cellNumTag = 34
    integer, parameter :: diagNumTag = 35
    integer, parameter :: gasRecOriginN = 36
    integer, parameter :: cellRecOriginN = 356
    integer, parameter :: diagRecOriginN = 388
    integer, parameter :: gasRecOriginC = 51
    integer, parameter :: cellRecOriginC = 179
    integer, parameter :: diagRecOriginC = 243
    integer, parameter :: rpGasOriginN = 425
    integer, parameter :: rpGasOriginC = 102
    integer, parameter :: fccGasOriginN = 110
    integer, parameter :: gasRecLeapN   = 5
    integer, parameter :: gasRecLeapC   = 2
    integer, parameter :: cellRecLeapN  = 1
    integer, parameter :: cellRecLeapC  = 2
    integer, parameter :: diagRecLeapN  = 1
    integer, parameter :: diagRecLeapC  = 2
    integer, parameter :: rpGasLeapN    = 11
    integer, parameter :: rpGasLeapC    = 3
    integer, parameter :: fccGasLeapN   = 6
    !> END GENERATED ProjectRecordOrigins
    logical :: EPPrjNTagFound(Npn)
    logical :: EPPrjCTagFound(Npc)
    type (Numerical) :: EPPrjNTags(Npn)
    type (Text) :: EPPrjCTags(Npc)
    !> BEGIN GENERATED EPPrjNTags - edit gen_project_tags.py, not this block
    data EPPrjNTags(1)%Label / 'binary_nbytes' / &
         EPPrjNTags(2)%Label / 'binary_hnlines' / &
         EPPrjNTags(3)%Label / 'col_ts' / &
         EPPrjNTags(4)%Label / 'col_co2' / &
         EPPrjNTags(5)%Label / 'col_h2o' / &
         EPPrjNTags(6)%Label / 'col_ch4' / &
         EPPrjNTags(7)%Label / 'col_gas4' / &
         EPPrjNTags(8)%Label / 'col_cell_t' / &
         EPPrjNTags(9)%Label / 'col_int_t_1' / &
         EPPrjNTags(10)%Label / 'col_int_t_2' / &
         EPPrjNTags(11)%Label / 'col_int_p' / &
         EPPrjNTags(12)%Label / 'col_air_t' / &
         EPPrjNTags(13)%Label / 'col_air_p' / &
         EPPrjNTags(14)%Label / 'col_diag_72' / &
         EPPrjNTags(15)%Label / 'col_diag_75' / &
         EPPrjNTags(16)%Label / 'col_diag_77' / &
         EPPrjNTags(17)%Label / 'gas_diff' / &
         EPPrjNTags(18)%Label / 'gas_mw' / &
         EPPrjNTags(19)%Label / 'sonic_output_rate' / &
         EPPrjNTags(20)%Label / 'col_diag_anem' / &
         EPPrjNTags(21)%Label / 'col_diag_staa' / &
         EPPrjNTags(22)%Label / 'col_diag_stad' / &
         EPPrjNTags(23)%Label / 'ru_its_meth' / &
         EPPrjNTags(24)%Label / 'ru_meth' / &
         EPPrjNTags(25)%Label / 'ru_tlag_max' / &
         EPPrjNTags(26)%Label / 'cec_h' / &
         EPPrjNTags(27)%Label / 'cec_min_o1_o2' / &
         EPPrjNTags(28)%Label / 'cec_min_octant' / &
         EPPrjNTags(29)%Label / 'cec_min_valid' / &
         EPPrjNTags(30)%Label / 'cec_signal_strength' / &
         EPPrjNTags(31)%Label / 'cec_max_gap_fill' / &
         EPPrjNTags(32)%Label / 'cec_max_stationarity' / &
         EPPrjNTags(33)%Label / 'gas_num' / &
         EPPrjNTags(34)%Label / 'cell_num' / &
         EPPrjNTags(35)%Label / 'diag_num' / &
         EPPrjNTags(36)%Label / 'gas_1_col' / &
         EPPrjNTags(37)%Label / 'gas_1_moist' / &
         EPPrjNTags(38)%Label / 'gas_1_cell' / &
         EPPrjNTags(39)%Label / 'gas_1_mw' / &
         EPPrjNTags(40)%Label / 'gas_1_diff' / &
         EPPrjNTags(41)%Label / 'gas_2_col' / &
         EPPrjNTags(42)%Label / 'gas_2_moist' / &
         EPPrjNTags(43)%Label / 'gas_2_cell' / &
         EPPrjNTags(44)%Label / 'gas_2_mw' / &
         EPPrjNTags(45)%Label / 'gas_2_diff' / &
         EPPrjNTags(46)%Label / 'gas_3_col' / &
         EPPrjNTags(47)%Label / 'gas_3_moist' / &
         EPPrjNTags(48)%Label / 'gas_3_cell' / &
         EPPrjNTags(49)%Label / 'gas_3_mw' / &
         EPPrjNTags(50)%Label / 'gas_3_diff' / &
         EPPrjNTags(51)%Label / 'gas_4_col' / &
         EPPrjNTags(52)%Label / 'gas_4_moist' / &
         EPPrjNTags(53)%Label / 'gas_4_cell' / &
         EPPrjNTags(54)%Label / 'gas_4_mw' / &
         EPPrjNTags(55)%Label / 'gas_4_diff' / &
         EPPrjNTags(56)%Label / 'gas_5_col' / &
         EPPrjNTags(57)%Label / 'gas_5_moist' / &
         EPPrjNTags(58)%Label / 'gas_5_cell' / &
         EPPrjNTags(59)%Label / 'gas_5_mw' / &
         EPPrjNTags(60)%Label / 'gas_5_diff' / &
         EPPrjNTags(61)%Label / 'gas_6_col' / &
         EPPrjNTags(62)%Label / 'gas_6_moist' / &
         EPPrjNTags(63)%Label / 'gas_6_cell' / &
         EPPrjNTags(64)%Label / 'gas_6_mw' / &
         EPPrjNTags(65)%Label / 'gas_6_diff' / &
         EPPrjNTags(66)%Label / 'gas_7_col' / &
         EPPrjNTags(67)%Label / 'gas_7_moist' / &
         EPPrjNTags(68)%Label / 'gas_7_cell' / &
         EPPrjNTags(69)%Label / 'gas_7_mw' / &
         EPPrjNTags(70)%Label / 'gas_7_diff' / &
         EPPrjNTags(71)%Label / 'gas_8_col' / &
         EPPrjNTags(72)%Label / 'gas_8_moist' / &
         EPPrjNTags(73)%Label / 'gas_8_cell' / &
         EPPrjNTags(74)%Label / 'gas_8_mw' / &
         EPPrjNTags(75)%Label / 'gas_8_diff' / &
         EPPrjNTags(76)%Label / 'gas_9_col' / &
         EPPrjNTags(77)%Label / 'gas_9_moist' / &
         EPPrjNTags(78)%Label / 'gas_9_cell' / &
         EPPrjNTags(79)%Label / 'gas_9_mw' / &
         EPPrjNTags(80)%Label / 'gas_9_diff' / &
         EPPrjNTags(81)%Label / 'gas_10_col' / &
         EPPrjNTags(82)%Label / 'gas_10_moist' / &
         EPPrjNTags(83)%Label / 'gas_10_cell' / &
         EPPrjNTags(84)%Label / 'gas_10_mw' / &
         EPPrjNTags(85)%Label / 'gas_10_diff' / &
         EPPrjNTags(86)%Label / 'gas_11_col' / &
         EPPrjNTags(87)%Label / 'gas_11_moist' / &
         EPPrjNTags(88)%Label / 'gas_11_cell' / &
         EPPrjNTags(89)%Label / 'gas_11_mw' / &
         EPPrjNTags(90)%Label / 'gas_11_diff' / &
         EPPrjNTags(91)%Label / 'gas_12_col' / &
         EPPrjNTags(92)%Label / 'gas_12_moist' / &
         EPPrjNTags(93)%Label / 'gas_12_cell' / &
         EPPrjNTags(94)%Label / 'gas_12_mw' / &
         EPPrjNTags(95)%Label / 'gas_12_diff' / &
         EPPrjNTags(96)%Label / 'gas_13_col' / &
         EPPrjNTags(97)%Label / 'gas_13_moist' / &
         EPPrjNTags(98)%Label / 'gas_13_cell' / &
         EPPrjNTags(99)%Label / 'gas_13_mw' / &
         EPPrjNTags(100)%Label / 'gas_13_diff' / &
         EPPrjNTags(101)%Label / 'gas_14_col' / &
         EPPrjNTags(102)%Label / 'gas_14_moist' / &
         EPPrjNTags(103)%Label / 'gas_14_cell' / &
         EPPrjNTags(104)%Label / 'gas_14_mw' / &
         EPPrjNTags(105)%Label / 'gas_14_diff' / &
         EPPrjNTags(106)%Label / 'gas_15_col' / &
         EPPrjNTags(107)%Label / 'gas_15_moist' / &
         EPPrjNTags(108)%Label / 'gas_15_cell' / &
         EPPrjNTags(109)%Label / 'gas_15_mw' / &
         EPPrjNTags(110)%Label / 'gas_15_diff' / &
         EPPrjNTags(111)%Label / 'gas_16_col' / &
         EPPrjNTags(112)%Label / 'gas_16_moist' / &
         EPPrjNTags(113)%Label / 'gas_16_cell' / &
         EPPrjNTags(114)%Label / 'gas_16_mw' / &
         EPPrjNTags(115)%Label / 'gas_16_diff' / &
         EPPrjNTags(116)%Label / 'gas_17_col' / &
         EPPrjNTags(117)%Label / 'gas_17_moist' / &
         EPPrjNTags(118)%Label / 'gas_17_cell' / &
         EPPrjNTags(119)%Label / 'gas_17_mw' / &
         EPPrjNTags(120)%Label / 'gas_17_diff' / &
         EPPrjNTags(121)%Label / 'gas_18_col' / &
         EPPrjNTags(122)%Label / 'gas_18_moist' / &
         EPPrjNTags(123)%Label / 'gas_18_cell' / &
         EPPrjNTags(124)%Label / 'gas_18_mw' / &
         EPPrjNTags(125)%Label / 'gas_18_diff' / &
         EPPrjNTags(126)%Label / 'gas_19_col' / &
         EPPrjNTags(127)%Label / 'gas_19_moist' / &
         EPPrjNTags(128)%Label / 'gas_19_cell' / &
         EPPrjNTags(129)%Label / 'gas_19_mw' / &
         EPPrjNTags(130)%Label / 'gas_19_diff' / &
         EPPrjNTags(131)%Label / 'gas_20_col' / &
         EPPrjNTags(132)%Label / 'gas_20_moist' / &
         EPPrjNTags(133)%Label / 'gas_20_cell' / &
         EPPrjNTags(134)%Label / 'gas_20_mw' / &
         EPPrjNTags(135)%Label / 'gas_20_diff' / &
         EPPrjNTags(136)%Label / 'gas_21_col' / &
         EPPrjNTags(137)%Label / 'gas_21_moist' / &
         EPPrjNTags(138)%Label / 'gas_21_cell' / &
         EPPrjNTags(139)%Label / 'gas_21_mw' / &
         EPPrjNTags(140)%Label / 'gas_21_diff' / &
         EPPrjNTags(141)%Label / 'gas_22_col' / &
         EPPrjNTags(142)%Label / 'gas_22_moist' / &
         EPPrjNTags(143)%Label / 'gas_22_cell' / &
         EPPrjNTags(144)%Label / 'gas_22_mw' / &
         EPPrjNTags(145)%Label / 'gas_22_diff' / &
         EPPrjNTags(146)%Label / 'gas_23_col' / &
         EPPrjNTags(147)%Label / 'gas_23_moist' / &
         EPPrjNTags(148)%Label / 'gas_23_cell' / &
         EPPrjNTags(149)%Label / 'gas_23_mw' / &
         EPPrjNTags(150)%Label / 'gas_23_diff' / &
         EPPrjNTags(151)%Label / 'gas_24_col' / &
         EPPrjNTags(152)%Label / 'gas_24_moist' / &
         EPPrjNTags(153)%Label / 'gas_24_cell' / &
         EPPrjNTags(154)%Label / 'gas_24_mw' / &
         EPPrjNTags(155)%Label / 'gas_24_diff' / &
         EPPrjNTags(156)%Label / 'gas_25_col' / &
         EPPrjNTags(157)%Label / 'gas_25_moist' / &
         EPPrjNTags(158)%Label / 'gas_25_cell' / &
         EPPrjNTags(159)%Label / 'gas_25_mw' / &
         EPPrjNTags(160)%Label / 'gas_25_diff' / &
         EPPrjNTags(161)%Label / 'gas_26_col' / &
         EPPrjNTags(162)%Label / 'gas_26_moist' / &
         EPPrjNTags(163)%Label / 'gas_26_cell' / &
         EPPrjNTags(164)%Label / 'gas_26_mw' / &
         EPPrjNTags(165)%Label / 'gas_26_diff' / &
         EPPrjNTags(166)%Label / 'gas_27_col' / &
         EPPrjNTags(167)%Label / 'gas_27_moist' / &
         EPPrjNTags(168)%Label / 'gas_27_cell' / &
         EPPrjNTags(169)%Label / 'gas_27_mw' / &
         EPPrjNTags(170)%Label / 'gas_27_diff' / &
         EPPrjNTags(171)%Label / 'gas_28_col' / &
         EPPrjNTags(172)%Label / 'gas_28_moist' / &
         EPPrjNTags(173)%Label / 'gas_28_cell' / &
         EPPrjNTags(174)%Label / 'gas_28_mw' / &
         EPPrjNTags(175)%Label / 'gas_28_diff' / &
         EPPrjNTags(176)%Label / 'gas_29_col' / &
         EPPrjNTags(177)%Label / 'gas_29_moist' / &
         EPPrjNTags(178)%Label / 'gas_29_cell' / &
         EPPrjNTags(179)%Label / 'gas_29_mw' / &
         EPPrjNTags(180)%Label / 'gas_29_diff' / &
         EPPrjNTags(181)%Label / 'gas_30_col' / &
         EPPrjNTags(182)%Label / 'gas_30_moist' / &
         EPPrjNTags(183)%Label / 'gas_30_cell' / &
         EPPrjNTags(184)%Label / 'gas_30_mw' / &
         EPPrjNTags(185)%Label / 'gas_30_diff' / &
         EPPrjNTags(186)%Label / 'gas_31_col' / &
         EPPrjNTags(187)%Label / 'gas_31_moist' / &
         EPPrjNTags(188)%Label / 'gas_31_cell' / &
         EPPrjNTags(189)%Label / 'gas_31_mw' / &
         EPPrjNTags(190)%Label / 'gas_31_diff' / &
         EPPrjNTags(191)%Label / 'gas_32_col' / &
         EPPrjNTags(192)%Label / 'gas_32_moist' / &
         EPPrjNTags(193)%Label / 'gas_32_cell' / &
         EPPrjNTags(194)%Label / 'gas_32_mw' / &
         EPPrjNTags(195)%Label / 'gas_32_diff' / &
         EPPrjNTags(196)%Label / 'gas_33_col' / &
         EPPrjNTags(197)%Label / 'gas_33_moist' / &
         EPPrjNTags(198)%Label / 'gas_33_cell' / &
         EPPrjNTags(199)%Label / 'gas_33_mw' / &
         EPPrjNTags(200)%Label / 'gas_33_diff' /
    data EPPrjNTags(201)%Label / 'gas_34_col' / &
         EPPrjNTags(202)%Label / 'gas_34_moist' / &
         EPPrjNTags(203)%Label / 'gas_34_cell' / &
         EPPrjNTags(204)%Label / 'gas_34_mw' / &
         EPPrjNTags(205)%Label / 'gas_34_diff' / &
         EPPrjNTags(206)%Label / 'gas_35_col' / &
         EPPrjNTags(207)%Label / 'gas_35_moist' / &
         EPPrjNTags(208)%Label / 'gas_35_cell' / &
         EPPrjNTags(209)%Label / 'gas_35_mw' / &
         EPPrjNTags(210)%Label / 'gas_35_diff' / &
         EPPrjNTags(211)%Label / 'gas_36_col' / &
         EPPrjNTags(212)%Label / 'gas_36_moist' / &
         EPPrjNTags(213)%Label / 'gas_36_cell' / &
         EPPrjNTags(214)%Label / 'gas_36_mw' / &
         EPPrjNTags(215)%Label / 'gas_36_diff' / &
         EPPrjNTags(216)%Label / 'gas_37_col' / &
         EPPrjNTags(217)%Label / 'gas_37_moist' / &
         EPPrjNTags(218)%Label / 'gas_37_cell' / &
         EPPrjNTags(219)%Label / 'gas_37_mw' / &
         EPPrjNTags(220)%Label / 'gas_37_diff' / &
         EPPrjNTags(221)%Label / 'gas_38_col' / &
         EPPrjNTags(222)%Label / 'gas_38_moist' / &
         EPPrjNTags(223)%Label / 'gas_38_cell' / &
         EPPrjNTags(224)%Label / 'gas_38_mw' / &
         EPPrjNTags(225)%Label / 'gas_38_diff' / &
         EPPrjNTags(226)%Label / 'gas_39_col' / &
         EPPrjNTags(227)%Label / 'gas_39_moist' / &
         EPPrjNTags(228)%Label / 'gas_39_cell' / &
         EPPrjNTags(229)%Label / 'gas_39_mw' / &
         EPPrjNTags(230)%Label / 'gas_39_diff' / &
         EPPrjNTags(231)%Label / 'gas_40_col' / &
         EPPrjNTags(232)%Label / 'gas_40_moist' / &
         EPPrjNTags(233)%Label / 'gas_40_cell' / &
         EPPrjNTags(234)%Label / 'gas_40_mw' / &
         EPPrjNTags(235)%Label / 'gas_40_diff' / &
         EPPrjNTags(236)%Label / 'gas_41_col' / &
         EPPrjNTags(237)%Label / 'gas_41_moist' / &
         EPPrjNTags(238)%Label / 'gas_41_cell' / &
         EPPrjNTags(239)%Label / 'gas_41_mw' / &
         EPPrjNTags(240)%Label / 'gas_41_diff' / &
         EPPrjNTags(241)%Label / 'gas_42_col' / &
         EPPrjNTags(242)%Label / 'gas_42_moist' / &
         EPPrjNTags(243)%Label / 'gas_42_cell' / &
         EPPrjNTags(244)%Label / 'gas_42_mw' / &
         EPPrjNTags(245)%Label / 'gas_42_diff' / &
         EPPrjNTags(246)%Label / 'gas_43_col' / &
         EPPrjNTags(247)%Label / 'gas_43_moist' / &
         EPPrjNTags(248)%Label / 'gas_43_cell' / &
         EPPrjNTags(249)%Label / 'gas_43_mw' / &
         EPPrjNTags(250)%Label / 'gas_43_diff' / &
         EPPrjNTags(251)%Label / 'gas_44_col' / &
         EPPrjNTags(252)%Label / 'gas_44_moist' / &
         EPPrjNTags(253)%Label / 'gas_44_cell' / &
         EPPrjNTags(254)%Label / 'gas_44_mw' / &
         EPPrjNTags(255)%Label / 'gas_44_diff' / &
         EPPrjNTags(256)%Label / 'gas_45_col' / &
         EPPrjNTags(257)%Label / 'gas_45_moist' / &
         EPPrjNTags(258)%Label / 'gas_45_cell' / &
         EPPrjNTags(259)%Label / 'gas_45_mw' / &
         EPPrjNTags(260)%Label / 'gas_45_diff' / &
         EPPrjNTags(261)%Label / 'gas_46_col' / &
         EPPrjNTags(262)%Label / 'gas_46_moist' / &
         EPPrjNTags(263)%Label / 'gas_46_cell' / &
         EPPrjNTags(264)%Label / 'gas_46_mw' / &
         EPPrjNTags(265)%Label / 'gas_46_diff' / &
         EPPrjNTags(266)%Label / 'gas_47_col' / &
         EPPrjNTags(267)%Label / 'gas_47_moist' / &
         EPPrjNTags(268)%Label / 'gas_47_cell' / &
         EPPrjNTags(269)%Label / 'gas_47_mw' / &
         EPPrjNTags(270)%Label / 'gas_47_diff' / &
         EPPrjNTags(271)%Label / 'gas_48_col' / &
         EPPrjNTags(272)%Label / 'gas_48_moist' / &
         EPPrjNTags(273)%Label / 'gas_48_cell' / &
         EPPrjNTags(274)%Label / 'gas_48_mw' / &
         EPPrjNTags(275)%Label / 'gas_48_diff' / &
         EPPrjNTags(276)%Label / 'gas_49_col' / &
         EPPrjNTags(277)%Label / 'gas_49_moist' / &
         EPPrjNTags(278)%Label / 'gas_49_cell' / &
         EPPrjNTags(279)%Label / 'gas_49_mw' / &
         EPPrjNTags(280)%Label / 'gas_49_diff' / &
         EPPrjNTags(281)%Label / 'gas_50_col' / &
         EPPrjNTags(282)%Label / 'gas_50_moist' / &
         EPPrjNTags(283)%Label / 'gas_50_cell' / &
         EPPrjNTags(284)%Label / 'gas_50_mw' / &
         EPPrjNTags(285)%Label / 'gas_50_diff' / &
         EPPrjNTags(286)%Label / 'gas_51_col' / &
         EPPrjNTags(287)%Label / 'gas_51_moist' / &
         EPPrjNTags(288)%Label / 'gas_51_cell' / &
         EPPrjNTags(289)%Label / 'gas_51_mw' / &
         EPPrjNTags(290)%Label / 'gas_51_diff' / &
         EPPrjNTags(291)%Label / 'gas_52_col' / &
         EPPrjNTags(292)%Label / 'gas_52_moist' / &
         EPPrjNTags(293)%Label / 'gas_52_cell' / &
         EPPrjNTags(294)%Label / 'gas_52_mw' / &
         EPPrjNTags(295)%Label / 'gas_52_diff' / &
         EPPrjNTags(296)%Label / 'gas_53_col' / &
         EPPrjNTags(297)%Label / 'gas_53_moist' / &
         EPPrjNTags(298)%Label / 'gas_53_cell' / &
         EPPrjNTags(299)%Label / 'gas_53_mw' / &
         EPPrjNTags(300)%Label / 'gas_53_diff' / &
         EPPrjNTags(301)%Label / 'gas_54_col' / &
         EPPrjNTags(302)%Label / 'gas_54_moist' / &
         EPPrjNTags(303)%Label / 'gas_54_cell' / &
         EPPrjNTags(304)%Label / 'gas_54_mw' / &
         EPPrjNTags(305)%Label / 'gas_54_diff' / &
         EPPrjNTags(306)%Label / 'gas_55_col' / &
         EPPrjNTags(307)%Label / 'gas_55_moist' / &
         EPPrjNTags(308)%Label / 'gas_55_cell' / &
         EPPrjNTags(309)%Label / 'gas_55_mw' / &
         EPPrjNTags(310)%Label / 'gas_55_diff' / &
         EPPrjNTags(311)%Label / 'gas_56_col' / &
         EPPrjNTags(312)%Label / 'gas_56_moist' / &
         EPPrjNTags(313)%Label / 'gas_56_cell' / &
         EPPrjNTags(314)%Label / 'gas_56_mw' / &
         EPPrjNTags(315)%Label / 'gas_56_diff' / &
         EPPrjNTags(316)%Label / 'gas_57_col' / &
         EPPrjNTags(317)%Label / 'gas_57_moist' / &
         EPPrjNTags(318)%Label / 'gas_57_cell' / &
         EPPrjNTags(319)%Label / 'gas_57_mw' / &
         EPPrjNTags(320)%Label / 'gas_57_diff' / &
         EPPrjNTags(321)%Label / 'gas_58_col' / &
         EPPrjNTags(322)%Label / 'gas_58_moist' / &
         EPPrjNTags(323)%Label / 'gas_58_cell' / &
         EPPrjNTags(324)%Label / 'gas_58_mw' / &
         EPPrjNTags(325)%Label / 'gas_58_diff' / &
         EPPrjNTags(326)%Label / 'gas_59_col' / &
         EPPrjNTags(327)%Label / 'gas_59_moist' / &
         EPPrjNTags(328)%Label / 'gas_59_cell' / &
         EPPrjNTags(329)%Label / 'gas_59_mw' / &
         EPPrjNTags(330)%Label / 'gas_59_diff' / &
         EPPrjNTags(331)%Label / 'gas_60_col' / &
         EPPrjNTags(332)%Label / 'gas_60_moist' / &
         EPPrjNTags(333)%Label / 'gas_60_cell' / &
         EPPrjNTags(334)%Label / 'gas_60_mw' / &
         EPPrjNTags(335)%Label / 'gas_60_diff' / &
         EPPrjNTags(336)%Label / 'gas_61_col' / &
         EPPrjNTags(337)%Label / 'gas_61_moist' / &
         EPPrjNTags(338)%Label / 'gas_61_cell' / &
         EPPrjNTags(339)%Label / 'gas_61_mw' / &
         EPPrjNTags(340)%Label / 'gas_61_diff' / &
         EPPrjNTags(341)%Label / 'gas_62_col' / &
         EPPrjNTags(342)%Label / 'gas_62_moist' / &
         EPPrjNTags(343)%Label / 'gas_62_cell' / &
         EPPrjNTags(344)%Label / 'gas_62_mw' / &
         EPPrjNTags(345)%Label / 'gas_62_diff' / &
         EPPrjNTags(346)%Label / 'gas_63_col' / &
         EPPrjNTags(347)%Label / 'gas_63_moist' / &
         EPPrjNTags(348)%Label / 'gas_63_cell' / &
         EPPrjNTags(349)%Label / 'gas_63_mw' / &
         EPPrjNTags(350)%Label / 'gas_63_diff' / &
         EPPrjNTags(351)%Label / 'gas_64_col' / &
         EPPrjNTags(352)%Label / 'gas_64_moist' / &
         EPPrjNTags(353)%Label / 'gas_64_cell' / &
         EPPrjNTags(354)%Label / 'gas_64_mw' / &
         EPPrjNTags(355)%Label / 'gas_64_diff' / &
         EPPrjNTags(356)%Label / 'cell_1_col' / &
         EPPrjNTags(357)%Label / 'cell_2_col' / &
         EPPrjNTags(358)%Label / 'cell_3_col' / &
         EPPrjNTags(359)%Label / 'cell_4_col' / &
         EPPrjNTags(360)%Label / 'cell_5_col' / &
         EPPrjNTags(361)%Label / 'cell_6_col' / &
         EPPrjNTags(362)%Label / 'cell_7_col' / &
         EPPrjNTags(363)%Label / 'cell_8_col' / &
         EPPrjNTags(364)%Label / 'cell_9_col' / &
         EPPrjNTags(365)%Label / 'cell_10_col' / &
         EPPrjNTags(366)%Label / 'cell_11_col' / &
         EPPrjNTags(367)%Label / 'cell_12_col' / &
         EPPrjNTags(368)%Label / 'cell_13_col' / &
         EPPrjNTags(369)%Label / 'cell_14_col' / &
         EPPrjNTags(370)%Label / 'cell_15_col' / &
         EPPrjNTags(371)%Label / 'cell_16_col' / &
         EPPrjNTags(372)%Label / 'cell_17_col' / &
         EPPrjNTags(373)%Label / 'cell_18_col' / &
         EPPrjNTags(374)%Label / 'cell_19_col' / &
         EPPrjNTags(375)%Label / 'cell_20_col' / &
         EPPrjNTags(376)%Label / 'cell_21_col' / &
         EPPrjNTags(377)%Label / 'cell_22_col' / &
         EPPrjNTags(378)%Label / 'cell_23_col' / &
         EPPrjNTags(379)%Label / 'cell_24_col' / &
         EPPrjNTags(380)%Label / 'cell_25_col' / &
         EPPrjNTags(381)%Label / 'cell_26_col' / &
         EPPrjNTags(382)%Label / 'cell_27_col' / &
         EPPrjNTags(383)%Label / 'cell_28_col' / &
         EPPrjNTags(384)%Label / 'cell_29_col' / &
         EPPrjNTags(385)%Label / 'cell_30_col' / &
         EPPrjNTags(386)%Label / 'cell_31_col' / &
         EPPrjNTags(387)%Label / 'cell_32_col' / &
         EPPrjNTags(388)%Label / 'diag_1_col' / &
         EPPrjNTags(389)%Label / 'diag_2_col' / &
         EPPrjNTags(390)%Label / 'diag_3_col' / &
         EPPrjNTags(391)%Label / 'diag_4_col' / &
         EPPrjNTags(392)%Label / 'diag_5_col' / &
         EPPrjNTags(393)%Label / 'diag_6_col' / &
         EPPrjNTags(394)%Label / 'diag_7_col' / &
         EPPrjNTags(395)%Label / 'diag_8_col' / &
         EPPrjNTags(396)%Label / 'diag_9_col' / &
         EPPrjNTags(397)%Label / 'diag_10_col' / &
         EPPrjNTags(398)%Label / 'diag_11_col' / &
         EPPrjNTags(399)%Label / 'diag_12_col' / &
         EPPrjNTags(400)%Label / 'diag_13_col' /
    data EPPrjNTags(401)%Label / 'diag_14_col' / &
         EPPrjNTags(402)%Label / 'diag_15_col' / &
         EPPrjNTags(403)%Label / 'diag_16_col' /
    !> END GENERATED EPPrjNTags

    !> BEGIN GENERATED EPPrjCTags - edit gen_project_tags.py, not this block
    data EPPrjCTags(1)%Label / 'sw_version' / &
         EPPrjCTags(2)%Label / 'ini_version' / &
         EPPrjCTags(3)%Label / 'file_name' / &
         EPPrjCTags(4)%Label / 'project_title' / &
         EPPrjCTags(5)%Label / 'project_id' / &
         EPPrjCTags(6)%Label / 'file_type' / &
         EPPrjCTags(7)%Label / 'file_prototype' / &
         EPPrjCTags(8)%Label / 'cfg_file' / &
         EPPrjCTags(9)%Label / 'use_pfile' / &
         EPPrjCTags(10)%Label / 'proj_file' / &
         EPPrjCTags(11)%Label / 'use_dyn_md_file' / &
         EPPrjCTags(12)%Label / 'dyn_metadata_file' / &
         EPPrjCTags(13)%Label / 'binary_eol' / &
         EPPrjCTags(14)%Label / 'binary_little_end' / &
         EPPrjCTags(15)%Label / 'master_sonic' / &
         EPPrjCTags(16)%Label / 'run_mode' / &
         EPPrjCTags(17)%Label / 'use_biom' / &
         EPPrjCTags(18)%Label / 'biom_file' / &
         EPPrjCTags(21)%Label / 'out_rich' / &
         EPPrjCTags(22)%Label / 'lf_meth' / &
         EPPrjCTags(23)%Label / 'hf_meth' / &
         EPPrjCTags(24)%Label / 'make_dataset' / &
         EPPrjCTags(25)%Label / 'pr_start_date' / &
         EPPrjCTags(26)%Label / 'pr_start_time' / &
         EPPrjCTags(27)%Label / 'pr_end_date' / &
         EPPrjCTags(28)%Label / 'pr_end_time' / &
         EPPrjCTags(29)%Label / 'biom_dir' / &
         EPPrjCTags(30)%Label / 'biom_ext' / &
         EPPrjCTags(31)%Label / 'biom_rec' / &
         EPPrjCTags(32)%Label / 'tob1_format' / &
         EPPrjCTags(33)%Label / 'wpl_meth' / &
         EPPrjCTags(34)%Label / 'foot_meth' / &
         EPPrjCTags(35)%Label / 'out_path' / &
         EPPrjCTags(36)%Label / 'err_label' / &
         EPPrjCTags(37)%Label / 'fix_out_format' / &
         EPPrjCTags(38)%Label / 'qc_meth' / &
         EPPrjCTags(39)%Label / 'out_metadata' / &
         EPPrjCTags(40)%Label / 'pr_subset' / &
         EPPrjCTags(41)%Label / 'out_mean_cosp' / &
         EPPrjCTags(42)%Label / 'out_biomet' / &
         EPPrjCTags(43)%Label / 'out_mean_spec' / &
         EPPrjCTags(44)%Label / 'bin_sp_avail' / &
         EPPrjCTags(45)%Label / 'full_sp_avail' / &
         EPPrjCTags(46)%Label / 'hf_correct_ghg_ba' / &
         EPPrjCTags(47)%Label / 'hf_correct_ghg_zoh' / &
         EPPrjCTags(48)%Label / 'fluxnet_standardize_biomet' / &
         EPPrjCTags(49)%Label / 'fluxnet_err_label' / &
         EPPrjCTags(50)%Label / 'cec_meth' / &
         EPPrjCTags(51)%Label / 'gas_1_var' / &
         EPPrjCTags(52)%Label / 'gas_1_instr' / &
         EPPrjCTags(53)%Label / 'gas_2_var' / &
         EPPrjCTags(54)%Label / 'gas_2_instr' / &
         EPPrjCTags(55)%Label / 'gas_3_var' / &
         EPPrjCTags(56)%Label / 'gas_3_instr' / &
         EPPrjCTags(57)%Label / 'gas_4_var' / &
         EPPrjCTags(58)%Label / 'gas_4_instr' / &
         EPPrjCTags(59)%Label / 'gas_5_var' / &
         EPPrjCTags(60)%Label / 'gas_5_instr' / &
         EPPrjCTags(61)%Label / 'gas_6_var' / &
         EPPrjCTags(62)%Label / 'gas_6_instr' / &
         EPPrjCTags(63)%Label / 'gas_7_var' / &
         EPPrjCTags(64)%Label / 'gas_7_instr' / &
         EPPrjCTags(65)%Label / 'gas_8_var' / &
         EPPrjCTags(66)%Label / 'gas_8_instr' / &
         EPPrjCTags(67)%Label / 'gas_9_var' / &
         EPPrjCTags(68)%Label / 'gas_9_instr' / &
         EPPrjCTags(69)%Label / 'gas_10_var' / &
         EPPrjCTags(70)%Label / 'gas_10_instr' / &
         EPPrjCTags(71)%Label / 'gas_11_var' / &
         EPPrjCTags(72)%Label / 'gas_11_instr' / &
         EPPrjCTags(73)%Label / 'gas_12_var' / &
         EPPrjCTags(74)%Label / 'gas_12_instr' / &
         EPPrjCTags(75)%Label / 'gas_13_var' / &
         EPPrjCTags(76)%Label / 'gas_13_instr' / &
         EPPrjCTags(77)%Label / 'gas_14_var' / &
         EPPrjCTags(78)%Label / 'gas_14_instr' / &
         EPPrjCTags(79)%Label / 'gas_15_var' / &
         EPPrjCTags(80)%Label / 'gas_15_instr' / &
         EPPrjCTags(81)%Label / 'gas_16_var' / &
         EPPrjCTags(82)%Label / 'gas_16_instr' / &
         EPPrjCTags(83)%Label / 'gas_17_var' / &
         EPPrjCTags(84)%Label / 'gas_17_instr' / &
         EPPrjCTags(85)%Label / 'gas_18_var' / &
         EPPrjCTags(86)%Label / 'gas_18_instr' / &
         EPPrjCTags(87)%Label / 'gas_19_var' / &
         EPPrjCTags(88)%Label / 'gas_19_instr' / &
         EPPrjCTags(89)%Label / 'gas_20_var' / &
         EPPrjCTags(90)%Label / 'gas_20_instr' / &
         EPPrjCTags(91)%Label / 'gas_21_var' / &
         EPPrjCTags(92)%Label / 'gas_21_instr' / &
         EPPrjCTags(93)%Label / 'gas_22_var' / &
         EPPrjCTags(94)%Label / 'gas_22_instr' / &
         EPPrjCTags(95)%Label / 'gas_23_var' / &
         EPPrjCTags(96)%Label / 'gas_23_instr' / &
         EPPrjCTags(97)%Label / 'gas_24_var' / &
         EPPrjCTags(98)%Label / 'gas_24_instr' / &
         EPPrjCTags(99)%Label / 'gas_25_var' / &
         EPPrjCTags(100)%Label / 'gas_25_instr' / &
         EPPrjCTags(101)%Label / 'gas_26_var' / &
         EPPrjCTags(102)%Label / 'gas_26_instr' / &
         EPPrjCTags(103)%Label / 'gas_27_var' / &
         EPPrjCTags(104)%Label / 'gas_27_instr' / &
         EPPrjCTags(105)%Label / 'gas_28_var' / &
         EPPrjCTags(106)%Label / 'gas_28_instr' / &
         EPPrjCTags(107)%Label / 'gas_29_var' / &
         EPPrjCTags(108)%Label / 'gas_29_instr' / &
         EPPrjCTags(109)%Label / 'gas_30_var' / &
         EPPrjCTags(110)%Label / 'gas_30_instr' / &
         EPPrjCTags(111)%Label / 'gas_31_var' / &
         EPPrjCTags(112)%Label / 'gas_31_instr' / &
         EPPrjCTags(113)%Label / 'gas_32_var' / &
         EPPrjCTags(114)%Label / 'gas_32_instr' / &
         EPPrjCTags(115)%Label / 'gas_33_var' / &
         EPPrjCTags(116)%Label / 'gas_33_instr' / &
         EPPrjCTags(117)%Label / 'gas_34_var' / &
         EPPrjCTags(118)%Label / 'gas_34_instr' / &
         EPPrjCTags(119)%Label / 'gas_35_var' / &
         EPPrjCTags(120)%Label / 'gas_35_instr' / &
         EPPrjCTags(121)%Label / 'gas_36_var' / &
         EPPrjCTags(122)%Label / 'gas_36_instr' / &
         EPPrjCTags(123)%Label / 'gas_37_var' / &
         EPPrjCTags(124)%Label / 'gas_37_instr' / &
         EPPrjCTags(125)%Label / 'gas_38_var' / &
         EPPrjCTags(126)%Label / 'gas_38_instr' / &
         EPPrjCTags(127)%Label / 'gas_39_var' / &
         EPPrjCTags(128)%Label / 'gas_39_instr' / &
         EPPrjCTags(129)%Label / 'gas_40_var' / &
         EPPrjCTags(130)%Label / 'gas_40_instr' / &
         EPPrjCTags(131)%Label / 'gas_41_var' / &
         EPPrjCTags(132)%Label / 'gas_41_instr' / &
         EPPrjCTags(133)%Label / 'gas_42_var' / &
         EPPrjCTags(134)%Label / 'gas_42_instr' / &
         EPPrjCTags(135)%Label / 'gas_43_var' / &
         EPPrjCTags(136)%Label / 'gas_43_instr' / &
         EPPrjCTags(137)%Label / 'gas_44_var' / &
         EPPrjCTags(138)%Label / 'gas_44_instr' / &
         EPPrjCTags(139)%Label / 'gas_45_var' / &
         EPPrjCTags(140)%Label / 'gas_45_instr' / &
         EPPrjCTags(141)%Label / 'gas_46_var' / &
         EPPrjCTags(142)%Label / 'gas_46_instr' / &
         EPPrjCTags(143)%Label / 'gas_47_var' / &
         EPPrjCTags(144)%Label / 'gas_47_instr' / &
         EPPrjCTags(145)%Label / 'gas_48_var' / &
         EPPrjCTags(146)%Label / 'gas_48_instr' / &
         EPPrjCTags(147)%Label / 'gas_49_var' / &
         EPPrjCTags(148)%Label / 'gas_49_instr' / &
         EPPrjCTags(149)%Label / 'gas_50_var' / &
         EPPrjCTags(150)%Label / 'gas_50_instr' / &
         EPPrjCTags(151)%Label / 'gas_51_var' / &
         EPPrjCTags(152)%Label / 'gas_51_instr' / &
         EPPrjCTags(153)%Label / 'gas_52_var' / &
         EPPrjCTags(154)%Label / 'gas_52_instr' / &
         EPPrjCTags(155)%Label / 'gas_53_var' / &
         EPPrjCTags(156)%Label / 'gas_53_instr' / &
         EPPrjCTags(157)%Label / 'gas_54_var' / &
         EPPrjCTags(158)%Label / 'gas_54_instr' / &
         EPPrjCTags(159)%Label / 'gas_55_var' / &
         EPPrjCTags(160)%Label / 'gas_55_instr' / &
         EPPrjCTags(161)%Label / 'gas_56_var' / &
         EPPrjCTags(162)%Label / 'gas_56_instr' / &
         EPPrjCTags(163)%Label / 'gas_57_var' / &
         EPPrjCTags(164)%Label / 'gas_57_instr' / &
         EPPrjCTags(165)%Label / 'gas_58_var' / &
         EPPrjCTags(166)%Label / 'gas_58_instr' / &
         EPPrjCTags(167)%Label / 'gas_59_var' / &
         EPPrjCTags(168)%Label / 'gas_59_instr' / &
         EPPrjCTags(169)%Label / 'gas_60_var' / &
         EPPrjCTags(170)%Label / 'gas_60_instr' / &
         EPPrjCTags(171)%Label / 'gas_61_var' / &
         EPPrjCTags(172)%Label / 'gas_61_instr' / &
         EPPrjCTags(173)%Label / 'gas_62_var' / &
         EPPrjCTags(174)%Label / 'gas_62_instr' / &
         EPPrjCTags(175)%Label / 'gas_63_var' / &
         EPPrjCTags(176)%Label / 'gas_63_instr' / &
         EPPrjCTags(177)%Label / 'gas_64_var' / &
         EPPrjCTags(178)%Label / 'gas_64_instr' / &
         EPPrjCTags(179)%Label / 'cell_1_var' / &
         EPPrjCTags(180)%Label / 'cell_1_instr' / &
         EPPrjCTags(181)%Label / 'cell_2_var' / &
         EPPrjCTags(182)%Label / 'cell_2_instr' / &
         EPPrjCTags(183)%Label / 'cell_3_var' / &
         EPPrjCTags(184)%Label / 'cell_3_instr' / &
         EPPrjCTags(185)%Label / 'cell_4_var' / &
         EPPrjCTags(186)%Label / 'cell_4_instr' / &
         EPPrjCTags(187)%Label / 'cell_5_var' / &
         EPPrjCTags(188)%Label / 'cell_5_instr' / &
         EPPrjCTags(189)%Label / 'cell_6_var' / &
         EPPrjCTags(190)%Label / 'cell_6_instr' / &
         EPPrjCTags(191)%Label / 'cell_7_var' / &
         EPPrjCTags(192)%Label / 'cell_7_instr' / &
         EPPrjCTags(193)%Label / 'cell_8_var' / &
         EPPrjCTags(194)%Label / 'cell_8_instr' / &
         EPPrjCTags(195)%Label / 'cell_9_var' / &
         EPPrjCTags(196)%Label / 'cell_9_instr' / &
         EPPrjCTags(197)%Label / 'cell_10_var' / &
         EPPrjCTags(198)%Label / 'cell_10_instr' / &
         EPPrjCTags(199)%Label / 'cell_11_var' / &
         EPPrjCTags(200)%Label / 'cell_11_instr' / &
         EPPrjCTags(201)%Label / 'cell_12_var' / &
         EPPrjCTags(202)%Label / 'cell_12_instr' /
    data EPPrjCTags(203)%Label / 'cell_13_var' / &
         EPPrjCTags(204)%Label / 'cell_13_instr' / &
         EPPrjCTags(205)%Label / 'cell_14_var' / &
         EPPrjCTags(206)%Label / 'cell_14_instr' / &
         EPPrjCTags(207)%Label / 'cell_15_var' / &
         EPPrjCTags(208)%Label / 'cell_15_instr' / &
         EPPrjCTags(209)%Label / 'cell_16_var' / &
         EPPrjCTags(210)%Label / 'cell_16_instr' / &
         EPPrjCTags(211)%Label / 'cell_17_var' / &
         EPPrjCTags(212)%Label / 'cell_17_instr' / &
         EPPrjCTags(213)%Label / 'cell_18_var' / &
         EPPrjCTags(214)%Label / 'cell_18_instr' / &
         EPPrjCTags(215)%Label / 'cell_19_var' / &
         EPPrjCTags(216)%Label / 'cell_19_instr' / &
         EPPrjCTags(217)%Label / 'cell_20_var' / &
         EPPrjCTags(218)%Label / 'cell_20_instr' / &
         EPPrjCTags(219)%Label / 'cell_21_var' / &
         EPPrjCTags(220)%Label / 'cell_21_instr' / &
         EPPrjCTags(221)%Label / 'cell_22_var' / &
         EPPrjCTags(222)%Label / 'cell_22_instr' / &
         EPPrjCTags(223)%Label / 'cell_23_var' / &
         EPPrjCTags(224)%Label / 'cell_23_instr' / &
         EPPrjCTags(225)%Label / 'cell_24_var' / &
         EPPrjCTags(226)%Label / 'cell_24_instr' / &
         EPPrjCTags(227)%Label / 'cell_25_var' / &
         EPPrjCTags(228)%Label / 'cell_25_instr' / &
         EPPrjCTags(229)%Label / 'cell_26_var' / &
         EPPrjCTags(230)%Label / 'cell_26_instr' / &
         EPPrjCTags(231)%Label / 'cell_27_var' / &
         EPPrjCTags(232)%Label / 'cell_27_instr' / &
         EPPrjCTags(233)%Label / 'cell_28_var' / &
         EPPrjCTags(234)%Label / 'cell_28_instr' / &
         EPPrjCTags(235)%Label / 'cell_29_var' / &
         EPPrjCTags(236)%Label / 'cell_29_instr' / &
         EPPrjCTags(237)%Label / 'cell_30_var' / &
         EPPrjCTags(238)%Label / 'cell_30_instr' / &
         EPPrjCTags(239)%Label / 'cell_31_var' / &
         EPPrjCTags(240)%Label / 'cell_31_instr' / &
         EPPrjCTags(241)%Label / 'cell_32_var' / &
         EPPrjCTags(242)%Label / 'cell_32_instr' / &
         EPPrjCTags(243)%Label / 'diag_1_var' / &
         EPPrjCTags(244)%Label / 'diag_1_instr' / &
         EPPrjCTags(245)%Label / 'diag_2_var' / &
         EPPrjCTags(246)%Label / 'diag_2_instr' / &
         EPPrjCTags(247)%Label / 'diag_3_var' / &
         EPPrjCTags(248)%Label / 'diag_3_instr' / &
         EPPrjCTags(249)%Label / 'diag_4_var' / &
         EPPrjCTags(250)%Label / 'diag_4_instr' / &
         EPPrjCTags(251)%Label / 'diag_5_var' / &
         EPPrjCTags(252)%Label / 'diag_5_instr' / &
         EPPrjCTags(253)%Label / 'diag_6_var' / &
         EPPrjCTags(254)%Label / 'diag_6_instr' / &
         EPPrjCTags(255)%Label / 'diag_7_var' / &
         EPPrjCTags(256)%Label / 'diag_7_instr' / &
         EPPrjCTags(257)%Label / 'diag_8_var' / &
         EPPrjCTags(258)%Label / 'diag_8_instr' / &
         EPPrjCTags(259)%Label / 'diag_9_var' / &
         EPPrjCTags(260)%Label / 'diag_9_instr' / &
         EPPrjCTags(261)%Label / 'diag_10_var' / &
         EPPrjCTags(262)%Label / 'diag_10_instr' / &
         EPPrjCTags(263)%Label / 'diag_11_var' / &
         EPPrjCTags(264)%Label / 'diag_11_instr' / &
         EPPrjCTags(265)%Label / 'diag_12_var' / &
         EPPrjCTags(266)%Label / 'diag_12_instr' / &
         EPPrjCTags(267)%Label / 'diag_13_var' / &
         EPPrjCTags(268)%Label / 'diag_13_instr' / &
         EPPrjCTags(269)%Label / 'diag_14_var' / &
         EPPrjCTags(270)%Label / 'diag_14_instr' / &
         EPPrjCTags(271)%Label / 'diag_15_var' / &
         EPPrjCTags(272)%Label / 'diag_15_instr' / &
         EPPrjCTags(273)%Label / 'diag_16_var' / &
         EPPrjCTags(274)%Label / 'diag_16_instr' /
    !> END GENERATED EPPrjCTags

    !> tags of the metadata file created by GHG software
    integer, parameter :: Nan = 929
    integer, parameter :: Nac = 790
    logical :: ANTagFound(Nan)
    logical :: ACTagFound(Nac)
    type (Numerical) :: ANTags(Nan)
    type (Text) :: ACTags(Nac)

    !> BEGIN GENERATED ANTags - edit gen_metadata_tags.py, not this block
    data ANTags(1)%Label / 'altitude' / &
         ANTags(2)%Label / 'latitude' / &
         ANTags(3)%Label / 'longitude' / &
         ANTags(4)%Label / 'canopy_height' / &
         ANTags(5)%Label / 'displacement_height' / &
         ANTags(6)%Label / 'roughness_length' / &
         ANTags(7)%Label / 'acquisition_frequency' / &
         ANTags(8)%Label / 'file_duration' / &
         ANTags(9)%Label / 'header_rows' / &
         ANTags(10)%Label / 'instr_1_height' / &
         ANTags(11)%Label / 'instr_1_north_offset' / &
         ANTags(12)%Label / 'instr_1_northward_separation' / &
         ANTags(13)%Label / 'instr_1_eastward_separation' / &
         ANTags(14)%Label / 'instr_1_vertical_separation' / &
         ANTags(15)%Label / 'instr_1_tube_diameter' / &
         ANTags(16)%Label / 'instr_1_tube_length' / &
         ANTags(17)%Label / 'instr_1_tube_flowrate' / &
         ANTags(18)%Label / 'instr_1_hpath_length' / &
         ANTags(19)%Label / 'instr_1_vpath_length' / &
         ANTags(20)%Label / 'instr_1_tau' / &
         ANTags(21)%Label / 'instr_1_kw' / &
         ANTags(22)%Label / 'instr_1_ko' / &
         ANTags(23)%Label / 'instr_1_void1' / &
         ANTags(24)%Label / 'instr_1_void2' / &
         ANTags(25)%Label / 'instr_2_height' / &
         ANTags(26)%Label / 'instr_2_north_offset' / &
         ANTags(27)%Label / 'instr_2_northward_separation' / &
         ANTags(28)%Label / 'instr_2_eastward_separation' / &
         ANTags(29)%Label / 'instr_2_vertical_separation' / &
         ANTags(30)%Label / 'instr_2_tube_diameter' / &
         ANTags(31)%Label / 'instr_2_tube_length' / &
         ANTags(32)%Label / 'instr_2_tube_flowrate' / &
         ANTags(33)%Label / 'instr_2_hpath_length' / &
         ANTags(34)%Label / 'instr_2_vpath_length' / &
         ANTags(35)%Label / 'instr_2_tau' / &
         ANTags(36)%Label / 'instr_2_kw' / &
         ANTags(37)%Label / 'instr_2_ko' / &
         ANTags(38)%Label / 'instr_2_void1' / &
         ANTags(39)%Label / 'instr_2_void2' / &
         ANTags(40)%Label / 'instr_3_height' / &
         ANTags(41)%Label / 'instr_3_north_offset' / &
         ANTags(42)%Label / 'instr_3_northward_separation' / &
         ANTags(43)%Label / 'instr_3_eastward_separation' / &
         ANTags(44)%Label / 'instr_3_vertical_separation' / &
         ANTags(45)%Label / 'instr_3_tube_diameter' / &
         ANTags(46)%Label / 'instr_3_tube_length' / &
         ANTags(47)%Label / 'instr_3_tube_flowrate' / &
         ANTags(48)%Label / 'instr_3_hpath_length' / &
         ANTags(49)%Label / 'instr_3_vpath_length' / &
         ANTags(50)%Label / 'instr_3_tau' / &
         ANTags(51)%Label / 'instr_3_kw' / &
         ANTags(52)%Label / 'instr_3_ko' / &
         ANTags(53)%Label / 'instr_3_void1' / &
         ANTags(54)%Label / 'instr_3_void2' / &
         ANTags(55)%Label / 'instr_4_height' / &
         ANTags(56)%Label / 'instr_4_north_offset' / &
         ANTags(57)%Label / 'instr_4_northward_separation' / &
         ANTags(58)%Label / 'instr_4_eastward_separation' / &
         ANTags(59)%Label / 'instr_4_vertical_separation' / &
         ANTags(60)%Label / 'instr_4_tube_diameter' / &
         ANTags(61)%Label / 'instr_4_tube_length' / &
         ANTags(62)%Label / 'instr_4_tube_flowrate' / &
         ANTags(63)%Label / 'instr_4_hpath_length' / &
         ANTags(64)%Label / 'instr_4_vpath_length' / &
         ANTags(65)%Label / 'instr_4_tau' / &
         ANTags(66)%Label / 'instr_4_kw' / &
         ANTags(67)%Label / 'instr_4_ko' / &
         ANTags(68)%Label / 'instr_4_void1' / &
         ANTags(69)%Label / 'instr_4_void2' / &
         ANTags(70)%Label / 'instr_5_height' / &
         ANTags(71)%Label / 'instr_5_north_offset' / &
         ANTags(72)%Label / 'instr_5_northward_separation' / &
         ANTags(73)%Label / 'instr_5_eastward_separation' / &
         ANTags(74)%Label / 'instr_5_vertical_separation' / &
         ANTags(75)%Label / 'instr_5_tube_diameter' / &
         ANTags(76)%Label / 'instr_5_tube_length' / &
         ANTags(77)%Label / 'instr_5_tube_flowrate' / &
         ANTags(78)%Label / 'instr_5_hpath_length' / &
         ANTags(79)%Label / 'instr_5_vpath_length' / &
         ANTags(80)%Label / 'instr_5_tau' / &
         ANTags(81)%Label / 'instr_5_kw' / &
         ANTags(82)%Label / 'instr_5_ko' / &
         ANTags(83)%Label / 'instr_5_void1' / &
         ANTags(84)%Label / 'instr_5_void2' / &
         ANTags(85)%Label / 'instr_6_height' / &
         ANTags(86)%Label / 'instr_6_north_offset' / &
         ANTags(87)%Label / 'instr_6_northward_separation' / &
         ANTags(88)%Label / 'instr_6_eastward_separation' / &
         ANTags(89)%Label / 'instr_6_vertical_separation' / &
         ANTags(90)%Label / 'instr_6_tube_diameter' / &
         ANTags(91)%Label / 'instr_6_tube_length' / &
         ANTags(92)%Label / 'instr_6_tube_flowrate' / &
         ANTags(93)%Label / 'instr_6_hpath_length' / &
         ANTags(94)%Label / 'instr_6_vpath_length' / &
         ANTags(95)%Label / 'instr_6_tau' / &
         ANTags(96)%Label / 'instr_6_kw' / &
         ANTags(97)%Label / 'instr_6_ko' / &
         ANTags(98)%Label / 'instr_6_void1' / &
         ANTags(99)%Label / 'instr_6_void2' / &
         ANTags(100)%Label / 'instr_7_height' / &
         ANTags(101)%Label / 'instr_7_north_offset' / &
         ANTags(102)%Label / 'instr_7_northward_separation' / &
         ANTags(103)%Label / 'instr_7_eastward_separation' / &
         ANTags(104)%Label / 'instr_7_vertical_separation' / &
         ANTags(105)%Label / 'instr_7_tube_diameter' / &
         ANTags(106)%Label / 'instr_7_tube_length' / &
         ANTags(107)%Label / 'instr_7_tube_flowrate' / &
         ANTags(108)%Label / 'instr_7_hpath_length' / &
         ANTags(109)%Label / 'instr_7_vpath_length' / &
         ANTags(110)%Label / 'instr_7_tau' / &
         ANTags(111)%Label / 'instr_7_kw' / &
         ANTags(112)%Label / 'instr_7_ko' / &
         ANTags(113)%Label / 'instr_7_void1' / &
         ANTags(114)%Label / 'instr_7_void2' / &
         ANTags(115)%Label / 'instr_8_height' / &
         ANTags(116)%Label / 'instr_8_north_offset' / &
         ANTags(117)%Label / 'instr_8_northward_separation' / &
         ANTags(118)%Label / 'instr_8_eastward_separation' / &
         ANTags(119)%Label / 'instr_8_vertical_separation' / &
         ANTags(120)%Label / 'instr_8_tube_diameter' / &
         ANTags(121)%Label / 'instr_8_tube_length' / &
         ANTags(122)%Label / 'instr_8_tube_flowrate' / &
         ANTags(123)%Label / 'instr_8_hpath_length' / &
         ANTags(124)%Label / 'instr_8_vpath_length' / &
         ANTags(125)%Label / 'instr_8_tau' / &
         ANTags(126)%Label / 'instr_8_kw' / &
         ANTags(127)%Label / 'instr_8_ko' / &
         ANTags(128)%Label / 'instr_8_void1' / &
         ANTags(129)%Label / 'instr_8_void2' / &
         ANTags(130)%Label / 'col_1_min_value' / &
         ANTags(131)%Label / 'col_1_max_value' / &
         ANTags(132)%Label / 'col_1_a_value' / &
         ANTags(133)%Label / 'col_1_b_value' / &
         ANTags(134)%Label / 'col_1_nom_timelag' / &
         ANTags(135)%Label / 'col_1_min_timelag' / &
         ANTags(136)%Label / 'col_1_max_timelag' / &
         ANTags(137)%Label / 'col_1_flag_threshold' / &
         ANTags(138)%Label / 'col_2_min_value' / &
         ANTags(139)%Label / 'col_2_max_value' / &
         ANTags(140)%Label / 'col_2_a_value' / &
         ANTags(141)%Label / 'col_2_b_value' / &
         ANTags(142)%Label / 'col_2_nom_timelag' / &
         ANTags(143)%Label / 'col_2_min_timelag' / &
         ANTags(144)%Label / 'col_2_max_timelag' / &
         ANTags(145)%Label / 'col_2_flag_threshold' / &
         ANTags(146)%Label / 'col_3_min_value' / &
         ANTags(147)%Label / 'col_3_max_value' / &
         ANTags(148)%Label / 'col_3_a_value' / &
         ANTags(149)%Label / 'col_3_b_value' / &
         ANTags(150)%Label / 'col_3_nom_timelag' / &
         ANTags(151)%Label / 'col_3_min_timelag' / &
         ANTags(152)%Label / 'col_3_max_timelag' / &
         ANTags(153)%Label / 'col_3_flag_threshold' / &
         ANTags(154)%Label / 'col_4_min_value' / &
         ANTags(155)%Label / 'col_4_max_value' / &
         ANTags(156)%Label / 'col_4_a_value' / &
         ANTags(157)%Label / 'col_4_b_value' / &
         ANTags(158)%Label / 'col_4_nom_timelag' / &
         ANTags(159)%Label / 'col_4_min_timelag' / &
         ANTags(160)%Label / 'col_4_max_timelag' / &
         ANTags(161)%Label / 'col_4_flag_threshold' / &
         ANTags(162)%Label / 'col_5_min_value' / &
         ANTags(163)%Label / 'col_5_max_value' / &
         ANTags(164)%Label / 'col_5_a_value' / &
         ANTags(165)%Label / 'col_5_b_value' / &
         ANTags(166)%Label / 'col_5_nom_timelag' / &
         ANTags(167)%Label / 'col_5_min_timelag' / &
         ANTags(168)%Label / 'col_5_max_timelag' / &
         ANTags(169)%Label / 'col_5_flag_threshold' / &
         ANTags(170)%Label / 'col_6_min_value' / &
         ANTags(171)%Label / 'col_6_max_value' / &
         ANTags(172)%Label / 'col_6_a_value' / &
         ANTags(173)%Label / 'col_6_b_value' / &
         ANTags(174)%Label / 'col_6_nom_timelag' / &
         ANTags(175)%Label / 'col_6_min_timelag' / &
         ANTags(176)%Label / 'col_6_max_timelag' / &
         ANTags(177)%Label / 'col_6_flag_threshold' / &
         ANTags(178)%Label / 'col_7_min_value' / &
         ANTags(179)%Label / 'col_7_max_value' / &
         ANTags(180)%Label / 'col_7_a_value' / &
         ANTags(181)%Label / 'col_7_b_value' / &
         ANTags(182)%Label / 'col_7_nom_timelag' / &
         ANTags(183)%Label / 'col_7_min_timelag' / &
         ANTags(184)%Label / 'col_7_max_timelag' / &
         ANTags(185)%Label / 'col_7_flag_threshold' / &
         ANTags(186)%Label / 'col_8_min_value' / &
         ANTags(187)%Label / 'col_8_max_value' / &
         ANTags(188)%Label / 'col_8_a_value' / &
         ANTags(189)%Label / 'col_8_b_value' / &
         ANTags(190)%Label / 'col_8_nom_timelag' / &
         ANTags(191)%Label / 'col_8_min_timelag' / &
         ANTags(192)%Label / 'col_8_max_timelag' / &
         ANTags(193)%Label / 'col_8_flag_threshold' / &
         ANTags(194)%Label / 'col_9_min_value' / &
         ANTags(195)%Label / 'col_9_max_value' / &
         ANTags(196)%Label / 'col_9_a_value' / &
         ANTags(197)%Label / 'col_9_b_value' / &
         ANTags(198)%Label / 'col_9_nom_timelag' / &
         ANTags(199)%Label / 'col_9_min_timelag' / &
         ANTags(200)%Label / 'col_9_max_timelag' /
    data ANTags(201)%Label / 'col_9_flag_threshold' / &
         ANTags(202)%Label / 'col_10_min_value' / &
         ANTags(203)%Label / 'col_10_max_value' / &
         ANTags(204)%Label / 'col_10_a_value' / &
         ANTags(205)%Label / 'col_10_b_value' / &
         ANTags(206)%Label / 'col_10_nom_timelag' / &
         ANTags(207)%Label / 'col_10_min_timelag' / &
         ANTags(208)%Label / 'col_10_max_timelag' / &
         ANTags(209)%Label / 'col_10_flag_threshold' / &
         ANTags(210)%Label / 'col_11_min_value' / &
         ANTags(211)%Label / 'col_11_max_value' / &
         ANTags(212)%Label / 'col_11_a_value' / &
         ANTags(213)%Label / 'col_11_b_value' / &
         ANTags(214)%Label / 'col_11_nom_timelag' / &
         ANTags(215)%Label / 'col_11_min_timelag' / &
         ANTags(216)%Label / 'col_11_max_timelag' / &
         ANTags(217)%Label / 'col_11_flag_threshold' / &
         ANTags(218)%Label / 'col_12_min_value' / &
         ANTags(219)%Label / 'col_12_max_value' / &
         ANTags(220)%Label / 'col_12_a_value' / &
         ANTags(221)%Label / 'col_12_b_value' / &
         ANTags(222)%Label / 'col_12_nom_timelag' / &
         ANTags(223)%Label / 'col_12_min_timelag' / &
         ANTags(224)%Label / 'col_12_max_timelag' / &
         ANTags(225)%Label / 'col_12_flag_threshold' / &
         ANTags(226)%Label / 'col_13_min_value' / &
         ANTags(227)%Label / 'col_13_max_value' / &
         ANTags(228)%Label / 'col_13_a_value' / &
         ANTags(229)%Label / 'col_13_b_value' / &
         ANTags(230)%Label / 'col_13_nom_timelag' / &
         ANTags(231)%Label / 'col_13_min_timelag' / &
         ANTags(232)%Label / 'col_13_max_timelag' / &
         ANTags(233)%Label / 'col_13_flag_threshold' / &
         ANTags(234)%Label / 'col_14_min_value' / &
         ANTags(235)%Label / 'col_14_max_value' / &
         ANTags(236)%Label / 'col_14_a_value' / &
         ANTags(237)%Label / 'col_14_b_value' / &
         ANTags(238)%Label / 'col_14_nom_timelag' / &
         ANTags(239)%Label / 'col_14_min_timelag' / &
         ANTags(240)%Label / 'col_14_max_timelag' / &
         ANTags(241)%Label / 'col_14_flag_threshold' / &
         ANTags(242)%Label / 'col_15_min_value' / &
         ANTags(243)%Label / 'col_15_max_value' / &
         ANTags(244)%Label / 'col_15_a_value' / &
         ANTags(245)%Label / 'col_15_b_value' / &
         ANTags(246)%Label / 'col_15_nom_timelag' / &
         ANTags(247)%Label / 'col_15_min_timelag' / &
         ANTags(248)%Label / 'col_15_max_timelag' / &
         ANTags(249)%Label / 'col_15_flag_threshold' / &
         ANTags(250)%Label / 'col_16_min_value' / &
         ANTags(251)%Label / 'col_16_max_value' / &
         ANTags(252)%Label / 'col_16_a_value' / &
         ANTags(253)%Label / 'col_16_b_value' / &
         ANTags(254)%Label / 'col_16_nom_timelag' / &
         ANTags(255)%Label / 'col_16_min_timelag' / &
         ANTags(256)%Label / 'col_16_max_timelag' / &
         ANTags(257)%Label / 'col_16_flag_threshold' / &
         ANTags(258)%Label / 'col_17_min_value' / &
         ANTags(259)%Label / 'col_17_max_value' / &
         ANTags(260)%Label / 'col_17_a_value' / &
         ANTags(261)%Label / 'col_17_b_value' / &
         ANTags(262)%Label / 'col_17_nom_timelag' / &
         ANTags(263)%Label / 'col_17_min_timelag' / &
         ANTags(264)%Label / 'col_17_max_timelag' / &
         ANTags(265)%Label / 'col_17_flag_threshold' / &
         ANTags(266)%Label / 'col_18_min_value' / &
         ANTags(267)%Label / 'col_18_max_value' / &
         ANTags(268)%Label / 'col_18_a_value' / &
         ANTags(269)%Label / 'col_18_b_value' / &
         ANTags(270)%Label / 'col_18_nom_timelag' / &
         ANTags(271)%Label / 'col_18_min_timelag' / &
         ANTags(272)%Label / 'col_18_max_timelag' / &
         ANTags(273)%Label / 'col_18_flag_threshold' / &
         ANTags(274)%Label / 'col_19_min_value' / &
         ANTags(275)%Label / 'col_19_max_value' / &
         ANTags(276)%Label / 'col_19_a_value' / &
         ANTags(277)%Label / 'col_19_b_value' / &
         ANTags(278)%Label / 'col_19_nom_timelag' / &
         ANTags(279)%Label / 'col_19_min_timelag' / &
         ANTags(280)%Label / 'col_19_max_timelag' / &
         ANTags(281)%Label / 'col_19_flag_threshold' / &
         ANTags(282)%Label / 'col_20_min_value' / &
         ANTags(283)%Label / 'col_20_max_value' / &
         ANTags(284)%Label / 'col_20_a_value' / &
         ANTags(285)%Label / 'col_20_b_value' / &
         ANTags(286)%Label / 'col_20_nom_timelag' / &
         ANTags(287)%Label / 'col_20_min_timelag' / &
         ANTags(288)%Label / 'col_20_max_timelag' / &
         ANTags(289)%Label / 'col_20_flag_threshold' / &
         ANTags(290)%Label / 'col_21_min_value' / &
         ANTags(291)%Label / 'col_21_max_value' / &
         ANTags(292)%Label / 'col_21_a_value' / &
         ANTags(293)%Label / 'col_21_b_value' / &
         ANTags(294)%Label / 'col_21_nom_timelag' / &
         ANTags(295)%Label / 'col_21_min_timelag' / &
         ANTags(296)%Label / 'col_21_max_timelag' / &
         ANTags(297)%Label / 'col_21_flag_threshold' / &
         ANTags(298)%Label / 'col_22_min_value' / &
         ANTags(299)%Label / 'col_22_max_value' / &
         ANTags(300)%Label / 'col_22_a_value' / &
         ANTags(301)%Label / 'col_22_b_value' / &
         ANTags(302)%Label / 'col_22_nom_timelag' / &
         ANTags(303)%Label / 'col_22_min_timelag' / &
         ANTags(304)%Label / 'col_22_max_timelag' / &
         ANTags(305)%Label / 'col_22_flag_threshold' / &
         ANTags(306)%Label / 'col_23_min_value' / &
         ANTags(307)%Label / 'col_23_max_value' / &
         ANTags(308)%Label / 'col_23_a_value' / &
         ANTags(309)%Label / 'col_23_b_value' / &
         ANTags(310)%Label / 'col_23_nom_timelag' / &
         ANTags(311)%Label / 'col_23_min_timelag' / &
         ANTags(312)%Label / 'col_23_max_timelag' / &
         ANTags(313)%Label / 'col_23_flag_threshold' / &
         ANTags(314)%Label / 'col_24_min_value' / &
         ANTags(315)%Label / 'col_24_max_value' / &
         ANTags(316)%Label / 'col_24_a_value' / &
         ANTags(317)%Label / 'col_24_b_value' / &
         ANTags(318)%Label / 'col_24_nom_timelag' / &
         ANTags(319)%Label / 'col_24_min_timelag' / &
         ANTags(320)%Label / 'col_24_max_timelag' / &
         ANTags(321)%Label / 'col_24_flag_threshold' / &
         ANTags(322)%Label / 'col_25_min_value' / &
         ANTags(323)%Label / 'col_25_max_value' / &
         ANTags(324)%Label / 'col_25_a_value' / &
         ANTags(325)%Label / 'col_25_b_value' / &
         ANTags(326)%Label / 'col_25_nom_timelag' / &
         ANTags(327)%Label / 'col_25_min_timelag' / &
         ANTags(328)%Label / 'col_25_max_timelag' / &
         ANTags(329)%Label / 'col_25_flag_threshold' / &
         ANTags(330)%Label / 'col_26_min_value' / &
         ANTags(331)%Label / 'col_26_max_value' / &
         ANTags(332)%Label / 'col_26_a_value' / &
         ANTags(333)%Label / 'col_26_b_value' / &
         ANTags(334)%Label / 'col_26_nom_timelag' / &
         ANTags(335)%Label / 'col_26_min_timelag' / &
         ANTags(336)%Label / 'col_26_max_timelag' / &
         ANTags(337)%Label / 'col_26_flag_threshold' / &
         ANTags(338)%Label / 'col_27_min_value' / &
         ANTags(339)%Label / 'col_27_max_value' / &
         ANTags(340)%Label / 'col_27_a_value' / &
         ANTags(341)%Label / 'col_27_b_value' / &
         ANTags(342)%Label / 'col_27_nom_timelag' / &
         ANTags(343)%Label / 'col_27_min_timelag' / &
         ANTags(344)%Label / 'col_27_max_timelag' / &
         ANTags(345)%Label / 'col_27_flag_threshold' / &
         ANTags(346)%Label / 'col_28_min_value' / &
         ANTags(347)%Label / 'col_28_max_value' / &
         ANTags(348)%Label / 'col_28_a_value' / &
         ANTags(349)%Label / 'col_28_b_value' / &
         ANTags(350)%Label / 'col_28_nom_timelag' / &
         ANTags(351)%Label / 'col_28_min_timelag' / &
         ANTags(352)%Label / 'col_28_max_timelag' / &
         ANTags(353)%Label / 'col_28_flag_threshold' / &
         ANTags(354)%Label / 'col_29_min_value' / &
         ANTags(355)%Label / 'col_29_max_value' / &
         ANTags(356)%Label / 'col_29_a_value' / &
         ANTags(357)%Label / 'col_29_b_value' / &
         ANTags(358)%Label / 'col_29_nom_timelag' / &
         ANTags(359)%Label / 'col_29_min_timelag' / &
         ANTags(360)%Label / 'col_29_max_timelag' / &
         ANTags(361)%Label / 'col_29_flag_threshold' / &
         ANTags(362)%Label / 'col_30_min_value' / &
         ANTags(363)%Label / 'col_30_max_value' / &
         ANTags(364)%Label / 'col_30_a_value' / &
         ANTags(365)%Label / 'col_30_b_value' / &
         ANTags(366)%Label / 'col_30_nom_timelag' / &
         ANTags(367)%Label / 'col_30_min_timelag' / &
         ANTags(368)%Label / 'col_30_max_timelag' / &
         ANTags(369)%Label / 'col_30_flag_threshold' / &
         ANTags(370)%Label / 'col_31_min_value' / &
         ANTags(371)%Label / 'col_31_max_value' / &
         ANTags(372)%Label / 'col_31_a_value' / &
         ANTags(373)%Label / 'col_31_b_value' / &
         ANTags(374)%Label / 'col_31_nom_timelag' / &
         ANTags(375)%Label / 'col_31_min_timelag' / &
         ANTags(376)%Label / 'col_31_max_timelag' / &
         ANTags(377)%Label / 'col_31_flag_threshold' / &
         ANTags(378)%Label / 'col_32_min_value' / &
         ANTags(379)%Label / 'col_32_max_value' / &
         ANTags(380)%Label / 'col_32_a_value' / &
         ANTags(381)%Label / 'col_32_b_value' / &
         ANTags(382)%Label / 'col_32_nom_timelag' / &
         ANTags(383)%Label / 'col_32_min_timelag' / &
         ANTags(384)%Label / 'col_32_max_timelag' / &
         ANTags(385)%Label / 'col_32_flag_threshold' / &
         ANTags(386)%Label / 'col_33_min_value' / &
         ANTags(387)%Label / 'col_33_max_value' / &
         ANTags(388)%Label / 'col_33_a_value' / &
         ANTags(389)%Label / 'col_33_b_value' / &
         ANTags(390)%Label / 'col_33_nom_timelag' / &
         ANTags(391)%Label / 'col_33_min_timelag' / &
         ANTags(392)%Label / 'col_33_max_timelag' / &
         ANTags(393)%Label / 'col_33_flag_threshold' / &
         ANTags(394)%Label / 'col_34_min_value' / &
         ANTags(395)%Label / 'col_34_max_value' / &
         ANTags(396)%Label / 'col_34_a_value' / &
         ANTags(397)%Label / 'col_34_b_value' / &
         ANTags(398)%Label / 'col_34_nom_timelag' / &
         ANTags(399)%Label / 'col_34_min_timelag' / &
         ANTags(400)%Label / 'col_34_max_timelag' /
    data ANTags(401)%Label / 'col_34_flag_threshold' / &
         ANTags(402)%Label / 'col_35_min_value' / &
         ANTags(403)%Label / 'col_35_max_value' / &
         ANTags(404)%Label / 'col_35_a_value' / &
         ANTags(405)%Label / 'col_35_b_value' / &
         ANTags(406)%Label / 'col_35_nom_timelag' / &
         ANTags(407)%Label / 'col_35_min_timelag' / &
         ANTags(408)%Label / 'col_35_max_timelag' / &
         ANTags(409)%Label / 'col_35_flag_threshold' / &
         ANTags(410)%Label / 'col_36_min_value' / &
         ANTags(411)%Label / 'col_36_max_value' / &
         ANTags(412)%Label / 'col_36_a_value' / &
         ANTags(413)%Label / 'col_36_b_value' / &
         ANTags(414)%Label / 'col_36_nom_timelag' / &
         ANTags(415)%Label / 'col_36_min_timelag' / &
         ANTags(416)%Label / 'col_36_max_timelag' / &
         ANTags(417)%Label / 'col_36_flag_threshold' / &
         ANTags(418)%Label / 'col_37_min_value' / &
         ANTags(419)%Label / 'col_37_max_value' / &
         ANTags(420)%Label / 'col_37_a_value' / &
         ANTags(421)%Label / 'col_37_b_value' / &
         ANTags(422)%Label / 'col_37_nom_timelag' / &
         ANTags(423)%Label / 'col_37_min_timelag' / &
         ANTags(424)%Label / 'col_37_max_timelag' / &
         ANTags(425)%Label / 'col_37_flag_threshold' / &
         ANTags(426)%Label / 'col_38_min_value' / &
         ANTags(427)%Label / 'col_38_max_value' / &
         ANTags(428)%Label / 'col_38_a_value' / &
         ANTags(429)%Label / 'col_38_b_value' / &
         ANTags(430)%Label / 'col_38_nom_timelag' / &
         ANTags(431)%Label / 'col_38_min_timelag' / &
         ANTags(432)%Label / 'col_38_max_timelag' / &
         ANTags(433)%Label / 'col_38_flag_threshold' / &
         ANTags(434)%Label / 'col_39_min_value' / &
         ANTags(435)%Label / 'col_39_max_value' / &
         ANTags(436)%Label / 'col_39_a_value' / &
         ANTags(437)%Label / 'col_39_b_value' / &
         ANTags(438)%Label / 'col_39_nom_timelag' / &
         ANTags(439)%Label / 'col_39_min_timelag' / &
         ANTags(440)%Label / 'col_39_max_timelag' / &
         ANTags(441)%Label / 'col_39_flag_threshold' / &
         ANTags(442)%Label / 'col_40_min_value' / &
         ANTags(443)%Label / 'col_40_max_value' / &
         ANTags(444)%Label / 'col_40_a_value' / &
         ANTags(445)%Label / 'col_40_b_value' / &
         ANTags(446)%Label / 'col_40_nom_timelag' / &
         ANTags(447)%Label / 'col_40_min_timelag' / &
         ANTags(448)%Label / 'col_40_max_timelag' / &
         ANTags(449)%Label / 'col_40_flag_threshold' / &
         ANTags(450)%Label / 'col_41_min_value' / &
         ANTags(451)%Label / 'col_41_max_value' / &
         ANTags(452)%Label / 'col_41_a_value' / &
         ANTags(453)%Label / 'col_41_b_value' / &
         ANTags(454)%Label / 'col_41_nom_timelag' / &
         ANTags(455)%Label / 'col_41_min_timelag' / &
         ANTags(456)%Label / 'col_41_max_timelag' / &
         ANTags(457)%Label / 'col_41_flag_threshold' / &
         ANTags(458)%Label / 'col_42_min_value' / &
         ANTags(459)%Label / 'col_42_max_value' / &
         ANTags(460)%Label / 'col_42_a_value' / &
         ANTags(461)%Label / 'col_42_b_value' / &
         ANTags(462)%Label / 'col_42_nom_timelag' / &
         ANTags(463)%Label / 'col_42_min_timelag' / &
         ANTags(464)%Label / 'col_42_max_timelag' / &
         ANTags(465)%Label / 'col_42_flag_threshold' / &
         ANTags(466)%Label / 'col_43_min_value' / &
         ANTags(467)%Label / 'col_43_max_value' / &
         ANTags(468)%Label / 'col_43_a_value' / &
         ANTags(469)%Label / 'col_43_b_value' / &
         ANTags(470)%Label / 'col_43_nom_timelag' / &
         ANTags(471)%Label / 'col_43_min_timelag' / &
         ANTags(472)%Label / 'col_43_max_timelag' / &
         ANTags(473)%Label / 'col_43_flag_threshold' / &
         ANTags(474)%Label / 'col_44_min_value' / &
         ANTags(475)%Label / 'col_44_max_value' / &
         ANTags(476)%Label / 'col_44_a_value' / &
         ANTags(477)%Label / 'col_44_b_value' / &
         ANTags(478)%Label / 'col_44_nom_timelag' / &
         ANTags(479)%Label / 'col_44_min_timelag' / &
         ANTags(480)%Label / 'col_44_max_timelag' / &
         ANTags(481)%Label / 'col_44_flag_threshold' / &
         ANTags(482)%Label / 'col_45_min_value' / &
         ANTags(483)%Label / 'col_45_max_value' / &
         ANTags(484)%Label / 'col_45_a_value' / &
         ANTags(485)%Label / 'col_45_b_value' / &
         ANTags(486)%Label / 'col_45_nom_timelag' / &
         ANTags(487)%Label / 'col_45_min_timelag' / &
         ANTags(488)%Label / 'col_45_max_timelag' / &
         ANTags(489)%Label / 'col_45_flag_threshold' / &
         ANTags(490)%Label / 'col_46_min_value' / &
         ANTags(491)%Label / 'col_46_max_value' / &
         ANTags(492)%Label / 'col_46_a_value' / &
         ANTags(493)%Label / 'col_46_b_value' / &
         ANTags(494)%Label / 'col_46_nom_timelag' / &
         ANTags(495)%Label / 'col_46_min_timelag' / &
         ANTags(496)%Label / 'col_46_max_timelag' / &
         ANTags(497)%Label / 'col_46_flag_threshold' / &
         ANTags(498)%Label / 'col_47_min_value' / &
         ANTags(499)%Label / 'col_47_max_value' / &
         ANTags(500)%Label / 'col_47_a_value' / &
         ANTags(501)%Label / 'col_47_b_value' / &
         ANTags(502)%Label / 'col_47_nom_timelag' / &
         ANTags(503)%Label / 'col_47_min_timelag' / &
         ANTags(504)%Label / 'col_47_max_timelag' / &
         ANTags(505)%Label / 'col_47_flag_threshold' / &
         ANTags(506)%Label / 'col_48_min_value' / &
         ANTags(507)%Label / 'col_48_max_value' / &
         ANTags(508)%Label / 'col_48_a_value' / &
         ANTags(509)%Label / 'col_48_b_value' / &
         ANTags(510)%Label / 'col_48_nom_timelag' / &
         ANTags(511)%Label / 'col_48_min_timelag' / &
         ANTags(512)%Label / 'col_48_max_timelag' / &
         ANTags(513)%Label / 'col_48_flag_threshold' / &
         ANTags(514)%Label / 'col_49_min_value' / &
         ANTags(515)%Label / 'col_49_max_value' / &
         ANTags(516)%Label / 'col_49_a_value' / &
         ANTags(517)%Label / 'col_49_b_value' / &
         ANTags(518)%Label / 'col_49_nom_timelag' / &
         ANTags(519)%Label / 'col_49_min_timelag' / &
         ANTags(520)%Label / 'col_49_max_timelag' / &
         ANTags(521)%Label / 'col_49_flag_threshold' / &
         ANTags(522)%Label / 'col_50_min_value' / &
         ANTags(523)%Label / 'col_50_max_value' / &
         ANTags(524)%Label / 'col_50_a_value' / &
         ANTags(525)%Label / 'col_50_b_value' / &
         ANTags(526)%Label / 'col_50_nom_timelag' / &
         ANTags(527)%Label / 'col_50_min_timelag' / &
         ANTags(528)%Label / 'col_50_max_timelag' / &
         ANTags(529)%Label / 'col_50_flag_threshold' / &
         ANTags(530)%Label / 'col_51_min_value' / &
         ANTags(531)%Label / 'col_51_max_value' / &
         ANTags(532)%Label / 'col_51_a_value' / &
         ANTags(533)%Label / 'col_51_b_value' / &
         ANTags(534)%Label / 'col_51_nom_timelag' / &
         ANTags(535)%Label / 'col_51_min_timelag' / &
         ANTags(536)%Label / 'col_51_max_timelag' / &
         ANTags(537)%Label / 'col_51_flag_threshold' / &
         ANTags(538)%Label / 'col_52_min_value' / &
         ANTags(539)%Label / 'col_52_max_value' / &
         ANTags(540)%Label / 'col_52_a_value' / &
         ANTags(541)%Label / 'col_52_b_value' / &
         ANTags(542)%Label / 'col_52_nom_timelag' / &
         ANTags(543)%Label / 'col_52_min_timelag' / &
         ANTags(544)%Label / 'col_52_max_timelag' / &
         ANTags(545)%Label / 'col_52_flag_threshold' / &
         ANTags(546)%Label / 'col_53_min_value' / &
         ANTags(547)%Label / 'col_53_max_value' / &
         ANTags(548)%Label / 'col_53_a_value' / &
         ANTags(549)%Label / 'col_53_b_value' / &
         ANTags(550)%Label / 'col_53_nom_timelag' / &
         ANTags(551)%Label / 'col_53_min_timelag' / &
         ANTags(552)%Label / 'col_53_max_timelag' / &
         ANTags(553)%Label / 'col_53_flag_threshold' / &
         ANTags(554)%Label / 'col_54_min_value' / &
         ANTags(555)%Label / 'col_54_max_value' / &
         ANTags(556)%Label / 'col_54_a_value' / &
         ANTags(557)%Label / 'col_54_b_value' / &
         ANTags(558)%Label / 'col_54_nom_timelag' / &
         ANTags(559)%Label / 'col_54_min_timelag' / &
         ANTags(560)%Label / 'col_54_max_timelag' / &
         ANTags(561)%Label / 'col_54_flag_threshold' / &
         ANTags(562)%Label / 'col_55_min_value' / &
         ANTags(563)%Label / 'col_55_max_value' / &
         ANTags(564)%Label / 'col_55_a_value' / &
         ANTags(565)%Label / 'col_55_b_value' / &
         ANTags(566)%Label / 'col_55_nom_timelag' / &
         ANTags(567)%Label / 'col_55_min_timelag' / &
         ANTags(568)%Label / 'col_55_max_timelag' / &
         ANTags(569)%Label / 'col_55_flag_threshold' / &
         ANTags(570)%Label / 'col_56_min_value' / &
         ANTags(571)%Label / 'col_56_max_value' / &
         ANTags(572)%Label / 'col_56_a_value' / &
         ANTags(573)%Label / 'col_56_b_value' / &
         ANTags(574)%Label / 'col_56_nom_timelag' / &
         ANTags(575)%Label / 'col_56_min_timelag' / &
         ANTags(576)%Label / 'col_56_max_timelag' / &
         ANTags(577)%Label / 'col_56_flag_threshold' / &
         ANTags(578)%Label / 'col_57_min_value' / &
         ANTags(579)%Label / 'col_57_max_value' / &
         ANTags(580)%Label / 'col_57_a_value' / &
         ANTags(581)%Label / 'col_57_b_value' / &
         ANTags(582)%Label / 'col_57_nom_timelag' / &
         ANTags(583)%Label / 'col_57_min_timelag' / &
         ANTags(584)%Label / 'col_57_max_timelag' / &
         ANTags(585)%Label / 'col_57_flag_threshold' / &
         ANTags(586)%Label / 'col_58_min_value' / &
         ANTags(587)%Label / 'col_58_max_value' / &
         ANTags(588)%Label / 'col_58_a_value' / &
         ANTags(589)%Label / 'col_58_b_value' / &
         ANTags(590)%Label / 'col_58_nom_timelag' / &
         ANTags(591)%Label / 'col_58_min_timelag' / &
         ANTags(592)%Label / 'col_58_max_timelag' / &
         ANTags(593)%Label / 'col_58_flag_threshold' / &
         ANTags(594)%Label / 'col_59_min_value' / &
         ANTags(595)%Label / 'col_59_max_value' / &
         ANTags(596)%Label / 'col_59_a_value' / &
         ANTags(597)%Label / 'col_59_b_value' / &
         ANTags(598)%Label / 'col_59_nom_timelag' / &
         ANTags(599)%Label / 'col_59_min_timelag' / &
         ANTags(600)%Label / 'col_59_max_timelag' /
    data ANTags(601)%Label / 'col_59_flag_threshold' / &
         ANTags(602)%Label / 'col_60_min_value' / &
         ANTags(603)%Label / 'col_60_max_value' / &
         ANTags(604)%Label / 'col_60_a_value' / &
         ANTags(605)%Label / 'col_60_b_value' / &
         ANTags(606)%Label / 'col_60_nom_timelag' / &
         ANTags(607)%Label / 'col_60_min_timelag' / &
         ANTags(608)%Label / 'col_60_max_timelag' / &
         ANTags(609)%Label / 'col_60_flag_threshold' / &
         ANTags(610)%Label / 'col_61_min_value' / &
         ANTags(611)%Label / 'col_61_max_value' / &
         ANTags(612)%Label / 'col_61_a_value' / &
         ANTags(613)%Label / 'col_61_b_value' / &
         ANTags(614)%Label / 'col_61_nom_timelag' / &
         ANTags(615)%Label / 'col_61_min_timelag' / &
         ANTags(616)%Label / 'col_61_max_timelag' / &
         ANTags(617)%Label / 'col_61_flag_threshold' / &
         ANTags(618)%Label / 'col_62_min_value' / &
         ANTags(619)%Label / 'col_62_max_value' / &
         ANTags(620)%Label / 'col_62_a_value' / &
         ANTags(621)%Label / 'col_62_b_value' / &
         ANTags(622)%Label / 'col_62_nom_timelag' / &
         ANTags(623)%Label / 'col_62_min_timelag' / &
         ANTags(624)%Label / 'col_62_max_timelag' / &
         ANTags(625)%Label / 'col_62_flag_threshold' / &
         ANTags(626)%Label / 'col_63_min_value' / &
         ANTags(627)%Label / 'col_63_max_value' / &
         ANTags(628)%Label / 'col_63_a_value' / &
         ANTags(629)%Label / 'col_63_b_value' / &
         ANTags(630)%Label / 'col_63_nom_timelag' / &
         ANTags(631)%Label / 'col_63_min_timelag' / &
         ANTags(632)%Label / 'col_63_max_timelag' / &
         ANTags(633)%Label / 'col_63_flag_threshold' / &
         ANTags(634)%Label / 'col_64_min_value' / &
         ANTags(635)%Label / 'col_64_max_value' / &
         ANTags(636)%Label / 'col_64_a_value' / &
         ANTags(637)%Label / 'col_64_b_value' / &
         ANTags(638)%Label / 'col_64_nom_timelag' / &
         ANTags(639)%Label / 'col_64_min_timelag' / &
         ANTags(640)%Label / 'col_64_max_timelag' / &
         ANTags(641)%Label / 'col_64_flag_threshold' / &
         ANTags(642)%Label / 'col_65_min_value' / &
         ANTags(643)%Label / 'col_65_max_value' / &
         ANTags(644)%Label / 'col_65_a_value' / &
         ANTags(645)%Label / 'col_65_b_value' / &
         ANTags(646)%Label / 'col_65_nom_timelag' / &
         ANTags(647)%Label / 'col_65_min_timelag' / &
         ANTags(648)%Label / 'col_65_max_timelag' / &
         ANTags(649)%Label / 'col_65_flag_threshold' / &
         ANTags(650)%Label / 'col_66_min_value' / &
         ANTags(651)%Label / 'col_66_max_value' / &
         ANTags(652)%Label / 'col_66_a_value' / &
         ANTags(653)%Label / 'col_66_b_value' / &
         ANTags(654)%Label / 'col_66_nom_timelag' / &
         ANTags(655)%Label / 'col_66_min_timelag' / &
         ANTags(656)%Label / 'col_66_max_timelag' / &
         ANTags(657)%Label / 'col_66_flag_threshold' / &
         ANTags(658)%Label / 'col_67_min_value' / &
         ANTags(659)%Label / 'col_67_max_value' / &
         ANTags(660)%Label / 'col_67_a_value' / &
         ANTags(661)%Label / 'col_67_b_value' / &
         ANTags(662)%Label / 'col_67_nom_timelag' / &
         ANTags(663)%Label / 'col_67_min_timelag' / &
         ANTags(664)%Label / 'col_67_max_timelag' / &
         ANTags(665)%Label / 'col_67_flag_threshold' / &
         ANTags(666)%Label / 'col_68_min_value' / &
         ANTags(667)%Label / 'col_68_max_value' / &
         ANTags(668)%Label / 'col_68_a_value' / &
         ANTags(669)%Label / 'col_68_b_value' / &
         ANTags(670)%Label / 'col_68_nom_timelag' / &
         ANTags(671)%Label / 'col_68_min_timelag' / &
         ANTags(672)%Label / 'col_68_max_timelag' / &
         ANTags(673)%Label / 'col_68_flag_threshold' / &
         ANTags(674)%Label / 'col_69_min_value' / &
         ANTags(675)%Label / 'col_69_max_value' / &
         ANTags(676)%Label / 'col_69_a_value' / &
         ANTags(677)%Label / 'col_69_b_value' / &
         ANTags(678)%Label / 'col_69_nom_timelag' / &
         ANTags(679)%Label / 'col_69_min_timelag' / &
         ANTags(680)%Label / 'col_69_max_timelag' / &
         ANTags(681)%Label / 'col_69_flag_threshold' / &
         ANTags(682)%Label / 'col_70_min_value' / &
         ANTags(683)%Label / 'col_70_max_value' / &
         ANTags(684)%Label / 'col_70_a_value' / &
         ANTags(685)%Label / 'col_70_b_value' / &
         ANTags(686)%Label / 'col_70_nom_timelag' / &
         ANTags(687)%Label / 'col_70_min_timelag' / &
         ANTags(688)%Label / 'col_70_max_timelag' / &
         ANTags(689)%Label / 'col_70_flag_threshold' / &
         ANTags(690)%Label / 'col_71_min_value' / &
         ANTags(691)%Label / 'col_71_max_value' / &
         ANTags(692)%Label / 'col_71_a_value' / &
         ANTags(693)%Label / 'col_71_b_value' / &
         ANTags(694)%Label / 'col_71_nom_timelag' / &
         ANTags(695)%Label / 'col_71_min_timelag' / &
         ANTags(696)%Label / 'col_71_max_timelag' / &
         ANTags(697)%Label / 'col_71_flag_threshold' / &
         ANTags(698)%Label / 'col_72_min_value' / &
         ANTags(699)%Label / 'col_72_max_value' / &
         ANTags(700)%Label / 'col_72_a_value' / &
         ANTags(701)%Label / 'col_72_b_value' / &
         ANTags(702)%Label / 'col_72_nom_timelag' / &
         ANTags(703)%Label / 'col_72_min_timelag' / &
         ANTags(704)%Label / 'col_72_max_timelag' / &
         ANTags(705)%Label / 'col_72_flag_threshold' / &
         ANTags(706)%Label / 'col_73_min_value' / &
         ANTags(707)%Label / 'col_73_max_value' / &
         ANTags(708)%Label / 'col_73_a_value' / &
         ANTags(709)%Label / 'col_73_b_value' / &
         ANTags(710)%Label / 'col_73_nom_timelag' / &
         ANTags(711)%Label / 'col_73_min_timelag' / &
         ANTags(712)%Label / 'col_73_max_timelag' / &
         ANTags(713)%Label / 'col_73_flag_threshold' / &
         ANTags(714)%Label / 'col_74_min_value' / &
         ANTags(715)%Label / 'col_74_max_value' / &
         ANTags(716)%Label / 'col_74_a_value' / &
         ANTags(717)%Label / 'col_74_b_value' / &
         ANTags(718)%Label / 'col_74_nom_timelag' / &
         ANTags(719)%Label / 'col_74_min_timelag' / &
         ANTags(720)%Label / 'col_74_max_timelag' / &
         ANTags(721)%Label / 'col_74_flag_threshold' / &
         ANTags(722)%Label / 'col_75_min_value' / &
         ANTags(723)%Label / 'col_75_max_value' / &
         ANTags(724)%Label / 'col_75_a_value' / &
         ANTags(725)%Label / 'col_75_b_value' / &
         ANTags(726)%Label / 'col_75_nom_timelag' / &
         ANTags(727)%Label / 'col_75_min_timelag' / &
         ANTags(728)%Label / 'col_75_max_timelag' / &
         ANTags(729)%Label / 'col_75_flag_threshold' / &
         ANTags(730)%Label / 'col_76_min_value' / &
         ANTags(731)%Label / 'col_76_max_value' / &
         ANTags(732)%Label / 'col_76_a_value' / &
         ANTags(733)%Label / 'col_76_b_value' / &
         ANTags(734)%Label / 'col_76_nom_timelag' / &
         ANTags(735)%Label / 'col_76_min_timelag' / &
         ANTags(736)%Label / 'col_76_max_timelag' / &
         ANTags(737)%Label / 'col_76_flag_threshold' / &
         ANTags(738)%Label / 'col_77_min_value' / &
         ANTags(739)%Label / 'col_77_max_value' / &
         ANTags(740)%Label / 'col_77_a_value' / &
         ANTags(741)%Label / 'col_77_b_value' / &
         ANTags(742)%Label / 'col_77_nom_timelag' / &
         ANTags(743)%Label / 'col_77_min_timelag' / &
         ANTags(744)%Label / 'col_77_max_timelag' / &
         ANTags(745)%Label / 'col_77_flag_threshold' / &
         ANTags(746)%Label / 'col_78_min_value' / &
         ANTags(747)%Label / 'col_78_max_value' / &
         ANTags(748)%Label / 'col_78_a_value' / &
         ANTags(749)%Label / 'col_78_b_value' / &
         ANTags(750)%Label / 'col_78_nom_timelag' / &
         ANTags(751)%Label / 'col_78_min_timelag' / &
         ANTags(752)%Label / 'col_78_max_timelag' / &
         ANTags(753)%Label / 'col_78_flag_threshold' / &
         ANTags(754)%Label / 'col_79_min_value' / &
         ANTags(755)%Label / 'col_79_max_value' / &
         ANTags(756)%Label / 'col_79_a_value' / &
         ANTags(757)%Label / 'col_79_b_value' / &
         ANTags(758)%Label / 'col_79_nom_timelag' / &
         ANTags(759)%Label / 'col_79_min_timelag' / &
         ANTags(760)%Label / 'col_79_max_timelag' / &
         ANTags(761)%Label / 'col_79_flag_threshold' / &
         ANTags(762)%Label / 'col_80_min_value' / &
         ANTags(763)%Label / 'col_80_max_value' / &
         ANTags(764)%Label / 'col_80_a_value' / &
         ANTags(765)%Label / 'col_80_b_value' / &
         ANTags(766)%Label / 'col_80_nom_timelag' / &
         ANTags(767)%Label / 'col_80_min_timelag' / &
         ANTags(768)%Label / 'col_80_max_timelag' / &
         ANTags(769)%Label / 'col_80_flag_threshold' / &
         ANTags(770)%Label / 'col_81_min_value' / &
         ANTags(771)%Label / 'col_81_max_value' / &
         ANTags(772)%Label / 'col_81_a_value' / &
         ANTags(773)%Label / 'col_81_b_value' / &
         ANTags(774)%Label / 'col_81_nom_timelag' / &
         ANTags(775)%Label / 'col_81_min_timelag' / &
         ANTags(776)%Label / 'col_81_max_timelag' / &
         ANTags(777)%Label / 'col_81_flag_threshold' / &
         ANTags(778)%Label / 'col_82_min_value' / &
         ANTags(779)%Label / 'col_82_max_value' / &
         ANTags(780)%Label / 'col_82_a_value' / &
         ANTags(781)%Label / 'col_82_b_value' / &
         ANTags(782)%Label / 'col_82_nom_timelag' / &
         ANTags(783)%Label / 'col_82_min_timelag' / &
         ANTags(784)%Label / 'col_82_max_timelag' / &
         ANTags(785)%Label / 'col_82_flag_threshold' / &
         ANTags(786)%Label / 'col_83_min_value' / &
         ANTags(787)%Label / 'col_83_max_value' / &
         ANTags(788)%Label / 'col_83_a_value' / &
         ANTags(789)%Label / 'col_83_b_value' / &
         ANTags(790)%Label / 'col_83_nom_timelag' / &
         ANTags(791)%Label / 'col_83_min_timelag' / &
         ANTags(792)%Label / 'col_83_max_timelag' / &
         ANTags(793)%Label / 'col_83_flag_threshold' / &
         ANTags(794)%Label / 'col_84_min_value' / &
         ANTags(795)%Label / 'col_84_max_value' / &
         ANTags(796)%Label / 'col_84_a_value' / &
         ANTags(797)%Label / 'col_84_b_value' / &
         ANTags(798)%Label / 'col_84_nom_timelag' / &
         ANTags(799)%Label / 'col_84_min_timelag' / &
         ANTags(800)%Label / 'col_84_max_timelag' /
    data ANTags(801)%Label / 'col_84_flag_threshold' / &
         ANTags(802)%Label / 'col_85_min_value' / &
         ANTags(803)%Label / 'col_85_max_value' / &
         ANTags(804)%Label / 'col_85_a_value' / &
         ANTags(805)%Label / 'col_85_b_value' / &
         ANTags(806)%Label / 'col_85_nom_timelag' / &
         ANTags(807)%Label / 'col_85_min_timelag' / &
         ANTags(808)%Label / 'col_85_max_timelag' / &
         ANTags(809)%Label / 'col_85_flag_threshold' / &
         ANTags(810)%Label / 'col_86_min_value' / &
         ANTags(811)%Label / 'col_86_max_value' / &
         ANTags(812)%Label / 'col_86_a_value' / &
         ANTags(813)%Label / 'col_86_b_value' / &
         ANTags(814)%Label / 'col_86_nom_timelag' / &
         ANTags(815)%Label / 'col_86_min_timelag' / &
         ANTags(816)%Label / 'col_86_max_timelag' / &
         ANTags(817)%Label / 'col_86_flag_threshold' / &
         ANTags(818)%Label / 'col_87_min_value' / &
         ANTags(819)%Label / 'col_87_max_value' / &
         ANTags(820)%Label / 'col_87_a_value' / &
         ANTags(821)%Label / 'col_87_b_value' / &
         ANTags(822)%Label / 'col_87_nom_timelag' / &
         ANTags(823)%Label / 'col_87_min_timelag' / &
         ANTags(824)%Label / 'col_87_max_timelag' / &
         ANTags(825)%Label / 'col_87_flag_threshold' / &
         ANTags(826)%Label / 'col_88_min_value' / &
         ANTags(827)%Label / 'col_88_max_value' / &
         ANTags(828)%Label / 'col_88_a_value' / &
         ANTags(829)%Label / 'col_88_b_value' / &
         ANTags(830)%Label / 'col_88_nom_timelag' / &
         ANTags(831)%Label / 'col_88_min_timelag' / &
         ANTags(832)%Label / 'col_88_max_timelag' / &
         ANTags(833)%Label / 'col_88_flag_threshold' / &
         ANTags(834)%Label / 'col_89_min_value' / &
         ANTags(835)%Label / 'col_89_max_value' / &
         ANTags(836)%Label / 'col_89_a_value' / &
         ANTags(837)%Label / 'col_89_b_value' / &
         ANTags(838)%Label / 'col_89_nom_timelag' / &
         ANTags(839)%Label / 'col_89_min_timelag' / &
         ANTags(840)%Label / 'col_89_max_timelag' / &
         ANTags(841)%Label / 'col_89_flag_threshold' / &
         ANTags(842)%Label / 'col_90_min_value' / &
         ANTags(843)%Label / 'col_90_max_value' / &
         ANTags(844)%Label / 'col_90_a_value' / &
         ANTags(845)%Label / 'col_90_b_value' / &
         ANTags(846)%Label / 'col_90_nom_timelag' / &
         ANTags(847)%Label / 'col_90_min_timelag' / &
         ANTags(848)%Label / 'col_90_max_timelag' / &
         ANTags(849)%Label / 'col_90_flag_threshold' / &
         ANTags(850)%Label / 'col_91_min_value' / &
         ANTags(851)%Label / 'col_91_max_value' / &
         ANTags(852)%Label / 'col_91_a_value' / &
         ANTags(853)%Label / 'col_91_b_value' / &
         ANTags(854)%Label / 'col_91_nom_timelag' / &
         ANTags(855)%Label / 'col_91_min_timelag' / &
         ANTags(856)%Label / 'col_91_max_timelag' / &
         ANTags(857)%Label / 'col_91_flag_threshold' / &
         ANTags(858)%Label / 'col_92_min_value' / &
         ANTags(859)%Label / 'col_92_max_value' / &
         ANTags(860)%Label / 'col_92_a_value' / &
         ANTags(861)%Label / 'col_92_b_value' / &
         ANTags(862)%Label / 'col_92_nom_timelag' / &
         ANTags(863)%Label / 'col_92_min_timelag' / &
         ANTags(864)%Label / 'col_92_max_timelag' / &
         ANTags(865)%Label / 'col_92_flag_threshold' / &
         ANTags(866)%Label / 'col_93_min_value' / &
         ANTags(867)%Label / 'col_93_max_value' / &
         ANTags(868)%Label / 'col_93_a_value' / &
         ANTags(869)%Label / 'col_93_b_value' / &
         ANTags(870)%Label / 'col_93_nom_timelag' / &
         ANTags(871)%Label / 'col_93_min_timelag' / &
         ANTags(872)%Label / 'col_93_max_timelag' / &
         ANTags(873)%Label / 'col_93_flag_threshold' / &
         ANTags(874)%Label / 'col_94_min_value' / &
         ANTags(875)%Label / 'col_94_max_value' / &
         ANTags(876)%Label / 'col_94_a_value' / &
         ANTags(877)%Label / 'col_94_b_value' / &
         ANTags(878)%Label / 'col_94_nom_timelag' / &
         ANTags(879)%Label / 'col_94_min_timelag' / &
         ANTags(880)%Label / 'col_94_max_timelag' / &
         ANTags(881)%Label / 'col_94_flag_threshold' / &
         ANTags(882)%Label / 'col_95_min_value' / &
         ANTags(883)%Label / 'col_95_max_value' / &
         ANTags(884)%Label / 'col_95_a_value' / &
         ANTags(885)%Label / 'col_95_b_value' / &
         ANTags(886)%Label / 'col_95_nom_timelag' / &
         ANTags(887)%Label / 'col_95_min_timelag' / &
         ANTags(888)%Label / 'col_95_max_timelag' / &
         ANTags(889)%Label / 'col_95_flag_threshold' / &
         ANTags(890)%Label / 'col_96_min_value' / &
         ANTags(891)%Label / 'col_96_max_value' / &
         ANTags(892)%Label / 'col_96_a_value' / &
         ANTags(893)%Label / 'col_96_b_value' / &
         ANTags(894)%Label / 'col_96_nom_timelag' / &
         ANTags(895)%Label / 'col_96_min_timelag' / &
         ANTags(896)%Label / 'col_96_max_timelag' / &
         ANTags(897)%Label / 'col_96_flag_threshold' / &
         ANTags(898)%Label / 'col_97_min_value' / &
         ANTags(899)%Label / 'col_97_max_value' / &
         ANTags(900)%Label / 'col_97_a_value' / &
         ANTags(901)%Label / 'col_97_b_value' / &
         ANTags(902)%Label / 'col_97_nom_timelag' / &
         ANTags(903)%Label / 'col_97_min_timelag' / &
         ANTags(904)%Label / 'col_97_max_timelag' / &
         ANTags(905)%Label / 'col_97_flag_threshold' / &
         ANTags(906)%Label / 'col_98_min_value' / &
         ANTags(907)%Label / 'col_98_max_value' / &
         ANTags(908)%Label / 'col_98_a_value' / &
         ANTags(909)%Label / 'col_98_b_value' / &
         ANTags(910)%Label / 'col_98_nom_timelag' / &
         ANTags(911)%Label / 'col_98_min_timelag' / &
         ANTags(912)%Label / 'col_98_max_timelag' / &
         ANTags(913)%Label / 'col_98_flag_threshold' / &
         ANTags(914)%Label / 'col_99_min_value' / &
         ANTags(915)%Label / 'col_99_max_value' / &
         ANTags(916)%Label / 'col_99_a_value' / &
         ANTags(917)%Label / 'col_99_b_value' / &
         ANTags(918)%Label / 'col_99_nom_timelag' / &
         ANTags(919)%Label / 'col_99_min_timelag' / &
         ANTags(920)%Label / 'col_99_max_timelag' / &
         ANTags(921)%Label / 'col_99_flag_threshold' / &
         ANTags(922)%Label / 'col_100_min_value' / &
         ANTags(923)%Label / 'col_100_max_value' / &
         ANTags(924)%Label / 'col_100_a_value' / &
         ANTags(925)%Label / 'col_100_b_value' / &
         ANTags(926)%Label / 'col_100_nom_timelag' / &
         ANTags(927)%Label / 'col_100_min_timelag' / &
         ANTags(928)%Label / 'col_100_max_timelag' / &
         ANTags(929)%Label / 'col_100_flag_threshold' /
    !> END GENERATED ANTags

    !> BEGIN GENERATED ACTags - edit gen_metadata_tags.py, not this block
    data ACTags(1)%Label / 'logger_sw_version' / &
         ACTags(2)%Label / 'title' / &
         ACTags(3)%Label / 'creation_date' / &
         ACTags(4)%Label / 'start_date' / &
         ACTags(5)%Label / 'end_date' / &
         ACTags(6)%Label / 'file_name' / &
         ACTags(7)%Label / 'data_path' / &
         ACTags(8)%Label / 'project_notes' / &
         ACTags(9)%Label / 'site_name' / &
         ACTags(10)%Label / 'site_id' / &
         ACTags(11)%Label / 'site_notes' / &
         ACTags(12)%Label / 'station_name' / &
         ACTags(13)%Label / 'station_id' / &
         ACTags(14)%Label / 'pc_time_settings' / &
         ACTags(15)%Label / 'timing_notes' / &
         ACTags(16)%Label / 'saved_native' / &
         ACTags(17)%Label / 'timestamp' / &
         ACTags(18)%Label / 'enable_processing' / &
         ACTags(19)%Label / 'iso_format' / &
         ACTags(20)%Label / 'tstamp_end' / &
         ACTags(21)%Label / 'native_format' / &
         ACTags(22)%Label / 'head_corr' / &
         ACTags(23)%Label / 'separator' / &
         ACTags(24)%Label / 'flag_discards_if_above' / &
         ACTags(25)%Label / 'instr_1_manufacturer' / &
         ACTags(26)%Label / 'instr_1_sw_version' / &
         ACTags(27)%Label / 'instr_1_model' / &
         ACTags(28)%Label / 'instr_1_sn' / &
         ACTags(29)%Label / 'instr_1_id' / &
         ACTags(30)%Label / 'instr_1_wformat' / &
         ACTags(31)%Label / 'instr_1_wref' / &
         ACTags(32)%Label / 'instr_1_head_corr' / &
         ACTags(33)%Label / 'instr_2_manufacturer' / &
         ACTags(34)%Label / 'instr_2_sw_version' / &
         ACTags(35)%Label / 'instr_2_model' / &
         ACTags(36)%Label / 'instr_2_sn' / &
         ACTags(37)%Label / 'instr_2_id' / &
         ACTags(38)%Label / 'instr_2_wformat' / &
         ACTags(39)%Label / 'instr_2_wref' / &
         ACTags(40)%Label / 'instr_2_head_corr' / &
         ACTags(41)%Label / 'instr_3_manufacturer' / &
         ACTags(42)%Label / 'instr_3_sw_version' / &
         ACTags(43)%Label / 'instr_3_model' / &
         ACTags(44)%Label / 'instr_3_sn' / &
         ACTags(45)%Label / 'instr_3_id' / &
         ACTags(46)%Label / 'instr_3_wformat' / &
         ACTags(47)%Label / 'instr_3_wref' / &
         ACTags(48)%Label / 'instr_3_head_corr' / &
         ACTags(49)%Label / 'instr_4_manufacturer' / &
         ACTags(50)%Label / 'instr_4_sw_version' / &
         ACTags(51)%Label / 'instr_4_model' / &
         ACTags(52)%Label / 'instr_4_sn' / &
         ACTags(53)%Label / 'instr_4_id' / &
         ACTags(54)%Label / 'instr_4_wformat' / &
         ACTags(55)%Label / 'instr_4_wref' / &
         ACTags(56)%Label / 'instr_4_head_corr' / &
         ACTags(57)%Label / 'instr_5_manufacturer' / &
         ACTags(58)%Label / 'instr_5_sw_version' / &
         ACTags(59)%Label / 'instr_5_model' / &
         ACTags(60)%Label / 'instr_5_sn' / &
         ACTags(61)%Label / 'instr_5_id' / &
         ACTags(62)%Label / 'instr_5_wformat' / &
         ACTags(63)%Label / 'instr_5_wref' / &
         ACTags(64)%Label / 'instr_5_head_corr' / &
         ACTags(65)%Label / 'instr_6_manufacturer' / &
         ACTags(66)%Label / 'instr_6_sw_version' / &
         ACTags(67)%Label / 'instr_6_model' / &
         ACTags(68)%Label / 'instr_6_sn' / &
         ACTags(69)%Label / 'instr_6_id' / &
         ACTags(70)%Label / 'instr_6_wformat' / &
         ACTags(71)%Label / 'instr_6_wref' / &
         ACTags(72)%Label / 'instr_6_head_corr' / &
         ACTags(73)%Label / 'instr_7_manufacturer' / &
         ACTags(74)%Label / 'instr_7_sw_version' / &
         ACTags(75)%Label / 'instr_7_model' / &
         ACTags(76)%Label / 'instr_7_sn' / &
         ACTags(77)%Label / 'instr_7_id' / &
         ACTags(78)%Label / 'instr_7_wformat' / &
         ACTags(79)%Label / 'instr_7_wref' / &
         ACTags(80)%Label / 'instr_7_head_corr' / &
         ACTags(81)%Label / 'instr_8_manufacturer' / &
         ACTags(82)%Label / 'instr_8_sw_version' / &
         ACTags(83)%Label / 'instr_8_model' / &
         ACTags(84)%Label / 'instr_8_sn' / &
         ACTags(85)%Label / 'instr_8_id' / &
         ACTags(86)%Label / 'instr_8_wformat' / &
         ACTags(87)%Label / 'instr_8_wref' / &
         ACTags(88)%Label / 'instr_8_head_corr' / &
         ACTags(89)%Label / 'data_label' / &
         ACTags(90)%Label / 'col_1_variable' / &
         ACTags(91)%Label / 'col_1_useit' / &
         ACTags(92)%Label / 'col_1_measure_type' / &
         ACTags(93)%Label / 'col_1_instrument' / &
         ACTags(94)%Label / 'col_1_unit_in' / &
         ACTags(95)%Label / 'col_1_conversion' / &
         ACTags(96)%Label / 'col_1_unit_out' / &
         ACTags(97)%Label / 'col_2_variable' / &
         ACTags(98)%Label / 'col_2_useit' / &
         ACTags(99)%Label / 'col_2_measure_type' / &
         ACTags(100)%Label / 'col_2_instrument' / &
         ACTags(101)%Label / 'col_2_unit_in' / &
         ACTags(102)%Label / 'col_2_conversion' / &
         ACTags(103)%Label / 'col_2_unit_out' / &
         ACTags(104)%Label / 'col_3_variable' / &
         ACTags(105)%Label / 'col_3_useit' / &
         ACTags(106)%Label / 'col_3_measure_type' / &
         ACTags(107)%Label / 'col_3_instrument' / &
         ACTags(108)%Label / 'col_3_unit_in' / &
         ACTags(109)%Label / 'col_3_conversion' / &
         ACTags(110)%Label / 'col_3_unit_out' / &
         ACTags(111)%Label / 'col_4_variable' / &
         ACTags(112)%Label / 'col_4_useit' / &
         ACTags(113)%Label / 'col_4_measure_type' / &
         ACTags(114)%Label / 'col_4_instrument' / &
         ACTags(115)%Label / 'col_4_unit_in' / &
         ACTags(116)%Label / 'col_4_conversion' / &
         ACTags(117)%Label / 'col_4_unit_out' / &
         ACTags(118)%Label / 'col_5_variable' / &
         ACTags(119)%Label / 'col_5_useit' / &
         ACTags(120)%Label / 'col_5_measure_type' / &
         ACTags(121)%Label / 'col_5_instrument' / &
         ACTags(122)%Label / 'col_5_unit_in' / &
         ACTags(123)%Label / 'col_5_conversion' / &
         ACTags(124)%Label / 'col_5_unit_out' / &
         ACTags(125)%Label / 'col_6_variable' / &
         ACTags(126)%Label / 'col_6_useit' / &
         ACTags(127)%Label / 'col_6_measure_type' / &
         ACTags(128)%Label / 'col_6_instrument' / &
         ACTags(129)%Label / 'col_6_unit_in' / &
         ACTags(130)%Label / 'col_6_conversion' / &
         ACTags(131)%Label / 'col_6_unit_out' / &
         ACTags(132)%Label / 'col_7_variable' / &
         ACTags(133)%Label / 'col_7_useit' / &
         ACTags(134)%Label / 'col_7_measure_type' / &
         ACTags(135)%Label / 'col_7_instrument' / &
         ACTags(136)%Label / 'col_7_unit_in' / &
         ACTags(137)%Label / 'col_7_conversion' / &
         ACTags(138)%Label / 'col_7_unit_out' / &
         ACTags(139)%Label / 'col_8_variable' / &
         ACTags(140)%Label / 'col_8_useit' / &
         ACTags(141)%Label / 'col_8_measure_type' / &
         ACTags(142)%Label / 'col_8_instrument' / &
         ACTags(143)%Label / 'col_8_unit_in' / &
         ACTags(144)%Label / 'col_8_conversion' / &
         ACTags(145)%Label / 'col_8_unit_out' / &
         ACTags(146)%Label / 'col_9_variable' / &
         ACTags(147)%Label / 'col_9_useit' / &
         ACTags(148)%Label / 'col_9_measure_type' / &
         ACTags(149)%Label / 'col_9_instrument' / &
         ACTags(150)%Label / 'col_9_unit_in' / &
         ACTags(151)%Label / 'col_9_conversion' / &
         ACTags(152)%Label / 'col_9_unit_out' / &
         ACTags(153)%Label / 'col_10_variable' / &
         ACTags(154)%Label / 'col_10_useit' / &
         ACTags(155)%Label / 'col_10_measure_type' / &
         ACTags(156)%Label / 'col_10_instrument' / &
         ACTags(157)%Label / 'col_10_unit_in' / &
         ACTags(158)%Label / 'col_10_conversion' / &
         ACTags(159)%Label / 'col_10_unit_out' / &
         ACTags(160)%Label / 'col_11_variable' / &
         ACTags(161)%Label / 'col_11_useit' / &
         ACTags(162)%Label / 'col_11_measure_type' / &
         ACTags(163)%Label / 'col_11_instrument' / &
         ACTags(164)%Label / 'col_11_unit_in' / &
         ACTags(165)%Label / 'col_11_conversion' / &
         ACTags(166)%Label / 'col_11_unit_out' / &
         ACTags(167)%Label / 'col_12_variable' / &
         ACTags(168)%Label / 'col_12_useit' / &
         ACTags(169)%Label / 'col_12_measure_type' / &
         ACTags(170)%Label / 'col_12_instrument' / &
         ACTags(171)%Label / 'col_12_unit_in' / &
         ACTags(172)%Label / 'col_12_conversion' / &
         ACTags(173)%Label / 'col_12_unit_out' / &
         ACTags(174)%Label / 'col_13_variable' / &
         ACTags(175)%Label / 'col_13_useit' / &
         ACTags(176)%Label / 'col_13_measure_type' / &
         ACTags(177)%Label / 'col_13_instrument' / &
         ACTags(178)%Label / 'col_13_unit_in' / &
         ACTags(179)%Label / 'col_13_conversion' / &
         ACTags(180)%Label / 'col_13_unit_out' / &
         ACTags(181)%Label / 'col_14_variable' / &
         ACTags(182)%Label / 'col_14_useit' / &
         ACTags(183)%Label / 'col_14_measure_type' / &
         ACTags(184)%Label / 'col_14_instrument' / &
         ACTags(185)%Label / 'col_14_unit_in' / &
         ACTags(186)%Label / 'col_14_conversion' / &
         ACTags(187)%Label / 'col_14_unit_out' / &
         ACTags(188)%Label / 'col_15_variable' / &
         ACTags(189)%Label / 'col_15_useit' / &
         ACTags(190)%Label / 'col_15_measure_type' / &
         ACTags(191)%Label / 'col_15_instrument' / &
         ACTags(192)%Label / 'col_15_unit_in' / &
         ACTags(193)%Label / 'col_15_conversion' / &
         ACTags(194)%Label / 'col_15_unit_out' / &
         ACTags(195)%Label / 'col_16_variable' / &
         ACTags(196)%Label / 'col_16_useit' / &
         ACTags(197)%Label / 'col_16_measure_type' / &
         ACTags(198)%Label / 'col_16_instrument' / &
         ACTags(199)%Label / 'col_16_unit_in' / &
         ACTags(200)%Label / 'col_16_conversion' /
    data ACTags(201)%Label / 'col_16_unit_out' / &
         ACTags(202)%Label / 'col_17_variable' / &
         ACTags(203)%Label / 'col_17_useit' / &
         ACTags(204)%Label / 'col_17_measure_type' / &
         ACTags(205)%Label / 'col_17_instrument' / &
         ACTags(206)%Label / 'col_17_unit_in' / &
         ACTags(207)%Label / 'col_17_conversion' / &
         ACTags(208)%Label / 'col_17_unit_out' / &
         ACTags(209)%Label / 'col_18_variable' / &
         ACTags(210)%Label / 'col_18_useit' / &
         ACTags(211)%Label / 'col_18_measure_type' / &
         ACTags(212)%Label / 'col_18_instrument' / &
         ACTags(213)%Label / 'col_18_unit_in' / &
         ACTags(214)%Label / 'col_18_conversion' / &
         ACTags(215)%Label / 'col_18_unit_out' / &
         ACTags(216)%Label / 'col_19_variable' / &
         ACTags(217)%Label / 'col_19_useit' / &
         ACTags(218)%Label / 'col_19_measure_type' / &
         ACTags(219)%Label / 'col_19_instrument' / &
         ACTags(220)%Label / 'col_19_unit_in' / &
         ACTags(221)%Label / 'col_19_conversion' / &
         ACTags(222)%Label / 'col_19_unit_out' / &
         ACTags(223)%Label / 'col_20_variable' / &
         ACTags(224)%Label / 'col_20_useit' / &
         ACTags(225)%Label / 'col_20_measure_type' / &
         ACTags(226)%Label / 'col_20_instrument' / &
         ACTags(227)%Label / 'col_20_unit_in' / &
         ACTags(228)%Label / 'col_20_conversion' / &
         ACTags(229)%Label / 'col_20_unit_out' / &
         ACTags(230)%Label / 'col_21_variable' / &
         ACTags(231)%Label / 'col_21_useit' / &
         ACTags(232)%Label / 'col_21_measure_type' / &
         ACTags(233)%Label / 'col_21_instrument' / &
         ACTags(234)%Label / 'col_21_unit_in' / &
         ACTags(235)%Label / 'col_21_conversion' / &
         ACTags(236)%Label / 'col_21_unit_out' / &
         ACTags(237)%Label / 'col_22_variable' / &
         ACTags(238)%Label / 'col_22_useit' / &
         ACTags(239)%Label / 'col_22_measure_type' / &
         ACTags(240)%Label / 'col_22_instrument' / &
         ACTags(241)%Label / 'col_22_unit_in' / &
         ACTags(242)%Label / 'col_22_conversion' / &
         ACTags(243)%Label / 'col_22_unit_out' / &
         ACTags(244)%Label / 'col_23_variable' / &
         ACTags(245)%Label / 'col_23_useit' / &
         ACTags(246)%Label / 'col_23_measure_type' / &
         ACTags(247)%Label / 'col_23_instrument' / &
         ACTags(248)%Label / 'col_23_unit_in' / &
         ACTags(249)%Label / 'col_23_conversion' / &
         ACTags(250)%Label / 'col_23_unit_out' / &
         ACTags(251)%Label / 'col_24_variable' / &
         ACTags(252)%Label / 'col_24_useit' / &
         ACTags(253)%Label / 'col_24_measure_type' / &
         ACTags(254)%Label / 'col_24_instrument' / &
         ACTags(255)%Label / 'col_24_unit_in' / &
         ACTags(256)%Label / 'col_24_conversion' / &
         ACTags(257)%Label / 'col_24_unit_out' / &
         ACTags(258)%Label / 'col_25_variable' / &
         ACTags(259)%Label / 'col_25_useit' / &
         ACTags(260)%Label / 'col_25_measure_type' / &
         ACTags(261)%Label / 'col_25_instrument' / &
         ACTags(262)%Label / 'col_25_unit_in' / &
         ACTags(263)%Label / 'col_25_conversion' / &
         ACTags(264)%Label / 'col_25_unit_out' / &
         ACTags(265)%Label / 'col_26_variable' / &
         ACTags(266)%Label / 'col_26_useit' / &
         ACTags(267)%Label / 'col_26_measure_type' / &
         ACTags(268)%Label / 'col_26_instrument' / &
         ACTags(269)%Label / 'col_26_unit_in' / &
         ACTags(270)%Label / 'col_26_conversion' / &
         ACTags(271)%Label / 'col_26_unit_out' / &
         ACTags(272)%Label / 'col_27_variable' / &
         ACTags(273)%Label / 'col_27_useit' / &
         ACTags(274)%Label / 'col_27_measure_type' / &
         ACTags(275)%Label / 'col_27_instrument' / &
         ACTags(276)%Label / 'col_27_unit_in' / &
         ACTags(277)%Label / 'col_27_conversion' / &
         ACTags(278)%Label / 'col_27_unit_out' / &
         ACTags(279)%Label / 'col_28_variable' / &
         ACTags(280)%Label / 'col_28_useit' / &
         ACTags(281)%Label / 'col_28_measure_type' / &
         ACTags(282)%Label / 'col_28_instrument' / &
         ACTags(283)%Label / 'col_28_unit_in' / &
         ACTags(284)%Label / 'col_28_conversion' / &
         ACTags(285)%Label / 'col_28_unit_out' / &
         ACTags(286)%Label / 'col_29_variable' / &
         ACTags(287)%Label / 'col_29_useit' / &
         ACTags(288)%Label / 'col_29_measure_type' / &
         ACTags(289)%Label / 'col_29_instrument' / &
         ACTags(290)%Label / 'col_29_unit_in' / &
         ACTags(291)%Label / 'col_29_conversion' / &
         ACTags(292)%Label / 'col_29_unit_out' / &
         ACTags(293)%Label / 'col_30_variable' / &
         ACTags(294)%Label / 'col_30_useit' / &
         ACTags(295)%Label / 'col_30_measure_type' / &
         ACTags(296)%Label / 'col_30_instrument' / &
         ACTags(297)%Label / 'col_30_unit_in' / &
         ACTags(298)%Label / 'col_30_conversion' / &
         ACTags(299)%Label / 'col_30_unit_out' / &
         ACTags(300)%Label / 'col_31_variable' / &
         ACTags(301)%Label / 'col_31_useit' / &
         ACTags(302)%Label / 'col_31_measure_type' / &
         ACTags(303)%Label / 'col_31_instrument' / &
         ACTags(304)%Label / 'col_31_unit_in' / &
         ACTags(305)%Label / 'col_31_conversion' / &
         ACTags(306)%Label / 'col_31_unit_out' / &
         ACTags(307)%Label / 'col_32_variable' / &
         ACTags(308)%Label / 'col_32_useit' / &
         ACTags(309)%Label / 'col_32_measure_type' / &
         ACTags(310)%Label / 'col_32_instrument' / &
         ACTags(311)%Label / 'col_32_unit_in' / &
         ACTags(312)%Label / 'col_32_conversion' / &
         ACTags(313)%Label / 'col_32_unit_out' / &
         ACTags(314)%Label / 'col_33_variable' / &
         ACTags(315)%Label / 'col_33_useit' / &
         ACTags(316)%Label / 'col_33_measure_type' / &
         ACTags(317)%Label / 'col_33_instrument' / &
         ACTags(318)%Label / 'col_33_unit_in' / &
         ACTags(319)%Label / 'col_33_conversion' / &
         ACTags(320)%Label / 'col_33_unit_out' / &
         ACTags(321)%Label / 'col_34_variable' / &
         ACTags(322)%Label / 'col_34_useit' / &
         ACTags(323)%Label / 'col_34_measure_type' / &
         ACTags(324)%Label / 'col_34_instrument' / &
         ACTags(325)%Label / 'col_34_unit_in' / &
         ACTags(326)%Label / 'col_34_conversion' / &
         ACTags(327)%Label / 'col_34_unit_out' / &
         ACTags(328)%Label / 'col_35_variable' / &
         ACTags(329)%Label / 'col_35_useit' / &
         ACTags(330)%Label / 'col_35_measure_type' / &
         ACTags(331)%Label / 'col_35_instrument' / &
         ACTags(332)%Label / 'col_35_unit_in' / &
         ACTags(333)%Label / 'col_35_conversion' / &
         ACTags(334)%Label / 'col_35_unit_out' / &
         ACTags(335)%Label / 'col_36_variable' / &
         ACTags(336)%Label / 'col_36_useit' / &
         ACTags(337)%Label / 'col_36_measure_type' / &
         ACTags(338)%Label / 'col_36_instrument' / &
         ACTags(339)%Label / 'col_36_unit_in' / &
         ACTags(340)%Label / 'col_36_conversion' / &
         ACTags(341)%Label / 'col_36_unit_out' / &
         ACTags(342)%Label / 'col_37_variable' / &
         ACTags(343)%Label / 'col_37_useit' / &
         ACTags(344)%Label / 'col_37_measure_type' / &
         ACTags(345)%Label / 'col_37_instrument' / &
         ACTags(346)%Label / 'col_37_unit_in' / &
         ACTags(347)%Label / 'col_37_conversion' / &
         ACTags(348)%Label / 'col_37_unit_out' / &
         ACTags(349)%Label / 'col_38_variable' / &
         ACTags(350)%Label / 'col_38_useit' / &
         ACTags(351)%Label / 'col_38_measure_type' / &
         ACTags(352)%Label / 'col_38_instrument' / &
         ACTags(353)%Label / 'col_38_unit_in' / &
         ACTags(354)%Label / 'col_38_conversion' / &
         ACTags(355)%Label / 'col_38_unit_out' / &
         ACTags(356)%Label / 'col_39_variable' / &
         ACTags(357)%Label / 'col_39_useit' / &
         ACTags(358)%Label / 'col_39_measure_type' / &
         ACTags(359)%Label / 'col_39_instrument' / &
         ACTags(360)%Label / 'col_39_unit_in' / &
         ACTags(361)%Label / 'col_39_conversion' / &
         ACTags(362)%Label / 'col_39_unit_out' / &
         ACTags(363)%Label / 'col_40_variable' / &
         ACTags(364)%Label / 'col_40_useit' / &
         ACTags(365)%Label / 'col_40_measure_type' / &
         ACTags(366)%Label / 'col_40_instrument' / &
         ACTags(367)%Label / 'col_40_unit_in' / &
         ACTags(368)%Label / 'col_40_conversion' / &
         ACTags(369)%Label / 'col_40_unit_out' / &
         ACTags(370)%Label / 'col_41_variable' / &
         ACTags(371)%Label / 'col_41_useit' / &
         ACTags(372)%Label / 'col_41_measure_type' / &
         ACTags(373)%Label / 'col_41_instrument' / &
         ACTags(374)%Label / 'col_41_unit_in' / &
         ACTags(375)%Label / 'col_41_conversion' / &
         ACTags(376)%Label / 'col_41_unit_out' / &
         ACTags(377)%Label / 'col_42_variable' / &
         ACTags(378)%Label / 'col_42_useit' / &
         ACTags(379)%Label / 'col_42_measure_type' / &
         ACTags(380)%Label / 'col_42_instrument' / &
         ACTags(381)%Label / 'col_42_unit_in' / &
         ACTags(382)%Label / 'col_42_conversion' / &
         ACTags(383)%Label / 'col_42_unit_out' / &
         ACTags(384)%Label / 'col_43_variable' / &
         ACTags(385)%Label / 'col_43_useit' / &
         ACTags(386)%Label / 'col_43_measure_type' / &
         ACTags(387)%Label / 'col_43_instrument' / &
         ACTags(388)%Label / 'col_43_unit_in' / &
         ACTags(389)%Label / 'col_43_conversion' / &
         ACTags(390)%Label / 'col_43_unit_out' / &
         ACTags(391)%Label / 'col_44_variable' / &
         ACTags(392)%Label / 'col_44_useit' / &
         ACTags(393)%Label / 'col_44_measure_type' / &
         ACTags(394)%Label / 'col_44_instrument' / &
         ACTags(395)%Label / 'col_44_unit_in' / &
         ACTags(396)%Label / 'col_44_conversion' / &
         ACTags(397)%Label / 'col_44_unit_out' / &
         ACTags(398)%Label / 'col_45_variable' / &
         ACTags(399)%Label / 'col_45_useit' / &
         ACTags(400)%Label / 'col_45_measure_type' /
    data ACTags(401)%Label / 'col_45_instrument' / &
         ACTags(402)%Label / 'col_45_unit_in' / &
         ACTags(403)%Label / 'col_45_conversion' / &
         ACTags(404)%Label / 'col_45_unit_out' / &
         ACTags(405)%Label / 'col_46_variable' / &
         ACTags(406)%Label / 'col_46_useit' / &
         ACTags(407)%Label / 'col_46_measure_type' / &
         ACTags(408)%Label / 'col_46_instrument' / &
         ACTags(409)%Label / 'col_46_unit_in' / &
         ACTags(410)%Label / 'col_46_conversion' / &
         ACTags(411)%Label / 'col_46_unit_out' / &
         ACTags(412)%Label / 'col_47_variable' / &
         ACTags(413)%Label / 'col_47_useit' / &
         ACTags(414)%Label / 'col_47_measure_type' / &
         ACTags(415)%Label / 'col_47_instrument' / &
         ACTags(416)%Label / 'col_47_unit_in' / &
         ACTags(417)%Label / 'col_47_conversion' / &
         ACTags(418)%Label / 'col_47_unit_out' / &
         ACTags(419)%Label / 'col_48_variable' / &
         ACTags(420)%Label / 'col_48_useit' / &
         ACTags(421)%Label / 'col_48_measure_type' / &
         ACTags(422)%Label / 'col_48_instrument' / &
         ACTags(423)%Label / 'col_48_unit_in' / &
         ACTags(424)%Label / 'col_48_conversion' / &
         ACTags(425)%Label / 'col_48_unit_out' / &
         ACTags(426)%Label / 'col_49_variable' / &
         ACTags(427)%Label / 'col_49_useit' / &
         ACTags(428)%Label / 'col_49_measure_type' / &
         ACTags(429)%Label / 'col_49_instrument' / &
         ACTags(430)%Label / 'col_49_unit_in' / &
         ACTags(431)%Label / 'col_49_conversion' / &
         ACTags(432)%Label / 'col_49_unit_out' / &
         ACTags(433)%Label / 'col_50_variable' / &
         ACTags(434)%Label / 'col_50_useit' / &
         ACTags(435)%Label / 'col_50_measure_type' / &
         ACTags(436)%Label / 'col_50_instrument' / &
         ACTags(437)%Label / 'col_50_unit_in' / &
         ACTags(438)%Label / 'col_50_conversion' / &
         ACTags(439)%Label / 'col_50_unit_out' / &
         ACTags(440)%Label / 'col_51_variable' / &
         ACTags(441)%Label / 'col_51_useit' / &
         ACTags(442)%Label / 'col_51_measure_type' / &
         ACTags(443)%Label / 'col_51_instrument' / &
         ACTags(444)%Label / 'col_51_unit_in' / &
         ACTags(445)%Label / 'col_51_conversion' / &
         ACTags(446)%Label / 'col_51_unit_out' / &
         ACTags(447)%Label / 'col_52_variable' / &
         ACTags(448)%Label / 'col_52_useit' / &
         ACTags(449)%Label / 'col_52_measure_type' / &
         ACTags(450)%Label / 'col_52_instrument' / &
         ACTags(451)%Label / 'col_52_unit_in' / &
         ACTags(452)%Label / 'col_52_conversion' / &
         ACTags(453)%Label / 'col_52_unit_out' / &
         ACTags(454)%Label / 'col_53_variable' / &
         ACTags(455)%Label / 'col_53_useit' / &
         ACTags(456)%Label / 'col_53_measure_type' / &
         ACTags(457)%Label / 'col_53_instrument' / &
         ACTags(458)%Label / 'col_53_unit_in' / &
         ACTags(459)%Label / 'col_53_conversion' / &
         ACTags(460)%Label / 'col_53_unit_out' / &
         ACTags(461)%Label / 'col_54_variable' / &
         ACTags(462)%Label / 'col_54_useit' / &
         ACTags(463)%Label / 'col_54_measure_type' / &
         ACTags(464)%Label / 'col_54_instrument' / &
         ACTags(465)%Label / 'col_54_unit_in' / &
         ACTags(466)%Label / 'col_54_conversion' / &
         ACTags(467)%Label / 'col_54_unit_out' / &
         ACTags(468)%Label / 'col_55_variable' / &
         ACTags(469)%Label / 'col_55_useit' / &
         ACTags(470)%Label / 'col_55_measure_type' / &
         ACTags(471)%Label / 'col_55_instrument' / &
         ACTags(472)%Label / 'col_55_unit_in' / &
         ACTags(473)%Label / 'col_55_conversion' / &
         ACTags(474)%Label / 'col_55_unit_out' / &
         ACTags(475)%Label / 'col_56_variable' / &
         ACTags(476)%Label / 'col_56_useit' / &
         ACTags(477)%Label / 'col_56_measure_type' / &
         ACTags(478)%Label / 'col_56_instrument' / &
         ACTags(479)%Label / 'col_56_unit_in' / &
         ACTags(480)%Label / 'col_56_conversion' / &
         ACTags(481)%Label / 'col_56_unit_out' / &
         ACTags(482)%Label / 'col_57_variable' / &
         ACTags(483)%Label / 'col_57_useit' / &
         ACTags(484)%Label / 'col_57_measure_type' / &
         ACTags(485)%Label / 'col_57_instrument' / &
         ACTags(486)%Label / 'col_57_unit_in' / &
         ACTags(487)%Label / 'col_57_conversion' / &
         ACTags(488)%Label / 'col_57_unit_out' / &
         ACTags(489)%Label / 'col_58_variable' / &
         ACTags(490)%Label / 'col_58_useit' / &
         ACTags(491)%Label / 'col_58_measure_type' / &
         ACTags(492)%Label / 'col_58_instrument' / &
         ACTags(493)%Label / 'col_58_unit_in' / &
         ACTags(494)%Label / 'col_58_conversion' / &
         ACTags(495)%Label / 'col_58_unit_out' / &
         ACTags(496)%Label / 'col_59_variable' / &
         ACTags(497)%Label / 'col_59_useit' / &
         ACTags(498)%Label / 'col_59_measure_type' / &
         ACTags(499)%Label / 'col_59_instrument' / &
         ACTags(500)%Label / 'col_59_unit_in' / &
         ACTags(501)%Label / 'col_59_conversion' / &
         ACTags(502)%Label / 'col_59_unit_out' / &
         ACTags(503)%Label / 'col_60_variable' / &
         ACTags(504)%Label / 'col_60_useit' / &
         ACTags(505)%Label / 'col_60_measure_type' / &
         ACTags(506)%Label / 'col_60_instrument' / &
         ACTags(507)%Label / 'col_60_unit_in' / &
         ACTags(508)%Label / 'col_60_conversion' / &
         ACTags(509)%Label / 'col_60_unit_out' / &
         ACTags(510)%Label / 'col_61_variable' / &
         ACTags(511)%Label / 'col_61_useit' / &
         ACTags(512)%Label / 'col_61_measure_type' / &
         ACTags(513)%Label / 'col_61_instrument' / &
         ACTags(514)%Label / 'col_61_unit_in' / &
         ACTags(515)%Label / 'col_61_conversion' / &
         ACTags(516)%Label / 'col_61_unit_out' / &
         ACTags(517)%Label / 'col_62_variable' / &
         ACTags(518)%Label / 'col_62_useit' / &
         ACTags(519)%Label / 'col_62_measure_type' / &
         ACTags(520)%Label / 'col_62_instrument' / &
         ACTags(521)%Label / 'col_62_unit_in' / &
         ACTags(522)%Label / 'col_62_conversion' / &
         ACTags(523)%Label / 'col_62_unit_out' / &
         ACTags(524)%Label / 'col_63_variable' / &
         ACTags(525)%Label / 'col_63_useit' / &
         ACTags(526)%Label / 'col_63_measure_type' / &
         ACTags(527)%Label / 'col_63_instrument' / &
         ACTags(528)%Label / 'col_63_unit_in' / &
         ACTags(529)%Label / 'col_63_conversion' / &
         ACTags(530)%Label / 'col_63_unit_out' / &
         ACTags(531)%Label / 'col_64_variable' / &
         ACTags(532)%Label / 'col_64_useit' / &
         ACTags(533)%Label / 'col_64_measure_type' / &
         ACTags(534)%Label / 'col_64_instrument' / &
         ACTags(535)%Label / 'col_64_unit_in' / &
         ACTags(536)%Label / 'col_64_conversion' / &
         ACTags(537)%Label / 'col_64_unit_out' / &
         ACTags(538)%Label / 'col_65_variable' / &
         ACTags(539)%Label / 'col_65_useit' / &
         ACTags(540)%Label / 'col_65_measure_type' / &
         ACTags(541)%Label / 'col_65_instrument' / &
         ACTags(542)%Label / 'col_65_unit_in' / &
         ACTags(543)%Label / 'col_65_conversion' / &
         ACTags(544)%Label / 'col_65_unit_out' / &
         ACTags(545)%Label / 'col_66_variable' / &
         ACTags(546)%Label / 'col_66_useit' / &
         ACTags(547)%Label / 'col_66_measure_type' / &
         ACTags(548)%Label / 'col_66_instrument' / &
         ACTags(549)%Label / 'col_66_unit_in' / &
         ACTags(550)%Label / 'col_66_conversion' / &
         ACTags(551)%Label / 'col_66_unit_out' / &
         ACTags(552)%Label / 'col_67_variable' / &
         ACTags(553)%Label / 'col_67_useit' / &
         ACTags(554)%Label / 'col_67_measure_type' / &
         ACTags(555)%Label / 'col_67_instrument' / &
         ACTags(556)%Label / 'col_67_unit_in' / &
         ACTags(557)%Label / 'col_67_conversion' / &
         ACTags(558)%Label / 'col_67_unit_out' / &
         ACTags(559)%Label / 'col_68_variable' / &
         ACTags(560)%Label / 'col_68_useit' / &
         ACTags(561)%Label / 'col_68_measure_type' / &
         ACTags(562)%Label / 'col_68_instrument' / &
         ACTags(563)%Label / 'col_68_unit_in' / &
         ACTags(564)%Label / 'col_68_conversion' / &
         ACTags(565)%Label / 'col_68_unit_out' / &
         ACTags(566)%Label / 'col_69_variable' / &
         ACTags(567)%Label / 'col_69_useit' / &
         ACTags(568)%Label / 'col_69_measure_type' / &
         ACTags(569)%Label / 'col_69_instrument' / &
         ACTags(570)%Label / 'col_69_unit_in' / &
         ACTags(571)%Label / 'col_69_conversion' / &
         ACTags(572)%Label / 'col_69_unit_out' / &
         ACTags(573)%Label / 'col_70_variable' / &
         ACTags(574)%Label / 'col_70_useit' / &
         ACTags(575)%Label / 'col_70_measure_type' / &
         ACTags(576)%Label / 'col_70_instrument' / &
         ACTags(577)%Label / 'col_70_unit_in' / &
         ACTags(578)%Label / 'col_70_conversion' / &
         ACTags(579)%Label / 'col_70_unit_out' / &
         ACTags(580)%Label / 'col_71_variable' / &
         ACTags(581)%Label / 'col_71_useit' / &
         ACTags(582)%Label / 'col_71_measure_type' / &
         ACTags(583)%Label / 'col_71_instrument' / &
         ACTags(584)%Label / 'col_71_unit_in' / &
         ACTags(585)%Label / 'col_71_conversion' / &
         ACTags(586)%Label / 'col_71_unit_out' / &
         ACTags(587)%Label / 'col_72_variable' / &
         ACTags(588)%Label / 'col_72_useit' / &
         ACTags(589)%Label / 'col_72_measure_type' / &
         ACTags(590)%Label / 'col_72_instrument' / &
         ACTags(591)%Label / 'col_72_unit_in' / &
         ACTags(592)%Label / 'col_72_conversion' / &
         ACTags(593)%Label / 'col_72_unit_out' / &
         ACTags(594)%Label / 'col_73_variable' / &
         ACTags(595)%Label / 'col_73_useit' / &
         ACTags(596)%Label / 'col_73_measure_type' / &
         ACTags(597)%Label / 'col_73_instrument' / &
         ACTags(598)%Label / 'col_73_unit_in' / &
         ACTags(599)%Label / 'col_73_conversion' / &
         ACTags(600)%Label / 'col_73_unit_out' /
    data ACTags(601)%Label / 'col_74_variable' / &
         ACTags(602)%Label / 'col_74_useit' / &
         ACTags(603)%Label / 'col_74_measure_type' / &
         ACTags(604)%Label / 'col_74_instrument' / &
         ACTags(605)%Label / 'col_74_unit_in' / &
         ACTags(606)%Label / 'col_74_conversion' / &
         ACTags(607)%Label / 'col_74_unit_out' / &
         ACTags(608)%Label / 'col_75_variable' / &
         ACTags(609)%Label / 'col_75_useit' / &
         ACTags(610)%Label / 'col_75_measure_type' / &
         ACTags(611)%Label / 'col_75_instrument' / &
         ACTags(612)%Label / 'col_75_unit_in' / &
         ACTags(613)%Label / 'col_75_conversion' / &
         ACTags(614)%Label / 'col_75_unit_out' / &
         ACTags(615)%Label / 'col_76_variable' / &
         ACTags(616)%Label / 'col_76_useit' / &
         ACTags(617)%Label / 'col_76_measure_type' / &
         ACTags(618)%Label / 'col_76_instrument' / &
         ACTags(619)%Label / 'col_76_unit_in' / &
         ACTags(620)%Label / 'col_76_conversion' / &
         ACTags(621)%Label / 'col_76_unit_out' / &
         ACTags(622)%Label / 'col_77_variable' / &
         ACTags(623)%Label / 'col_77_useit' / &
         ACTags(624)%Label / 'col_77_measure_type' / &
         ACTags(625)%Label / 'col_77_instrument' / &
         ACTags(626)%Label / 'col_77_unit_in' / &
         ACTags(627)%Label / 'col_77_conversion' / &
         ACTags(628)%Label / 'col_77_unit_out' / &
         ACTags(629)%Label / 'col_78_variable' / &
         ACTags(630)%Label / 'col_78_useit' / &
         ACTags(631)%Label / 'col_78_measure_type' / &
         ACTags(632)%Label / 'col_78_instrument' / &
         ACTags(633)%Label / 'col_78_unit_in' / &
         ACTags(634)%Label / 'col_78_conversion' / &
         ACTags(635)%Label / 'col_78_unit_out' / &
         ACTags(636)%Label / 'col_79_variable' / &
         ACTags(637)%Label / 'col_79_useit' / &
         ACTags(638)%Label / 'col_79_measure_type' / &
         ACTags(639)%Label / 'col_79_instrument' / &
         ACTags(640)%Label / 'col_79_unit_in' / &
         ACTags(641)%Label / 'col_79_conversion' / &
         ACTags(642)%Label / 'col_79_unit_out' / &
         ACTags(643)%Label / 'col_80_variable' / &
         ACTags(644)%Label / 'col_80_useit' / &
         ACTags(645)%Label / 'col_80_measure_type' / &
         ACTags(646)%Label / 'col_80_instrument' / &
         ACTags(647)%Label / 'col_80_unit_in' / &
         ACTags(648)%Label / 'col_80_conversion' / &
         ACTags(649)%Label / 'col_80_unit_out' / &
         ACTags(650)%Label / 'col_81_variable' / &
         ACTags(651)%Label / 'col_81_useit' / &
         ACTags(652)%Label / 'col_81_measure_type' / &
         ACTags(653)%Label / 'col_81_instrument' / &
         ACTags(654)%Label / 'col_81_unit_in' / &
         ACTags(655)%Label / 'col_81_conversion' / &
         ACTags(656)%Label / 'col_81_unit_out' / &
         ACTags(657)%Label / 'col_82_variable' / &
         ACTags(658)%Label / 'col_82_useit' / &
         ACTags(659)%Label / 'col_82_measure_type' / &
         ACTags(660)%Label / 'col_82_instrument' / &
         ACTags(661)%Label / 'col_82_unit_in' / &
         ACTags(662)%Label / 'col_82_conversion' / &
         ACTags(663)%Label / 'col_82_unit_out' / &
         ACTags(664)%Label / 'col_83_variable' / &
         ACTags(665)%Label / 'col_83_useit' / &
         ACTags(666)%Label / 'col_83_measure_type' / &
         ACTags(667)%Label / 'col_83_instrument' / &
         ACTags(668)%Label / 'col_83_unit_in' / &
         ACTags(669)%Label / 'col_83_conversion' / &
         ACTags(670)%Label / 'col_83_unit_out' / &
         ACTags(671)%Label / 'col_84_variable' / &
         ACTags(672)%Label / 'col_84_useit' / &
         ACTags(673)%Label / 'col_84_measure_type' / &
         ACTags(674)%Label / 'col_84_instrument' / &
         ACTags(675)%Label / 'col_84_unit_in' / &
         ACTags(676)%Label / 'col_84_conversion' / &
         ACTags(677)%Label / 'col_84_unit_out' / &
         ACTags(678)%Label / 'col_85_variable' / &
         ACTags(679)%Label / 'col_85_useit' / &
         ACTags(680)%Label / 'col_85_measure_type' / &
         ACTags(681)%Label / 'col_85_instrument' / &
         ACTags(682)%Label / 'col_85_unit_in' / &
         ACTags(683)%Label / 'col_85_conversion' / &
         ACTags(684)%Label / 'col_85_unit_out' / &
         ACTags(685)%Label / 'col_86_variable' / &
         ACTags(686)%Label / 'col_86_useit' / &
         ACTags(687)%Label / 'col_86_measure_type' / &
         ACTags(688)%Label / 'col_86_instrument' / &
         ACTags(689)%Label / 'col_86_unit_in' / &
         ACTags(690)%Label / 'col_86_conversion' / &
         ACTags(691)%Label / 'col_86_unit_out' / &
         ACTags(692)%Label / 'col_87_variable' / &
         ACTags(693)%Label / 'col_87_useit' / &
         ACTags(694)%Label / 'col_87_measure_type' / &
         ACTags(695)%Label / 'col_87_instrument' / &
         ACTags(696)%Label / 'col_87_unit_in' / &
         ACTags(697)%Label / 'col_87_conversion' / &
         ACTags(698)%Label / 'col_87_unit_out' / &
         ACTags(699)%Label / 'col_88_variable' / &
         ACTags(700)%Label / 'col_88_useit' / &
         ACTags(701)%Label / 'col_88_measure_type' / &
         ACTags(702)%Label / 'col_88_instrument' / &
         ACTags(703)%Label / 'col_88_unit_in' / &
         ACTags(704)%Label / 'col_88_conversion' / &
         ACTags(705)%Label / 'col_88_unit_out' / &
         ACTags(706)%Label / 'col_89_variable' / &
         ACTags(707)%Label / 'col_89_useit' / &
         ACTags(708)%Label / 'col_89_measure_type' / &
         ACTags(709)%Label / 'col_89_instrument' / &
         ACTags(710)%Label / 'col_89_unit_in' / &
         ACTags(711)%Label / 'col_89_conversion' / &
         ACTags(712)%Label / 'col_89_unit_out' / &
         ACTags(713)%Label / 'col_90_variable' / &
         ACTags(714)%Label / 'col_90_useit' / &
         ACTags(715)%Label / 'col_90_measure_type' / &
         ACTags(716)%Label / 'col_90_instrument' / &
         ACTags(717)%Label / 'col_90_unit_in' / &
         ACTags(718)%Label / 'col_90_conversion' / &
         ACTags(719)%Label / 'col_90_unit_out' / &
         ACTags(720)%Label / 'col_91_variable' / &
         ACTags(721)%Label / 'col_91_useit' / &
         ACTags(722)%Label / 'col_91_measure_type' / &
         ACTags(723)%Label / 'col_91_instrument' / &
         ACTags(724)%Label / 'col_91_unit_in' / &
         ACTags(725)%Label / 'col_91_conversion' / &
         ACTags(726)%Label / 'col_91_unit_out' / &
         ACTags(727)%Label / 'col_92_variable' / &
         ACTags(728)%Label / 'col_92_useit' / &
         ACTags(729)%Label / 'col_92_measure_type' / &
         ACTags(730)%Label / 'col_92_instrument' / &
         ACTags(731)%Label / 'col_92_unit_in' / &
         ACTags(732)%Label / 'col_92_conversion' / &
         ACTags(733)%Label / 'col_92_unit_out' / &
         ACTags(734)%Label / 'col_93_variable' / &
         ACTags(735)%Label / 'col_93_useit' / &
         ACTags(736)%Label / 'col_93_measure_type' / &
         ACTags(737)%Label / 'col_93_instrument' / &
         ACTags(738)%Label / 'col_93_unit_in' / &
         ACTags(739)%Label / 'col_93_conversion' / &
         ACTags(740)%Label / 'col_93_unit_out' / &
         ACTags(741)%Label / 'col_94_variable' / &
         ACTags(742)%Label / 'col_94_useit' / &
         ACTags(743)%Label / 'col_94_measure_type' / &
         ACTags(744)%Label / 'col_94_instrument' / &
         ACTags(745)%Label / 'col_94_unit_in' / &
         ACTags(746)%Label / 'col_94_conversion' / &
         ACTags(747)%Label / 'col_94_unit_out' / &
         ACTags(748)%Label / 'col_95_variable' / &
         ACTags(749)%Label / 'col_95_useit' / &
         ACTags(750)%Label / 'col_95_measure_type' / &
         ACTags(751)%Label / 'col_95_instrument' / &
         ACTags(752)%Label / 'col_95_unit_in' / &
         ACTags(753)%Label / 'col_95_conversion' / &
         ACTags(754)%Label / 'col_95_unit_out' / &
         ACTags(755)%Label / 'col_96_variable' / &
         ACTags(756)%Label / 'col_96_useit' / &
         ACTags(757)%Label / 'col_96_measure_type' / &
         ACTags(758)%Label / 'col_96_instrument' / &
         ACTags(759)%Label / 'col_96_unit_in' / &
         ACTags(760)%Label / 'col_96_conversion' / &
         ACTags(761)%Label / 'col_96_unit_out' / &
         ACTags(762)%Label / 'col_97_variable' / &
         ACTags(763)%Label / 'col_97_useit' / &
         ACTags(764)%Label / 'col_97_measure_type' / &
         ACTags(765)%Label / 'col_97_instrument' / &
         ACTags(766)%Label / 'col_97_unit_in' / &
         ACTags(767)%Label / 'col_97_conversion' / &
         ACTags(768)%Label / 'col_97_unit_out' / &
         ACTags(769)%Label / 'col_98_variable' / &
         ACTags(770)%Label / 'col_98_useit' / &
         ACTags(771)%Label / 'col_98_measure_type' / &
         ACTags(772)%Label / 'col_98_instrument' / &
         ACTags(773)%Label / 'col_98_unit_in' / &
         ACTags(774)%Label / 'col_98_conversion' / &
         ACTags(775)%Label / 'col_98_unit_out' / &
         ACTags(776)%Label / 'col_99_variable' / &
         ACTags(777)%Label / 'col_99_useit' / &
         ACTags(778)%Label / 'col_99_measure_type' / &
         ACTags(779)%Label / 'col_99_instrument' / &
         ACTags(780)%Label / 'col_99_unit_in' / &
         ACTags(781)%Label / 'col_99_conversion' / &
         ACTags(782)%Label / 'col_99_unit_out' / &
         ACTags(783)%Label / 'col_100_variable' / &
         ACTags(784)%Label / 'col_100_useit' / &
         ACTags(785)%Label / 'col_100_measure_type' / &
         ACTags(786)%Label / 'col_100_instrument' / &
         ACTags(787)%Label / 'col_100_unit_in' / &
         ACTags(788)%Label / 'col_100_conversion' / &
         ACTags(789)%Label / 'col_100_unit_out' / &
         ACTags(790)%Label / 'instr_9_manufacturer' /
    !> END GENERATED ACTags
end module m_common_global_var
