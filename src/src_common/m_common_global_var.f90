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
    use m_log
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
    !> Raw column holding each gas's calibration reference, indexed by gas
    !> slot; 0 where the project supplies none. Was a single scalar that could
    !> only ever calibrate the fourth slot, so a site running a calibration
    !> reference against any other gas had it silently ignored.
    integer :: GasCalRefCol(GHGNumVar)
    !> Gas slot each user column calibrates, or 0. The 'cal-ref' marker on
    !> UserCol says only *that* a column is a calibration reference; which gas
    !> it belongs to used to be implicit, because there could be only one and
    !> it was always the fourth slot.
    integer :: UserCalRefSlot(MaxUserVar)

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
    character(4),  parameter :: Log_FilePadding         = '_log'
    character(11), parameter :: QCdetails_FilePadding  = '_qc_details'
    character(14), parameter :: Legend_FilePadding     = '_column_legend'
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
    !> The half-hourly PWB table. It is also the re-readable cache, so a run
    !> that wrote it can be replayed exactly from it.
    character(12), parameter :: PwbTimelag_FilePadding  = '_pwb_timelag'
    !> _pwb_diagnostics and _pwb_summary stood here. The first repeated the
    !> per-period cache one column apart and is folded into it; the second
    !> held per-gas tallies the run log already prints. PWB writes two files
    !> now: the half-hourly table (which is also the cache) and the aggregate
    !> _pwb_timelag_opt.
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
    !> Slot *defaults*, not species constants. Superseded per record by
    !> WriteProcessingProjectVariables, which writes each gas's own diffusivity
    !> from the project - including for slots one to four. Nothing may read
    !> Dc(co2) and mean "carbon dioxide": on a project that orders its records
    !> differently, slot five holds something else.
    data (Dc(mmm), mmm = histGas1, histGas4) / 0.00001381d0, 0.00002178d0, 0.00001952d0, 0.00001436d0/ !--> Massman (1998, Atm Env, Table 2)
    real(kind = sgl) :: MW(E2NumVar) !< Molecular weights
    !> Defaults keyed by legacy slot position, not by species. A record that
    !> carries its own mw supersedes the entry for its slot - including slots
    !> one to four - so nothing may read MW(co2) and mean "carbon dioxide" or
    !> MW(h2o) and mean "water". Use MW_H2O for the latter.
    data (MW(mmm), mmm = histGas1, histGas4) / 44.0095e-3, 18.0153e-3, 16.0425e-3, 44.0128e-3/
    !> The molecular weight of water as a physical constant, for the latent
    !> heat, evapotranspiration and vapour-density terms - all of which are
    !> about water itself and not about whatever gas occupies a slot. These
    !> read MW(h2o), which is the same number only for as long as record two
    !> happens to be water: a project declaring N2O there made every LE, E and
    !> ET wrong by a factor of 2.44.
    !>
    !> Deliberately `sgl` with the identical literal, so the promotion to
    !> double is bit-for-bit what the data statement above produced.
    real(kind = sgl), parameter :: MW_H2O = 18.0153e-3
    real(kind = dbl), parameter :: h2o_to_ET =  0.0648d0  !< To convert between H2O flux [mmol m-2 s-1] and ET flux (mm  hour-1)
    real(kind = dbl), parameter :: p = 3.14159265358979323846d0 !< Greek pi
    real(kind = dbl), parameter :: StdVair = 0.02245d0  !< gas molar volume at 25 �C and 101.325 kPa
    real(kind = dbl), parameter :: vk = 0.41d0 !< Von Karman constant
    !> Dry air and water, and the ratio the WPL correction runs on.
    !>
    !> Both were rounded to four figures, and their errors happened to cancel:
    !> 0.02897/0.01802 gave mu = 1.60766 against a true 1.60777, low by 0.007%.
    !> Refining water alone would have moved it to 1.60807 - high by 0.019%,
    !> worse than the pair it replaced. They only make sense together.
    real(kind = dbl), parameter :: Md = 0.0289647d0 !< molecular weight of dry air [kg_d/mol_d]
    real(kind = dbl), parameter :: mu = Md / 18.0153d-3
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

    !> "Derive this time-lag search window from the instrument geometry."
    !>
    !> AdjustTimelagOptSettings derives a window from tube transit time or from
    !> the sensor separations when a gas declares none, and takes the declared
    !> one otherwise. The two cases were told apart by a bare -1000 written out
    !> in one place and a project-file value of -1000.1 written in another, with
    !> nothing tying them together - so a gas left at the zero default read as
    !> "the user asked for [0, 0]" and searched nothing at all.
    !>
    !> The test is `< TlagDeriveThreshold`, so the default sits strictly below
    !> it. Real windows cannot reach here: the widest an open-path gas derives
    !> is twice its sensor separation, and a separation would have to exceed
    !> about five hundred metres.
    real(kind = dbl), parameter :: TlagDeriveThreshold = -1000d0
    real(kind = dbl), parameter :: TlagDeriveWindow = -1000.1d0
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
    character(4), parameter  :: LogExt                  = '.log'

    !> logging variables and parameters
    character(10) :: LogInteger
    logical :: LogAll = .false. !< working variable, for debug only
    logical :: co2_new_sw_ver = .false.

    integer, parameter :: ErrLab1 = 2

    !> gU..gPe were a second, seventeen-wide numbering of the same variables
    !> the E2Col slots already number, with four gases wired in at 8..11.
    !> Nothing read any of the seventeen. A duplicate numbering is how the
    !> instrument-role index came to cap analysers at four, so this one goes
    !> rather than waiting to be picked up.

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

    !> The site and sonic fields only. Fifty-six per-analyser entries used to
    !> follow - fourteen each for co2, h2o, ch4 and the fourth gas - and every
    !> one of the index constants naming them was unread: the gas fields are
    !> matched per record through DynMDGasFieldNames and GasSlotFromDynMDTag,
    !> which reach every analyser rather than the first four. A header column
    !> such as co2_irga_model still resolves; it just resolves on that pass.
    integer, parameter :: NumStdDynMDVars = 19
    character(64) :: StdDynMDVars(NumStdDynMDVars)
    data (StdDynMDVars(mmm), mmm = 1, NumStdDynMDVars) /'date', 'time', 'latitude', 'longitude', 'altitude',&
                'file_length', 'acquisition_frequency', &
                'canopy_height', 'displacement_height', 'roughness_length', &
                'master_sonic_manufacturer', 'master_sonic_model', 'master_sonic_height', &
                'master_sonic_wformat', 'master_sonic_wref', 'master_sonic_north_offset', &
                'master_sonic_hpath_length', 'master_sonic_vpath_length', 'master_sonic_tau' /


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
    !> 5.0.0 is the record format: gases, cell measurements and diagnostics are
    !> described by indexed records rather than one column per fixed role. The
    !> GUI migrates a 4.x file on open and saves it in this format.
    character(5), parameter :: MaxSupportedIniVer = '5.1.0'
    type(DateType), parameter :: &
        tsNull= DateType(0, 0, 0, 0, 0)
    type(InstrumentType), parameter :: &
        NullInstrument = InstrumentType('none', 'none', 'none', 'none', 'none', &
        error, error, error, error, &
        error, error, error, error, error, error, error, error, error, error, &
        error, .false., 0, 'none', 'none', &
        'none', 'none', .false., .false., errSwVer)
    type(RawFlagType), parameter :: NullRawFlag = RawFlagType(0, error, .false.)
    type(ColType), parameter :: &
        NullCol = Coltype('none', 'none', 'none', 'none', '', '', '', &
        0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, error, NullInstrument, NullRawFlag, &
        .false., .false., error, &
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
            tau = error, ustar = error, L = error, zL = error, &
            E_at = error, LE_at = error, ET_at = error, H_at = error, &
            tau_at = error, L_at = error, zL_at = error)

    !> Whether the project names a biomet relative humidity column.
    !>
    !> ResolveGasRef needs it - a gas may name the biomet as its moisture
    !> source, and automatic resolution falls back to it - and that function
    !> lives in src_common, where bSetup is not visible. Set once from the ini
    !> rather than asked of biomet%val, which is a per-period value and is not
    !> yet retrieved when DefineE2Set runs.
    logical :: BiometRhConfigured = .false.

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
    !> Where a slower column's real samples sit, and where the spectral window
    !> starts, both needed to rebuild that column's spectrum on its own grid.
    !>
    !> SpecPhase is the offset of the first non-error sample within the
    !> column's own sampling interval, measured in E2Primes rows and recorded
    !> before FixDatasetForSpectra interpolates the gaps away - after that
    !> there is nothing left to measure it from. SpecRowOffset is the row of
    !> E2Primes the spectral set starts at, so the two can be combined.
    !>
    !> Sampling an interpolated column at the wrong phase does not fail, it
    !> quietly attenuates: every sample becomes a linear blend of two real ones.
    integer :: SpecPhase(E2NumVar) = 0
    integer :: SpecRowOffset = 0
    !> Whether the "this gas has no absolute limits" warning has been said.
    !> TestAbsoluteLimits runs once per averaging period, and the condition is a
    !> property of the project rather than of the period, so it is said once.
    logical :: AlLimitsWarned = .false.
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
    integer, parameter :: Npn = 509
    integer, parameter :: Npc = 314
    !> BEGIN GENERATED ProjectRecordOrigins - edit gen_project_tags.py, not this block
    !> Slot origins for the appended gas/cell/diag records. The value
    !> is the index of the FIRST field of record 1, so record i field f
    !> is <origin> + (i-1)*<leap> + f, with f zero-based.
    integer, parameter :: gasNumTag = 33
    integer, parameter :: cellNumTag = 34
    integer, parameter :: diagNumTag = 35
    integer, parameter :: gasRecOriginN = 36
    integer, parameter :: cellRecOriginN = 420
    integer, parameter :: diagRecOriginN = 452
    integer, parameter :: cecNumTag = 468
    integer, parameter :: cecRecOriginN = 469
    integer, parameter :: agcNumTag = 493
    integer, parameter :: agcRecOriginN = 494
    integer, parameter :: gasRecOriginC = 51
    integer, parameter :: cellRecOriginC = 179
    integer, parameter :: diagRecOriginC = 243
    integer, parameter :: cecRecOriginC = 275
    integer, parameter :: agcRecOriginC = 283
    integer, parameter :: rpGasOriginN = 425
    integer, parameter :: rpInstrMaxLackN = 357
    integer, parameter :: rpGasOriginC = 102
    integer, parameter :: fccGasOriginN = 110
    integer, parameter :: fccGasOriginC = 28
    integer, parameter :: gasRecLeapN   = 6
    integer, parameter :: gasRecLeapC   = 2
    integer, parameter :: cellRecLeapN  = 1
    integer, parameter :: cellRecLeapC  = 2
    integer, parameter :: diagRecLeapN  = 1
    integer, parameter :: diagRecLeapC  = 2
    integer, parameter :: rpGasLeapN    = 25
    integer, parameter :: rpGasLeapC    = 3
    integer, parameter :: cecRecLeapN   = 3
    integer, parameter :: cecRecLeapC   = 1
    integer, parameter :: agcRecLeapN   = 1
    integer, parameter :: agcRecLeapC   = 2
    integer, parameter :: fccGasLeapN   = 6
    integer, parameter :: fccGasLeapC   = 1
    !> END GENERATED ProjectRecordOrigins
    logical :: EPPrjNTagFound(Npn)
    logical :: EPPrjCTagFound(Npc)
    type (Numerical) :: EPPrjNTags(Npn)
    type (Text) :: EPPrjCTags(Npc)
    !> BEGIN GENERATED EPPrjNTags - edit gen_project_tags.py, not this block
    data EPPrjNTags(1)%Label / 'binary_nbytes' / &
         EPPrjNTags(2)%Label / 'binary_hnlines' / &
         EPPrjNTags(3)%Label / 'col_ts' / &
         EPPrjNTags(4)%Label / 'cec_singular_band' / &
         EPPrjNTags(5)%Label / 'cec_stationarity_mode' / &
         EPPrjNTags(6)%Label / 'cec_min_flux_sigma' / &
         EPPrjNTags(7)%Label / '' / &
         EPPrjNTags(8)%Label / '' / &
         EPPrjNTags(9)%Label / '' / &
         EPPrjNTags(10)%Label / '' / &
         EPPrjNTags(11)%Label / '' / &
         EPPrjNTags(12)%Label / 'col_air_t' / &
         EPPrjNTags(13)%Label / 'col_air_p' / &
         EPPrjNTags(14)%Label / '' / &
         EPPrjNTags(15)%Label / '' / &
         EPPrjNTags(16)%Label / '' / &
         EPPrjNTags(17)%Label / '' / &
         EPPrjNTags(18)%Label / '' / &
         EPPrjNTags(19)%Label / 'sonic_output_rate' / &
         EPPrjNTags(20)%Label / '' / &
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
         EPPrjNTags(41)%Label / 'gas_1_fluxnet_default' / &
         EPPrjNTags(42)%Label / 'gas_2_col' / &
         EPPrjNTags(43)%Label / 'gas_2_moist' / &
         EPPrjNTags(44)%Label / 'gas_2_cell' / &
         EPPrjNTags(45)%Label / 'gas_2_mw' / &
         EPPrjNTags(46)%Label / 'gas_2_diff' / &
         EPPrjNTags(47)%Label / 'gas_2_fluxnet_default' / &
         EPPrjNTags(48)%Label / 'gas_3_col' / &
         EPPrjNTags(49)%Label / 'gas_3_moist' / &
         EPPrjNTags(50)%Label / 'gas_3_cell' / &
         EPPrjNTags(51)%Label / 'gas_3_mw' / &
         EPPrjNTags(52)%Label / 'gas_3_diff' / &
         EPPrjNTags(53)%Label / 'gas_3_fluxnet_default' / &
         EPPrjNTags(54)%Label / 'gas_4_col' / &
         EPPrjNTags(55)%Label / 'gas_4_moist' / &
         EPPrjNTags(56)%Label / 'gas_4_cell' / &
         EPPrjNTags(57)%Label / 'gas_4_mw' / &
         EPPrjNTags(58)%Label / 'gas_4_diff' / &
         EPPrjNTags(59)%Label / 'gas_4_fluxnet_default' / &
         EPPrjNTags(60)%Label / 'gas_5_col' / &
         EPPrjNTags(61)%Label / 'gas_5_moist' / &
         EPPrjNTags(62)%Label / 'gas_5_cell' / &
         EPPrjNTags(63)%Label / 'gas_5_mw' / &
         EPPrjNTags(64)%Label / 'gas_5_diff' / &
         EPPrjNTags(65)%Label / 'gas_5_fluxnet_default' / &
         EPPrjNTags(66)%Label / 'gas_6_col' / &
         EPPrjNTags(67)%Label / 'gas_6_moist' / &
         EPPrjNTags(68)%Label / 'gas_6_cell' / &
         EPPrjNTags(69)%Label / 'gas_6_mw' / &
         EPPrjNTags(70)%Label / 'gas_6_diff' / &
         EPPrjNTags(71)%Label / 'gas_6_fluxnet_default' / &
         EPPrjNTags(72)%Label / 'gas_7_col' / &
         EPPrjNTags(73)%Label / 'gas_7_moist' / &
         EPPrjNTags(74)%Label / 'gas_7_cell' / &
         EPPrjNTags(75)%Label / 'gas_7_mw' / &
         EPPrjNTags(76)%Label / 'gas_7_diff' / &
         EPPrjNTags(77)%Label / 'gas_7_fluxnet_default' / &
         EPPrjNTags(78)%Label / 'gas_8_col' / &
         EPPrjNTags(79)%Label / 'gas_8_moist' / &
         EPPrjNTags(80)%Label / 'gas_8_cell' / &
         EPPrjNTags(81)%Label / 'gas_8_mw' / &
         EPPrjNTags(82)%Label / 'gas_8_diff' / &
         EPPrjNTags(83)%Label / 'gas_8_fluxnet_default' / &
         EPPrjNTags(84)%Label / 'gas_9_col' / &
         EPPrjNTags(85)%Label / 'gas_9_moist' / &
         EPPrjNTags(86)%Label / 'gas_9_cell' / &
         EPPrjNTags(87)%Label / 'gas_9_mw' / &
         EPPrjNTags(88)%Label / 'gas_9_diff' / &
         EPPrjNTags(89)%Label / 'gas_9_fluxnet_default' / &
         EPPrjNTags(90)%Label / 'gas_10_col' / &
         EPPrjNTags(91)%Label / 'gas_10_moist' / &
         EPPrjNTags(92)%Label / 'gas_10_cell' / &
         EPPrjNTags(93)%Label / 'gas_10_mw' / &
         EPPrjNTags(94)%Label / 'gas_10_diff' / &
         EPPrjNTags(95)%Label / 'gas_10_fluxnet_default' / &
         EPPrjNTags(96)%Label / 'gas_11_col' / &
         EPPrjNTags(97)%Label / 'gas_11_moist' / &
         EPPrjNTags(98)%Label / 'gas_11_cell' / &
         EPPrjNTags(99)%Label / 'gas_11_mw' / &
         EPPrjNTags(100)%Label / 'gas_11_diff' / &
         EPPrjNTags(101)%Label / 'gas_11_fluxnet_default' / &
         EPPrjNTags(102)%Label / 'gas_12_col' / &
         EPPrjNTags(103)%Label / 'gas_12_moist' / &
         EPPrjNTags(104)%Label / 'gas_12_cell' / &
         EPPrjNTags(105)%Label / 'gas_12_mw' / &
         EPPrjNTags(106)%Label / 'gas_12_diff' / &
         EPPrjNTags(107)%Label / 'gas_12_fluxnet_default' / &
         EPPrjNTags(108)%Label / 'gas_13_col' / &
         EPPrjNTags(109)%Label / 'gas_13_moist' / &
         EPPrjNTags(110)%Label / 'gas_13_cell' / &
         EPPrjNTags(111)%Label / 'gas_13_mw' / &
         EPPrjNTags(112)%Label / 'gas_13_diff' / &
         EPPrjNTags(113)%Label / 'gas_13_fluxnet_default' / &
         EPPrjNTags(114)%Label / 'gas_14_col' / &
         EPPrjNTags(115)%Label / 'gas_14_moist' / &
         EPPrjNTags(116)%Label / 'gas_14_cell' / &
         EPPrjNTags(117)%Label / 'gas_14_mw' / &
         EPPrjNTags(118)%Label / 'gas_14_diff' / &
         EPPrjNTags(119)%Label / 'gas_14_fluxnet_default' / &
         EPPrjNTags(120)%Label / 'gas_15_col' / &
         EPPrjNTags(121)%Label / 'gas_15_moist' / &
         EPPrjNTags(122)%Label / 'gas_15_cell' / &
         EPPrjNTags(123)%Label / 'gas_15_mw' / &
         EPPrjNTags(124)%Label / 'gas_15_diff' / &
         EPPrjNTags(125)%Label / 'gas_15_fluxnet_default' / &
         EPPrjNTags(126)%Label / 'gas_16_col' / &
         EPPrjNTags(127)%Label / 'gas_16_moist' / &
         EPPrjNTags(128)%Label / 'gas_16_cell' / &
         EPPrjNTags(129)%Label / 'gas_16_mw' / &
         EPPrjNTags(130)%Label / 'gas_16_diff' / &
         EPPrjNTags(131)%Label / 'gas_16_fluxnet_default' / &
         EPPrjNTags(132)%Label / 'gas_17_col' / &
         EPPrjNTags(133)%Label / 'gas_17_moist' / &
         EPPrjNTags(134)%Label / 'gas_17_cell' / &
         EPPrjNTags(135)%Label / 'gas_17_mw' / &
         EPPrjNTags(136)%Label / 'gas_17_diff' / &
         EPPrjNTags(137)%Label / 'gas_17_fluxnet_default' / &
         EPPrjNTags(138)%Label / 'gas_18_col' / &
         EPPrjNTags(139)%Label / 'gas_18_moist' / &
         EPPrjNTags(140)%Label / 'gas_18_cell' / &
         EPPrjNTags(141)%Label / 'gas_18_mw' / &
         EPPrjNTags(142)%Label / 'gas_18_diff' / &
         EPPrjNTags(143)%Label / 'gas_18_fluxnet_default' / &
         EPPrjNTags(144)%Label / 'gas_19_col' / &
         EPPrjNTags(145)%Label / 'gas_19_moist' / &
         EPPrjNTags(146)%Label / 'gas_19_cell' / &
         EPPrjNTags(147)%Label / 'gas_19_mw' / &
         EPPrjNTags(148)%Label / 'gas_19_diff' / &
         EPPrjNTags(149)%Label / 'gas_19_fluxnet_default' / &
         EPPrjNTags(150)%Label / 'gas_20_col' / &
         EPPrjNTags(151)%Label / 'gas_20_moist' / &
         EPPrjNTags(152)%Label / 'gas_20_cell' / &
         EPPrjNTags(153)%Label / 'gas_20_mw' / &
         EPPrjNTags(154)%Label / 'gas_20_diff' / &
         EPPrjNTags(155)%Label / 'gas_20_fluxnet_default' / &
         EPPrjNTags(156)%Label / 'gas_21_col' / &
         EPPrjNTags(157)%Label / 'gas_21_moist' / &
         EPPrjNTags(158)%Label / 'gas_21_cell' / &
         EPPrjNTags(159)%Label / 'gas_21_mw' / &
         EPPrjNTags(160)%Label / 'gas_21_diff' / &
         EPPrjNTags(161)%Label / 'gas_21_fluxnet_default' / &
         EPPrjNTags(162)%Label / 'gas_22_col' / &
         EPPrjNTags(163)%Label / 'gas_22_moist' / &
         EPPrjNTags(164)%Label / 'gas_22_cell' / &
         EPPrjNTags(165)%Label / 'gas_22_mw' / &
         EPPrjNTags(166)%Label / 'gas_22_diff' / &
         EPPrjNTags(167)%Label / 'gas_22_fluxnet_default' / &
         EPPrjNTags(168)%Label / 'gas_23_col' / &
         EPPrjNTags(169)%Label / 'gas_23_moist' / &
         EPPrjNTags(170)%Label / 'gas_23_cell' / &
         EPPrjNTags(171)%Label / 'gas_23_mw' / &
         EPPrjNTags(172)%Label / 'gas_23_diff' / &
         EPPrjNTags(173)%Label / 'gas_23_fluxnet_default' / &
         EPPrjNTags(174)%Label / 'gas_24_col' / &
         EPPrjNTags(175)%Label / 'gas_24_moist' / &
         EPPrjNTags(176)%Label / 'gas_24_cell' / &
         EPPrjNTags(177)%Label / 'gas_24_mw' / &
         EPPrjNTags(178)%Label / 'gas_24_diff' / &
         EPPrjNTags(179)%Label / 'gas_24_fluxnet_default' / &
         EPPrjNTags(180)%Label / 'gas_25_col' / &
         EPPrjNTags(181)%Label / 'gas_25_moist' / &
         EPPrjNTags(182)%Label / 'gas_25_cell' / &
         EPPrjNTags(183)%Label / 'gas_25_mw' / &
         EPPrjNTags(184)%Label / 'gas_25_diff' / &
         EPPrjNTags(185)%Label / 'gas_25_fluxnet_default' / &
         EPPrjNTags(186)%Label / 'gas_26_col' / &
         EPPrjNTags(187)%Label / 'gas_26_moist' / &
         EPPrjNTags(188)%Label / 'gas_26_cell' / &
         EPPrjNTags(189)%Label / 'gas_26_mw' / &
         EPPrjNTags(190)%Label / 'gas_26_diff' / &
         EPPrjNTags(191)%Label / 'gas_26_fluxnet_default' / &
         EPPrjNTags(192)%Label / 'gas_27_col' / &
         EPPrjNTags(193)%Label / 'gas_27_moist' / &
         EPPrjNTags(194)%Label / 'gas_27_cell' / &
         EPPrjNTags(195)%Label / 'gas_27_mw' / &
         EPPrjNTags(196)%Label / 'gas_27_diff' / &
         EPPrjNTags(197)%Label / 'gas_27_fluxnet_default' / &
         EPPrjNTags(198)%Label / 'gas_28_col' / &
         EPPrjNTags(199)%Label / 'gas_28_moist' / &
         EPPrjNTags(200)%Label / 'gas_28_cell' /
    data EPPrjNTags(201)%Label / 'gas_28_mw' / &
         EPPrjNTags(202)%Label / 'gas_28_diff' / &
         EPPrjNTags(203)%Label / 'gas_28_fluxnet_default' / &
         EPPrjNTags(204)%Label / 'gas_29_col' / &
         EPPrjNTags(205)%Label / 'gas_29_moist' / &
         EPPrjNTags(206)%Label / 'gas_29_cell' / &
         EPPrjNTags(207)%Label / 'gas_29_mw' / &
         EPPrjNTags(208)%Label / 'gas_29_diff' / &
         EPPrjNTags(209)%Label / 'gas_29_fluxnet_default' / &
         EPPrjNTags(210)%Label / 'gas_30_col' / &
         EPPrjNTags(211)%Label / 'gas_30_moist' / &
         EPPrjNTags(212)%Label / 'gas_30_cell' / &
         EPPrjNTags(213)%Label / 'gas_30_mw' / &
         EPPrjNTags(214)%Label / 'gas_30_diff' / &
         EPPrjNTags(215)%Label / 'gas_30_fluxnet_default' / &
         EPPrjNTags(216)%Label / 'gas_31_col' / &
         EPPrjNTags(217)%Label / 'gas_31_moist' / &
         EPPrjNTags(218)%Label / 'gas_31_cell' / &
         EPPrjNTags(219)%Label / 'gas_31_mw' / &
         EPPrjNTags(220)%Label / 'gas_31_diff' / &
         EPPrjNTags(221)%Label / 'gas_31_fluxnet_default' / &
         EPPrjNTags(222)%Label / 'gas_32_col' / &
         EPPrjNTags(223)%Label / 'gas_32_moist' / &
         EPPrjNTags(224)%Label / 'gas_32_cell' / &
         EPPrjNTags(225)%Label / 'gas_32_mw' / &
         EPPrjNTags(226)%Label / 'gas_32_diff' / &
         EPPrjNTags(227)%Label / 'gas_32_fluxnet_default' / &
         EPPrjNTags(228)%Label / 'gas_33_col' / &
         EPPrjNTags(229)%Label / 'gas_33_moist' / &
         EPPrjNTags(230)%Label / 'gas_33_cell' / &
         EPPrjNTags(231)%Label / 'gas_33_mw' / &
         EPPrjNTags(232)%Label / 'gas_33_diff' / &
         EPPrjNTags(233)%Label / 'gas_33_fluxnet_default' / &
         EPPrjNTags(234)%Label / 'gas_34_col' / &
         EPPrjNTags(235)%Label / 'gas_34_moist' / &
         EPPrjNTags(236)%Label / 'gas_34_cell' / &
         EPPrjNTags(237)%Label / 'gas_34_mw' / &
         EPPrjNTags(238)%Label / 'gas_34_diff' / &
         EPPrjNTags(239)%Label / 'gas_34_fluxnet_default' / &
         EPPrjNTags(240)%Label / 'gas_35_col' / &
         EPPrjNTags(241)%Label / 'gas_35_moist' / &
         EPPrjNTags(242)%Label / 'gas_35_cell' / &
         EPPrjNTags(243)%Label / 'gas_35_mw' / &
         EPPrjNTags(244)%Label / 'gas_35_diff' / &
         EPPrjNTags(245)%Label / 'gas_35_fluxnet_default' / &
         EPPrjNTags(246)%Label / 'gas_36_col' / &
         EPPrjNTags(247)%Label / 'gas_36_moist' / &
         EPPrjNTags(248)%Label / 'gas_36_cell' / &
         EPPrjNTags(249)%Label / 'gas_36_mw' / &
         EPPrjNTags(250)%Label / 'gas_36_diff' / &
         EPPrjNTags(251)%Label / 'gas_36_fluxnet_default' / &
         EPPrjNTags(252)%Label / 'gas_37_col' / &
         EPPrjNTags(253)%Label / 'gas_37_moist' / &
         EPPrjNTags(254)%Label / 'gas_37_cell' / &
         EPPrjNTags(255)%Label / 'gas_37_mw' / &
         EPPrjNTags(256)%Label / 'gas_37_diff' / &
         EPPrjNTags(257)%Label / 'gas_37_fluxnet_default' / &
         EPPrjNTags(258)%Label / 'gas_38_col' / &
         EPPrjNTags(259)%Label / 'gas_38_moist' / &
         EPPrjNTags(260)%Label / 'gas_38_cell' / &
         EPPrjNTags(261)%Label / 'gas_38_mw' / &
         EPPrjNTags(262)%Label / 'gas_38_diff' / &
         EPPrjNTags(263)%Label / 'gas_38_fluxnet_default' / &
         EPPrjNTags(264)%Label / 'gas_39_col' / &
         EPPrjNTags(265)%Label / 'gas_39_moist' / &
         EPPrjNTags(266)%Label / 'gas_39_cell' / &
         EPPrjNTags(267)%Label / 'gas_39_mw' / &
         EPPrjNTags(268)%Label / 'gas_39_diff' / &
         EPPrjNTags(269)%Label / 'gas_39_fluxnet_default' / &
         EPPrjNTags(270)%Label / 'gas_40_col' / &
         EPPrjNTags(271)%Label / 'gas_40_moist' / &
         EPPrjNTags(272)%Label / 'gas_40_cell' / &
         EPPrjNTags(273)%Label / 'gas_40_mw' / &
         EPPrjNTags(274)%Label / 'gas_40_diff' / &
         EPPrjNTags(275)%Label / 'gas_40_fluxnet_default' / &
         EPPrjNTags(276)%Label / 'gas_41_col' / &
         EPPrjNTags(277)%Label / 'gas_41_moist' / &
         EPPrjNTags(278)%Label / 'gas_41_cell' / &
         EPPrjNTags(279)%Label / 'gas_41_mw' / &
         EPPrjNTags(280)%Label / 'gas_41_diff' / &
         EPPrjNTags(281)%Label / 'gas_41_fluxnet_default' / &
         EPPrjNTags(282)%Label / 'gas_42_col' / &
         EPPrjNTags(283)%Label / 'gas_42_moist' / &
         EPPrjNTags(284)%Label / 'gas_42_cell' / &
         EPPrjNTags(285)%Label / 'gas_42_mw' / &
         EPPrjNTags(286)%Label / 'gas_42_diff' / &
         EPPrjNTags(287)%Label / 'gas_42_fluxnet_default' / &
         EPPrjNTags(288)%Label / 'gas_43_col' / &
         EPPrjNTags(289)%Label / 'gas_43_moist' / &
         EPPrjNTags(290)%Label / 'gas_43_cell' / &
         EPPrjNTags(291)%Label / 'gas_43_mw' / &
         EPPrjNTags(292)%Label / 'gas_43_diff' / &
         EPPrjNTags(293)%Label / 'gas_43_fluxnet_default' / &
         EPPrjNTags(294)%Label / 'gas_44_col' / &
         EPPrjNTags(295)%Label / 'gas_44_moist' / &
         EPPrjNTags(296)%Label / 'gas_44_cell' / &
         EPPrjNTags(297)%Label / 'gas_44_mw' / &
         EPPrjNTags(298)%Label / 'gas_44_diff' / &
         EPPrjNTags(299)%Label / 'gas_44_fluxnet_default' / &
         EPPrjNTags(300)%Label / 'gas_45_col' / &
         EPPrjNTags(301)%Label / 'gas_45_moist' / &
         EPPrjNTags(302)%Label / 'gas_45_cell' / &
         EPPrjNTags(303)%Label / 'gas_45_mw' / &
         EPPrjNTags(304)%Label / 'gas_45_diff' / &
         EPPrjNTags(305)%Label / 'gas_45_fluxnet_default' / &
         EPPrjNTags(306)%Label / 'gas_46_col' / &
         EPPrjNTags(307)%Label / 'gas_46_moist' / &
         EPPrjNTags(308)%Label / 'gas_46_cell' / &
         EPPrjNTags(309)%Label / 'gas_46_mw' / &
         EPPrjNTags(310)%Label / 'gas_46_diff' / &
         EPPrjNTags(311)%Label / 'gas_46_fluxnet_default' / &
         EPPrjNTags(312)%Label / 'gas_47_col' / &
         EPPrjNTags(313)%Label / 'gas_47_moist' / &
         EPPrjNTags(314)%Label / 'gas_47_cell' / &
         EPPrjNTags(315)%Label / 'gas_47_mw' / &
         EPPrjNTags(316)%Label / 'gas_47_diff' / &
         EPPrjNTags(317)%Label / 'gas_47_fluxnet_default' / &
         EPPrjNTags(318)%Label / 'gas_48_col' / &
         EPPrjNTags(319)%Label / 'gas_48_moist' / &
         EPPrjNTags(320)%Label / 'gas_48_cell' / &
         EPPrjNTags(321)%Label / 'gas_48_mw' / &
         EPPrjNTags(322)%Label / 'gas_48_diff' / &
         EPPrjNTags(323)%Label / 'gas_48_fluxnet_default' / &
         EPPrjNTags(324)%Label / 'gas_49_col' / &
         EPPrjNTags(325)%Label / 'gas_49_moist' / &
         EPPrjNTags(326)%Label / 'gas_49_cell' / &
         EPPrjNTags(327)%Label / 'gas_49_mw' / &
         EPPrjNTags(328)%Label / 'gas_49_diff' / &
         EPPrjNTags(329)%Label / 'gas_49_fluxnet_default' / &
         EPPrjNTags(330)%Label / 'gas_50_col' / &
         EPPrjNTags(331)%Label / 'gas_50_moist' / &
         EPPrjNTags(332)%Label / 'gas_50_cell' / &
         EPPrjNTags(333)%Label / 'gas_50_mw' / &
         EPPrjNTags(334)%Label / 'gas_50_diff' / &
         EPPrjNTags(335)%Label / 'gas_50_fluxnet_default' / &
         EPPrjNTags(336)%Label / 'gas_51_col' / &
         EPPrjNTags(337)%Label / 'gas_51_moist' / &
         EPPrjNTags(338)%Label / 'gas_51_cell' / &
         EPPrjNTags(339)%Label / 'gas_51_mw' / &
         EPPrjNTags(340)%Label / 'gas_51_diff' / &
         EPPrjNTags(341)%Label / 'gas_51_fluxnet_default' / &
         EPPrjNTags(342)%Label / 'gas_52_col' / &
         EPPrjNTags(343)%Label / 'gas_52_moist' / &
         EPPrjNTags(344)%Label / 'gas_52_cell' / &
         EPPrjNTags(345)%Label / 'gas_52_mw' / &
         EPPrjNTags(346)%Label / 'gas_52_diff' / &
         EPPrjNTags(347)%Label / 'gas_52_fluxnet_default' / &
         EPPrjNTags(348)%Label / 'gas_53_col' / &
         EPPrjNTags(349)%Label / 'gas_53_moist' / &
         EPPrjNTags(350)%Label / 'gas_53_cell' / &
         EPPrjNTags(351)%Label / 'gas_53_mw' / &
         EPPrjNTags(352)%Label / 'gas_53_diff' / &
         EPPrjNTags(353)%Label / 'gas_53_fluxnet_default' / &
         EPPrjNTags(354)%Label / 'gas_54_col' / &
         EPPrjNTags(355)%Label / 'gas_54_moist' / &
         EPPrjNTags(356)%Label / 'gas_54_cell' / &
         EPPrjNTags(357)%Label / 'gas_54_mw' / &
         EPPrjNTags(358)%Label / 'gas_54_diff' / &
         EPPrjNTags(359)%Label / 'gas_54_fluxnet_default' / &
         EPPrjNTags(360)%Label / 'gas_55_col' / &
         EPPrjNTags(361)%Label / 'gas_55_moist' / &
         EPPrjNTags(362)%Label / 'gas_55_cell' / &
         EPPrjNTags(363)%Label / 'gas_55_mw' / &
         EPPrjNTags(364)%Label / 'gas_55_diff' / &
         EPPrjNTags(365)%Label / 'gas_55_fluxnet_default' / &
         EPPrjNTags(366)%Label / 'gas_56_col' / &
         EPPrjNTags(367)%Label / 'gas_56_moist' / &
         EPPrjNTags(368)%Label / 'gas_56_cell' / &
         EPPrjNTags(369)%Label / 'gas_56_mw' / &
         EPPrjNTags(370)%Label / 'gas_56_diff' / &
         EPPrjNTags(371)%Label / 'gas_56_fluxnet_default' / &
         EPPrjNTags(372)%Label / 'gas_57_col' / &
         EPPrjNTags(373)%Label / 'gas_57_moist' / &
         EPPrjNTags(374)%Label / 'gas_57_cell' / &
         EPPrjNTags(375)%Label / 'gas_57_mw' / &
         EPPrjNTags(376)%Label / 'gas_57_diff' / &
         EPPrjNTags(377)%Label / 'gas_57_fluxnet_default' / &
         EPPrjNTags(378)%Label / 'gas_58_col' / &
         EPPrjNTags(379)%Label / 'gas_58_moist' / &
         EPPrjNTags(380)%Label / 'gas_58_cell' / &
         EPPrjNTags(381)%Label / 'gas_58_mw' / &
         EPPrjNTags(382)%Label / 'gas_58_diff' / &
         EPPrjNTags(383)%Label / 'gas_58_fluxnet_default' / &
         EPPrjNTags(384)%Label / 'gas_59_col' / &
         EPPrjNTags(385)%Label / 'gas_59_moist' / &
         EPPrjNTags(386)%Label / 'gas_59_cell' / &
         EPPrjNTags(387)%Label / 'gas_59_mw' / &
         EPPrjNTags(388)%Label / 'gas_59_diff' / &
         EPPrjNTags(389)%Label / 'gas_59_fluxnet_default' / &
         EPPrjNTags(390)%Label / 'gas_60_col' / &
         EPPrjNTags(391)%Label / 'gas_60_moist' / &
         EPPrjNTags(392)%Label / 'gas_60_cell' / &
         EPPrjNTags(393)%Label / 'gas_60_mw' / &
         EPPrjNTags(394)%Label / 'gas_60_diff' / &
         EPPrjNTags(395)%Label / 'gas_60_fluxnet_default' / &
         EPPrjNTags(396)%Label / 'gas_61_col' / &
         EPPrjNTags(397)%Label / 'gas_61_moist' / &
         EPPrjNTags(398)%Label / 'gas_61_cell' / &
         EPPrjNTags(399)%Label / 'gas_61_mw' / &
         EPPrjNTags(400)%Label / 'gas_61_diff' /
    data EPPrjNTags(401)%Label / 'gas_61_fluxnet_default' / &
         EPPrjNTags(402)%Label / 'gas_62_col' / &
         EPPrjNTags(403)%Label / 'gas_62_moist' / &
         EPPrjNTags(404)%Label / 'gas_62_cell' / &
         EPPrjNTags(405)%Label / 'gas_62_mw' / &
         EPPrjNTags(406)%Label / 'gas_62_diff' / &
         EPPrjNTags(407)%Label / 'gas_62_fluxnet_default' / &
         EPPrjNTags(408)%Label / 'gas_63_col' / &
         EPPrjNTags(409)%Label / 'gas_63_moist' / &
         EPPrjNTags(410)%Label / 'gas_63_cell' / &
         EPPrjNTags(411)%Label / 'gas_63_mw' / &
         EPPrjNTags(412)%Label / 'gas_63_diff' / &
         EPPrjNTags(413)%Label / 'gas_63_fluxnet_default' / &
         EPPrjNTags(414)%Label / 'gas_64_col' / &
         EPPrjNTags(415)%Label / 'gas_64_moist' / &
         EPPrjNTags(416)%Label / 'gas_64_cell' / &
         EPPrjNTags(417)%Label / 'gas_64_mw' / &
         EPPrjNTags(418)%Label / 'gas_64_diff' / &
         EPPrjNTags(419)%Label / 'gas_64_fluxnet_default' / &
         EPPrjNTags(420)%Label / 'cell_1_col' / &
         EPPrjNTags(421)%Label / 'cell_2_col' / &
         EPPrjNTags(422)%Label / 'cell_3_col' / &
         EPPrjNTags(423)%Label / 'cell_4_col' / &
         EPPrjNTags(424)%Label / 'cell_5_col' / &
         EPPrjNTags(425)%Label / 'cell_6_col' / &
         EPPrjNTags(426)%Label / 'cell_7_col' / &
         EPPrjNTags(427)%Label / 'cell_8_col' / &
         EPPrjNTags(428)%Label / 'cell_9_col' / &
         EPPrjNTags(429)%Label / 'cell_10_col' / &
         EPPrjNTags(430)%Label / 'cell_11_col' / &
         EPPrjNTags(431)%Label / 'cell_12_col' / &
         EPPrjNTags(432)%Label / 'cell_13_col' / &
         EPPrjNTags(433)%Label / 'cell_14_col' / &
         EPPrjNTags(434)%Label / 'cell_15_col' / &
         EPPrjNTags(435)%Label / 'cell_16_col' / &
         EPPrjNTags(436)%Label / 'cell_17_col' / &
         EPPrjNTags(437)%Label / 'cell_18_col' / &
         EPPrjNTags(438)%Label / 'cell_19_col' / &
         EPPrjNTags(439)%Label / 'cell_20_col' / &
         EPPrjNTags(440)%Label / 'cell_21_col' / &
         EPPrjNTags(441)%Label / 'cell_22_col' / &
         EPPrjNTags(442)%Label / 'cell_23_col' / &
         EPPrjNTags(443)%Label / 'cell_24_col' / &
         EPPrjNTags(444)%Label / 'cell_25_col' / &
         EPPrjNTags(445)%Label / 'cell_26_col' / &
         EPPrjNTags(446)%Label / 'cell_27_col' / &
         EPPrjNTags(447)%Label / 'cell_28_col' / &
         EPPrjNTags(448)%Label / 'cell_29_col' / &
         EPPrjNTags(449)%Label / 'cell_30_col' / &
         EPPrjNTags(450)%Label / 'cell_31_col' / &
         EPPrjNTags(451)%Label / 'cell_32_col' / &
         EPPrjNTags(452)%Label / 'diag_1_col' / &
         EPPrjNTags(453)%Label / 'diag_2_col' / &
         EPPrjNTags(454)%Label / 'diag_3_col' / &
         EPPrjNTags(455)%Label / 'diag_4_col' / &
         EPPrjNTags(456)%Label / 'diag_5_col' / &
         EPPrjNTags(457)%Label / 'diag_6_col' / &
         EPPrjNTags(458)%Label / 'diag_7_col' / &
         EPPrjNTags(459)%Label / 'diag_8_col' / &
         EPPrjNTags(460)%Label / 'diag_9_col' / &
         EPPrjNTags(461)%Label / 'diag_10_col' / &
         EPPrjNTags(462)%Label / 'diag_11_col' / &
         EPPrjNTags(463)%Label / 'diag_12_col' / &
         EPPrjNTags(464)%Label / 'diag_13_col' / &
         EPPrjNTags(465)%Label / 'diag_14_col' / &
         EPPrjNTags(466)%Label / 'diag_15_col' / &
         EPPrjNTags(467)%Label / 'diag_16_col' / &
         EPPrjNTags(468)%Label / 'cec_num' / &
         EPPrjNTags(469)%Label / 'cec_1_meth' / &
         EPPrjNTags(470)%Label / 'cec_1_co2' / &
         EPPrjNTags(471)%Label / 'cec_1_h2o' / &
         EPPrjNTags(472)%Label / 'cec_2_meth' / &
         EPPrjNTags(473)%Label / 'cec_2_co2' / &
         EPPrjNTags(474)%Label / 'cec_2_h2o' / &
         EPPrjNTags(475)%Label / 'cec_3_meth' / &
         EPPrjNTags(476)%Label / 'cec_3_co2' / &
         EPPrjNTags(477)%Label / 'cec_3_h2o' / &
         EPPrjNTags(478)%Label / 'cec_4_meth' / &
         EPPrjNTags(479)%Label / 'cec_4_co2' / &
         EPPrjNTags(480)%Label / 'cec_4_h2o' / &
         EPPrjNTags(481)%Label / 'cec_5_meth' / &
         EPPrjNTags(482)%Label / 'cec_5_co2' / &
         EPPrjNTags(483)%Label / 'cec_5_h2o' / &
         EPPrjNTags(484)%Label / 'cec_6_meth' / &
         EPPrjNTags(485)%Label / 'cec_6_co2' / &
         EPPrjNTags(486)%Label / 'cec_6_h2o' / &
         EPPrjNTags(487)%Label / 'cec_7_meth' / &
         EPPrjNTags(488)%Label / 'cec_7_co2' / &
         EPPrjNTags(489)%Label / 'cec_7_h2o' / &
         EPPrjNTags(490)%Label / 'cec_8_meth' / &
         EPPrjNTags(491)%Label / 'cec_8_co2' / &
         EPPrjNTags(492)%Label / 'cec_8_h2o' / &
         EPPrjNTags(493)%Label / 'agc_num' / &
         EPPrjNTags(494)%Label / 'agc_1_col' / &
         EPPrjNTags(495)%Label / 'agc_2_col' / &
         EPPrjNTags(496)%Label / 'agc_3_col' / &
         EPPrjNTags(497)%Label / 'agc_4_col' / &
         EPPrjNTags(498)%Label / 'agc_5_col' / &
         EPPrjNTags(499)%Label / 'agc_6_col' / &
         EPPrjNTags(500)%Label / 'agc_7_col' / &
         EPPrjNTags(501)%Label / 'agc_8_col' / &
         EPPrjNTags(502)%Label / 'agc_9_col' / &
         EPPrjNTags(503)%Label / 'agc_10_col' / &
         EPPrjNTags(504)%Label / 'agc_11_col' / &
         EPPrjNTags(505)%Label / 'agc_12_col' / &
         EPPrjNTags(506)%Label / 'agc_13_col' / &
         EPPrjNTags(507)%Label / 'agc_14_col' / &
         EPPrjNTags(508)%Label / 'agc_15_col' / &
         EPPrjNTags(509)%Label / 'agc_16_col' /
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
         EPPrjCTags(37)%Label / '' / &
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
         EPPrjCTags(274)%Label / 'diag_16_instr' / &
         EPPrjCTags(275)%Label / 'cec_1_extra' / &
         EPPrjCTags(276)%Label / 'cec_2_extra' / &
         EPPrjCTags(277)%Label / 'cec_3_extra' / &
         EPPrjCTags(278)%Label / 'cec_4_extra' / &
         EPPrjCTags(279)%Label / 'cec_5_extra' / &
         EPPrjCTags(280)%Label / 'cec_6_extra' / &
         EPPrjCTags(281)%Label / 'cec_7_extra' / &
         EPPrjCTags(282)%Label / 'cec_8_extra' / &
         EPPrjCTags(283)%Label / 'agc_1_var' / &
         EPPrjCTags(284)%Label / 'agc_1_instr' / &
         EPPrjCTags(285)%Label / 'agc_2_var' / &
         EPPrjCTags(286)%Label / 'agc_2_instr' / &
         EPPrjCTags(287)%Label / 'agc_3_var' / &
         EPPrjCTags(288)%Label / 'agc_3_instr' / &
         EPPrjCTags(289)%Label / 'agc_4_var' / &
         EPPrjCTags(290)%Label / 'agc_4_instr' / &
         EPPrjCTags(291)%Label / 'agc_5_var' / &
         EPPrjCTags(292)%Label / 'agc_5_instr' / &
         EPPrjCTags(293)%Label / 'agc_6_var' / &
         EPPrjCTags(294)%Label / 'agc_6_instr' / &
         EPPrjCTags(295)%Label / 'agc_7_var' / &
         EPPrjCTags(296)%Label / 'agc_7_instr' / &
         EPPrjCTags(297)%Label / 'agc_8_var' / &
         EPPrjCTags(298)%Label / 'agc_8_instr' / &
         EPPrjCTags(299)%Label / 'agc_9_var' / &
         EPPrjCTags(300)%Label / 'agc_9_instr' / &
         EPPrjCTags(301)%Label / 'agc_10_var' / &
         EPPrjCTags(302)%Label / 'agc_10_instr' / &
         EPPrjCTags(303)%Label / 'agc_11_var' / &
         EPPrjCTags(304)%Label / 'agc_11_instr' / &
         EPPrjCTags(305)%Label / 'agc_12_var' / &
         EPPrjCTags(306)%Label / 'agc_12_instr' / &
         EPPrjCTags(307)%Label / 'agc_13_var' / &
         EPPrjCTags(308)%Label / 'agc_13_instr' / &
         EPPrjCTags(309)%Label / 'agc_14_var' / &
         EPPrjCTags(310)%Label / 'agc_14_instr' / &
         EPPrjCTags(311)%Label / 'agc_15_var' / &
         EPPrjCTags(312)%Label / 'agc_15_instr' / &
         EPPrjCTags(313)%Label / 'agc_16_var' / &
         EPPrjCTags(314)%Label / 'agc_16_instr' /
    !> END GENERATED EPPrjCTags

    !> tags of the metadata file created by GHG software
    integer, parameter :: Nan = 1929
    integer, parameter :: Nac = 1490
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
         ANTags(23)%Label / 'instr_1_ac_freq' / &
         ANTags(24)%Label / 'instr_1_integrates' / &
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
         ANTags(38)%Label / 'instr_2_ac_freq' / &
         ANTags(39)%Label / 'instr_2_integrates' / &
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
         ANTags(53)%Label / 'instr_3_ac_freq' / &
         ANTags(54)%Label / 'instr_3_integrates' / &
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
         ANTags(68)%Label / 'instr_4_ac_freq' / &
         ANTags(69)%Label / 'instr_4_integrates' / &
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
         ANTags(83)%Label / 'instr_5_ac_freq' / &
         ANTags(84)%Label / 'instr_5_integrates' / &
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
         ANTags(98)%Label / 'instr_6_ac_freq' / &
         ANTags(99)%Label / 'instr_6_integrates' / &
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
         ANTags(113)%Label / 'instr_7_ac_freq' / &
         ANTags(114)%Label / 'instr_7_integrates' / &
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
         ANTags(128)%Label / 'instr_8_ac_freq' / &
         ANTags(129)%Label / 'instr_8_integrates' / &
         ANTags(130)%Label / 'col_1_min_value' / &
         ANTags(131)%Label / 'col_1_max_value' / &
         ANTags(132)%Label / 'col_1_a_value' / &
         ANTags(133)%Label / 'col_1_b_value' / &
         ANTags(134)%Label / 'col_1_nom_timelag' / &
         ANTags(135)%Label / 'col_1_min_timelag' / &
         ANTags(136)%Label / 'col_1_max_timelag' / &
         ANTags(137)%Label / 'col_1_flag_threshold' / &
         ANTags(138)%Label / 'col_1_error_value' / &
         ANTags(139)%Label / 'col_2_min_value' / &
         ANTags(140)%Label / 'col_2_max_value' / &
         ANTags(141)%Label / 'col_2_a_value' / &
         ANTags(142)%Label / 'col_2_b_value' / &
         ANTags(143)%Label / 'col_2_nom_timelag' / &
         ANTags(144)%Label / 'col_2_min_timelag' / &
         ANTags(145)%Label / 'col_2_max_timelag' / &
         ANTags(146)%Label / 'col_2_flag_threshold' / &
         ANTags(147)%Label / 'col_2_error_value' / &
         ANTags(148)%Label / 'col_3_min_value' / &
         ANTags(149)%Label / 'col_3_max_value' / &
         ANTags(150)%Label / 'col_3_a_value' / &
         ANTags(151)%Label / 'col_3_b_value' / &
         ANTags(152)%Label / 'col_3_nom_timelag' / &
         ANTags(153)%Label / 'col_3_min_timelag' / &
         ANTags(154)%Label / 'col_3_max_timelag' / &
         ANTags(155)%Label / 'col_3_flag_threshold' / &
         ANTags(156)%Label / 'col_3_error_value' / &
         ANTags(157)%Label / 'col_4_min_value' / &
         ANTags(158)%Label / 'col_4_max_value' / &
         ANTags(159)%Label / 'col_4_a_value' / &
         ANTags(160)%Label / 'col_4_b_value' / &
         ANTags(161)%Label / 'col_4_nom_timelag' / &
         ANTags(162)%Label / 'col_4_min_timelag' / &
         ANTags(163)%Label / 'col_4_max_timelag' / &
         ANTags(164)%Label / 'col_4_flag_threshold' / &
         ANTags(165)%Label / 'col_4_error_value' / &
         ANTags(166)%Label / 'col_5_min_value' / &
         ANTags(167)%Label / 'col_5_max_value' / &
         ANTags(168)%Label / 'col_5_a_value' / &
         ANTags(169)%Label / 'col_5_b_value' / &
         ANTags(170)%Label / 'col_5_nom_timelag' / &
         ANTags(171)%Label / 'col_5_min_timelag' / &
         ANTags(172)%Label / 'col_5_max_timelag' / &
         ANTags(173)%Label / 'col_5_flag_threshold' / &
         ANTags(174)%Label / 'col_5_error_value' / &
         ANTags(175)%Label / 'col_6_min_value' / &
         ANTags(176)%Label / 'col_6_max_value' / &
         ANTags(177)%Label / 'col_6_a_value' / &
         ANTags(178)%Label / 'col_6_b_value' / &
         ANTags(179)%Label / 'col_6_nom_timelag' / &
         ANTags(180)%Label / 'col_6_min_timelag' / &
         ANTags(181)%Label / 'col_6_max_timelag' / &
         ANTags(182)%Label / 'col_6_flag_threshold' / &
         ANTags(183)%Label / 'col_6_error_value' / &
         ANTags(184)%Label / 'col_7_min_value' / &
         ANTags(185)%Label / 'col_7_max_value' / &
         ANTags(186)%Label / 'col_7_a_value' / &
         ANTags(187)%Label / 'col_7_b_value' / &
         ANTags(188)%Label / 'col_7_nom_timelag' / &
         ANTags(189)%Label / 'col_7_min_timelag' / &
         ANTags(190)%Label / 'col_7_max_timelag' / &
         ANTags(191)%Label / 'col_7_flag_threshold' / &
         ANTags(192)%Label / 'col_7_error_value' / &
         ANTags(193)%Label / 'col_8_min_value' / &
         ANTags(194)%Label / 'col_8_max_value' / &
         ANTags(195)%Label / 'col_8_a_value' / &
         ANTags(196)%Label / 'col_8_b_value' / &
         ANTags(197)%Label / 'col_8_nom_timelag' / &
         ANTags(198)%Label / 'col_8_min_timelag' / &
         ANTags(199)%Label / 'col_8_max_timelag' / &
         ANTags(200)%Label / 'col_8_flag_threshold' /
    data ANTags(201)%Label / 'col_8_error_value' / &
         ANTags(202)%Label / 'col_9_min_value' / &
         ANTags(203)%Label / 'col_9_max_value' / &
         ANTags(204)%Label / 'col_9_a_value' / &
         ANTags(205)%Label / 'col_9_b_value' / &
         ANTags(206)%Label / 'col_9_nom_timelag' / &
         ANTags(207)%Label / 'col_9_min_timelag' / &
         ANTags(208)%Label / 'col_9_max_timelag' / &
         ANTags(209)%Label / 'col_9_flag_threshold' / &
         ANTags(210)%Label / 'col_9_error_value' / &
         ANTags(211)%Label / 'col_10_min_value' / &
         ANTags(212)%Label / 'col_10_max_value' / &
         ANTags(213)%Label / 'col_10_a_value' / &
         ANTags(214)%Label / 'col_10_b_value' / &
         ANTags(215)%Label / 'col_10_nom_timelag' / &
         ANTags(216)%Label / 'col_10_min_timelag' / &
         ANTags(217)%Label / 'col_10_max_timelag' / &
         ANTags(218)%Label / 'col_10_flag_threshold' / &
         ANTags(219)%Label / 'col_10_error_value' / &
         ANTags(220)%Label / 'col_11_min_value' / &
         ANTags(221)%Label / 'col_11_max_value' / &
         ANTags(222)%Label / 'col_11_a_value' / &
         ANTags(223)%Label / 'col_11_b_value' / &
         ANTags(224)%Label / 'col_11_nom_timelag' / &
         ANTags(225)%Label / 'col_11_min_timelag' / &
         ANTags(226)%Label / 'col_11_max_timelag' / &
         ANTags(227)%Label / 'col_11_flag_threshold' / &
         ANTags(228)%Label / 'col_11_error_value' / &
         ANTags(229)%Label / 'col_12_min_value' / &
         ANTags(230)%Label / 'col_12_max_value' / &
         ANTags(231)%Label / 'col_12_a_value' / &
         ANTags(232)%Label / 'col_12_b_value' / &
         ANTags(233)%Label / 'col_12_nom_timelag' / &
         ANTags(234)%Label / 'col_12_min_timelag' / &
         ANTags(235)%Label / 'col_12_max_timelag' / &
         ANTags(236)%Label / 'col_12_flag_threshold' / &
         ANTags(237)%Label / 'col_12_error_value' / &
         ANTags(238)%Label / 'col_13_min_value' / &
         ANTags(239)%Label / 'col_13_max_value' / &
         ANTags(240)%Label / 'col_13_a_value' / &
         ANTags(241)%Label / 'col_13_b_value' / &
         ANTags(242)%Label / 'col_13_nom_timelag' / &
         ANTags(243)%Label / 'col_13_min_timelag' / &
         ANTags(244)%Label / 'col_13_max_timelag' / &
         ANTags(245)%Label / 'col_13_flag_threshold' / &
         ANTags(246)%Label / 'col_13_error_value' / &
         ANTags(247)%Label / 'col_14_min_value' / &
         ANTags(248)%Label / 'col_14_max_value' / &
         ANTags(249)%Label / 'col_14_a_value' / &
         ANTags(250)%Label / 'col_14_b_value' / &
         ANTags(251)%Label / 'col_14_nom_timelag' / &
         ANTags(252)%Label / 'col_14_min_timelag' / &
         ANTags(253)%Label / 'col_14_max_timelag' / &
         ANTags(254)%Label / 'col_14_flag_threshold' / &
         ANTags(255)%Label / 'col_14_error_value' / &
         ANTags(256)%Label / 'col_15_min_value' / &
         ANTags(257)%Label / 'col_15_max_value' / &
         ANTags(258)%Label / 'col_15_a_value' / &
         ANTags(259)%Label / 'col_15_b_value' / &
         ANTags(260)%Label / 'col_15_nom_timelag' / &
         ANTags(261)%Label / 'col_15_min_timelag' / &
         ANTags(262)%Label / 'col_15_max_timelag' / &
         ANTags(263)%Label / 'col_15_flag_threshold' / &
         ANTags(264)%Label / 'col_15_error_value' / &
         ANTags(265)%Label / 'col_16_min_value' / &
         ANTags(266)%Label / 'col_16_max_value' / &
         ANTags(267)%Label / 'col_16_a_value' / &
         ANTags(268)%Label / 'col_16_b_value' / &
         ANTags(269)%Label / 'col_16_nom_timelag' / &
         ANTags(270)%Label / 'col_16_min_timelag' / &
         ANTags(271)%Label / 'col_16_max_timelag' / &
         ANTags(272)%Label / 'col_16_flag_threshold' / &
         ANTags(273)%Label / 'col_16_error_value' / &
         ANTags(274)%Label / 'col_17_min_value' / &
         ANTags(275)%Label / 'col_17_max_value' / &
         ANTags(276)%Label / 'col_17_a_value' / &
         ANTags(277)%Label / 'col_17_b_value' / &
         ANTags(278)%Label / 'col_17_nom_timelag' / &
         ANTags(279)%Label / 'col_17_min_timelag' / &
         ANTags(280)%Label / 'col_17_max_timelag' / &
         ANTags(281)%Label / 'col_17_flag_threshold' / &
         ANTags(282)%Label / 'col_17_error_value' / &
         ANTags(283)%Label / 'col_18_min_value' / &
         ANTags(284)%Label / 'col_18_max_value' / &
         ANTags(285)%Label / 'col_18_a_value' / &
         ANTags(286)%Label / 'col_18_b_value' / &
         ANTags(287)%Label / 'col_18_nom_timelag' / &
         ANTags(288)%Label / 'col_18_min_timelag' / &
         ANTags(289)%Label / 'col_18_max_timelag' / &
         ANTags(290)%Label / 'col_18_flag_threshold' / &
         ANTags(291)%Label / 'col_18_error_value' / &
         ANTags(292)%Label / 'col_19_min_value' / &
         ANTags(293)%Label / 'col_19_max_value' / &
         ANTags(294)%Label / 'col_19_a_value' / &
         ANTags(295)%Label / 'col_19_b_value' / &
         ANTags(296)%Label / 'col_19_nom_timelag' / &
         ANTags(297)%Label / 'col_19_min_timelag' / &
         ANTags(298)%Label / 'col_19_max_timelag' / &
         ANTags(299)%Label / 'col_19_flag_threshold' / &
         ANTags(300)%Label / 'col_19_error_value' / &
         ANTags(301)%Label / 'col_20_min_value' / &
         ANTags(302)%Label / 'col_20_max_value' / &
         ANTags(303)%Label / 'col_20_a_value' / &
         ANTags(304)%Label / 'col_20_b_value' / &
         ANTags(305)%Label / 'col_20_nom_timelag' / &
         ANTags(306)%Label / 'col_20_min_timelag' / &
         ANTags(307)%Label / 'col_20_max_timelag' / &
         ANTags(308)%Label / 'col_20_flag_threshold' / &
         ANTags(309)%Label / 'col_20_error_value' / &
         ANTags(310)%Label / 'col_21_min_value' / &
         ANTags(311)%Label / 'col_21_max_value' / &
         ANTags(312)%Label / 'col_21_a_value' / &
         ANTags(313)%Label / 'col_21_b_value' / &
         ANTags(314)%Label / 'col_21_nom_timelag' / &
         ANTags(315)%Label / 'col_21_min_timelag' / &
         ANTags(316)%Label / 'col_21_max_timelag' / &
         ANTags(317)%Label / 'col_21_flag_threshold' / &
         ANTags(318)%Label / 'col_21_error_value' / &
         ANTags(319)%Label / 'col_22_min_value' / &
         ANTags(320)%Label / 'col_22_max_value' / &
         ANTags(321)%Label / 'col_22_a_value' / &
         ANTags(322)%Label / 'col_22_b_value' / &
         ANTags(323)%Label / 'col_22_nom_timelag' / &
         ANTags(324)%Label / 'col_22_min_timelag' / &
         ANTags(325)%Label / 'col_22_max_timelag' / &
         ANTags(326)%Label / 'col_22_flag_threshold' / &
         ANTags(327)%Label / 'col_22_error_value' / &
         ANTags(328)%Label / 'col_23_min_value' / &
         ANTags(329)%Label / 'col_23_max_value' / &
         ANTags(330)%Label / 'col_23_a_value' / &
         ANTags(331)%Label / 'col_23_b_value' / &
         ANTags(332)%Label / 'col_23_nom_timelag' / &
         ANTags(333)%Label / 'col_23_min_timelag' / &
         ANTags(334)%Label / 'col_23_max_timelag' / &
         ANTags(335)%Label / 'col_23_flag_threshold' / &
         ANTags(336)%Label / 'col_23_error_value' / &
         ANTags(337)%Label / 'col_24_min_value' / &
         ANTags(338)%Label / 'col_24_max_value' / &
         ANTags(339)%Label / 'col_24_a_value' / &
         ANTags(340)%Label / 'col_24_b_value' / &
         ANTags(341)%Label / 'col_24_nom_timelag' / &
         ANTags(342)%Label / 'col_24_min_timelag' / &
         ANTags(343)%Label / 'col_24_max_timelag' / &
         ANTags(344)%Label / 'col_24_flag_threshold' / &
         ANTags(345)%Label / 'col_24_error_value' / &
         ANTags(346)%Label / 'col_25_min_value' / &
         ANTags(347)%Label / 'col_25_max_value' / &
         ANTags(348)%Label / 'col_25_a_value' / &
         ANTags(349)%Label / 'col_25_b_value' / &
         ANTags(350)%Label / 'col_25_nom_timelag' / &
         ANTags(351)%Label / 'col_25_min_timelag' / &
         ANTags(352)%Label / 'col_25_max_timelag' / &
         ANTags(353)%Label / 'col_25_flag_threshold' / &
         ANTags(354)%Label / 'col_25_error_value' / &
         ANTags(355)%Label / 'col_26_min_value' / &
         ANTags(356)%Label / 'col_26_max_value' / &
         ANTags(357)%Label / 'col_26_a_value' / &
         ANTags(358)%Label / 'col_26_b_value' / &
         ANTags(359)%Label / 'col_26_nom_timelag' / &
         ANTags(360)%Label / 'col_26_min_timelag' / &
         ANTags(361)%Label / 'col_26_max_timelag' / &
         ANTags(362)%Label / 'col_26_flag_threshold' / &
         ANTags(363)%Label / 'col_26_error_value' / &
         ANTags(364)%Label / 'col_27_min_value' / &
         ANTags(365)%Label / 'col_27_max_value' / &
         ANTags(366)%Label / 'col_27_a_value' / &
         ANTags(367)%Label / 'col_27_b_value' / &
         ANTags(368)%Label / 'col_27_nom_timelag' / &
         ANTags(369)%Label / 'col_27_min_timelag' / &
         ANTags(370)%Label / 'col_27_max_timelag' / &
         ANTags(371)%Label / 'col_27_flag_threshold' / &
         ANTags(372)%Label / 'col_27_error_value' / &
         ANTags(373)%Label / 'col_28_min_value' / &
         ANTags(374)%Label / 'col_28_max_value' / &
         ANTags(375)%Label / 'col_28_a_value' / &
         ANTags(376)%Label / 'col_28_b_value' / &
         ANTags(377)%Label / 'col_28_nom_timelag' / &
         ANTags(378)%Label / 'col_28_min_timelag' / &
         ANTags(379)%Label / 'col_28_max_timelag' / &
         ANTags(380)%Label / 'col_28_flag_threshold' / &
         ANTags(381)%Label / 'col_28_error_value' / &
         ANTags(382)%Label / 'col_29_min_value' / &
         ANTags(383)%Label / 'col_29_max_value' / &
         ANTags(384)%Label / 'col_29_a_value' / &
         ANTags(385)%Label / 'col_29_b_value' / &
         ANTags(386)%Label / 'col_29_nom_timelag' / &
         ANTags(387)%Label / 'col_29_min_timelag' / &
         ANTags(388)%Label / 'col_29_max_timelag' / &
         ANTags(389)%Label / 'col_29_flag_threshold' / &
         ANTags(390)%Label / 'col_29_error_value' / &
         ANTags(391)%Label / 'col_30_min_value' / &
         ANTags(392)%Label / 'col_30_max_value' / &
         ANTags(393)%Label / 'col_30_a_value' / &
         ANTags(394)%Label / 'col_30_b_value' / &
         ANTags(395)%Label / 'col_30_nom_timelag' / &
         ANTags(396)%Label / 'col_30_min_timelag' / &
         ANTags(397)%Label / 'col_30_max_timelag' / &
         ANTags(398)%Label / 'col_30_flag_threshold' / &
         ANTags(399)%Label / 'col_30_error_value' / &
         ANTags(400)%Label / 'col_31_min_value' /
    data ANTags(401)%Label / 'col_31_max_value' / &
         ANTags(402)%Label / 'col_31_a_value' / &
         ANTags(403)%Label / 'col_31_b_value' / &
         ANTags(404)%Label / 'col_31_nom_timelag' / &
         ANTags(405)%Label / 'col_31_min_timelag' / &
         ANTags(406)%Label / 'col_31_max_timelag' / &
         ANTags(407)%Label / 'col_31_flag_threshold' / &
         ANTags(408)%Label / 'col_31_error_value' / &
         ANTags(409)%Label / 'col_32_min_value' / &
         ANTags(410)%Label / 'col_32_max_value' / &
         ANTags(411)%Label / 'col_32_a_value' / &
         ANTags(412)%Label / 'col_32_b_value' / &
         ANTags(413)%Label / 'col_32_nom_timelag' / &
         ANTags(414)%Label / 'col_32_min_timelag' / &
         ANTags(415)%Label / 'col_32_max_timelag' / &
         ANTags(416)%Label / 'col_32_flag_threshold' / &
         ANTags(417)%Label / 'col_32_error_value' / &
         ANTags(418)%Label / 'col_33_min_value' / &
         ANTags(419)%Label / 'col_33_max_value' / &
         ANTags(420)%Label / 'col_33_a_value' / &
         ANTags(421)%Label / 'col_33_b_value' / &
         ANTags(422)%Label / 'col_33_nom_timelag' / &
         ANTags(423)%Label / 'col_33_min_timelag' / &
         ANTags(424)%Label / 'col_33_max_timelag' / &
         ANTags(425)%Label / 'col_33_flag_threshold' / &
         ANTags(426)%Label / 'col_33_error_value' / &
         ANTags(427)%Label / 'col_34_min_value' / &
         ANTags(428)%Label / 'col_34_max_value' / &
         ANTags(429)%Label / 'col_34_a_value' / &
         ANTags(430)%Label / 'col_34_b_value' / &
         ANTags(431)%Label / 'col_34_nom_timelag' / &
         ANTags(432)%Label / 'col_34_min_timelag' / &
         ANTags(433)%Label / 'col_34_max_timelag' / &
         ANTags(434)%Label / 'col_34_flag_threshold' / &
         ANTags(435)%Label / 'col_34_error_value' / &
         ANTags(436)%Label / 'col_35_min_value' / &
         ANTags(437)%Label / 'col_35_max_value' / &
         ANTags(438)%Label / 'col_35_a_value' / &
         ANTags(439)%Label / 'col_35_b_value' / &
         ANTags(440)%Label / 'col_35_nom_timelag' / &
         ANTags(441)%Label / 'col_35_min_timelag' / &
         ANTags(442)%Label / 'col_35_max_timelag' / &
         ANTags(443)%Label / 'col_35_flag_threshold' / &
         ANTags(444)%Label / 'col_35_error_value' / &
         ANTags(445)%Label / 'col_36_min_value' / &
         ANTags(446)%Label / 'col_36_max_value' / &
         ANTags(447)%Label / 'col_36_a_value' / &
         ANTags(448)%Label / 'col_36_b_value' / &
         ANTags(449)%Label / 'col_36_nom_timelag' / &
         ANTags(450)%Label / 'col_36_min_timelag' / &
         ANTags(451)%Label / 'col_36_max_timelag' / &
         ANTags(452)%Label / 'col_36_flag_threshold' / &
         ANTags(453)%Label / 'col_36_error_value' / &
         ANTags(454)%Label / 'col_37_min_value' / &
         ANTags(455)%Label / 'col_37_max_value' / &
         ANTags(456)%Label / 'col_37_a_value' / &
         ANTags(457)%Label / 'col_37_b_value' / &
         ANTags(458)%Label / 'col_37_nom_timelag' / &
         ANTags(459)%Label / 'col_37_min_timelag' / &
         ANTags(460)%Label / 'col_37_max_timelag' / &
         ANTags(461)%Label / 'col_37_flag_threshold' / &
         ANTags(462)%Label / 'col_37_error_value' / &
         ANTags(463)%Label / 'col_38_min_value' / &
         ANTags(464)%Label / 'col_38_max_value' / &
         ANTags(465)%Label / 'col_38_a_value' / &
         ANTags(466)%Label / 'col_38_b_value' / &
         ANTags(467)%Label / 'col_38_nom_timelag' / &
         ANTags(468)%Label / 'col_38_min_timelag' / &
         ANTags(469)%Label / 'col_38_max_timelag' / &
         ANTags(470)%Label / 'col_38_flag_threshold' / &
         ANTags(471)%Label / 'col_38_error_value' / &
         ANTags(472)%Label / 'col_39_min_value' / &
         ANTags(473)%Label / 'col_39_max_value' / &
         ANTags(474)%Label / 'col_39_a_value' / &
         ANTags(475)%Label / 'col_39_b_value' / &
         ANTags(476)%Label / 'col_39_nom_timelag' / &
         ANTags(477)%Label / 'col_39_min_timelag' / &
         ANTags(478)%Label / 'col_39_max_timelag' / &
         ANTags(479)%Label / 'col_39_flag_threshold' / &
         ANTags(480)%Label / 'col_39_error_value' / &
         ANTags(481)%Label / 'col_40_min_value' / &
         ANTags(482)%Label / 'col_40_max_value' / &
         ANTags(483)%Label / 'col_40_a_value' / &
         ANTags(484)%Label / 'col_40_b_value' / &
         ANTags(485)%Label / 'col_40_nom_timelag' / &
         ANTags(486)%Label / 'col_40_min_timelag' / &
         ANTags(487)%Label / 'col_40_max_timelag' / &
         ANTags(488)%Label / 'col_40_flag_threshold' / &
         ANTags(489)%Label / 'col_40_error_value' / &
         ANTags(490)%Label / 'col_41_min_value' / &
         ANTags(491)%Label / 'col_41_max_value' / &
         ANTags(492)%Label / 'col_41_a_value' / &
         ANTags(493)%Label / 'col_41_b_value' / &
         ANTags(494)%Label / 'col_41_nom_timelag' / &
         ANTags(495)%Label / 'col_41_min_timelag' / &
         ANTags(496)%Label / 'col_41_max_timelag' / &
         ANTags(497)%Label / 'col_41_flag_threshold' / &
         ANTags(498)%Label / 'col_41_error_value' / &
         ANTags(499)%Label / 'col_42_min_value' / &
         ANTags(500)%Label / 'col_42_max_value' / &
         ANTags(501)%Label / 'col_42_a_value' / &
         ANTags(502)%Label / 'col_42_b_value' / &
         ANTags(503)%Label / 'col_42_nom_timelag' / &
         ANTags(504)%Label / 'col_42_min_timelag' / &
         ANTags(505)%Label / 'col_42_max_timelag' / &
         ANTags(506)%Label / 'col_42_flag_threshold' / &
         ANTags(507)%Label / 'col_42_error_value' / &
         ANTags(508)%Label / 'col_43_min_value' / &
         ANTags(509)%Label / 'col_43_max_value' / &
         ANTags(510)%Label / 'col_43_a_value' / &
         ANTags(511)%Label / 'col_43_b_value' / &
         ANTags(512)%Label / 'col_43_nom_timelag' / &
         ANTags(513)%Label / 'col_43_min_timelag' / &
         ANTags(514)%Label / 'col_43_max_timelag' / &
         ANTags(515)%Label / 'col_43_flag_threshold' / &
         ANTags(516)%Label / 'col_43_error_value' / &
         ANTags(517)%Label / 'col_44_min_value' / &
         ANTags(518)%Label / 'col_44_max_value' / &
         ANTags(519)%Label / 'col_44_a_value' / &
         ANTags(520)%Label / 'col_44_b_value' / &
         ANTags(521)%Label / 'col_44_nom_timelag' / &
         ANTags(522)%Label / 'col_44_min_timelag' / &
         ANTags(523)%Label / 'col_44_max_timelag' / &
         ANTags(524)%Label / 'col_44_flag_threshold' / &
         ANTags(525)%Label / 'col_44_error_value' / &
         ANTags(526)%Label / 'col_45_min_value' / &
         ANTags(527)%Label / 'col_45_max_value' / &
         ANTags(528)%Label / 'col_45_a_value' / &
         ANTags(529)%Label / 'col_45_b_value' / &
         ANTags(530)%Label / 'col_45_nom_timelag' / &
         ANTags(531)%Label / 'col_45_min_timelag' / &
         ANTags(532)%Label / 'col_45_max_timelag' / &
         ANTags(533)%Label / 'col_45_flag_threshold' / &
         ANTags(534)%Label / 'col_45_error_value' / &
         ANTags(535)%Label / 'col_46_min_value' / &
         ANTags(536)%Label / 'col_46_max_value' / &
         ANTags(537)%Label / 'col_46_a_value' / &
         ANTags(538)%Label / 'col_46_b_value' / &
         ANTags(539)%Label / 'col_46_nom_timelag' / &
         ANTags(540)%Label / 'col_46_min_timelag' / &
         ANTags(541)%Label / 'col_46_max_timelag' / &
         ANTags(542)%Label / 'col_46_flag_threshold' / &
         ANTags(543)%Label / 'col_46_error_value' / &
         ANTags(544)%Label / 'col_47_min_value' / &
         ANTags(545)%Label / 'col_47_max_value' / &
         ANTags(546)%Label / 'col_47_a_value' / &
         ANTags(547)%Label / 'col_47_b_value' / &
         ANTags(548)%Label / 'col_47_nom_timelag' / &
         ANTags(549)%Label / 'col_47_min_timelag' / &
         ANTags(550)%Label / 'col_47_max_timelag' / &
         ANTags(551)%Label / 'col_47_flag_threshold' / &
         ANTags(552)%Label / 'col_47_error_value' / &
         ANTags(553)%Label / 'col_48_min_value' / &
         ANTags(554)%Label / 'col_48_max_value' / &
         ANTags(555)%Label / 'col_48_a_value' / &
         ANTags(556)%Label / 'col_48_b_value' / &
         ANTags(557)%Label / 'col_48_nom_timelag' / &
         ANTags(558)%Label / 'col_48_min_timelag' / &
         ANTags(559)%Label / 'col_48_max_timelag' / &
         ANTags(560)%Label / 'col_48_flag_threshold' / &
         ANTags(561)%Label / 'col_48_error_value' / &
         ANTags(562)%Label / 'col_49_min_value' / &
         ANTags(563)%Label / 'col_49_max_value' / &
         ANTags(564)%Label / 'col_49_a_value' / &
         ANTags(565)%Label / 'col_49_b_value' / &
         ANTags(566)%Label / 'col_49_nom_timelag' / &
         ANTags(567)%Label / 'col_49_min_timelag' / &
         ANTags(568)%Label / 'col_49_max_timelag' / &
         ANTags(569)%Label / 'col_49_flag_threshold' / &
         ANTags(570)%Label / 'col_49_error_value' / &
         ANTags(571)%Label / 'col_50_min_value' / &
         ANTags(572)%Label / 'col_50_max_value' / &
         ANTags(573)%Label / 'col_50_a_value' / &
         ANTags(574)%Label / 'col_50_b_value' / &
         ANTags(575)%Label / 'col_50_nom_timelag' / &
         ANTags(576)%Label / 'col_50_min_timelag' / &
         ANTags(577)%Label / 'col_50_max_timelag' / &
         ANTags(578)%Label / 'col_50_flag_threshold' / &
         ANTags(579)%Label / 'col_50_error_value' / &
         ANTags(580)%Label / 'col_51_min_value' / &
         ANTags(581)%Label / 'col_51_max_value' / &
         ANTags(582)%Label / 'col_51_a_value' / &
         ANTags(583)%Label / 'col_51_b_value' / &
         ANTags(584)%Label / 'col_51_nom_timelag' / &
         ANTags(585)%Label / 'col_51_min_timelag' / &
         ANTags(586)%Label / 'col_51_max_timelag' / &
         ANTags(587)%Label / 'col_51_flag_threshold' / &
         ANTags(588)%Label / 'col_51_error_value' / &
         ANTags(589)%Label / 'col_52_min_value' / &
         ANTags(590)%Label / 'col_52_max_value' / &
         ANTags(591)%Label / 'col_52_a_value' / &
         ANTags(592)%Label / 'col_52_b_value' / &
         ANTags(593)%Label / 'col_52_nom_timelag' / &
         ANTags(594)%Label / 'col_52_min_timelag' / &
         ANTags(595)%Label / 'col_52_max_timelag' / &
         ANTags(596)%Label / 'col_52_flag_threshold' / &
         ANTags(597)%Label / 'col_52_error_value' / &
         ANTags(598)%Label / 'col_53_min_value' / &
         ANTags(599)%Label / 'col_53_max_value' / &
         ANTags(600)%Label / 'col_53_a_value' /
    data ANTags(601)%Label / 'col_53_b_value' / &
         ANTags(602)%Label / 'col_53_nom_timelag' / &
         ANTags(603)%Label / 'col_53_min_timelag' / &
         ANTags(604)%Label / 'col_53_max_timelag' / &
         ANTags(605)%Label / 'col_53_flag_threshold' / &
         ANTags(606)%Label / 'col_53_error_value' / &
         ANTags(607)%Label / 'col_54_min_value' / &
         ANTags(608)%Label / 'col_54_max_value' / &
         ANTags(609)%Label / 'col_54_a_value' / &
         ANTags(610)%Label / 'col_54_b_value' / &
         ANTags(611)%Label / 'col_54_nom_timelag' / &
         ANTags(612)%Label / 'col_54_min_timelag' / &
         ANTags(613)%Label / 'col_54_max_timelag' / &
         ANTags(614)%Label / 'col_54_flag_threshold' / &
         ANTags(615)%Label / 'col_54_error_value' / &
         ANTags(616)%Label / 'col_55_min_value' / &
         ANTags(617)%Label / 'col_55_max_value' / &
         ANTags(618)%Label / 'col_55_a_value' / &
         ANTags(619)%Label / 'col_55_b_value' / &
         ANTags(620)%Label / 'col_55_nom_timelag' / &
         ANTags(621)%Label / 'col_55_min_timelag' / &
         ANTags(622)%Label / 'col_55_max_timelag' / &
         ANTags(623)%Label / 'col_55_flag_threshold' / &
         ANTags(624)%Label / 'col_55_error_value' / &
         ANTags(625)%Label / 'col_56_min_value' / &
         ANTags(626)%Label / 'col_56_max_value' / &
         ANTags(627)%Label / 'col_56_a_value' / &
         ANTags(628)%Label / 'col_56_b_value' / &
         ANTags(629)%Label / 'col_56_nom_timelag' / &
         ANTags(630)%Label / 'col_56_min_timelag' / &
         ANTags(631)%Label / 'col_56_max_timelag' / &
         ANTags(632)%Label / 'col_56_flag_threshold' / &
         ANTags(633)%Label / 'col_56_error_value' / &
         ANTags(634)%Label / 'col_57_min_value' / &
         ANTags(635)%Label / 'col_57_max_value' / &
         ANTags(636)%Label / 'col_57_a_value' / &
         ANTags(637)%Label / 'col_57_b_value' / &
         ANTags(638)%Label / 'col_57_nom_timelag' / &
         ANTags(639)%Label / 'col_57_min_timelag' / &
         ANTags(640)%Label / 'col_57_max_timelag' / &
         ANTags(641)%Label / 'col_57_flag_threshold' / &
         ANTags(642)%Label / 'col_57_error_value' / &
         ANTags(643)%Label / 'col_58_min_value' / &
         ANTags(644)%Label / 'col_58_max_value' / &
         ANTags(645)%Label / 'col_58_a_value' / &
         ANTags(646)%Label / 'col_58_b_value' / &
         ANTags(647)%Label / 'col_58_nom_timelag' / &
         ANTags(648)%Label / 'col_58_min_timelag' / &
         ANTags(649)%Label / 'col_58_max_timelag' / &
         ANTags(650)%Label / 'col_58_flag_threshold' / &
         ANTags(651)%Label / 'col_58_error_value' / &
         ANTags(652)%Label / 'col_59_min_value' / &
         ANTags(653)%Label / 'col_59_max_value' / &
         ANTags(654)%Label / 'col_59_a_value' / &
         ANTags(655)%Label / 'col_59_b_value' / &
         ANTags(656)%Label / 'col_59_nom_timelag' / &
         ANTags(657)%Label / 'col_59_min_timelag' / &
         ANTags(658)%Label / 'col_59_max_timelag' / &
         ANTags(659)%Label / 'col_59_flag_threshold' / &
         ANTags(660)%Label / 'col_59_error_value' / &
         ANTags(661)%Label / 'col_60_min_value' / &
         ANTags(662)%Label / 'col_60_max_value' / &
         ANTags(663)%Label / 'col_60_a_value' / &
         ANTags(664)%Label / 'col_60_b_value' / &
         ANTags(665)%Label / 'col_60_nom_timelag' / &
         ANTags(666)%Label / 'col_60_min_timelag' / &
         ANTags(667)%Label / 'col_60_max_timelag' / &
         ANTags(668)%Label / 'col_60_flag_threshold' / &
         ANTags(669)%Label / 'col_60_error_value' / &
         ANTags(670)%Label / 'col_61_min_value' / &
         ANTags(671)%Label / 'col_61_max_value' / &
         ANTags(672)%Label / 'col_61_a_value' / &
         ANTags(673)%Label / 'col_61_b_value' / &
         ANTags(674)%Label / 'col_61_nom_timelag' / &
         ANTags(675)%Label / 'col_61_min_timelag' / &
         ANTags(676)%Label / 'col_61_max_timelag' / &
         ANTags(677)%Label / 'col_61_flag_threshold' / &
         ANTags(678)%Label / 'col_61_error_value' / &
         ANTags(679)%Label / 'col_62_min_value' / &
         ANTags(680)%Label / 'col_62_max_value' / &
         ANTags(681)%Label / 'col_62_a_value' / &
         ANTags(682)%Label / 'col_62_b_value' / &
         ANTags(683)%Label / 'col_62_nom_timelag' / &
         ANTags(684)%Label / 'col_62_min_timelag' / &
         ANTags(685)%Label / 'col_62_max_timelag' / &
         ANTags(686)%Label / 'col_62_flag_threshold' / &
         ANTags(687)%Label / 'col_62_error_value' / &
         ANTags(688)%Label / 'col_63_min_value' / &
         ANTags(689)%Label / 'col_63_max_value' / &
         ANTags(690)%Label / 'col_63_a_value' / &
         ANTags(691)%Label / 'col_63_b_value' / &
         ANTags(692)%Label / 'col_63_nom_timelag' / &
         ANTags(693)%Label / 'col_63_min_timelag' / &
         ANTags(694)%Label / 'col_63_max_timelag' / &
         ANTags(695)%Label / 'col_63_flag_threshold' / &
         ANTags(696)%Label / 'col_63_error_value' / &
         ANTags(697)%Label / 'col_64_min_value' / &
         ANTags(698)%Label / 'col_64_max_value' / &
         ANTags(699)%Label / 'col_64_a_value' / &
         ANTags(700)%Label / 'col_64_b_value' / &
         ANTags(701)%Label / 'col_64_nom_timelag' / &
         ANTags(702)%Label / 'col_64_min_timelag' / &
         ANTags(703)%Label / 'col_64_max_timelag' / &
         ANTags(704)%Label / 'col_64_flag_threshold' / &
         ANTags(705)%Label / 'col_64_error_value' / &
         ANTags(706)%Label / 'col_65_min_value' / &
         ANTags(707)%Label / 'col_65_max_value' / &
         ANTags(708)%Label / 'col_65_a_value' / &
         ANTags(709)%Label / 'col_65_b_value' / &
         ANTags(710)%Label / 'col_65_nom_timelag' / &
         ANTags(711)%Label / 'col_65_min_timelag' / &
         ANTags(712)%Label / 'col_65_max_timelag' / &
         ANTags(713)%Label / 'col_65_flag_threshold' / &
         ANTags(714)%Label / 'col_65_error_value' / &
         ANTags(715)%Label / 'col_66_min_value' / &
         ANTags(716)%Label / 'col_66_max_value' / &
         ANTags(717)%Label / 'col_66_a_value' / &
         ANTags(718)%Label / 'col_66_b_value' / &
         ANTags(719)%Label / 'col_66_nom_timelag' / &
         ANTags(720)%Label / 'col_66_min_timelag' / &
         ANTags(721)%Label / 'col_66_max_timelag' / &
         ANTags(722)%Label / 'col_66_flag_threshold' / &
         ANTags(723)%Label / 'col_66_error_value' / &
         ANTags(724)%Label / 'col_67_min_value' / &
         ANTags(725)%Label / 'col_67_max_value' / &
         ANTags(726)%Label / 'col_67_a_value' / &
         ANTags(727)%Label / 'col_67_b_value' / &
         ANTags(728)%Label / 'col_67_nom_timelag' / &
         ANTags(729)%Label / 'col_67_min_timelag' / &
         ANTags(730)%Label / 'col_67_max_timelag' / &
         ANTags(731)%Label / 'col_67_flag_threshold' / &
         ANTags(732)%Label / 'col_67_error_value' / &
         ANTags(733)%Label / 'col_68_min_value' / &
         ANTags(734)%Label / 'col_68_max_value' / &
         ANTags(735)%Label / 'col_68_a_value' / &
         ANTags(736)%Label / 'col_68_b_value' / &
         ANTags(737)%Label / 'col_68_nom_timelag' / &
         ANTags(738)%Label / 'col_68_min_timelag' / &
         ANTags(739)%Label / 'col_68_max_timelag' / &
         ANTags(740)%Label / 'col_68_flag_threshold' / &
         ANTags(741)%Label / 'col_68_error_value' / &
         ANTags(742)%Label / 'col_69_min_value' / &
         ANTags(743)%Label / 'col_69_max_value' / &
         ANTags(744)%Label / 'col_69_a_value' / &
         ANTags(745)%Label / 'col_69_b_value' / &
         ANTags(746)%Label / 'col_69_nom_timelag' / &
         ANTags(747)%Label / 'col_69_min_timelag' / &
         ANTags(748)%Label / 'col_69_max_timelag' / &
         ANTags(749)%Label / 'col_69_flag_threshold' / &
         ANTags(750)%Label / 'col_69_error_value' / &
         ANTags(751)%Label / 'col_70_min_value' / &
         ANTags(752)%Label / 'col_70_max_value' / &
         ANTags(753)%Label / 'col_70_a_value' / &
         ANTags(754)%Label / 'col_70_b_value' / &
         ANTags(755)%Label / 'col_70_nom_timelag' / &
         ANTags(756)%Label / 'col_70_min_timelag' / &
         ANTags(757)%Label / 'col_70_max_timelag' / &
         ANTags(758)%Label / 'col_70_flag_threshold' / &
         ANTags(759)%Label / 'col_70_error_value' / &
         ANTags(760)%Label / 'col_71_min_value' / &
         ANTags(761)%Label / 'col_71_max_value' / &
         ANTags(762)%Label / 'col_71_a_value' / &
         ANTags(763)%Label / 'col_71_b_value' / &
         ANTags(764)%Label / 'col_71_nom_timelag' / &
         ANTags(765)%Label / 'col_71_min_timelag' / &
         ANTags(766)%Label / 'col_71_max_timelag' / &
         ANTags(767)%Label / 'col_71_flag_threshold' / &
         ANTags(768)%Label / 'col_71_error_value' / &
         ANTags(769)%Label / 'col_72_min_value' / &
         ANTags(770)%Label / 'col_72_max_value' / &
         ANTags(771)%Label / 'col_72_a_value' / &
         ANTags(772)%Label / 'col_72_b_value' / &
         ANTags(773)%Label / 'col_72_nom_timelag' / &
         ANTags(774)%Label / 'col_72_min_timelag' / &
         ANTags(775)%Label / 'col_72_max_timelag' / &
         ANTags(776)%Label / 'col_72_flag_threshold' / &
         ANTags(777)%Label / 'col_72_error_value' / &
         ANTags(778)%Label / 'col_73_min_value' / &
         ANTags(779)%Label / 'col_73_max_value' / &
         ANTags(780)%Label / 'col_73_a_value' / &
         ANTags(781)%Label / 'col_73_b_value' / &
         ANTags(782)%Label / 'col_73_nom_timelag' / &
         ANTags(783)%Label / 'col_73_min_timelag' / &
         ANTags(784)%Label / 'col_73_max_timelag' / &
         ANTags(785)%Label / 'col_73_flag_threshold' / &
         ANTags(786)%Label / 'col_73_error_value' / &
         ANTags(787)%Label / 'col_74_min_value' / &
         ANTags(788)%Label / 'col_74_max_value' / &
         ANTags(789)%Label / 'col_74_a_value' / &
         ANTags(790)%Label / 'col_74_b_value' / &
         ANTags(791)%Label / 'col_74_nom_timelag' / &
         ANTags(792)%Label / 'col_74_min_timelag' / &
         ANTags(793)%Label / 'col_74_max_timelag' / &
         ANTags(794)%Label / 'col_74_flag_threshold' / &
         ANTags(795)%Label / 'col_74_error_value' / &
         ANTags(796)%Label / 'col_75_min_value' / &
         ANTags(797)%Label / 'col_75_max_value' / &
         ANTags(798)%Label / 'col_75_a_value' / &
         ANTags(799)%Label / 'col_75_b_value' / &
         ANTags(800)%Label / 'col_75_nom_timelag' /
    data ANTags(801)%Label / 'col_75_min_timelag' / &
         ANTags(802)%Label / 'col_75_max_timelag' / &
         ANTags(803)%Label / 'col_75_flag_threshold' / &
         ANTags(804)%Label / 'col_75_error_value' / &
         ANTags(805)%Label / 'col_76_min_value' / &
         ANTags(806)%Label / 'col_76_max_value' / &
         ANTags(807)%Label / 'col_76_a_value' / &
         ANTags(808)%Label / 'col_76_b_value' / &
         ANTags(809)%Label / 'col_76_nom_timelag' / &
         ANTags(810)%Label / 'col_76_min_timelag' / &
         ANTags(811)%Label / 'col_76_max_timelag' / &
         ANTags(812)%Label / 'col_76_flag_threshold' / &
         ANTags(813)%Label / 'col_76_error_value' / &
         ANTags(814)%Label / 'col_77_min_value' / &
         ANTags(815)%Label / 'col_77_max_value' / &
         ANTags(816)%Label / 'col_77_a_value' / &
         ANTags(817)%Label / 'col_77_b_value' / &
         ANTags(818)%Label / 'col_77_nom_timelag' / &
         ANTags(819)%Label / 'col_77_min_timelag' / &
         ANTags(820)%Label / 'col_77_max_timelag' / &
         ANTags(821)%Label / 'col_77_flag_threshold' / &
         ANTags(822)%Label / 'col_77_error_value' / &
         ANTags(823)%Label / 'col_78_min_value' / &
         ANTags(824)%Label / 'col_78_max_value' / &
         ANTags(825)%Label / 'col_78_a_value' / &
         ANTags(826)%Label / 'col_78_b_value' / &
         ANTags(827)%Label / 'col_78_nom_timelag' / &
         ANTags(828)%Label / 'col_78_min_timelag' / &
         ANTags(829)%Label / 'col_78_max_timelag' / &
         ANTags(830)%Label / 'col_78_flag_threshold' / &
         ANTags(831)%Label / 'col_78_error_value' / &
         ANTags(832)%Label / 'col_79_min_value' / &
         ANTags(833)%Label / 'col_79_max_value' / &
         ANTags(834)%Label / 'col_79_a_value' / &
         ANTags(835)%Label / 'col_79_b_value' / &
         ANTags(836)%Label / 'col_79_nom_timelag' / &
         ANTags(837)%Label / 'col_79_min_timelag' / &
         ANTags(838)%Label / 'col_79_max_timelag' / &
         ANTags(839)%Label / 'col_79_flag_threshold' / &
         ANTags(840)%Label / 'col_79_error_value' / &
         ANTags(841)%Label / 'col_80_min_value' / &
         ANTags(842)%Label / 'col_80_max_value' / &
         ANTags(843)%Label / 'col_80_a_value' / &
         ANTags(844)%Label / 'col_80_b_value' / &
         ANTags(845)%Label / 'col_80_nom_timelag' / &
         ANTags(846)%Label / 'col_80_min_timelag' / &
         ANTags(847)%Label / 'col_80_max_timelag' / &
         ANTags(848)%Label / 'col_80_flag_threshold' / &
         ANTags(849)%Label / 'col_80_error_value' / &
         ANTags(850)%Label / 'col_81_min_value' / &
         ANTags(851)%Label / 'col_81_max_value' / &
         ANTags(852)%Label / 'col_81_a_value' / &
         ANTags(853)%Label / 'col_81_b_value' / &
         ANTags(854)%Label / 'col_81_nom_timelag' / &
         ANTags(855)%Label / 'col_81_min_timelag' / &
         ANTags(856)%Label / 'col_81_max_timelag' / &
         ANTags(857)%Label / 'col_81_flag_threshold' / &
         ANTags(858)%Label / 'col_81_error_value' / &
         ANTags(859)%Label / 'col_82_min_value' / &
         ANTags(860)%Label / 'col_82_max_value' / &
         ANTags(861)%Label / 'col_82_a_value' / &
         ANTags(862)%Label / 'col_82_b_value' / &
         ANTags(863)%Label / 'col_82_nom_timelag' / &
         ANTags(864)%Label / 'col_82_min_timelag' / &
         ANTags(865)%Label / 'col_82_max_timelag' / &
         ANTags(866)%Label / 'col_82_flag_threshold' / &
         ANTags(867)%Label / 'col_82_error_value' / &
         ANTags(868)%Label / 'col_83_min_value' / &
         ANTags(869)%Label / 'col_83_max_value' / &
         ANTags(870)%Label / 'col_83_a_value' / &
         ANTags(871)%Label / 'col_83_b_value' / &
         ANTags(872)%Label / 'col_83_nom_timelag' / &
         ANTags(873)%Label / 'col_83_min_timelag' / &
         ANTags(874)%Label / 'col_83_max_timelag' / &
         ANTags(875)%Label / 'col_83_flag_threshold' / &
         ANTags(876)%Label / 'col_83_error_value' / &
         ANTags(877)%Label / 'col_84_min_value' / &
         ANTags(878)%Label / 'col_84_max_value' / &
         ANTags(879)%Label / 'col_84_a_value' / &
         ANTags(880)%Label / 'col_84_b_value' / &
         ANTags(881)%Label / 'col_84_nom_timelag' / &
         ANTags(882)%Label / 'col_84_min_timelag' / &
         ANTags(883)%Label / 'col_84_max_timelag' / &
         ANTags(884)%Label / 'col_84_flag_threshold' / &
         ANTags(885)%Label / 'col_84_error_value' / &
         ANTags(886)%Label / 'col_85_min_value' / &
         ANTags(887)%Label / 'col_85_max_value' / &
         ANTags(888)%Label / 'col_85_a_value' / &
         ANTags(889)%Label / 'col_85_b_value' / &
         ANTags(890)%Label / 'col_85_nom_timelag' / &
         ANTags(891)%Label / 'col_85_min_timelag' / &
         ANTags(892)%Label / 'col_85_max_timelag' / &
         ANTags(893)%Label / 'col_85_flag_threshold' / &
         ANTags(894)%Label / 'col_85_error_value' / &
         ANTags(895)%Label / 'col_86_min_value' / &
         ANTags(896)%Label / 'col_86_max_value' / &
         ANTags(897)%Label / 'col_86_a_value' / &
         ANTags(898)%Label / 'col_86_b_value' / &
         ANTags(899)%Label / 'col_86_nom_timelag' / &
         ANTags(900)%Label / 'col_86_min_timelag' / &
         ANTags(901)%Label / 'col_86_max_timelag' / &
         ANTags(902)%Label / 'col_86_flag_threshold' / &
         ANTags(903)%Label / 'col_86_error_value' / &
         ANTags(904)%Label / 'col_87_min_value' / &
         ANTags(905)%Label / 'col_87_max_value' / &
         ANTags(906)%Label / 'col_87_a_value' / &
         ANTags(907)%Label / 'col_87_b_value' / &
         ANTags(908)%Label / 'col_87_nom_timelag' / &
         ANTags(909)%Label / 'col_87_min_timelag' / &
         ANTags(910)%Label / 'col_87_max_timelag' / &
         ANTags(911)%Label / 'col_87_flag_threshold' / &
         ANTags(912)%Label / 'col_87_error_value' / &
         ANTags(913)%Label / 'col_88_min_value' / &
         ANTags(914)%Label / 'col_88_max_value' / &
         ANTags(915)%Label / 'col_88_a_value' / &
         ANTags(916)%Label / 'col_88_b_value' / &
         ANTags(917)%Label / 'col_88_nom_timelag' / &
         ANTags(918)%Label / 'col_88_min_timelag' / &
         ANTags(919)%Label / 'col_88_max_timelag' / &
         ANTags(920)%Label / 'col_88_flag_threshold' / &
         ANTags(921)%Label / 'col_88_error_value' / &
         ANTags(922)%Label / 'col_89_min_value' / &
         ANTags(923)%Label / 'col_89_max_value' / &
         ANTags(924)%Label / 'col_89_a_value' / &
         ANTags(925)%Label / 'col_89_b_value' / &
         ANTags(926)%Label / 'col_89_nom_timelag' / &
         ANTags(927)%Label / 'col_89_min_timelag' / &
         ANTags(928)%Label / 'col_89_max_timelag' / &
         ANTags(929)%Label / 'col_89_flag_threshold' / &
         ANTags(930)%Label / 'col_89_error_value' / &
         ANTags(931)%Label / 'col_90_min_value' / &
         ANTags(932)%Label / 'col_90_max_value' / &
         ANTags(933)%Label / 'col_90_a_value' / &
         ANTags(934)%Label / 'col_90_b_value' / &
         ANTags(935)%Label / 'col_90_nom_timelag' / &
         ANTags(936)%Label / 'col_90_min_timelag' / &
         ANTags(937)%Label / 'col_90_max_timelag' / &
         ANTags(938)%Label / 'col_90_flag_threshold' / &
         ANTags(939)%Label / 'col_90_error_value' / &
         ANTags(940)%Label / 'col_91_min_value' / &
         ANTags(941)%Label / 'col_91_max_value' / &
         ANTags(942)%Label / 'col_91_a_value' / &
         ANTags(943)%Label / 'col_91_b_value' / &
         ANTags(944)%Label / 'col_91_nom_timelag' / &
         ANTags(945)%Label / 'col_91_min_timelag' / &
         ANTags(946)%Label / 'col_91_max_timelag' / &
         ANTags(947)%Label / 'col_91_flag_threshold' / &
         ANTags(948)%Label / 'col_91_error_value' / &
         ANTags(949)%Label / 'col_92_min_value' / &
         ANTags(950)%Label / 'col_92_max_value' / &
         ANTags(951)%Label / 'col_92_a_value' / &
         ANTags(952)%Label / 'col_92_b_value' / &
         ANTags(953)%Label / 'col_92_nom_timelag' / &
         ANTags(954)%Label / 'col_92_min_timelag' / &
         ANTags(955)%Label / 'col_92_max_timelag' / &
         ANTags(956)%Label / 'col_92_flag_threshold' / &
         ANTags(957)%Label / 'col_92_error_value' / &
         ANTags(958)%Label / 'col_93_min_value' / &
         ANTags(959)%Label / 'col_93_max_value' / &
         ANTags(960)%Label / 'col_93_a_value' / &
         ANTags(961)%Label / 'col_93_b_value' / &
         ANTags(962)%Label / 'col_93_nom_timelag' / &
         ANTags(963)%Label / 'col_93_min_timelag' / &
         ANTags(964)%Label / 'col_93_max_timelag' / &
         ANTags(965)%Label / 'col_93_flag_threshold' / &
         ANTags(966)%Label / 'col_93_error_value' / &
         ANTags(967)%Label / 'col_94_min_value' / &
         ANTags(968)%Label / 'col_94_max_value' / &
         ANTags(969)%Label / 'col_94_a_value' / &
         ANTags(970)%Label / 'col_94_b_value' / &
         ANTags(971)%Label / 'col_94_nom_timelag' / &
         ANTags(972)%Label / 'col_94_min_timelag' / &
         ANTags(973)%Label / 'col_94_max_timelag' / &
         ANTags(974)%Label / 'col_94_flag_threshold' / &
         ANTags(975)%Label / 'col_94_error_value' / &
         ANTags(976)%Label / 'col_95_min_value' / &
         ANTags(977)%Label / 'col_95_max_value' / &
         ANTags(978)%Label / 'col_95_a_value' / &
         ANTags(979)%Label / 'col_95_b_value' / &
         ANTags(980)%Label / 'col_95_nom_timelag' / &
         ANTags(981)%Label / 'col_95_min_timelag' / &
         ANTags(982)%Label / 'col_95_max_timelag' / &
         ANTags(983)%Label / 'col_95_flag_threshold' / &
         ANTags(984)%Label / 'col_95_error_value' / &
         ANTags(985)%Label / 'col_96_min_value' / &
         ANTags(986)%Label / 'col_96_max_value' / &
         ANTags(987)%Label / 'col_96_a_value' / &
         ANTags(988)%Label / 'col_96_b_value' / &
         ANTags(989)%Label / 'col_96_nom_timelag' / &
         ANTags(990)%Label / 'col_96_min_timelag' / &
         ANTags(991)%Label / 'col_96_max_timelag' / &
         ANTags(992)%Label / 'col_96_flag_threshold' / &
         ANTags(993)%Label / 'col_96_error_value' / &
         ANTags(994)%Label / 'col_97_min_value' / &
         ANTags(995)%Label / 'col_97_max_value' / &
         ANTags(996)%Label / 'col_97_a_value' / &
         ANTags(997)%Label / 'col_97_b_value' / &
         ANTags(998)%Label / 'col_97_nom_timelag' / &
         ANTags(999)%Label / 'col_97_min_timelag' / &
         ANTags(1000)%Label / 'col_97_max_timelag' /
    data ANTags(1001)%Label / 'col_97_flag_threshold' / &
         ANTags(1002)%Label / 'col_97_error_value' / &
         ANTags(1003)%Label / 'col_98_min_value' / &
         ANTags(1004)%Label / 'col_98_max_value' / &
         ANTags(1005)%Label / 'col_98_a_value' / &
         ANTags(1006)%Label / 'col_98_b_value' / &
         ANTags(1007)%Label / 'col_98_nom_timelag' / &
         ANTags(1008)%Label / 'col_98_min_timelag' / &
         ANTags(1009)%Label / 'col_98_max_timelag' / &
         ANTags(1010)%Label / 'col_98_flag_threshold' / &
         ANTags(1011)%Label / 'col_98_error_value' / &
         ANTags(1012)%Label / 'col_99_min_value' / &
         ANTags(1013)%Label / 'col_99_max_value' / &
         ANTags(1014)%Label / 'col_99_a_value' / &
         ANTags(1015)%Label / 'col_99_b_value' / &
         ANTags(1016)%Label / 'col_99_nom_timelag' / &
         ANTags(1017)%Label / 'col_99_min_timelag' / &
         ANTags(1018)%Label / 'col_99_max_timelag' / &
         ANTags(1019)%Label / 'col_99_flag_threshold' / &
         ANTags(1020)%Label / 'col_99_error_value' / &
         ANTags(1021)%Label / 'col_100_min_value' / &
         ANTags(1022)%Label / 'col_100_max_value' / &
         ANTags(1023)%Label / 'col_100_a_value' / &
         ANTags(1024)%Label / 'col_100_b_value' / &
         ANTags(1025)%Label / 'col_100_nom_timelag' / &
         ANTags(1026)%Label / 'col_100_min_timelag' / &
         ANTags(1027)%Label / 'col_100_max_timelag' / &
         ANTags(1028)%Label / 'col_100_flag_threshold' / &
         ANTags(1029)%Label / 'col_100_error_value' / &
         ANTags(1030)%Label / 'col_101_min_value' / &
         ANTags(1031)%Label / 'col_101_max_value' / &
         ANTags(1032)%Label / 'col_101_a_value' / &
         ANTags(1033)%Label / 'col_101_b_value' / &
         ANTags(1034)%Label / 'col_101_nom_timelag' / &
         ANTags(1035)%Label / 'col_101_min_timelag' / &
         ANTags(1036)%Label / 'col_101_max_timelag' / &
         ANTags(1037)%Label / 'col_101_flag_threshold' / &
         ANTags(1038)%Label / 'col_101_error_value' / &
         ANTags(1039)%Label / 'col_102_min_value' / &
         ANTags(1040)%Label / 'col_102_max_value' / &
         ANTags(1041)%Label / 'col_102_a_value' / &
         ANTags(1042)%Label / 'col_102_b_value' / &
         ANTags(1043)%Label / 'col_102_nom_timelag' / &
         ANTags(1044)%Label / 'col_102_min_timelag' / &
         ANTags(1045)%Label / 'col_102_max_timelag' / &
         ANTags(1046)%Label / 'col_102_flag_threshold' / &
         ANTags(1047)%Label / 'col_102_error_value' / &
         ANTags(1048)%Label / 'col_103_min_value' / &
         ANTags(1049)%Label / 'col_103_max_value' / &
         ANTags(1050)%Label / 'col_103_a_value' / &
         ANTags(1051)%Label / 'col_103_b_value' / &
         ANTags(1052)%Label / 'col_103_nom_timelag' / &
         ANTags(1053)%Label / 'col_103_min_timelag' / &
         ANTags(1054)%Label / 'col_103_max_timelag' / &
         ANTags(1055)%Label / 'col_103_flag_threshold' / &
         ANTags(1056)%Label / 'col_103_error_value' / &
         ANTags(1057)%Label / 'col_104_min_value' / &
         ANTags(1058)%Label / 'col_104_max_value' / &
         ANTags(1059)%Label / 'col_104_a_value' / &
         ANTags(1060)%Label / 'col_104_b_value' / &
         ANTags(1061)%Label / 'col_104_nom_timelag' / &
         ANTags(1062)%Label / 'col_104_min_timelag' / &
         ANTags(1063)%Label / 'col_104_max_timelag' / &
         ANTags(1064)%Label / 'col_104_flag_threshold' / &
         ANTags(1065)%Label / 'col_104_error_value' / &
         ANTags(1066)%Label / 'col_105_min_value' / &
         ANTags(1067)%Label / 'col_105_max_value' / &
         ANTags(1068)%Label / 'col_105_a_value' / &
         ANTags(1069)%Label / 'col_105_b_value' / &
         ANTags(1070)%Label / 'col_105_nom_timelag' / &
         ANTags(1071)%Label / 'col_105_min_timelag' / &
         ANTags(1072)%Label / 'col_105_max_timelag' / &
         ANTags(1073)%Label / 'col_105_flag_threshold' / &
         ANTags(1074)%Label / 'col_105_error_value' / &
         ANTags(1075)%Label / 'col_106_min_value' / &
         ANTags(1076)%Label / 'col_106_max_value' / &
         ANTags(1077)%Label / 'col_106_a_value' / &
         ANTags(1078)%Label / 'col_106_b_value' / &
         ANTags(1079)%Label / 'col_106_nom_timelag' / &
         ANTags(1080)%Label / 'col_106_min_timelag' / &
         ANTags(1081)%Label / 'col_106_max_timelag' / &
         ANTags(1082)%Label / 'col_106_flag_threshold' / &
         ANTags(1083)%Label / 'col_106_error_value' / &
         ANTags(1084)%Label / 'col_107_min_value' / &
         ANTags(1085)%Label / 'col_107_max_value' / &
         ANTags(1086)%Label / 'col_107_a_value' / &
         ANTags(1087)%Label / 'col_107_b_value' / &
         ANTags(1088)%Label / 'col_107_nom_timelag' / &
         ANTags(1089)%Label / 'col_107_min_timelag' / &
         ANTags(1090)%Label / 'col_107_max_timelag' / &
         ANTags(1091)%Label / 'col_107_flag_threshold' / &
         ANTags(1092)%Label / 'col_107_error_value' / &
         ANTags(1093)%Label / 'col_108_min_value' / &
         ANTags(1094)%Label / 'col_108_max_value' / &
         ANTags(1095)%Label / 'col_108_a_value' / &
         ANTags(1096)%Label / 'col_108_b_value' / &
         ANTags(1097)%Label / 'col_108_nom_timelag' / &
         ANTags(1098)%Label / 'col_108_min_timelag' / &
         ANTags(1099)%Label / 'col_108_max_timelag' / &
         ANTags(1100)%Label / 'col_108_flag_threshold' / &
         ANTags(1101)%Label / 'col_108_error_value' / &
         ANTags(1102)%Label / 'col_109_min_value' / &
         ANTags(1103)%Label / 'col_109_max_value' / &
         ANTags(1104)%Label / 'col_109_a_value' / &
         ANTags(1105)%Label / 'col_109_b_value' / &
         ANTags(1106)%Label / 'col_109_nom_timelag' / &
         ANTags(1107)%Label / 'col_109_min_timelag' / &
         ANTags(1108)%Label / 'col_109_max_timelag' / &
         ANTags(1109)%Label / 'col_109_flag_threshold' / &
         ANTags(1110)%Label / 'col_109_error_value' / &
         ANTags(1111)%Label / 'col_110_min_value' / &
         ANTags(1112)%Label / 'col_110_max_value' / &
         ANTags(1113)%Label / 'col_110_a_value' / &
         ANTags(1114)%Label / 'col_110_b_value' / &
         ANTags(1115)%Label / 'col_110_nom_timelag' / &
         ANTags(1116)%Label / 'col_110_min_timelag' / &
         ANTags(1117)%Label / 'col_110_max_timelag' / &
         ANTags(1118)%Label / 'col_110_flag_threshold' / &
         ANTags(1119)%Label / 'col_110_error_value' / &
         ANTags(1120)%Label / 'col_111_min_value' / &
         ANTags(1121)%Label / 'col_111_max_value' / &
         ANTags(1122)%Label / 'col_111_a_value' / &
         ANTags(1123)%Label / 'col_111_b_value' / &
         ANTags(1124)%Label / 'col_111_nom_timelag' / &
         ANTags(1125)%Label / 'col_111_min_timelag' / &
         ANTags(1126)%Label / 'col_111_max_timelag' / &
         ANTags(1127)%Label / 'col_111_flag_threshold' / &
         ANTags(1128)%Label / 'col_111_error_value' / &
         ANTags(1129)%Label / 'col_112_min_value' / &
         ANTags(1130)%Label / 'col_112_max_value' / &
         ANTags(1131)%Label / 'col_112_a_value' / &
         ANTags(1132)%Label / 'col_112_b_value' / &
         ANTags(1133)%Label / 'col_112_nom_timelag' / &
         ANTags(1134)%Label / 'col_112_min_timelag' / &
         ANTags(1135)%Label / 'col_112_max_timelag' / &
         ANTags(1136)%Label / 'col_112_flag_threshold' / &
         ANTags(1137)%Label / 'col_112_error_value' / &
         ANTags(1138)%Label / 'col_113_min_value' / &
         ANTags(1139)%Label / 'col_113_max_value' / &
         ANTags(1140)%Label / 'col_113_a_value' / &
         ANTags(1141)%Label / 'col_113_b_value' / &
         ANTags(1142)%Label / 'col_113_nom_timelag' / &
         ANTags(1143)%Label / 'col_113_min_timelag' / &
         ANTags(1144)%Label / 'col_113_max_timelag' / &
         ANTags(1145)%Label / 'col_113_flag_threshold' / &
         ANTags(1146)%Label / 'col_113_error_value' / &
         ANTags(1147)%Label / 'col_114_min_value' / &
         ANTags(1148)%Label / 'col_114_max_value' / &
         ANTags(1149)%Label / 'col_114_a_value' / &
         ANTags(1150)%Label / 'col_114_b_value' / &
         ANTags(1151)%Label / 'col_114_nom_timelag' / &
         ANTags(1152)%Label / 'col_114_min_timelag' / &
         ANTags(1153)%Label / 'col_114_max_timelag' / &
         ANTags(1154)%Label / 'col_114_flag_threshold' / &
         ANTags(1155)%Label / 'col_114_error_value' / &
         ANTags(1156)%Label / 'col_115_min_value' / &
         ANTags(1157)%Label / 'col_115_max_value' / &
         ANTags(1158)%Label / 'col_115_a_value' / &
         ANTags(1159)%Label / 'col_115_b_value' / &
         ANTags(1160)%Label / 'col_115_nom_timelag' / &
         ANTags(1161)%Label / 'col_115_min_timelag' / &
         ANTags(1162)%Label / 'col_115_max_timelag' / &
         ANTags(1163)%Label / 'col_115_flag_threshold' / &
         ANTags(1164)%Label / 'col_115_error_value' / &
         ANTags(1165)%Label / 'col_116_min_value' / &
         ANTags(1166)%Label / 'col_116_max_value' / &
         ANTags(1167)%Label / 'col_116_a_value' / &
         ANTags(1168)%Label / 'col_116_b_value' / &
         ANTags(1169)%Label / 'col_116_nom_timelag' / &
         ANTags(1170)%Label / 'col_116_min_timelag' / &
         ANTags(1171)%Label / 'col_116_max_timelag' / &
         ANTags(1172)%Label / 'col_116_flag_threshold' / &
         ANTags(1173)%Label / 'col_116_error_value' / &
         ANTags(1174)%Label / 'col_117_min_value' / &
         ANTags(1175)%Label / 'col_117_max_value' / &
         ANTags(1176)%Label / 'col_117_a_value' / &
         ANTags(1177)%Label / 'col_117_b_value' / &
         ANTags(1178)%Label / 'col_117_nom_timelag' / &
         ANTags(1179)%Label / 'col_117_min_timelag' / &
         ANTags(1180)%Label / 'col_117_max_timelag' / &
         ANTags(1181)%Label / 'col_117_flag_threshold' / &
         ANTags(1182)%Label / 'col_117_error_value' / &
         ANTags(1183)%Label / 'col_118_min_value' / &
         ANTags(1184)%Label / 'col_118_max_value' / &
         ANTags(1185)%Label / 'col_118_a_value' / &
         ANTags(1186)%Label / 'col_118_b_value' / &
         ANTags(1187)%Label / 'col_118_nom_timelag' / &
         ANTags(1188)%Label / 'col_118_min_timelag' / &
         ANTags(1189)%Label / 'col_118_max_timelag' / &
         ANTags(1190)%Label / 'col_118_flag_threshold' / &
         ANTags(1191)%Label / 'col_118_error_value' / &
         ANTags(1192)%Label / 'col_119_min_value' / &
         ANTags(1193)%Label / 'col_119_max_value' / &
         ANTags(1194)%Label / 'col_119_a_value' / &
         ANTags(1195)%Label / 'col_119_b_value' / &
         ANTags(1196)%Label / 'col_119_nom_timelag' / &
         ANTags(1197)%Label / 'col_119_min_timelag' / &
         ANTags(1198)%Label / 'col_119_max_timelag' / &
         ANTags(1199)%Label / 'col_119_flag_threshold' / &
         ANTags(1200)%Label / 'col_119_error_value' /
    data ANTags(1201)%Label / 'col_120_min_value' / &
         ANTags(1202)%Label / 'col_120_max_value' / &
         ANTags(1203)%Label / 'col_120_a_value' / &
         ANTags(1204)%Label / 'col_120_b_value' / &
         ANTags(1205)%Label / 'col_120_nom_timelag' / &
         ANTags(1206)%Label / 'col_120_min_timelag' / &
         ANTags(1207)%Label / 'col_120_max_timelag' / &
         ANTags(1208)%Label / 'col_120_flag_threshold' / &
         ANTags(1209)%Label / 'col_120_error_value' / &
         ANTags(1210)%Label / 'col_121_min_value' / &
         ANTags(1211)%Label / 'col_121_max_value' / &
         ANTags(1212)%Label / 'col_121_a_value' / &
         ANTags(1213)%Label / 'col_121_b_value' / &
         ANTags(1214)%Label / 'col_121_nom_timelag' / &
         ANTags(1215)%Label / 'col_121_min_timelag' / &
         ANTags(1216)%Label / 'col_121_max_timelag' / &
         ANTags(1217)%Label / 'col_121_flag_threshold' / &
         ANTags(1218)%Label / 'col_121_error_value' / &
         ANTags(1219)%Label / 'col_122_min_value' / &
         ANTags(1220)%Label / 'col_122_max_value' / &
         ANTags(1221)%Label / 'col_122_a_value' / &
         ANTags(1222)%Label / 'col_122_b_value' / &
         ANTags(1223)%Label / 'col_122_nom_timelag' / &
         ANTags(1224)%Label / 'col_122_min_timelag' / &
         ANTags(1225)%Label / 'col_122_max_timelag' / &
         ANTags(1226)%Label / 'col_122_flag_threshold' / &
         ANTags(1227)%Label / 'col_122_error_value' / &
         ANTags(1228)%Label / 'col_123_min_value' / &
         ANTags(1229)%Label / 'col_123_max_value' / &
         ANTags(1230)%Label / 'col_123_a_value' / &
         ANTags(1231)%Label / 'col_123_b_value' / &
         ANTags(1232)%Label / 'col_123_nom_timelag' / &
         ANTags(1233)%Label / 'col_123_min_timelag' / &
         ANTags(1234)%Label / 'col_123_max_timelag' / &
         ANTags(1235)%Label / 'col_123_flag_threshold' / &
         ANTags(1236)%Label / 'col_123_error_value' / &
         ANTags(1237)%Label / 'col_124_min_value' / &
         ANTags(1238)%Label / 'col_124_max_value' / &
         ANTags(1239)%Label / 'col_124_a_value' / &
         ANTags(1240)%Label / 'col_124_b_value' / &
         ANTags(1241)%Label / 'col_124_nom_timelag' / &
         ANTags(1242)%Label / 'col_124_min_timelag' / &
         ANTags(1243)%Label / 'col_124_max_timelag' / &
         ANTags(1244)%Label / 'col_124_flag_threshold' / &
         ANTags(1245)%Label / 'col_124_error_value' / &
         ANTags(1246)%Label / 'col_125_min_value' / &
         ANTags(1247)%Label / 'col_125_max_value' / &
         ANTags(1248)%Label / 'col_125_a_value' / &
         ANTags(1249)%Label / 'col_125_b_value' / &
         ANTags(1250)%Label / 'col_125_nom_timelag' / &
         ANTags(1251)%Label / 'col_125_min_timelag' / &
         ANTags(1252)%Label / 'col_125_max_timelag' / &
         ANTags(1253)%Label / 'col_125_flag_threshold' / &
         ANTags(1254)%Label / 'col_125_error_value' / &
         ANTags(1255)%Label / 'col_126_min_value' / &
         ANTags(1256)%Label / 'col_126_max_value' / &
         ANTags(1257)%Label / 'col_126_a_value' / &
         ANTags(1258)%Label / 'col_126_b_value' / &
         ANTags(1259)%Label / 'col_126_nom_timelag' / &
         ANTags(1260)%Label / 'col_126_min_timelag' / &
         ANTags(1261)%Label / 'col_126_max_timelag' / &
         ANTags(1262)%Label / 'col_126_flag_threshold' / &
         ANTags(1263)%Label / 'col_126_error_value' / &
         ANTags(1264)%Label / 'col_127_min_value' / &
         ANTags(1265)%Label / 'col_127_max_value' / &
         ANTags(1266)%Label / 'col_127_a_value' / &
         ANTags(1267)%Label / 'col_127_b_value' / &
         ANTags(1268)%Label / 'col_127_nom_timelag' / &
         ANTags(1269)%Label / 'col_127_min_timelag' / &
         ANTags(1270)%Label / 'col_127_max_timelag' / &
         ANTags(1271)%Label / 'col_127_flag_threshold' / &
         ANTags(1272)%Label / 'col_127_error_value' / &
         ANTags(1273)%Label / 'col_128_min_value' / &
         ANTags(1274)%Label / 'col_128_max_value' / &
         ANTags(1275)%Label / 'col_128_a_value' / &
         ANTags(1276)%Label / 'col_128_b_value' / &
         ANTags(1277)%Label / 'col_128_nom_timelag' / &
         ANTags(1278)%Label / 'col_128_min_timelag' / &
         ANTags(1279)%Label / 'col_128_max_timelag' / &
         ANTags(1280)%Label / 'col_128_flag_threshold' / &
         ANTags(1281)%Label / 'col_128_error_value' / &
         ANTags(1282)%Label / 'col_129_min_value' / &
         ANTags(1283)%Label / 'col_129_max_value' / &
         ANTags(1284)%Label / 'col_129_a_value' / &
         ANTags(1285)%Label / 'col_129_b_value' / &
         ANTags(1286)%Label / 'col_129_nom_timelag' / &
         ANTags(1287)%Label / 'col_129_min_timelag' / &
         ANTags(1288)%Label / 'col_129_max_timelag' / &
         ANTags(1289)%Label / 'col_129_flag_threshold' / &
         ANTags(1290)%Label / 'col_129_error_value' / &
         ANTags(1291)%Label / 'col_130_min_value' / &
         ANTags(1292)%Label / 'col_130_max_value' / &
         ANTags(1293)%Label / 'col_130_a_value' / &
         ANTags(1294)%Label / 'col_130_b_value' / &
         ANTags(1295)%Label / 'col_130_nom_timelag' / &
         ANTags(1296)%Label / 'col_130_min_timelag' / &
         ANTags(1297)%Label / 'col_130_max_timelag' / &
         ANTags(1298)%Label / 'col_130_flag_threshold' / &
         ANTags(1299)%Label / 'col_130_error_value' / &
         ANTags(1300)%Label / 'col_131_min_value' / &
         ANTags(1301)%Label / 'col_131_max_value' / &
         ANTags(1302)%Label / 'col_131_a_value' / &
         ANTags(1303)%Label / 'col_131_b_value' / &
         ANTags(1304)%Label / 'col_131_nom_timelag' / &
         ANTags(1305)%Label / 'col_131_min_timelag' / &
         ANTags(1306)%Label / 'col_131_max_timelag' / &
         ANTags(1307)%Label / 'col_131_flag_threshold' / &
         ANTags(1308)%Label / 'col_131_error_value' / &
         ANTags(1309)%Label / 'col_132_min_value' / &
         ANTags(1310)%Label / 'col_132_max_value' / &
         ANTags(1311)%Label / 'col_132_a_value' / &
         ANTags(1312)%Label / 'col_132_b_value' / &
         ANTags(1313)%Label / 'col_132_nom_timelag' / &
         ANTags(1314)%Label / 'col_132_min_timelag' / &
         ANTags(1315)%Label / 'col_132_max_timelag' / &
         ANTags(1316)%Label / 'col_132_flag_threshold' / &
         ANTags(1317)%Label / 'col_132_error_value' / &
         ANTags(1318)%Label / 'col_133_min_value' / &
         ANTags(1319)%Label / 'col_133_max_value' / &
         ANTags(1320)%Label / 'col_133_a_value' / &
         ANTags(1321)%Label / 'col_133_b_value' / &
         ANTags(1322)%Label / 'col_133_nom_timelag' / &
         ANTags(1323)%Label / 'col_133_min_timelag' / &
         ANTags(1324)%Label / 'col_133_max_timelag' / &
         ANTags(1325)%Label / 'col_133_flag_threshold' / &
         ANTags(1326)%Label / 'col_133_error_value' / &
         ANTags(1327)%Label / 'col_134_min_value' / &
         ANTags(1328)%Label / 'col_134_max_value' / &
         ANTags(1329)%Label / 'col_134_a_value' / &
         ANTags(1330)%Label / 'col_134_b_value' / &
         ANTags(1331)%Label / 'col_134_nom_timelag' / &
         ANTags(1332)%Label / 'col_134_min_timelag' / &
         ANTags(1333)%Label / 'col_134_max_timelag' / &
         ANTags(1334)%Label / 'col_134_flag_threshold' / &
         ANTags(1335)%Label / 'col_134_error_value' / &
         ANTags(1336)%Label / 'col_135_min_value' / &
         ANTags(1337)%Label / 'col_135_max_value' / &
         ANTags(1338)%Label / 'col_135_a_value' / &
         ANTags(1339)%Label / 'col_135_b_value' / &
         ANTags(1340)%Label / 'col_135_nom_timelag' / &
         ANTags(1341)%Label / 'col_135_min_timelag' / &
         ANTags(1342)%Label / 'col_135_max_timelag' / &
         ANTags(1343)%Label / 'col_135_flag_threshold' / &
         ANTags(1344)%Label / 'col_135_error_value' / &
         ANTags(1345)%Label / 'col_136_min_value' / &
         ANTags(1346)%Label / 'col_136_max_value' / &
         ANTags(1347)%Label / 'col_136_a_value' / &
         ANTags(1348)%Label / 'col_136_b_value' / &
         ANTags(1349)%Label / 'col_136_nom_timelag' / &
         ANTags(1350)%Label / 'col_136_min_timelag' / &
         ANTags(1351)%Label / 'col_136_max_timelag' / &
         ANTags(1352)%Label / 'col_136_flag_threshold' / &
         ANTags(1353)%Label / 'col_136_error_value' / &
         ANTags(1354)%Label / 'col_137_min_value' / &
         ANTags(1355)%Label / 'col_137_max_value' / &
         ANTags(1356)%Label / 'col_137_a_value' / &
         ANTags(1357)%Label / 'col_137_b_value' / &
         ANTags(1358)%Label / 'col_137_nom_timelag' / &
         ANTags(1359)%Label / 'col_137_min_timelag' / &
         ANTags(1360)%Label / 'col_137_max_timelag' / &
         ANTags(1361)%Label / 'col_137_flag_threshold' / &
         ANTags(1362)%Label / 'col_137_error_value' / &
         ANTags(1363)%Label / 'col_138_min_value' / &
         ANTags(1364)%Label / 'col_138_max_value' / &
         ANTags(1365)%Label / 'col_138_a_value' / &
         ANTags(1366)%Label / 'col_138_b_value' / &
         ANTags(1367)%Label / 'col_138_nom_timelag' / &
         ANTags(1368)%Label / 'col_138_min_timelag' / &
         ANTags(1369)%Label / 'col_138_max_timelag' / &
         ANTags(1370)%Label / 'col_138_flag_threshold' / &
         ANTags(1371)%Label / 'col_138_error_value' / &
         ANTags(1372)%Label / 'col_139_min_value' / &
         ANTags(1373)%Label / 'col_139_max_value' / &
         ANTags(1374)%Label / 'col_139_a_value' / &
         ANTags(1375)%Label / 'col_139_b_value' / &
         ANTags(1376)%Label / 'col_139_nom_timelag' / &
         ANTags(1377)%Label / 'col_139_min_timelag' / &
         ANTags(1378)%Label / 'col_139_max_timelag' / &
         ANTags(1379)%Label / 'col_139_flag_threshold' / &
         ANTags(1380)%Label / 'col_139_error_value' / &
         ANTags(1381)%Label / 'col_140_min_value' / &
         ANTags(1382)%Label / 'col_140_max_value' / &
         ANTags(1383)%Label / 'col_140_a_value' / &
         ANTags(1384)%Label / 'col_140_b_value' / &
         ANTags(1385)%Label / 'col_140_nom_timelag' / &
         ANTags(1386)%Label / 'col_140_min_timelag' / &
         ANTags(1387)%Label / 'col_140_max_timelag' / &
         ANTags(1388)%Label / 'col_140_flag_threshold' / &
         ANTags(1389)%Label / 'col_140_error_value' / &
         ANTags(1390)%Label / 'col_141_min_value' / &
         ANTags(1391)%Label / 'col_141_max_value' / &
         ANTags(1392)%Label / 'col_141_a_value' / &
         ANTags(1393)%Label / 'col_141_b_value' / &
         ANTags(1394)%Label / 'col_141_nom_timelag' / &
         ANTags(1395)%Label / 'col_141_min_timelag' / &
         ANTags(1396)%Label / 'col_141_max_timelag' / &
         ANTags(1397)%Label / 'col_141_flag_threshold' / &
         ANTags(1398)%Label / 'col_141_error_value' / &
         ANTags(1399)%Label / 'col_142_min_value' / &
         ANTags(1400)%Label / 'col_142_max_value' /
    data ANTags(1401)%Label / 'col_142_a_value' / &
         ANTags(1402)%Label / 'col_142_b_value' / &
         ANTags(1403)%Label / 'col_142_nom_timelag' / &
         ANTags(1404)%Label / 'col_142_min_timelag' / &
         ANTags(1405)%Label / 'col_142_max_timelag' / &
         ANTags(1406)%Label / 'col_142_flag_threshold' / &
         ANTags(1407)%Label / 'col_142_error_value' / &
         ANTags(1408)%Label / 'col_143_min_value' / &
         ANTags(1409)%Label / 'col_143_max_value' / &
         ANTags(1410)%Label / 'col_143_a_value' / &
         ANTags(1411)%Label / 'col_143_b_value' / &
         ANTags(1412)%Label / 'col_143_nom_timelag' / &
         ANTags(1413)%Label / 'col_143_min_timelag' / &
         ANTags(1414)%Label / 'col_143_max_timelag' / &
         ANTags(1415)%Label / 'col_143_flag_threshold' / &
         ANTags(1416)%Label / 'col_143_error_value' / &
         ANTags(1417)%Label / 'col_144_min_value' / &
         ANTags(1418)%Label / 'col_144_max_value' / &
         ANTags(1419)%Label / 'col_144_a_value' / &
         ANTags(1420)%Label / 'col_144_b_value' / &
         ANTags(1421)%Label / 'col_144_nom_timelag' / &
         ANTags(1422)%Label / 'col_144_min_timelag' / &
         ANTags(1423)%Label / 'col_144_max_timelag' / &
         ANTags(1424)%Label / 'col_144_flag_threshold' / &
         ANTags(1425)%Label / 'col_144_error_value' / &
         ANTags(1426)%Label / 'col_145_min_value' / &
         ANTags(1427)%Label / 'col_145_max_value' / &
         ANTags(1428)%Label / 'col_145_a_value' / &
         ANTags(1429)%Label / 'col_145_b_value' / &
         ANTags(1430)%Label / 'col_145_nom_timelag' / &
         ANTags(1431)%Label / 'col_145_min_timelag' / &
         ANTags(1432)%Label / 'col_145_max_timelag' / &
         ANTags(1433)%Label / 'col_145_flag_threshold' / &
         ANTags(1434)%Label / 'col_145_error_value' / &
         ANTags(1435)%Label / 'col_146_min_value' / &
         ANTags(1436)%Label / 'col_146_max_value' / &
         ANTags(1437)%Label / 'col_146_a_value' / &
         ANTags(1438)%Label / 'col_146_b_value' / &
         ANTags(1439)%Label / 'col_146_nom_timelag' / &
         ANTags(1440)%Label / 'col_146_min_timelag' / &
         ANTags(1441)%Label / 'col_146_max_timelag' / &
         ANTags(1442)%Label / 'col_146_flag_threshold' / &
         ANTags(1443)%Label / 'col_146_error_value' / &
         ANTags(1444)%Label / 'col_147_min_value' / &
         ANTags(1445)%Label / 'col_147_max_value' / &
         ANTags(1446)%Label / 'col_147_a_value' / &
         ANTags(1447)%Label / 'col_147_b_value' / &
         ANTags(1448)%Label / 'col_147_nom_timelag' / &
         ANTags(1449)%Label / 'col_147_min_timelag' / &
         ANTags(1450)%Label / 'col_147_max_timelag' / &
         ANTags(1451)%Label / 'col_147_flag_threshold' / &
         ANTags(1452)%Label / 'col_147_error_value' / &
         ANTags(1453)%Label / 'col_148_min_value' / &
         ANTags(1454)%Label / 'col_148_max_value' / &
         ANTags(1455)%Label / 'col_148_a_value' / &
         ANTags(1456)%Label / 'col_148_b_value' / &
         ANTags(1457)%Label / 'col_148_nom_timelag' / &
         ANTags(1458)%Label / 'col_148_min_timelag' / &
         ANTags(1459)%Label / 'col_148_max_timelag' / &
         ANTags(1460)%Label / 'col_148_flag_threshold' / &
         ANTags(1461)%Label / 'col_148_error_value' / &
         ANTags(1462)%Label / 'col_149_min_value' / &
         ANTags(1463)%Label / 'col_149_max_value' / &
         ANTags(1464)%Label / 'col_149_a_value' / &
         ANTags(1465)%Label / 'col_149_b_value' / &
         ANTags(1466)%Label / 'col_149_nom_timelag' / &
         ANTags(1467)%Label / 'col_149_min_timelag' / &
         ANTags(1468)%Label / 'col_149_max_timelag' / &
         ANTags(1469)%Label / 'col_149_flag_threshold' / &
         ANTags(1470)%Label / 'col_149_error_value' / &
         ANTags(1471)%Label / 'col_150_min_value' / &
         ANTags(1472)%Label / 'col_150_max_value' / &
         ANTags(1473)%Label / 'col_150_a_value' / &
         ANTags(1474)%Label / 'col_150_b_value' / &
         ANTags(1475)%Label / 'col_150_nom_timelag' / &
         ANTags(1476)%Label / 'col_150_min_timelag' / &
         ANTags(1477)%Label / 'col_150_max_timelag' / &
         ANTags(1478)%Label / 'col_150_flag_threshold' / &
         ANTags(1479)%Label / 'col_150_error_value' / &
         ANTags(1480)%Label / 'col_151_min_value' / &
         ANTags(1481)%Label / 'col_151_max_value' / &
         ANTags(1482)%Label / 'col_151_a_value' / &
         ANTags(1483)%Label / 'col_151_b_value' / &
         ANTags(1484)%Label / 'col_151_nom_timelag' / &
         ANTags(1485)%Label / 'col_151_min_timelag' / &
         ANTags(1486)%Label / 'col_151_max_timelag' / &
         ANTags(1487)%Label / 'col_151_flag_threshold' / &
         ANTags(1488)%Label / 'col_151_error_value' / &
         ANTags(1489)%Label / 'col_152_min_value' / &
         ANTags(1490)%Label / 'col_152_max_value' / &
         ANTags(1491)%Label / 'col_152_a_value' / &
         ANTags(1492)%Label / 'col_152_b_value' / &
         ANTags(1493)%Label / 'col_152_nom_timelag' / &
         ANTags(1494)%Label / 'col_152_min_timelag' / &
         ANTags(1495)%Label / 'col_152_max_timelag' / &
         ANTags(1496)%Label / 'col_152_flag_threshold' / &
         ANTags(1497)%Label / 'col_152_error_value' / &
         ANTags(1498)%Label / 'col_153_min_value' / &
         ANTags(1499)%Label / 'col_153_max_value' / &
         ANTags(1500)%Label / 'col_153_a_value' / &
         ANTags(1501)%Label / 'col_153_b_value' / &
         ANTags(1502)%Label / 'col_153_nom_timelag' / &
         ANTags(1503)%Label / 'col_153_min_timelag' / &
         ANTags(1504)%Label / 'col_153_max_timelag' / &
         ANTags(1505)%Label / 'col_153_flag_threshold' / &
         ANTags(1506)%Label / 'col_153_error_value' / &
         ANTags(1507)%Label / 'col_154_min_value' / &
         ANTags(1508)%Label / 'col_154_max_value' / &
         ANTags(1509)%Label / 'col_154_a_value' / &
         ANTags(1510)%Label / 'col_154_b_value' / &
         ANTags(1511)%Label / 'col_154_nom_timelag' / &
         ANTags(1512)%Label / 'col_154_min_timelag' / &
         ANTags(1513)%Label / 'col_154_max_timelag' / &
         ANTags(1514)%Label / 'col_154_flag_threshold' / &
         ANTags(1515)%Label / 'col_154_error_value' / &
         ANTags(1516)%Label / 'col_155_min_value' / &
         ANTags(1517)%Label / 'col_155_max_value' / &
         ANTags(1518)%Label / 'col_155_a_value' / &
         ANTags(1519)%Label / 'col_155_b_value' / &
         ANTags(1520)%Label / 'col_155_nom_timelag' / &
         ANTags(1521)%Label / 'col_155_min_timelag' / &
         ANTags(1522)%Label / 'col_155_max_timelag' / &
         ANTags(1523)%Label / 'col_155_flag_threshold' / &
         ANTags(1524)%Label / 'col_155_error_value' / &
         ANTags(1525)%Label / 'col_156_min_value' / &
         ANTags(1526)%Label / 'col_156_max_value' / &
         ANTags(1527)%Label / 'col_156_a_value' / &
         ANTags(1528)%Label / 'col_156_b_value' / &
         ANTags(1529)%Label / 'col_156_nom_timelag' / &
         ANTags(1530)%Label / 'col_156_min_timelag' / &
         ANTags(1531)%Label / 'col_156_max_timelag' / &
         ANTags(1532)%Label / 'col_156_flag_threshold' / &
         ANTags(1533)%Label / 'col_156_error_value' / &
         ANTags(1534)%Label / 'col_157_min_value' / &
         ANTags(1535)%Label / 'col_157_max_value' / &
         ANTags(1536)%Label / 'col_157_a_value' / &
         ANTags(1537)%Label / 'col_157_b_value' / &
         ANTags(1538)%Label / 'col_157_nom_timelag' / &
         ANTags(1539)%Label / 'col_157_min_timelag' / &
         ANTags(1540)%Label / 'col_157_max_timelag' / &
         ANTags(1541)%Label / 'col_157_flag_threshold' / &
         ANTags(1542)%Label / 'col_157_error_value' / &
         ANTags(1543)%Label / 'col_158_min_value' / &
         ANTags(1544)%Label / 'col_158_max_value' / &
         ANTags(1545)%Label / 'col_158_a_value' / &
         ANTags(1546)%Label / 'col_158_b_value' / &
         ANTags(1547)%Label / 'col_158_nom_timelag' / &
         ANTags(1548)%Label / 'col_158_min_timelag' / &
         ANTags(1549)%Label / 'col_158_max_timelag' / &
         ANTags(1550)%Label / 'col_158_flag_threshold' / &
         ANTags(1551)%Label / 'col_158_error_value' / &
         ANTags(1552)%Label / 'col_159_min_value' / &
         ANTags(1553)%Label / 'col_159_max_value' / &
         ANTags(1554)%Label / 'col_159_a_value' / &
         ANTags(1555)%Label / 'col_159_b_value' / &
         ANTags(1556)%Label / 'col_159_nom_timelag' / &
         ANTags(1557)%Label / 'col_159_min_timelag' / &
         ANTags(1558)%Label / 'col_159_max_timelag' / &
         ANTags(1559)%Label / 'col_159_flag_threshold' / &
         ANTags(1560)%Label / 'col_159_error_value' / &
         ANTags(1561)%Label / 'col_160_min_value' / &
         ANTags(1562)%Label / 'col_160_max_value' / &
         ANTags(1563)%Label / 'col_160_a_value' / &
         ANTags(1564)%Label / 'col_160_b_value' / &
         ANTags(1565)%Label / 'col_160_nom_timelag' / &
         ANTags(1566)%Label / 'col_160_min_timelag' / &
         ANTags(1567)%Label / 'col_160_max_timelag' / &
         ANTags(1568)%Label / 'col_160_flag_threshold' / &
         ANTags(1569)%Label / 'col_160_error_value' / &
         ANTags(1570)%Label / 'col_161_min_value' / &
         ANTags(1571)%Label / 'col_161_max_value' / &
         ANTags(1572)%Label / 'col_161_a_value' / &
         ANTags(1573)%Label / 'col_161_b_value' / &
         ANTags(1574)%Label / 'col_161_nom_timelag' / &
         ANTags(1575)%Label / 'col_161_min_timelag' / &
         ANTags(1576)%Label / 'col_161_max_timelag' / &
         ANTags(1577)%Label / 'col_161_flag_threshold' / &
         ANTags(1578)%Label / 'col_161_error_value' / &
         ANTags(1579)%Label / 'col_162_min_value' / &
         ANTags(1580)%Label / 'col_162_max_value' / &
         ANTags(1581)%Label / 'col_162_a_value' / &
         ANTags(1582)%Label / 'col_162_b_value' / &
         ANTags(1583)%Label / 'col_162_nom_timelag' / &
         ANTags(1584)%Label / 'col_162_min_timelag' / &
         ANTags(1585)%Label / 'col_162_max_timelag' / &
         ANTags(1586)%Label / 'col_162_flag_threshold' / &
         ANTags(1587)%Label / 'col_162_error_value' / &
         ANTags(1588)%Label / 'col_163_min_value' / &
         ANTags(1589)%Label / 'col_163_max_value' / &
         ANTags(1590)%Label / 'col_163_a_value' / &
         ANTags(1591)%Label / 'col_163_b_value' / &
         ANTags(1592)%Label / 'col_163_nom_timelag' / &
         ANTags(1593)%Label / 'col_163_min_timelag' / &
         ANTags(1594)%Label / 'col_163_max_timelag' / &
         ANTags(1595)%Label / 'col_163_flag_threshold' / &
         ANTags(1596)%Label / 'col_163_error_value' / &
         ANTags(1597)%Label / 'col_164_min_value' / &
         ANTags(1598)%Label / 'col_164_max_value' / &
         ANTags(1599)%Label / 'col_164_a_value' / &
         ANTags(1600)%Label / 'col_164_b_value' /
    data ANTags(1601)%Label / 'col_164_nom_timelag' / &
         ANTags(1602)%Label / 'col_164_min_timelag' / &
         ANTags(1603)%Label / 'col_164_max_timelag' / &
         ANTags(1604)%Label / 'col_164_flag_threshold' / &
         ANTags(1605)%Label / 'col_164_error_value' / &
         ANTags(1606)%Label / 'col_165_min_value' / &
         ANTags(1607)%Label / 'col_165_max_value' / &
         ANTags(1608)%Label / 'col_165_a_value' / &
         ANTags(1609)%Label / 'col_165_b_value' / &
         ANTags(1610)%Label / 'col_165_nom_timelag' / &
         ANTags(1611)%Label / 'col_165_min_timelag' / &
         ANTags(1612)%Label / 'col_165_max_timelag' / &
         ANTags(1613)%Label / 'col_165_flag_threshold' / &
         ANTags(1614)%Label / 'col_165_error_value' / &
         ANTags(1615)%Label / 'col_166_min_value' / &
         ANTags(1616)%Label / 'col_166_max_value' / &
         ANTags(1617)%Label / 'col_166_a_value' / &
         ANTags(1618)%Label / 'col_166_b_value' / &
         ANTags(1619)%Label / 'col_166_nom_timelag' / &
         ANTags(1620)%Label / 'col_166_min_timelag' / &
         ANTags(1621)%Label / 'col_166_max_timelag' / &
         ANTags(1622)%Label / 'col_166_flag_threshold' / &
         ANTags(1623)%Label / 'col_166_error_value' / &
         ANTags(1624)%Label / 'col_167_min_value' / &
         ANTags(1625)%Label / 'col_167_max_value' / &
         ANTags(1626)%Label / 'col_167_a_value' / &
         ANTags(1627)%Label / 'col_167_b_value' / &
         ANTags(1628)%Label / 'col_167_nom_timelag' / &
         ANTags(1629)%Label / 'col_167_min_timelag' / &
         ANTags(1630)%Label / 'col_167_max_timelag' / &
         ANTags(1631)%Label / 'col_167_flag_threshold' / &
         ANTags(1632)%Label / 'col_167_error_value' / &
         ANTags(1633)%Label / 'col_168_min_value' / &
         ANTags(1634)%Label / 'col_168_max_value' / &
         ANTags(1635)%Label / 'col_168_a_value' / &
         ANTags(1636)%Label / 'col_168_b_value' / &
         ANTags(1637)%Label / 'col_168_nom_timelag' / &
         ANTags(1638)%Label / 'col_168_min_timelag' / &
         ANTags(1639)%Label / 'col_168_max_timelag' / &
         ANTags(1640)%Label / 'col_168_flag_threshold' / &
         ANTags(1641)%Label / 'col_168_error_value' / &
         ANTags(1642)%Label / 'col_169_min_value' / &
         ANTags(1643)%Label / 'col_169_max_value' / &
         ANTags(1644)%Label / 'col_169_a_value' / &
         ANTags(1645)%Label / 'col_169_b_value' / &
         ANTags(1646)%Label / 'col_169_nom_timelag' / &
         ANTags(1647)%Label / 'col_169_min_timelag' / &
         ANTags(1648)%Label / 'col_169_max_timelag' / &
         ANTags(1649)%Label / 'col_169_flag_threshold' / &
         ANTags(1650)%Label / 'col_169_error_value' / &
         ANTags(1651)%Label / 'col_170_min_value' / &
         ANTags(1652)%Label / 'col_170_max_value' / &
         ANTags(1653)%Label / 'col_170_a_value' / &
         ANTags(1654)%Label / 'col_170_b_value' / &
         ANTags(1655)%Label / 'col_170_nom_timelag' / &
         ANTags(1656)%Label / 'col_170_min_timelag' / &
         ANTags(1657)%Label / 'col_170_max_timelag' / &
         ANTags(1658)%Label / 'col_170_flag_threshold' / &
         ANTags(1659)%Label / 'col_170_error_value' / &
         ANTags(1660)%Label / 'col_171_min_value' / &
         ANTags(1661)%Label / 'col_171_max_value' / &
         ANTags(1662)%Label / 'col_171_a_value' / &
         ANTags(1663)%Label / 'col_171_b_value' / &
         ANTags(1664)%Label / 'col_171_nom_timelag' / &
         ANTags(1665)%Label / 'col_171_min_timelag' / &
         ANTags(1666)%Label / 'col_171_max_timelag' / &
         ANTags(1667)%Label / 'col_171_flag_threshold' / &
         ANTags(1668)%Label / 'col_171_error_value' / &
         ANTags(1669)%Label / 'col_172_min_value' / &
         ANTags(1670)%Label / 'col_172_max_value' / &
         ANTags(1671)%Label / 'col_172_a_value' / &
         ANTags(1672)%Label / 'col_172_b_value' / &
         ANTags(1673)%Label / 'col_172_nom_timelag' / &
         ANTags(1674)%Label / 'col_172_min_timelag' / &
         ANTags(1675)%Label / 'col_172_max_timelag' / &
         ANTags(1676)%Label / 'col_172_flag_threshold' / &
         ANTags(1677)%Label / 'col_172_error_value' / &
         ANTags(1678)%Label / 'col_173_min_value' / &
         ANTags(1679)%Label / 'col_173_max_value' / &
         ANTags(1680)%Label / 'col_173_a_value' / &
         ANTags(1681)%Label / 'col_173_b_value' / &
         ANTags(1682)%Label / 'col_173_nom_timelag' / &
         ANTags(1683)%Label / 'col_173_min_timelag' / &
         ANTags(1684)%Label / 'col_173_max_timelag' / &
         ANTags(1685)%Label / 'col_173_flag_threshold' / &
         ANTags(1686)%Label / 'col_173_error_value' / &
         ANTags(1687)%Label / 'col_174_min_value' / &
         ANTags(1688)%Label / 'col_174_max_value' / &
         ANTags(1689)%Label / 'col_174_a_value' / &
         ANTags(1690)%Label / 'col_174_b_value' / &
         ANTags(1691)%Label / 'col_174_nom_timelag' / &
         ANTags(1692)%Label / 'col_174_min_timelag' / &
         ANTags(1693)%Label / 'col_174_max_timelag' / &
         ANTags(1694)%Label / 'col_174_flag_threshold' / &
         ANTags(1695)%Label / 'col_174_error_value' / &
         ANTags(1696)%Label / 'col_175_min_value' / &
         ANTags(1697)%Label / 'col_175_max_value' / &
         ANTags(1698)%Label / 'col_175_a_value' / &
         ANTags(1699)%Label / 'col_175_b_value' / &
         ANTags(1700)%Label / 'col_175_nom_timelag' / &
         ANTags(1701)%Label / 'col_175_min_timelag' / &
         ANTags(1702)%Label / 'col_175_max_timelag' / &
         ANTags(1703)%Label / 'col_175_flag_threshold' / &
         ANTags(1704)%Label / 'col_175_error_value' / &
         ANTags(1705)%Label / 'col_176_min_value' / &
         ANTags(1706)%Label / 'col_176_max_value' / &
         ANTags(1707)%Label / 'col_176_a_value' / &
         ANTags(1708)%Label / 'col_176_b_value' / &
         ANTags(1709)%Label / 'col_176_nom_timelag' / &
         ANTags(1710)%Label / 'col_176_min_timelag' / &
         ANTags(1711)%Label / 'col_176_max_timelag' / &
         ANTags(1712)%Label / 'col_176_flag_threshold' / &
         ANTags(1713)%Label / 'col_176_error_value' / &
         ANTags(1714)%Label / 'col_177_min_value' / &
         ANTags(1715)%Label / 'col_177_max_value' / &
         ANTags(1716)%Label / 'col_177_a_value' / &
         ANTags(1717)%Label / 'col_177_b_value' / &
         ANTags(1718)%Label / 'col_177_nom_timelag' / &
         ANTags(1719)%Label / 'col_177_min_timelag' / &
         ANTags(1720)%Label / 'col_177_max_timelag' / &
         ANTags(1721)%Label / 'col_177_flag_threshold' / &
         ANTags(1722)%Label / 'col_177_error_value' / &
         ANTags(1723)%Label / 'col_178_min_value' / &
         ANTags(1724)%Label / 'col_178_max_value' / &
         ANTags(1725)%Label / 'col_178_a_value' / &
         ANTags(1726)%Label / 'col_178_b_value' / &
         ANTags(1727)%Label / 'col_178_nom_timelag' / &
         ANTags(1728)%Label / 'col_178_min_timelag' / &
         ANTags(1729)%Label / 'col_178_max_timelag' / &
         ANTags(1730)%Label / 'col_178_flag_threshold' / &
         ANTags(1731)%Label / 'col_178_error_value' / &
         ANTags(1732)%Label / 'col_179_min_value' / &
         ANTags(1733)%Label / 'col_179_max_value' / &
         ANTags(1734)%Label / 'col_179_a_value' / &
         ANTags(1735)%Label / 'col_179_b_value' / &
         ANTags(1736)%Label / 'col_179_nom_timelag' / &
         ANTags(1737)%Label / 'col_179_min_timelag' / &
         ANTags(1738)%Label / 'col_179_max_timelag' / &
         ANTags(1739)%Label / 'col_179_flag_threshold' / &
         ANTags(1740)%Label / 'col_179_error_value' / &
         ANTags(1741)%Label / 'col_180_min_value' / &
         ANTags(1742)%Label / 'col_180_max_value' / &
         ANTags(1743)%Label / 'col_180_a_value' / &
         ANTags(1744)%Label / 'col_180_b_value' / &
         ANTags(1745)%Label / 'col_180_nom_timelag' / &
         ANTags(1746)%Label / 'col_180_min_timelag' / &
         ANTags(1747)%Label / 'col_180_max_timelag' / &
         ANTags(1748)%Label / 'col_180_flag_threshold' / &
         ANTags(1749)%Label / 'col_180_error_value' / &
         ANTags(1750)%Label / 'col_181_min_value' / &
         ANTags(1751)%Label / 'col_181_max_value' / &
         ANTags(1752)%Label / 'col_181_a_value' / &
         ANTags(1753)%Label / 'col_181_b_value' / &
         ANTags(1754)%Label / 'col_181_nom_timelag' / &
         ANTags(1755)%Label / 'col_181_min_timelag' / &
         ANTags(1756)%Label / 'col_181_max_timelag' / &
         ANTags(1757)%Label / 'col_181_flag_threshold' / &
         ANTags(1758)%Label / 'col_181_error_value' / &
         ANTags(1759)%Label / 'col_182_min_value' / &
         ANTags(1760)%Label / 'col_182_max_value' / &
         ANTags(1761)%Label / 'col_182_a_value' / &
         ANTags(1762)%Label / 'col_182_b_value' / &
         ANTags(1763)%Label / 'col_182_nom_timelag' / &
         ANTags(1764)%Label / 'col_182_min_timelag' / &
         ANTags(1765)%Label / 'col_182_max_timelag' / &
         ANTags(1766)%Label / 'col_182_flag_threshold' / &
         ANTags(1767)%Label / 'col_182_error_value' / &
         ANTags(1768)%Label / 'col_183_min_value' / &
         ANTags(1769)%Label / 'col_183_max_value' / &
         ANTags(1770)%Label / 'col_183_a_value' / &
         ANTags(1771)%Label / 'col_183_b_value' / &
         ANTags(1772)%Label / 'col_183_nom_timelag' / &
         ANTags(1773)%Label / 'col_183_min_timelag' / &
         ANTags(1774)%Label / 'col_183_max_timelag' / &
         ANTags(1775)%Label / 'col_183_flag_threshold' / &
         ANTags(1776)%Label / 'col_183_error_value' / &
         ANTags(1777)%Label / 'col_184_min_value' / &
         ANTags(1778)%Label / 'col_184_max_value' / &
         ANTags(1779)%Label / 'col_184_a_value' / &
         ANTags(1780)%Label / 'col_184_b_value' / &
         ANTags(1781)%Label / 'col_184_nom_timelag' / &
         ANTags(1782)%Label / 'col_184_min_timelag' / &
         ANTags(1783)%Label / 'col_184_max_timelag' / &
         ANTags(1784)%Label / 'col_184_flag_threshold' / &
         ANTags(1785)%Label / 'col_184_error_value' / &
         ANTags(1786)%Label / 'col_185_min_value' / &
         ANTags(1787)%Label / 'col_185_max_value' / &
         ANTags(1788)%Label / 'col_185_a_value' / &
         ANTags(1789)%Label / 'col_185_b_value' / &
         ANTags(1790)%Label / 'col_185_nom_timelag' / &
         ANTags(1791)%Label / 'col_185_min_timelag' / &
         ANTags(1792)%Label / 'col_185_max_timelag' / &
         ANTags(1793)%Label / 'col_185_flag_threshold' / &
         ANTags(1794)%Label / 'col_185_error_value' / &
         ANTags(1795)%Label / 'col_186_min_value' / &
         ANTags(1796)%Label / 'col_186_max_value' / &
         ANTags(1797)%Label / 'col_186_a_value' / &
         ANTags(1798)%Label / 'col_186_b_value' / &
         ANTags(1799)%Label / 'col_186_nom_timelag' / &
         ANTags(1800)%Label / 'col_186_min_timelag' /
    data ANTags(1801)%Label / 'col_186_max_timelag' / &
         ANTags(1802)%Label / 'col_186_flag_threshold' / &
         ANTags(1803)%Label / 'col_186_error_value' / &
         ANTags(1804)%Label / 'col_187_min_value' / &
         ANTags(1805)%Label / 'col_187_max_value' / &
         ANTags(1806)%Label / 'col_187_a_value' / &
         ANTags(1807)%Label / 'col_187_b_value' / &
         ANTags(1808)%Label / 'col_187_nom_timelag' / &
         ANTags(1809)%Label / 'col_187_min_timelag' / &
         ANTags(1810)%Label / 'col_187_max_timelag' / &
         ANTags(1811)%Label / 'col_187_flag_threshold' / &
         ANTags(1812)%Label / 'col_187_error_value' / &
         ANTags(1813)%Label / 'col_188_min_value' / &
         ANTags(1814)%Label / 'col_188_max_value' / &
         ANTags(1815)%Label / 'col_188_a_value' / &
         ANTags(1816)%Label / 'col_188_b_value' / &
         ANTags(1817)%Label / 'col_188_nom_timelag' / &
         ANTags(1818)%Label / 'col_188_min_timelag' / &
         ANTags(1819)%Label / 'col_188_max_timelag' / &
         ANTags(1820)%Label / 'col_188_flag_threshold' / &
         ANTags(1821)%Label / 'col_188_error_value' / &
         ANTags(1822)%Label / 'col_189_min_value' / &
         ANTags(1823)%Label / 'col_189_max_value' / &
         ANTags(1824)%Label / 'col_189_a_value' / &
         ANTags(1825)%Label / 'col_189_b_value' / &
         ANTags(1826)%Label / 'col_189_nom_timelag' / &
         ANTags(1827)%Label / 'col_189_min_timelag' / &
         ANTags(1828)%Label / 'col_189_max_timelag' / &
         ANTags(1829)%Label / 'col_189_flag_threshold' / &
         ANTags(1830)%Label / 'col_189_error_value' / &
         ANTags(1831)%Label / 'col_190_min_value' / &
         ANTags(1832)%Label / 'col_190_max_value' / &
         ANTags(1833)%Label / 'col_190_a_value' / &
         ANTags(1834)%Label / 'col_190_b_value' / &
         ANTags(1835)%Label / 'col_190_nom_timelag' / &
         ANTags(1836)%Label / 'col_190_min_timelag' / &
         ANTags(1837)%Label / 'col_190_max_timelag' / &
         ANTags(1838)%Label / 'col_190_flag_threshold' / &
         ANTags(1839)%Label / 'col_190_error_value' / &
         ANTags(1840)%Label / 'col_191_min_value' / &
         ANTags(1841)%Label / 'col_191_max_value' / &
         ANTags(1842)%Label / 'col_191_a_value' / &
         ANTags(1843)%Label / 'col_191_b_value' / &
         ANTags(1844)%Label / 'col_191_nom_timelag' / &
         ANTags(1845)%Label / 'col_191_min_timelag' / &
         ANTags(1846)%Label / 'col_191_max_timelag' / &
         ANTags(1847)%Label / 'col_191_flag_threshold' / &
         ANTags(1848)%Label / 'col_191_error_value' / &
         ANTags(1849)%Label / 'col_192_min_value' / &
         ANTags(1850)%Label / 'col_192_max_value' / &
         ANTags(1851)%Label / 'col_192_a_value' / &
         ANTags(1852)%Label / 'col_192_b_value' / &
         ANTags(1853)%Label / 'col_192_nom_timelag' / &
         ANTags(1854)%Label / 'col_192_min_timelag' / &
         ANTags(1855)%Label / 'col_192_max_timelag' / &
         ANTags(1856)%Label / 'col_192_flag_threshold' / &
         ANTags(1857)%Label / 'col_192_error_value' / &
         ANTags(1858)%Label / 'col_193_min_value' / &
         ANTags(1859)%Label / 'col_193_max_value' / &
         ANTags(1860)%Label / 'col_193_a_value' / &
         ANTags(1861)%Label / 'col_193_b_value' / &
         ANTags(1862)%Label / 'col_193_nom_timelag' / &
         ANTags(1863)%Label / 'col_193_min_timelag' / &
         ANTags(1864)%Label / 'col_193_max_timelag' / &
         ANTags(1865)%Label / 'col_193_flag_threshold' / &
         ANTags(1866)%Label / 'col_193_error_value' / &
         ANTags(1867)%Label / 'col_194_min_value' / &
         ANTags(1868)%Label / 'col_194_max_value' / &
         ANTags(1869)%Label / 'col_194_a_value' / &
         ANTags(1870)%Label / 'col_194_b_value' / &
         ANTags(1871)%Label / 'col_194_nom_timelag' / &
         ANTags(1872)%Label / 'col_194_min_timelag' / &
         ANTags(1873)%Label / 'col_194_max_timelag' / &
         ANTags(1874)%Label / 'col_194_flag_threshold' / &
         ANTags(1875)%Label / 'col_194_error_value' / &
         ANTags(1876)%Label / 'col_195_min_value' / &
         ANTags(1877)%Label / 'col_195_max_value' / &
         ANTags(1878)%Label / 'col_195_a_value' / &
         ANTags(1879)%Label / 'col_195_b_value' / &
         ANTags(1880)%Label / 'col_195_nom_timelag' / &
         ANTags(1881)%Label / 'col_195_min_timelag' / &
         ANTags(1882)%Label / 'col_195_max_timelag' / &
         ANTags(1883)%Label / 'col_195_flag_threshold' / &
         ANTags(1884)%Label / 'col_195_error_value' / &
         ANTags(1885)%Label / 'col_196_min_value' / &
         ANTags(1886)%Label / 'col_196_max_value' / &
         ANTags(1887)%Label / 'col_196_a_value' / &
         ANTags(1888)%Label / 'col_196_b_value' / &
         ANTags(1889)%Label / 'col_196_nom_timelag' / &
         ANTags(1890)%Label / 'col_196_min_timelag' / &
         ANTags(1891)%Label / 'col_196_max_timelag' / &
         ANTags(1892)%Label / 'col_196_flag_threshold' / &
         ANTags(1893)%Label / 'col_196_error_value' / &
         ANTags(1894)%Label / 'col_197_min_value' / &
         ANTags(1895)%Label / 'col_197_max_value' / &
         ANTags(1896)%Label / 'col_197_a_value' / &
         ANTags(1897)%Label / 'col_197_b_value' / &
         ANTags(1898)%Label / 'col_197_nom_timelag' / &
         ANTags(1899)%Label / 'col_197_min_timelag' / &
         ANTags(1900)%Label / 'col_197_max_timelag' / &
         ANTags(1901)%Label / 'col_197_flag_threshold' / &
         ANTags(1902)%Label / 'col_197_error_value' / &
         ANTags(1903)%Label / 'col_198_min_value' / &
         ANTags(1904)%Label / 'col_198_max_value' / &
         ANTags(1905)%Label / 'col_198_a_value' / &
         ANTags(1906)%Label / 'col_198_b_value' / &
         ANTags(1907)%Label / 'col_198_nom_timelag' / &
         ANTags(1908)%Label / 'col_198_min_timelag' / &
         ANTags(1909)%Label / 'col_198_max_timelag' / &
         ANTags(1910)%Label / 'col_198_flag_threshold' / &
         ANTags(1911)%Label / 'col_198_error_value' / &
         ANTags(1912)%Label / 'col_199_min_value' / &
         ANTags(1913)%Label / 'col_199_max_value' / &
         ANTags(1914)%Label / 'col_199_a_value' / &
         ANTags(1915)%Label / 'col_199_b_value' / &
         ANTags(1916)%Label / 'col_199_nom_timelag' / &
         ANTags(1917)%Label / 'col_199_min_timelag' / &
         ANTags(1918)%Label / 'col_199_max_timelag' / &
         ANTags(1919)%Label / 'col_199_flag_threshold' / &
         ANTags(1920)%Label / 'col_199_error_value' / &
         ANTags(1921)%Label / 'col_200_min_value' / &
         ANTags(1922)%Label / 'col_200_max_value' / &
         ANTags(1923)%Label / 'col_200_a_value' / &
         ANTags(1924)%Label / 'col_200_b_value' / &
         ANTags(1925)%Label / 'col_200_nom_timelag' / &
         ANTags(1926)%Label / 'col_200_min_timelag' / &
         ANTags(1927)%Label / 'col_200_max_timelag' / &
         ANTags(1928)%Label / 'col_200_flag_threshold' / &
         ANTags(1929)%Label / 'col_200_error_value' /
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
         ACTags(790)%Label / 'col_101_variable' / &
         ACTags(791)%Label / 'col_101_useit' / &
         ACTags(792)%Label / 'col_101_measure_type' / &
         ACTags(793)%Label / 'col_101_instrument' / &
         ACTags(794)%Label / 'col_101_unit_in' / &
         ACTags(795)%Label / 'col_101_conversion' / &
         ACTags(796)%Label / 'col_101_unit_out' / &
         ACTags(797)%Label / 'col_102_variable' / &
         ACTags(798)%Label / 'col_102_useit' / &
         ACTags(799)%Label / 'col_102_measure_type' / &
         ACTags(800)%Label / 'col_102_instrument' /
    data ACTags(801)%Label / 'col_102_unit_in' / &
         ACTags(802)%Label / 'col_102_conversion' / &
         ACTags(803)%Label / 'col_102_unit_out' / &
         ACTags(804)%Label / 'col_103_variable' / &
         ACTags(805)%Label / 'col_103_useit' / &
         ACTags(806)%Label / 'col_103_measure_type' / &
         ACTags(807)%Label / 'col_103_instrument' / &
         ACTags(808)%Label / 'col_103_unit_in' / &
         ACTags(809)%Label / 'col_103_conversion' / &
         ACTags(810)%Label / 'col_103_unit_out' / &
         ACTags(811)%Label / 'col_104_variable' / &
         ACTags(812)%Label / 'col_104_useit' / &
         ACTags(813)%Label / 'col_104_measure_type' / &
         ACTags(814)%Label / 'col_104_instrument' / &
         ACTags(815)%Label / 'col_104_unit_in' / &
         ACTags(816)%Label / 'col_104_conversion' / &
         ACTags(817)%Label / 'col_104_unit_out' / &
         ACTags(818)%Label / 'col_105_variable' / &
         ACTags(819)%Label / 'col_105_useit' / &
         ACTags(820)%Label / 'col_105_measure_type' / &
         ACTags(821)%Label / 'col_105_instrument' / &
         ACTags(822)%Label / 'col_105_unit_in' / &
         ACTags(823)%Label / 'col_105_conversion' / &
         ACTags(824)%Label / 'col_105_unit_out' / &
         ACTags(825)%Label / 'col_106_variable' / &
         ACTags(826)%Label / 'col_106_useit' / &
         ACTags(827)%Label / 'col_106_measure_type' / &
         ACTags(828)%Label / 'col_106_instrument' / &
         ACTags(829)%Label / 'col_106_unit_in' / &
         ACTags(830)%Label / 'col_106_conversion' / &
         ACTags(831)%Label / 'col_106_unit_out' / &
         ACTags(832)%Label / 'col_107_variable' / &
         ACTags(833)%Label / 'col_107_useit' / &
         ACTags(834)%Label / 'col_107_measure_type' / &
         ACTags(835)%Label / 'col_107_instrument' / &
         ACTags(836)%Label / 'col_107_unit_in' / &
         ACTags(837)%Label / 'col_107_conversion' / &
         ACTags(838)%Label / 'col_107_unit_out' / &
         ACTags(839)%Label / 'col_108_variable' / &
         ACTags(840)%Label / 'col_108_useit' / &
         ACTags(841)%Label / 'col_108_measure_type' / &
         ACTags(842)%Label / 'col_108_instrument' / &
         ACTags(843)%Label / 'col_108_unit_in' / &
         ACTags(844)%Label / 'col_108_conversion' / &
         ACTags(845)%Label / 'col_108_unit_out' / &
         ACTags(846)%Label / 'col_109_variable' / &
         ACTags(847)%Label / 'col_109_useit' / &
         ACTags(848)%Label / 'col_109_measure_type' / &
         ACTags(849)%Label / 'col_109_instrument' / &
         ACTags(850)%Label / 'col_109_unit_in' / &
         ACTags(851)%Label / 'col_109_conversion' / &
         ACTags(852)%Label / 'col_109_unit_out' / &
         ACTags(853)%Label / 'col_110_variable' / &
         ACTags(854)%Label / 'col_110_useit' / &
         ACTags(855)%Label / 'col_110_measure_type' / &
         ACTags(856)%Label / 'col_110_instrument' / &
         ACTags(857)%Label / 'col_110_unit_in' / &
         ACTags(858)%Label / 'col_110_conversion' / &
         ACTags(859)%Label / 'col_110_unit_out' / &
         ACTags(860)%Label / 'col_111_variable' / &
         ACTags(861)%Label / 'col_111_useit' / &
         ACTags(862)%Label / 'col_111_measure_type' / &
         ACTags(863)%Label / 'col_111_instrument' / &
         ACTags(864)%Label / 'col_111_unit_in' / &
         ACTags(865)%Label / 'col_111_conversion' / &
         ACTags(866)%Label / 'col_111_unit_out' / &
         ACTags(867)%Label / 'col_112_variable' / &
         ACTags(868)%Label / 'col_112_useit' / &
         ACTags(869)%Label / 'col_112_measure_type' / &
         ACTags(870)%Label / 'col_112_instrument' / &
         ACTags(871)%Label / 'col_112_unit_in' / &
         ACTags(872)%Label / 'col_112_conversion' / &
         ACTags(873)%Label / 'col_112_unit_out' / &
         ACTags(874)%Label / 'col_113_variable' / &
         ACTags(875)%Label / 'col_113_useit' / &
         ACTags(876)%Label / 'col_113_measure_type' / &
         ACTags(877)%Label / 'col_113_instrument' / &
         ACTags(878)%Label / 'col_113_unit_in' / &
         ACTags(879)%Label / 'col_113_conversion' / &
         ACTags(880)%Label / 'col_113_unit_out' / &
         ACTags(881)%Label / 'col_114_variable' / &
         ACTags(882)%Label / 'col_114_useit' / &
         ACTags(883)%Label / 'col_114_measure_type' / &
         ACTags(884)%Label / 'col_114_instrument' / &
         ACTags(885)%Label / 'col_114_unit_in' / &
         ACTags(886)%Label / 'col_114_conversion' / &
         ACTags(887)%Label / 'col_114_unit_out' / &
         ACTags(888)%Label / 'col_115_variable' / &
         ACTags(889)%Label / 'col_115_useit' / &
         ACTags(890)%Label / 'col_115_measure_type' / &
         ACTags(891)%Label / 'col_115_instrument' / &
         ACTags(892)%Label / 'col_115_unit_in' / &
         ACTags(893)%Label / 'col_115_conversion' / &
         ACTags(894)%Label / 'col_115_unit_out' / &
         ACTags(895)%Label / 'col_116_variable' / &
         ACTags(896)%Label / 'col_116_useit' / &
         ACTags(897)%Label / 'col_116_measure_type' / &
         ACTags(898)%Label / 'col_116_instrument' / &
         ACTags(899)%Label / 'col_116_unit_in' / &
         ACTags(900)%Label / 'col_116_conversion' / &
         ACTags(901)%Label / 'col_116_unit_out' / &
         ACTags(902)%Label / 'col_117_variable' / &
         ACTags(903)%Label / 'col_117_useit' / &
         ACTags(904)%Label / 'col_117_measure_type' / &
         ACTags(905)%Label / 'col_117_instrument' / &
         ACTags(906)%Label / 'col_117_unit_in' / &
         ACTags(907)%Label / 'col_117_conversion' / &
         ACTags(908)%Label / 'col_117_unit_out' / &
         ACTags(909)%Label / 'col_118_variable' / &
         ACTags(910)%Label / 'col_118_useit' / &
         ACTags(911)%Label / 'col_118_measure_type' / &
         ACTags(912)%Label / 'col_118_instrument' / &
         ACTags(913)%Label / 'col_118_unit_in' / &
         ACTags(914)%Label / 'col_118_conversion' / &
         ACTags(915)%Label / 'col_118_unit_out' / &
         ACTags(916)%Label / 'col_119_variable' / &
         ACTags(917)%Label / 'col_119_useit' / &
         ACTags(918)%Label / 'col_119_measure_type' / &
         ACTags(919)%Label / 'col_119_instrument' / &
         ACTags(920)%Label / 'col_119_unit_in' / &
         ACTags(921)%Label / 'col_119_conversion' / &
         ACTags(922)%Label / 'col_119_unit_out' / &
         ACTags(923)%Label / 'col_120_variable' / &
         ACTags(924)%Label / 'col_120_useit' / &
         ACTags(925)%Label / 'col_120_measure_type' / &
         ACTags(926)%Label / 'col_120_instrument' / &
         ACTags(927)%Label / 'col_120_unit_in' / &
         ACTags(928)%Label / 'col_120_conversion' / &
         ACTags(929)%Label / 'col_120_unit_out' / &
         ACTags(930)%Label / 'col_121_variable' / &
         ACTags(931)%Label / 'col_121_useit' / &
         ACTags(932)%Label / 'col_121_measure_type' / &
         ACTags(933)%Label / 'col_121_instrument' / &
         ACTags(934)%Label / 'col_121_unit_in' / &
         ACTags(935)%Label / 'col_121_conversion' / &
         ACTags(936)%Label / 'col_121_unit_out' / &
         ACTags(937)%Label / 'col_122_variable' / &
         ACTags(938)%Label / 'col_122_useit' / &
         ACTags(939)%Label / 'col_122_measure_type' / &
         ACTags(940)%Label / 'col_122_instrument' / &
         ACTags(941)%Label / 'col_122_unit_in' / &
         ACTags(942)%Label / 'col_122_conversion' / &
         ACTags(943)%Label / 'col_122_unit_out' / &
         ACTags(944)%Label / 'col_123_variable' / &
         ACTags(945)%Label / 'col_123_useit' / &
         ACTags(946)%Label / 'col_123_measure_type' / &
         ACTags(947)%Label / 'col_123_instrument' / &
         ACTags(948)%Label / 'col_123_unit_in' / &
         ACTags(949)%Label / 'col_123_conversion' / &
         ACTags(950)%Label / 'col_123_unit_out' / &
         ACTags(951)%Label / 'col_124_variable' / &
         ACTags(952)%Label / 'col_124_useit' / &
         ACTags(953)%Label / 'col_124_measure_type' / &
         ACTags(954)%Label / 'col_124_instrument' / &
         ACTags(955)%Label / 'col_124_unit_in' / &
         ACTags(956)%Label / 'col_124_conversion' / &
         ACTags(957)%Label / 'col_124_unit_out' / &
         ACTags(958)%Label / 'col_125_variable' / &
         ACTags(959)%Label / 'col_125_useit' / &
         ACTags(960)%Label / 'col_125_measure_type' / &
         ACTags(961)%Label / 'col_125_instrument' / &
         ACTags(962)%Label / 'col_125_unit_in' / &
         ACTags(963)%Label / 'col_125_conversion' / &
         ACTags(964)%Label / 'col_125_unit_out' / &
         ACTags(965)%Label / 'col_126_variable' / &
         ACTags(966)%Label / 'col_126_useit' / &
         ACTags(967)%Label / 'col_126_measure_type' / &
         ACTags(968)%Label / 'col_126_instrument' / &
         ACTags(969)%Label / 'col_126_unit_in' / &
         ACTags(970)%Label / 'col_126_conversion' / &
         ACTags(971)%Label / 'col_126_unit_out' / &
         ACTags(972)%Label / 'col_127_variable' / &
         ACTags(973)%Label / 'col_127_useit' / &
         ACTags(974)%Label / 'col_127_measure_type' / &
         ACTags(975)%Label / 'col_127_instrument' / &
         ACTags(976)%Label / 'col_127_unit_in' / &
         ACTags(977)%Label / 'col_127_conversion' / &
         ACTags(978)%Label / 'col_127_unit_out' / &
         ACTags(979)%Label / 'col_128_variable' / &
         ACTags(980)%Label / 'col_128_useit' / &
         ACTags(981)%Label / 'col_128_measure_type' / &
         ACTags(982)%Label / 'col_128_instrument' / &
         ACTags(983)%Label / 'col_128_unit_in' / &
         ACTags(984)%Label / 'col_128_conversion' / &
         ACTags(985)%Label / 'col_128_unit_out' / &
         ACTags(986)%Label / 'col_129_variable' / &
         ACTags(987)%Label / 'col_129_useit' / &
         ACTags(988)%Label / 'col_129_measure_type' / &
         ACTags(989)%Label / 'col_129_instrument' / &
         ACTags(990)%Label / 'col_129_unit_in' / &
         ACTags(991)%Label / 'col_129_conversion' / &
         ACTags(992)%Label / 'col_129_unit_out' / &
         ACTags(993)%Label / 'col_130_variable' / &
         ACTags(994)%Label / 'col_130_useit' / &
         ACTags(995)%Label / 'col_130_measure_type' / &
         ACTags(996)%Label / 'col_130_instrument' / &
         ACTags(997)%Label / 'col_130_unit_in' / &
         ACTags(998)%Label / 'col_130_conversion' / &
         ACTags(999)%Label / 'col_130_unit_out' / &
         ACTags(1000)%Label / 'col_131_variable' /
    data ACTags(1001)%Label / 'col_131_useit' / &
         ACTags(1002)%Label / 'col_131_measure_type' / &
         ACTags(1003)%Label / 'col_131_instrument' / &
         ACTags(1004)%Label / 'col_131_unit_in' / &
         ACTags(1005)%Label / 'col_131_conversion' / &
         ACTags(1006)%Label / 'col_131_unit_out' / &
         ACTags(1007)%Label / 'col_132_variable' / &
         ACTags(1008)%Label / 'col_132_useit' / &
         ACTags(1009)%Label / 'col_132_measure_type' / &
         ACTags(1010)%Label / 'col_132_instrument' / &
         ACTags(1011)%Label / 'col_132_unit_in' / &
         ACTags(1012)%Label / 'col_132_conversion' / &
         ACTags(1013)%Label / 'col_132_unit_out' / &
         ACTags(1014)%Label / 'col_133_variable' / &
         ACTags(1015)%Label / 'col_133_useit' / &
         ACTags(1016)%Label / 'col_133_measure_type' / &
         ACTags(1017)%Label / 'col_133_instrument' / &
         ACTags(1018)%Label / 'col_133_unit_in' / &
         ACTags(1019)%Label / 'col_133_conversion' / &
         ACTags(1020)%Label / 'col_133_unit_out' / &
         ACTags(1021)%Label / 'col_134_variable' / &
         ACTags(1022)%Label / 'col_134_useit' / &
         ACTags(1023)%Label / 'col_134_measure_type' / &
         ACTags(1024)%Label / 'col_134_instrument' / &
         ACTags(1025)%Label / 'col_134_unit_in' / &
         ACTags(1026)%Label / 'col_134_conversion' / &
         ACTags(1027)%Label / 'col_134_unit_out' / &
         ACTags(1028)%Label / 'col_135_variable' / &
         ACTags(1029)%Label / 'col_135_useit' / &
         ACTags(1030)%Label / 'col_135_measure_type' / &
         ACTags(1031)%Label / 'col_135_instrument' / &
         ACTags(1032)%Label / 'col_135_unit_in' / &
         ACTags(1033)%Label / 'col_135_conversion' / &
         ACTags(1034)%Label / 'col_135_unit_out' / &
         ACTags(1035)%Label / 'col_136_variable' / &
         ACTags(1036)%Label / 'col_136_useit' / &
         ACTags(1037)%Label / 'col_136_measure_type' / &
         ACTags(1038)%Label / 'col_136_instrument' / &
         ACTags(1039)%Label / 'col_136_unit_in' / &
         ACTags(1040)%Label / 'col_136_conversion' / &
         ACTags(1041)%Label / 'col_136_unit_out' / &
         ACTags(1042)%Label / 'col_137_variable' / &
         ACTags(1043)%Label / 'col_137_useit' / &
         ACTags(1044)%Label / 'col_137_measure_type' / &
         ACTags(1045)%Label / 'col_137_instrument' / &
         ACTags(1046)%Label / 'col_137_unit_in' / &
         ACTags(1047)%Label / 'col_137_conversion' / &
         ACTags(1048)%Label / 'col_137_unit_out' / &
         ACTags(1049)%Label / 'col_138_variable' / &
         ACTags(1050)%Label / 'col_138_useit' / &
         ACTags(1051)%Label / 'col_138_measure_type' / &
         ACTags(1052)%Label / 'col_138_instrument' / &
         ACTags(1053)%Label / 'col_138_unit_in' / &
         ACTags(1054)%Label / 'col_138_conversion' / &
         ACTags(1055)%Label / 'col_138_unit_out' / &
         ACTags(1056)%Label / 'col_139_variable' / &
         ACTags(1057)%Label / 'col_139_useit' / &
         ACTags(1058)%Label / 'col_139_measure_type' / &
         ACTags(1059)%Label / 'col_139_instrument' / &
         ACTags(1060)%Label / 'col_139_unit_in' / &
         ACTags(1061)%Label / 'col_139_conversion' / &
         ACTags(1062)%Label / 'col_139_unit_out' / &
         ACTags(1063)%Label / 'col_140_variable' / &
         ACTags(1064)%Label / 'col_140_useit' / &
         ACTags(1065)%Label / 'col_140_measure_type' / &
         ACTags(1066)%Label / 'col_140_instrument' / &
         ACTags(1067)%Label / 'col_140_unit_in' / &
         ACTags(1068)%Label / 'col_140_conversion' / &
         ACTags(1069)%Label / 'col_140_unit_out' / &
         ACTags(1070)%Label / 'col_141_variable' / &
         ACTags(1071)%Label / 'col_141_useit' / &
         ACTags(1072)%Label / 'col_141_measure_type' / &
         ACTags(1073)%Label / 'col_141_instrument' / &
         ACTags(1074)%Label / 'col_141_unit_in' / &
         ACTags(1075)%Label / 'col_141_conversion' / &
         ACTags(1076)%Label / 'col_141_unit_out' / &
         ACTags(1077)%Label / 'col_142_variable' / &
         ACTags(1078)%Label / 'col_142_useit' / &
         ACTags(1079)%Label / 'col_142_measure_type' / &
         ACTags(1080)%Label / 'col_142_instrument' / &
         ACTags(1081)%Label / 'col_142_unit_in' / &
         ACTags(1082)%Label / 'col_142_conversion' / &
         ACTags(1083)%Label / 'col_142_unit_out' / &
         ACTags(1084)%Label / 'col_143_variable' / &
         ACTags(1085)%Label / 'col_143_useit' / &
         ACTags(1086)%Label / 'col_143_measure_type' / &
         ACTags(1087)%Label / 'col_143_instrument' / &
         ACTags(1088)%Label / 'col_143_unit_in' / &
         ACTags(1089)%Label / 'col_143_conversion' / &
         ACTags(1090)%Label / 'col_143_unit_out' / &
         ACTags(1091)%Label / 'col_144_variable' / &
         ACTags(1092)%Label / 'col_144_useit' / &
         ACTags(1093)%Label / 'col_144_measure_type' / &
         ACTags(1094)%Label / 'col_144_instrument' / &
         ACTags(1095)%Label / 'col_144_unit_in' / &
         ACTags(1096)%Label / 'col_144_conversion' / &
         ACTags(1097)%Label / 'col_144_unit_out' / &
         ACTags(1098)%Label / 'col_145_variable' / &
         ACTags(1099)%Label / 'col_145_useit' / &
         ACTags(1100)%Label / 'col_145_measure_type' / &
         ACTags(1101)%Label / 'col_145_instrument' / &
         ACTags(1102)%Label / 'col_145_unit_in' / &
         ACTags(1103)%Label / 'col_145_conversion' / &
         ACTags(1104)%Label / 'col_145_unit_out' / &
         ACTags(1105)%Label / 'col_146_variable' / &
         ACTags(1106)%Label / 'col_146_useit' / &
         ACTags(1107)%Label / 'col_146_measure_type' / &
         ACTags(1108)%Label / 'col_146_instrument' / &
         ACTags(1109)%Label / 'col_146_unit_in' / &
         ACTags(1110)%Label / 'col_146_conversion' / &
         ACTags(1111)%Label / 'col_146_unit_out' / &
         ACTags(1112)%Label / 'col_147_variable' / &
         ACTags(1113)%Label / 'col_147_useit' / &
         ACTags(1114)%Label / 'col_147_measure_type' / &
         ACTags(1115)%Label / 'col_147_instrument' / &
         ACTags(1116)%Label / 'col_147_unit_in' / &
         ACTags(1117)%Label / 'col_147_conversion' / &
         ACTags(1118)%Label / 'col_147_unit_out' / &
         ACTags(1119)%Label / 'col_148_variable' / &
         ACTags(1120)%Label / 'col_148_useit' / &
         ACTags(1121)%Label / 'col_148_measure_type' / &
         ACTags(1122)%Label / 'col_148_instrument' / &
         ACTags(1123)%Label / 'col_148_unit_in' / &
         ACTags(1124)%Label / 'col_148_conversion' / &
         ACTags(1125)%Label / 'col_148_unit_out' / &
         ACTags(1126)%Label / 'col_149_variable' / &
         ACTags(1127)%Label / 'col_149_useit' / &
         ACTags(1128)%Label / 'col_149_measure_type' / &
         ACTags(1129)%Label / 'col_149_instrument' / &
         ACTags(1130)%Label / 'col_149_unit_in' / &
         ACTags(1131)%Label / 'col_149_conversion' / &
         ACTags(1132)%Label / 'col_149_unit_out' / &
         ACTags(1133)%Label / 'col_150_variable' / &
         ACTags(1134)%Label / 'col_150_useit' / &
         ACTags(1135)%Label / 'col_150_measure_type' / &
         ACTags(1136)%Label / 'col_150_instrument' / &
         ACTags(1137)%Label / 'col_150_unit_in' / &
         ACTags(1138)%Label / 'col_150_conversion' / &
         ACTags(1139)%Label / 'col_150_unit_out' / &
         ACTags(1140)%Label / 'col_151_variable' / &
         ACTags(1141)%Label / 'col_151_useit' / &
         ACTags(1142)%Label / 'col_151_measure_type' / &
         ACTags(1143)%Label / 'col_151_instrument' / &
         ACTags(1144)%Label / 'col_151_unit_in' / &
         ACTags(1145)%Label / 'col_151_conversion' / &
         ACTags(1146)%Label / 'col_151_unit_out' / &
         ACTags(1147)%Label / 'col_152_variable' / &
         ACTags(1148)%Label / 'col_152_useit' / &
         ACTags(1149)%Label / 'col_152_measure_type' / &
         ACTags(1150)%Label / 'col_152_instrument' / &
         ACTags(1151)%Label / 'col_152_unit_in' / &
         ACTags(1152)%Label / 'col_152_conversion' / &
         ACTags(1153)%Label / 'col_152_unit_out' / &
         ACTags(1154)%Label / 'col_153_variable' / &
         ACTags(1155)%Label / 'col_153_useit' / &
         ACTags(1156)%Label / 'col_153_measure_type' / &
         ACTags(1157)%Label / 'col_153_instrument' / &
         ACTags(1158)%Label / 'col_153_unit_in' / &
         ACTags(1159)%Label / 'col_153_conversion' / &
         ACTags(1160)%Label / 'col_153_unit_out' / &
         ACTags(1161)%Label / 'col_154_variable' / &
         ACTags(1162)%Label / 'col_154_useit' / &
         ACTags(1163)%Label / 'col_154_measure_type' / &
         ACTags(1164)%Label / 'col_154_instrument' / &
         ACTags(1165)%Label / 'col_154_unit_in' / &
         ACTags(1166)%Label / 'col_154_conversion' / &
         ACTags(1167)%Label / 'col_154_unit_out' / &
         ACTags(1168)%Label / 'col_155_variable' / &
         ACTags(1169)%Label / 'col_155_useit' / &
         ACTags(1170)%Label / 'col_155_measure_type' / &
         ACTags(1171)%Label / 'col_155_instrument' / &
         ACTags(1172)%Label / 'col_155_unit_in' / &
         ACTags(1173)%Label / 'col_155_conversion' / &
         ACTags(1174)%Label / 'col_155_unit_out' / &
         ACTags(1175)%Label / 'col_156_variable' / &
         ACTags(1176)%Label / 'col_156_useit' / &
         ACTags(1177)%Label / 'col_156_measure_type' / &
         ACTags(1178)%Label / 'col_156_instrument' / &
         ACTags(1179)%Label / 'col_156_unit_in' / &
         ACTags(1180)%Label / 'col_156_conversion' / &
         ACTags(1181)%Label / 'col_156_unit_out' / &
         ACTags(1182)%Label / 'col_157_variable' / &
         ACTags(1183)%Label / 'col_157_useit' / &
         ACTags(1184)%Label / 'col_157_measure_type' / &
         ACTags(1185)%Label / 'col_157_instrument' / &
         ACTags(1186)%Label / 'col_157_unit_in' / &
         ACTags(1187)%Label / 'col_157_conversion' / &
         ACTags(1188)%Label / 'col_157_unit_out' / &
         ACTags(1189)%Label / 'col_158_variable' / &
         ACTags(1190)%Label / 'col_158_useit' / &
         ACTags(1191)%Label / 'col_158_measure_type' / &
         ACTags(1192)%Label / 'col_158_instrument' / &
         ACTags(1193)%Label / 'col_158_unit_in' / &
         ACTags(1194)%Label / 'col_158_conversion' / &
         ACTags(1195)%Label / 'col_158_unit_out' / &
         ACTags(1196)%Label / 'col_159_variable' / &
         ACTags(1197)%Label / 'col_159_useit' / &
         ACTags(1198)%Label / 'col_159_measure_type' / &
         ACTags(1199)%Label / 'col_159_instrument' / &
         ACTags(1200)%Label / 'col_159_unit_in' /
    data ACTags(1201)%Label / 'col_159_conversion' / &
         ACTags(1202)%Label / 'col_159_unit_out' / &
         ACTags(1203)%Label / 'col_160_variable' / &
         ACTags(1204)%Label / 'col_160_useit' / &
         ACTags(1205)%Label / 'col_160_measure_type' / &
         ACTags(1206)%Label / 'col_160_instrument' / &
         ACTags(1207)%Label / 'col_160_unit_in' / &
         ACTags(1208)%Label / 'col_160_conversion' / &
         ACTags(1209)%Label / 'col_160_unit_out' / &
         ACTags(1210)%Label / 'col_161_variable' / &
         ACTags(1211)%Label / 'col_161_useit' / &
         ACTags(1212)%Label / 'col_161_measure_type' / &
         ACTags(1213)%Label / 'col_161_instrument' / &
         ACTags(1214)%Label / 'col_161_unit_in' / &
         ACTags(1215)%Label / 'col_161_conversion' / &
         ACTags(1216)%Label / 'col_161_unit_out' / &
         ACTags(1217)%Label / 'col_162_variable' / &
         ACTags(1218)%Label / 'col_162_useit' / &
         ACTags(1219)%Label / 'col_162_measure_type' / &
         ACTags(1220)%Label / 'col_162_instrument' / &
         ACTags(1221)%Label / 'col_162_unit_in' / &
         ACTags(1222)%Label / 'col_162_conversion' / &
         ACTags(1223)%Label / 'col_162_unit_out' / &
         ACTags(1224)%Label / 'col_163_variable' / &
         ACTags(1225)%Label / 'col_163_useit' / &
         ACTags(1226)%Label / 'col_163_measure_type' / &
         ACTags(1227)%Label / 'col_163_instrument' / &
         ACTags(1228)%Label / 'col_163_unit_in' / &
         ACTags(1229)%Label / 'col_163_conversion' / &
         ACTags(1230)%Label / 'col_163_unit_out' / &
         ACTags(1231)%Label / 'col_164_variable' / &
         ACTags(1232)%Label / 'col_164_useit' / &
         ACTags(1233)%Label / 'col_164_measure_type' / &
         ACTags(1234)%Label / 'col_164_instrument' / &
         ACTags(1235)%Label / 'col_164_unit_in' / &
         ACTags(1236)%Label / 'col_164_conversion' / &
         ACTags(1237)%Label / 'col_164_unit_out' / &
         ACTags(1238)%Label / 'col_165_variable' / &
         ACTags(1239)%Label / 'col_165_useit' / &
         ACTags(1240)%Label / 'col_165_measure_type' / &
         ACTags(1241)%Label / 'col_165_instrument' / &
         ACTags(1242)%Label / 'col_165_unit_in' / &
         ACTags(1243)%Label / 'col_165_conversion' / &
         ACTags(1244)%Label / 'col_165_unit_out' / &
         ACTags(1245)%Label / 'col_166_variable' / &
         ACTags(1246)%Label / 'col_166_useit' / &
         ACTags(1247)%Label / 'col_166_measure_type' / &
         ACTags(1248)%Label / 'col_166_instrument' / &
         ACTags(1249)%Label / 'col_166_unit_in' / &
         ACTags(1250)%Label / 'col_166_conversion' / &
         ACTags(1251)%Label / 'col_166_unit_out' / &
         ACTags(1252)%Label / 'col_167_variable' / &
         ACTags(1253)%Label / 'col_167_useit' / &
         ACTags(1254)%Label / 'col_167_measure_type' / &
         ACTags(1255)%Label / 'col_167_instrument' / &
         ACTags(1256)%Label / 'col_167_unit_in' / &
         ACTags(1257)%Label / 'col_167_conversion' / &
         ACTags(1258)%Label / 'col_167_unit_out' / &
         ACTags(1259)%Label / 'col_168_variable' / &
         ACTags(1260)%Label / 'col_168_useit' / &
         ACTags(1261)%Label / 'col_168_measure_type' / &
         ACTags(1262)%Label / 'col_168_instrument' / &
         ACTags(1263)%Label / 'col_168_unit_in' / &
         ACTags(1264)%Label / 'col_168_conversion' / &
         ACTags(1265)%Label / 'col_168_unit_out' / &
         ACTags(1266)%Label / 'col_169_variable' / &
         ACTags(1267)%Label / 'col_169_useit' / &
         ACTags(1268)%Label / 'col_169_measure_type' / &
         ACTags(1269)%Label / 'col_169_instrument' / &
         ACTags(1270)%Label / 'col_169_unit_in' / &
         ACTags(1271)%Label / 'col_169_conversion' / &
         ACTags(1272)%Label / 'col_169_unit_out' / &
         ACTags(1273)%Label / 'col_170_variable' / &
         ACTags(1274)%Label / 'col_170_useit' / &
         ACTags(1275)%Label / 'col_170_measure_type' / &
         ACTags(1276)%Label / 'col_170_instrument' / &
         ACTags(1277)%Label / 'col_170_unit_in' / &
         ACTags(1278)%Label / 'col_170_conversion' / &
         ACTags(1279)%Label / 'col_170_unit_out' / &
         ACTags(1280)%Label / 'col_171_variable' / &
         ACTags(1281)%Label / 'col_171_useit' / &
         ACTags(1282)%Label / 'col_171_measure_type' / &
         ACTags(1283)%Label / 'col_171_instrument' / &
         ACTags(1284)%Label / 'col_171_unit_in' / &
         ACTags(1285)%Label / 'col_171_conversion' / &
         ACTags(1286)%Label / 'col_171_unit_out' / &
         ACTags(1287)%Label / 'col_172_variable' / &
         ACTags(1288)%Label / 'col_172_useit' / &
         ACTags(1289)%Label / 'col_172_measure_type' / &
         ACTags(1290)%Label / 'col_172_instrument' / &
         ACTags(1291)%Label / 'col_172_unit_in' / &
         ACTags(1292)%Label / 'col_172_conversion' / &
         ACTags(1293)%Label / 'col_172_unit_out' / &
         ACTags(1294)%Label / 'col_173_variable' / &
         ACTags(1295)%Label / 'col_173_useit' / &
         ACTags(1296)%Label / 'col_173_measure_type' / &
         ACTags(1297)%Label / 'col_173_instrument' / &
         ACTags(1298)%Label / 'col_173_unit_in' / &
         ACTags(1299)%Label / 'col_173_conversion' / &
         ACTags(1300)%Label / 'col_173_unit_out' / &
         ACTags(1301)%Label / 'col_174_variable' / &
         ACTags(1302)%Label / 'col_174_useit' / &
         ACTags(1303)%Label / 'col_174_measure_type' / &
         ACTags(1304)%Label / 'col_174_instrument' / &
         ACTags(1305)%Label / 'col_174_unit_in' / &
         ACTags(1306)%Label / 'col_174_conversion' / &
         ACTags(1307)%Label / 'col_174_unit_out' / &
         ACTags(1308)%Label / 'col_175_variable' / &
         ACTags(1309)%Label / 'col_175_useit' / &
         ACTags(1310)%Label / 'col_175_measure_type' / &
         ACTags(1311)%Label / 'col_175_instrument' / &
         ACTags(1312)%Label / 'col_175_unit_in' / &
         ACTags(1313)%Label / 'col_175_conversion' / &
         ACTags(1314)%Label / 'col_175_unit_out' / &
         ACTags(1315)%Label / 'col_176_variable' / &
         ACTags(1316)%Label / 'col_176_useit' / &
         ACTags(1317)%Label / 'col_176_measure_type' / &
         ACTags(1318)%Label / 'col_176_instrument' / &
         ACTags(1319)%Label / 'col_176_unit_in' / &
         ACTags(1320)%Label / 'col_176_conversion' / &
         ACTags(1321)%Label / 'col_176_unit_out' / &
         ACTags(1322)%Label / 'col_177_variable' / &
         ACTags(1323)%Label / 'col_177_useit' / &
         ACTags(1324)%Label / 'col_177_measure_type' / &
         ACTags(1325)%Label / 'col_177_instrument' / &
         ACTags(1326)%Label / 'col_177_unit_in' / &
         ACTags(1327)%Label / 'col_177_conversion' / &
         ACTags(1328)%Label / 'col_177_unit_out' / &
         ACTags(1329)%Label / 'col_178_variable' / &
         ACTags(1330)%Label / 'col_178_useit' / &
         ACTags(1331)%Label / 'col_178_measure_type' / &
         ACTags(1332)%Label / 'col_178_instrument' / &
         ACTags(1333)%Label / 'col_178_unit_in' / &
         ACTags(1334)%Label / 'col_178_conversion' / &
         ACTags(1335)%Label / 'col_178_unit_out' / &
         ACTags(1336)%Label / 'col_179_variable' / &
         ACTags(1337)%Label / 'col_179_useit' / &
         ACTags(1338)%Label / 'col_179_measure_type' / &
         ACTags(1339)%Label / 'col_179_instrument' / &
         ACTags(1340)%Label / 'col_179_unit_in' / &
         ACTags(1341)%Label / 'col_179_conversion' / &
         ACTags(1342)%Label / 'col_179_unit_out' / &
         ACTags(1343)%Label / 'col_180_variable' / &
         ACTags(1344)%Label / 'col_180_useit' / &
         ACTags(1345)%Label / 'col_180_measure_type' / &
         ACTags(1346)%Label / 'col_180_instrument' / &
         ACTags(1347)%Label / 'col_180_unit_in' / &
         ACTags(1348)%Label / 'col_180_conversion' / &
         ACTags(1349)%Label / 'col_180_unit_out' / &
         ACTags(1350)%Label / 'col_181_variable' / &
         ACTags(1351)%Label / 'col_181_useit' / &
         ACTags(1352)%Label / 'col_181_measure_type' / &
         ACTags(1353)%Label / 'col_181_instrument' / &
         ACTags(1354)%Label / 'col_181_unit_in' / &
         ACTags(1355)%Label / 'col_181_conversion' / &
         ACTags(1356)%Label / 'col_181_unit_out' / &
         ACTags(1357)%Label / 'col_182_variable' / &
         ACTags(1358)%Label / 'col_182_useit' / &
         ACTags(1359)%Label / 'col_182_measure_type' / &
         ACTags(1360)%Label / 'col_182_instrument' / &
         ACTags(1361)%Label / 'col_182_unit_in' / &
         ACTags(1362)%Label / 'col_182_conversion' / &
         ACTags(1363)%Label / 'col_182_unit_out' / &
         ACTags(1364)%Label / 'col_183_variable' / &
         ACTags(1365)%Label / 'col_183_useit' / &
         ACTags(1366)%Label / 'col_183_measure_type' / &
         ACTags(1367)%Label / 'col_183_instrument' / &
         ACTags(1368)%Label / 'col_183_unit_in' / &
         ACTags(1369)%Label / 'col_183_conversion' / &
         ACTags(1370)%Label / 'col_183_unit_out' / &
         ACTags(1371)%Label / 'col_184_variable' / &
         ACTags(1372)%Label / 'col_184_useit' / &
         ACTags(1373)%Label / 'col_184_measure_type' / &
         ACTags(1374)%Label / 'col_184_instrument' / &
         ACTags(1375)%Label / 'col_184_unit_in' / &
         ACTags(1376)%Label / 'col_184_conversion' / &
         ACTags(1377)%Label / 'col_184_unit_out' / &
         ACTags(1378)%Label / 'col_185_variable' / &
         ACTags(1379)%Label / 'col_185_useit' / &
         ACTags(1380)%Label / 'col_185_measure_type' / &
         ACTags(1381)%Label / 'col_185_instrument' / &
         ACTags(1382)%Label / 'col_185_unit_in' / &
         ACTags(1383)%Label / 'col_185_conversion' / &
         ACTags(1384)%Label / 'col_185_unit_out' / &
         ACTags(1385)%Label / 'col_186_variable' / &
         ACTags(1386)%Label / 'col_186_useit' / &
         ACTags(1387)%Label / 'col_186_measure_type' / &
         ACTags(1388)%Label / 'col_186_instrument' / &
         ACTags(1389)%Label / 'col_186_unit_in' / &
         ACTags(1390)%Label / 'col_186_conversion' / &
         ACTags(1391)%Label / 'col_186_unit_out' / &
         ACTags(1392)%Label / 'col_187_variable' / &
         ACTags(1393)%Label / 'col_187_useit' / &
         ACTags(1394)%Label / 'col_187_measure_type' / &
         ACTags(1395)%Label / 'col_187_instrument' / &
         ACTags(1396)%Label / 'col_187_unit_in' / &
         ACTags(1397)%Label / 'col_187_conversion' / &
         ACTags(1398)%Label / 'col_187_unit_out' / &
         ACTags(1399)%Label / 'col_188_variable' / &
         ACTags(1400)%Label / 'col_188_useit' /
    data ACTags(1401)%Label / 'col_188_measure_type' / &
         ACTags(1402)%Label / 'col_188_instrument' / &
         ACTags(1403)%Label / 'col_188_unit_in' / &
         ACTags(1404)%Label / 'col_188_conversion' / &
         ACTags(1405)%Label / 'col_188_unit_out' / &
         ACTags(1406)%Label / 'col_189_variable' / &
         ACTags(1407)%Label / 'col_189_useit' / &
         ACTags(1408)%Label / 'col_189_measure_type' / &
         ACTags(1409)%Label / 'col_189_instrument' / &
         ACTags(1410)%Label / 'col_189_unit_in' / &
         ACTags(1411)%Label / 'col_189_conversion' / &
         ACTags(1412)%Label / 'col_189_unit_out' / &
         ACTags(1413)%Label / 'col_190_variable' / &
         ACTags(1414)%Label / 'col_190_useit' / &
         ACTags(1415)%Label / 'col_190_measure_type' / &
         ACTags(1416)%Label / 'col_190_instrument' / &
         ACTags(1417)%Label / 'col_190_unit_in' / &
         ACTags(1418)%Label / 'col_190_conversion' / &
         ACTags(1419)%Label / 'col_190_unit_out' / &
         ACTags(1420)%Label / 'col_191_variable' / &
         ACTags(1421)%Label / 'col_191_useit' / &
         ACTags(1422)%Label / 'col_191_measure_type' / &
         ACTags(1423)%Label / 'col_191_instrument' / &
         ACTags(1424)%Label / 'col_191_unit_in' / &
         ACTags(1425)%Label / 'col_191_conversion' / &
         ACTags(1426)%Label / 'col_191_unit_out' / &
         ACTags(1427)%Label / 'col_192_variable' / &
         ACTags(1428)%Label / 'col_192_useit' / &
         ACTags(1429)%Label / 'col_192_measure_type' / &
         ACTags(1430)%Label / 'col_192_instrument' / &
         ACTags(1431)%Label / 'col_192_unit_in' / &
         ACTags(1432)%Label / 'col_192_conversion' / &
         ACTags(1433)%Label / 'col_192_unit_out' / &
         ACTags(1434)%Label / 'col_193_variable' / &
         ACTags(1435)%Label / 'col_193_useit' / &
         ACTags(1436)%Label / 'col_193_measure_type' / &
         ACTags(1437)%Label / 'col_193_instrument' / &
         ACTags(1438)%Label / 'col_193_unit_in' / &
         ACTags(1439)%Label / 'col_193_conversion' / &
         ACTags(1440)%Label / 'col_193_unit_out' / &
         ACTags(1441)%Label / 'col_194_variable' / &
         ACTags(1442)%Label / 'col_194_useit' / &
         ACTags(1443)%Label / 'col_194_measure_type' / &
         ACTags(1444)%Label / 'col_194_instrument' / &
         ACTags(1445)%Label / 'col_194_unit_in' / &
         ACTags(1446)%Label / 'col_194_conversion' / &
         ACTags(1447)%Label / 'col_194_unit_out' / &
         ACTags(1448)%Label / 'col_195_variable' / &
         ACTags(1449)%Label / 'col_195_useit' / &
         ACTags(1450)%Label / 'col_195_measure_type' / &
         ACTags(1451)%Label / 'col_195_instrument' / &
         ACTags(1452)%Label / 'col_195_unit_in' / &
         ACTags(1453)%Label / 'col_195_conversion' / &
         ACTags(1454)%Label / 'col_195_unit_out' / &
         ACTags(1455)%Label / 'col_196_variable' / &
         ACTags(1456)%Label / 'col_196_useit' / &
         ACTags(1457)%Label / 'col_196_measure_type' / &
         ACTags(1458)%Label / 'col_196_instrument' / &
         ACTags(1459)%Label / 'col_196_unit_in' / &
         ACTags(1460)%Label / 'col_196_conversion' / &
         ACTags(1461)%Label / 'col_196_unit_out' / &
         ACTags(1462)%Label / 'col_197_variable' / &
         ACTags(1463)%Label / 'col_197_useit' / &
         ACTags(1464)%Label / 'col_197_measure_type' / &
         ACTags(1465)%Label / 'col_197_instrument' / &
         ACTags(1466)%Label / 'col_197_unit_in' / &
         ACTags(1467)%Label / 'col_197_conversion' / &
         ACTags(1468)%Label / 'col_197_unit_out' / &
         ACTags(1469)%Label / 'col_198_variable' / &
         ACTags(1470)%Label / 'col_198_useit' / &
         ACTags(1471)%Label / 'col_198_measure_type' / &
         ACTags(1472)%Label / 'col_198_instrument' / &
         ACTags(1473)%Label / 'col_198_unit_in' / &
         ACTags(1474)%Label / 'col_198_conversion' / &
         ACTags(1475)%Label / 'col_198_unit_out' / &
         ACTags(1476)%Label / 'col_199_variable' / &
         ACTags(1477)%Label / 'col_199_useit' / &
         ACTags(1478)%Label / 'col_199_measure_type' / &
         ACTags(1479)%Label / 'col_199_instrument' / &
         ACTags(1480)%Label / 'col_199_unit_in' / &
         ACTags(1481)%Label / 'col_199_conversion' / &
         ACTags(1482)%Label / 'col_199_unit_out' / &
         ACTags(1483)%Label / 'col_200_variable' / &
         ACTags(1484)%Label / 'col_200_useit' / &
         ACTags(1485)%Label / 'col_200_measure_type' / &
         ACTags(1486)%Label / 'col_200_instrument' / &
         ACTags(1487)%Label / 'col_200_unit_in' / &
         ACTags(1488)%Label / 'col_200_conversion' / &
         ACTags(1489)%Label / 'col_200_unit_out' / &
         ACTags(1490)%Label / 'instr_9_manufacturer' /
    !> END GENERATED ACTags
end module m_common_global_var
