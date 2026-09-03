!***************************************************************************
! eddyflow-rp_main.f90
! -------------------
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
! \brief       Program main
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
program EddyFlowRP
    use m_rp_global_var
    use m_cec
    use m_pwb_timelag, only: ResetPwbDiagnostics, ReportPwbDiagnostics, InitPwbTimelagCache, &
        ReadPwbTimelagCache, WritePwbTimelagCache, SetPwbPeriodTimestamp, &
        PostProcessPwbTimelagCache, &
        RecordPwbTimelagOptPeriod, RebuildPwbTimelagOptFromCache, &
        ResetPwbAggregateSummary, AddPwbTimelagSummaryDataset, ResolvePwbAggregateSummary
    use m_ghg_prefetch, only: GhgPrefetchCleanup
    use m_prepass_parallel, only: PlanPrepassBatches, PrepassSlice, &
        StartPrepassBatches, WaitPrepassBatches, &
        WriteTlagBatchDump, MergeTlagBatchDumps, &
        WritePwbBatchDump, MergePwbBatchDumps, &
        WritePfBatchDump, MergePfBatchDumps
    !use netcdf
    !use iso_c_binding
    !use iso_fortran_env
    implicit none

    !> Local variables
    integer :: NumberOfOkPeriods = 0
    integer :: PeriodRecords = 0
    integer :: N2 = 0
    integer :: pcount = 0
    integer :: LatestRawFileIndx = 0
    integer :: NumRawFiles = 0
    integer :: NumberOfPeriods
    integer :: i
    integer :: j
    !> The hygrometer a gas is corrected with, from its own record.
    integer :: msl
    !> Mole fraction of the water a gas names, from whichever source.
    real(kind = dbl) :: chi_moist
    integer :: cec_p
    integer :: cec_k
    integer :: cec_ntarget
    integer :: cec_slots(MaxNumCecTargets)
    real(kind = dbl) :: cec_totals(MaxNumCecTargets)
    real(kind = dbl) :: cec_errors(MaxNumCecTargets)
    real(kind = dbl), allocatable :: CecPrimes(:, :)
    logical :: cec_ok
    integer :: SpecRow
    integer :: Nmax
    integer :: Nmin
    integer :: max_nsmpl
    integer :: pfn
    !> How the two assessment pre-passes were split across worker processes.
    !> 1 for a serial run, which is every run that did not ask for workers,
    !> every run too short to be worth splitting, and every PWB cache pre-pass.
    integer :: toWorkers
    integer :: pfWorkers
    integer :: sliceStart
    integer :: sliceEnd
    logical :: toParallel
    logical :: pfParallel
    !> Iterative correction: the pass counter, how many passes this project
    !> asked for, the worst relative change at the last comparison, and the
    !> previous pass's gas fluxes to compare against.
    integer :: corr_pass
    integer :: corr_passes
    real(kind = dbl) :: iter_dev
    real(kind = dbl) :: prev_gas_flux(GHGNumVar)
    real(kind = dbl), external :: WorstRelativeChange
    integer :: err_cnt1
    integer :: sec
    integer :: faulty_col
    integer :: STFlg(GHGNumVar)
    integer :: DTFlg(GHGNumVar)
    integer :: month
    integer :: day
    integer :: PeriodActualRecords
    integer :: ton
    integer :: PwbTimelagN
    integer :: int_doy
    integer :: bLastFile
    integer :: bLastRec
    integer :: TotNumFile
    integer :: NumFileNoRecurse
    integer :: pfStartTimestampIndx
    integer :: pfEndTimestampIndx
    integer :: toStartTimestampIndx
    integer :: toEndTimestampIndx
    integer :: rpStartTimestampIndx
    integer :: rpEndTimestampIndx
    integer :: tlagn(E2NumVar)
    integer :: MaxNumFileRecords
    integer :: NextRawFileIndx
    integer :: InitGasCalRefCol(GHGNumVar)
    integer :: nCalibEvents
    integer :: NumDynRecords
    integer :: clean
    integer :: dirty
    integer :: latestCleaning
    integer :: NumBiometFiles
    integer :: mkdir_status
    integer :: del_status

    integer, allocatable :: toH2On(:)
    integer, allocatable :: pfNumElem(:)

    real(kind = dbl) :: MissingRecords
    real(kind = dbl) :: float_doy
    real(kind = dbl) :: PP(3, 3)
    real(kind = dbl) :: Mat(3, 3)
    real(kind = dbl) :: Mat2d(2, 2)
    real(kind = dbl) :: pfVec(3)
    real(kind = dbl) :: pfVec2d(2)
    real(kind = dbl) :: PFb2d(2, MaxNumWSect) = 0.d0

    real(kind = dbl), allocatable :: bf(:)
    real(kind = sgl), allocatable :: Raw(:, :)
    real(kind = dbl), allocatable :: E2Set(:, :)
    real(kind = dbl), allocatable :: E2Primes(:, :)
    real(kind = dbl), allocatable :: UserSet(:, :)
    real(kind = dbl), allocatable :: UserPrimes(:, :)
    real(kind = dbl), allocatable :: DiagSet(:, :)
    real(kind = dbl), allocatable :: SpecSet(:, :)
    real(kind = dbl), allocatable :: pfWindBySect(:, :, :)
    real(kind = dbl), allocatable :: pfWind(:, :)

    character(10) :: loggedDate
    character(10) :: date
    character(5) :: time
    character(PathLen) :: suffixOutString
    character(64) :: TmpString1
    character(32) :: char_doy
    character(10) :: tmpDate
    character(5) ::  tmpTime
    character(128) ::  PeriodSkipMessage

    logical :: skip_period
    !> Whether any configured gas needs the FCC-only spectral path.
    logical :: has_fcc_only_gas
    logical :: passed(32)
    logical :: MetaIsNeeded = .true.
    logical :: EmbBiometDataExist = .false.
    logical :: AddUserStatsHeader = .true.
    logical :: IniFileNotFound
    logical :: initialize
    logical :: initializeBiometOut
    logical :: initializeFluxnetOut
    logical :: InitializeStorage
    logical :: InitOutVarPresence
    logical :: SingMat
    logical :: make_dataset_common
    logical :: make_dataset_rp
    logical :: FilterWhat(E2NumVar)
    logical :: FileEndReached
    logical :: toInit
    logical :: BiometDataFound
    logical :: AssessmentOnly
    logical :: FakeGoPlanarFit(1)
    logical :: PwbCacheRecognized
    logical :: PwbCacheValid

    logical, allocatable :: GoPlanarFit(:)

    type (FileListType), allocatable :: RawFileList(:)
    type (FileListType), allocatable :: bFileList(:)
    type (DateType), allocatable :: RawTimeSeries(:)
    type (DateType), allocatable :: MasterTimeSeries(:)
    type (DateType) :: tsStart, tsEnd
    type (DateType) :: LastMetadataTimestamp
    type (DateType) :: tsDatasetStart
    type (DateType) :: tsDatasetEnd
    type (DateType) :: auxStartTimestamp
    type (DateType) :: auxEndTimestamp
    type (DateType) :: SelectedStartTimestamp
    type (DateType) :: SelectedEndTimestamp
    type (StatsType) :: PrevStats
    type (AmbientStateType) :: prevAmbient
    type (QCType) :: StDiff
    type (QCType) :: DtDiff
    type (ColType) :: BypassCol(MaxNumCol)
    type (TestType) :: auxTest
    type(TimeLagOptType), allocatable :: TimelagOpt(:)
    type(TimeLagOptType), allocatable :: PwbTimelagOpt(:)
    type(TimeLagDatasetType), allocatable :: toSet(:)
    integer :: TimelagOptSize = 0
    integer :: PwbTimelagOptSize = 0

    integer, external :: NumOfPeriods
    integer, external :: NumberOfFilesInSubperiod
    real(kind = dbl), external :: LaggedCovarianceNoError
    real(kind = dbl) , external :: Poly6
    integer, external :: CreateDir
    include '../src_common/interfaces.inc'


    !***************************************************************************
    !***************************************************************************
    !****** INITIALIZATION PART COMMON TO ALL SW COMPONENTS ********************
    !***************************************************************************
    !***************************************************************************
    !> Connect the log before anything is said. It has no name yet - the
    !> output folder is not known until the project is read - so it starts on a
    !> scratch file, and LogInit copies it across once there is somewhere to
    !> put it. Late instead, and every line up to that point would write to a
    !> stray fort.163 in whatever directory the run started in.
    call LogStart()

    call LogSay('')
    call LogSay(' *******************')
    call LogSay('  Executing EddyFlow ')
    call LogSay(' *******************')
    call LogSay('')

    app = rp_app

    !> Initialize environment
    call InitEnv()

    !> By detault, create FLUXNET output
    EddyFlowProj%out_fluxnet = .true.

    !> Read setup file
    call ReadIniRP('RawProcess')
    AssessmentOnly = RPsetup%pf_assessment_only .or. &
        RPsetup%tlag_assessment_only
    allocate(bf(Meth%spec%nbins + 1))

    !> A project with no analyser is a legitimate configuration - an
    !> anemometer still measures momentum, sensible heat and stability - but
    !> the gas half of every output file is error codes, and that is worth
    !> saying before the user goes looking for the cause. Here rather than in
    !> the period loop, so it is said once.
    if (EddyFlowProj%gas_num <= 0) call ExceptionHandler(108)

    !> Add run-mode tag to Timestamp_FilePadding
    call TagRunMode()
    call ResetPwbDiagnostics()
    call ResetPwbAggregateSummary()
    pwb_last_optimal_lag = error
    pwb_last_optimal_origin = 0
    pwb_has_previous = .false.

    !> EddyFlow Express settings
    if (EddyFlowProj%run_mode == 'express') call ConfigureForExpress()
    if (EddyFlowProj%run_mode == 'md_retrieval') call ConfigureForMdRetrieval()
    if (EddyFlowProj%fluxnet_mode) call ConfigureForFluxnet()

    !> Define message for skipped periods
    if (EddyFlowProj%run_mode /= 'md_retrieval') then
        PeriodSkipMessage = '   Flux averaging period processing time: '
    else
        PeriodSkipMessage = '  Metadata retrieving time: '
    end if

    !> Selects which datasets should be filled with error codes,
    !> based on user selection
    make_dataset_common = EddyFlowProj%make_dataset
    make_dataset_rp     = EddyFlowProj%make_dataset

    !> Selects which files to output, considering the selected
    !> spectral correction method
    if (EddyFlowProj%out_avrg_cosp &
        .or. EddyFlowProj%out_avrg_spec &
        .or. (EddyFlowProj%hf_meth /= 'none' &
        .and. EddyFlowProj%hf_meth /= 'moncrieff_97' &
        .and. EddyFlowProj%hf_meth /= 'massman_00')) then
        !> in this cases, passage is needed to FCC, so:
        !> don't output files, don't create dataset
        !> don't output metadata
        !> don't calculate spectral correction
        !> don't calculate fluxes 2/3
        !> don't calculate footprint
        EddyFlowProj%fcc_follows     = .true.
        !> FCC owns the public full output in this flow. RP still writes the
        !> parent FLUXNET essentials file that FCC uses, including gas4 values.
        EddyFlowProj%out_full        = .false.
        EddyFlowProj%out_md          = .false.
        make_dataset_common         = .false.
    else
        !> in this cases, does what selected by user
        EddyFlowProj%fcc_follows  = .false.
    end if

    !> If running in embedded mode, override some settings
    if (EddyFlowProj%run_env == 'embedded') call ConfigureForEmbedded()
    if (EddyFlowProj%run_env == 'embedded') RPsetup%out_st = .false.

    !> Create output directory if it does not exist, otherwise is silent
    mkdir_status = CreateDir('"' //trim(adjustl(Dir%main_out)) // '"')

    !> The run log, as soon as there is a folder to put it in. Everything said
    !> before this - the banner, the project file, any exception raised while
    !> reading it - was buffered and is flushed here.
    call InitRunLog()

    !> Check on filename template
    call tsValidateTemplate(EddyFlowProj%fname_template)

    !> Detect number of raw files and allocate RawFileList
    call NumberOfFilesInDir(Dir%main_in, '.'//EddyFlowProj%fext, .true., &
        EddyFlowProj%fname_template, TotNumFile, NumFileNoRecurse)

    if (RPsetup%recurse) then
        NumRawFiles = TotNumFile
    else
        NumRawFiles = NumFileNoRecurse
    end if
    allocate(RawFileList(NumRawFiles))

     !> Store names of data files in RawFileList
    call FileListByExt(Dir%main_in, '.'//EddyFlowProj%fext, .true., .true., &
        EddyFlowProj%fname_template, EddyFlowLog%iso_format, .true., &
        RPsetup%recurse, RawFileList, size(RawFileList), .true., indent0)

    if (EddyFlowProj%use_extmd_file) then
        !> If requested, read external metadata file \n
        !> This is the case with non-GHG files or with GHG files if user \n
        !> explicitly selects an alternative metadata file
        write(*,'(a)', advance = 'no') ' Reading alternative metadata file: "' &
            // AuxFile%metadata(1:len_trim(AuxFile%metadata)) // '"..'
        write(ulog,'(a)', advance = 'no') ' Reading alternative metadata file: "' &
            // AuxFile%metadata(1:len_trim(AuxFile%metadata)) // '"..'
        call ReadMetadataFile(Col, AuxFile%metadata, IniFileNotFound, .true.)
        if (IniFileNotFound) then
            write(*, *)
            write(ulog, *)
            call ExceptionHandler(22)
        end if
        !> Retrieve variables to be used (from EddyFlow project file) \n
        !> and define user-variables
        call DefineUsedVariables(Col)
        MetaIsNeeded = .false.
        call MetadataFileValidation(Col, passed, faulty_col)
        if (.not. passed(1)) then
            write(*, *)
            write(ulog, *)
            call InformOfMetadataProblem(passed, faulty_col)
            call ExceptionHandler(23)
        end if
        call LogSay(' Done.')
    else
        !> In case of standard GHG processing, without alternative metadata \n
        !> file one GHG file must be opened to read the metadata content for \n
        !> importing information that is necessary before looping on all \n
        !> GHG files.
        allocate(Raw(1, 1))
        BypassCol = NullCol
        if (EddyFlowProj%ftype == 'licor_ghg') then
            i = 1
            do while (i <= NumRawFiles)
                !> The preamble reads one file to learn the columns and
                !> stops, so there is nothing to fetch ahead.
                call ReadLicorGhgArchive(RawFileList(i)%path, -1, -1, Col, &
                    BypassCol, .true., .false., .false., &
                    EddyFlowProj%run_mode /= 'md_retrieval', &
                    Raw, size(Raw, 1), size(Raw, 2), skip_period, passed, &
                    faulty_col, PeriodRecords, FileEndReached, .false., '')
                if (.not. skip_period .and. passed(1)) exit
                i = i + 1
                call InformOfMetadataProblem(passed, faulty_col)
            end do
            if (skip_period .or. (.not. passed(1))) call ExceptionHandler(32)
        end if
        deallocate(Raw)
    end if

    !> MasterSonic-related settings
    call DetectMasterSonic(Col, NumCol)

    !> Override/adjust settings related to the MasterSonic
    call OverrideMasterSonicRelatedSettings()

    !> Now that metadata are read, can set avrg_len in case user didn't
    if (RPsetup%avrg_len <= 0) RPsetup%avrg_len = nint(Metadata%file_length)
    if (EddyFlowProj%run_mode == 'md_retrieval') &
        RPsetup%avrg_len = nint(Metadata%file_length)

    !> Adjust time constant for planar fit if needed
    if (Meth%det == 'ld') then
        !> If time constant is larger than flux averaging interval,
        !> limit time constant to flux averaging interval and notify
        if (RPsetup%Tconst > RPsetup%avrg_len) then
            call ExceptionHandler(91)
            RPsetup%Tconst = nint(RPsetup%avrg_len * 6d1)
        end if
        !> Default to avrg_len anyway
        if (RPsetup%Tconst <= 0) RPsetup%Tconst = nint(RPsetup%avrg_len * 6d1)
    end if

    !> Some convenient variables
    DatafileDateStep = DateType(0, 0, 0, 0, nint(Metadata%file_length))
    DateStep         = DateType(0, 0, 0, 0, RPsetup%avrg_len)
    MaxNumFileRecords   = nint(Metadata%file_length * 60d0 * Metadata%ac_freq)
    MaxPeriodNumRecords = nint(RPsetup%avrg_len     * 60d0 * Metadata%ac_freq)

    !> Remember bypass columns (or columns detected
    !> from reading a sample GHG file)
    BypassCol = Col

    !> Initialize external biomet data
    if (index(EddyFlowProj%biomet_data, 'ext_') /= 0) then
        if (EddyFlowProj%biomet_data == 'ext_dir') then
            call LogSay(' Reading external biomet file(s) from:')
            write(*,'(a)') '  ' // trim(adjustl(Dir%biomet))
            write(ulog,'(a)') '  ' // trim(adjustl(Dir%biomet))
            call NumberOfFilesInDir(Dir%biomet, &
                trim(adjustl(EddyFlowProj%biomet_tail)), &
                .false., 'none', TotNumFile, NumFileNoRecurse)
            if (EddyFlowProj%biomet_recurse) then
                NumBiometFiles = TotNumFile
            else
                NumBiometFiles = NumFileNoRecurse
            end if
            call LogSay(' Done.')
        else
            NumBiometFiles = 1
        end if
        if (.not. allocated(bFileList)) allocate(bFileList(NumBiometFiles))
        call InitExternalBiomet(bFileList, size(bFileList))
    else
        allocate(bFileList(1))
    end if

    !> Open biomet output file
    if (index(EddyFlowProj%biomet_data, 'ext_') /= 0 .and. nbVars > 0) &
        call InitBiometOut()

    !> Initialize dynamic metadata by reading the file
    !> and figuring out available variables
    if (EddyFlowProj%use_dynmd_file) call InitDynamicMetadata(NumDynRecords)

    !> Determine potential radiation, based on lat/long info from metadata file
    PotRad = PotentialRadiation(Metadata%lat)

    !> Initialize output files for "user" variables (non-sensitive variables)
    !> if at least one such variable exists
    if (NumUserVar > 0) call InitUserOutFiles()

    !> Retrieve timestamp array in chronological order and
    !> order RawFileList, also in chronological order
    call FilesInChronologicalOrder(RawFileList, size(RawFileList), &
        tsDatasetStart, tsDatasetEnd, '')

    !> Adjust Start/End timestamps to define the boundaries of the
    !> RawTimeSeries. Retrieve the beginning time of first file
    !> and end time of last file.
    if (EddyFlowLog%tstamp_end) then
        tsDatasetStart = tsDatasetStart - DatafileDateStep
    else
        tsDatasetEnd = tsDatasetEnd + DatafileDateStep
    end if
    call tsRoundToMinute(tsDatasetStart, RPsetup%avrg_len, 'earlier')

    !> Retrieve NumberOfPeriods and allocate RawTimeSeries
    NumberOfPeriods = NumOfPeriods(tsDatasetStart, tsDatasetEnd, DateStep)
    allocate(RawTimeSeries(NumberOfPeriods + 1))

    !> Create timestamp array for full dataset
    call CreateTimeSeries(tsDatasetStart, tsDatasetEnd, DateStep, &
        RawTimeSeries, size(RawTimeSeries), .true.)

    !> Check the dynamic metadata file for calibration data.
    !> If found, builds up time series of absorptance drifts
    if (DriftCorr%method /= 'none') then
        allocate(tsDrifts(NumberOfPeriods + 1))
        allocate(Calib(0:NumDynRecords))  !< elem. 0 is to alloc. start of period
        allocate(tmpCalib(0:NumDynRecords))
        if (EddyFlowProj%use_dynmd_file) &
            call driftRetrieveCalibrationEvents(nCalibEvents)
    end if

    !> Define exp-binned frequencies extending \n
    !> from f_min = 1/(Flux avrg length) Hz
    !> to f_max = AcFreq/2 (= Nyquist frequency) Hz
    !> It is defined a priori, constant for each data period regardless
    !> of the actual length of the averaging periods
    if (EddyFlowProj%run_mode /=  'md_retrieval') &
        call BinnedFrequencyVector(bf, Meth%spec%nbins, &
            RPsetup%avrg_len, Metadata%ac_freq)

    !> Allocate array containing all potential data:
    !> rows: all rows potentially needed for current period
    !> columns: all except ignored ones and flag columns
    if (.not. allocated(Raw)) allocate(Raw(MaxPeriodNumRecords, NumAllVar))

    !***************************************************************************
    !***************************************************************************
    !******* TIME LAG OPTIMIZATION IF REQUESTED ********************************
    !***************************************************************************
    !***************************************************************************

    if (AssessmentOnly) then
        call LogSay(' Auxiliary assessment-only session requested.')
        if (RPsetup%tlag_assessment_only) &
            call LogSay('  Time-lag optimization will be created.')
        if (RPsetup%pf_assessment_only) &
            call LogSay('  Planar-fit file will be created.')
        call LogSay('')
    end if

    !> Method 5 uses a per-period cache when an existing time-lag file was
    !> selected.  A legacy aggregate optimizer file keeps its established path.
    if (trim(adjustl(Meth%tlag)) == 'pwb' .and. PwbCacheUpdateRequested) then
        call ReadPwbTimelagCache(AuxFile%to, PwbCacheRecognized, PwbCacheValid)
        if (PwbCacheRecognized) then
            if (.not. PwbCacheValid) &
                error stop 'PWB half-hourly time-lag table could not be used; see the message above.'
            call LogSay(' PWB half-hourly time-lag table found, retrieving content..')
            call LogSay(' PWB mode: exact per-period reuse; missing periods will be detected.')
        else
            TimeLagOptSelected = .true.
            Meth%tlag = 'tlag_opt'
            call LogSay(' PWB mode: aggregate/RH-class time-lag reuse; PWB detection is disabled.')
        end if
    elseif (PwbCacheGenerate) then
        call InitPwbTimelagCache()
        call LogSay(' PWB mode: whole-run time-lag pre-pass, then production processing from its table.')
    elseif (trim(adjustl(Meth%tlag)) == 'pwb') then
        call InitPwbTimelagCache()
        call LogSay(' PWB mode: live detection during production processing.')
    end if
    if (Meth%tlag == 'pwb') then
        PwbTimelagOptSize = size(RawTimeSeries) - 1
        allocate(PwbTimelagOpt(PwbTimelagOptSize))
        PwbTimelagN = 0
    end if

    if ((trim(adjustl(Meth%tlag)) == 'tlag_opt' .or. PwbCacheGenerate) .and. &
        (.not. AssessmentOnly .or. RPsetup%tlag_assessment_only)) then
        if (.not. RPsetup%to_onthefly) then
            call ReadTimelagOptFile(TOSetup%h2o_nclass)
            if (TOSetup%h2o_nclass > 1) &
                TOSetup%h2o_class_size = floor(100d0 / TOSetup%h2o_nclass)
        else
            call LogSay(' Performing time-lag optimization:')

            if (PwbCacheGenerate .and. EddyFlowProj%subperiod) then
                !> PWB caches cover the actual requested output range, not the
                !> optional aggregate-optimizer assessment subperiod.
                call DateTimeToDateType(EddyFlowProj%start_date, EddyFlowProj%start_time, auxStartTimestamp)
                call DateTimeToDateType(EddyFlowProj%end_date, EddyFlowProj%end_time, auxEndTimestamp)
                call tsExtractSubperiodIndexes(RawTimeSeries, size(RawTimeSeries), auxStartTimestamp, &
                    auxEndTimestamp, toStartTimestampIndx, toEndTimestampIndx)
                toEndTimestampIndx = toEndTimestampIndx + 1
                if (toStartTimestampIndx == nint(error) .or. toEndTimestampIndx == nint(error)) &
                    call ExceptionHandler(49)
            elseif (TOSetup%subperiod .and. .not. PwbCacheGenerate) then
                !> Timestamps of start and end of time-lag optimization period
                call DateTimeToDateType(TOSetup%start_date, TOSetup%start_time, auxStartTimestamp)
                call DateTimeToDateType(TOSetup%end_date, TOSetup%end_time, auxEndTimestamp)

                !> In RawTimeSeries, detect indices of first and last files
                !> relevant to time-lag optimization
                call tsExtractSubperiodIndexes(RawTimeSeries, &
                    size(RawTimeSeries), auxStartTimestamp, auxEndTimestamp, &
                    toStartTimestampIndx, toEndTimestampIndx)
                toEndTimestampIndx = toEndTimestampIndx + 1

                if (toStartTimestampIndx == nint(error) &
                    .or. toEndTimestampIndx == nint(error)) &
                    call ExceptionHandler(49)
            else
                toStartTimestampIndx = 1
                toEndTimestampIndx = size(RawTimeSeries)
            end if

            !> Count maximum number of periods for timelag optimization
            write(TmpString1, '(i7)') toEndTimestampIndx - toStartTimestampIndx
            write(*, '(a)') '  Maximum number of flux averaging periods &
                &available for time-lag optimization: ' &
                // trim(adjustl(TmpString1))
            write(ulog, '(a)') '  Maximum number of flux averaging periods &
                &available for time-lag optimization: ' &
                // trim(adjustl(TmpString1))

            !> Allocate variables that depend upon maximum number of periods
            if (.not. PwbCacheGenerate) then
                TimelagOptSize = toEndTimestampIndx - toStartTimestampIndx
                allocate(TimelagOpt(TimelagOptSize))
            end if

            !> Loop on selected files and calculate relevant statistics
            ton = 0
            month = 0
            day   = 0
            pcount = toStartTimestampIndx - 1
            LatestRawFileIndx = 1
            bLastFile = 1
            bLastRec = 0
            DynamicMetadata = ErrDynamicMetadata
            LastMetadataTimestamp = DateType(0, 0, 0, 0, 0)
            toInit = .true.

            !> Every period in the range is read, reduced, and turned into one
            !> record that depends on no other period's. So the range can be
            !> cut into slices and each slice run in a copy of this program.
            !> PlanPrepassBatches decides whether that is worth doing - it
            !> never is for a worker, which would otherwise spawn workers of
            !> its own, nor for a range too short to pay for the processes.
            !> A PWB cache pre-pass may be split too, now. It could not while
            !> the streaming classifier's verdict reached the output: that
            !> verdict depends on the last settled detection before a period,
            !> and a slice starting cold has none. Three things carried it -
            !> the terminal fallback's lag, the aggregate dataset's membership
            !> and the donor tally - and each is now taken from the settled
            !> table, which is built once, by the parent, over every slice.
            !> What a worker produces is evidence, and evidence does not
            !> depend on where the walk began.
            call PlanPrepassBatches(toEndTimestampIndx - toStartTimestampIndx, &
                .true., toWorkers)

            !> A worker was handed its slice on the command line. Narrowed
            !> here rather than where the range was computed, so the arrays
            !> above are still allocated for the whole range and the same
            !> bounds checks hold in parent and worker alike.
            if (BatchIndex > 0) then
                toStartTimestampIndx = BatchSliceStart
                toEndTimestampIndx = BatchSliceEnd
                pcount = toStartTimestampIndx - 1
            end if

            !> The workers go first so they are already reading raw data while
            !> this process works through the first slice. The parent takes a
            !> slice rather than waiting: the code after the loop reads state
            !> the loop establishes, and a parent that had skipped it would
            !> reach that code with the state unset.
            toParallel = toWorkers > 1
            if (toParallel) then
                call StartPrepassBatches('to', toStartTimestampIndx, &
                    toEndTimestampIndx, toWorkers)
                call PrepassSlice(toStartTimestampIndx, toEndTimestampIndx, &
                    toWorkers, 1, sliceStart, sliceEnd)
                toStartTimestampIndx = sliceStart
                toEndTimestampIndx = sliceEnd
                pcount = toStartTimestampIndx - 1
            end if

            to_periods_loop: do
                pcount = pcount + 1

                !> If embedded metadata are to be used,
                !> reinitialize column information to null
                if (EddyFlowProj%use_extmd_file) then
                    Col = BypassCol
                else
                    Col = NullCol
                end if

                !> Normal exit instruction: either the last period was
                !> dealt with, or raw files are finished
                if (LatestRawFileIndx > NumRawFiles &
                    .or. pcount >= toEndTimestampIndx) exit to_periods_loop

                !> Define initial/final timestamps
                !> of current period, say [8:00 - 8:30)
                tsStart = RawTimeSeries(pcount)
                tsEnd   = RawTimeSeries(pcount + 1)
                call DateTypeToDateTime(tsEnd, tmpDate, tmpTime)
                call SetPwbPeriodTimestamp(tmpDate, tmpTime)

                !> Search file containing data starting from the
                !> time closest to tsStart. Searches only from most current
                !> file onward, to avoid wasting time
                call FirstFileOfCurrentPeriod(tsStart, tsEnd, RawFileList, &
                    NumRawFiles, LatestRawFileIndx, NextRawFileIndx, skip_period)

                !> Averaging period advancement
                 if (day /= 0) then
                    if (EddyFlowProj%caller == 'console') then
                        call LogSayNoAdv('#')
                    else
                        call DisplayProgress('avrg_interval', &
                            '   another small step to the time-lag: ', &
                            tsStart, 'yes')
                    end if
                end if

                !> Daily advancement
                if (day /= tsStart%day &
                    .or. month /= tsStart%month) then
                    month = tsStart%month
                    day   = tsStart%day
                    if (EddyFlowProj%caller == 'console') then
                        write(*, '(a)')
                        write(ulog, '(a)')
                        call DisplayProgress('daily','  Importing data for ', &
                            tsStart, 'no')
                    else
                        call DisplayProgress('daily','  Importing data for ', &
                            tsStart, 'yes')
                    end if
                end if

                if (skip_period) cycle to_periods_loop

                !> Import dataset for current period. If using embedded biomet,
                !> also read biomet data. On entrance, NextRawFileIndx contains
                !> the index of the file to start the current period with
                !> On exit, LatestRawFileIndx contains index of latest file used
                call ImportCurrentPeriod(tsStart, tsEnd, &
                    RawFileList, NumRawFiles, NextRawFileIndx, BypassCol, &
                    MaxNumFileRecords, MetaIsNeeded, &
                    EddyFlowProj%biomet_data == 'embedded', .false., &
                    Raw, size(Raw, 1), size(Raw, 2), PeriodRecords, &
                    EmbBiometDataExist, skip_period, LatestRawFileIndx, Col, &
                    .false.)

                if (skip_period) cycle to_periods_loop

                !> Period skip control with message
                MissingRecords = dfloat(MaxPeriodNumRecords - PeriodRecords) &
                    / dfloat(MaxPeriodNumRecords) * 100d0
                if (PeriodRecords > 0 .and. MissingRecords > RPsetup%max_lack) &
                    cycle to_periods_loop

                !> Filter raw data for user-defined flags
                if (RPsetup%filter_by_raw_flags) &
                    call FilterDatasetForFlags(Col, Raw, &
                        size(Raw, 1), size(Raw, 2))

                !***************************************************************
                !**** RAW FILE IMPORT FINISHES HERE ****************************
                !**** NOW STARTS DATASET DEFINITION ****************************
                !***************************************************************

                !> Allocate arrays for actual data processing
                if (.not. allocated(E2Set)) &
                    allocate(E2Set(PeriodRecords, E2NumVar))
                if (.not. allocated(E2Primes)) &
                    allocate(E2Primes(PeriodRecords, E2NumVar))
                if (.not. allocated(DiagSet)) &
                    allocate(DiagSet(PeriodRecords, MaxNumDiag))

                !> Define EddyFlow set of variables for the following processing
                call DefineE2Set(Col, Raw,   size(Raw, 1),     Size(Raw, 2), &
                                    E2Set,   size(E2Set, 1),   Size(E2Set, 2), &
                                    DiagSet, size(DiagSet, 1), Size(DiagSet, 2))

                !> If H2O instrument path type is 'open', doesn't make sense
                !> to use RH classes so set it to 1.
                if (toInit) then
                    if (E2Col(PrimaryWaterOutSlot())%instr%path_type == 'open') then
                        TOSetup%h2o_nclass = 1
                        toInit = .false.
                    end if
                end if

                !> Clean up E2Set, eliminating values that are clearly unphysical
                call CleanUpE2Set(E2Set, size(E2Set, 1), size(E2Set, 2))

                !> Define as not present, variables for which
                !> too many values are out-ranged
                call EliminateCorruptedVariables(E2Set, size(E2Set, 1), &
                    size(E2Set, 2), skip_period, .false.)

                !> If either u, v or w have been eliminated, stops processing this period
                if (skip_period) then
                    if(allocated(E2Set)) deallocate(E2Set)
                    if(allocated(E2Primes)) deallocate(E2Primes)
                    if(allocated(DiagSet)) deallocate(DiagSet)
                    cycle to_periods_loop
                end if

                !> Every configured gas, not the four historical slots: a
                !> project whose gases all sit past the fourth record would
                !> otherwise have every period discarded here as gas-less.
                if (.not. any(E2Col(firstGas:lastGas)%present)) then
                    if(allocated(E2Set)) deallocate(E2Set)
                    if(allocated(E2Primes)) deallocate(E2Primes)
                    if(allocated(DiagSet)) deallocate(DiagSet)
                    cycle to_periods_loop
                end if

                !> Update metadata if dynamic metadata are to be used
                if (EddyFlowProj%use_dynmd_file) &
                    call RetrieveDynamicMetadata(tsEnd, E2Col, size(E2Col))

                !> Retrieve biomet data if they exist
                if (index(EddyFlowProj%biomet_data, 'ext_') /= 0) then
                    call BiometRetrieveExternalData(bFileList, size(bFileList), &
                        bLastFile, bLastRec, tsStart, &
                        tsEnd, BiometDataFound, .false.)
                elseif (EddyFlowProj%biomet_data == 'embedded') then
                    call BiometRetrieveEmbeddedData(EmbBiometDataExist, .false.)
                end if

                !> Calculate relative separations between
                !> the analyzers and the anemometer used
                call DefineRelativeSeparations()

                !> Override users choices if needed
                call OverrideSettings()

                !***************************************************************
                !**** DATASET DEFINITION FINISHES HERE *************************
                !**** NOW STARTS RAW DATA REDUCTION ****************************
                !***************************************************************
                !> Interpret diagnostics and filter accordingly
                if (NumDiag > 0) then
                    call InterpretLicorDiagnostics(DiagSet, &
                        size(DiagSet, 1), size(DiagSet, 2))
                    call FilterDatasetForDiagnostics(E2Set, &
                        size(E2Set, 1), size(E2Set, 2), &
                        DiagSet, size(DiagSet, 1), size(DiagSet, 2), &
                        DiagAnemometer, .true.)
                end if
                if(allocated(DiagSet)) deallocate(DiagSet)

                !> Adjust coordinate systems if the case
                call AdjustSonicCoordinates(E2Set, size(E2Set, 1), size(E2Set, 2))

                !> Filter for wind direction if requested
                if (RPSetup%apply_wdf) &
                    call FilterDatasetForWindDirection(E2Set, size(E2Set, 1), size(E2Set, 2))

                !> Calculate basic stats
                call BasicStats(E2Set, size(E2Set, 1), size(E2Set, 2), 1, .false.)
                Stats1 = Stats

                !> Calculate raw screening flags and despike data if requeste
                auxTest = TestType(.true., .false., .false., .true., .false., &
                    .false., .false., .false., .false., .false., .false.)
                call StatisticalScreening(E2Set, &
                    size(E2Set, 1), size(E2Set, 2), auxTest, .false.)

                !> Define as not present, variables for which
                !> too many values are out-ranged
                call EliminateCorruptedVariables(E2Set, &
                    size(E2Set, 1), size(E2Set, 2), skip_period, .false.)

                !> If either u, v or w have been eliminated,
                !> stops processing this period
                if (skip_period) then
                    if(allocated(E2Set)) deallocate(E2Set)
                    if(allocated(E2Primes)) deallocate(E2Primes)
                    cycle to_periods_loop
                end if

                !> Calculate basic stats
                call BasicStats(E2Set, &
                    size(E2Set, 1), size(E2Set, 2), 2, .false.)
                Stats2 = Stats

                !> Apply raw-level cross wind correction
                !> (after Liu et al. 2001), if requested
                if (RPsetup%calib_cw) &
                    call CrossWindCorr(E2Col(u), E2Set, &
                        size(E2Set, 1), size(E2Set, 2), .false.)

                !> Calculate basic stats and output them as requested
                call BasicStats(E2Set, &
                    size(E2Set, 1), size(E2Set, 2), 3, .false.)
                Stats3 = Stats

                !> Angle-of-attack calibration
                call AoaCalibration(E2Set, size(E2Set, 1), size(E2Set, 2))

                !> Gill WindMaster w-boost
                if (RPsetup%calib_wboost) &
                    call ApplyGillWmWBoost(E2Set, size(E2Set, 1), size(E2Set, 2))

                !> Metek head correction, the same one the flux loop applies.
                !> Without it here, a lag or a plane would be worked out from
                !> a wind the fluxes never see. Its companion, the
                !> inclinometer correction, cannot follow: it reads its angles
                !> from the custom columns, and DefineUserSet has not run yet
                !> at this point in the program. Recorded rather than worked
                !> around, because moving DefineUserSet is a layout change and
                !> this correction has no effect on any current dataset.
                call MetekHeadCorrection(E2Set, size(E2Set, 1), size(E2Set, 2))

                !> Calculate basic stats
                call BasicStats(E2Set, &
                    size(E2Set, 1), size(E2Set, 2), 4, .false.)
                Stats4 = Stats

                !> Apply rotations for tilt correction, if requested
                FakeGoPlanarFit = .false.
                call TiltCorrection('double_rotation', FakeGoPlanarFit, E2Set, &
                    size(E2Set, 1), size(E2Set, 2), 1, Essentials%yaw, &
                    Essentials%pitch, Essentials%roll, .false.)

                !> Calculate basic stats
                call BasicStats(E2Set, &
                    size(E2Set, 1), size(E2Set, 2), 5, .false.)
                Stats5 = Stats

                !> PWB detection on rotated, pre-WPL data. This used to be
                !> optional; both alternatives ran on rotated 20 Hz data, so
                !> the choice was only ever whether the gas series had been
                !> through the pointwise mixing-ratio conversion first. That
                !> conversion runs before time-lag compensation, so after it
                !> the cell temperature and water signals sit in the gas
                !> series at the wrong relative lag - which is the one series
                !> being cross-correlated here.
                if (PwbCacheGenerate) then
                    call RetrieveSensorParams()
                    call SetTimelags()
                    pwb_detect_only_mode = .true.
                    call TimeLagHandle('pwb', E2Set, size(E2Set, 1), size(E2Set, 2), &
                        pwb_raw_ActTLag, pwb_raw_TLag, pwb_raw_DefTlagUsed, .false.)
                    pwb_raw_Result = PWBResult
                    pwb_raw_detection_done = .true.
                end if

                !> Convert to mixing ratios (if requested and if the case)
                if (EddyFlowProj%wpl) &
                    call PointByPointToMixingRatio(E2Set, &
                        size(E2Set, 1), size(E2Set, 2), .false.)

                !> Initialize analytic spectral corrections,
                !> retrieving sensor parameters
                call RetrieveSensorParams()

                if (PwbCacheGenerate) then
                    call SetTimelags()
                    call TimeLagHandle('pwb', E2Set, &
                        size(E2Set, 1), size(E2Set, 2), Essentials%actual_timelag, &
                        Essentials%used_timelag, Essentials%def_tlag, .true.)
                else
                    !> Adjust min/max time-lags associated to columns, to fit
                    !> user settings in the Time lag optimizer dialog
                    call AdjustTimelagOptSettings()
                    call TimeLagHandle('maxcov', E2Set, &
                        size(E2Set, 1), size(E2Set, 2), Essentials%actual_timelag, &
                        Essentials%used_timelag, Essentials%def_tlag, .true.)
                end if

                !> Calculate basic stats
                call BasicStats(E2Set, &
                    size(E2Set, 1), size(E2Set, 2), 6, .false.)
                Stats6 = Stats

                !> Calculate air and cell parameters
                call AirAndCellParameters()

                !> Apply filter for absolute limits test, if the case
                FilterWhat = .false.
                FilterWhat(firstGas:lastGas) = .true.
                call FilterDatasetForPhysicalThresholds(E2Set, &
                    size(E2Set, 1), size(E2Set, 2), FilterWhat)

                !> Define as not present, variables for which
                !> too many values are out-ranged
                call EliminateCorruptedVariables(E2Set, &
                    size(E2Set, 1), size(E2Set, 2), skip_period, .false.)

                !> If either u, v or w have been eliminated,
                !> stops processing this period
                if (skip_period) then
                    if(allocated(E2Set)) deallocate(E2Set)
                    if(allocated(E2Primes)) deallocate(E2Primes)
                    cycle to_periods_loop
                end if

                call Fluctuations(E2Set, E2Primes, size(E2Set, 1), &
                    size(E2Set, 2), RPsetup%Tconst, Stats6, E2Col)
                if (allocated(E2Set)) deallocate(E2Set)

                !> Calculate basic stats
                call BasicStats(E2Primes, &
                    size(E2Primes, 1), size(E2Primes, 2), 7, .false.)
                if (allocated(E2Primes)) deallocate(E2Primes)
                Stats7 = Stats

                !***************************************************************
                !**** RAW DATA REDUCTION FINISHES HERE. ************************
                !**** TENTATIVE FLUX CALCULATION        ************************
                !***************************************************************

                !> Average mole fractions in [umol mol_a-1] and [mmol mol_a-1]
                call MoleFractionsAndMixingRatios()

                !> Calculate parameters for flux computation
                call FluxParams(.false.)

                !> Calculate fluxes at Level 0
                call Fluxes0_rp(.false.)

                if (PwbCacheGenerate) then
                    PwbTimelagN = PwbTimelagN + 1
                    if (.not. allocated(PwbTimelagOpt) .or. PwbTimelagOptSize <= 0 &
                        .or. PwbTimelagN > PwbTimelagOptSize) &
                        error stop 'PWB time-lag optimization dataset is not allocated safely.'
                    !> Only what the table cannot say - the humidity and
                    !> which period this is. The lags come from the settled
                    !> table afterwards, not from the streaming guess that
                    !> this period has just been given.
                    call RecordPwbTimelagOptPeriod(PwbTimelagOpt, &
                        PwbTimelagOptSize, PwbTimelagN)
                else
                    !> Store values only for aggregate time-lag optimization.
                    ton = ton + 1
                    if (.not. allocated(TimelagOpt) .or. TimelagOptSize <= 0 &
                        .or. ton > TimelagOptSize) &
                        error stop 'Time-lag optimization dataset is not allocated safely.'
                    call AddToTimelagOptDataset(TimelagOpt, TimelagOptSize, ton)
                end if

            end do to_periods_loop
            write(*, '(a)')
            write(ulog, '(a)')
            call LogSay(' Done.')

            !> Now collect what the other slices produced and append them to
            !> this one, which leaves the dataset in period order - the order
            !> a single loop over the whole range would have built it in.
            if (toParallel) then
                call WaitPrepassBatches('to', toWorkers)
                if (PwbCacheGenerate) then
                    call MergePwbBatchDumps('to', toWorkers, PwbTimelagOpt, &
                        PwbTimelagOptSize, PwbTimelagN)
                else
                    call MergeTlagBatchDumps('to', toWorkers, TimelagOpt, &
                        TimelagOptSize, ton)
                end if
            end if

            !> A worker's job ends here: it hands back the records its slice
            !> produced and stops. The fit itself - OptimizeTimelags, the
            !> cache post-processing - is done once, by the parent, over every
            !> slice concatenated, so a worker doing it too would be both
            !> wasted work and a second answer nobody reads.
            if (BatchIndex > 0) then
                if (PwbCacheGenerate) then
                    call WritePwbBatchDump(PwbTimelagOpt, &
                        PwbTimelagOptSize, PwbTimelagN)
                else
                    call WriteTlagBatchDump(TimelagOpt, TimelagOptSize, ton)
                end if
                call LogSay(' Time-lag pre-pass slice finished.')
                stop ''
            end if

            !*******************************************************************
            !**** RAW DATA REDUCTION FINISHES HERE.    *************************
            !**** NOW STARTS TIME LAG OPT CALCULATIONS *************************
            !*******************************************************************

            if (PwbCacheGenerate) then
                !> Every period has been read, so the gap-filling can look
                !> forwards as well as back. This is what the streaming
                !> classifier in timelag_handle cannot do.
                call PostProcessPwbTimelagCache()
                !> The aggregate dataset is built here, from the finished
                !> table, rather than accumulated as the walk went.
                call RebuildPwbTimelagOptFromCache(PwbTimelagOpt, &
                    PwbTimelagOptSize, PwbTimelagN)
                call WritePwbTimelagCache()
                if (PwbTimelagN > 0) then
                    if (.not. allocated(PwbTimelagOpt) .or. PwbTimelagOptSize <= 0 &
                        .or. PwbTimelagN > PwbTimelagOptSize) &
                        error stop 'PWB time-lag optimization dataset is not allocated safely.'
                    allocate(toSet(PwbTimelagN))
                    call FixTimelagOptDataset(PwbTimelagOpt, PwbTimelagOptSize, &
                        toSet, size(toSet), tlagn, size(tlagn))
                    allocate(toH2On(TOSetup%h2o_nclass))
                    call OptimizeTimelags(toSet, size(toSet), tlagn, E2NumVar, toH2On, &
                        TOSetup%h2o_nclass, TOSetup%h2o_class_size)
                    call ResolvePwbAggregateSummary(tlagn)
                    PwbAggregateSummary = .true.
                    call WriteOutTimelagOptimization(tlagn, E2NumVar, toH2On, &
                        TOSetup%h2o_nclass, TOSetup%h2o_class_size)
                    PwbAggregateSummary = .false.
                    deallocate(toSet)
                    deallocate(toH2On)
                    PwbTimelagN = 0
                    call ResetPwbAggregateSummary()
                end if
                call ReportPwbDiagnostics()
                call ResetPwbDiagnostics()
                call LogSay(' PWB time-lag pre-pass finished.')
            else
            !> Adjust time-lag opt dataset to eliminate errors,
            !> so that it's easier to treat them later
            allocate (toSet(ton))
            call FixTimelagOptDataset(TimelagOpt, TimelagOptSize, &
                toSet, size(toSet), tlagn, size(tlagn))
            if (allocated(TimelagOpt)) deallocate(TimelagOpt)
            TimelagOptSize = 0

            allocate(toH2On(TOSetup%h2o_nclass))

            !> Optimize time-lags                                        ******* Improve readability of this subroutine interface
            call OptimizeTimelags(toSet, size(toSet), tlagn, E2NumVar, toH2On, & 
                TOSetup%h2o_nclass, TOSetup%h2o_class_size)

            !> Write time-lag optimization results on output file
            if (.not. (Meth%tlag == 'maxcov')) &
                call WriteOutTimelagOptimization(tlagn, E2NumVar, &
                    toH2On, TOSetup%h2o_nclass, TOSetup%h2o_class_size)

            if (allocated(toH2On)) deallocate(toH2On)
            call LogSay(' Time-lag optimization session terminated.')
            end if
            write(*,'(a)')
            write(ulog,'(a)')
        end if
    end if

    !***************************************************************************
    !***************************************************************************
    !********************** PLANAR FIT IF REQUESTED ****************************
    !***************************************************************************
    !***************************************************************************
    if (index(Meth%rot(1:len_trim(Meth%rot)), 'planar_fit') /= 0 .and. &
        (.not. AssessmentOnly .or. RPsetup%pf_assessment_only)) then
        if (.not. RPsetup%pf_onthefly) then
            call ReadPlanarFitFile()
            if (.not. allocated(GoPlanarFit)) &
                allocate(GoPlanarFit(PFSetup%num_sec))
            GoPlanarFit = .true.
            secloop2: do sec = 1, PFSetup%num_sec
                do i = 1, 3
                    do j = 1, 3
                        if (PFMat(i, j, sec) == error) then
                            GoPlanarFit(sec) = .false.
                            cycle secloop2
                        end if
                    end do
                end do
            end do secloop2
        else
            call LogSay(' Performing planar-fit assessment:')

            !> If zero sectors were selected, set to 1 sector by
            !> default and inform
            if (PFSetup%num_sec == 0) then
                call ExceptionHandler(38)
                PFSetup%num_sec = 1
            end if

            !> Allocate variables depending upon number of sectors
            if (.not. allocated(pfNumElem))  &
                allocate(pfNumElem(PFSetup%num_sec))

            if (PFSetup%subperiod) then
                !> Timestamps of start and end of planar fit period
                call DateTimeToDateType(PFSetup%start_date, PFSetup%start_time, &
                    auxStartTimestamp)
                call DateTimeToDateType(PFSetup%end_date, PFSetup%end_time, &
                    auxEndTimestamp)

                !> In RawTimeSeries, detect indexes of first and last files
                !> relevant to planar fit
                call tsExtractSubperiodIndexes(RawTimeSeries, &
                    size(RawTimeSeries), auxStartTimestamp, auxEndTimestamp, &
                    pfStartTimestampIndx, pfEndTimestampIndx)
                pfEndTimestampIndx = pfEndTimestampIndx + 1

                if (pfStartTimestampIndx == nint(error) &
                    .or. pfEndTimestampIndx == nint(error)) &
                    call ExceptionHandler(48)
            else
                pfStartTimestampIndx = 1
                pfEndTimestampIndx = size(RawTimeSeries)
            end if

            !> Count maximum number of periods for planar fit
            write(TmpString1, '(i7)') &
                pfEndTimestampIndx - pfStartTimestampIndx
            write(*, '(a)') '  Maximum number of &
                &flux averaging periods available for planar-fit: ' &
                // trim(adjustl(TmpString1))
            write(ulog, '(a)') '  Maximum number of &
                &flux averaging periods available for planar-fit: ' &
                // trim(adjustl(TmpString1))

            !> Allocate variables that depend upon maximum number of
            !> periods for planar fit
            allocate(pfWind(pfEndTimestampIndx - pfStartTimestampIndx, 3))

            !> Loop on selected files and calculate relevant statistics
            pcount = pfStartTimestampIndx - 1
            pfn = 0
            LatestRawFileIndx = 1
            month = 0
            day   = 0
            DynamicMetadata = ErrDynamicMetadata
            LastMetadataTimestamp = DateType(0, 0, 0, 0, 0)

            !> Same split as the time-lag pre-pass above, and simpler: a
            !> planar-fit period contributes three numbers - the mean wind -
            !> and nothing it computes depends on the period before it, so
            !> there is no warm-up to read and the concatenation is exact.
            call PlanPrepassBatches(pfEndTimestampIndx - pfStartTimestampIndx, &
                .true., pfWorkers)

            if (BatchIndex > 0) then
                pfStartTimestampIndx = BatchSliceStart
                pfEndTimestampIndx = BatchSliceEnd
                pcount = pfStartTimestampIndx - 1
            end if

            pfParallel = pfWorkers > 1
            if (pfParallel) then
                call StartPrepassBatches('pf', pfStartTimestampIndx, &
                    pfEndTimestampIndx, pfWorkers)
                call PrepassSlice(pfStartTimestampIndx, pfEndTimestampIndx, &
                    pfWorkers, 1, sliceStart, sliceEnd)
                pfStartTimestampIndx = sliceStart
                pfEndTimestampIndx = sliceEnd
                pcount = pfStartTimestampIndx - 1
            end if

            pf_periods_loop: do
                pcount = pcount + 1

                !> If embedded metadata are to be used,
                !> reinitialize column information to null
                if (EddyFlowProj%use_extmd_file) then
                    Col = BypassCol
                else
                    Col = NullCol

                end if

                !> Normal exit instruction: either the last period
                !> was dealt with, or raw files are finished
                if (LatestRawFileIndx > NumRawFiles &
                    .or. pcount >= pfEndTimestampIndx) exit pf_periods_loop

                !> Define initial/final timestamps of
                !> current period, say [8:00 - 8:30)
                tsStart = RawTimeSeries(pcount)
                tsEnd   = RawTimeSeries(pcount + 1)

                !> Search file containing data starting from the time closest to
                !> tsStart. Searches only from most current file
                !> onward, to avoid wasting time
                call FirstFileOfCurrentPeriod(tsStart, tsEnd, &
                    RawFileList, NumRawFiles, LatestRawFileIndx, &
                    NextRawFileIndx, skip_period)

                !> Averaging period advancement
                if (day /= 0) then
                    if (EddyFlowProj%caller == 'console') then
                        call LogSayNoAdv('#')
                    else
                        call DisplayProgress('avrg_interval', &
                            '   another small step to the planar-fit: ', &
                                tsStart, 'yes')
                    end if
                end if

                !> Daily advancement
                if (day /= tsStart%day &
                    .or. month /= tsStart%month) then
                    month = tsStart%month
                    day   = tsStart%day
                    if (EddyFlowProj%caller == 'console') then
                        write(*, '(a)')
                        write(ulog, '(a)')
                        call DisplayProgress('daily', &
                            '  Importing wind data for ', tsStart, 'no')
                    else
                        call DisplayProgress('daily', &
                            '  Importing wind data for ', tsStart, 'yes')
                    end if
                end if

                if (skip_period) cycle pf_periods_loop

                !> Import dataset for current period. If using embedded biomet,
                !> also read biomet data
                !> On entrance, NextRawFileIndx contains the index of the file
                !> to start the current period with.
                !> On exit, LatestRawFileIndx contains the index of
                !> the latest file used
                call ImportCurrentPeriod(tsStart, tsEnd, &
                    RawFileList, NumRawFiles, NextRawFileIndx, BypassCol,  &
                    MaxNumFileRecords, MetaIsNeeded, &
                    .false., .false., &
                    Raw, size(Raw, 1), size(Raw, 2), PeriodRecords, &
                    EmbBiometDataExist, skip_period, LatestRawFileIndx, Col, &
                    .false.)
                if (skip_period) cycle pf_periods_loop

                !> Period skip control with message
                MissingRecords = dfloat(MaxPeriodNumRecords - PeriodRecords) &
                    / dfloat(MaxPeriodNumRecords) * 100d0
                if (PeriodRecords > 0 .and. MissingRecords > RPsetup%max_lack) &
                    cycle pf_periods_loop

                !> Filter raw data for user-defined flags
                if (RPsetup%filter_by_raw_flags) &
                    call FilterDatasetForFlags(Col, &
                        Raw, size(Raw, 1), size(Raw, 2))

                !***************************************************************
                !**** RAW FILE IMPORT FINISHES HERE. ***************************
                !**** NOW STARTS DATASET DEFINITION  ***************************
                !***************************************************************

                !> Allocate arrays for actual data processing
                if (.not. allocated(E2Set)) &
                    allocate(E2Set(PeriodRecords, E2NumVar))
                if (.not. allocated(DiagSet)) &
                    allocate(DiagSet(PeriodRecords, MaxNumDiag))

                !> Define EddyFlow set of variables for the following processing
                call DefineE2Set(Col, Raw,   size(Raw, 1),     Size(Raw, 2), &
                                    E2Set,   size(E2Set, 1),   Size(E2Set, 2), &
                                    DiagSet, size(DiagSet, 1), Size(DiagSet, 2))
                !if (allocated(DiagSet))  deallocate(DiagSet)

                !> Clean up E2Set, eliminating values that are clearly un-physical
                call CleanUpE2Set(E2Set, size(E2Set, 1), size(E2Set, 2))

                !> Define as not present, variables for which
                !> too many values are outranged
                call EliminateCorruptedVariables(E2Set, &
                    size(E2Set, 1), size(E2Set, 2), skip_period, .false.)

                !> If either u, v or w have been eliminated,
                !> stops processing this period
                if (skip_period) then
                    if(allocated(E2Set)) deallocate(E2Set)
                    if(allocated(DiagSet)) deallocate(DiagSet)
                    cycle pf_periods_loop
                end if

                !> Update metadata if dynamic metadata are to be used
                if (EddyFlowProj%use_dynmd_file) &
                    call RetrieveDynamicMetadata(tsEnd, &
                        E2Col, size(E2Col))

                !> Override users choices if needed
                call OverrideSettings()

                !***************************************************************
                !**** DATASET DEFINITION FINISHES HERE. ************************
                !**** NOW STARTS RAW DATA REDUCTION     ************************
                !***************************************************************
                !> Filter only for sonic diagnostics (IRGA is irrelevant in
                !> planar fit)
                if (NumDiag > 0) then
                    call FilterDatasetForDiagnostics(E2Set, size(E2Set, 1), &
                        size(E2Set, 2), DiagSet, &
                        size(DiagSet, 1), size(DiagSet, 2), &
                        DiagAnemometer, .false.)
                end if
                if(allocated(DiagSet)) deallocate(DiagSet)

                !> Adjust coordinate systems if the case
                call AdjustSonicCoordinates(E2Set, &
                    size(E2Set, 1), size(E2Set, 2))

                !> Filter for wind direction if requested
                if (RPSetup%apply_wdf) &
                    call FilterDatasetForWindDirection(E2Set, size(E2Set, 1), size(E2Set, 2))

                !> Calculate basic stats
                call BasicStats(E2Set, size(E2Set, 1), size(E2Set, 2), &
                    1, .false.)
                Stats1 = Stats

                !> Calculate raw screening flags and despike data if requeste
                auxTest = TestType(.true., .false., .false., .true., .false., &
                    .false., .false., .false., .false., .false., .false.)
                call StatisticalScreening(E2Set, &
                    size(E2Set, 1), size(E2Set, 2), auxTest, .false.)

                !> Define as not present, variables for which too
                !> many values are outranged
                call EliminateCorruptedVariables(E2Set, &
                    size(E2Set, 1), size(E2Set, 2), skip_period, .false.)

                !> If either u, v or w have been eliminated,
                !> stops processing this period
                if (skip_period) then
                    if(allocated(E2Set)) deallocate(E2Set)
                    cycle pf_periods_loop
                end if

                !> Calculate basic stats
                call BasicStats(E2Set, size(E2Set, 1), size(E2Set, 2), &
                    2, .false.)
                Stats2 = Stats
                Stats3 = Stats

                !> Angle-of-attack calibration
                call AoaCalibration(E2Set, size(E2Set, 1), size(E2Set, 2))

                !> Gill WindMaster w-boost
                if (RPsetup%calib_wboost) &
                    call ApplyGillWmWBoost(E2Set, size(E2Set, 1), size(E2Set, 2))

                !> Metek head correction, the same one the flux loop applies.
                !> Without it here, a lag or a plane would be worked out from
                !> a wind the fluxes never see. Its companion, the
                !> inclinometer correction, cannot follow: it reads its angles
                !> from the custom columns, and DefineUserSet has not run yet
                !> at this point in the program. Recorded rather than worked
                !> around, because moving DefineUserSet is a layout change and
                !> this correction has no effect on any current dataset.
                call MetekHeadCorrection(E2Set, size(E2Set, 1), size(E2Set, 2))

                !> Calculate basic stats
                call BasicStats(E2Set, size(E2Set, 1), size(E2Set, 2), &
                    4, .false.)
                Stats4 = Stats
                if (allocated(E2Set)) deallocate(E2Set)

                !> Store statistics needed for planar fit calcuations
                pfn = pfn + 1
                pfWind(pfn, u) = Stats4%Mean(u)
                pfWind(pfn, v) = Stats4%Mean(v)
                pfWind(pfn, w) = Stats4%Mean(w)
            end do pf_periods_loop
            write(*, '(a)')
            write(ulog, '(a)')
            call LogSay(' Done.')

            if (pfParallel) then
                call WaitPrepassBatches('pf', pfWorkers)
                call MergePfBatchDumps(pfWorkers, pfWind, size(pfWind, 1), pfn)
            end if

            !> As above: a worker hands back its slice of the wind means and
            !> stops. The sector regressions are the parent's job.
            if (BatchIndex > 0) then
                call WritePfBatchDump(pfWind, size(pfWind, 1), pfn)
                call LogSay(' Planar-fit pre-pass slice finished.')
                stop ''
            end if

            !*******************************************************************
            !**** RAW DATA REDUCTION FINISHES HERE.  ***************************
            !**** NOW STARTS PLANAR FIT CALCULATIONS ***************************
            !*******************************************************************

            !> Allocate sector-wise wind array
            if (.not. allocated(pfWindBySect)) &
                allocate(pfWindBySect(pfn, 3, PFSetup%num_sec))
            if (.not. allocated(GoPlanarFit)) &
                allocate(GoPlanarFit(PFSetup%num_sec))

            !> Check if wind components are within specified limits
            where (dsqrt(pfWind(1:pfn, u)**2 + pfWind(1:pfn, v)**2) < PFSetup%u_min &
                   .or. dsqrt(pfWind(1:pfn, u)**2 + pfWind(1:pfn, v)**2) > 20d0 &
                   .or. pfWind(1:pfn, w) > PFSetup%w_max)
                pfWind(1:pfn, u) = error
                pfWind(1:pfn, v) = error
                pfWind(1:pfn, w) = error
            end where

            !> Counts and log excluded stats lines
            err_cnt1 = 0
            do i = 1, pfn
                if(pfWind(i, u) == error .or. pfWind(i, v) == error &
                    .or. pfWind(i, w) == error) err_cnt1 = err_cnt1 + 1
            end do
            !if (err_cnt1 /= 0) then
                !Insert here call to ExceptionHandle and notify
                !that at least 1 wind data was excluded
            !end if

            !> Sort wind data according to wind sector
            !> (from pfWind to pfWindBySect)
            call SortWindBySector(pfWind(1:pfn, u:w), pfn, &
                pfNumElem, pfWindBySect)
            deallocate(pfWind)

            !> Some logging
            write(LogInteger, '(i6)') PFSetup%num_sec
            write(*, '(a)') ' Calculating planar fit rotation matrices for ' &
                                  // trim(adjustl(LogInteger)) // ' sector(s).'
            write(ulog, '(a)') ' Calculating planar fit rotation matrices for ' &
                                  // trim(adjustl(LogInteger)) // ' sector(s).'

            !> Loop over wind sectors
            GoPlanarFit = .true.
            secloop: do sec = 1, PFSetup%num_sec
                write(*, '(a, i2, a)', advance = 'no') '  Sector n.', sec, '..'
                write(ulog, '(a, i2, a)', advance = 'no') '  Sector n.', sec, '..'
                if (PFSetup%wsect_exclude(sec)) then
                    GoPlanarFit(sec) = .false.
                    PFb(:, sec) = error
                    PFMat(:, :, sec) = error
                    call ExceptionHandler(41)
                    cycle secloop
                end if
                if(pfNumElem(sec) < PFSetup%min_per_sec) then
                    GoPlanarFit(sec) = .false.
                    PFb(:, sec) = error
                    PFMat(:, :, sec) = error
                    call ExceptionHandler(33)
                    cycle secloop
                end if

                allocate(pfWind(pfNumElem(sec), 3))
                do i = 1, pfNumElem(sec)
                    pfWind(i, :) = pfWindBySect(i, :, sec)
                end do

                !> Calculate auxiliary variables (see, e.g., Van Dijk et al.2004)
                call PlanarFitAuxParams(pfWind, pfNumElem(sec), Mat, pfVec)
                deallocate(pfWind)

                if (Meth%rot(1:len_trim(Meth%rot)) == 'planar_fit') then
                    !> Invert matrix --> Mat^(-1)
                    call MatrixInversion(Mat, 3, SingMat)

                    !> If singular matrix is found, set results to error
                    if (SingMat) then
                        GoPlanarFit(sec) = .false.
                        PFb(:, sec) = error
                        PFMat(:, :, sec) = error
                        call ExceptionHandler(34)
                        cycle secloop
                    end if

                    !> Calculate plane coefficients: PFb = Mat^(-1) * pfVec
                    PFb(:, sec) = 0d0
                    do i = u, w
                        PFb(:, sec) = PFb(:, sec) + dble(Mat(:, i)) * pfVec(i)
                    end do
                elseif (Meth%rot(1:len_trim(Meth%rot)) &
                        == 'planar_fit_no_bias') then
                    !> Define tensors of 2 elements (out of the 3-elements ones)
                    Mat2d(1,1:2) = Mat(2, 2:3)
                    Mat2d(2,1:2) = Mat(3, 2:3)
                    pfVec2d(1:2) = pfVec(2:3)

                    !> Invert matrix --> Mat^(-1)
                    call MatrixInversion(Mat2d, 2, SingMat)

                    !> If singular matrix is found, set results to error
                    if (SingMat) then
                        GoPlanarFit(sec) = .false.
                        PFb(:, sec) = error
                        PFMat(:, :, sec) = error
                        call ExceptionHandler(34)
                        cycle secloop
                    end if

                    !> Calculate plane coefficients: PFb = Mat^(-1) * pfVec
                    PFb2d(:, sec) = 0d0
                    do i = 1, 2
                        PFb2d(:, sec) = PFb2d(:, sec) &
                            + dble(Mat2d(:, i)) * pfVec2d(i)
                    end do
                    PFb(1, sec) = 0d0
                    PFb(2:3, sec) = PFb2d(1:2, sec)
                end if

                !> Calculate PP (PF rotation matrix, see
                !> Wilczak et al. 2001, BLM)
                call PlanarFitRotationMatrix(sec, PP)

                !> Update sector-wise rotation matrix
                PFMat(:, :, sec) = PP

                call LogSay(' Done.')
            end do secloop

            !> Fix sectors without calculations, using closest
            !> sector in the angular direction defined by user
            if (index(PFSetup%fix,'clockwise') /= 0) &
                call FixPlanarfitSectors(GoPlanarFit, size(GoPlanarFit))

            !> Write planar fit results on output file
            if (PFSetup%num_sec >= 1) &
                call WriteOutPlanarFit(pfNumElem, PFSetup%num_sec)

            if (allocated (pfWindBySect)) deallocate(pfWindBySect)
            if (allocated (pfNumElem)) deallocate(pfNumElem)
            call LogSay(' Planar Fit session terminated.')
            write(*,'(a)')
            write(ulog,'(a)')
        end if
    elseif (.not. AssessmentOnly) then
        if (.not. allocated(GoPlanarFit)) allocate(GoPlanarFit(PFSetup%num_sec))
    end if

    if (AssessmentOnly) then
        if (allocated(bf)) deallocate(bf)
        if (allocated(Raw)) deallocate(Raw)
        call LogSay(' Auxiliary assessment-only session completed.')
        call LogSay(' Normal raw-data processing and flux calculation were not run.')
        stop ''
    end if

    !***************************************************************************
    !***************************************************************************
    !********************** PART COMMON TO NEXT TWO BIG LOOPS ******************
    !***************************************************************************
    !***************************************************************************

    !> Create TimeSeries for actual raw data processing
    if (EddyFlowProj%subperiod) then
        !> If user selected a sub-period, create time series corresponding to
        !> that period and verify that there is any overlap with RawTimeSeries
        call DateTimeToDateType(EddyFlowProj%start_date, &
            EddyFlowProj%start_time, SelectedStartTimestamp)
        call DateTimeToDateType(EddyFlowProj%end_date, &
            EddyFlowProj%end_time, SelectedEndTimestamp)

        NumberOfPeriods = NumOfPeriods(SelectedStartTimestamp, &
            SelectedEndTimestamp, DateStep)
        allocate(MasterTimeSeries(NumberOfPeriods + 1))
        call CreateTimeSeries(SelectedStartTimestamp, SelectedEndTimestamp, &
            DateStep, MasterTimeSeries, size(MasterTimeSeries), .false.)

        !> Verify at least partial overlap
        if (MasterTimeSeries(1) > RawTimeSeries(size(RawTimeSeries)) &
            .or. MasterTimeSeries(size(MasterTimeSeries)) < RawTimeSeries(1)) &
            call ExceptionHandler(46)
    else
        allocate(MasterTimeSeries(size(RawTimeSeries)))
        MasterTimeSeries = RawTimeSeries
    end if

    !> ONLY TEMPORARY: you can now actually eliminate rpStartTimestampIndx and
    !> rpEndTimestampIndx
    rpStartTimestampIndx = 1
    rpEndTimestampIndx = size(MasterTimeSeries)

    !***************************************************************************
    !***************************************************************************
    !******************** DEFINITION OF CALIBRATION EVENTS *********************
    !***************************************************************************
    !***************************************************************************
    if (DriftCorr%method /= 'none' .and. nCalibEvents > 0) then
        call LogSay(' Elaborating IRGA calibration-check history..')

        !> Loop on periods to be processed
        pcount = rpStartTimestampIndx - 1
        LatestRawFileIndx = 1
        latestCleaning = 0
        loggedDate = 'none'
        drift_loop: do
            pcount = pcount + 1

            !> If embedded metadata are to be used,
            !> reinitialize column information to null
            if (EddyFlowProj%use_extmd_file) then
                Col = BypassCol
            else
                Col = NullCol
            end if

            !> Normal exit instruction: either the last period was
            !> dealt with, or raw files are finished
            if (pcount > rpEndTimestampIndx - 1) exit drift_loop

            !> Normal exit instruction: if all cleaning
            !> events have been processed
            if (latestCleaning >= nCalibEvents) exit drift_loop

            !> Define initial/final timestamps of
            !> current period (say, [8:00 to 8:30))
            tsStart = MasterTimeSeries(pcount)
            tsEnd   = MasterTimeSeries(pcount + 1)

            !> If files are finished, keep going until the end of the selected
            !> period
            if (LatestRawFileIndx > NumRawFiles) cycle drift_loop

            !> Search file containing data starting from the time
            !> closest to tsStart
            !> Searches only from most current file onward, to avoid wasting time
            call FirstFileOfCurrentPeriod(tsStart, &
                tsEnd, RawFileList, NumRawFiles, &
                LatestRawFileIndx, NextRawFileIndx, skip_period)
            if (skip_period) cycle drift_loop

            !> Identify if current period is a cleaning event. If not, skip it
            call tsRelaxedMatch(tsStart, &
                Calib(latestCleaning + 1: nCalibEvents)%ts, &
                nCalibEvents - latestCleaning, &
                datetype(0, 0, 0, 3, 0), 'strictly later', clean)

            !> Identify if current period is closest possible to a cleaning event
            call tsRelaxedMatch(tsStart, &
                Calib(latestCleaning + 1: nCalibEvents)%ts, &
                nCalibEvents - latestCleaning, &
                datetype(0, 0, 0, 3, 0), 'strictly before', dirty)

            !> Cycle if file is not relevant to anything
            if (pcount /= rpStartTimestampIndx &
                .and. clean <= 0 .and. dirty <= 0) then
                LatestRawFileIndx = LatestRawFileIndx + 1
                cycle drift_loop
            end if

            !> Log out if the case
            if (pcount /= rpStartTimestampIndx) then
                call DateTypeToDateTime(tsStart, date, time)
                if (date /= loggedDate) then
                    write(*, '(a)') &
                        '  Calibration-check data found on: ' // date(1:10)
                    write(ulog, '(a)') &
                        '  Calibration-check data found on: ' // date(1:10)
                    loggedDate = date
                end if
            end if

            !> Import dataset for current period. If using embedded biomet,
            !> also read biomet data. On entrance, NextRawFileIndx contains
            !> the index of the file to start the current period with
            !> On exit, LatestRawFileIndx contains the index of the
            !> latest file used
            call ImportCurrentPeriod(tsStart, tsEnd, &
                RawFileList, NumRawFiles, NextRawFileIndx, BypassCol, &
                MaxNumFileRecords, MetaIsNeeded, &
                EddyFlowProj%biomet_data == 'embedded', .false., Raw, &
                size(Raw, 1), size(Raw, 2), PeriodRecords, EmbBiometDataExist, &
                skip_period, LatestRawFileIndx, Col, .false.)

            !> Period skip control
            if (skip_period) cycle drift_loop

            !> Period skip control with message
            MissingRecords = dfloat(MaxPeriodNumRecords - PeriodRecords) &
                / dfloat(MaxPeriodNumRecords) * 100d0
            if (PeriodRecords > 0 &
                .and. MissingRecords > RPsetup%max_lack) cycle drift_loop

            !> Calculate reference counts
            call ReferenceCounts(dble(Raw), size(Raw, 1), size(Raw, 2))

            !> The whole gas block, not the first two slots.
            !>
            !> Every other step of the drift chain spans firstGas:lastGas -
            !> DriftRetrieveCalibrationEvents fills %offset and %ref over it,
            !> DriftCorrection consumes it - and only these three assignments
            !> were (co2:h2o). Calib is initialised to error, so a third gas's
            !> %ri and %rf stayed error, the correction masked it out, and a
            !> gas given its own reference and offset columns was silently
            !> never drift-corrected at all.
            !>
            !> Special case of first file in the dataset: used to initialize
            !> drift history assuming cleaned instrument at the beginning
            if (pcount == rpStartTimestampIndx) then
                Calib(0)%ts = MasterTimeSeries(rpStartTimestampIndx)
                call DateTypeToDateTime(Calib(0)%ts, Calib(0)%date, Calib(0)%time)
                Calib(0)%ri(firstGas:lastGas) = refCounts(firstGas:lastGas)
                cycle drift_loop
            end if

            !> Case of cleaning event
            !> Assign relevant ri to current Calib dataset
            if (clean > 0) then
                Calib(latestCleaning + clean)%ri(firstGas:lastGas) = &
                    refCounts(firstGas:lastGas)
                latestCleaning = latestCleaning + clean
            end if
            !> Case of most dirty file (right before next cleaning event)
            !> Calculate and assign relevant quantities to current Calib dataset
            if (dirty > 0) then
                Calib(latestCleaning)%rf(firstGas:lastGas) = &
                    refCounts(firstGas:lastGas)
            end if
        end do drift_loop

        !> Rearrange so that all data needed between
        !> t1 and t2 are stored in Calib(t2)
        tmpCalib = Calib
        do i = 0, nCalibEvents - 1
            !> Calculate number of periods between
            !> each calibration-check event
            Calib(i+1)%numPeriods = &
                NumOfPeriods(Calib(i)%ts, Calib(i+1)%ts, DateStep)

            Calib(i+1)%ri = tmpCalib(i)%ri
            Calib(i+1)%rf = tmpCalib(i)%rf
        end do
        !> For Calib(0) (beginning of dataset), set at clean instrument
        Calib(0)%offset = 0d0
        Calib(0)%ri = error
        Calib(0)%rf = error

!> Only needed and valid for ICOS dataset, where H2O does not
!> start from "clean" but with some offset on day 1, Aug. 3.
!> Artificially set initial ri to the mean value at July 23,
!> when H2O signal was actually "clean", i.e. gives
!> same concentration of LI-7000.
!Calib(1)%ri(h2o) = 34703.78d0

    end if

    !***************************************************************************
    !***************************************************************************
    !***************************** RAW DATA PROCESSING *************************
    !***************************************************************************
    !***************************************************************************

    !> Start loop over all periods contained in the MasterTimeSeries
    NumberOfOkPeriods = 0
    pcount = rpStartTimestampIndx - 1
    LatestRawFileIndx = 1
    bLastFile = 1
    bLastRec = 0
    initialize = .true.
    initializeBiometOut = .true.
    initializeFluxnetOut = .true.
    InitializeStorage = .true.
    InitOutVarPresence = .true.
    DynamicMetadata = ErrDynamicMetadata
    InitGasCalRefCol = GasCalRefCol

    periods_loop: do
        GasCalRefCol = InitGasCalRefCol
        !> Reset CEC state at the start of every period.
        nCecPairs = 0
        do cec_p = 1, MaxNumCecPairs
            call ResetCecDescriptor(CECDescriptor(cec_p))
            call ResetCecFlux(CECFlux(cec_p))
        end do

        !***********************************************************************
        !**** RAW FILE IMPORT **************************************************
        !***********************************************************************

        if (pcount == rpStartTimestampIndx - 1) then
            !> Some log out
            if (EddyFlowProj%run_mode /=  'md_retrieval') then
                call hms_delta_print(' Start raw data processing: ', '')
            else
                call hms_delta_print(' Start metadata retrieving: ', '')
            end if
            call LogSay(' Processing time period:')
            call DateTypeToDateTime(MasterTimeSeries(rpStartTimestampIndx), &
                tmpDate, tmpTime)
            write(*, '(a)') '  Start: ' // tmpDate // ' ' // tmpTime
            write(ulog, '(a)') '  Start: ' // tmpDate // ' ' // tmpTime
            call DateTypeToDateTime(MasterTimeSeries(rpEndTimestampIndx - 1) &
                + DateStep , tmpDate, tmpTime)
            write(*, '(a)') '    End: ' // tmpDate // ' ' // tmpTime
            write(ulog, '(a)') '    End: ' // tmpDate // ' ' // tmpTime
            write(TmpString1, '(i7)') &
                rpEndTimestampIndx - rpStartTimestampIndx
            write(*, '(a)') '  Total number of flux averaging periods: ' &
                // trim(adjustl(TmpString1))
            write(ulog, '(a)') '  Total number of flux averaging periods: ' &
                // trim(adjustl(TmpString1))
            write(*, '(a)')
            write(ulog, '(a)')
        end if
        pcount = pcount + 1

        !> If embedded metadata are to be used, reinitialize
        !> column information to null
        if (EddyFlowProj%use_extmd_file) then
            Col = BypassCol
        else
            Col = NullCol
        end if

        !> Initialize biomet variables to be used in computations
        biomet%val = error

        !> Normal exit instruction: either the last period was dealt with,
        !> or raw files are finished
        if (pcount > rpEndTimestampIndx - 1) exit periods_loop

        !> Define initial/final timestamps of current period (say, 8:00 to 8:29)
        tsStart = MasterTimeSeries(pcount)
        tsEnd   = MasterTimeSeries(pcount + 1)

        !> Associate timestamp of end of the period to current Stats
        call DateTypeToDateTime(tsStart, date, time)
        call DateTypeToDateTime(tsStart, Stats%start_date, Stats%start_time)
        call DateTypeToDateTime(tsEnd, Stats%date, Stats%time)
        call SetPwbPeriodTimestamp(Stats%date, Stats%time)

        !> Some logging
        if (EddyFlowProj%run_mode /= 'md_retrieval') then
            write(*, '(a)')
            write(ulog, '(a)')
            call hms_current_print(' ',': processing new &
                &flux averaging period', .true.)
            write(*, '(a)') ' From: ' &
                // trim(date)   // ' ' // trim(time)
            write(ulog, '(a)') ' From: ' &
                // trim(date)   // ' ' // trim(time)
            write(*, '(a)') '   To: ' &
                // trim(Stats%date) // ' ' // trim(Stats%time)
            write(ulog, '(a)') '   To: ' &
                // trim(Stats%date) // ' ' // trim(Stats%time)
        end if

        !> Define initial part of each output string
        call DateTimeToDOY(Stats%date, Stats%time, int_doy, float_doy)
        write(char_doy, *) float_doy
        call ShrinkString(char_doy)
        suffixOutString =  trim(Stats%date) // ',' // trim(Stats%time) &
                   // ',' // char_doy(1: index(char_doy, '.')+ 4)

        !> Only for external biomet files: retrieve biomet data for current
        !> period and write on output. Even if current period is skipped,
        !> biomet data will be on output.
        if (index(EddyFlowProj%biomet_data, 'ext_') /= 0) then
            call BiometRetrieveExternalData(bFileList, size(bFileList), &
                bLastFile, bLastRec, tsStart, &
                tsEnd, BiometDataFound, .true.)
            call WriteOutBiomet(suffixOutString, .false.)
        end if

        !> If files are finished, keep going until the end of the selected
        !> period
        if (LatestRawFileIndx > NumRawFiles) then
            if (EddyFlowProj%run_mode /= 'md_retrieval') then
                call ExceptionHandler(53)
                if (EddyFlowProj%out_fluxnet) call WriteOutFluxnetOnlyBiomet()
            end if
            call hms_delta_print(PeriodSkipMessage,'')
            cycle periods_loop
        end if

        !> Update metadata if dynamic metadata are to be used
        if (EddyFlowProj%use_dynmd_file) &
        call RetrieveDynamicMetadata(tsEnd, E2Col, size(E2Col))

        MaxNumFileRecords   = nint(Metadata%file_length * 60d0 * Metadata%ac_freq)
        MaxPeriodNumRecords = nint(RPsetup%avrg_len     * 60d0 * Metadata%ac_freq)

        !> Search file containing data starting from the time
        !> closest to tsStart. Searches only from most current
        !> file onward, to avoid wasting time
        call FirstFileOfCurrentPeriod(tsStart, tsEnd, &
            RawFileList, NumRawFiles, LatestRawFileIndx, &
            NextRawFileIndx, skip_period)

        Essentials%fname = trim(adjustl(RawFileList(NextRawFileIndx)%name))

        suffixOutString =  trim(adjustl(RawFileList(NextRawFileIndx)%name)) &
            // ',' // suffixOutString

        !> Exception handling
        if (skip_period) then
            if (EddyFlowProj%run_mode /= 'md_retrieval') then
                call ExceptionHandler(53)
                if (EddyFlowProj%out_fluxnet) call WriteOutFluxnetOnlyBiomet()
            end if
            call hms_delta_print(PeriodSkipMessage,'')
            cycle periods_loop
        end if

        !> Import dataset for current period. If using embedded biomet, also
        !> read biomet data. On entrance, NextRawFileIndx contains index of
        !> file to start the current period with. On exit,
        !> LatestRawFileIndx contains the index of the latest file used
        call ImportCurrentPeriod(tsStart, tsEnd, RawFileList, &
            NumRawFiles, NextRawFileIndx, BypassCol, MaxNumFileRecords, &
            MetaIsNeeded, EddyFlowProj%biomet_data == 'embedded', .true., &
            Raw, size(Raw, 1), size(Raw, 2), PeriodRecords, &
            EmbBiometDataExist, skip_period, LatestRawFileIndx, Col, .true.)

        !> If it's running in metadata retriever mode,
        !> create a dummy dataset 1 minute long
        if (EddyFlowProj%run_mode == 'md_retrieval') then
            PeriodRecords = nint(Metadata%ac_freq * Metadata%file_length * 60d0)
            Raw = 1d0
            NumUserVar = 0
        else
            !> Retrieve biomet data for current period
            if (EddyFlowProj%biomet_data == 'embedded') then

                !> Retrieve biomet data from the already created bSet
                !> Basically, here only convert units and perform average
                !> over the averaging period
                call BiometRetrieveEmbeddedData(EmbBiometDataExist, .true.)

                !> Open biomet output file in case of embedded biomet files
                if(initializeBiometOut .and. nbVars > 0) then
                    call InitBiometOut()
                    initializeBiometOut  = .false.
                end if
                !> Write biomet output
                call WriteOutBiomet(suffixOutString, .true.)
            end if

            if (.not. allocated(UserCol)) &
                allocate(UserCol(NumUserVar))
            call DefineVars(Col, size(Raw, 2), NumUserVar)

            if (initializeFluxnetOut .and. EddyFlowProj%out_fluxnet) then
                call InitFluxnetFile_rp()
                initializeFluxnetOut  = .false.
            end if

            !> Period skip control
            if (skip_period) then
                if (EddyFlowProj%out_fluxnet) call WriteOutFluxnetOnlyBiomet()
                call hms_delta_print(PeriodSkipMessage,'')
                cycle periods_loop
            end if

            !> Number of valid records imported from raw files
            Essentials%n_in = &
                CountRecordsAndValues(dble(Raw), size(Raw, 1), size(Raw, 2))

            !> Some logging
            write(*, '(a, i6)') '  Number of valid records available for this period: ', Essentials%n_in
            write(ulog, '(a, i6)') '  Number of valid records available for this period: ', Essentials%n_in

            !> Period skip control
            MissingRecords = dfloat(MaxPeriodNumRecords - Essentials%n_in) &
                / dfloat(MaxPeriodNumRecords) * 100d0
            if (Essentials%n_in > 0 .and. MissingRecords > RPsetup%max_lack) then
                if (EddyFlowProj%out_fluxnet) call WriteOutFluxnetOnlyBiomet()
                call ExceptionHandler(58)
                call hms_delta_print(PeriodSkipMessage,'')
                cycle periods_loop
            end if

            !> Filter raw data for user-defined flags
            if (RPsetup%filter_by_raw_flags) &
                call FilterDatasetForFlags(Col, Raw, size(Raw, 1), size(Raw, 2))
            Essentials%n_after_custom_flags = &
                CountRecordsAndValues(dble(Raw), size(Raw, 1), size(Raw, 2))

            !> Period skip control
            MissingRecords = dfloat(MaxPeriodNumRecords - Essentials%n_after_custom_flags) &
                / dfloat(MaxPeriodNumRecords) * 100d0
            if (MissingRecords > RPsetup%max_lack) then
                if (EddyFlowProj%out_fluxnet) call WriteOutFluxnetOnlyBiomet()
                call ExceptionHandler(58)
                call hms_delta_print(PeriodSkipMessage,'')
                cycle periods_loop
            end if

            !> If drift correction is to be performed with signal strength
            !> proxy, calculate mean refCounts for current period
            if (DriftCorr%method == 'signal_strength') &
                call ReferenceCounts(dble(Raw), size(Raw, 1), size(Raw, 2))
        end if

        !***********************************************************************
        !**** RAW FILE IMPORT FINISHES HERE. NOW STARTS DATASET DEFINITION *****
        !***********************************************************************

        !> Allocate arrays for actual data processing
        if (.not. allocated(E2Set))    &
            allocate(E2Set(PeriodRecords, E2NumVar))
        if (.not. allocated(E2Primes)) &
            allocate(E2Primes(PeriodRecords, E2NumVar))
        if (.not. allocated(DiagSet))  &
            allocate(DiagSet(PeriodRecords, MaxNumDiag))

        !> Define EddyFlow set of variables for the following processing
        call DefineE2Set(Col, Raw,   size(Raw, 1),     Size(Raw, 2), &
                            E2Set,   size(E2Set, 1),   Size(E2Set, 2), &
                            DiagSet, size(DiagSet, 1), Size(DiagSet, 2))

        !> Some convenient variables
        if (InitOutVarPresence) then
            OutVarPresent(u:E2NumVar) = E2Col(u:E2NumVar)%present
            InitOutVarPresence = .false.
        end if

        !> Define User set of variables, for main statistics
            if (.not. allocated(UserSet)) &
                allocate(UserSet(PeriodRecords, NumUserVar))
            if (.not. allocated(UserCol)) &
                allocate(UserCol(NumUserVar))
            if (.not. allocated(UserPrimes)) &
                allocate(UserPrimes(PeriodRecords, NumUserVar))
            call DefineUserSet(Col, Raw, size(Raw, 1), size(Raw, 2), &
                UserSet, size(UserSet, 1), size(UserSet, 2))

        RowLags = 0
        if (EddyFlowProj%run_mode /= 'md_retrieval') then

            ! !> Update metadata if dynamic metadata are to be used
            if (EddyFlowProj%use_dynmd_file) &
                call RetrieveDynamicMetadata(tsEnd, E2Col, size(E2Col))

            !> Calculate relative separations between the analyzers
            !> and the anemometer used
            call DefineRelativeSeparations()

            !> Override users choices if needed
            call OverrideSettings()

            !> Determine whether it is day or night-time,
            call AssessDayTime(Stats%date, Stats%time)

            !*******************************************************************
            !**** DATASET DEFINITION FINISHES HERE. ****************************
            !**** STARTS RAW DATA REDUCTION         ****************************
            !*******************************************************************
            !> Interpret diagnostics and filter accordingly
            call InterpretLicorDiagnostics(DiagSet, &
                size(DiagSet, 1), size(DiagSet, 2))
            call FilterDatasetForDiagnostics(E2Set, size(E2Set, 1), &
                size(E2Set, 2), DiagSet, &
                size(DiagSet, 1), size(DiagSet, 2), &
                DiagAnemometer, .true.)
            if(allocated(DiagSet)) deallocate(DiagSet)

            !> Adjust coordinate systems if the case
            call AdjustSonicCoordinates(E2Set, size(E2Set, 1), size(E2Set, 2))

            !> Filter for wind direction if requested
            if (RPSetup%apply_wdf) &
                call FilterDatasetForWindDirection(E2Set, size(E2Set, 1), size(E2Set, 2))

            !> Number of valid records after filtering for wind direction
            Essentials%n_after_wdf = &
                CountRecordsAndValues(E2Set, size(E2Set, 1), size(E2Set, 2))
            PeriodActualRecords = Essentials%n_after_wdf
            
            !> Period skip control
            MissingRecords = dfloat(MaxPeriodNumRecords - Essentials%n_after_wdf) &
                / dfloat(MaxPeriodNumRecords) * 100d0
            if (MissingRecords > RPsetup%max_lack) then
                if (EddyFlowProj%out_fluxnet) call WriteOutFluxnetOnlyBiomet()
                if(allocated(E2Set)) deallocate(E2Set)
                if(allocated(E2Primes)) deallocate(E2Primes)
                if(allocated(UserSet)) deallocate(UserSet)
                if(allocated(UserPrimes)) deallocate(UserPrimes)
                call ExceptionHandler(58)
                call LogSayList('')
                call hms_delta_print(PeriodSkipMessage,'')
                cycle periods_loop
            end if

            !> Generate cell temperature dataset if the case, using either
            !> (1) native cell temperature, (2) weighted average of ti1 and ti2,
            !> (3) either ti1 or ti2 depending on availability
            call GenerateTcell(E2Set, size(E2Set, 1), size(E2Set, 2))

            !> Filter Tcell to simulate slower response temperature measurement
            !> for conversion to mixing ratio
            !if (RPsetup%tcell_filter_tconst /= 0) &
            !call CRA(E2Set, size(E2Set, 1), size(E2Set, 2), Metadata%ac_freq, &
            !    RPsetup%tcell_filter_tconst, tc)
        end if

        !> Now that variables have been properly assigned, can initialize
        !> main output files. This is done also if run is in
        !> metadata retriever mode
        if(initialize) then
            call InitOutFiles_rp()
            initialize = .false.
        end if

        !> Output first level of stats
        if (EddyFlowProj%run_mode /= 'md_retrieval') then
            pwb_raw_detection_done = .false.

            !> ===== 1. RAW DATA AS READ FROM FILES ============================
            !> Output raw dataset first level
            if (RPsetup%out_raw(1)) call OutRawData(Stats%date, Stats%time, &
                E2Set, size(E2Set, 1), size(E2Set, 2), 1)
            !> Calculate basic stats and output them as requested
            call BasicStats(E2Set, size(E2Set, 1), size(E2Set, 2), 1, .true.)
            Stats1 = Stats
            if (RPsetup%out_st(1)) &
                call WriteOutStats(ust1, Stats1, suffixOutString, PeriodRecords)
            if (NumUserVar > 0) then
                call UserBasicStats(UserSet, &
                    size(UserSet, 1), size(UserSet, 2), 1)
                if (RPsetup%out_st(1)) &
                    call WriteOutUserStats(u_user_st1, suffixOutString, &
                        PeriodRecords, AddUserStatsHeader)
                    AddUserStatsHeader = .false.
            end if

            !> Based on mean value, if sonic (or fast) temperature
            !> is out-ranged, search alternative one.
            if (Stats1%Mean(ts) < 220d0 .or. Stats1%Mean(ts) > 340d0) &
                call ReplaceSonicTemperature(E2Set, size(E2Set, 1), &
                    size(E2Set, 2), UserSet, size(UserSet, 1), size(UserSet, 2))

            !> ===== 2. STATISTICAL SCREENING ==================================
            !> Calculate raw screening flags and despike data if requested
            call StatisticalScreening(E2Set, &
                size(E2Set, 1), size(E2Set, 2), Test, .true.)
            if (NumUserVar > 0) call DespikeUserSet(UserSet, &
                size(UserSet, 1), size(UserSet, 2))

            !> Define as not present, variables for which
            !> too many values are out-ranged
            call EliminateCorruptedVariables(E2Set, &
                size(E2Set, 1), size(E2Set, 2), skip_period, .true.)

            !> If either u, v or w have been eliminated,
            !> stops processing this period
                if (skip_period) then
                if (EddyFlowProj%out_fluxnet) call WriteOutFluxnetOnlyBiomet()
                if(allocated(E2Set)) deallocate(E2Set)
                if(allocated(E2Primes)) deallocate(E2Primes)
                if(allocated(UserSet)) deallocate(UserSet)
                if(allocated(UserPrimes)) deallocate(UserPrimes)
                call ExceptionHandler(59)
                call LogSayList('')
                call hms_delta_print(PeriodSkipMessage,'')
                cycle periods_loop
            end if

            !> If got until here, incrase number of ok periods
            NumberOfOkPeriods = NumberOfOkPeriods + 1
            
            !> Count values available for each variable and value pairs 
            !> available for each main w-covariance
            !>> 
            Essentials%n = ierror
            Essentials%n_wcov = ierror
            !> Wind data
            Essentials%n(w) = &
                CountRecordsAndValues(E2Set, size(E2Set, 1), size(E2Set, 2), w)
            Essentials%n_wcov(u) = &
                CountRecordsAndValues(E2Set, size(E2Set, 1), size(E2Set, 2), w, u)
            !> Gas data
            do j = ts, lastGas
                if (E2Col(j)%present) then
                    Essentials%n(j) = &
                        CountRecordsAndValues(E2Set, size(E2Set, 1), size(E2Set, 2), j)
                    Essentials%n_wcov(j) = &
                        CountRecordsAndValues(E2Set, size(E2Set, 1), size(E2Set, 2), w, j)
                end if
            end do

            !> If a 4th gas calibration has to be done (using a 'cal-ref'
            !> column from UserCol) does so. Note that so far the calibration
            !> procedure is fully customized on the needs of a
            !> specific O3 analyzer
            call CalibrateGases(E2Set, size(E2Set, 1), size(E2Set, 2))

            !> Output raw dataset second level
            if (RPsetup%out_raw(2)) call OutRawData(Stats%date, Stats%time, &
                E2Set, size(E2Set, 1), size(E2Set, 2), 2)
            !> Calculate basic stats and output them as requested
            call BasicStats(E2Set, size(E2Set, 1), size(E2Set, 2), 2, .true.)
            Stats2 = Stats
            if (RPsetup%out_st(2)) &
                call WriteOutStats(ust2, Stats2, suffixOutString, PeriodRecords)
            if (NumUserVar > 0) then
                call UserBasicStats(UserSet, &
                    size(UserSet, 1), size(UserSet, 2), 2)
                if (RPsetup%out_st(2)) &
                    call WriteOutUserStats(u_user_st2, suffixOutString, &
                        PeriodRecords, AddUserStatsHeader)
                    AddUserStatsHeader = .false.
            end if

            !> ===== 3. CROSS-WIND CORRECTION ==================================
            !> Apply raw-level cross wind correction
            !> (after Liu et al. 2001), if requested
            if (RPsetup%calib_cw) then
                call CrossWindCorr(E2Col(u), E2Set, &
                    size(E2Set, 1), size(E2Set, 2), .true.)
            else
                write(*,'(a)') '  Cross-wind correction not requested &
                    &or not applicable'
                write(ulog,'(a)') '  Cross-wind correction not requested &
                    &or not applicable'
            end if

            !> Output raw dataset third level
            if (RPsetup%out_raw(3)) call OutRawData(Stats%date, Stats%time, &
                E2Set, size(E2Set, 1), size(E2Set, 2), 3)
            !> Calculate basic stats and output them as requested
            call BasicStats(E2Set, size(E2Set, 1), size(E2Set, 2), 3, .true.)
            Stats3 = Stats
            if (RPsetup%out_st(3)) &
                call WriteOutStats(ust3, Stats3, suffixOutString, PeriodRecords)
            if (NumUserVar > 0) then
                call UserBasicStats(UserSet, &
                    size(UserSet, 1), size(UserSet, 2), 3)
                if (RPsetup%out_st(3)) &
                    call WriteOutUserStats(u_user_st3, suffixOutString, &
                        PeriodRecords, AddUserStatsHeader)
                    AddUserStatsHeader = .false.
            end if

            !> ===== 4. ANGLE OF ATTACK CORRECTION =============================
            !> Angle-of-attack calibration
            call AoaCalibration(E2Set, size(E2Set, 1), size(E2Set, 2))

            !> Gill WindMaster w-boost
            if (RPsetup%calib_wboost) &
                call ApplyGillWmWBoost(E2Set, size(E2Set, 1), size(E2Set, 2))

            !> Output raw dataset forth level
            if (RPsetup%out_raw(4)) call OutRawData(Stats%date, Stats%time, &
                E2Set, size(E2Set, 1), size(E2Set, 2), 4)
            !> Calculate basic stats and output them as requested
            call BasicStats(E2Set, size(E2Set, 1), size(E2Set, 2), 4, .true.)
            Stats4 = Stats
            if (RPsetup%out_st(4)) &
                call WriteOutStats(ust4, Stats4, suffixOutString, PeriodRecords)
            if (NumUserVar > 0) then
                call UserBasicStats(UserSet, &
                    size(UserSet, 1), size(UserSet, 2), 4)
                if (RPsetup%out_st(4)) &
                    call WriteOutUserStats(u_user_st4, suffixOutString, &
                        PeriodRecords, AddUserStatsHeader)
                    AddUserStatsHeader = .false.
            end if

            !> ===== 4.1 CORRECTION OF CALIBRATION DRIFTS ======================
            if (DriftCorr%method /= 'none' .and. nCalibEvents /= 0) &
                call DriftCorrection(E2Set, size(E2Set, 1), size(E2Set, 2), &
                    E2Col, size(E2Col), nCalibEvents, tsStart)

            !> ===== 4.2 SONIC HARDWARE CORRECTIONS ===========================
            !> Both act on the raw wind in the sonic's own frame, before any
            !> rotation removes the mean tilt - the head correction because
            !> flow distortion is a property of the direction the wind came
            !> from relative to the instrument, and the inclinometer because
            !> it is putting the instrument's own frame right in the first
            !> place. The head correction goes first: it corrects what the
            !> transducers measured, and the inclinometer then says where
            !> those transducers were pointing.
            call MetekHeadCorrection(E2Set, size(E2Set, 1), size(E2Set, 2))
            if (NumUserVar > 0 .and. allocated(UserSet)) then
                call InclinometerTilt(E2Set, size(E2Set, 1), size(E2Set, 2), &
                    UserSet, size(UserSet, 1), size(UserSet, 2))
            else if (RPSetup%tilt_sensor_meth /= 'none') then
                call LogSay('  Inclinometer tilt correction asked for, but &
                    &the project describes no extra columns to read the &
                    &angles from - skipped.')
            end if

            !> ===== 5. TILT CORRECTION ========================================
            !> Apply rotations for tilt correction, if requested.
            !> NOTE: rotation is applied BEFORE the WPL mixing-ratio conversion so
            !> that pre-WPL PWB time-lag detection sees despiked + rotated
            !> concentrations, matching the Python/RFlux reference. TiltCorrection
            !> only touches wind (u,v,w) and PointByPointToMixingRatio only touches
            !> gas columns, so the two operations commute and the final E2Set is
            !> identical regardless of their order.
            call TiltCorrection(Meth%rot, GoPlanarFit, E2Set, &
                size(E2Set, 1), size(E2Set, 2), PFSetup%num_sec, &
                Essentials%yaw, Essentials%pitch, Essentials%roll, .true.)

            !> PWB time-lag detection on despiked + rotated, pre-WPL data.
            !> Lags are stored now and applied later at the normal timelag
            !> stage. Detecting before PointByPointToMixingRatio matters
            !> because that conversion runs before time-lag compensation: after
            !> it, cell temperature and water sit in the gas series at the
            !> wrong relative lag, and the gas series is what is being
            !> cross-correlated here.
            if (Meth%tlag == 'pwb') then
                call RetrieveSensorParams()
                call SetTimelags()
                pwb_detect_only_mode = .true.
                call TimeLagHandle('pwb', E2Set, size(E2Set, 1), size(E2Set, 2), &
                    pwb_raw_ActTLag, pwb_raw_TLag, pwb_raw_DefTlagUsed, .false.)
                pwb_raw_Result = PWBResult
                pwb_raw_detection_done = .true.
            end if

            !> Remove the spectroscopic effect of water vapour, before the
            !> conversion below and independently of it: what the analyser
            !> reported is biased whatever units it reported in, and the bias
            !> is there whether or not WPL was asked for.
            call SpectroscopicClosedPath(E2Set, &
                size(E2Set, 1), size(E2Set, 2), .true.)

            !> Convert to mixing ratios (if WPL requested, and if the case)
            if (EddyFlowProj%wpl) &
                call PointByPointToMixingRatio(E2Set, &
                    size(E2Set, 1), size(E2Set, 2), .true.)

            !> Output raw dataset fifth level
            if (RPsetup%out_raw(5)) call OutRawData(Stats%date, Stats%time, &
                E2Set, size(E2Set, 1), size(E2Set, 2), 5)
            !> Calculate basic stats and output them as requested
            call BasicStats(E2Set, size(E2Set, 1), size(E2Set, 2), 5, .true.)
            Stats5 = Stats
            if (RPsetup%out_st(5)) &
                call WriteOutStats(ust5, Stats5, suffixOutString, PeriodRecords)
            if (NumUserVar > 0) then
                call UserBasicStats(UserSet, &
                    size(UserSet, 1), size(UserSet, 2), 5)
                if (RPsetup%out_st(5)) &
                    call WriteOutUserStats(u_user_st5, suffixOutString, &
                        PeriodRecords, AddUserStatsHeader)
                    AddUserStatsHeader = .false.
            end if

            !> ===== 6. TIMELAG COMPENSATION  ==================================
            !> If available, for files others than GHG, replace instrument
            !> flow rates provided by user with mean values from raw files
            if (EddyFlowProj%ftype /= 'licor_ghg' &
                .or. EddyFlowProj%use_extmd_file) then
                !> Over every configured gas, not the first four. The measured
                !> flow rate drives tube velocity, Reynolds number and so the
                !> tube-attenuation transfer function; bounded at the fourth
                !> slot, a gas past it silently kept the flow rate declared in
                !> the metadata while its neighbours used the measured one -
                !> so the same analyser's gases were corrected with different
                !> flow rates, and moving a gas between slots changed its
                !> correction factor.
                do i = firstGas, lastGas
                    if (i - firstGas + 1 > &
                        min(EddyFlowProj%gas_num, MaxNumGases)) exit
                    if (NumUserVar > 0) then
                        do j = 1, NumUserVar
                            if (UserCol(j)%var == 'flowrate' &
                                .and. (UserCol(j)%instr_name == E2Col(i)%instr_name &
                                    .or. ((len_trim(UserCol(j)%instr_name) == 0 &
                                        .or. len_trim(E2Col(i)%instr_name) == 0) &
                                        .and. UserCol(j)%instr%model == E2Col(i)%instr%model)) &
                                .and. UserStats%Mean(j) /= 0d0 &
                                .and. UserStats%Mean(j) /= error) then
                                E2Col(i)%instr%tube_f = UserStats%Mean(j)
                                exit
                            end if
                        end do
                    end if
                end do
            end if

            !> Retrieving instruments parameters
            call RetrieveSensorParams()

            !> Defines plausible nominal timelags and timelag ranges
            call SetTimelags()

            !> Calculate and compensate time-lags
            if (TimeLagOptSelected) Meth%tlag = 'maxcov&default'
            call TimeLagHandle(Meth%tlag(1:len_trim(Meth%tlag)), E2Set, &
                size(E2Set, 1), size(E2Set, 2), Essentials%actual_timelag, &
                Essentials%used_timelag, Essentials%def_tlag, .false.)
            if (TimeLagOptSelected) Meth%tlag = 'tlag_opt'
            if (Meth%tlag == 'pwb') then
                PwbTimelagN = PwbTimelagN + 1
                if (.not. allocated(PwbTimelagOpt) .or. PwbTimelagOptSize <= 0 &
                    .or. PwbTimelagN > PwbTimelagOptSize) &
                    error stop 'PWB time-lag optimization dataset is not allocated safely.'
                call AddPwbTimelagSummaryDataset(PwbTimelagOpt, PwbTimelagOptSize, PwbTimelagN)
            end if

            !> ===== 6.1 FILTERING MOLAR DENSITY DATA FOR ABSOLUTE LIMITS TEST  ====================
            if (EddyFlowProj%run_mode /= 'md_retrieval') then
                !> Estimate temperatures, pressures and relevant
                !> air molar volumes
                call AirAndCellParameters()
                if (Test%al .and. RPsetup%filter_al) then
                    !> Apply filter for absolute limits test, if the case
                    FilterWhat = .false.
                    FilterWhat(firstGas:lastGas) = .true.
                    call FilterDatasetForPhysicalThresholds(E2Set, &
                        size(E2Set, 1), size(E2Set, 2), FilterWhat)
                    !> Define as not present, variables for which &
                    !> too many values are out-ranged
                    call EliminateCorruptedVariables(E2Set, &
                        size(E2Set, 1), size(E2Set, 2), skip_period, .true.)
                end if
            end if

            !> Output raw dataset sixth level
            if (RPsetup%out_raw(6)) call OutRawData(Stats%date, Stats%time, &
                E2Set, size(E2Set, 1), size(E2Set, 2), 6)
            !> Calculate basic stats and output them as requested
            call BasicStats(E2Set, size(E2Set, 1), size(E2Set, 2), 6, .true.)
            Stats6 = Stats
            if (RPsetup%out_st(6)) &
                call WriteOutStats(ust6, Stats6, suffixOutString, PeriodRecords)
            if (NumUserVar > 0) then
                call UserBasicStats(UserSet, &
                    size(UserSet, 1), size(UserSet, 2), 6)
                if (RPsetup%out_st(6)) &
                    call WriteOutUserStats(u_user_st6, suffixOutString, &
                        PeriodRecords, AddUserStatsHeader)
                    AddUserStatsHeader = .false.
            end if

            !> Quality check test for stationarity
            call StationarityTest(E2Set, size(E2Set, 1), size(E2Set, 2), StDiff)

            !> Calculate wind speed and maximum wind speed
            call MaxWindSpeed(E2Set, &
                size(E2Set, 1), size(E2Set, 2), Ambient%MWS)
            if (Stats6%mean(u) /= error .and. Stats6%mean(v) /= error &
                .and. Stats6%mean(w) /= error) then
                Ambient%WS = dsqrt(Stats6%mean(u)**2 &
                    + Stats6%mean(v)**2 + Stats6%mean(w)**2)
            else
                Ambient%WS = error
            end if

            !> ===== 6.2 QC tests =============================================
            !> Calculate Kurtosis Index on differenced variables
            call KID(E2Set(:, 1:GHGNumVar), size(E2Set, 1), GHGNumVar)

            !> Calculate Longest Gap Duration
            call LongestGapDuration(E2Set(:, 1:GHGNumVar), size(E2Set, 1), GHGNumVar)

            !> ===== 6.3 CONDITIONAL EDDY COVARIANCE ==========================
            !> Before the run's own detrending, because the partition screens
            !> the raw series on the analyser diagnostics and then detrends what
            !> survives - and E2Set is gone a few lines below. Screening after
            !> the fact would leave the trend fitted through the samples being
            !> rejected.
            !>
            !> Also, necessarily, before FixDatasetForSpectra: that interpolates
            !> E2Primes for the spectra, and a partition built on fabricated
            !> samples would pass a completeness gate it should not.
            nCecPairs = 0
            do cec_p = 1, MaxNumCecPairs
                call ResetCecDescriptor(CECDescriptor(cec_p))
                call ResetCecFlux(CECFlux(cec_p))
            end do
            if (EddyFlowProj%do_cec > 0) then
                call LogSayNoAdv('  Calculating CEC partitioning..')
                call CecPairs(CecPairList, nCecPairs)
                if (nCecPairs > 0) then
                    if (.not. allocated(CecPrimes)) &
                        allocate(CecPrimes(PeriodRecords, MaxNumCecTargets + 1))
                    do cec_p = 1, nCecPairs
                        call BuildCecPrimes(CecPairList(cec_p), E2Set, &
                            size(E2Set, 1), size(E2Set, 2), UserSet, NumUserVar, &
                            RPsetup%Tconst, EddyFlowProj%cec, CecPrimes, cec_ok)
                        if (.not. cec_ok) cycle
                        call ExtractCecDescriptor(CecPairList(cec_p), CecPrimes, &
                            size(CecPrimes, 1), &
                            StDiff%w_gas(CecPairList(cec_p)%carbon_slot), &
                            StDiff%w_gas(CecPairList(cec_p)%water_slot), &
                            CECDescriptor(cec_p), EddyFlowProj%cec)
                    end do
                    if (allocated(CecPrimes)) deallocate(CecPrimes)
                end if
                call LogSay(' Done.')
            end if

            !> ===== 7. DETRENDING =============================================
            !> Calculate fluctuations based on chosen detrending method
            call LogSayNoAdv('  Detrending..')
            call Fluctuations(E2Set, E2Primes, &
                size(E2Set, 1), size(E2Set, 2), RPsetup%Tconst, Stats, E2Col)
            if (allocated(E2Set)) deallocate(E2Set)
            call LogSay(' Done.')
            if (NumUserVar > 0) then
                call LogSayNoAdv('  Detrending user set..')
                call UserFluctuations(UserSet, UserPrimes, &
                    size(UserSet, 1), size(UserSet, 2), &
                    RPsetup%Tconst, UserStats, UserCol)
                call LogSay(' Done.')
            end if

            !> Output raw dataset seventh level
            if (RPsetup%out_raw(7)) &
                call OutRawData(Stats%date, Stats%time, E2Primes, &
                    size(E2Primes, 1), size(E2Primes, 2), 7)
            !> Calculate basic stats and output them as requested
            call BasicStats(E2Primes, &
                size(E2Primes, 1), size(E2Primes, 2), 7, .true.)
            Stats7 = Stats
            if (RPsetup%out_st(7)) &
                call WriteOutStats(ust7, Stats7, suffixOutString, PeriodRecords)
            if (NumUserVar > 0) then
                call UserBasicStats(UserPrimes, &
                    size(UserPrimes, 1), size(UserPrimes, 2), 7)
                if (RPsetup%out_st(7)) &
                    call WriteOutUserStats(u_user_st7, suffixOutString, &
                        PeriodRecords, AddUserStatsHeader)
                    AddUserStatsHeader = .false.
            end if
            if (allocated(UserPrimes)) deallocate(UserPrimes)

            !> ===== 7.1 QC tests =============================================
            !> Fisher's test
            !> ncol must describe the section actually passed, not E2Primes's
            !> full second dimension - that is E2NumVar, half as wide again as
            !> the GHGNumVar columns handed over, so the explicit-shape dummy
            !> inside read past the end of them. Matches the KID call above.
            call Fisher(E2Primes(:, 1:GHGNumVar), size(E2Primes, 1), GHGNumVar)

            !> Cross-correlation R^2 test for repeated values - informational
            !> only (see cross_corr_test.f90's own header), so gated the
            !> same as its sibling raw-signal diagnostics rather than
            !> always spending the two extra CCF passes per variable.
            if (Test%rf) &
                call CrossCorrTest(E2Primes(:, 1:GHGNumVar), size(E2Primes, 1), size(E2Primes, 2))

            !> Calculate Mahrt's random error and Nonstationarity ratio anyway.
            call RU_Mahrt_98(E2Primes, size(E2Primes, 1), size(E2Primes, 2))

            !> If requested, estimate random error
            call RandomUncertaintyHandle(E2Primes, size(E2Primes, 1), size(E2Primes, 2))

            if (allocated(UserSet)) deallocate(UserSet)

            !*******************************************************************
            !**** RAW DATA REDUCTION FINISHES HERE *****************************
            !**** CALCULATE AND OUTPUT CO-SPECTRA  *****************************
            !*******************************************************************
            if (RPsetup%do_spectral_analysis) then
                SpecCol = E2Col

                !> Replace gaps with linear interpolation of neighbouring data
                call FixDatasetForSpectra(E2Primes, &
                    size(E2Primes, 1), size(E2Primes, 2), N2)

                !> Set length of dataset by stripping
                !> trailing/leading error codes
                Nmax = maxval(RowLags)
                Nmin = minval(RowLags)
                max_nsmpl = N2 - (Nmax - Nmin)
                if(RPsetup%power_of_two) then
                    !> Calculate power-of-two closest to number of
                    !> available samples
                    call PowerOfTwo(max_nsmpl, SpecRow)
                else
                    !> use all samples
                    SpecRow = max_nsmpl
                end if

                Nmin = - Nmin
                !> Which row of E2Primes the spectral set starts at. A slower
                !> column's samples are located relative to E2Primes, and this
                !> is what carries that phase across the trim.
                SpecRowOffset = Nmin
                !> Recalculate basic statistics with PeriodRecords=SpecRow
                !> as number of observations
                allocate(SpecSet(SpecRow, E2NumVar))
                SpecSet(1: SpecRow, :) = E2Primes(Nmin + 1: Nmin + SpecRow, :)
                call BasicStats(SpecSet, &
                    size(SpecSet, 1), size(SpecSet, 2), 8, .true.)

                !> Calculate spectra and cospectra and output them all
                !> Pass the whole gas block, not just the four legacy slots:
                !> SpectralAnalysis loops u..GHGNumVar internally, so a
                !> narrower slice reads past the end of the array as soon as a
                !> project describes more than four gases.
                call SpectralAnalysis(Stats%date, Stats%time, bf, &
                    SpecSet(:, u:lastGas), size(SpecSet, 1), lastGas)
                if (allocated(SpecSet)) deallocate(SpecSet)

                !> Reset stats to Stats7, after the parenthesis
                !> of spectral analysis
                Stats = Stats7
            else
                Essentials%degH(:) = error
            end if
        end if
        !> E2Primes deallocation is deferred to after CEC computation below
        if (allocated(UserPrimes)) deallocate(UserPrimes)
        if (allocated(UserSet)) deallocate(UserSet)

        !***********************************************************************
        !**** (CO)SPECTRA CALCULATION FINISHES HERE  ***************************
        !**** NOW STARTS FLUX COMPUTATION/CORRECTION ***************************
        !***********************************************************************
        if (EddyFlowProj%run_mode /= 'md_retrieval') then

            !> Average mole fractions in [umol mol_a-1] and [mmol mol_a-1]
            call MoleFractionsAndMixingRatios()

            !> Calculate parameters for flux computation
            call FluxParams(.true.)

            !> Cleared every period, and for every gas.
            !>
            !> The loop below only assigns to gases on an LI-7700, so without
            !> this a gas that is not on one would read whatever the array
            !> last held - and, with a scalar, so would a period in which the
            !> loop found nothing at all. Only IsLi7700 gases ever consult it,
            !> but leaving stale numbers where a reader might look is how the
            !> next fault gets built.
            Mul7700(:)%A = error
            Mul7700(:)%B = error
            Mul7700(:)%C = error

            !> LI-7700 spectroscopic correction. It applies to whichever gas
            !> the LI-7700 measures, which is a question about the analyser -
            !> asked of slot seven, it both missed a 7700 sitting on any other
            !> record and, on a project whose seventh slot holds something
            !> else, would have scaled the wrong gas. The multipliers depend
            !> only on P, T and water, so they are computed once.
            do j = firstGas, lastGas
                if (E2Col(j)%Instr%model(1:max(1, &
                    len_trim(E2Col(j)%Instr%model) - 2)) /= 'li7700') cycle
                !> Calculate multipliers for LI-7700 spectroscopic correction,
                !> from the water this gas is corrected with.
                !>
                !> Eq. 6.13 is a spectroscopic property of the sample the
                !> analyser sees, so it wants that analyser's humidity. This
                !> read the site's, which on a two-analyser site corrects a
                !> 7700 with a hygrometer it does not share air with - and
                !> through the fallback variant, which names a trace gas when
                !> a project has no water at all.
                !> Whatever the gas names, including the biomet. A biomet RH
                !> is not the analyser's own air either, but it is a
                !> measurement of the air the analyser is sampling from,
                !> which another instrument's cell is not - and it is what
                !> the user asked for by naming it.
                msl = E2Col(j)%moist_ref
                if (msl == biometMoistRef) then
                    chi_moist = Ambient%chi_biomet
                elseif (msl >= firstGas .and. msl <= lastGas) then
                    chi_moist = Stats%chi(msl)
                else
                    cycle
                end if
                if (chi_moist == error) cycle
                call Multipliers7700(Stats%Pr, Ambient%Ta, &
                    chi_moist, &
                    Mul7700(j)%A, Mul7700(j)%B, Mul7700(j)%C)
                !> Modify mole fraction and mixing ratio to account for
                !> key(T,P), Eq. 6.13 of LI-7700 manual
                !> Uses multiplies A, because this is equal to key.
                Stats%chi(j) = Stats%chi(j) * Mul7700(j)%A
                Stats%r(j)   = Stats%r(j)   * Mul7700(j)%A
            end do

            !> Calculate LI-7500 surface heating correction if requested
            call BurbaTerms()

            !> Calculate fluxes at Level 0
            call Fluxes0_rp(.true.)

            !> As of now, still use CO2 analyzer software version as a proxy for
            !> Logger software version. However, the machinery is in place
            !> for using logger version from [Station], simply remove the
            !> following line.
!            if (E2Col(co2)%instr%sw_ver /= errSwVer) then
            !> From the first configured gas's analyser. Was E2Col(co2) -
            !> slot five - which is that analyser only when CO2 is record one.
            Metadata%logger_swver = E2Col(FirstConfiguredGasSlot())%instr%sw_ver
!            elseif (E2Col(h2o)%instr%sw_ver /= errSwVer) then
!                Metadata%logger_swver = E2Col(h2o)%instr%sw_ver
!            end if

            if (.not. EddyFlowProj%fcc_follows) then
                !> Spectral correction and the two flux levels, once or
                !> repeatedly.
                !>
                !> The three depend on each other in a circle: the analytic
                !> cospectrum is evaluated at z/L, z/L comes from the
                !> corrected sensible heat flux, and that flux is what the
                !> spectral correction produces. One pass leaves them
                !> disagreeing - the correction was computed at a stability
                !> the run then went on to revise.
                !>
                !> Iterating closes the circle. Nothing here accumulates:
                !> Fluxes1_rp rebuilds Flux1 from Flux0 and BPCF, and
                !> Fluxes23_rp rebuilds Flux2 and Flux3 from Flux1, so each
                !> pass is a fresh correction of the same raw covariances
                !> rather than a correction of a correction. Fluxes23_rp
                !> already recomputes Ambient%zL, which is what the next pass
                !> reads - the feedback path was there, only the repetition
                !> was missing.
                !>
                !> Off by default and a single pass then, which is exactly
                !> what this block did before.
                iter_dev = error
                corr_passes = 1
                if (EddyFlowProj%corr_iter_meth) &
                    corr_passes = EddyFlowProj%corr_iter_max
                do corr_pass = 1, corr_passes
                    if (corr_pass > 1) prev_gas_flux = Flux3%gas

                    !> Low-pass and high-pass spectral correction factors
                    call BandPassSpectralCorrections(E2Col(u)%Instr%height, &
                        Metadata%d, E2Col(u:GHGNumVar)%present, Ambient%WS, Ambient%Ta, &
                        Ambient%zL, Metadata%ac_freq, RPsetup%avrg_len, &
                        Metadata%logger_swver, Meth%det, &
                        RPsetup%Tconst, corr_pass == 1, E2Col(u:GHGNumVar)%instr, 1)

                    !> Calculate fluxes at Level 1
                    call Fluxes1_rp()

                    !> Calculate fluxes at Level 2 and Level 3
                    call Fluxes23_rp()

                    if (corr_pass > 1) then
                        iter_dev = WorstRelativeChange(prev_gas_flux, Flux3%gas)
                        !> A tolerance of zero never fires, which is EddyUH's
                        !> behaviour: it runs every pass and tests nothing.
                        if (EddyFlowProj%corr_iter_tol > 0d0 &
                            .and. iter_dev /= error &
                            .and. iter_dev < EddyFlowProj%corr_iter_tol) exit
                    end if
                end do
                Essentials%corr_iter_dev = iter_dev

                !> Footprint estimation
                foot_model_used = Meth%foot(1:len_trim(Meth%foot))
                call FootprintHandle(Stats%Cov(w, w), Ambient%us, &
                    Ambient%zL, Ambient%WS, Ambient%L, &
                    E2Col(u)%Instr%height, Metadata%d, Metadata%z0)
            else
                !> Whether any gas needs the FCC-only path. That is a species
                !> question - CO2, water and CH4 have in-situ spectral
                !> corrections and nothing else does - so it is asked of the
                !> records rather than of the fourth slot. A project with COS
                !> on record five and nothing on record four used to take the
                !> else arm and lose every flux for the period.
                has_fcc_only_gas = .false.
                do j = firstGas, lastGas
                    if (.not. OutVarPresent(j)) cycle
                    if (.not. HasInSituSpectralCorrection(j)) &
                        has_fcc_only_gas = .true.
                end do
                if (has_fcc_only_gas) then
                    !> Those gases cannot use in-situ spectral corrections
                    !> (FCC-only): compute with BPCF=1.0, then zero the ones
                    !> that can so FCC corrects those later.
                    call BandPassSpectralCorrections(E2Col(u)%Instr%height, &
                        Metadata%d, E2Col(u:GHGNumVar)%present, Ambient%WS, Ambient%Ta, &
                        Ambient%zL, Metadata%ac_freq, RPsetup%avrg_len, &
                        Metadata%logger_swver, Meth%det, &
                        RPsetup%Tconst, .true., E2Col(u:GHGNumVar)%instr, 1)
                    call Fluxes1_rp()
                    call Fluxes23_rp()
                    do j = firstGas, lastGas
                        if (.not. HasInSituSpectralCorrection(j)) cycle
                        Flux1%gas(j) = error
                        Flux2%gas(j) = error
                        Flux3%gas(j) = error
                    end do
                else
                    Flux1 = errFlux
                    Flux2 = errFlux
                    Flux3 = errFlux
                    BPCF = errBPCF
                end if
                Foot = errFootprint
                foot_model_used = 'none'
            end if

            !> RP applies the CEC descriptors only when it owns the final
            !> corrected totals. Otherwise FCC applies them to FCC's Flux3.
            if (EddyFlowProj%do_cec > 0 .and. .not. EddyFlowProj%fcc_follows) then
                do cec_p = 1, nCecPairs
                    call CecTargetSlots(CecPairList(cec_p), cec_slots, cec_ntarget)
                    cec_totals = error
                    cec_errors = error
                    do cec_k = 1, cec_ntarget
                        if (cec_slots(cec_k) < firstGas &
                            .or. cec_slots(cec_k) > lastGas) cycle
                        cec_totals(cec_k) = Flux3%gas(cec_slots(cec_k))
                        !> From the same slot as the total beside it. Reading
                        !> either from a different one would test a flux
                        !> against another gas's error and say nothing.
                        cec_errors(cec_k) = Essentials%rand_uncer(cec_slots(cec_k))
                    end do
                    call ApplyCecDescriptor(CECDescriptor(cec_p), cec_totals, &
                        cec_errors, EddyFlowProj%cec, CECFlux(cec_p))
                end do
            end if
            if (allocated(E2Primes)) deallocate(E2Primes)

            !> Calculate storage terms
            if(InitializeStorage) then
                Stor%H  = error
                Stor%LE = error
                Stor%of = error
                InitializeStorage = .false.
            else
                call Storage(PrevStats, prevAmbient)
            end if
            if (Test%stor_clean) call StoreStorCache(Stats%date, Stats%time)
            PrevStats = Stats
            prevAmbient = Ambient
            prevBiomet = biomet

            !> Well developed turbulence conditions test,
            !> after Foken et al. (2004, Handbook of Microm.)
            call DevelopedTurbulenceTest(DtDiff)

            !> flagging the file, after Foken et al. (2004, Handbook of Microm.)
            call QualityFlags(Flux2, StDiff, DtDiff, STFlg, DTFlg, QCFlag, .true., .true., Test%rf)

            !> Write details on output files if requested
            if(RPsetup%out_qc_details .and. Meth%qcflag /= 'none') &
                call WriteOutQCDetails(suffixOutString, StDiff, DtDiff, STFlg, DTFlg)

            !> Update values of AGC and RSSI as available
            call SetLicorDiagnostics(NumUserVar)
        end if

        !>Write out full output file (main express output)
        if (EddyFlowProj%out_full) &
            call WriteOutFull(suffixOutString, PeriodRecords, PeriodActualRecords)

        !>Write out full output file (main express output)
        if (EddyFlowProj%out_md) &
            call WriteOutMetadata(suffixOutString)
            if (EddyFlowProj%out_fluxnet) &
                call WriteOutFluxnet(StDiff, DtDiff, STFlg, DTFlg)

        if (EddyFlowProj%run_mode /= 'md_retrieval') then
            call hms_delta_print('  Flux averaging period processing time: ','')
        else
            call hms_delta_print('  Metadata retrieving time: ','')
        end if
        write(*, *)
        write(ulog, *)

        if (allocated(UserCol))  deallocate(UserCol)
        if (allocated(E2Set))    deallocate(E2Set)
        if (allocated(E2Primes)) deallocate(E2Primes)
        if (allocated(DiagSet))  deallocate(DiagSet)
        if (allocated(UserSet))  deallocate(UserSet)
    end do periods_loop
    if (Test%stor_clean .and. StorCacheN > 0) call PostProcessStorClean()
    if (allocated(bf)) deallocate(bf)
    if (Meth%tlag == 'pwb') call ReportPwbDiagnostics()
    if (Meth%tlag == 'pwb' .and. PwbCacheDirty) call WritePwbTimelagCache()
    if (Meth%tlag == 'pwb' .and. PwbTimelagN > 0) then
        if (.not. allocated(PwbTimelagOpt) .or. PwbTimelagOptSize <= 0 &
            .or. PwbTimelagN > PwbTimelagOptSize) &
            error stop 'PWB time-lag optimization dataset is not allocated safely.'
        allocate(toSet(PwbTimelagN))
        call FixTimelagOptDataset(PwbTimelagOpt, PwbTimelagOptSize, &
            toSet, size(toSet), tlagn, size(tlagn))
        allocate(toH2On(TOSetup%h2o_nclass))
        call OptimizeTimelags(toSet, size(toSet), tlagn, E2NumVar, toH2On, &
            TOSetup%h2o_nclass, TOSetup%h2o_class_size)
        call ResolvePwbAggregateSummary(tlagn)
        PwbAggregateSummary = .true.
        call WriteOutTimelagOptimization(tlagn, E2NumVar, toH2On, &
            TOSetup%h2o_nclass, TOSetup%h2o_class_size)
        PwbAggregateSummary = .false.
        deallocate(toSet)
        deallocate(toH2On)
        deallocate(PwbTimelagOpt)
    end if

    !***************************************************************************
    !**** FLUX COMPUTATION FINISHES HERE.                      *****************
    !**** NOW STARTS DATASET CREATION AND OUTPUT FILE HANDLING *****************
    !***************************************************************************
    close(ust1)
    close(ust2)
    close(ust3)
    close(ust4)
    close(ust5)
    close(ust6)
    close(ust7)
    close(u_user_st1)
    close(u_user_st2)
    close(u_user_st3)
    close(u_user_st4)
    close(u_user_st5)
    close(u_user_st6)
    close(u_user_st7)
    close(umd)
    close(uflx)
    close(ufnet_e)
    close(ufnet_b)
    close(uaflx)
    close(uex)
    close(uflxnt)
    close(ubiomet)
    close(uqc)


    !> If no averaging period was performed, return message and cancel tmp files
    if (NumberOfOkPeriods == 0 .and. EddyFlowProj%run_mode /= 'md_retrieval') then
        !> Delete files in output folder
        del_status = system(trim(comm_del) // ' "' // trim(adjustl(Dir%main_out)) &
            // '*' // Timestamp_FilePadding //'*"'  // comm_err_redirect)

        !> Alerting and closing run
        write(*,'(a)')
        write(ulog,'(a)')
        call ExceptionHandler(35)
    end if

    !> Creating datasets from output files
    write(*, '(a)')
    write(ulog, '(a)')
    write(*, '(a)') ' Raw data processing terminated. &
        &Creating continuous datasets if necessary..'
    write(ulog, '(a)') ' Raw data processing terminated. &
        &Creating continuous datasets if necessary..'

    if (make_dataset_common) then
        call CreateDatasetsCommon(MasterTimeSeries, size(MasterTimeSeries), &
            rpStartTimestampIndx, rpEndTimestampIndx, 'RP')
    else
        call RenameTmpFilesCommon()
    end if
    if (make_dataset_rp) then
        call CreateDatasetsRP(MasterTimeSeries, size(MasterTimeSeries), &
            rpStartTimestampIndx, rpEndTimestampIndx)
    else
        call RenameTmpFilesRP()
    end if
    call LogSay(' Done.')

    !> Edit .eddypro file updating path to ex_file
    call ForceSlash(FLUXNET_Path, .false.)
    call EditIniFile(trim(PrjPath), 'ex_file', &
        trim(FLUXNET_Path(1:index(FLUXNET_Path, '.tmp')-1)))

    if (EddyFlowProj%run_env /= 'embedded') &
        write(*, '(a)') ' FLUXNET file path: ' &
            // trim(FLUXNET_Path(1:index(FLUXNET_Path, '.tmp')-1))
        write(ulog, '(a)') ' FLUXNET file path: ' &
            // trim(FLUXNET_Path(1:index(FLUXNET_Path, '.tmp')-1))

    !> Copy ".eddypro" file into output folder
    if (.not. EddyFlowProj%fcc_follows) then
        call CopyFile(trim(adjustl(PrjPath)), &
        trim(adjustl(Dir%main_out)) // 'processing' &
        // Timestamp_FilePadding // '.eddyflow')
    end if

    !> Whatever the last prefetch left behind. Desktop mode removes the whole
    !> temporary directory below and would take it with it; embedded mode
    !> keeps that directory between runs and would not.
    call GhgPrefetchCleanup()

    !> Delete tmp folder if running in embedded mode
    if(EddyFlowProj%run_env == 'desktop') &
        del_status = system(trim(comm_rmdir) // ' "' &
        // trim(adjustl(TmpDir)) // '"')

    if (.not. EddyFlowProj%fcc_follows) then
        call LogSay('')
        call LogSay(' ****************************************************')
        call LogSay(' Program EddyFlow executed gracefully.')
        call LogSay(' Check results in the selected output directory.     ')
        call LogSay(' ****************************************************')
    end if
    stop ''
end program EddyFlowRP
