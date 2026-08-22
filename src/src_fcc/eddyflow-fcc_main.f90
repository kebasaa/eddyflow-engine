!***************************************************************************
! eddyflow-fcc_main.f90
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
Program EddyFlowFCC
    use m_fx_global_var
    use m_cec
    implicit none

    integer, external :: CreateDir
    integer :: i
    integer :: gas
    !> Iterative correction, twinned with eddyflow-rp_main.f90.
    integer :: corr_pass
    integer :: corr_passes
    real(kind = dbl) :: iter_dev
    real(kind = dbl) :: iter_zL
    real(kind = dbl) :: prev_gas_flux(GHGNumVar)
    real(kind = dbl), external :: WorstRelativeChange
    integer :: month
    integer :: nbins
    integer :: open_status
    integer :: nfit(GHGNumVar, 2)
    integer :: day
    integer :: STFlg(GHGNumVar)
    integer :: DTFlg(GHGNumVar)
    integer :: int_doy
    integer :: NumValidExRecords
    integer :: FirstValidRecord
    integer :: NumExRecords
    integer :: NumberOfPeriods
    integer :: fxStartTimestampIndx
    integer :: fxEndTimestampIndx
    integer :: saStartTimestampIndx
    integer :: saEndTimestampIndx
    integer :: NumFullFiles
    integer :: NumFullFilesNoRecurse
    integer :: NumBinnedFiles
    integer :: NumBinnedFilesNoRecurse
    integer :: fcount
    integer :: nrow_full
    integer :: del_status

    character(10) :: sDate, eDate
    character(5) :: sTime, eTime

    logical :: skip
    logical :: InitializeOuputFiles
    logical :: skip_spectra
    logical :: skip_cospectra
    logical :: ValidRecord
    logical :: EndOfFileReached
    logical :: exEndReached

    !> Derived type variables
    type(DateType) :: exStartTimestamp
    type(DateType) :: exEndTimestamp
    type(DateType) :: saStartTimestamp
    type(DateType) :: saEndTimestamp
    type(DateType) :: binStartTimestamp
    type(DateType) :: binEndTimestamp
    type(DateType) :: SelectedStartTimestamp
    type(DateType) :: SelectedEndTimestamp
    type(Datetype) :: CurrentTimestamp
    type(SpectraSetType) :: BinSpec(MaxNumBins)
    type(SpectraSetType) :: BinCosp(MaxNumBins)
    type(SpectraSetType) :: BinCospForStable(MaxNumBins)
    type(SpectraSetType) :: BinCospForUnstable(MaxNumBins)
    type(InstrumentType) :: AuxInstrument(GHGNumVar)
    type(QCType) :: StDiff
    type(QCType) :: DtDiff
    type(ExType) :: lEx
    integer :: cec_p
    integer :: cec_k
    integer :: cec_slot
    real(kind = dbl) :: cec_totals(MaxNumCecTargets)
    real(kind = dbl) :: cec_errors(MaxNumCecTargets)

    !> Allocatable variabled
    type(DateType), allocatable :: exTimeSeries(:)
    type(DateType), allocatable :: MasterTimeSeries(:)
    type(FitSpectraType), allocatable :: FitUnstable(:)
    type(FitSpectraType), allocatable :: FitStable(:)

    !> External functions
    integer, external :: NumOfPeriods
    integer, external :: NumberOfFilesInSubperiod
    real(kind = dbl), external :: func

    include '../src_common/interfaces.inc'

    !*******************************************************************************
    !*******************************************************************************

    !> Connect the log before anything is said - see the RP main.
    call LogStart()

    app = fcc_app

    !> Initialize environment
    write(*, '(a)')
    write(ulog, '(a)')
    call InitEnv()

    !> By detault, create FLUXNET output
    EddyFlowProj%out_fluxnet = .true.

    call LogSay('Starting flux computation and correction session..')
    write(*, '(a)')
    write(ulog, '(a)')

    !> Read ".eddypro" file for both spectral analysis and flux correction
    call ReadIniFCC('FluxCorrection')
    call ResetSpectralAssessmentDiagnostics()

    !> Add run-mode tag to Timestamp_FilePadding
    call TagRunMode()

    !> The run log. The output folder is RP's and already exists by the time
    !> FCC runs; the timestamp is FCC's own, so the two logs never collide.
    call InitRunLog()

    !> If running in embedded mode, override some settings
    if (EddyFlowProj%run_env == 'embedded') &
        call ConfigureForEmbedded('EddyFlow-FCC')

    if (EddyFlowProj%fluxnet_mode) call ConfigureForFluxnet()

    !> Before anything reads a record. InitExVars parses the whole file by
    !> field position, so a file from an older RP does not announce itself -
    !> every record simply fails and the run dies on "no valid data records",
    !> which tells the user nothing about why. The check used to sit after
    !> this and could therefore never fire.
    call CheckExFileVintageAt(AuxFile%ex)

    !> Preliminarily read essential files and retrieve a few information
    call InitExVars(exStartTimestamp, exEndTimestamp, &
        NumExRecords, NumValidExRecords, FirstValidRecord)

    call ReadExRecord(AuxFile%ex, udf, FirstValidRecord, lEx, ValidRecord, EndOfFileReached)

    !> If no good records are found stop execution
    if (NumValidExRecords <= 0) call ExceptionHandler(61)

    !> Retrieve NumberOfPeriods and allocate exTimeSeries
    NumberOfPeriods = NumOfPeriods(exStartTimestamp, exEndTimestamp, DateStep)
    allocate(exTimeSeries(NumberOfPeriods + 1))

    call CreateTimeSeries(exStartTimestamp, exEndTimestamp, &
        DateStep, exTimeSeries, size(exTimeSeries), .true.)

    !> Define MasterTimeSeries for the period to be considered
    if (EddyFlowProj%subperiod) then
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
        if (MasterTimeSeries(1) > exTimeSeries(size(exTimeSeries)) &
            .or. MasterTimeSeries(size(MasterTimeSeries)) < exTimeSeries(1)) &
            call ExceptionHandler(50)
        !> Set start/end indexes in MasterTimeSeries. Note that definition of
        !> fxStartTimestampIndx works, but isn't very nice. Should be improved.
        fxEndTimestampIndx = size(MasterTimeSeries)
        if (EddyFlowProj%make_dataset) then
            fxStartTimestampIndx = 1
        else
            fxStartTimestampIndx = 2
        end if
    else
        allocate(MasterTimeSeries(size(exTimeSeries)))
        MasterTimeSeries = exTimeSeries
        !> Set start/end indexes in MasterTimeSeries.
        fxStartTimestampIndx = 1
        fxEndTimestampIndx = size(MasterTimeSeries)
    end if

    !> Initialize import of full cospectra files
    if (FCCsetup%import_full_cospectra) then
        call NumberOfFilesInDir(Dir%full, '.csv', .false., '', &
            NumFullFiles, NumFullFilesNoRecurse)
        allocate(FullFileList(NumFullFiles))

        !> Read names of full cospectra files
        call FileListByExt(Dir%full, '.csv', .true., .false., &
            FullFilePrototype, .true., .true., .true., FullFilelist, &
            size(FullFilelist), .true., '')

        !> Retrieve length of full cospectra for later allocation
        call FullCospectraLength(FullFilelist(1)%path, nrow_full)
    end if

    !****************************************************************
    !****************************************************************
    !*********** SPECTRAL ASSESSMENT SECTION IF REQUESTED ***********
    !****************************************************************
    !****************************************************************

    !> If spectral analysis must be performed or co-spectral
    !> outputs are requested, start loop on cospectra files
    if (FCCsetup%pass_thru_spectral_assessment) then
        call ResetSpectralAssessmentDiagnostics()
        if (FCCsetup%do_spectral_assessment) write(*, '(a)') &
            ' Starting "spectral assessment" session..'
        if (FCCsetup%do_spectral_assessment) write(ulog, '(a)') &
            ' Starting "spectral assessment" session..'
        if (EddyFlowProj%out_avrg_cosp .or. EddyFlowProj%out_avrg_spec) then
            call LogSay(' Reading (co)spectra from:')
            write(*, '(a)') '  ' // trim(adjustl(Dir%binned))
            write(ulog, '(a)') '  ' // trim(adjustl(Dir%binned))
            call LogSay('')
        end if

        !> Convert start/end to timestamps
        call DateTimeToDateType(FCCsetup%SA%start_date, &
            FCCsetup%SA%start_time, saStartTimestamp)
        call DateTimeToDateType(FCCsetup%SA%end_date, &
            FCCsetup%SA%end_time, saEndTimestamp)

        !> Spectra hold the filestamp of the end of the period, so increase
        !> start timestamp by DateStep
        saStartTimestamp = saStartTimestamp + DateStep

        !> Read names of binned (co)spectra files
        call NumberOfFilesInDir(Dir%binned, '.csv', .false., '', &
            NumBinnedFiles, NumBinnedFilesNoRecurse)

        if(NumBinnedFiles <= 0) then
            !> Exit loop and set spectral correction to Moncrieff.
            !> Also, cannot create any ensemble spectral output
            EddyFlowProj%out_avrg_cosp = .false.
            EddyFlowProj%out_avrg_spec = .false.
            FCCsetup%do_spectral_assessment = .false.
            if (FCCsetup%SA%in_situ) then
                EddyFlowProj%hf_meth = 'moncrieff_97'
                FCCsetup%SA%in_situ = .false.
            end if
            call ExceptionHandler(89)
            goto 100
        end if
        allocate(BinnedFileList(NumBinnedFiles))

        !> Create list of binned spectra file names
        call FileListByExt(Dir%binned, '.csv', .true., .false., &
            BinnedFilePrototype, .true., .true., .true., &
            BinnedFileList, size(BinnedFileList), .true., ' ')

        !> Set file list in chronological order
        call FilesInChronologicalOrder(BinnedFileList, size(BinnedFileList), &
            binStartTimestamp, binEndTimestamp, ' ')

        !> Detect first and last binned files to be used for spectral \n
        !> assessment based on user's dates selection
        call tsExtractSubperiodIndexesFromFilelist(BinnedFileList, &
            size(BinnedFileList), saStartTimestamp, saEndTimestamp, &
            saStartTimestampIndx, saEndTimestampIndx)

        if(saStartTimestampIndx <= 0 .or. saEndTimestampIndx <= 0) then
            !> Exit loop and set spectral correction to Moncrieff.
            !> Also, cannot create any ensemble spectral output
            EddyFlowProj%out_avrg_cosp = .false.
            EddyFlowProj%out_avrg_spec = .false.
            FCCsetup%do_spectral_assessment = .false.
            if (FCCsetup%SA%in_situ) then
                EddyFlowProj%hf_meth = 'moncrieff_97'
                FCCsetup%SA%in_situ = .false.
            end if
            call ExceptionHandler(90)
            goto 100
        end if

        !> Some logging
        call DateTypeToDateTime(binStartTimestamp - DateStep, sDate, sTime)
        call DateTypeToDateTime(binEndTimestamp, eDate, eTime)
        call LogSay('')
        call LogSay('  Period covered by available binned (co)spectra files:')
        write(*, '(a)') '   Start: ' // sDate // ' ' // sTime
        write(ulog, '(a)') '   Start: ' // sDate // ' ' // sTime
        write(*, '(a)') '   End:   ' // eDate // ' ' // eTime
        write(ulog, '(a)') '   End:   ' // eDate // ' ' // eTime

        if (FCCsetup%SA%subperiod) then
            call DateTypeToDateTime(saStartTimestamp, sDate, sTime)
            call DateTypeToDateTime(saEndTimestamp + Datetype(0, 0, 0, 0, 1), eDate, eTime)
            call LogSay('')
            call LogSay('  Selected (co)spectra sub-period:')
            write(*, '(a)') '   Start: ' // sDate // ' ' // sTime
            write(ulog, '(a)') '   Start: ' // sDate // ' ' // sTime
            write(*, '(a)') '   End:   ' // eDate // ' ' // eTime
            write(ulog, '(a)') '   End:   ' // eDate // ' ' // eTime
        end if

        call LogSay('')
        write(LogInteger, '(i8)') saEndTimestampIndx - saStartTimestampIndx + 1
        write(*, '(a)') '  Importing, sorting and ensemble-averaging up to ' &
            // trim(adjustl(LogInteger)) // ' binned (co)spectra from files.. '
        write(ulog, '(a)') '  Importing, sorting and ensemble-averaging up to ' &
            // trim(adjustl(LogInteger)) // ' binned (co)spectra from files.. '

        !> Create an exponentially spaced frequency array in a range \n
        !> wide enough to accommodate any possible normalized frequency
        dkf(1) = 1d0 / (60d0 * 60d0 * 4d0)     !< 1 / (4 hours in seconds)
        dkf(ndkf + 1) = 200d0 / 2d0            !< Max normalized freq of 200 Hz
        do i = 2, ndkf
            dkf(i) = dkf(1) * dexp(dble(i - 1) &
                * (dlog(dkf(ndkf + 1)) - dlog(dkf(1))) / dble(ndkf))
        end do

        !> Open Ex file to keep it ready for reading
        !> and exit with error in case of problems opening the file
        open(uex, file = AuxFile%ex, status = 'old', iostat = open_status)
        if (open_status /= 0) call ExceptionHandler(60)

        !> Skip header in Ex file, after checking it was written by this
        !> version - see CheckExFileVintage.
        call CheckExFileVintage()

        !> Loop to import binned (co)spectra
        month = 0
        day   = 0
        nfit = 0
        allocate(FitStable(0))
        allocate(FitUnstable(0))
        fcount = saStartTimestampIndx - 1
        binned_loop: do
            !> Update file counter
            fcount = fcount + 1

            !> Normal exit instruction
            if (fcount > saEndTimestampIndx) exit binned_loop

            !> Read (co)spectra from file
            SADiagSelectedFiles = SADiagSelectedFiles + 1
            call ReadBinnedFile(BinnedFileList(fcount), BinSpec, BinCosp, &
                size(BinSpec), nbins, skip)
            if (skip) cycle binned_loop
            SADiagReadableFiles = SADiagReadableFiles + 1

            !> Show advancement
            if (day /= BinnedFileList(fcount)%timestamp%Day &
                .or. month /= BinnedFileList(fcount)%timestamp%Month) then
                month = BinnedFileList(fcount)%timestamp%Month
                day   = BinnedFileList(fcount)%timestamp%Day
                call DisplayProgress('daily', &
                    '  Importing binned (co)spectra for ', &
                    BinnedFileList(fcount)%timestamp, 'yes')
            end if

            !> Retrieve ex information for current spectra
            call RetrieveExVarsByTimestamp(uex, &
                BinnedFileList(fcount)%timestamp, lEx, exEndReached, skip)

            if (exEndReached) exit binned_loop
            if (skip) cycle binned_loop
            SADiagMatchedRecords = SADiagMatchedRecords + 1

            !> Allocate variables that depend upon nbins and initialize them
            if (.not. allocated(MeanBinSpec)) then
                allocate(MeanBinSpec(nbins, MaxGasClasses))
                allocate(dMeanBinSpec(nbins, MaxGasClasses))
                MeanBinSpec = NullMeanSpec
                dMeanBinSpec = NullMeanSpec
            end if
            if (.not. allocated(MeanBinCosp)) then
                allocate(MeanBinCosp(nbins, MaxGasClasses))
                MeanBinCosp = NullMeanSpec
                deallocate(FitUnstable)
                allocate(FitUnstable(nbins * &
                    (saEndTimestampIndx - saStartTimestampIndx + 1)))
                FitUnstable = NullFitCosp
                deallocate(FitStable)
                allocate(FitStable  (nbins * &
                    (saEndTimestampIndx - saStartTimestampIndx + 1)))
                FitStable   = NullFitCosp
            end if

            !> Eliminate (co)spectra based on user-selected quality criteria
            call CospectraQAQC(BinSpec, BinCosp, size(BinSpec), lEx, &
                BinCospForStable, BinCospForUnstable, &
                skip_spectra, skip_cospectra)

            !> Sort current spectra in relevant classes
            if (.not. skip_spectra) &
                call SpectraSortingAndAveraging(lEx, BinSpec, &
                    size(BinSpec), nbins)

            !> Sort current cospectra in time-slot classes
            if (EddyFlowProj%out_avrg_cosp .and. .not. skip_cospectra) then

                !> Add current cospectra to dataset for regression
                call AddToCospectraFitDataset(lEx, BinCospForStable, &
                    BinCospForUnstable, size(BinCospForStable), nfit, &
                    size(nfit, 1), size(nfit, 2), nbins, FitStable, &
                    FitUnstable, size(FitStable))

                !> Sort current cospectra in time slot classes
                call CospectraSortingAndAveraging(BinCosp, size(BinCosp), &
                    lEx%end_time, nbins)
            end if
        end do binned_loop
        close(uex)
        call LogSay('  Done.')

        !> Write number of imported spectra and cospectra on stdout
        if (EddyFlowProj%out_avrg_spec .or. FCCsetup%do_spectral_assessment) &
            call ReportImportedSpectra(nbins)

        if (EddyFlowProj%out_avrg_cosp .and. allocated(FitStable)) then
            !> If cospectra were found for fitting, fit Massman model
            call FitCospectralModel(nfit, size(nfit, 1), size(nfit, 2), &
                FitStable, FitUnstable, size(FitStable))

            !> Ensemble average cospectra in stable and unstable stratifications
            call EnsembleCospectraByStability(nfit, size(nfit, 1), &
                size(nfit, 2), FitStable, FitUnstable, size(FitStable))

        end if

        !> Normalize sums for obtaining mean spectra
        call NormalizeMeanSpectraCospectra(nbins)

        !> Detect which average spectra are available
        !> for each gas and each class
        call AvailableMeanSpectraCospectra(nbins)

        !> Determine low-pass TF cut-off frequencies, RH-sorted (H2O)
        !> and time-sorted (CO2/CH4/GAS4)
        call FitTFModels(nbins, FCCsetup%do_spectral_assessment)

        !> Spectral attenuation assessment
        if (FCCsetup%do_spectral_assessment) then

            !> Determine analytical relation fc/RH
            call FitRh2Fco()

            !> If necessary, calculate spectral correction factor models
            !> as from Ibrom et al. (2007)
            call CorrectionFactorModel(AuxFile%ex, NumExRecords)
        else
            !> If an in-situ method was chosen, and spectral
            !> assessment file is available, read file
            if (FCCsetup%SA%in_situ) call ReadSpectralAssessmentFile()
        end if

        call ReportSpectralAssessmentDiagnostics(skip_spectra)

        !> Write everything on output files
        call OutputSpectralAssessmentResults(nbins)
        write(*,'(a)')
        write(ulog,'(a)')


    else
        !> If an in-situ method was chosen, and spectral
        !> assessment file is available, read file
        if (FCCsetup%SA%in_situ) call ReadSpectralAssessmentFile()
    end if

100 continue

    !***************************************************************************
    !***************************************************************************
    !****** MAIN CYCLE ON RESULTS RECORDS RETRIEVED FROM ESSENTIALS FILE *******
    !***************************************************************************
    !***************************************************************************

    !> Establish present variables
    ! call EstablishPresentVariables()

    !> Open Ex file to keep it ready for reading
    !> and exit with error in case of problems opening the file
    open(uex, file = AuxFile%ex, status = 'old', iostat = open_status)
    if (open_status /= 0) call ExceptionHandler(60)
    !> The header was skipped unread. ReadExRecord parses this file by field
    !> position, so a file written by an older RP is not merely out of date -
    !> its per-gas moisture records are three fields wide where the reader now
    !> expects seven, and a list-directed read simply continues into the next
    !> gas's fields and returns plausible numbers for the wrong slot. Nothing
    !> downstream can notice. Checking one column name turns that into a stop.
    call CheckExFileVintage()



    month = 0
    day   = 0
    InitializeOuputFiles = .true.
    ex_loop: do i = 1, NumExRecords

        !> Read record from essentials file
        call ReadExRecord('', uex, -1, lEx, ValidRecord, EndOfFileReached)

        !> Initialize presence of key variables for outputting results
        ! if (InitializeOuputFiles) &
        !     fcc_var_present(u:GHGNumVar) = lEx%var_present(u:GHGNumVar)

        !> If end of file was reached, exit loop
        if (EndOfFileReached) exit ex_loop

        !> If invalid record was found, cycle loop
        if (.not. ValidRecord) cycle ex_loop

        !> Retrieve timestamp
        call DateTimeToDateType(lEx%end_date, lEx%end_time, CurrentTimestamp)

        !> If current timestamp is < start selected timestamp, cycle
        if (CurrentTimestamp < MasterTimeSeries(fxStartTimestampIndx)) &
            cycle ex_loop

        !> If current timestamp is > end selected range, exit
        if (CurrentTimestamp > MasterTimeSeries(fxEndTimestampIndx)) &
            exit ex_loop

        !> Show advancement
        call DateTimetoDOY(lEx%end_date, lEx%end_time, int_doy, float_doy)
        if (day /= CurrentTimestamp%day &
            .or. month /= CurrentTimestamp%month) then
            month = CurrentTimestamp%month
            day   = CurrentTimestamp%day
            call DisplayProgress('daily', '  Calculating fluxes for ', &
                CurrentTimestamp, 'yes')
        end if

        !> Band-pass spectral correction factors
        BPCF%of(:) = 1d0

        !> Create aux variables to pass to BandPassSpectralCorrections
        AuxInstrument = NullInstrument
        AuxInstrument(sonic) = lEx%instr(sonic)
        !> Indexed by gas *slot*, not by instrument role. lEx%instr is indexed
        !> by role (ico2..igas4) and so only ever reaches four gases; past that
        !> the role index addresses an unrelated instrument. lEx%gas_instr is
        !> the per-slot view, mirrored from those four after their unit
        !> conversions, so the historical slots are unchanged.
        do gas = firstGas, lastGas
            AuxInstrument(gas) = lEx%gas_instr(gas)
        end do
        if (.not. allocated(FullFileList)) allocate(FullFileList(1))

        !> Spectral correction and the two flux levels, once or repeatedly.
        !>
        !> Twinned with the loop in eddyflow-rp_main.f90 and for the same
        !> reason: the analytic cospectrum is evaluated at a stability the
        !> corrected heat flux itself determines. Off by default and a single
        !> pass then, which is exactly what this block did before.
        !>
        !> The stability is threaded through a local rather than read from
        !> lEx twice. BandPassSpectralCorrections is handed lEx%Flux0%zL - the
        !> Level-0 value the ex record carries - while Fluxes23 writes
        !> lEx%zL, so the two are different fields and the feedback would
        !> otherwise not close. RP has the same shape with one name,
        !> Ambient%zL, which is why the connection is easy to miss here.
        iter_dev = error
        corr_passes = 1
        if (EddyFlowProj%corr_iter_meth) corr_passes = EddyFlowProj%corr_iter_max
        iter_zL = lEx%Flux0%zL
        do corr_pass = 1, corr_passes
            if (corr_pass > 1) prev_gas_flux = Flux3%gas

            !> Bad pass spectral correction factors
            call BandPassSpectralCorrections(lEx%instr(sonic)%height, &
                lEx%disp_height, lEx%var_present, lEx%WS, lEx%Ta, iter_zL, &
                lEx%ac_freq, nint(lEx%avrg_length), lEx%logger_swver, &
                lEx%det_meth, nint(lEx%det_timec), .false., AuxInstrument, &
                size(FullFileList), FullFileList, nrow_full, lEx, FCCsetup)

            !> Calculate fluxes at Level 1
            call Fluxes1(lEx)

            !> Calculate fluxes at Level 2 and Level 3
            call Fluxes23(lEx)

            !> Fluxes23 leaves the previous value standing when it cannot
            !> form a new one, so this is either the refined stability or the
            !> one the pass started from - never an error code.
            iter_zL = lEx%zL

            if (corr_pass > 1) then
                iter_dev = WorstRelativeChange(prev_gas_flux, Flux3%gas)
                !> A tolerance of zero never fires, which is EddyUH's
                !> behaviour: it runs every pass and tests nothing.
                if (EddyFlowProj%corr_iter_tol > 0d0 &
                    .and. iter_dev /= error &
                    .and. iter_dev < EddyFlowProj%corr_iter_tol) exit
            end if
        end do
        lEx%corr_iter_dev = iter_dev

        !> Apply RP's high-frequency CEC descriptors to FCC's authoritative
        !> corrected totals - one per pairing.
        !>
        !> The target slots come from the descriptor, not from the project: the
        !> descriptor is what RP actually computed, and reading the totals for
        !> some other list would pair a ratio with a flux it does not describe.
        do cec_p = 1, MaxNumCecPairs
            call ResetCecFlux(CECFlux(cec_p))
        end do
        if (EddyFlowProj%do_cec > 0) then
            do cec_p = 1, min(lEx%n_cec, MaxNumCecPairs)
                cec_totals = error
                cec_errors = error
                do cec_k = 1, lEx%cec(cec_p)%n_target
                    cec_slot = lEx%cec(cec_p)%target(cec_k)%slot
                    if (cec_slot < firstGas .or. cec_slot > lastGas) cycle
                    cec_totals(cec_k) = Flux3%gas(cec_slot)
                    !> RP estimated it and the essentials row carried it here,
                    !> from the same slot as the total beside it.
                    cec_errors(cec_k) = lEx%rand_uncer(cec_slot)
                end do
                call ApplyCecDescriptor(lEx%cec(cec_p), cec_totals, &
                    cec_errors, EddyFlowProj%cec, CECFlux(cec_p))
            end do
        end if

        !> Calculate footprint estimation   
        foot_model_used = Meth%foot(1:len_trim(Meth%foot))
        call FootprintHandle(lEx%var(w), lEx%ustar, lEx%zL, lEx%WS, lEx%L, &
            lEx%instr(sonic)%height, lEx%disp_height, lEx%rough_length)

        !> Calculate quality flags
        StDiff%w_u    = nint(lEx%TAU_SS)
        StDiff%w_ts   = nint(lEx%H_SS)
        !> Every configured gas. lEx%F_SS is slot-indexed and read_ex_record
        !> already fills firstGas..firstGas+n_layout_gas-1 from the file, so
        !> the data was there all along - only this copy stopped at the fourth
        !> slot, leaving QualityFlags to make a flag out of unset memory for
        !> every gas past it. That reached the full output as qc_LE whenever
        !> the site's water sat further out.
        StDiff%w_gas = nint(error)
        do gas = firstGas, lastGas
            if (.not. fcc_var_present(gas)) cycle
            StDiff%w_gas(gas) = nint(lEx%F_SS(gas))
        end do
        DtDiff%u      = nint(lEx%U_ITC)
        DtDiff%w      = nint(lEx%W_ITC)
        DtDiff%ts     = nint(lEx%TS_ITC)
        call QualityFlags(Flux2, StDiff, DtDiff, STFlg, DTFlg, QCFlag, .false.)

        !> Initialize output files
        if (InitializeOuputFiles) then
            call InitOutFiles(lEx)
            InitializeOuputFiles = .false.
        end if

        if (EddyFlowProj%out_full .and. .not. lEx%not_enough_data) call WriteOutFullFcc(lEx)
        if (EddyFlowProj%out_md .and. .not. lEx%not_enough_data) call WriteOutMetadataFcc(lEx)
        if (EddyFlowProj%out_fluxnet) call WriteOutFluxnetFcc(lEx)

    end do ex_loop
    close(uex)
    close(uflx)
    close(ufnet_e)
    close(uaflx)
    close(umd)
    close(uflxnt)

    write(*,*)
    write(ulog,*)
    write(*,*)
    write(ulog,*)
    call sleep(1)

    !> Creating datasets from output files
    if (EddyFlowProj%make_dataset) then
        call CreateDatasetsCommon(MasterTimeSeries, size(MasterTimeSeries), &
            fxStartTimestampIndx, fxEndTimestampIndx, 'FCC')
    else
        call RenameTmpFilesCommon()
    end if

    !> Delete tmp folder if running in embedded mode
    if(EddyFlowProj%run_env == 'desktop') &
        del_status = system(trim(comm_rmdir) // ' "' &
        // trim(adjustl(TmpDir)) // '"')

    !> Delete parent fluxnet file
    if (.not. FCCsetup%keep_parent) then
        call system(comm_del // '"' // trim(adjustl(AuxFile%ex)) // '"' // comm_err_redirect)
    end if

    !> Copy ".eddypro" file into output folder
    call CopyFile(trim(adjustl(PrjPath)), &
    trim(adjustl(Dir%main_out)) // 'processing' &
    // Timestamp_FilePadding // '.eddyflow')
    call ApplyAutomaticSpectralConfiguration(trim(adjustl(Dir%main_out)) // 'processing' &
        // Timestamp_FilePadding // '.eddyflow')


    call LogSay('')
    call LogSay(' ****************************************************')
    call LogSay(' Program EddyFlow executed gracefully.')
    call LogSay(' Check results in the selected output directory.     ')
    call LogSay(' ****************************************************')
    stop ''

contains

!***************************************************************************
!
! \brief       Read past the ex-file header, refusing one this version cannot
!              parse.
! \author      Jonathan Muller
! \note        Both places that open the ex file used to `read(uex, *)` and
!              throw the header away. That was harmless while the format only
!              ever grew at the end - a reader stops when it has what it wants.
!              It stopped being harmless when the per-gas moisture records went
!              from three fields to seven: a list-directed read of seven values
!              from a three-field record does not fail, it continues into the
!              next gas's fields and returns numbers that look entirely
!              ordinary for the wrong slot, and every flux computed from them
!              is quietly wrong.
!
!              So one column name is checked. NUM_WATER_FLUX appears only in
!              files written by this version, and its absence means the
!              moisture records are the narrow ones. Re-running RP regenerates
!              the file; there is nothing to migrate.
!***************************************************************************
!> Read the header off an already-open essentials unit, and judge it. Doubles
!> as the header skip both callers need.
subroutine CheckExFileVintage()
    implicit none
    character(LongOutstringLen) :: header

    header = ''
    read(uex, '(a)', iostat = open_status) header
    if (open_status /= 0) call ExceptionHandler(60)
    call JudgeExHeader(header)
end subroutine CheckExFileVintage

!> The same judgement, on a file nothing has opened yet.
!>
!> Wanted separately because the only useful place to make it is before
!> InitExVars, which is before any unit is open - and after InitExVars is too
!> late, the records having already failed to parse.
subroutine CheckExFileVintageAt(path)
    implicit none
    character(*), intent(in) :: path

    integer :: unt
    integer :: ios
    character(LongOutstringLen) :: header

    open(newunit = unt, file = path, status = 'old', iostat = ios)
    if (ios /= 0) call ExceptionHandler(60)
    header = ''
    read(unt, '(a)', iostat = ios) header
    close(unt)
    if (ios /= 0) call ExceptionHandler(60)
    call JudgeExHeader(header)
end subroutine CheckExFileVintageAt

!> Which column names a file must carry to be readable by this version.
subroutine JudgeExHeader(header)
    implicit none
    character(*), intent(in) :: header

    !> Two markers, because the row grew twice. NUM_WATER_FLUX arrived with
    !> the per-hygrometer families and H2O_BIOMET_MOLE_FRACTION with the
    !> biomet triple, which sits in the *fixed* part - a file carrying the
    !> first but not the second parses three fields short from there on.
    if (index(header, 'NUM_WATER_FLUX') <= 0 &
        .or. index(header, 'H2O_BIOMET_MOLE_FRACTION') <= 0) &
        call ExceptionHandler(107)

    !> The third marker is conditional, because a project with the partition
    !> off writes no CEC block at all and a header without one is not old. But
    !> a header that HAS the block and lacks CEC_NS_ was written before the
    !> partition-stability statistic joined the per-target fields, and would
    !> parse one field short per target from there to the end of the row.
    if (index(header, 'CEC_METH') > 0 .and. index(header, 'CEC_NS_') <= 0) &
        call ExceptionHandler(107)
end subroutine JudgeExHeader

end program EddyFlowFCC
