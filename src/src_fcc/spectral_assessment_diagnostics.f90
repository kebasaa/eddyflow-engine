!*******************************************************************************
! spectral_assessment_diagnostics.f90
! ----------------------------------
! Summarise the readiness and outcome of FCC spectral assessment.
!*******************************************************************************
subroutine ResetSpectralAssessmentDiagnostics()
    use m_fx_global_var
    implicit none

    SADiagSelectedFiles = 0
    SADiagReadableFiles = 0
    SADiagMatchedRecords = 0
    SADiagUsableWT = 0
    SADiagRejectedUstar = 0
    SADiagRejectedVM = 0
    SADiagRejectedFoken = 0
    SADiagRejectedFlux = 0
    SADiagAccepted = 0
    SADiagDegradedUnstable = 0
    SADiagDegradedStable = 0
    SADiagFluxCandidateCount = 0
    SADiagFluxCandidateCapacity = 0
    SAAutoMin = error
    SAAutoMax = error
    SAAutoApplyMin = .false.
    SAAutoApplyMax = .false.
    SADiagFilePath = 'none'
    if (allocated(SADiagFluxCandidate)) deallocate(SADiagFluxCandidate)
    if (allocated(SADiagFluxCandidateClass)) deallocate(SADiagFluxCandidateClass)
end subroutine ResetSpectralAssessmentDiagnostics

!*******************************************************************************
subroutine RecordSpectralAssessmentFluxCandidate(gas, stability, flux, cls)
    use m_fx_global_var
    implicit none
    integer, intent(in) :: gas
    integer, intent(in) :: stability
    integer, intent(in) :: cls
    real(kind = dbl), intent(in) :: flux
    integer :: new_capacity
    real(kind = dbl), allocatable :: tmp_flux(:, :, :)
    integer, allocatable :: tmp_class(:, :, :)

    !> Every configured gas. Bounded at the fourth, the readiness report
    !> could only ever say flux=0 for a gas past it, whatever the data said.
    if (gas < firstGas .or. gas > lastGas) return
    if (stability < SADiagUnstable .or. stability > SADiagStable) return
    if (cls < 1 .or. cls > MaxGasClasses) return
    if (flux == error .or. flux < 0d0) return

    if (.not. allocated(SADiagFluxCandidate)) then
        SADiagFluxCandidateCapacity = 256
        allocate(SADiagFluxCandidate(SADiagFluxCandidateCapacity, 2, GHGNumVar))
        allocate(SADiagFluxCandidateClass(SADiagFluxCandidateCapacity, 2, GHGNumVar))
        SADiagFluxCandidate = error
        SADiagFluxCandidateClass = 0
    elseif (SADiagFluxCandidateCount(stability, gas) == SADiagFluxCandidateCapacity) then
        new_capacity = SADiagFluxCandidateCapacity * 2
        allocate(tmp_flux(new_capacity, 2, GHGNumVar))
        allocate(tmp_class(new_capacity, 2, GHGNumVar))
        tmp_flux = error
        tmp_class = 0
        tmp_flux(1:SADiagFluxCandidateCapacity, :, :) = SADiagFluxCandidate
        tmp_class(1:SADiagFluxCandidateCapacity, :, :) = SADiagFluxCandidateClass
        call move_alloc(tmp_flux, SADiagFluxCandidate)
        call move_alloc(tmp_class, SADiagFluxCandidateClass)
        SADiagFluxCandidateCapacity = new_capacity
    end if

    SADiagFluxCandidateCount(stability, gas) = SADiagFluxCandidateCount(stability, gas) + 1
    SADiagFluxCandidate(SADiagFluxCandidateCount(stability, gas), stability, gas) = flux
    SADiagFluxCandidateClass(SADiagFluxCandidateCount(stability, gas), stability, gas) = cls
end subroutine RecordSpectralAssessmentFluxCandidate

!*******************************************************************************
subroutine ReportSpectralAssessmentDiagnostics(assessment_ready)
    use m_fx_global_var
    implicit none
    logical, intent(out) :: assessment_ready
    integer, external :: CreateDir
    integer :: cls
    integer :: gas
    integer :: h2o_classes
    integer :: h2o_fits
    integer :: gas_classes
    !> The site's water record. These counters are water's, and they read the
    !> sixth slot whatever species that record held.
    integer :: wsl
    integer :: report_unit
    integer :: mkdir_status
    integer :: open_status
    character(PathLen) :: SpecDir
    character(PathLen) :: FilePath
    character(128) :: Filename
    character(32) :: method
    logical :: gas_ok
    logical :: corrected_ok
    integer :: donor
    integer :: n_insitu
    integer :: n_analytic
    character(32), external :: IntToText
    character(16), external :: GasName
    character(12), external :: StatusLabel
    character(12), external :: OutcomeLabel
    logical, external :: GasHasSpectralFit
    logical, external :: HygrometerHasRhFit
    logical, external :: GasInsituWithSource
    character(64), external :: DonorNote
    include '../src_common/interfaces_1.inc'

    assessment_ready = .true.
    !> The site's hygrometer, or none. With no water the fallback slot six
    !> holds another gas, and its classes were counted - and reported - as
    !> water's RH classes.
    wsl = PrimaryWaterSlot()
    method = trim(adjustl(EddyFlowProj%hf_meth))
    h2o_classes = 0
    h2o_fits = 0
    if (allocated(MeanBinSpec) .and. wsl >= firstGas) then
        do cls = RH10, RH90
            if (MeanBinSpecAvailable(cls, wsl)) h2o_classes = h2o_classes + 1
            if (RegPar(wsl, cls)%fc /= error) h2o_fits = h2o_fits + 1
        end do
    end if


    if (FCCsetup%do_spectral_assessment) then
        !> The file is worth writing when any gas got a fit, not when every
        !> one did. It used to be h2o_ok .and. (every gas ok), so a site whose
        !> hygrometers could not be fitted lost the blocks of the gases that
        !> had been - and the report called them PASS on the way past.
        n_insitu = 0
        n_analytic = 0
        do gas = firstGas, lastGas
            if (.not. fcc_var_present(gas)) cycle
            gas_ok = GasHasSpectralFit(gas)
            !> A hygrometer needs the RH regression on top of its own fit:
            !> that is what turns humidity into a cut-off frequency. Its own
            !> regression, or the project-wide one it falls back to - the
            !> pair the correction evaluates, not the primary's for all.
            if (GasSlotIsWater(gas)) then
                if (gas_ok) gas_ok = HygrometerHasRhFit(gas)
            end if
            if (gas_ok) then
                n_insitu = n_insitu + 1
            else
                n_analytic = n_analytic + 1
            end if
        end do
        assessment_ready = n_insitu > 0
    end if

    mkdir_status = CreateDir('"' // Dir%main_out(1:len_trim(Dir%main_out)) // '"')
    SpecDir = Dir%main_out(1:len_trim(Dir%main_out)) // SubDirSpecAn // slash
    mkdir_status = CreateDir('"' // SpecDir(1:len_trim(SpecDir)) // '"')
    Filename = EddyFlowProj%id(1:len_trim(EddyFlowProj%id)) // &
        '_spectral_correction_diagnostics' // Timestamp_FilePadding // TxtExt
    FilePath = SpecDir(1:len_trim(SpecDir)) // Filename(1:len_trim(Filename))
    SADiagFilePath = FilePath
    open(newunit = report_unit, file = FilePath, status = 'replace', &
        action = 'write', iostat = open_status)

    call EmitReportLine(report_unit, open_status, '')
    call EmitReportLine(report_unit, open_status, 'Spectral corrections readiness report')
    call EmitReportLine(report_unit, open_status, '-------------------------------------')
    call EmitReportLine(report_unit, open_status, 'Selected method: ' // trim(method) // '... ' // &
        trim(StatusLabel(FCCsetup%SA%in_situ, 'FAIL')))
    if (FCCsetup%do_spectral_assessment) then
        call EmitReportLine(report_unit, open_status, 'Assessment mode: on-the-fly... PASS')
    else
        call EmitReportLine(report_unit, open_status, 'Assessment mode: existing file or not requested... NOT APPLIED')
    end if
    call EmitReportLine(report_unit, open_status, 'Binned spectra input files: ' // &
        trim(IntToText(SADiagReadableFiles)) // '/' // trim(IntToText(SADiagSelectedFiles)) // '... ' // &
        trim(StatusLabel(SADiagReadableFiles > 0 .or. .not. FCCsetup%do_spectral_assessment, 'FAIL')))
    if (method == 'fratini_12') then
        call EmitReportLine(report_unit, open_status, 'Full cospectra (Fratini et al. 2012)... PASS')
    else
        call EmitReportLine(report_unit, open_status, 'Full cospectra (Fratini et al. 2012)... NOT APPLIED')
    end if
    call EmitReportLine(report_unit, open_status, 'Matching essentials records: ' // &
        trim(IntToText(SADiagMatchedRecords)) // '... ' // &
        trim(StatusLabel(SADiagMatchedRecords > 0 .or. .not. FCCsetup%do_spectral_assessment, 'FAIL')))
    call EmitReportLine(report_unit, open_status, 'Usable w/T spectra: ' // &
        trim(IntToText(SADiagUsableWT)) // '... ' // &
        trim(StatusLabel(SADiagUsableWT > 0 .or. .not. FCCsetup%do_spectral_assessment, 'FAIL')))
    call EmitReportLine(report_unit, open_status, 'u* filter rejections: ' // &
        trim(IntToText(SADiagRejectedUstar)) // '... ' // &
        trim(StatusLabel(SADiagRejectedUstar == 0, 'WARNING')))
    if (FCCsetup%SA%filter_cosp_by_vm_flags) then
        call EmitReportLine(report_unit, open_status, 'VM filtering... PASS (rejections reported per gas)')
    else
        call EmitReportLine(report_unit, open_status, 'VM filtering... NOT APPLIED')
    end if
    if (FCCsetup%SA%foken_lim >= 0) then
        call EmitReportLine(report_unit, open_status, 'Foken filtering (limit ' // &
            trim(IntToText(FCCsetup%SA%foken_lim)) // ')... PASS (rejections reported per gas)')
    else
        call EmitReportLine(report_unit, open_status, 'Foken filtering... NOT APPLIED')
    end if
    if (FCCsetup%SA%automatic_config) then
        call EmitReportLine(report_unit, open_status, &
            'Automatic spectral configuration: ENABLED - qualifying recommendations will be saved to the output project.')
    else
        call EmitReportLine(report_unit, open_status, 'Automatic spectral configuration: DISABLED')
    end if

    !> Only gases the project configures. Running the full slot range would
    !> emit sixty "unavailable" lines for slots that do not exist.
    do gas = firstGas, lastGas
        if (gas - firstGas + 1 > min(EddyFlowProj%gas_num, MaxNumGases)) exit
        gas_classes = 0
        if (allocated(MeanBinSpec)) then
            do cls = 1, MaxGasClasses
                if (MeanBinSpecAvailable(cls, gas)) gas_classes = gas_classes + 1
            end do
        end if
        !> What this gas will actually get. A count of ensemble classes is not
        !> the same question: the correction needs a fitted cut-off, and a
        !> hygrometer needs the RH regression on top of it.
        gas_ok = GasHasSpectralFit(gas)
        if (GasSlotIsWater(gas)) then
            if (gas_ok) gas_ok = HygrometerHasRhFit(gas)
        end if
        !> PASS/FAIL is this gas's own fit; the outcome is what corrects it,
        !> which may be a tube-mate's fit (SpectralSourceSlot) - in situ when
        !> a gas that failed borrows a good one, analytic when the donor it
        !> names has nothing.
        corrected_ok = gas_ok
        donor = 0
        if (.not. GasSlotIsWater(gas)) &
            corrected_ok = GasInsituWithSource(gas, donor)
        if (.not. fcc_var_present(gas)) then
            call EmitReportLine(report_unit, open_status, trim(GasName(gas)) // ': unavailable... NOT APPLIED')
        else
            call EmitReportLine(report_unit, open_status, trim(GasName(gas)) // ': accepted periods=' // &
                trim(IntToText(SADiagAccepted(gas))) // ', flux=' // trim(IntToText(SADiagRejectedFlux(gas))) // &
                ', VM=' // trim(IntToText(SADiagRejectedVM(gas))) // ', Foken=' // &
                trim(IntToText(SADiagRejectedFoken(gas))) // ', groups=' // &
                trim(IntToText(FCCsetup%SA%nclass(gas))) // ', valid classes=' // &
                trim(IntToText(gas_classes)) // '... ' // &
                trim(StatusLabel(gas_ok, 'FAIL')) // &
                '  -> corrected ' // trim(OutcomeLabel(corrected_ok)) // &
                trim(DonorNote(donor)))
            if (.not. gas_ok .and. SADiagRejectedFlux(gas) > 0) &
                call ReportFluxLimitSuggestions(report_unit, open_status, gas)
        end if
    end do
    if (wsl < firstGas) then
        call EmitReportLine(report_unit, open_status, 'H2O RH classes: no hygrometer configured... NOT APPLIED')
    else
        call EmitReportLine(report_unit, open_status, 'H2O RH classes: ' // trim(IntToText(h2o_classes)) // &
            '/9, fitted cut-offs=' // trim(IntToText(h2o_fits)) // ', minimum samples/class=' // &
            trim(IntToText(FCCsetup%SA%min_smpl)) // '... ' // &
            trim(StatusLabel(h2o_classes >= 1, 'FAIL')))
    end if
    if (h2o_classes > 0 .and. h2o_classes < 3) then
        call EmitReportLine(report_unit, open_status, &
            'H2O RH coverage: WARNING - one class can create a file; three or more are recommended.')
    end if
    call EmitReportLine(report_unit, open_status, 'Valid degraded wT covariance: unstable=' // &
        trim(IntToText(SADiagDegradedUnstable)) // ', stable=' // trim(IntToText(SADiagDegradedStable)) // '... ' // &
        trim(StatusLabel(SADiagDegradedUnstable + SADiagDegradedStable > 0, 'WARNING')))

    if (FCCsetup%do_spectral_assessment) then
        if (.not. assessment_ready) then
            call EmitReportLine(report_unit, open_status, &
                'Spectral assessment: FAILED - no gas was fitted, so no assessment file is created.')
            call EmitReportLine(report_unit, open_status, &
                'Standard EddyFlow errors and Moncrieff fallback will follow during output/correction.')
            call EmitReportLine(report_unit, open_status, &
                'Adjust the reported QA thresholds/filters or select more qualifying periods.')
        else if (n_analytic == 0) then
            call EmitReportLine(report_unit, open_status, &
                'Spectral assessment: SUCCESS - assessment file will be created, all ' // &
                trim(IntToText(n_insitu)) // ' gases corrected in situ.')
        else
            !> The case this report exists for. Naming both counts is the point:
            !> a file covering four gases out of six looks exactly like a file
            !> covering all six unless the difference is stated.
            call EmitReportLine(report_unit, open_status, &
                'Spectral assessment: PARTIAL - assessment file will be created. ' // &
                trim(IntToText(n_insitu)) // ' gas(es) corrected in situ, ' // &
                trim(IntToText(n_analytic)) // ' with the analytical method instead.')
            call EmitReportLine(report_unit, open_status, &
                'The per-gas lines above say which is which; a gas marked FAIL is corrected analytically.')
            call EmitReportLine(report_unit, open_status, &
                'That is not an error - it is the best available for those gases - but their')
            call EmitReportLine(report_unit, open_status, &
                'correction reads instrument geometry rather than their own measured spectra.')
        end if
    end if
    if (open_status == 0) close(report_unit)
    write(*, '(a)') ' Spectral-correction diagnostics written to: ' // trim(FilePath)
    write(ulog, '(a)') ' Spectral-correction diagnostics written to: ' // trim(FilePath)
end subroutine ReportSpectralAssessmentDiagnostics

!*******************************************************************************
!> The periods in which a gas the in-situ spectral correction was meant to take
!> was corrected with Moncrieff et al. (1997) instead, one line per gas with
!> the count, then each reason with the timestamps it applies to. Written at
!> the end of the FCC run, once the workers' periods are merged in; a run
!> without any writes nothing.
!>
!> The readiness report is written before the flux loop and can only say what
!> the assessment fitted. A gas loses its correction here period by period -
!> no humidity to find a hygrometer's cut-off from, a failed fit for this
!> period's stability - and its factor used to read as a correction.
subroutine ReportInsituFallbacks()
    use m_fx_global_var
    implicit none
    integer :: gas
    integer :: reason
    integer :: i
    integer :: n
    integer :: nr
    integer :: pos
    logical :: first_reason
    logical :: first_stamp
    character(:), allocatable :: line
    character(32), external :: IntToText
    character(16), external :: GasName
    character(20), parameter :: ReasonText(2) = &
        [character(20) :: 'no cut-off frequency', 'failed fit']

    if (InsituFallbackN < 1) return
    call LogSay('')
    call LogSay(' In-situ spectral correction not possible for some periods; ' // &
        'corrected with Moncrieff et al. (1997) instead:')
    do gas = firstGas, lastGas
        n = count(InsituFallbacks(1:InsituFallbackN)%gas == gas)
        if (n == 0) cycle
        !> Sized once: a run can hold thousands of such periods, and growing
        !> the line by concatenation would copy it for every one.
        allocate(character(64 + 2 * 24 + 18 * n) :: line)
        pos = 0
        call Put('  ' // trim(GasName(gas)) // ': ' // trim(IntToText(n)))
        if (n == 1) then
            call Put(' period')
        else
            call Put(' periods')
        end if
        first_reason = .true.
        do reason = InsituNoCutoff, InsituFitFailed
            nr = 0
            do i = 1, InsituFallbackN
                if (InsituFallbacks(i)%gas == gas .and. InsituFallbacks(i)%reason == reason) &
                    nr = nr + 1
            end do
            if (nr == 0) cycle
            if (first_reason) then
                call Put(' - ')
            else
                call Put('; ')
            end if
            first_reason = .false.
            call Put(trim(ReasonText(reason)) // ':')
            first_stamp = .true.
            do i = 1, InsituFallbackN
                if (InsituFallbacks(i)%gas /= gas .or. InsituFallbacks(i)%reason /= reason) cycle
                if (first_stamp) then
                    call Put(' ')
                else
                    call Put(', ')
                end if
                first_stamp = .false.
                call Put(InsituFallbacks(i)%date // ' ' // InsituFallbacks(i)%time)
            end do
        end do
        call LogSay(line(1:pos))
        deallocate(line)
    end do

contains

    subroutine Put(s)
        character(*), intent(in) :: s
        line(pos + 1:pos + len(s)) = s
        pos = pos + len(s)
    end subroutine Put
end subroutine ReportInsituFallbacks

!*******************************************************************************
subroutine ReportFluxLimitSuggestions(report_unit, open_status, gas)
    use m_fx_global_var
    implicit none
    integer, intent(in) :: report_unit
    integer, intent(in) :: open_status
    integer, intent(in) :: gas
    integer :: stability
    integer :: n
    integer :: i
    integer :: projected
    integer :: valid_classes
    integer :: current_valid_classes
    integer :: class_counts(MaxGasClasses)
    integer :: current_class_counts(MaxGasClasses)
    real(kind = dbl) :: current_min
    real(kind = dbl) :: current_max
    real(kind = dbl) :: suggested_min
    real(kind = dbl) :: suggested_max
    real(kind = dbl), allocatable :: values(:)
    logical :: upper_limited
    character(32) :: min_label
    character(32) :: max_label
    character(16) :: stability_label
    character(32) :: current_min_text
    character(32) :: current_max_text
    character(32) :: suggested_min_text
    character(32) :: suggested_max_text
    character(32), external :: IntToText
    real(kind = dbl), external :: SpectralFluxQuantile
    real(kind = dbl), external :: RoundFluxLimitDown
    real(kind = dbl), external :: RoundFluxLimitUp
    character(16), external :: GasName

    do stability = SADiagUnstable, SADiagStable
        n = SADiagFluxCandidateCount(stability, gas)
        if (stability == SADiagUnstable) then
            stability_label = 'Unstable'
        else
            stability_label = 'Stable'
        end if
        if (n == 0) then
            call EmitReportLine(report_unit, open_status, '  ' // trim(stability_label) // &
                ' flux-limit suggestion: NOT AVAILABLE - no records survived non-flux QA (u*, VM, or Foken).')
            cycle
        end if

        call SpectralFluxLimitSettings(gas, stability, current_min, current_max, min_label, max_label)
        allocate(values(n))
        values = SADiagFluxCandidate(1:n, stability, gas)
        call SortSpectralFluxValues(values, n)
        suggested_min = RoundFluxLimitDown(SpectralFluxQuantile(values, n, 0.10d0))
        upper_limited = any(values > current_max)
        suggested_max = current_max
        if (upper_limited) suggested_max = RoundFluxLimitUp(SpectralFluxQuantile(values, n, 0.99d0))

        class_counts = 0
        projected = 0
        do i = 1, n
            if (SADiagFluxCandidate(i, stability, gas) >= suggested_min .and. &
                SADiagFluxCandidate(i, stability, gas) <= suggested_max) then
                projected = projected + 1
                class_counts(SADiagFluxCandidateClass(i, stability, gas)) = &
                    class_counts(SADiagFluxCandidateClass(i, stability, gas)) + 1
            end if
        end do
        valid_classes = count(class_counts >= FCCsetup%SA%min_smpl)
        current_valid_classes = 0
        current_class_counts = 0
        do i = 1, n
            if (SADiagFluxCandidate(i, stability, gas) >= current_min .and. &
                SADiagFluxCandidate(i, stability, gas) <= current_max) then
                current_class_counts(SADiagFluxCandidateClass(i, stability, gas)) = &
                    current_class_counts(SADiagFluxCandidateClass(i, stability, gas)) + 1
            end if
        end do
        current_valid_classes = count(current_class_counts >= FCCsetup%SA%min_smpl)
        if (suggested_min < current_min .and. valid_classes >= 1) then
            SAAutoMin(stability, gas) = suggested_min
            SAAutoApplyMin(stability, gas) = .true.
        end if
        if (upper_limited .and. suggested_max > current_max .and. &
            valid_classes > current_valid_classes) then
            if (.not. SAAutoApplyMax(gas) .or. suggested_max > SAAutoMax(gas)) &
                SAAutoMax(gas) = suggested_max
            SAAutoApplyMax(gas) = .true.
        end if
        write(current_min_text, '(g0.6)') current_min
        write(current_max_text, '(g0.6)') current_max
        write(suggested_min_text, '(g0.6)') suggested_min
        write(suggested_max_text, '(g0.6)') suggested_max
        call EmitReportLine(report_unit, open_status, '  ' // trim(stability_label) // ' eligible flux records=' // &
            trim(IntToText(n)) // '; current ' // trim(min_label) // '=' // trim(current_min_text) // &
            ', ' // trim(max_label) // '=' // trim(current_max_text) // '.')
        call EmitReportLine(report_unit, open_status, '  Suggested ' // trim(min_label) // '=' // &
            trim(suggested_min_text) // ' (10th percentile; review before applying).')
        if (upper_limited) call EmitReportLine(report_unit, open_status, '  Suggested ' // trim(max_label) // '=' // &
            trim(suggested_max_text) // ' (99th percentile; current upper limit excludes eligible records).')
        !> Gathered for the end of the run as well, where they are the last
        !> thing printed instead of one block among thousands of period lines.
        call AddRunSuggestion('flux', trim(GasName(gas)) // ' (' // trim(stability_label) // &
            '): ' // trim(min_label) // '=' // trim(suggested_min_text) // ', now ' // &
            trim(current_min_text) // ' (10th percentile of eligible fluxes).')
        if (upper_limited) call AddRunSuggestion('flux', trim(GasName(gas)) // ' (' // &
            trim(stability_label) // '): ' // trim(max_label) // '=' // trim(suggested_max_text) // &
            ', now ' // trim(current_max_text) // ' (99th percentile; the current limit excludes eligible records).')
        call EmitReportLine(report_unit, open_status, '  Projected accepted records=' // trim(IntToText(projected)) // &
            ', valid classes=' // trim(IntToText(valid_classes)) // ' at sa_min_smpl=' // &
            trim(IntToText(FCCsetup%SA%min_smpl)) // '.')

        !> Two different failures wear the same report line, and the flux limit
        !> only answers one of them. A gas rejected for small fluxes is fixed by
        !> lowering the floor - that is the N2O and COS case, and the projection
        !> above shows classes appearing. A gas whose classes are simply too
        !> thinly populated is not: no floor produces a tenth sample in a class
        !> that only ever held four, and the suggested floor is derived as a
        !> percentile of what is already there, so applying it removes records.
        !> Saying which one this is saves the reader from trying the lever that
        !> cannot move.
        if (valid_classes == 0) then
            call EmitReportLine(report_unit, open_status, &
                '  Not a flux-limit problem: no class reaches sa_min_smpl=' // &
                trim(IntToText(FCCsetup%SA%min_smpl)) // ' at any flux limit. Best class holds ' // &
                trim(IntToText(maxval(class_counts))) // ' record(s).')
            call EmitReportLine(report_unit, open_status, &
                '  Widen the assessment period so the classes fill, or lower sa_min_smpl.')
        else if (suggested_min > current_min) then
            call EmitReportLine(report_unit, open_status, &
                '  Note: the suggested ' // trim(min_label) // ' is above the current one, so it')
            call EmitReportLine(report_unit, open_status, &
                '  discards records rather than admitting them. It raises the class count only')
            call EmitReportLine(report_unit, open_status, &
                '  because weak fluxes fit poorly; keep the current limit if records are scarce.')
        end if
        deallocate(values)
    end do
end subroutine ReportFluxLimitSuggestions

!*******************************************************************************
subroutine ApplyAutomaticSpectralConfiguration(output_project)
    use m_fx_global_var
    implicit none
    character(*), intent(in) :: output_project
    integer :: gas
    integer :: stability
    real(kind = dbl) :: current_min
    real(kind = dbl) :: current_max
    character(32) :: min_label
    character(32) :: max_label
    logical :: changes_written

    if (.not. FCCsetup%SA%automatic_config) return

    changes_written = any(SAAutoApplyMin) .or. any(SAAutoApplyMax) .or. any(SAAutoSource > 0)
    if (changes_written) then
        call LogSay(' Automatic spectral configuration: updating output project file.')
        do gas = firstGas, lastGas
            do stability = SADiagUnstable, SADiagStable
                if (.not. SAAutoApplyMin(stability, gas)) cycle
                call SpectralFluxLimitSettings(gas, stability, current_min, current_max, min_label, max_label)
                call WriteAutomaticSpectralSetting(output_project, trim(min_label), current_min, &
                    SAAutoMin(stability, gas), '10th-percentile eligible flux minimum')
            end do
            if (.not. SAAutoApplyMax(gas)) cycle
            call SpectralFluxLimitSettings(gas, SADiagUnstable, current_min, current_max, &
                min_label, max_label)
            call WriteAutomaticSpectralSetting(output_project, trim(max_label), current_max, &
                SAAutoMax(gas), &
                '99th-percentile eligible flux maximum')
        end do
        !> A gas left on automatic whose own transfer function looks
        !> implausible is given the tube-mate ReportSpectralSources suggested.
        !> Written as a record key, which an automatic gas does not have yet,
        !> so it has to be added to its section rather than rewritten.
        do gas = firstGas, lastGas
            if (SAAutoSource(gas) <= 0) cycle
            call WriteAutomaticSpectralSource(output_project, gas, SAAutoSource(gas))
        end do
        call EditIniFile(trim(output_project), 'automatic_spectra_config', '0')
        call AppendAutomaticSpectralConfigDiagnostic('Automatic spectral configuration: applied to output project: ' // &
            trim(output_project))
        call AppendAutomaticSpectralConfigDiagnostic( &
            'Automatic spectral configuration: current FCC results used the original settings; ' // &
            'rerun with this output project.')
        write(*, '(a)') ' Automatic spectral configuration saved to: ' // trim(output_project)
        write(ulog, '(a)') ' Automatic spectral configuration saved to: ' // trim(output_project)
        call LogSay(' Current FCC results used the original settings; rerun with this output project.')
    else
        call AppendAutomaticSpectralConfigDiagnostic( &
            'Automatic spectral configuration: no qualifying recommendations; output project unchanged.')
        call LogSay(' Automatic spectral configuration: no qualifying recommendations; output project unchanged.')
    end if
end subroutine ApplyAutomaticSpectralConfiguration

!*******************************************************************************
subroutine WriteAutomaticSpectralSetting(output_project, tag, original, replacement, reason)
    use m_fx_global_var
    implicit none
    character(*), intent(in) :: output_project
    character(*), intent(in) :: tag
    real(kind = dbl), intent(in) :: original
    real(kind = dbl), intent(in) :: replacement
    character(*), intent(in) :: reason
    character(32) :: original_text
    character(32) :: replacement_text
    character(512) :: line

    write(original_text, '(f0.6)') original
    write(replacement_text, '(f0.6)') replacement
    call EditIniFile(trim(output_project), trim(tag), trim(replacement_text))
    line = '  ' // trim(tag) // ': ' // trim(original_text) // ' -> ' // &
        trim(replacement_text) // ' (' // trim(reason) // ').'
    call AppendAutomaticSpectralConfigDiagnostic(trim(line))
    write(*, '(a)') trim(line)
    write(ulog, '(a)') trim(line)
end subroutine WriteAutomaticSpectralSetting

!*******************************************************************************
!> gas_<i>_sa_source=<j> in the output project: `gas` takes the transfer
!> function of `donor`, both gas slots.
!*******************************************************************************
subroutine WriteAutomaticSpectralSource(output_project, gas, donor)
    use m_fx_global_var
    implicit none
    character(*), intent(in) :: output_project
    integer, intent(in) :: gas
    integer, intent(in) :: donor
    character(32) :: tag
    character(16) :: value
    character(512) :: line
    character(32), external :: GasRecordLabel

    write(tag, '(a,i0,a)') 'gas_', gas - firstGas + 1, '_sa_source'
    write(value, '(i0)') donor - firstGas + 1
    call SetIniKeyInSection(trim(output_project), 'FluxCorrection_SpectralAnalysis_General', &
        trim(tag), trim(value))
    line = '  ' // trim(tag) // ': automatic -> ' // trim(value) // ' (' // &
        trim(GasRecordLabel(gas)) // ' takes ' // trim(GasRecordLabel(donor)) // &
        '''s transfer function; its own looks implausible).'
    call AppendAutomaticSpectralConfigDiagnostic(trim(line))
    write(*, '(a)') trim(line)
    write(ulog, '(a)') trim(line)
end subroutine WriteAutomaticSpectralSource

!*******************************************************************************
!> Write into the output project the transfer-function sources automatic
!> selection applied in this run (SAAppliedSource), as gas_<i>_sa_source.
!>
!> Unlike the suggestions ApplyAutomaticSpectralConfiguration writes, these
!> are not advice for the next run: they are what corrected this one, so they
!> are saved whatever automatic_spectra_config says. The project kept beside
!> the results then names every gas's source outright, and still reproduces
!> them if automatic selection is later switched off.
!*******************************************************************************
subroutine SaveAppliedSpectralSources(output_project)
    use m_fx_global_var
    implicit none
    character(*), intent(in) :: output_project
    integer :: gas

    if (.not. any(SAAppliedSource > 0)) return
    call LogSay(' Transfer-function sources chosen automatically in this run: ' // &
        'saved to the output project.')
    do gas = firstGas, lastGas
        if (SAAppliedSource(gas) <= 0) cycle
        call WriteAutomaticSpectralSource(output_project, gas, SAAppliedSource(gas))
    end do
end subroutine SaveAppliedSpectralSources

!*******************************************************************************
subroutine AppendAutomaticSpectralConfigDiagnostic(line)
    use m_fx_global_var
    implicit none
    character(*), intent(in) :: line
    integer :: report_unit
    integer :: open_status

    if (trim(SADiagFilePath) == 'none') return
    open(newunit = report_unit, file = trim(SADiagFilePath), status = 'old', &
        position = 'append', action = 'write', iostat = open_status)
    if (open_status /= 0) return
    write(report_unit, '(a)') trim(line)
    close(report_unit)
end subroutine AppendAutomaticSpectralConfigDiagnostic

!*******************************************************************************
subroutine SpectralFluxLimitSettings(gas, stability, minimum, maximum, min_label, max_label)
    use m_fx_global_var
    implicit none
    integer, intent(in) :: gas
    integer, intent(in) :: stability
    real(kind = dbl), intent(out) :: minimum
    real(kind = dbl), intent(out) :: maximum
    character(*), intent(out) :: min_label
    character(*), intent(out) :: max_label
    logical, external :: GasSlotIsWater

    !> Water is gated on the latent heat flux, which is a one-per-site
    !> quantity, so it keeps an arm of its own - resolved by species rather
    !> than read off the sixth slot.
    if (GasSlotIsWater(gas)) then
        if (stability == SADiagUnstable) then
            minimum = FCCsetup%SA%min_un_LE
            min_label = 'sa_min_un_le'
        else
            minimum = FCCsetup%SA%min_st_LE
            min_label = 'sa_min_st_le'
        end if
        maximum = FCCsetup%SA%max_LE
        max_label = 'sa_max_le'
        return
    end if

    !> Every other gas reads its own threshold. The four arms this replaces
    !> ended in a `case default` that read gas4's limits, so a fifth gas was
    !> filtered on the fourth slot's thresholds - even though ReadIniFCC had
    !> already loaded that gas's own value into min_un_gas(gas) from its
    !> record. The setting was written, read, and then not used.
    if (stability == SADiagUnstable) then
        minimum = FCCsetup%SA%min_un_gas(gas)
        min_label = SpectralLimitKey(gas, 'min_un')
    else
        minimum = FCCsetup%SA%min_st_gas(gas)
        min_label = SpectralLimitKey(gas, 'min_st')
    end if
    maximum = FCCsetup%SA%max_gas(gas)
    max_label = SpectralLimitKey(gas, 'max')

contains

!> The project-file key that carries this gas's limit.
!>
!> This must be the key ReadIniFCC actually reads, because the automatic
!> configuration does not merely name it in a report - it writes it, through
!> EditIniFile, into the project the next run reads.
!>
!> It used to compose `sa_min_un_co2` and friends for the first four slots and
!> `sa_min_un_gas_5_record` beyond them. Neither exists: the flat sa_min_*_<gas>
!> tags were retired with the record format, and the record spelling is
!> `gas_<i>_sa_min_un`. So automatic_spectra_config wrote settings nothing
!> reads, and reported changes that never took effect - for every gas,
!> including the first four.
character(32) function SpectralLimitKey(gas_slot, suffix)
    integer, intent(in) :: gas_slot
    character(*), intent(in) :: suffix

    write(SpectralLimitKey, '(a,i0,a)') &
        'gas_', gas_slot - firstGas + 1, '_sa_' // trim(suffix)
end function SpectralLimitKey

end subroutine SpectralFluxLimitSettings

!*******************************************************************************
real(kind = dbl) function SpectralFluxQuantile(values, n, probability)
    use m_common_global_var
    implicit none
    integer, intent(in) :: n
    real(kind = dbl), intent(in) :: values(n)
    real(kind = dbl), intent(in) :: probability
    integer :: rank

    rank = max(1, min(n, ceiling(probability * dfloat(n))))
    SpectralFluxQuantile = values(rank)
end function SpectralFluxQuantile

!*******************************************************************************
real(kind = dbl) function RoundFluxLimitDown(value)
    use m_common_global_var
    implicit none
    real(kind = dbl), intent(in) :: value
    real(kind = dbl) :: scale

    if (value <= 0d0) then
        RoundFluxLimitDown = value
    else
        scale = 10d0 ** (floor(log10(value)) - 2)
        RoundFluxLimitDown = floor(value / scale) * scale
    end if
end function RoundFluxLimitDown

!*******************************************************************************
real(kind = dbl) function RoundFluxLimitUp(value)
    use m_common_global_var
    implicit none
    real(kind = dbl), intent(in) :: value
    real(kind = dbl) :: scale

    if (value <= 0d0) then
        RoundFluxLimitUp = value
    else
        scale = 10d0 ** (floor(log10(value)) - 2)
        RoundFluxLimitUp = ceiling(value / scale) * scale
    end if
end function RoundFluxLimitUp

!*******************************************************************************
recursive subroutine SortSpectralFluxValues(values, n)
    use m_common_global_var
    implicit none
    integer, intent(in) :: n
    real(kind = dbl), intent(inout) :: values(n)
    integer :: i
    integer :: j
    real(kind = dbl) :: pivot
    real(kind = dbl) :: tmp

    if (n <= 1) return
    pivot = values((n + 1) / 2)
    i = 1
    j = n
    do
        do while (values(i) < pivot)
            i = i + 1
        end do
        do while (values(j) > pivot)
            j = j - 1
        end do
        if (i >= j) exit
        tmp = values(i)
        values(i) = values(j)
        values(j) = tmp
        i = i + 1
        j = j - 1
    end do
    call SortSpectralFluxValues(values(1:j), j)
    call SortSpectralFluxValues(values(j + 1:n), n - j)
end subroutine SortSpectralFluxValues

!*******************************************************************************
subroutine EmitReportLine(report_unit, open_status, line)
    use m_index_parameters
    use m_log
    implicit none
    integer, intent(in) :: report_unit
    integer, intent(in) :: open_status
    character(*), intent(in) :: line
    write(*, '(a)') '  ' // trim(line)
    write(ulog, '(a)') '  ' // trim(line)
    if (open_status == 0) write(report_unit, '(a)') trim(line)
end subroutine EmitReportLine

!*******************************************************************************
function IntToText(value) result(text)
    implicit none
    integer, intent(in) :: value
    character(32) :: text
    write(text, '(i0)') value
end function IntToText

!*******************************************************************************
function GasName(gas) result(name)
    use m_common_global_var
    implicit none
    integer, intent(in) :: gas
    character(16) :: name
    character(64) :: tags(GHGNumVar)
    include '../src_common/interfaces_1.inc'

    !> Every slot takes its record's species. Nothing here may assume a slot
    !> holds a given gas: slots are assigned by record order, so a project
    !> that declares its gases in a different order puts something other than
    !> CO2 in the first one. Pinning the first three named the wrong species,
    !> and calling everything past them "Gas 4" named five gases identically.
    call clearstr(name)
    if (gas >= firstGas .and. gas <= lastGas) then
        call SpectralGasNames(tags)
        if (len_trim(tags(gas)) > 0) then
            name = tags(gas)(1:min(len_trim(tags(gas)), len(name)))
            call uppercase(name)
        end if
    end if

    !> Only when the record names nothing at all. Positional, because at that
    !> point there is no species to report - but the number is the record's,
    !> not a fixed slot's.
    if (len_trim(name) == 0) write(name, '(a,i0)') 'Gas ', gas - firstGas + 1
end function GasName

!*******************************************************************************
function StatusLabel(passed, failed_label) result(label)
    implicit none
    logical, intent(in) :: passed
    character(*), intent(in) :: failed_label
    character(12) :: label
    if (passed) then
        label = 'PASS'
    else
        label = failed_label
    end if
end function StatusLabel

!*******************************************************************************
!> Which correction a gas ends up with, in the report's own words. Kept beside
!> StatusLabel so the two read as a pair: PASS/FAIL says whether the assessment
!> covered the gas, this says what that means for its flux.
!*******************************************************************************
character(12) function OutcomeLabel(in_situ)
    implicit none
    logical, intent(in) :: in_situ

    if (in_situ) then
        OutcomeLabel = 'in situ'
    else
        OutcomeLabel = 'analytically'
    end if
end function OutcomeLabel

!*******************************************************************************
!> Whether a non-water gas is corrected in situ in any month, and through
!> whose fit: `donor` is the first tube-mate SpectralSourceSlot hands it in a
!> month that has a cut-off, 0 when its own fit serves throughout.
!>
!> The readiness report's "corrected in situ" has to follow the source the
!> correction reads. Asked of the gas's own fit, a COS whose fit failed but
!> which takes CO2's would be reported as corrected analytically, and one
!> whose named donor has nothing would be reported as in situ.
!*******************************************************************************
logical function GasInsituWithSource(gas, donor)
    use m_fx_global_var
    implicit none
    integer, intent(in) :: gas
    integer, intent(out) :: donor
    integer :: month
    integer :: src
    integer :: cls
    integer, external :: SpectralSourceSlot

    GasInsituWithSource = .false.
    donor = 0
    do month = 1, 12
        src = SpectralSourceSlot(gas, month, FCCsetup%SA)
        cls = FCCsetup%SA%class(src, month)
        if (cls < 1 .or. cls > MaxGasClasses) cycle
        if (RegPar(src, cls)%fc == error) cycle
        GasInsituWithSource = .true.
        if (src /= gas .and. donor == 0) donor = src
    end do
end function GasInsituWithSource

!*******************************************************************************
!> " with CO2_2 (miro_mga1_5_2)'s transfer function", or nothing for 0.
!*******************************************************************************
character(64) function DonorNote(donor)
    use m_fx_global_var
    implicit none
    integer, intent(in) :: donor
    character(32), external :: GasRecordLabel

    DonorNote = ''
    if (donor <= 0) return
    DonorNote = ' with ' // trim(GasRecordLabel(donor)) // '''s transfer function'
end function DonorNote

!*******************************************************************************
!> Keep a suggestion for the end of the run. `section` is 'sources' or 'flux'.
!> Past MaxRunSuggestions the rest are dropped - they were printed where they
!> were made, so nothing is lost but the repeat.
!*******************************************************************************
subroutine AddRunSuggestion(section, line)
    use m_fx_global_var
    implicit none
    character(*), intent(in) :: section
    character(*), intent(in) :: line

    if (RunSuggestionN >= MaxRunSuggestions) return
    RunSuggestionN = RunSuggestionN + 1
    RunSuggestionSection(RunSuggestionN) = section
    RunSuggestions(RunSuggestionN) = line
end subroutine AddRunSuggestion

!*******************************************************************************
!> Every suggestion the run made, as its very last message before the banner.
!> Nothing when there is none.
!*******************************************************************************
subroutine ReportRunSuggestions()
    use m_fx_global_var
    implicit none
    integer :: i
    integer :: k
    character(8), parameter :: Sections(2) = [character(8) :: 'sources', 'flux']
    character(96), parameter :: Headings(2) = [character(96) :: &
        ' Spectral transfer functions that look implausible (gas_<i>_sa_source):', &
        ' Flux limits for the spectral assessment:']

    if (RunSuggestionN < 1) return
    call LogSay('')
    call LogSay(' Suggested settings (review before applying)')
    call LogSay(' -------------------------------------------')
    do k = 1, size(Sections)
        if (.not. any(RunSuggestionSection(1:RunSuggestionN) == Sections(k))) cycle
        call LogSay(trim(Headings(k)))
        do i = 1, RunSuggestionN
            if (RunSuggestionSection(i) /= Sections(k)) cycle
            call LogSay('   ' // trim(RunSuggestions(i)))
        end do
    end do
end subroutine ReportRunSuggestions

!*******************************************************************************
!> Which gases are corrected with a tube-mate's transfer function, and which
!> would be: one line per gas and run of months with the same answer, logged
!> now and - where it is a suggestion - kept for the end of the run.
!>
!> Called once by the parent, after the assessment is read or fitted and
!> before the context is handed to workers. Workers re-evaluate the same rule
!> per period from the same RegPar, so they need no report of their own.
!>
!> Also records, for ApplyAutomaticSpectralConfiguration, the donor a gas left
!> on automatic would take: the one suggested for the most months.
!*******************************************************************************
subroutine ReportSpectralSources()
    use m_fx_global_var
    implicit none
    integer :: gas
    integer :: month
    integer :: m0
    integer :: used
    integer :: sugg
    integer :: cur_used
    integer :: cur_sugg
    integer :: first_used
    integer :: counts(GHGNumVar)
    logical :: header
    logical :: all_same
    logical, external :: GasSlotIsWater
    integer, external :: SpectralSourceSlot
    integer, external :: SpectralSuggestedSource

    SAAutoSource = 0
    SAAppliedSource = 0
    if (.not. FCCsetup%SA%in_situ) return
    header = .false.
    do gas = firstGas, lastGas
        if (gas - firstGas + 1 > min(EddyFlowProj%gas_num, MaxNumGases)) exit
        if (GasSlotIsWater(gas)) cycle
        if (EddyFlowProj%gas(gas - firstGas + 1)%col <= 0) cycle
        counts = 0
        m0 = 1
        cur_used = SpectralSourceSlot(gas, 1, FCCsetup%SA)
        cur_sugg = SpectralSuggestedSource(gas, 1, FCCsetup%SA)
        if (cur_sugg > 0) counts(cur_sugg) = counts(cur_sugg) + 1
        first_used = cur_used
        all_same = .true.
        do month = 2, 12
            used = SpectralSourceSlot(gas, month, FCCsetup%SA)
            sugg = SpectralSuggestedSource(gas, month, FCCsetup%SA)
            if (sugg > 0) counts(sugg) = counts(sugg) + 1
            if (used /= first_used) all_same = .false.
            if (used /= cur_used .or. sugg /= cur_sugg) then
                call SourceRange(gas, m0, month - 1, cur_used, cur_sugg)
                m0 = month
                cur_used = used
                cur_sugg = sugg
            end if
        end do
        call SourceRange(gas, m0, 12, cur_used, cur_sugg)
        if (FCCsetup%SA%source(gas) == SASourceAuto .and. .not. FCCsetup%SA%auto_source &
            .and. maxval(counts) > 0) SAAutoSource(gas) = maxloc(counts, 1)
        !> What automatic selection did, when it did the same thing all year:
        !> one gas_<i>_sa_source states it exactly. A gas switched for some
        !> months only stays automatic in the output project - the key has no
        !> months - and the copy's sa_auto_source=1 repeats the same choice.
        if (FCCsetup%SA%source(gas) == SASourceAuto .and. FCCsetup%SA%auto_source &
            .and. all_same .and. first_used /= gas) SAAppliedSource(gas) = first_used
    end do

contains

    !> One line for months m1..m2 of `gas`, which take `used`'s fit and for
    !> which `sugg` is suggested (0: none). Silent when the gas takes its own
    !> fit and nothing better is suggested.
    subroutine SourceRange(gas, m1, m2, used, sugg)
        integer, intent(in) :: gas
        integer, intent(in) :: m1
        integer, intent(in) :: m2
        integer, intent(in) :: used
        integer, intent(in) :: sugg
        integer :: other
        character(256) :: line
        character(48) :: status
        character(32), external :: GasRecordLabel

        if (used == gas .and. sugg == 0) return
        if (used /= gas) then
            other = used
        else
            other = sugg
        end if

        if (used /= gas .and. FCCsetup%SA%source(gas) > 0) then
            status = 'applied: set by user'
        else if (used /= gas) then
            status = 'applied: automatic'
        else if (FCCsetup%SA%source(gas) == SASourceOwn) then
            status = 'not applied: "Own assessment" set by user'
        else
            status = 'not applied: automatic selection is off'
        end if

        line = trim(GasRecordLabel(gas)) // ' ' // trim(MonthSpan(m1, m2)) // &
            ': own ' // trim(FitText(gas, m1)) // ' -> ' // &
            trim(GasRecordLabel(other)) // ': ' // trim(FitText(other, m1)) // &
            '  [' // trim(status) // ']'

        if (.not. header) then
            call LogSay('')
            call LogSay(' Spectral transfer functions taken from tube-mates on the same analyser:')
            header = .true.
        end if
        call LogSay('  ' // trim(line))
        if (sugg > 0) call AddRunSuggestion('sources', trim(line))
    end subroutine SourceRange

    !> "Fn 0.036, fc 0.00199 Hz" for a gas's own class of `month`
    character(64) function FitText(slot, month)
        integer, intent(in) :: slot
        integer, intent(in) :: month
        integer :: cls

        cls = FCCsetup%SA%class(slot, month)
        FitText = 'no fit'
        if (cls < 1 .or. cls > MaxGasClasses) return
        if (RegPar(slot, cls)%fc == error) return
        FitText = 'Fn ' // trim(ShortReal(RegPar(slot, cls)%Fn)) // ', fc ' // &
            trim(ShortReal(RegPar(slot, cls)%fc)) // ' Hz'
    end function FitText

    !> Three significant figures or so, with the leading zero gfortran's F0.d
    !> leaves out
    character(24) function ShortReal(x)
        real(kind = dbl), intent(in) :: x
        character(24) :: t

        if (x == error) then
            ShortReal = 'n/a'
            return
        end if
        if (dabs(x) >= 0.1d0) then
            write(t, '(f0.3)') x
        else if (dabs(x) >= 0.01d0) then
            write(t, '(f0.4)') x
        else
            write(t, '(f0.5)') x
        end if
        t = adjustl(t)
        if (t(1:1) == '.') t = '0' // trim(t)
        if (t(1:2) == '-.') t = '-0' // trim(t(2:))
        ShortReal = t
    end function ShortReal

    !> "Jan-Dec", "Mar", "Apr-Sep"
    character(8) function MonthSpan(m1, m2)
        integer, intent(in) :: m1
        integer, intent(in) :: m2
        character(3), parameter :: Names(12) = [character(3) :: 'Jan', 'Feb', &
            'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

        if (m1 == m2) then
            MonthSpan = Names(m1)
        else
            MonthSpan = Names(m1) // '-' // Names(m2)
        end if
    end function MonthSpan

end subroutine ReportSpectralSources
