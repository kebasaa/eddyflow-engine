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

    if (gas < co2 .or. gas > gas4) return
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
    integer :: report_unit
    integer :: mkdir_status
    integer :: open_status
    character(PathLen) :: SpecDir
    character(PathLen) :: FilePath
    character(128) :: Filename
    character(32) :: method
    logical :: gas_ok
    logical :: h2o_ok
    character(32), external :: IntToText
    character(16), external :: GasName
    character(12), external :: StatusLabel
    include '../src_common/interfaces_1.inc'

    assessment_ready = .true.
    method = trim(adjustl(EddyFlowProj%hf_meth))
    h2o_classes = 0
    h2o_fits = 0
    if (allocated(MeanBinSpec)) then
        do cls = RH10, RH90
            if (MeanBinSpecAvailable(cls, h2o)) h2o_classes = h2o_classes + 1
            if (RegPar(h2o, cls)%fc /= error) h2o_fits = h2o_fits + 1
        end do
    end if

    if (FCCsetup%do_spectral_assessment) then
        h2o_ok = h2o_classes >= 1 .and. &
            RegPar(dum, dum)%e1 /= error .and. RegPar(dum, dum)%e2 /= error .and. &
            RegPar(dum, dum)%e3 /= error
        assessment_ready = h2o_ok
        do gas = co2, gas4
            if (gas == h2o .or. .not. fcc_var_present(gas)) cycle
            gas_ok = .false.
            do cls = 1, MaxGasClasses
                if (RegPar(gas, cls)%fc /= error) gas_ok = .true.
            end do
            assessment_ready = assessment_ready .and. gas_ok
        end do
    end if

    mkdir_status = CreateDir('"' // Dir%main_out(1:len_trim(Dir%main_out)) // '"')
    SpecDir = Dir%main_out(1:len_trim(Dir%main_out)) // SubDirSpecAn // slash
    mkdir_status = CreateDir('"' // SpecDir(1:len_trim(SpecDir)) // '"')
    Filename = EddyFlowProj%id(1:len_trim(EddyFlowProj%id)) // &
        '_spectral_correction_diagnostics' // Timestamp_FilePadding // TxtExt
    FilePath = SpecDir(1:len_trim(SpecDir)) // Filename(1:len_trim(Filename))
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

    do gas = co2, gas4
        gas_classes = 0
        if (allocated(MeanBinSpec)) then
            do cls = 1, MaxGasClasses
                if (MeanBinSpecAvailable(cls, gas)) gas_classes = gas_classes + 1
            end do
        end if
        if (.not. fcc_var_present(gas)) then
            call EmitReportLine(report_unit, open_status, trim(GasName(gas)) // ': unavailable... NOT APPLIED')
        else
            call EmitReportLine(report_unit, open_status, trim(GasName(gas)) // ': accepted periods=' // &
                trim(IntToText(SADiagAccepted(gas))) // ', flux=' // trim(IntToText(SADiagRejectedFlux(gas))) // &
                ', VM=' // trim(IntToText(SADiagRejectedVM(gas))) // ', Foken=' // &
                trim(IntToText(SADiagRejectedFoken(gas))) // ', valid classes=' // &
                trim(IntToText(gas_classes)) // '... ' // &
                trim(StatusLabel(gas_classes > 0, 'FAIL')))
            if (gas_classes == 0 .and. SADiagRejectedFlux(gas) > 0) &
                call ReportFluxLimitSuggestions(report_unit, open_status, gas)
        end if
    end do
    call EmitReportLine(report_unit, open_status, 'H2O RH classes: ' // trim(IntToText(h2o_classes)) // &
        '/9, fitted cut-offs=' // trim(IntToText(h2o_fits)) // ', minimum samples/class=' // &
        trim(IntToText(FCCsetup%SA%min_smpl)) // '... ' // &
        trim(StatusLabel(h2o_classes >= 1, 'FAIL')))
    if (h2o_classes > 0 .and. h2o_classes < 3) then
        call EmitReportLine(report_unit, open_status, &
            'H2O RH coverage: WARNING - one class can create a file; three or more are recommended.')
    end if
    call EmitReportLine(report_unit, open_status, 'Valid degraded wT covariance: unstable=' // &
        trim(IntToText(SADiagDegradedUnstable)) // ', stable=' // trim(IntToText(SADiagDegradedStable)) // '... ' // &
        trim(StatusLabel(SADiagDegradedUnstable + SADiagDegradedStable > 0, 'WARNING')))

    if (FCCsetup%do_spectral_assessment) then
        if (assessment_ready) then
            call EmitReportLine(report_unit, open_status, 'Spectral assessment: SUCCESS - assessment file will be created.')
        else
            call EmitReportLine(report_unit, open_status, 'Spectral assessment: FAILED - assessment file will not be created.')
            call EmitReportLine(report_unit, open_status, &
                'Standard EddyFlow errors and Moncrieff fallback will follow during output/correction.')
            call EmitReportLine(report_unit, open_status, &
                'Adjust the reported QA thresholds/filters or select more qualifying periods.')
        end if
    end if
    if (open_status == 0) close(report_unit)
    write(*, '(a)') ' Spectral-correction diagnostics written to: ' // trim(FilePath)
end subroutine ReportSpectralAssessmentDiagnostics

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
    integer :: class_counts(MaxGasClasses)
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
        call EmitReportLine(report_unit, open_status, '  Projected accepted records=' // trim(IntToText(projected)) // &
            ', valid classes=' // trim(IntToText(valid_classes)) // ' at sa_min_smpl=' // &
            trim(IntToText(FCCsetup%SA%min_smpl)) // '.')
        deallocate(values)
    end do
end subroutine ReportFluxLimitSuggestions

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

    select case (gas)
        case (h2o)
            if (stability == SADiagUnstable) then
                minimum = FCCsetup%SA%min_un_LE
                min_label = 'sa_min_un_le'
            else
                minimum = FCCsetup%SA%min_st_LE
                min_label = 'sa_min_st_le'
            end if
            maximum = FCCsetup%SA%max_LE
            max_label = 'sa_max_le'
        case (co2)
            if (stability == SADiagUnstable) then
                minimum = FCCsetup%SA%min_un_co2
                min_label = 'sa_min_un_co2'
            else
                minimum = FCCsetup%SA%min_st_co2
                min_label = 'sa_min_st_co2'
            end if
            maximum = FCCsetup%SA%max_co2
            max_label = 'sa_max_co2'
        case (ch4)
            if (stability == SADiagUnstable) then
                minimum = FCCsetup%SA%min_un_ch4
                min_label = 'sa_min_un_ch4'
            else
                minimum = FCCsetup%SA%min_st_ch4
                min_label = 'sa_min_st_ch4'
            end if
            maximum = FCCsetup%SA%max_ch4
            max_label = 'sa_max_ch4'
        case default
            if (stability == SADiagUnstable) then
                minimum = FCCsetup%SA%min_un_gas4
                min_label = 'sa_min_un_gas4'
            else
                minimum = FCCsetup%SA%min_st_gas4
                min_label = 'sa_min_st_gas4'
            end if
            maximum = FCCsetup%SA%max_gas4
            max_label = 'sa_max_gas4'
    end select
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
    implicit none
    integer, intent(in) :: report_unit
    integer, intent(in) :: open_status
    character(*), intent(in) :: line
    write(*, '(a)') '  ' // trim(line)
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
    select case (gas)
        case (co2)
            name = 'CO2'
        case (h2o)
            name = 'H2O'
        case (ch4)
            name = 'CH4'
        case default
            name = 'Gas 4'
    end select
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
