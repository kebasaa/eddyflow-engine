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
end subroutine ResetSpectralAssessmentDiagnostics

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
