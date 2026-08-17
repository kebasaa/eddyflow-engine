!***************************************************************************
! init_ex_vars.f90
! ----------------
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
! \brief       Reads essentials file, retrieving all information that might \n
!              be useful to other programs
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine InitExVars(StartTimestamp, EndTimestamp, NumRecords, NumValidRecords, FirstValidRecord)
    use m_fx_global_var
    implicit none
    !> In/out variables
    integer, intent(out) :: NumRecords
    integer, intent(out) :: NumValidRecords
    integer, intent(out) :: FirstValidRecord
    type(DateType), intent(out) :: StartTimestamp
    type(DateType), intent(out) :: EndTimestamp
    !> local variables
    integer :: open_status
    integer :: j
    integer :: k
    integer :: gas
    integer :: field_start
    integer :: field_end
    integer :: field_count
    integer :: marker_custom
    integer :: marker_biomet
    character(128) :: custom_field
    character(64) :: custom_label
    logical :: label_has_alpha
    logical :: ValidRecord
    logical :: EndOfFileReached
    logical :: InitializationPerformed
    type (ExType) :: lEX
    include '../src_common/interfaces_1.inc'

    write(*,'(a)') &
        ' Initializing retrieval of EddyFlow-RP results from file: '
    write(ulog,'(a)') &
        ' Initializing retrieval of EddyFlow-RP results from file: '
    write(*,'(a)') '  "' // trim(adjustl(AuxFile%ex)) // '"..'
    write(ulog,'(a)') '  "' // trim(adjustl(AuxFile%ex)) // '"..'

    !> Open EX file
    open(udf, file = AuxFile%ex, status = 'old', iostat = open_status)

    !> Exit with error in case of problems opening the file
    if (open_status /= 0) call ExceptionHandler(60)

    call LogSay('  File found, importing content..')

    !> Store header to string, for writing it on output
    read(udf, '(a)') fluxnet_header

    !> The fourth gas's label used to be recovered here by finding ',FCH4,'
    !> in the header and taking the column after it, then stripping one
    !> character to drop the 'F' flux prefix. That assumed a gas named
    !> methane exists, that the fourth gas is written immediately after it,
    !> and that its flux tag is one character longer than its name - none of
    !> which holds once records assign species to slots by declaration order.
    !> On a project without methane `index` returned 0 and the label became a
    !> slice of the first columns of the header. Every consumer now names its
    !> gas from the project's own records, via SpectralGasNames.


    UserVarHeader = ''
    marker_custom = index(fluxnet_header, 'NUM_CUSTOM_VARS')
    marker_biomet = index(fluxnet_header, 'NUM_BIOMET_VARS')
    if (marker_custom > 0 .and. marker_biomet > marker_custom) then
        field_start = marker_custom + len('NUM_CUSTOM_VARS') + 1
        field_count = 0
        do while (field_start < marker_biomet .and. field_count < MaxUserVar)
            field_end = field_start + index(fluxnet_header(field_start:), ',') - 2
            if (field_end < field_start) exit
            field_count = field_count + 1
            call clearstr(custom_field)
            call clearstr(custom_label)
            custom_field = fluxnet_header(field_start:field_end)
            custom_label = replace2(custom_field, 'CUSTOM_', '')
            call lowercase(custom_label)
            label_has_alpha = index(custom_label, 'a') > 0 .or. index(custom_label, 'b') > 0 &
                .or. index(custom_label, 'c') > 0 .or. index(custom_label, 'd') > 0 &
                .or. index(custom_label, 'e') > 0 .or. index(custom_label, 'f') > 0 &
                .or. index(custom_label, 'g') > 0 .or. index(custom_label, 'h') > 0 &
                .or. index(custom_label, 'i') > 0 .or. index(custom_label, 'j') > 0 &
                .or. index(custom_label, 'k') > 0 .or. index(custom_label, 'l') > 0 &
                .or. index(custom_label, 'm') > 0 .or. index(custom_label, 'n') > 0 &
                .or. index(custom_label, 'o') > 0 .or. index(custom_label, 'p') > 0 &
                .or. index(custom_label, 'q') > 0 .or. index(custom_label, 'r') > 0 &
                .or. index(custom_label, 's') > 0 .or. index(custom_label, 't') > 0 &
                .or. index(custom_label, 'u') > 0 .or. index(custom_label, 'v') > 0 &
                .or. index(custom_label, 'w') > 0 .or. index(custom_label, 'x') > 0 &
                .or. index(custom_label, 'y') > 0 .or. index(custom_label, 'z') > 0
            if (label_has_alpha) then
                if (len_trim(custom_label) > 5 &
                    .and. custom_label(len_trim(custom_label) - 4:len_trim(custom_label)) == '_mean') then
                    custom_label = custom_label(1:len_trim(custom_label) - 5)
                end if
                if (len_trim(custom_label) <= len(custom_label) - 5) &
                    custom_label = custom_label(1:len_trim(custom_label)) // '_mean'
                UserVarHeader(field_count) = custom_label
            end if
            field_start = field_end + 2
        end do
    end if

    !> Initialize variables that are determined for the whole
    !> dataset (presence of certain variables)
    Diag7200%present = .false.
    Diag7500%present = .false.
    Diag7700%present = .false.
    fcc_var_present = .false.
    FCCMetadata%ru = .false.
    FCCMetadata%ac_freq = -1
    FCCMetadata%GasAcFreq = error
    DateStep = DateType(0, 0, 0, 0, ierror)

    !> Cycle on all records
    NumRecords = 0
    NumValidRecords = 0
    InitializationPerformed = .false.

    do
        !> Read essentials record
        call ReadExRecord('', udf, -1, lEx, ValidRecord, EndOfFileReached)
        if (EndOfFileReached) exit

        !> Counts
        NumRecords = NumRecords + 1

        if (NumValidRecords == 0 .and. ValidRecord) FirstValidRecord = NumRecords

        if (ValidRecord) NumValidRecords = NumValidRecords + 1

        !> Handles dates
        if (ValidRecord .and. NumValidRecords == 1) &
            call DateTimeToDateType(lEx%end_date, lEX%end_time, StartTimestamp)
        if (ValidRecord) &
            call DateTimeToDateType(lEx%end_date, lEX%end_time, EndTimestamp)

        !> Initializations
        if (ValidRecord .and. .not. InitializationPerformed) then

            !> Look for variable presence, over every gas slot the project
            !> configures rather than the first four. This gate is what every
            !> FCC output loop tests, so leaving it four-bounded made gases 5+
            !> absent from the full output no matter how wide the loops were.
            if (lEx%WS /= error) fcc_var_present(u:w) = .true.
            if (lEx%Ts /= error) fcc_var_present(ts)  = .true.
            !> A gas is present because the project names a column for it, not
            !> because this particular record carried a value. Asking the
            !> essentials record instead - its measure type is the error code
            !> whenever the gas was filtered away - made a gas that lost its
            !> data lose its columns as well, so COS and a second analyser's
            !> CO2 and H2O were absent from the full output altogether rather
            !> than present and empty. A filtered gas is written as the error
            !> label; a missing column is not the way to say "no data".
            !>
            !> Only as far as the project configures, and only records that
            !> name a column: an unconfigured slot would otherwise emit a
            !> column family per empty slot, and a record with no column - a
            !> species the site does not measure - has nothing to report here.
            do k = 1, min(EddyFlowProj%gas_num, MaxNumGases)
                gas = firstGas + k - 1
                if (gas > lastGas) exit
                if (EddyFlowProj%gas(k)%col > 0) fcc_var_present(gas) = .true.
            end do
                
            !> Determine whether LI-COR's flags are available
            if (.not. Diag7200%present) then
                do j = 1, 9
                    if (lEx%licor_flags(j) /= error) then
                        Diag7200%present = .true.
                        exit
                    end if
                end do
            end if

            if (.not. Diag7500%present) then
                do j = 10, 13
                    if (lEx%licor_flags(j) /= error) then
                        Diag7500%present = .true.
                        exit
                    end if
                end do
            end if

            if (.not. Diag7700%present) then
                do j = 14, 29
                    if (lEx%licor_flags(j) /= error) then
                        Diag7700%present = .true.
                        exit
                    end if
                end do
            end if

            !> Reads DateStep
            if (DateStep == DateType(0, 0, 0, 0, ierror)) DateStep = DateType(0, 0, 0, 0, nint(lEx%avrg_length))

            !> Define whether random uncertainty was calculated by
            !> looking at only 1 value (if one value is -6999d0, all
            !> of them are the same)
            if (lEx%rand_uncer(u) == aflx_error) FCCMetadata%ru = .false.

            !> Acquisition frequency and gas analyser path type for H2O
            if (FCCMetadata%ac_freq <= 0) FCCMetadata%ac_freq = lEx%ac_freq
            !> From the site's water record. Was lEx%instr(ih2o), the water
            !> role of the retired five-wide instrument numbering.
            FCCMetadata%H2oPathType = &
                lEx%gas_instr(PrimaryWaterOutSlot())%path_type
            !> And every slot's own, for the routines that must not assume the
            !> primary's - FitRh2Fco fits one hygrometer per slot now.
            do gas = firstGas, lastGas
                FCCMetadata%GasPathType(gas) = &
                    lEx%gas_instr(gas)%path_type
                !> And each analyser's own rate, for the checks that must not
                !> apply the station's Nyquist to a slower instrument.
                FCCMetadata%GasAcFreq(gas) = &
                    lEx%gas_instr(gas)%ac_freq
            end do
        end if

        if (all(fcc_var_present) .and. Diag7200%present .and. Diag7500%present .and. Diag7700%present .and. &
           FCCMetadata%ac_freq > 0 .and. FCCMetadata%ru .and. DateStep /= DateType(0, 0, 0, 0, ierror)) then
            InitializationPerformed = .true.
        end if
    end do
    close(udf)

    !> Adjust start timestamp so that Start/End define the whole period
    !> From beginning of first period to end of last period
    StartTimestamp = StartTimestamp - DateStep
    call LogSay(' Done.')
end subroutine InitExVars
