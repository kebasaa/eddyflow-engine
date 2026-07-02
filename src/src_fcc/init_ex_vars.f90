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
    integer :: st
    integer :: en
    integer :: j
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
    write(*,'(a)') '  "' // trim(adjustl(AuxFile%ex)) // '"..'

    !> Open EX file
    open(udf, file = AuxFile%ex, status = 'old', iostat = open_status)

    !> Exit with error in case of problems opening the file
    if (open_status /= 0) call ExceptionHandler(60)

    write(*, '(a)') '  File found, importing content..'

    !> Store header to string, for writing it on output
    read(udf, '(a)') fluxnet_header

    st = index(fluxnet_header, ',FCH4,') + 6
    en = st + index(fluxnet_header(st:), ',') - 2
    g4lab = fluxnet_header(st+1:en)
    call lowercase(g4lab)
    g4l = len_trim(g4lab)
    

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

            !> Look for variable presence (u thru GS4)
            if (lEx%WS /= error) fcc_var_present(u:w) = .true.
            if (lEx%Ts /= error) fcc_var_present(ts)  = .true.
            do gas = co2, gas4
                fcc_var_present(gas) = lEx%measure_type_int(gas) /= ierror .or. fcc_var_present(gas)  
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
            FCCMetadata%H2oPathType = lEx%instr(ih2o)%path_type
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
    write(*,'(a)') ' Done.'
end subroutine InitExVars
