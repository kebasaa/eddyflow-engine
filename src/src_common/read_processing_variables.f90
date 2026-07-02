!***************************************************************************
! read_processing_variables.f90
! -----------------------------
! Reads the [ProcessingVariables] project-file group.
!***************************************************************************
subroutine ReadProcessingVariables(IniFile)
    use m_common_global_var
    implicit none
    character(*), intent(in) :: IniFile
    integer :: uini, io_status, eq_pos, count_value, file_row, enabled_count, n_tags
    logical :: in_group, found_group
    character(ShortInstringLen) :: dataline
    type(Text) :: tags(MaxNLinesIni)
    integer :: GetProcessingInt
    logical :: GetProcessingLogical
    character(iniLabelLen) :: RowKey

    call ResetProcessingVariables()
    uini = udf
    open(uini, file = IniFile, status = 'old', iostat = io_status)
    if (io_status /= 0) then
        call ExceptionHandler(7)
        return
    end if

    in_group = .false.
    found_group = .false.
    n_tags = 0
    do
        read(uini, '(a)', iostat = io_status) dataline
        if (io_status /= 0) exit
        call stripstr(dataline)
        if (len_trim(dataline) == 0) cycle
        if (dataline(1:1) == ';') cycle
        if (dataline(1:1) == '[') then
            if (trim(adjustl(dataline)) == '[ProcessingVariables]') then
                in_group = .true.
                found_group = .true.
                cycle
            elseif (in_group) then
                exit
            else
                cycle
            end if
        end if
        if (.not. in_group) cycle
        eq_pos = index(dataline, '=')
        if (eq_pos <= 1) cycle
        n_tags = n_tags + 1
        if (n_tags > MaxNLinesIni) call ProcessingVariablesError('Too many ProcessingVariables tags')
        tags(n_tags)%Label = dataline(1:eq_pos - 1)
        tags(n_tags)%Value = dataline(eq_pos + 1:len_trim(dataline))
    end do
    close(uini)

    if (.not. found_group) then
        call ProcessingVariablesError('[ProcessingVariables] group is required')
        return
    end if

    EddyFlowProj%processing%has_processing_variables = .true.
    count_value = GetProcessingInt(tags, n_tags, 'count', 0)
    if (count_value <= 0) call ProcessingVariablesError('ProcessingVariables count must be positive')
    if (count_value > MaxProcessingVariables) &
        call ProcessingVariablesError('ProcessingVariables count exceeds MaxProcessingVariables')
    EddyFlowProj%processing%file_count = count_value

    enabled_count = 0
    do file_row = 1, count_value
        if (GetProcessingLogical(tags, n_tags, RowKey(file_row, 'enabled'), .true.)) then
            enabled_count = enabled_count + 1
            call ReadProcessingVariableRow(tags, n_tags, file_row, EddyFlowProj%processing%rows(enabled_count))
        end if
    end do
    EddyFlowProj%processing%count = enabled_count
    call ValidateProcessingVariables()
    call LogProcessingVariables()
end subroutine ReadProcessingVariables

subroutine ResetProcessingVariables()
    use m_common_global_var
    implicit none
    integer :: i
    EddyFlowProj%processing = GasCollectionType()
    do i = 1, MaxProcessingVariables
        EddyFlowProj%processing%rows(i) = ProcessingVariableType()
        EddyFlowProj%processing%results(i) = GasResultType()
        EddyFlowProj%processing%stats(i) = GasStatsType()
    end do
end subroutine ResetProcessingVariables

subroutine ReadProcessingVariableRow(tags, n_tags, file_row, row)
    use m_common_global_var
    implicit none
    integer, intent(in) :: n_tags, file_row
    type(Text), intent(in) :: tags(MaxNLinesIni)
    type(ProcessingVariableType), intent(out) :: row
    integer :: GetProcessingInt
    real(kind = dbl) :: GetProcessingReal
    character(iniLabelLen) :: RowKey
    character(iniValueLen) :: GetProcessingText
    character(64) :: ExpectedProcessingId
    character(64) :: configured_id, expected_id

    row = ProcessingVariableType()
    row%enabled = .true.
    row%file_row = file_row
    row%gas_col = GetProcessingInt(tags, n_tags, RowKey(file_row, 'gas_col'), -1)
    row%gas_name = trim(adjustl(GetProcessingText(tags, n_tags, RowKey(file_row, 'gas'), '')))
    call lowercase(row%gas_name)
    row%irga_id = trim(adjustl(GetProcessingText(tags, n_tags, RowKey(file_row, 'irga'), '')))
    row%irga_index = GetProcessingInt(tags, n_tags, RowKey(file_row, 'irga_index'), 0)
    row%gas_instance_index = GetProcessingInt(tags, n_tags, RowKey(file_row, 'gas_index'), 0)
    row%molecular_weight = GetProcessingReal(tags, n_tags, RowKey(file_row, 'mw'), error)
    row%molecular_diffusivity = GetProcessingReal(tags, n_tags, RowKey(file_row, 'diff'), error)
    if (row%molecular_weight > 0d0) row%molecular_weight = row%molecular_weight * 1d-3
    if (row%molecular_diffusivity > 0d0) row%molecular_diffusivity = row%molecular_diffusivity * 1d-4
    row%reference_h2o_id = trim(adjustl(GetProcessingText(tags, n_tags, RowKey(file_row, 'h2o_ref'), '')))
    row%col_cell_t = GetProcessingInt(tags, n_tags, RowKey(file_row, 'cell_t'), -1)
    row%col_int_t_1 = GetProcessingInt(tags, n_tags, RowKey(file_row, 'int_t_1'), -1)
    row%col_int_t_2 = GetProcessingInt(tags, n_tags, RowKey(file_row, 'int_t_2'), -1)
    row%col_int_p = GetProcessingInt(tags, n_tags, RowKey(file_row, 'int_p'), -1)
    row%col_air_t = GetProcessingInt(tags, n_tags, RowKey(file_row, 'air_t'), -1)
    row%col_air_p = GetProcessingInt(tags, n_tags, RowKey(file_row, 'air_p'), -1)
    row%col_diag = GetProcessingInt(tags, n_tags, RowKey(file_row, 'diag'), -1)
    configured_id = trim(adjustl(GetProcessingText(tags, n_tags, RowKey(file_row, 'id'), '')))
    expected_id = ExpectedProcessingId(row%gas_name, row%irga_index, row%gas_instance_index)
    if (len_trim(configured_id) == 0) then
        row%processing_id = expected_id
    else
        row%processing_id = configured_id
    end if
end subroutine ReadProcessingVariableRow

subroutine ValidateProcessingVariables()
    use m_common_global_var
    implicit none
    integer :: i, j, h2o_count, only_h2o, ref_index, FindProcessingVariableById
    character(64) :: ExpectedProcessingId
    character(64) :: expected_id
    h2o_count = 0
    only_h2o = 0
    if (EddyFlowProj%processing%count <= 0) &
        call ProcessingVariablesError('ProcessingVariables must contain at least one enabled row')
    do i = 1, EddyFlowProj%processing%count
        if (len_trim(EddyFlowProj%processing%rows(i)%processing_id) == 0) &
            call ProcessingVariablesError('Enabled ProcessingVariables row has empty id')
        if (len_trim(EddyFlowProj%processing%rows(i)%gas_name) == 0) &
            call ProcessingVariablesError('Enabled ProcessingVariables row has empty gas')
        if (EddyFlowProj%processing%rows(i)%gas_col <= 0) &
            call ProcessingVariablesError('Enabled ProcessingVariables row has non-positive gas_col')
        if (EddyFlowProj%processing%rows(i)%irga_index <= 0) &
            call ProcessingVariablesError('Enabled ProcessingVariables row has non-positive irga_index')
        if (EddyFlowProj%processing%rows(i)%gas_instance_index <= 0) &
            call ProcessingVariablesError('Enabled ProcessingVariables row has non-positive gas_index')
        expected_id = ExpectedProcessingId(EddyFlowProj%processing%rows(i)%gas_name, &
            EddyFlowProj%processing%rows(i)%irga_index, &
            EddyFlowProj%processing%rows(i)%gas_instance_index)
        if (trim(EddyFlowProj%processing%rows(i)%processing_id) /= trim(expected_id)) &
            call ProcessingVariablesError('ProcessingVariables id must match gas_irgaIndex_gasIndex: ' // &
                trim(EddyFlowProj%processing%rows(i)%processing_id) // ' /= ' // trim(expected_id))
        do j = i + 1, EddyFlowProj%processing%count
            if (trim(EddyFlowProj%processing%rows(i)%processing_id) == trim(EddyFlowProj%processing%rows(j)%processing_id)) &
                call ProcessingVariablesError('Duplicate ProcessingVariables id: ' // &
                    trim(EddyFlowProj%processing%rows(i)%processing_id))
        end do
        if (trim(EddyFlowProj%processing%rows(i)%gas_name) == 'h2o') then
            h2o_count = h2o_count + 1
            only_h2o = i
        end if
    end do
    do i = 1, EddyFlowProj%processing%count
        if (trim(EddyFlowProj%processing%rows(i)%gas_name) == 'h2o') then
            if (len_trim(EddyFlowProj%processing%rows(i)%reference_h2o_id) == 0) &
                EddyFlowProj%processing%rows(i)%reference_h2o_id = EddyFlowProj%processing%rows(i)%processing_id
            if (trim(EddyFlowProj%processing%rows(i)%reference_h2o_id) /= &
                trim(EddyFlowProj%processing%rows(i)%processing_id)) &
                call ProcessingVariablesError('H2O row may only self-reference: ' // &
                    trim(EddyFlowProj%processing%rows(i)%processing_id))
            EddyFlowProj%processing%rows(i)%h2o_ref_index = i
            cycle
        end if
        if (len_trim(EddyFlowProj%processing%rows(i)%reference_h2o_id) == 0) then
            if (h2o_count == 1) then
                EddyFlowProj%processing%rows(i)%reference_h2o_id = EddyFlowProj%processing%rows(only_h2o)%processing_id
            else
                call ProcessingVariablesError('Ambiguous or missing H2O reference for ' // &
                    trim(EddyFlowProj%processing%rows(i)%processing_id))
            end if
        end if
        ref_index = FindProcessingVariableById(EddyFlowProj%processing%rows(i)%reference_h2o_id)
        if (ref_index <= 0) call ProcessingVariablesError('Invalid H2O reference for ' // &
            trim(EddyFlowProj%processing%rows(i)%processing_id))
        if (trim(EddyFlowProj%processing%rows(ref_index)%gas_name) /= 'h2o') &
            call ProcessingVariablesError('H2O reference is not an H2O row')
        EddyFlowProj%processing%rows(i)%h2o_ref_index = ref_index
    end do
end subroutine ValidateProcessingVariables

subroutine LogProcessingVariables()
    use m_common_global_var
    implicit none
    integer :: i
    if (EddyFlowProj%processing%count <= 0) return
    write(*,'(a,i0)') ' Processing variable rows enabled: ', EddyFlowProj%processing%count
    do i = 1, EddyFlowProj%processing%count
        write(*,'(a,a,a,a,a,i0,a,a,a,f8.4,a,f8.5)') '  ', trim(EddyFlowProj%processing%rows(i)%processing_id), &
            ' gas=', trim(EddyFlowProj%processing%rows(i)%gas_name), ' col=', EddyFlowProj%processing%rows(i)%gas_col, &
            ' h2o_ref=', trim(EddyFlowProj%processing%rows(i)%reference_h2o_id), &
            ' mw=', EddyFlowProj%processing%rows(i)%molecular_weight, &
            ' diff=', EddyFlowProj%processing%rows(i)%molecular_diffusivity
    end do
end subroutine LogProcessingVariables

integer function FindProcessingVariableById(processing_id)
    use m_common_global_var
    implicit none
    character(*), intent(in) :: processing_id
    integer :: i
    FindProcessingVariableById = 0
    do i = 1, EddyFlowProj%processing%count
        if (trim(EddyFlowProj%processing%rows(i)%processing_id) == trim(processing_id)) then
            FindProcessingVariableById = i
            return
        end if
    end do
end function FindProcessingVariableById

character(64) function ExpectedProcessingId(gas_name, irga_index, gas_index)
    implicit none
    character(*), intent(in) :: gas_name
    integer, intent(in) :: irga_index, gas_index
    character(32) :: normalized_gas
    normalized_gas = trim(adjustl(gas_name))
    call lowercase(normalized_gas)
    write(ExpectedProcessingId, '(a,a,i0,a,i0)') trim(normalized_gas), '_', irga_index, '_', gas_index
end function ExpectedProcessingId

character(iniLabelLen) function RowKey(row_number, suffix)
    use m_typedef
    implicit none
    integer, intent(in) :: row_number
    character(*), intent(in) :: suffix
    write(RowKey, '(a,i0,a,a)') 'row_', row_number, '_', trim(suffix)
end function RowKey

character(iniValueLen) function GetProcessingText(tags, n_tags, key, default_value)
    use m_common_global_var
    implicit none
    integer, intent(in) :: n_tags
    type(Text), intent(in) :: tags(MaxNLinesIni)
    character(*), intent(in) :: key, default_value
    integer :: i
    GetProcessingText = default_value
    do i = 1, n_tags
        if (trim(tags(i)%Label) == trim(key)) then
            GetProcessingText = tags(i)%Value
            return
        end if
    end do
end function GetProcessingText

integer function GetProcessingInt(tags, n_tags, key, default_value)
    use m_common_global_var
    implicit none
    integer, intent(in) :: n_tags, default_value
    type(Text), intent(in) :: tags(MaxNLinesIni)
    character(*), intent(in) :: key
    integer :: read_status
    character(iniValueLen) :: raw, GetProcessingText
    raw = GetProcessingText(tags, n_tags, key, '')
    if (len_trim(raw) == 0) then
        GetProcessingInt = default_value
    else
        read(raw, *, iostat = read_status) GetProcessingInt
        if (read_status /= 0) GetProcessingInt = default_value
    end if
end function GetProcessingInt

real(kind = dbl) function GetProcessingReal(tags, n_tags, key, default_value)
    use m_common_global_var
    implicit none
    integer, intent(in) :: n_tags
    real(kind = dbl), intent(in) :: default_value
    type(Text), intent(in) :: tags(MaxNLinesIni)
    character(*), intent(in) :: key
    integer :: read_status
    character(iniValueLen) :: raw, GetProcessingText
    raw = GetProcessingText(tags, n_tags, key, '')
    if (len_trim(raw) == 0) then
        GetProcessingReal = default_value
    else
        read(raw, *, iostat = read_status) GetProcessingReal
        if (read_status /= 0) GetProcessingReal = default_value
    end if
end function GetProcessingReal

logical function GetProcessingLogical(tags, n_tags, key, default_value)
    use m_common_global_var
    implicit none
    integer, intent(in) :: n_tags
    logical, intent(in) :: default_value
    type(Text), intent(in) :: tags(MaxNLinesIni)
    character(*), intent(in) :: key
    character(iniValueLen) :: raw, GetProcessingText
    raw = trim(adjustl(GetProcessingText(tags, n_tags, key, '')))
    if (len_trim(raw) == 0) then
        GetProcessingLogical = default_value
    else
        GetProcessingLogical = raw(1:1) == '1' .or. raw(1:1) == 't' .or. raw(1:1) == 'T'
    end if
end function GetProcessingLogical

subroutine ProcessingVariablesError(message)
    implicit none
    character(*), intent(in) :: message
    write(*,'(a)') ' EddyFlow project configuration error: ' // trim(message)
    stop 1
end subroutine ProcessingVariablesError

