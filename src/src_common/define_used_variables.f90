!***************************************************************************
! define_used_variables.f90
! -------------------------
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
! \brief       Determine whether to use individual data column or not, based either on the
!              fact that it's coming from a master sonic, or on the value assigned
!              to %useit (for variables other than sonics).
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine DefineUsedVariables(LocCol)
    use m_common_global_var
    implicit none
    !> in/out variables
    type(ColType), intent(inout) :: LocCol(MaxNumCol)
    !> local variables
    integer :: i
    integer :: slot
    integer :: selected_ts_col
    logical :: ts_found
    integer, external :: GasSlotFromDynMDTag


    NumUserVar = 0

    !> Associate the "flag" properties to columns selected as such by user
    do i = 1, NumRawFlags
        LocCol(RawFlag(i)%col)%flag = RawFlag(i)
    end do

    !> Associate the "master_sonic" property to the relevant instrument
    !> This automatically associate the master_sonic property to the
    !> relevant sonic variables
    LocCol%Instr%master_sonic = .false.
    do i = 1, NumCol
        if (index(EddyFlowProj%master_sonic, &
            trim(adjustl(LocCol(i)%Instr%model))) /= 0) then
            LocCol(i)%Instr%master_sonic = .true.
        end if
    end do

    LocCol%useit = .false.
    !> Information in EddyFlow project file (user explicitly selects which
    !> variables are to be used)
    where (EddyFlowProj%Col(firstGas:E2NumVar) > 0)
        LocCol(EddyFlowProj%Col(firstGas:E2NumVar))%useit = .true.
    endwhere

    where (EddyFlowProj%Col(E2NumVar + diag72 :E2NumVar + diagAnem) > 0)
        LocCol(EddyFlowProj%Col(E2NumVar + diag72 :E2NumVar + diagAnem))%useit = .true.
    endwhere

    !> Columns named by gas/cell/diagnostic records must be marked here too.
    !> This runs while LocCol is still indexed by .metadata column number, and
    !> unused columns are dropped straight afterwards - a record column left
    !> unmarked never reaches DefineE2Set, so the record silently selects
    !> nothing.
    do i = 1, min(EddyFlowProj%gas_num, MaxNumGases)
        if (EddyFlowProj%gas(i)%col > 0 .and. &
            EddyFlowProj%gas(i)%col <= MaxNumCol) &
            LocCol(EddyFlowProj%gas(i)%col)%useit = .true.
    end do
    do i = 1, min(EddyFlowProj%cell_num, MaxNumCellCols)
        if (EddyFlowProj%cell(i)%col > 0 .and. &
            EddyFlowProj%cell(i)%col <= MaxNumCol) &
            LocCol(EddyFlowProj%cell(i)%col)%useit = .true.
    end do
    do i = 1, min(EddyFlowProj%diag_num, MaxNumDiagCols)
        if (EddyFlowProj%diag(i)%col > 0 .and. &
            EddyFlowProj%diag(i)%col <= MaxNumCol) &
            LocCol(EddyFlowProj%diag(i)%col)%useit = .true.
    end do

    !> The fourth gas's column used to be renamed to 'n2o' here so the rest of
    !> the engine would treat it as that species. It was gated on col_gas4,
    !> which is retired, so the rename could never fire - and a record names
    !> its own species, which is what ApplyGasRecords resolves. Renaming a
    !> column to n2o regardless of what it measured is the assumption this
    !> whole effort removes.

    !> Diagnostic flags
    NumDiag = 0
    Diag7200%present = .false.
    Diag7500%present = .false.
    Diag7700%present = .false.
    DiagAnemometer%binary_flag_present = .false.
    DiagAnemometer%staa_present = .false.
    if (EddyFlowProj%Col(E2NumVar + diag72) > 0) then
        LocCol(EddyFlowProj%Col(E2NumVar + diag72))%useit = .true.
        NumDiag = NumDiag + 1
        Diag7200%present = .true.
    end if
    if (EddyFlowProj%Col(E2NumVar + diag75) > 0) then
        LocCol(EddyFlowProj%Col(E2NumVar + diag75))%useit = .true.
        NumDiag = NumDiag + 1
        Diag7500%present = .true.
    end if
    if (EddyFlowProj%Col(E2NumVar + diag77) > 0) then
        LocCol(EddyFlowProj%Col(E2NumVar + diag77))%useit = .true.
        NumDiag = NumDiag + 1
        Diag7700%present = .true.
    end if
    if (EddyFlowProj%Col(E2NumVar + diagAnem) > 0) then
        LocCol(EddyFlowProj%Col(E2NumVar + diagAnem))%useit = .true.
        NumDiag = NumDiag + 1
        DiagAnemometer%binary_flag_present = .true.
    end if
    if (EddyFlowProj%Col(E2NumVar + diagStaA) > 0) then
        LocCol(EddyFlowProj%Col(E2NumVar + diagStaA))%useit = .true.
        NumDiag = NumDiag + 1
        DiagAnemometer%staa_present = .true.
    end if

    !> Loop on the actual number of columns and determine
    !> whether to use them or not
    GasCalRefCol = 0
    do i = 1, NumCol
        !> Variables from the master_sonic are to be used
        if (LocCol(i)%instr%master_sonic) then
            LocCol(i)%useit = .true.
            cycle
        end if
        !> Count users variables, made up of: sonic variables from a non-master
        !> sonic; irga variables without the property "use_it", and
        !> those with a custom label
        if (IsCustomOutputColumn(LocCol(i)) .and. NumUserVar < MaxUserVar - 1) &
            NumUserVar = NumUserVar + 1

        !> Detect whether a gas calibration data column is available.
        !>
        !> `<gas>_cal-ref` names the gas it calibrates, resolved the same way
        !> the drift subsystem resolves `<gas>_ref`. A bare `cal-ref` with no
        !> prefix keeps calibrating the fourth slot, which is what every
        !> metadata file written before this says and means.
        if (index(LocCol(i)%var, 'cal-ref') /= 0) then
            slot = GasSlotFromDynMDTag(LocCol(i)%var, '_cal-ref')
            if (slot <= 0) slot = gas4
            GasCalRefCol(slot) = i
        end if
    end do

    !> If user selects a different temperature reading
    !> (instead of sonic temperature) for sensible heat flux, redefine Ts
    !> as that column, but changes instrument category to "fast_t_sensor"
    !> to remember that it does not need water vapor correction and
    !> (sonic-specific) spectral corrections.
    if (EddyFlowProj%Col(ts) > 0) then
        selected_ts_col = EddyFlowProj%Col(ts)
        LocCol(selected_ts_col)%useit = .true.
        LocCol(selected_ts_col)%var = 'ts'
        LocCol(selected_ts_col)%instr%category = 'fast_t_sensor'
        !> Search Ts or SoS from master sonic and change property in
        !> "don't use it", so now it will fall into the "non sensitive"
        !> variables group. Note that the total number of User Variables
        !> did not change
        ts_found = .false.
        do i = 1, NumCol
            if (i /= selected_ts_col &
                .and. LocCol(i)%instr%master_sonic &
                .and. (LocCol(i)%var == 'ts' .or. LocCol(i)%var == 'sos' ) &
                .and. LocCol(i)%useit) then
                LocCol(i)%useit = .false.
                ts_found = .true.
                !exit
            end if
        end do
        !> If Ts was not there instead, the number of user variables must
        !> be reduced by one, because one was used
        !> as a fast temperature
        if (.not. ts_found) NumUserVar = NumUserVar - 1
    end if

contains

logical function IsCustomOutputColumn(col)
    type(ColType), intent(in) :: col
    character(32) :: var

    IsCustomOutputColumn = .false.
    if (col%useit) return

    var = col%var
    call lowercase(var)
    if (len_trim(var) == 0) return
    select case (trim(var))
        case ('ignore', 'not_numeric', 'none', 'flag_1', 'flag_2', &
              'agc', 'rssi')
            return
        case default
            IsCustomOutputColumn = .true.
    end select
end function IsCustomOutputColumn
end subroutine DefineUsedVariables
