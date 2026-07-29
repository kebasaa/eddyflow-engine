!***************************************************************************
! define_e2_set.f90
! -----------------
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
! \brief       Define "E2Set", the pre-defined set of variables needed for any following \n
!              processing. Variables are: u, v, w, ts, co2, h2o, ch4, gas4, tc, ti1, ti2, pi, te, pe
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine DefineE2Set(LocCol, Raw, nrow, ncol, E2Set, e2nrow, e2ncol, DiagSet, dnrow, dncol)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    integer, intent(in) :: e2nrow, e2ncol
    integer, intent(in) :: dnrow, dncol
    type(ColType), intent(in) :: LocCol(MaxNumCol)
    real(kind = sgl), intent(in)  :: Raw(nrow, ncol)
    real(kind = dbl), intent(out) :: E2Set(e2nrow, e2ncol)
    real(kind = dbl), intent(out) :: DiagSet(dnrow, dncol)
    !> local variables
    integer :: j
    !> Cell-record instruments seen so far, in first-seen order. Host
    !> variables so that both ApplyCellDiagRecords and CellInstrIndex reach
    !> them; an internal procedure cannot see a sibling's locals.
    integer :: nCellInstr
    integer :: i
    character(32) :: cellInstr(MaxNumInstruments)


    E2Col = NullCol
    E2Set = error
    DiagSet = error
    !> Initialised here, not inside ApplyCellDiagRecords: a project with no
    !> cell records never enters that routine, and the resolution loop below
    !> would then read an uninitialised count.
    nCellInstr = 0
    cellInstr = ''
    !> First, master sonic wind components are handled
    do j = 1, ncol
        !> u-component of wind vector
        if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) == 'u' .and. LocCol(j)%Instr%master_sonic) then
            E2Col(u) = LocCol(j)
            E2Col(u)%present = .true.
            E2Set(1:e2nrow, u) = Raw(1:e2nrow, j)
            cycle
        end if
        !> v-component of wind vector
        if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) == 'v' .and. LocCol(j)%Instr%master_sonic) then
            E2Col(v) = LocCol(j)
            E2Col(v)%present = .true.
            E2Set(1:e2nrow, v) = Raw(1:e2nrow, j)
            cycle
        end if
        !> w-component of wind vector
        if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) == 'w' .and. LocCol(j)%Instr%master_sonic) then
            E2Col(w) = LocCol(j)
            E2Col(w)%present = .true.
            E2Set(1:e2nrow, w) = Raw(1:e2nrow, j)
            cycle
        end if
        !> sonic/fast temperature
        if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) == 'ts' .and. LocCol(j)%useit) then
            E2Col(ts) = LocCol(j)
            E2Col(ts)%present = .true.
            E2Set(1:e2nrow, ts) = Raw(1:e2nrow, j)
            cycle
        end if
        !> speed of sound (already converted into sonic temperature)
        if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) == 'sos' .and. LocCol(j)%Instr%master_sonic) then
            E2Col(ts) = LocCol(j)
            E2Col(ts)%present = .true.
            E2Set(1:e2nrow, ts) = Raw(1:e2nrow, j)
            cycle
        end if
    end do

    !> Now, all remaining EddyFlow standard variables (concentrations, temperatures and pressure)
    do j = 1, ncol
        !> master co2 concentration
        if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) == 'co2' .and. LocCol(j)%useit) then
            E2Col(co2) = LocCol(j)
            E2Col(co2)%present = .true.
            E2Set(1:e2nrow, co2) = Raw(1:e2nrow, j)
            cycle
        end if
        !> master h2o concentration
        if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) == 'h2o' .and. LocCol(j)%useit) then
            E2Col(h2o) = LocCol(j)
            E2Col(h2o)%present = .true.
            E2Set(1:e2nrow, h2o) = Raw(1:e2nrow, j)
            cycle
        end if
        !> master ch4 concentration
        if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) == 'ch4' .and. LocCol(j)%useit) then
            E2Col(ch4) = LocCol(j)
            E2Col(ch4)%present = .true.
            E2Set(1:e2nrow, ch4) = Raw(1:e2nrow, j)
            cycle
        end if
        !> master 4th gas concentration
        if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) == 'n2o' .and. LocCol(j)%useit) then
            E2Col(gas4) = LocCol(j)
            E2Col(gas4)%present = .true.
            E2Set(1:e2nrow, gas4) = Raw(1:e2nrow, j)
            cycle
        end if
        !> master cell temperature
        if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) == 'cell_t' .and. LocCol(j)%useit) then
            E2Col(tc) = LocCol(j)
            E2Col(tc)%present = .true.
            E2Set(1:e2nrow, tc) = Raw(1:e2nrow, j)
            cycle
        end if
        !> master internal temperature 1
        if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) == 'int_t_1' .and. LocCol(j)%useit) then
            E2Col(ti1) = LocCol(j)
            E2Col(ti1)%present = .true.
            E2Set(1:e2nrow, ti1) = Raw(1:e2nrow, j)
            cycle
        end if
        !> master internal temperature 2
        if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) == 'int_t_2' .and. LocCol(j)%useit) then
            E2Col(ti2) = LocCol(j)
            E2Col(ti2)%present = .true.
            E2Set(1:e2nrow, ti2) = Raw(1:e2nrow, j)
            cycle
        end if
        !> master internal pressure
        if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) == 'int_p' .and. LocCol(j)%useit) then
            E2Col(pi) = LocCol(j)
            E2Col(pi)%present = .true.
            E2Set(1:e2nrow, pi) = Raw(1:e2nrow, j)
            cycle
        end if
        !> master air temperature
        if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) == 'air_t' .and. LocCol(j)%useit) then
            E2Col(te) = LocCol(j)
            E2Col(te)%present = .true.
            E2Set(1:e2nrow, te) = Raw(1:e2nrow, j)
            cycle
        end if
        !> master air pressure
        if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) == 'air_p' .and. LocCol(j)%useit) then
            E2Col(pe) = LocCol(j)
            E2Col(pe)%present = .true.
            E2Set(1:e2nrow, pe) = Raw(1:e2nrow, j)
            cycle
        end if
    end do

    do j = 1, ncol
        !> Diagnostic flags set
        if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) == 'diag_72' .and. LocCol(j)%useit) then
            DiagSet(1:dnrow, diag72) = Raw(1:dnrow, j)
            cycle
        end if
        if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) == 'diag_75' .and. LocCol(j)%useit) then
            DiagSet(1:dnrow, diag75) = Raw(1:dnrow, j)
            cycle
        end if
        if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) == 'diag_77' .and. LocCol(j)%useit) then
            DiagSet(1:dnrow, diag77) = Raw(1:dnrow, j)
            cycle
        end if
        if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) == 'anemometer_diagnostic' .and. LocCol(j)%useit) then
            DiagSet(1:dnrow, diagAnem) = Raw(1:dnrow, j)
            cycle
        end if
        if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) == 'Gill_StaA' .and. LocCol(j)%useit) then
            DiagSet(1:dnrow, diagStaA) = Raw(1:dnrow, j)
            cycle
        end if
        if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) == 'Gill_StaD' .and. LocCol(j)%useit) then
            DiagSet(1:dnrow, diagStaD) = Raw(1:dnrow, j)
            cycle
        end if
    end do

    !> If the project file describes gas records, they take over the gas slots.
    !> The name matching above can only fill the four legacy slots, one column
    !> per species; records address a slot per measurement, so the same species
    !> can appear more than once. With no records nothing below runs and the
    !> legacy selection stands.
    call ApplyGasRecords(LocCol, Raw, nrow, ncol, E2Set, e2nrow, e2ncol)
    call ApplyCellDiagRecords(LocCol, Raw, nrow, ncol, E2Set, e2nrow, e2ncol, &
                              DiagSet, dnrow, dncol)

    !> Point every gas at the cell block of the analyser that measured it.
    !>
    !> cell_ref holds the *base* slot of that instrument's block, so the four
    !> quantities are cell_ref + 0..3 in the order cell_t, int_t_1, int_t_2,
    !> int_p. Gases whose instrument has no cell record fall back to the first
    !> block, which is where a single-analyser project's data has always been.
    do j = firstGas, lastGas
        if (.not. E2Col(j)%present) cycle
        E2Col(j)%cell_ref = firstCell
        if (nCellInstr <= 0) cycle
        do i = 1, nCellInstr
            if (trim(cellInstr(i)) == trim(E2Col(j)%instr%model) .or. &
                trim(cellInstr(i)) == trim(E2Col(j)%Instr%ep_label)) then
                E2Col(j)%cell_ref = firstCell + (i - 1) * NumCellPerInstr
                exit
            end if
        end do
    end do

    !> Default every unresolved gas to the H2O in the fixed slot. The legacy
    !> selection has exactly one H2O, so this is what the single-analyser case
    !> has always used implicitly; making it explicit means both paths describe
    !> their moisture the same way instead of the legacy one leaving it blank.
    if (E2Col(h2o)%present) then
        do j = firstGas, lastGas
            if (.not. E2Col(j)%present) cycle
            if (E2Col(j)%moist_ref < firstGas .or. E2Col(j)%moist_ref > lastGas) &
                E2Col(j)%moist_ref = h2o
        end do
    end if

contains

!***************************************************************************
!> Position of `instr` in the cell-instrument list, adding it if new.
!>
!> Returns 0 once the list is full, dropping the extra records rather than
!> letting them wrap onto another instrument's slots.
integer function CellInstrIndex(instr) result(idx)
    implicit none
    character(*), intent(in) :: instr
    integer :: j

    do j = 1, nCellInstr
        if (trim(cellInstr(j)) == trim(instr)) then
            idx = j
            return
        end if
    end do
    if (nCellInstr >= MaxNumInstruments) then
        idx = 0
        return
    end if
    nCellInstr = nCellInstr + 1
    cellInstr(nCellInstr) = instr
    idx = nCellInstr
end function CellInstrIndex

!***************************************************************************
!> Fill the cell temperature/pressure and diagnostic slots from records.
!>
!> LIMITATION: E2Col has one slot each for cell_t / int_t_1 / int_t_2 / int_p,
!> so only one instrument's cell measurements can be held at a time - a second
!> record naming the same quantity overwrites the first. Genuine per-instrument
!> cell T/P needs the cell slots widened the way the gas slots were, and until
!> then GasRecordType%cell stays unresolved (cell_ref = 0) and the historical
!> global slots apply.
!***************************************************************************
subroutine ApplyCellDiagRecords(LocCol, Raw, nrow, ncol, E2Set, e2nrow, e2ncol, &
                                DiagSet, dnrow, dncol)
    implicit none
    integer, intent(in) :: nrow, ncol, e2nrow, e2ncol, dnrow, dncol
    type(ColType), intent(in) :: LocCol(MaxNumCol)
    real(kind = sgl), intent(in) :: Raw(nrow, ncol)
    real(kind = dbl), intent(inout) :: E2Set(e2nrow, e2ncol)
    real(kind = dbl), intent(inout) :: DiagSet(dnrow, dncol)
    integer :: i, slot, src, offset, k

    if (EddyFlowProj%cell_num > 0) then
        do slot = firstCell, lastCell
            E2Col(slot) = NullCol
            E2Set(1:e2nrow, slot) = error
        end do
        do i = 1, min(EddyFlowProj%cell_num, MaxNumCellCols)
            src = LocColByOrigCol(LocCol, ncol, EddyFlowProj%cell(i)%col)
            if (src <= 0) cycle
            select case (trim(EddyFlowProj%cell(i)%var))
                case ('cell_t');  offset = 0
                case ('int_t_1'); offset = 1
                case ('int_t_2'); offset = 2
                case ('int_p');   offset = 3
                case default;     cycle
            end select
            !> Each instrument gets its own set of cell slots. Keyed on the
            !> record's instrument rather than on arrival order, so the two
            !> halves of one analyser's cell data land together; records used
            !> to share a single set, and a second analyser silently
            !> overwrote the first.
            k = CellInstrIndex(trim(EddyFlowProj%cell(i)%instr))
            if (k <= 0) cycle
            slot = firstCell + (k - 1) * NumCellPerInstr + offset
            if (slot > lastCell) cycle
            E2Col(slot) = LocCol(src)
            E2Col(slot)%present = .true.
            E2Set(1:e2nrow, slot) = Raw(1:e2nrow, src)
        end do
    end if

    if (EddyFlowProj%diag_num > 0) then
        DiagSet(1:dnrow, diag72)   = error
        DiagSet(1:dnrow, diag75)   = error
        DiagSet(1:dnrow, diag77)   = error
        DiagSet(1:dnrow, diagAnem) = error
        do i = 1, min(EddyFlowProj%diag_num, MaxNumDiagCols)
            src = LocColByOrigCol(LocCol, ncol, EddyFlowProj%diag(i)%col)
            if (src <= 0) cycle
            select case (trim(EddyFlowProj%diag(i)%var))
                case ('diag_72'); slot = diag72
                case ('diag_75'); slot = diag75
                case ('diag_77'); slot = diag77
                case ('anemometer_diagnostic'); slot = diagAnem
                case default;     cycle
            end select
            DiagSet(1:dnrow, slot) = Raw(1:dnrow, src)
        end do
    end if
end subroutine ApplyCellDiagRecords

!***************************************************************************
!> Fill E2Col/E2Set gas slots from EddyFlowProj%gas, then resolve each gas's
!> moisture reference to an E2Col slot.
!***************************************************************************
subroutine ApplyGasRecords(LocCol, Raw, nrow, ncol, E2Set, e2nrow, e2ncol)
    implicit none
    integer, intent(in) :: nrow, ncol, e2nrow, e2ncol
    type(ColType), intent(in) :: LocCol(MaxNumCol)
    real(kind = sgl), intent(in) :: Raw(nrow, ncol)
    real(kind = dbl), intent(inout) :: E2Set(e2nrow, e2ncol)
    integer :: i, slot, src

    !> Records are the only way a project names its gases. A file without
    !> them is a pre-5.0.0 project that has not been through the interface;
    !> processing would silently produce no gas fluxes at all, so it is
    !> refused instead. The GUI migrates such a file on open.
    if (EddyFlowProj%gas_num <= 0) then
        write(*, '(a)') '  Fatal error(99)> Project file describes no gas &
            &records (gas_num). Open and save it in the EddyFlow interface &
            &to bring it up to the current format.'
        call ExceptionHandler(99)
    end if

    !> A record set replaces the gas selection wholesale, so clear the slots
    !> the name matching filled before re-populating from records.
    do slot = firstGas, lastGas
        E2Col(slot) = NullCol
        E2Set(1:e2nrow, slot) = error
    end do

    do i = 1, min(EddyFlowProj%gas_num, MaxNumGases)
        src = LocColByOrigCol(LocCol, ncol, EddyFlowProj%gas(i)%col)
        if (src <= 0) cycle
        slot = firstGas + i - 1
        E2Col(slot) = LocCol(src)
        E2Col(slot)%present = .true.
        E2Set(1:e2nrow, slot) = Raw(1:e2nrow, src)
    end do

    !> Resolve references only once every slot is populated, so that "same
    !> instrument" is evaluated against what is actually present.
    do i = 1, min(EddyFlowProj%gas_num, MaxNumGases)
        slot = firstGas + i - 1
        if (.not. E2Col(slot)%present) cycle
        E2Col(slot)%moist_ref = ResolveGasRef(i, EddyFlowProj%gas(i)%moist, 'h2o')
        E2Col(slot)%cell_ref = 0
    end do
end subroutine ApplyGasRecords

!***************************************************************************
!> Translate a .metadata column number into its position in LocCol/Raw.
!>
!> Columns the project does not use are dropped before this point, so the two
!> numbering schemes diverge. The project file stores metadata numbers, so a
!> record must be matched on %orig_col rather than indexed directly - doing
!> the latter silently selects a different variable.
!***************************************************************************
integer function LocColByOrigCol(LocCol, ncol, wanted)
    implicit none
    type(ColType), intent(in) :: LocCol(MaxNumCol)
    integer, intent(in) :: ncol
    integer, intent(in) :: wanted
    integer :: k

    LocColByOrigCol = 0
    if (wanted < 1) return
    do k = 1, ncol
        if (LocCol(k)%orig_col == wanted) then
            LocColByOrigCol = k
            return
        end if
    end do
end function LocColByOrigCol

!***************************************************************************
!> Map a record-level reference onto an E2Col slot.
!>
!> An explicit reference (>0) is a 1-based index into the gas records. Zero
!> means "auto": prefer a record of the wanted species on the same instrument,
!> falling back to the first of that species anywhere. The GUI shows the same
!> default, so the two rules must not drift - change them together.
!***************************************************************************
integer function ResolveGasRef(gasIdx, ref, wantedVar)
    implicit none
    integer, intent(in) :: gasIdx
    integer, intent(in) :: ref
    character(*), intent(in) :: wantedVar
    integer :: k, slot

    ResolveGasRef = 0

    if (ref > 0 .and. ref <= EddyFlowProj%gas_num) then
        slot = firstGas + ref - 1
        if (E2Col(slot)%present) ResolveGasRef = slot
        return
    end if

    !> 1. same instrument. 'other' and 'none' are not identities: many
    !> unrelated variables carry them, so matching on those would pair gases
    !> with arbitrary partners.
    if (trim(EddyFlowProj%gas(gasIdx)%instr) /= 'other' .and. &
        trim(EddyFlowProj%gas(gasIdx)%instr) /= 'none') then
        do k = 1, min(EddyFlowProj%gas_num, MaxNumGases)
            slot = firstGas + k - 1
            if (.not. E2Col(slot)%present) cycle
            if (trim(EddyFlowProj%gas(k)%var) /= trim(wantedVar)) cycle
            if (trim(EddyFlowProj%gas(k)%instr) == &
                trim(EddyFlowProj%gas(gasIdx)%instr)) then
                ResolveGasRef = slot
                return
            end if
        end do
    end if

    !> 2. first record of the wanted species, whichever instrument
    do k = 1, min(EddyFlowProj%gas_num, MaxNumGases)
        slot = firstGas + k - 1
        if (.not. E2Col(slot)%present) cycle
        if (trim(EddyFlowProj%gas(k)%var) == trim(wantedVar)) then
            ResolveGasRef = slot
            return
        end if
    end do
end function ResolveGasRef

end subroutine DefineE2Set
