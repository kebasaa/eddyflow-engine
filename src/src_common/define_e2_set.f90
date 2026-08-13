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
    integer :: wsl
    !> The hygrometer a gas is corrected with, when reporting the ones that
    !> borrow across analysers.
    integer :: msl
    !> Gases already reported as borrowing. Saved, because this routine runs
    !> once per averaging period and the report belongs to the run.
    logical, save :: crossWaterWarned(GHGNumVar) = .false.
    !> A gas record's index and the cell record it names, for resolving an
    !> explicit gas_<i>_cell.
    integer :: gasrec
    integer :: cellrec
    character(32) :: cellInstr(MaxNumInstruments)
    !> Whether any cell record names an analyser. Decides whether a gas that
    !> matches none of them falls back to the first block or has no cell.
    logical :: cellInstrNamed
    integer, external :: PrimaryWaterSlot
    logical, external :: GasSlotIsWater
    character(32), external :: GasOutputLabel
    !> Declared in the host so the internal procedures inherit it. The lookup
    !> moved to file scope for DefineVars, which asks the same question on the
    !> earlier pass.
    integer, external :: LocColByOrigCol


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

    !> The remaining standard variables: cell temperatures and pressures, and
    !> ambient temperature and pressure.
    !>
    !> The gases are not here. Four branches used to match a column's species
    !> name onto a fixed slot - co2 to five, h2o to six and so on - and
    !> ApplyGasRecords then cleared firstGas..lastGas and refilled it from the
    !> records a few lines below, so every one of those assignments was
    !> overwritten before anything read it. What the branches could express
    !> was also strictly less: one column per species, and nothing at all past
    !> the fourth.
    do j = 1, ncol
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
    !>
    !> A record's own gas_<i>_cell wins. It names the cell record whose
    !> analyser holds this gas, and it exists for the cases the instrument
    !> name cannot express: two analysers reported under one model string, or
    !> a gas plumbed through a cell that is not its own instrument's. It was
    !> parsed and then discarded - ApplyGasRecords set cell_ref to zero and
    !> this loop overwrote it from the name match - because the cell slots
    !> were one global set and there was nothing for an index to select. They
    !> have been per-instrument since, so the field can mean what it says.
    !> Whether the cell records name analysers at all.
    !>
    !> They either do or they do not, and the answer decides what an unmatched
    !> gas means. Records that name nobody describe the site's single cell,
    !> which is every pre-record project, and every gas uses it. Records that
    !> name analysers are making a claim about which cell is whose, and a gas
    !> none of them names has no cell record - not the first one's.
    cellInstrNamed = .false.
    do i = 1, nCellInstr
        if (len_trim(cellInstr(i)) > 0) cellInstrNamed = .true.
    end do

    do j = firstGas, lastGas
        if (.not. E2Col(j)%present) cycle
        !> Unresolved until this gas's own cell is found. It defaulted to
        !> firstCell, which reads as "every gas has a cell" and on a site whose
        !> analysers do not all have one handed the first analyser's cell
        !> temperature and pressure to gases measured somewhere else. Va, the
        !> point-by-point dilution and the cell-temperature covariance all read
        !> this field, so one default silently fed three corrections.
        E2Col(j)%cell_ref = 0

        gasrec = j - firstGas + 1
        cellrec = 0
        if (gasrec >= 1 .and. gasrec <= min(EddyFlowProj%gas_num, MaxNumGases)) &
            cellrec = EddyFlowProj%gas(gasrec)%cell
        if (cellrec >= 1 .and. &
            cellrec <= min(EddyFlowProj%cell_num, MaxNumCellCols)) then
            i = CellInstrIndex(trim(EddyFlowProj%cell(cellrec)%instr))
            if (i >= 1) then
                E2Col(j)%cell_ref = firstCell + (i - 1) * NumCellPerInstr
                cycle
            end if
        end if

        !> Undescribed cells belong to everyone; named cells belong to the
        !> analyser they name.
        if (nCellInstr <= 0 .or. .not. cellInstrNamed) then
            if (nCellInstr > 0) E2Col(j)%cell_ref = firstCell
            cycle
        end if
        do i = 1, nCellInstr
            if (trim(cellInstr(i)) == trim(E2Col(j)%instr%model) .or. &
                trim(cellInstr(i)) == trim(E2Col(j)%Instr%ep_label)) then
                E2Col(j)%cell_ref = firstCell + (i - 1) * NumCellPerInstr
                exit
            end if
        end do
    end do

    !> Default every unresolved gas to the site's primary water. The legacy
    !> selection has exactly one H2O, so this is what the single-analyser case
    !> has always used implicitly; making it explicit means both paths describe
    !> their moisture the same way instead of the legacy one leaving it blank.
    !>
    !> Resolved by species rather than taken from the fixed slot. That slot is
    !> record two, which is water only by convention - a project declaring its
    !> water elsewhere had every unresolved gas pointed at whatever species
    !> record two happened to hold, and the WPL correction then ran against a
    !> trace gas's density as though it were humidity.
    !> A biomet reference is not unresolved and must survive this. Left to the
    !> bounds test alone it would not: biometMoistRef is outside firstGas..
    !> lastGas by construction, so every gas the user pointed at the biomet
    !> would be quietly re-pointed at the primary hygrometer here, and the
    !> selection in the interface would do nothing. The bug would show up as a
    !> WPL correction against the wrong water, which is the one failure mode
    !> this file has spent the most comments on.
    wsl = PrimaryWaterSlot()
    if (wsl >= firstGas) then
        if (E2Col(wsl)%present) then
            do j = firstGas, lastGas
                if (.not. E2Col(j)%present) cycle
                if (E2Col(j)%moist_ref == biometMoistRef) cycle
                if (E2Col(j)%moist_ref < firstGas .or. E2Col(j)%moist_ref > lastGas) &
                    E2Col(j)%moist_ref = wsl
            end do
        end if
    end if

    !> Say which gases are corrected with another analyser's water.
    !>
    !> A legitimate configuration - it is what a site with one hygrometer and
    !> two analysers has - but a compromise, and the compromise used to be
    !> invisible. The dilution and the water-flux covariance both declined the
    !> pairing outright while the mean WPL terms took it, so the gas came out
    !> corrected by a water it was held not to share, with nothing in the
    !> output to show for it. Both honour it now, so the one thing left to do
    !> is say so.
    !>
    !> Once per gas for the whole run. This routine is called for every
    !> averaging period, so an unguarded message is one per period per gas -
    !> ninety-seven of them on a single day of CH-LAE, which buries the log
    !> the warning exists to be read in.
    do j = firstGas, lastGas
        if (.not. E2Col(j)%present) cycle
        if (GasSlotIsWater(j)) cycle
        msl = E2Col(j)%moist_ref
        if (msl < firstGas .or. msl > lastGas) cycle
        if (j == msl .or. .not. E2Col(msl)%present) cycle
        if (E2Col(j)%instr%model == E2Col(msl)%instr%model) cycle
        if (crossWaterWarned(j)) cycle
        crossWaterWarned(j) = .true.
        write(*, '(a)') '  Warning(106)> ' // trim(GasOutputLabel(j)) // ' on ' &
            // trim(E2Col(j)%instr%model) // ' is corrected with the water on ' &
            // trim(E2Col(msl)%instr%model) // '.'
        call ExceptionHandler(106)
    end do

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
!> Each instrument owns a block of four slots - cell_t, int_t_1, int_t_2,
!> int_p - keyed on the record's instrument rather than on arrival order, so a
!> site with two analysers keeps both. This carried a LIMITATION note saying
!> the slots were one global set and a second record overwrote the first, and
!> that GasRecordType%cell therefore stayed unresolved; the slots were widened
!> and the note outlived the limitation it described. DefineE2Set honours
!> gas_<i>_cell now.
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

    !> Records are the only way a project names its gases, so a file that
    !> states no gas count at all is a pre-5.0.0 project that has not been
    !> through the interface. Processing it would silently produce no gas
    !> fluxes, so it is refused; the GUI migrates such a file on open.
    !>
    !> A file that states `gas_num=0` is a different thing entirely: a site
    !> with an anemometer and no analyser. It measures wind and sonic
    !> temperature, so it has a heat flux, and it is processed normally. This
    !> test used to be `gas_num <= 0`, which refused it.
    if (.not. EddyFlowProj%gas_num_stated) then
        write(*, '(a)') '  Fatal error(99)> Project file states no gas count &
            &(gas_num). Open and save it in the EddyFlow interface to bring &
            &it up to the current format.'
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
        !> Cleared, not resolved: the cell blocks are not filled until
        !> ApplyCellDiagRecords has run, so DefineE2Set resolves this once both
        !> are populated - honouring gas_<i>_cell first and falling back on the
        !> instrument name.
        E2Col(slot)%cell_ref = 0
    end do
end subroutine ApplyGasRecords

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
    character(32) :: candidate, wanted

    !> Species matching is case-insensitive, as it is in GasSlotIsWater.
    !> The interface normalises a slug to lower case when it creates it, but
    !> writes it to the file verbatim and re-normalises nothing on the way
    !> back in - so a hand-edited project naming H2O rather than h2o would
    !> have resolved no moisture reference at all, silently.
    wanted = wantedVar
    call uppercase(wanted)

    ResolveGasRef = 0

    !> 0. the biomet, named explicitly. Honoured whether or not the project
    !> also has a hygrometer - the point of naming it is to say which of the
    !> two you want. Only if the project has no biomet RH column at all does
    !> this fall through to the automatic rules, since the alternative is a
    !> reference to a measurement that does not exist.
    if (ref == biometMoistRef) then
        if (BiometRhConfigured) then
            ResolveGasRef = biometMoistRef
            return
        end if
    elseif (ref > 0 .and. ref <= EddyFlowProj%gas_num) then
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
            candidate = EddyFlowProj%gas(k)%var
            call uppercase(candidate)
            if (trim(adjustl(candidate)) /= trim(wanted)) cycle
            if (trim(EddyFlowProj%gas(k)%instr) == &
                trim(EddyFlowProj%gas(gasIdx)%instr)) then
                ResolveGasRef = slot
                return
            end if
        end do
    end if

    !> 2. the biomet relative humidity, for water.
    !>
    !> This step used to be "the first record of the wanted species, whichever
    !> instrument" - so a gas whose own analyser carried no hygrometer silently
    !> borrowed another analyser's water, taken through a different cell at a
    !> different time lag. That is the compromise Warning(106) exists to
    !> announce, and announcing it was the best that could be done while it was
    !> the only fallback there was. A site RH sensor is the better answer: it
    !> measures the air rather than the inside of an unrelated instrument.
    !>
    !> Borrowing is still reachable, but only by asking for it - a user who
    !> names another analyser's H2O in the interface gets it, and gets the
    !> warning. What has gone is arriving there without being asked.
    !>
    !> Water only. The function takes a species because it was written to be
    !> general, but the biomet has a humidity and nothing else, so this arm
    !> cannot answer for any other gas.
    if (trim(wanted) == 'H2O' .and. BiometRhConfigured) then
        ResolveGasRef = biometMoistRef
        return
    end if

    !> Nothing resolved - and the caller does not leave it there. The
    !> default-to-primary block in DefineE2Set points such a gas at the site's
    !> primary hygrometer, so on a project with no biomet RH a gas whose own
    !> analyser has none is still corrected with another analyser's water, and
    !> still says so through Warning(106).
    !>
    !> So dropping the species search above narrows what happens automatically
    !> rather than abolishing it: with a biomet RH the gas takes that instead
    !> of an unrelated cell, and without one the outcome is what it always was.
end function ResolveGasRef

end subroutine DefineE2Set

!***************************************************************************
!
! \brief       Cell-pressure slot of the analyser that measured `gas`.
! \author      Jonathan Muller
! \note        Offset 3 of that instrument's cell block, matching the layout
!              ApplyCellDiagRecords writes and AirAndCellParameters reads.
!              External rather than contained, because the flux code and the
!              FLUXNET writer must agree on it: the writer puts this
!              covariance in the file and the flux code consumes it, so a
!              second copy of the arithmetic is a silent mismatch waiting to
!              happen. Falls back to instrument 1's `pi`, which is where a
!              single-analyser project's pressure has always been.
!***************************************************************************
integer function cellPressureSlot(gas) result(slot)
    use m_common_global_var
    implicit none
    integer, intent(in) :: gas

    slot = pi
    if (gas < firstGas .or. gas > lastGas) return
    if (E2Col(gas)%cell_ref < firstCell .or. E2Col(gas)%cell_ref > lastCell) return
    slot = E2Col(gas)%cell_ref + 3
end function cellPressureSlot

!***************************************************************************
!> Translate a .metadata column number into its position in LocCol/Raw.
!>
!> Columns the project does not use are dropped before this point, so the two
!> numbering schemes diverge. The project file stores metadata numbers, so a
!> record must be matched on %orig_col rather than indexed directly - doing
!> the latter silently selects a different variable.
!>
!> At file scope rather than contained in DefineE2Set, because DefineVars asks
!> the same question of the same records on the earlier pass and had no way to
!> reach it. Two copies of "which column is metadata column N" is the class of
!> drift this effort keeps unwinding.
!***************************************************************************
integer function LocColByOrigCol(LocCol, ncol, wanted)
    use m_common_global_var
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
