!***************************************************************************
! define_all_var_set.f90
! ----------------------
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
! \brief       Creates a dataset with with all variables in phisical and \n
!              standardized units
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine DefineAllVarSet(LocCol, fRaw, nrow, ncol, N)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow
    integer, intent(in) :: ncol
    integer, intent(in) :: N
    type(ColType), intent(inout) :: LocCol(MaxNumCol)
    real(kind = sgl), intent(inout) :: fRaw(nrow, ncol)
    !> local variables
    integer :: j
    integer :: gasSlot
    real(kind = sgl) :: DumVec(N)
    integer, external :: RecordGasSlot
    integer, external :: HistoricGasSlot
    logical, external :: IsHistoricGasVar

    !> Converts units, according to information in the metadata file
    !> Physical, non-standard units are detected and converted into standard
    !> Volt (or any other "non-physical") are converted as well, using conversion parameters
    do j = 1, NumAllVar
        !> A gas named by a record is a trace gas whatever its species, but the
        !> select below can only match literal names. Without this, a column
        !> holding COS - or any gas in a slot past the fourth - falls through
        !> to `case default` and receives no unit conversion at all: a ppb
        !> reading is then used as though it were ppm, a silent factor of 1000.
        !> Legacy projects escaped this only because DefineUsedVariables
        !> renamed the fourth gas's column to 'n2o' so it borrowed N2O's arm.
        gasSlot = RecordGasSlot(LocCol(j)%orig_col)
        if (gasSlot > 0 .and. .not. IsHistoricGasVar(LocCol(j)%var)) then
            call ConvertTraceGasUnits(LocCol, fRaw, nrow, ncol, N, j, gasSlot)
            cycle
        end if

        select case (LocCol(j)%var)
            !> Wind components, taken to [m s-1]
            case('u', 'v', 'w', 'sos')
                if (LocCol(j)%conversion_type == 'none') then
                    select case(LocCol(j)%unit_in(1:len_trim(LocCol(j)%unit_in)))
                        case ('m_sec')
                            if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) /= 'sos') cycle
                        case ('mm_sec')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) * 1e-3
                            end where
                            if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) /= 'sos') cycle
                        case ('cm_sec')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) * 1e-2
                            end where
                            if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) /= 'sos') cycle
                        case default
                            cycle
                    end select
                else
                    DumVec(1:N) = fRaw(1:N, j)
                    call LinearConversion(LocCol, DumVec(1:N), N, j)
                    fRaw(1:N, j) = DumVec(1:N)
                    select case(LocCol(j)%unit_out(1:len_trim(LocCol(j)%unit_out)))
                        case ('m_sec')
                            if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) /= 'sos') cycle
                        case ('mm_sec')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) * 1e-3
                            end where
                            if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) /= 'sos') cycle
                        case ('cm_sec')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) * 1e-2
                            end where
                            if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) /= 'sos') cycle
                        case default
                            if (LocCol(j)%var(1:len_trim(LocCol(j)%var)) /= 'sos') cycle
                        end select
                end if
                !> if speed-of-sound, it is taken from [m s-1] to sonic temperature [K]
                if(LocCol(j)%var(1:len_trim(LocCol(j)%var)) == 'sos') then
                    where(fRaw(1:N, j) /= error)
                        fRaw(1:N, j) = (fRaw(1:N, j))**2 / 403.
                    end where
                    cycle
                end if

            !> Temperatures (K)
            case('ts', 'cell_t', 'int_t_1', 'int_t_2', 'air_t')
                if (LocCol(j)%conversion_type == 'none') then
                    select case(LocCol(j)%unit_in(1:len_trim(LocCol(j)%unit_in)))
                        case ('kelvin')
                            cycle
                        case ('ckelvin')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) * 1e-2
                            end where
                            cycle
                        case ('celsius')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) + 273.15
                            end where
                        case ('ccelsius')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) * 1e-2 + 273.15
                            end where
                            cycle
                        case default
                            cycle
                    end select
                else
                    DumVec(1:N) = fRaw(1:N, j)
                    call LinearConversion(LocCol, DumVec(1:N), N, j)
                    fRaw(1:N, j) = DumVec(1:N)
                    select case(LocCol(j)%unit_out(1:len_trim(LocCol(j)%unit_out)))
                        case ('kelvin')
                            cycle
                        case ('ckelvin')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) * 1e-2
                            end where
                            cycle
                        case ('celsius')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) + 273.15
                            end where
                        case ('ccelsius')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) * 1e-2 + 273.15
                            end where
                            cycle
                        case default
                            cycle
                    end select
                end if

            !> Pressures (Pa)
            case('int_p', 'air_p')
                if (LocCol(j)%conversion_type == 'none') then
                    select case(LocCol(j)%unit_in(1:len_trim(LocCol(j)%unit_in)))
                        case ('pa')
                            cycle
                        case ('hpa')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) * 1e2
                            end where
                            cycle
                        case ('kpa')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) * 1e3
                            end where
                            cycle
                    end select
                else
                    DumVec(1:N) = fRaw(1:N, j)
                    call LinearConversion(LocCol, DumVec(1:N), N, j)
                    fRaw(1:N, j) = DumVec(1:N)
                    select case(LocCol(j)%unit_out(1:len_trim(LocCol(j)%unit_out)))
                        case ('pa')
                            cycle
                        case ('hpa')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) * 1e2
                            end where
                            cycle
                        case ('kpa')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) * 1e3
                            end where
                            cycle
                        case default
                            cycle
                    end select
                end if

            !> Concentrations of trace gases
            !> Concentrations of trace gases. The conversion is identical for
            !> every species once the molecular weight is known, so it lives in
            !> one place and is keyed on the gas slot; a record-named gas is
            !> routed to the same code above.
            case('co2', 'ch4', 'n2o')
                call ConvertTraceGasUnits(LocCol, fRaw, nrow, ncol, N, j, &
                                          HistoricGasSlot(LocCol(j)%var))
                cycle


            !> Concentrations of water vapour
            case('h2o')
                if (LocCol(j)%conversion_type == 'none') then
                    select case(LocCol(j)%unit_in(1:len_trim(LocCol(j)%unit_in)))
                        case ('mmol_m3', 'ppt')
                            cycle
                        case ('umol_m3', 'ppm')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) * 1e-3
                            end where
                        case ('ppb')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) * 1e-6
                            end where
                        case ('g_m3')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) / MW(h2o)
                            end where
                        case ('mg_m3')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) / MW(h2o) * 1e-3
                            end where
                        case ('ug_m3')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) / MW(h2o) * 1e-6
                            end where
                        case default
                            cycle
                    end select
                else
                    DumVec(1:N) = fRaw(1:N, j)
                    call LinearConversion(LocCol, DumVec(1:N), N, j)
                    fRaw(1:N, j) = DumVec(1:N)
                    select case(LocCol(j)%unit_out(1:len_trim(LocCol(j)%unit_out)))
                        case ('umol_m3', 'ppm')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) * 1e-3
                            end where
                        case ('ppb')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) * 1e-6
                            end where
                        case ('g_m3')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) / MW(h2o)
                            end where
                        case ('mg_m3')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) / MW(h2o) * 1e-3
                            end where
                        case ('ug_m3')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) / MW(h2o) * 1e-6
                            end where
                        case default
                            cycle
                    end select
                end if

            !> Flow rates (m+3s-1)
            case('flowrate')
                if (LocCol(j)%conversion_type == 'none') then
                    select case(LocCol(j)%unit_in(1:len_trim(LocCol(j)%unit_in)))
                        case ('m3_s')
                            cycle
                        case ('cm3_s')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) * 1e-6
                            end where
                            cycle
                        case ('lit_m')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) *  1.66666667e-5
                            end where
                            cycle
                        case ('ft3_s')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) * 0.028316846592e0
                            end where
                            cycle
                        case ('in3_s')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) * 1.6387064e-5
                            end where
                            cycle
                        case default
                            cycle
                    end select
                else
                    DumVec(1:N) = fRaw(1:N, j)
                    call LinearConversion(LocCol, DumVec(1:N), N, j)
                    fRaw(1:N, j) = DumVec(1:N)
                    select case(LocCol(j)%unit_out(1:len_trim(LocCol(j)%unit_out)))
                        case ('m3_s')
                            cycle
                        case ('cm3_s')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) * 1e-6
                            end where
                            cycle
                        case ('lit_m')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) *  1.66666667e-5
                            end where
                            cycle
                        case ('ft3_s')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) * 0.028316846592
                            end where
                            cycle
                        case ('in3_s')
                            where(fRaw(1:N, j) /= error)
                                fRaw(1:N, j) = fRaw(1:N, j) * 1.6387064e-5
                            end where
                            cycle
                    end select
                end if

            !> All other variables (custom variables)
            case default
                DumVec(1:N) = fRaw(1:N, j)
                call LinearConversion(LocCol, DumVec(1:N), N, j)
                fRaw(1:N, j) = DumVec(1:N)
        end select
    end do
end subroutine DefineAllVarSet

!***************************************************************************
!
! \brief       Performs a linear rescaling as either gain/offset or min/max.
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine LinearConversion(LocCol, Vec, nrow, j)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow
    integer, intent(in) :: j
    type(ColType), intent(in) :: LocCol(MaxNumCol)
    real(kind = sgl), intent(inout) :: Vec(nrow)

    select case (LocCol(j)%conversion_type(1:len_trim(LocCol(j)%conversion_type)))
        case ('zero_fullscale')
            where(Vec(:) /= error)
                Vec(:) = sngl(((LocCol(j)%b - LocCol(j)%a) / &
                    (LocCol(j)%max - LocCol(j)%min)) * (dble(Vec(:)) - LocCol(j)%min) + LocCol(j)%a)
            end where
        case ('gain_offset')
            where(Vec(:) /= error)
                Vec(:) = sngl(LocCol(j)%a * dble(Vec(:)) + LocCol(j)%b)
            end where
        case default
            continue
    end select
end subroutine LinearConversion

!***************************************************************************
!
! \brief       Gas slot a metadata column is assigned to by the project's
!              gas records, or 0 if no record names it.
! \author      Jonathan Muller
! \note        Runs while LocCol is still indexed by metadata column number,
!              which is exactly what a record stores, so no translation is
!              needed here.
!***************************************************************************
integer function RecordGasSlot(orig_col)
    use m_common_global_var
    implicit none
    integer, intent(in) :: orig_col
    integer :: i

    RecordGasSlot = 0
    if (orig_col <= 0) return
    do i = 1, min(EddyFlowProj%gas_num, MaxNumGases)
        if (EddyFlowProj%gas(i)%col /= orig_col) cycle
        if (firstGas + i - 1 > lastGas) return
        RecordGasSlot = firstGas + i - 1
        return
    end do
end function RecordGasSlot

!***************************************************************************
!
! \brief       Whether a column name is one the unit conversion already
!              dispatches on by name.
! \author      Jonathan Muller
! \note        h2o is included so a record naming water vapour keeps using
!              the water arm, whose conversions differ from the trace gases'.
!***************************************************************************
logical function IsHistoricGasVar(var)
    implicit none
    character(*), intent(in) :: var

    select case (trim(adjustl(var)))
        case ('co2', 'h2o', 'ch4', 'n2o'); IsHistoricGasVar = .true.
        case default;                      IsHistoricGasVar = .false.
    end select
end function IsHistoricGasVar

!***************************************************************************
!
! \brief       Gas slot behind one of the historically named gas columns.
! \author      Jonathan Muller
! \note        'n2o' maps to gas4: the fourth slot has always been N2O's by
!              default, and legacy projects renamed whatever gas occupied it
!              to 'n2o' so that it would reach this mapping at all.
!***************************************************************************
integer function HistoricGasSlot(var)
    use m_common_global_var
    implicit none
    character(*), intent(in) :: var

    select case (trim(adjustl(var)))
        case ('co2'); HistoricGasSlot = co2
        case ('h2o'); HistoricGasSlot = h2o
        case ('ch4'); HistoricGasSlot = ch4
        case ('n2o'); HistoricGasSlot = gas4
        case default; HistoricGasSlot = gas4
    end select
end function HistoricGasSlot

!***************************************************************************
!
! \brief       Converts one trace-gas column to the standard mixing-ratio
!              units, using the molecular weight of the slot it occupies.
! \author      Gerardo Fratini, generalised by Jonathan Muller
! \note        Extracted verbatim from the co2/ch4/n2o arm of DefineAllVarSet
!              so that the named gases and the record-named gases cannot
!              drift apart. The only change is that the molecular weight is
!              taken from the gas slot instead of a three-way switch on the
!              column name, which is what allows any species to be handled.
!***************************************************************************
subroutine ConvertTraceGasUnits(LocCol, fRaw, nrow, ncol, N, j, gas_slot)
    use m_common_global_var
    implicit none
    integer, intent(in) :: nrow, ncol, N, j, gas_slot
    type(ColType), intent(in) :: LocCol(MaxNumCol)
    real(kind = sgl), intent(inout) :: fRaw(nrow, ncol)
    real(kind = sgl) :: DumVec(N)

    if (LocCol(j)%conversion_type == 'none') then
        select case(LocCol(j)%unit_in(1:len_trim(LocCol(j)%unit_in)))
            case ('mmol_m3', 'ppm')
                return
            !> The three scalings below are deliberately left unguarded, as
            !> they were in the code this was extracted from. Adding a
            !> `where (/= error)` here would change every existing gas
            !> result, so it is a separate decision, not a tidy-up.
            case ('ppt')
                fRaw(1:N, j) = fRaw(1:N, j) * 1e3
            case ('umol_m3', 'ppb')
                fRaw(1:N, j) = fRaw(1:N, j) * 1e-3
            case ('pmol_mol')
                fRaw(1:N, j) = fRaw(1:N, j) * 1e-6
            case ('g_m3')
                where(fRaw(1:N, j) /= error)
                    fRaw(1:N, j) = fRaw(1:N, j) / MW(gas_slot)
                end where
            case ('mg_m3')
                where(fRaw(1:N, j) /= error)
                    fRaw(1:N, j) = fRaw(1:N, j) / MW(gas_slot) * 1e-3
                end where
            case ('ug_m3')
                where(fRaw(1:N, j) /= error)
                    fRaw(1:N, j) = fRaw(1:N, j) / MW(gas_slot) * 1e-6
                end where
            case default
                return
        end select
    else
        DumVec(1:N) = fRaw(1:N, j)
        call LinearConversion(LocCol, DumVec(1:N), N, j)
        fRaw(1:N, j) = DumVec(1:N)
        select case(LocCol(j)%unit_out(1:len_trim(LocCol(j)%unit_out)))
            case ('ppt')
                where(fRaw(1:N, j) /= error)
                    fRaw(1:N, j) = fRaw(1:N, j) * 1e3
                end where
            case ('umol_m3', 'ppb')
                where(fRaw(1:N, j) /= error)
                    fRaw(1:N, j) = fRaw(1:N, j) * 1e-3
                end where
            case ('pmol_mol')
                where(fRaw(1:N, j) /= error)
                    fRaw(1:N, j) = fRaw(1:N, j) * 1e-6
                end where
            case ('g_m3')
                where(fRaw(1:N, j) /= error)
                    fRaw(1:N, j) = fRaw(1:N, j) / MW(gas_slot)
                end where
            case ('mg_m3')
                where(fRaw(1:N, j) /= error)
                    fRaw(1:N, j) = fRaw(1:N, j) / MW(gas_slot) * 1e-3
                end where
            case ('ug_m3')
                where(fRaw(1:N, j) /= error)
                    fRaw(1:N, j) = fRaw(1:N, j) / MW(gas_slot) * 1e-6
                end where
        end select
    end if
end subroutine ConvertTraceGasUnits

!***************************************************************************
!
! \brief       Label to use for a gas slot in output headers.
! \author      Jonathan Muller
! \note        Headers are written before the first data file is read, so
!              E2Col is still empty in a record project - ApplyGasRecords
!              fills it per file. Legacy projects had the label by then
!              because the retired col_gas4 tag pointed straight at the
!              metadata column. The record names the same column, so this
!              resolves it the same way and yields the same label.
!
!              Takes the slot rather than assuming the fourth: the full
!              output carries a column set per configured gas, so every slot
!              needs a name, and deriving one only for slot four is what left
!              gases 5+ out of that file entirely.
!***************************************************************************
function GasOutputLabel(gas_slot) result(label)
    use m_common_global_var
    implicit none
    integer, intent(in) :: gas_slot
    character(32) :: label
    integer :: rec4
    integer :: i

    rec4 = gas_slot - firstGas + 1

    call clearstr(label)
    !> The record is the authority, exactly as it is for the input unit. It
    !> was tempting to prefer E2Col when populated, on the grounds that the
    !> per-file path has "already resolved" the slot - but E2Col is also
    !> filled by name matching before the records are applied, so with five
    !> gases configured it could hand back the species of a different slot.
    !> That produced a FLUXNET row with two sets of N2O_* columns and no COS
    !> columns at all.
    if (EddyFlowProj%gas_num >= rec4 .and. EddyFlowProj%gas(rec4)%col > 0) then
        !> Prefer the metadata column's own label, which is what the legacy
        !> path used, so an upgraded project keeps its column names.
        do i = 1, MaxNumCol
            if (Col(i)%orig_col /= EddyFlowProj%gas(rec4)%col) cycle
            if (len_trim(Col(i)%label) == 0 .or. &
                trim(Col(i)%label) == 'none') exit
            label = trim(Col(i)%label)
            return
        end do
        !> No usable label: fall back to the species the record names.
        if (len_trim(EddyFlowProj%gas(rec4)%var) > 0 .and. &
            trim(EddyFlowProj%gas(rec4)%var) /= 'none') then
            label = trim(EddyFlowProj%gas(rec4)%var)
            return
        end if
    end if

    !> No record for this slot: fall back to whatever the per-file path
    !> resolved, then to the slot name. 'gas4' is kept for the fourth slot so
    !> a project that resolves nothing still produces the historical name.
    if (len_trim(E2Col(gas_slot)%label) > 0 .and. &
        trim(E2Col(gas_slot)%label) /= 'none') then
        label = trim(E2Col(gas_slot)%label)
        return
    end if
    if (gas_slot == gas4) then
        label = 'gas4'
    else
        write(label, '(a,i0)') 'gas', gas_slot - firstGas + 1
    end if
end function GasOutputLabel

!***************************************************************************
!
! \brief       Input unit of a gas slot, for the output-units decision.
! \author      Jonathan Muller
! \note        Same problem, and same resolution, as GasOutputLabel: the unit
!              decides whether the full output is written in nmol or umol,
!              and it is consulted where E2Col has not been filled from the
!              records yet. Getting it wrong does not corrupt the numbers -
!              the label and the scaling move together - but it silently
!              changes the units an upgraded project reports in.
!***************************************************************************
function GasUnitIn(gas_slot) result(unit_in)
    use m_common_global_var
    implicit none
    integer, intent(in) :: gas_slot
    character(32) :: unit_in
    integer :: rec4
    integer :: i

    rec4 = gas_slot - firstGas + 1
    call clearstr(unit_in)
    !> The metadata column is the authority here: unit_in describes what the
    !> data file contains, which no amount of processing changes. E2Col's copy
    !> is only a fallback, for a project that names no records.
    if (EddyFlowProj%gas_num >= rec4 .and. EddyFlowProj%gas(rec4)%col > 0) then
        do i = 1, MaxNumCol
            if (Col(i)%orig_col /= EddyFlowProj%gas(rec4)%col) cycle
            if (len_trim(Col(i)%unit_in) > 0) unit_in = trim(Col(i)%unit_in)
            return
        end do
    end if

    if (len_trim(E2Col(gas_slot)%unit_in) > 0 .and. &
        trim(E2Col(gas_slot)%unit_in) /= 'none') unit_in = trim(E2Col(gas_slot)%unit_in)
end function GasUnitIn
