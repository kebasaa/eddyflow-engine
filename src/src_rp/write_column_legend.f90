!***************************************************************************
! write_column_legend.f90
! -----------------------
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
! \brief       Write the small file saying which column suffix is which
!              instrument.
! \author      Jonathan Muller
! \note        A species measured twice is `co2_1` and `co2_2`, and which is
!              which follows the order the project declares its gas records
!              in. That is a serviceable rule and completely opaque from
!              outside: nothing in `h2o_2_flux` says LI-7200. The answer was
!              recoverable only from the project file, or from the per-period
!              `_metadata` output - which does carry `irga_model` under the
!              same tags, but is optional, is one row per averaging period,
!              and so states the mapping tens of thousands of times.
!
!              This states it once. Written unconditionally: it is a handful
!              of lines, and a legend that might not be there is one a reader
!              cannot rely on.
!
!              Every tag comes from the helpers that name the columns -
!              FullOutputGasTags, FluxnetLayoutTags as SelectFluxnetGasSlots
!              filled it, WaterOutSlots. Rebuilding the names here would put
!              the header/row divergence this codebase keeps closing into the
!              one file whose whole job is to explain the headers, where it
!              would be the last place anyone thought to look.
!
!              RP only, though both executables produce output. FCC has no
!              E2Col and no FluxnetLayoutTags, and it cannot run without an ex
!              file that RP wrote from this same project - so the legend beside
!              that file is already correct, and an FCC-only run needs no
!              second copy.
!***************************************************************************
subroutine WriteColumnLegend()
    use m_rp_global_var
    implicit none
    !> local variables
    integer :: open_status
    integer :: dot
    integer :: j
    integer :: gas
    integer :: rec
    character(PathLen) :: Legend_Path
    character(LongOutstringLen) :: csv_row
    character(32) :: gas_tag(GHGNumVar)
    integer :: waterSlots(GHGNumVar)
    character(8) :: waterTags(GHGNumVar)
    integer :: nWater
    integer :: designatedWater
    integer :: designatedCarbon
    character(32) :: species
    character(32), external :: GasOutputLabel
    integer, external :: PrimaryWaterSlot
    integer, external :: PrimaryCarbonSlot

    Legend_Path = Dir%main_out(1:len_trim(Dir%main_out)) &
                // EddyFlowProj%id(1:len_trim(EddyFlowProj%id)) &
                // Legend_FilePadding // Timestamp_FilePadding // CsvExt
    open(uleg, file = Legend_Path, iostat = open_status, encoding = 'utf-8')
    if (open_status /= 0) return

    write(uleg, '(a)') 'columns,file,species,instrument,raw_column,primary'

    designatedWater = PrimaryWaterSlot()
    designatedCarbon = PrimaryCarbonSlot()

    !> The full-output stems, in the order the full output writes them.
    call FullOutputGasTags(gas_tag)
    do gas = firstGas, lastGas
        rec = gas - firstGas + 1
        if (rec > min(EddyFlowProj%gas_num, MaxNumGases)) exit
        if (len_trim(gas_tag(gas)) == 0) cycle

        call clearstr(csv_row)
        !> The stem carries a trailing underscore for concatenation; the
        !> column family is named without it.
        call AddDatum(csv_row, &
            gas_tag(gas)(1:len_trim(gas_tag(gas)) - 1), separator)
        call AddDatum(csv_row, 'full_output', separator)
        species = GasOutputLabel(gas)
        call uppercase(species)
        call AddDatum(csv_row, trim(species), separator)
        call AddDatum(csv_row, trim(EddyFlowProj%gas(rec)%instr), separator)
        call AddIntDatumToDataline(EddyFlowProj%gas(rec)%col, csv_row, &
            EddyFlowProj%err_label)
        call AddDatum(csv_row, &
            trim(YesNo(gas == designatedWater .or. gas == designatedCarbon)), &
            separator)
        write(uleg, '(a)') csv_row(1:len_trim(csv_row) - 1)
    end do

    !> The FLUXNET spellings, which follow a different rule: a required
    !> species keeps its bare standard name for the designated record where
    !> the full output numbers both occurrences. Present only once
    !> SelectFluxnetGasSlots has filled the tags, which it does when the
    !> FLUXNET output is on.
    do j = 1, nFluxnetLayoutSlots
        gas = FluxnetLayoutSlots(j)
        rec = gas - firstGas + 1
        call clearstr(csv_row)
        call AddDatum(csv_row, trim(FluxnetLayoutTags(j)), separator)
        call AddDatum(csv_row, 'fluxnet', separator)
        if (rec >= 1 .and. rec <= min(EddyFlowProj%gas_num, MaxNumGases)) then
            species = GasOutputLabel(gas)
            call uppercase(species)
            call AddDatum(csv_row, trim(species), separator)
            call AddDatum(csv_row, trim(EddyFlowProj%gas(rec)%instr), separator)
            call AddIntDatumToDataline(EddyFlowProj%gas(rec)%col, csv_row, &
                EddyFlowProj%err_label)
            call AddDatum(csv_row, &
                trim(YesNo(gas == designatedWater &
                           .or. gas == designatedCarbon)), separator)
        else
            !> A required variable no record names. The column exists because
            !> FLUXNET requires CO2, H2O and CH4 whatever the site measures,
            !> and carries the error label throughout. Saying so is more use
            !> than leaving the row out and letting a reader wonder.
            call AddDatum(csv_row, trim(FluxnetLayoutTags(j)), separator)
            call AddDatum(csv_row, 'not measured', separator)
            call AddDatum(csv_row, '', separator)
            call AddDatum(csv_row, 'no', separator)
        end if
        write(uleg, '(a)') csv_row(1:len_trim(csv_row) - 1)
    end do

    !> The per-hygrometer flux families. Six column names share one suffix,
    !> so one row states it for all six.
    call WaterOutSlots(waterSlots, waterTags, nWater)
    do j = 1, nWater
        gas = waterSlots(j)
        rec = gas - firstGas + 1
        if (rec < 1 .or. rec > min(EddyFlowProj%gas_num, MaxNumGases)) cycle
        call clearstr(csv_row)
        call AddDatum(csv_row, 'H' // trim(waterTags(j)) // ' LE' &
            // trim(waterTags(j)) // ' ET' // trim(waterTags(j)) // ' TAU' &
            // trim(waterTags(j)) // ' MO_LENGTH' // trim(waterTags(j)) &
            // ' ZL' // trim(waterTags(j)), separator)
        call AddDatum(csv_row, 'both', separator)
        call AddDatum(csv_row, 'H2O', separator)
        call AddDatum(csv_row, trim(EddyFlowProj%gas(rec)%instr), separator)
        call AddIntDatumToDataline(EddyFlowProj%gas(rec)%col, csv_row, &
            EddyFlowProj%err_label)
        call AddDatum(csv_row, trim(YesNo(gas == designatedWater)), separator)
        write(uleg, '(a)') csv_row(1:len_trim(csv_row) - 1)
    end do

    close(uleg)

contains

    !> Spelled out rather than 1/0: the file is meant to be read by a person
    !> who is already puzzled about what a suffix means.
    function YesNo(flag) result(word)
        logical, intent(in) :: flag
        character(3) :: word
        if (flag) then
            word = 'yes'
        else
            word = 'no'
        end if
    end function YesNo

end subroutine WriteColumnLegend
