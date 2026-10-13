!***************************************************************************
! write_out_metadata_fcc.f90
! --------------------------
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
! \brief       Write results on output files
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine WriteOutMetadataFcc(lEx)
    use m_fx_global_var
    implicit none
    !> in/out variables
    Type(ExType), intent(in) :: lEx
    character(LongOutstringLen) :: csv_row

    !> local variables
    integer :: gas
    character(DatumLen) :: field_val
    include '../src_common/interfaces_1.inc'


    call clearstr(csv_row)
    !> Preliminary timestmap information
    write(field_val, *) lEx%fname(1:len_trim(lEx%fname))
    call AddDatum(csv_row, field_val, separator)
    write(field_val, *) lEx%end_date(1:10)
    call AddDatum(csv_row, field_val, separator)
    write(field_val, *) lEx%end_time(1:5)
    call AddDatum(csv_row, field_val, separator)
    !> DOY as RP writes it in its own metadata file: four decimals, cut
    !> rather than rounded. The two files are one file - which program
    !> writes it depends only on whether FCC follows - so they share a
    !> layout; without this every column after time was off by one.
    write(field_val, *) float_doy
    call ShrinkString(field_val)
    call AddDatum(csv_row, field_val(1:index(field_val, '.') + 4), separator)

    !> Site location and characteristics
    write(field_val, *) lEx%lat
    call AddDatum(csv_row, field_val, separator)
    write(field_val, *) lEx%lon
    call AddDatum(csv_row, field_val, separator)
    write(field_val, *) lEx%alt
    call AddDatum(csv_row, field_val, separator)
    write(field_val, *) lEx%canopy_height
    call AddDatum(csv_row, field_val, separator)
    write(field_val, *) lEx%disp_height
    call AddDatum(csv_row, field_val, separator)
    write(field_val, *) lEx%rough_length
    call AddDatum(csv_row, field_val, separator)

    !> Acquisition setup
    write(field_val, *) lEx%file_length
    call AddDatum(csv_row, field_val, separator)
    write(field_val, *) lEx%ac_freq
    call AddDatum(csv_row, field_val, separator)
    !> Master sonic height and north offset
    write(field_val, *) lEx%instr(sonic)%firm(1:len_trim(lEx%instr(sonic)%firm))
    call AddDatum(csv_row, field_val, separator)
    write(field_val, *) lEx%instr(sonic)%model(1:len_trim(lEx%instr(sonic)%model))
    call AddDatum(csv_row, field_val, separator)
    write(field_val, *) lEx%instr(sonic)%height
    call AddDatum(csv_row, field_val, separator)
    write(field_val, *) lEx%instr(sonic)%wformat
    call AddDatum(csv_row, field_val, separator)
    write(field_val, *) lEx%instr(sonic)%wref
    call AddDatum(csv_row, field_val, separator)
    write(field_val, *) lEx%instr(sonic)%north_offset
    call AddDatum(csv_row, field_val, separator)
    !> Lengths in the units RP writes and the .metadata file states: the
    !> ex reader converted them to SI, so they go back here. They used to
    !> be written as read, so the same project's metadata said 0.25 when
    !> FCC wrote it and 25.0 when RP did.
    write(field_val, *) OutUnits(lEx%instr(sonic)%hpath_length, 1d2)
    call AddDatum(csv_row, field_val, separator)
    write(field_val, *) OutUnits(lEx%instr(sonic)%vpath_length, 1d2)
    call AddDatum(csv_row, field_val, separator)
    write(field_val, *) lEx%instr(sonic)%tau
    call AddDatum(csv_row, field_val, separator)
    !> irgas
    !> One block per configured gas, matching the header loop in InitOutFiles.
    !>
    !> This indexed lEx%instr, which has one entry per instrument *role* -
    !> co2, h2o, ch4, the fourth gas, the sonic - so it could never describe
    !> more than four analysers, and past the fourth the role index addresses
    !> an unrelated instrument. lEx%gas_instr is per gas *slot* and the reader
    !> mirrors the four historical analysers into it after their unit
    !> conversions, so the values here are unchanged for a four-gas project.
    do gas = firstGas, lastGas
        if (fcc_var_present(gas)) then
            write(field_val, *) lEx%gas_instr(gas)%firm(1:len_trim(lEx%gas_instr(gas)%firm))
            call AddDatum(csv_row, field_val, separator)
            write(field_val, *) lEx%gas_instr(gas)%model(1:len_trim(lEx%gas_instr(gas)%model))
            call AddDatum(csv_row, field_val, separator)
            write(field_val, *) lEx%measure_type(gas)
            call AddDatum(csv_row, field_val, separator)
            write(field_val, *) OutUnits(lEx%gas_instr(gas)%nsep, 1d2)
            call AddDatum(csv_row, field_val, separator)
            write(field_val, *) OutUnits(lEx%gas_instr(gas)%esep, 1d2)
            call AddDatum(csv_row, field_val, separator)
            write(field_val, *) OutUnits(lEx%gas_instr(gas)%vsep, 1d2)
            call AddDatum(csv_row, field_val, separator)
            write(field_val, *) OutUnits(lEx%gas_instr(gas)%tube_l, 1d2)
            call AddDatum(csv_row, field_val, separator)
            write(field_val, *) OutUnits(lEx%gas_instr(gas)%tube_d, 1d3)
            call AddDatum(csv_row, field_val, separator)
            write(field_val, *) OutUnits(lEx%gas_instr(gas)%tube_f, 6d4)
            call AddDatum(csv_row, field_val, separator)
            write(field_val, *) lEx%gas_instr(gas)%kw
            call AddDatum(csv_row, field_val, separator)
            write(field_val, *) lEx%gas_instr(gas)%ko
            call AddDatum(csv_row, field_val, separator)
            write(field_val, *) OutUnits(lEx%gas_instr(gas)%hpath_length, 1d2)
            call AddDatum(csv_row, field_val, separator)
            write(field_val, *) OutUnits(lEx%gas_instr(gas)%vpath_length, 1d2)
            call AddDatum(csv_row, field_val, separator)
            write(field_val, *) lEx%gas_instr(gas)%tau
            call AddDatum(csv_row, field_val, separator)
        end if
    end do
    !> '(a)', as RP writes it: list-directed output is free to lead with a
    !> blank, which then becomes part of the filename field.
    write(umd, '(a)') csv_row(1:len_trim(csv_row) - 1)

contains

    !> SI to the file's units, the error code passing through unscaled.
    !> The same function as WriteOutMetadata's; the two move together.
    real(kind = dbl) function OutUnits(si, gain)
        real(kind = dbl), intent(in) :: si
        real(kind = dbl), intent(in) :: gain

        OutUnits = merge(si * gain, si, si /= error)
    end function OutUnits
end subroutine WriteOutMetadataFcc
