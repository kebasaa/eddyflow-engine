!***************************************************************************
! write_out_metadata.f90
! ----------------------
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
! \brief       Write all results on (temporary) output files
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine WriteOutMetadata(init_string)
    use m_rp_global_var
    implicit none
    !> in/out variables
    character(*), intent(in) :: init_string
    !> local variables
    integer :: gas
!    integer :: prof
    character(LongOutstringLen) :: csv_row
    character(DatumLen) :: field_val
    include '../src_common/interfaces.inc'


    call clearstr(csv_row)
    call AddDatum(csv_row, init_string, separator)

    !> Site location and characteristics
    write(field_val, *) Metadata%lat
    call AddDatum(csv_row, field_val, separator)
    write(field_val, *) Metadata%lon
    call AddDatum(csv_row, field_val, separator)
    write(field_val, *) Metadata%alt
    call AddDatum(csv_row, field_val, separator)
    write(field_val, *) Metadata%canopy_height
    call AddDatum(csv_row, field_val, separator)
    write(field_val, *) Metadata%d
    call AddDatum(csv_row, field_val, separator)
    write(field_val, *) Metadata%z0
    call AddDatum(csv_row, field_val, separator)
    !> Acquisition setup
    write(field_val, *) Metadata%file_length
    call AddDatum(csv_row, field_val, separator)
    write(field_val, *) Metadata%ac_freq
    call AddDatum(csv_row, field_val, separator)
    !> Master sonic height and north offset
    write(field_val, *) E2Col(u)%instr%firm(1:len_trim(E2Col(u)%Instr%firm))
    call AddDatum(csv_row, field_val, separator)
    write(field_val, *) E2Col(u)%instr%model(1:len_trim(E2Col(u)%Instr%model))
    call AddDatum(csv_row, field_val, separator)
    write(field_val, *) E2Col(u)%instr%height
    call AddDatum(csv_row, field_val, separator)
    write(field_val, *) E2Col(u)%instr%wformat
    call AddDatum(csv_row, field_val, separator)
    write(field_val, *) E2Col(u)%instr%wref
    call AddDatum(csv_row, field_val, separator)
    write(field_val, *) E2Col(u)%instr%north_offset
    call AddDatum(csv_row, field_val, separator)
    write(field_val, *) E2Col(u)%instr%hpath_length * 1d2
    call AddDatum(csv_row, field_val, separator)
    write(field_val, *) E2Col(u)%instr%vpath_length  * 1d2
    call AddDatum(csv_row, field_val, separator)
    write(field_val, *) E2Col(u)%instr%tau
    call AddDatum(csv_row, field_val, separator)
    !> irgas
    do gas = co2, gas4
        if (OutVarPresent(gas)) then
            write(field_val, *) E2Col(gas)%instr%firm(1:len_trim(E2Col(gas)%Instr%firm))
            call AddDatum(csv_row, field_val, separator)
            write(field_val, *) E2Col(gas)%instr%model(1:len_trim(E2Col(gas)%Instr%model))
            call AddDatum(csv_row, field_val, separator)
            write(field_val, *) E2Col(gas)%measure_type
            call AddDatum(csv_row, field_val, separator)
            write(field_val, *) E2Col(gas)%instr%nsep * 1d2
            call AddDatum(csv_row, field_val, separator)
            write(field_val, *) E2Col(gas)%instr%esep * 1d2
            call AddDatum(csv_row, field_val, separator)
            write(field_val, *) E2Col(gas)%instr%vsep * 1d2
            call AddDatum(csv_row, field_val, separator)
            write(field_val, *) E2Col(gas)%instr%tube_l * 1d2
            call AddDatum(csv_row, field_val, separator)
            write(field_val, *) E2Col(gas)%instr%tube_d * 1d3
            call AddDatum(csv_row, field_val, separator)
            write(field_val, *) E2Col(gas)%instr%tube_f * 6d4
            call AddDatum(csv_row, field_val, separator)
            write(field_val, *) E2Col(gas)%instr%kw
            call AddDatum(csv_row, field_val, separator)
            write(field_val, *) E2Col(gas)%instr%ko
            call AddDatum(csv_row, field_val, separator)
            write(field_val, *) E2Col(gas)%instr%hpath_length * 1d2
            call AddDatum(csv_row, field_val, separator)
            write(field_val, *) E2Col(gas)%instr%vpath_length * 1d2
            call AddDatum(csv_row, field_val, separator)
            write(field_val, *) E2Col(gas)%instr%tau
            call AddDatum(csv_row, field_val, separator)
        end if
    end do

    write(umd, '(a)') csv_row(1:len_trim(csv_row) - 1)

end subroutine WriteOutMetadata
