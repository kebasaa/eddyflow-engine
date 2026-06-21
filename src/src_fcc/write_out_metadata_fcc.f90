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
    character(16000) :: csv_row

    !> local variables
    integer :: igas
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
    write(field_val, *) lEx%instr(sonic)%hpath_length
    call AddDatum(csv_row, field_val, separator)
    write(field_val, *) lEx%instr(sonic)%vpath_length
    call AddDatum(csv_row, field_val, separator)
    write(field_val, *) lEx%instr(sonic)%tau
    call AddDatum(csv_row, field_val, separator)
    !> irgas
    do igas = ico2, igas4
        if (fcc_var_present(3 + igas)) then
            write(field_val, *) lEx%instr(igas)%firm(1:len_trim(lEx%instr(igas)%firm))
            call AddDatum(csv_row, field_val, separator)
            write(field_val, *) lEx%instr(igas)%model(1:len_trim(lEx%instr(igas)%model))
            call AddDatum(csv_row, field_val, separator)
            write(field_val, *) lEx%measure_type(3 + igas)
            call AddDatum(csv_row, field_val, separator)
            write(field_val, *) lEx%instr(igas)%nsep
            call AddDatum(csv_row, field_val, separator)
            write(field_val, *) lEx%instr(igas)%esep
            call AddDatum(csv_row, field_val, separator)
            write(field_val, *) lEx%instr(igas)%vsep
            call AddDatum(csv_row, field_val, separator)
            write(field_val, *) lEx%instr(igas)%tube_l
            call AddDatum(csv_row, field_val, separator)
            write(field_val, *) lEx%instr(igas)%tube_d
            call AddDatum(csv_row, field_val, separator)
            write(field_val, *) lEx%instr(igas)%tube_f
            call AddDatum(csv_row, field_val, separator)
            write(field_val, *) lEx%instr(igas)%kw
            call AddDatum(csv_row, field_val, separator)
            write(field_val, *) lEx%instr(igas)%ko
            call AddDatum(csv_row, field_val, separator)
            write(field_val, *) lEx%instr(igas)%hpath_length
            call AddDatum(csv_row, field_val, separator)
            write(field_val, *) lEx%instr(igas)%vpath_length
            call AddDatum(csv_row, field_val, separator)
            write(field_val, *) lEx%instr(igas)%tau
            call AddDatum(csv_row, field_val, separator)
        end if
    end do
    write(umd,*) csv_row(1:len_trim(csv_row) - 1)
end subroutine WriteOutMetadataFcc
