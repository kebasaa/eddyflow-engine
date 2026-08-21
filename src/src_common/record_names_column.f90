!***************************************************************************
! record_names_column.f90
! -----------------------
! Copyright (c) 2026-    , ETH Zurich, Jonathan Muller
!
! This file is part of EddyFlow(R).
!
! EddyFlow (TM) is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version. You should have received a copy
! of the GNU General Public License along with EddyFlow (R). If not,
! see <http://www.gnu.org/licenses/>.
!
! EddyFlow(R) contains additional Open Source Components. The licenses
! and/or notices these Components can be found in the file LIBRARIES.txt.
!
! EddyFlow(R) is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
!***************************************************************************
!
! \brief       Whether a measurement record still describes the column it
!              names.
! \author      Jonathan Muller
! \note        A record is a pair: a raw column number, and what that column
!              measures. The metadata is the authority on the second half, and
!              the two can disagree - re-declare a column in the Raw File
!              Description and the record naming it is left behind, pointing
!              at a column that now measures something else entirely.
!
!              Such a record has to be inert. Honoured, two things followed. A
!              leftover diag_72 record on a column since re-declared AGC still
!              counted as a diagnostic, so a project with a real diagnostic
!              elsewhere had two records competing for the one slot and
!              MetadataFileValidation refused the run - over a record the
!              interface does not show and the user cannot remove. And where it
!              did not collide it won the slot instead, and the engine decoded
!              a signal-strength percentage as a bitfield.
!
!              `ignore` and `not_numeric` are the metadata saying the column
!              holds nothing usable, and they are the case this generalises:
!              they were already exempt, for exactly this reason.
!
!              Names are compared case-insensitively, so a metadata file
!              written by hand or by another tool is read the way it was meant.
! \sa          ColumnIsSelectable in define_used_variables.f90
!***************************************************************************
logical function RecordNamesColumn(col_var, record_var)
    implicit none
    !> in/out variables
    character(*), intent(in) :: col_var
    character(*), intent(in) :: record_var
    !> local variables
    character(32) :: col
    character(32) :: rec

    RecordNamesColumn = .false.

    col = adjustl(col_var)
    rec = adjustl(record_var)
    call lowercase(col)
    call lowercase(rec)

    if (len_trim(col) == 0 .or. len_trim(rec) == 0) return
    if (trim(col) == 'ignore' .or. trim(col) == 'not_numeric') return

    !> The anemometer diagnostic is the one measurement whose record slug is
    !> not the metadata's own spelling of it: the interface writes `diag_anem`
    !> where the Raw File Description says `anemometer_diagnostic`. Both are
    !> accepted on both sides rather than one being declared correct, because
    !> project files carrying either are already in the field.
    if (trim(col) == 'diag_anem') col = 'anemometer_diagnostic'
    if (trim(rec) == 'diag_anem') rec = 'anemometer_diagnostic'

    RecordNamesColumn = trim(col) == trim(rec)
end function RecordNamesColumn
