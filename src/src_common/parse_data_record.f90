!***************************************************************************
! parse_data_record.f90
! ---------------------
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
! \brief       Split one data record into its hot columns, a field at a time.
! \author      Gerardo Fratini
! \note        Extracted from ImportAsciiWithText, which has always done this;
!              ImportAscii now falls back to it when its whole-record read
!              fails, so the two cannot drift apart.
!
!              The point of parsing a field at a time is what happens to a field
!              that will not parse: it becomes the error code and the rest of the
!              record is kept. A single list-directed read over the whole record
!              cannot do that - one `NA` in one gas column fails the statement,
!              and the row is discarded entire, taking the wind data with it.
! \sa          import_ascii.f90, import_ascii_with_text.f90
!***************************************************************************
subroutine ParseDataRecord(dataline, LocCol, Vec, ncol)
    use m_common_global_var
    implicit none
    !> in/out variables
    !> One value per RAW file column, in file order, ignored columns included -
    !> the same shape the whole-record read fills, so a caller can swap one for
    !> the other and compact afterwards exactly as it already does. Filling hot
    !> slots here instead silently shifts every column past the first ignored
    !> one, which is what it did until this was measured.
    integer, intent(in) :: ncol
    character(*), intent(in) :: dataline
    type(ColType), intent(in) :: LocCol(MaxNumCol)
    real(kind = sgl), intent(out) :: Vec(ncol)
    !> local variables
    integer :: j
    integer :: intsep
    integer :: io_status
    character(len(dataline)) :: rest
    character(DatumLen) :: datum

    Vec(1:ncol) = error
    rest = dataline

    !> Multiple separators collapse only for spaces, as in the caller this came
    !> from: with commas, two in a row is an empty field and means "no reading",
    !> which is not the same as one separator.
    if (FileInterpreter%separator == ' ') &
        call StripConsecutiveChar(rest, FileInterpreter%separator)

    do j = 1, min(NumCol, ncol)
        intsep = index(rest, FileInterpreter%separator)
        if (intsep == 0) intsep = len_trim(rest) + 1
        if (len_trim(rest) == 0) exit
        datum = rest(1:intsep - 1)
        rest = rest(intsep + 1: len_trim(rest))
        if (LocCol(j)%var == 'ignore' .or. LocCol(j)%var == 'not_numeric') cycle

        !> Two Gill diagnostic words are hexadecimal, and the metadata declares
        !> the column numeric. Translated here rather than rejected.
        if (EddyFlowProj%col(E2NumVar + DiagAnem) == j &
            .and. (index(EddyFlowProj%master_sonic, 'wm') /= 0 &
            .or. index(EddyFlowProj%master_sonic, 'hs') /= 0)) then
            if (trim(datum) == '0A') datum = '10'
            if (trim(datum) == '0B') datum = '11'
        end if

        !> An empty field and an unparseable one mean the same thing, and both
        !> land here: NA, N/A, missing, or whatever else a logger writes in
        !> words. One field's worth of data is lost, not the record.
        read(datum, *, iostat = io_status) Vec(j)
        if (io_status /= 0) Vec(j) = error
    end do
end subroutine ParseDataRecord
