!***************************************************************************
! retrieve_ex_var_by_timestamp.f90
! --------------------------------
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
! \brief       Retrieve "essentials" information from file, based on
!              timestamp provided on input
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine RetrieveExVarsByTimestamp(unt, Timestamp, lEx, endReached, skip)
    use m_fx_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: unt
    type(DateType), intent(in) :: Timestamp
    logical, intent(out) :: endReached
    logical, intent(out) :: skip
    type(ExType), intent(out) :: lEx
    !> Local variables
    type(DateType) :: ExTimestamp
    logical :: EndOfFileReached
    logical :: ValidRecord

    skip = .false.
    endReached = .false.
    do
        call ReadExRecord('', unt, -1, lEx, ValidRecord, EndOfFileReached)

        !> If end of files was reached, exit routine with error flag
        if (EndOfFileReached) then
            endReached = .true.
            skip = .true.
            return
        end if

        !> If timestamp matches, exit routine (with error flag if the case)
        call DateTimeToDateType(lEx%end_date, lEx%end_time, ExTimestamp)
        if (ExTimestamp == Timestamp) then
            if (.not. ValidRecord) skip = .true.
            return
        end if

        !> If timestamp exceeds the one looked for, backwards essentials unit
        !> and exit with error code
        if (ExTimestamp >= Timestamp) then
            backspace(unt)
            skip = .true.
            return
        end if
    end do
end subroutine RetrieveExVarsByTimestamp
