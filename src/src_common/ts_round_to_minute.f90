!***************************************************************************
! ts_round_to_minute.f90
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
! \brief       Rounds the given timestamp to the provided "precision"
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine tsRoundToMinute(tstamp, approx, where_to)
    use m_common_global_var
    implicit none
    !> In/out variables
    integer, intent(in) :: approx
    character(*), intent(in) :: where_to
    type (DateType), intent(inout) :: tstamp
    !> Local variables
    integer :: base
    integer :: off


    base = approx * (tstamp%minute / approx)
    off  = tstamp%minute - base
    if (off == 0) return

    if (where_to == 'later') then
        tstamp = tstamp + datetype(0, 0, 0, 0, approx - off)
    elseif (where_to == 'earlier') then
        tstamp = tstamp - datetype(0, 0, 0, 0, off)
    elseif(where_to == 'closest') then
        if (off >= approx / 2) then
            tstamp = tstamp + datetype(0, 0, 0, 0, approx - off)
        else
            tstamp = tstamp - datetype(0, 0, 0, 0, off)
        end if
    end if
end subroutine
