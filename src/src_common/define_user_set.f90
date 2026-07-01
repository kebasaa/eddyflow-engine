!***************************************************************************
! define_user_set.f90
! -------------------
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
! \brief       Define "UserSet", the pre-defined set of variables  \n
!              needed for any following processing. \n
!              Variables are: u, v, w, ts, co2, h2o, ch4, gas4, tc, tc, \n
!              ti1, ti2, pi, te, pe
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine DefineUserSet(LocCol, Raw, nrow, ncol, UserSet, unrow, uncol)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    integer, intent(in) :: unrow, uncol
    type(ColType), intent(in) :: LocCol(MaxNumCol)
    real(kind = sgl), intent(in) :: Raw(nrow, ncol)
    real(kind = dbl), intent(out) :: UserSet(unrow, uncol)
    !> local variables
    integer :: j
    integer :: jj
    character(len(LocCol%label)), external :: replace


    UserCol = NullCol
    UserSet = error
    jj = 0
    do j = 1, ncol
        if (.not. IsCustomOutputColumn(LocCol(j))) cycle
        jj = jj + 1
        if (jj > uncol) then
            jj = uncol
            exit
        end if
        UserCol(jj) = LocCol(j)
        UserCol(jj)%present = .true.
        UserSet(1:unrow, jj) = Raw(1:unrow, j)
        !> Replace spaces with underscores
        UserCol(jj)%label = replace(UserCol(jj)%label, &
            ' ', '_', len(UserCol(jj)%label))
        !> Special case of 4th gas calibration reference
        if (j == Gas4CalRefCol) UserCol(jj)%var = 'cal-ref'
    end do
    NumUserVar = jj

contains

logical function IsCustomOutputColumn(col)
    type(ColType), intent(in) :: col
    character(32) :: var

    IsCustomOutputColumn = .false.
    if (col%useit) return

    var = col%var
    call lowercase(var)
    if (len_trim(var) == 0) return
    select case (trim(var))
        case ('ignore', 'not_numeric', 'none', 'flag_1', 'flag_2', &
              'agc', 'rssi')
            return
        case default
            IsCustomOutputColumn = .true.
    end select
end function IsCustomOutputColumn
end subroutine DefineUserSet
