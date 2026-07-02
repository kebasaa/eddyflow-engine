!***************************************************************************
! define_vars.f90
! ---------------
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
! \brief       Map instrument columns to E2Col/UserCol arrays
! \author      Jonathan Muller, ETH Zurich
!
!***************************************************************************
subroutine DefineVars(LocCol, ncol, uncol)
    use m_rp_global_var
    implicit none
    integer, intent(in) :: ncol, uncol
    type(ColType), intent(in) :: LocCol(MaxNumCol)
    integer :: idx, usr_cnt
    character(len(LocCol%label)), external :: replace

    E2Col = NullCol

    ! Pass 1: sonic wind components (master sonic only)
    do idx = 1, ncol
        if (.not. LocCol(idx)%Instr%master_sonic) cycle
        select case (trim(LocCol(idx)%var))
            case ('u');   E2Col(u)  = LocCol(idx);  E2Col(u)%present  = .true.
            case ('v');   E2Col(v)  = LocCol(idx);  E2Col(v)%present  = .true.
            case ('w');   E2Col(w)  = LocCol(idx);  E2Col(w)%present  = .true.
            case ('ts');  E2Col(ts) = LocCol(idx);  E2Col(ts)%present = .true.
            case ('sos'); E2Col(ts) = LocCol(idx);  E2Col(ts)%present = .true.
        end select
    end do

    ! Pass 2: gas species and auxiliary met (useit required)
    do idx = 1, ncol
        if (.not. LocCol(idx)%useit) cycle
        select case (trim(LocCol(idx)%var))
            case ('co2');     E2Col(co2)  = LocCol(idx);  E2Col(co2)%present  = .true.
            case ('h2o');     E2Col(h2o)  = LocCol(idx);  E2Col(h2o)%present  = .true.
            case ('ch4');     E2Col(ch4)  = LocCol(idx);  E2Col(ch4)%present  = .true.
            case ('n2o');     E2Col(gas4) = LocCol(idx);  E2Col(gas4)%present = .true.
            case ('ts');      E2Col(ts)   = LocCol(idx);  E2Col(ts)%present   = .true.
            case ('cell_t');  E2Col(tc)   = LocCol(idx);  E2Col(tc)%present   = .true.
            case ('int_t_1'); E2Col(ti1)  = LocCol(idx);  E2Col(ti1)%present  = .true.
            case ('int_t_2'); E2Col(ti2)  = LocCol(idx);  E2Col(ti2)%present  = .true.
            case ('int_p');   E2Col(pi)   = LocCol(idx);  E2Col(pi)%present   = .true.
            case ('air_t');   E2Col(te)   = LocCol(idx);  E2Col(te)%present   = .true.
            case ('air_p');   E2Col(pe)   = LocCol(idx);  E2Col(pe)%present   = .true.
        end select
    end do

    ! Pass 3: collect user columns (not mapped to standard slots above)
    UserCol = NullCol
    usr_cnt = 0
    do idx = 1, ncol
        if (.not. IsCustomOutputColumn(LocCol(idx))) cycle
        usr_cnt = usr_cnt + 1
        if (usr_cnt > uncol) then
            usr_cnt = uncol
            exit
        end if
        UserCol(usr_cnt) = LocCol(idx)
        UserCol(usr_cnt)%present = .true.
        UserCol(usr_cnt)%label = replace(UserCol(usr_cnt)%label, &
            ' ', '_', len(UserCol(usr_cnt)%label))
        if (idx == Gas4CalRefCol) UserCol(usr_cnt)%var = 'cal-ref'
    end do
    NumUserVar = usr_cnt

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
end subroutine DefineVars
