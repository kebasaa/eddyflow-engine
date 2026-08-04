!***************************************************************************
! int_2_flag.f90
! --------------
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
! \brief       Converts integer flags into characher flags
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine Int2Flags(len)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: len

    !> sr, ar, do, al, sk, ds and tl carry one digit per variable. They are
    !> built directly as strings by the tests themselves (see PackFlagString),
    !> because the old base-10 integer packing overflows a 32-bit integer once
    !> the variable count grows. tl was the last of them to go through the
    !> integer route, and being four digits wide is what bounded the test that
    !> fills it.
    !>
    !> aa and ns are single whole-run outcomes, not per-variable strings, so
    !> they stay.
    call int2char(IntHF%aa, CharHF%aa, len)
    call int2char(IntHF%ns, CharHF%ns, len)
end subroutine Int2Flags
