!***************************************************************************
! longest_gap_duration.f90
! ------------------------
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
! \brief       Longest consecutive gap per variable
! \author      Jonathan Muller, ETH Zurich
!
!***************************************************************************
subroutine LongestGapDuration(Set, nrow, ncol)
    use m_rp_global_var
    implicit none
    integer, intent(in) :: nrow, ncol
    real(kind = dbl), intent(in) :: Set(nrow, ncol)
    integer :: icol
    integer, external :: LongestVariableGap

    do icol = u, GHGNumVar
        if (E2Col(icol)%present) then
            Essentials%LGD(icol) = LongestVariableGap(Set(:, icol), nrow) / Metadata%ac_freq
        else
            Essentials%LGD(icol) = error
        end if
    end do
end subroutine LongestGapDuration


integer function LongestVariableGap(arr, nrow)
    use m_rp_global_var
    implicit none
    integer, intent(in) :: nrow
    real(kind = dbl), intent(in) :: arr(nrow)
    integer :: pos, run_len

    LongestVariableGap = 0
    pos = 1
    do while (pos <= nrow)
        if (arr(pos) == error) then
            run_len = 0
            do while (pos <= nrow)
                if (arr(pos) /= error) exit
                run_len = run_len + 1
                pos = pos + 1
            end do
            if (run_len > LongestVariableGap) LongestVariableGap = run_len
        else
            pos = pos + 1
        end if
    end do
end function LongestVariableGap
