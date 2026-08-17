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
    integer :: spacing
    integer, external :: LongestVariableGap
    real(kind = dbl), external :: ColumnAcFreq

    do icol = u, GHGNumVar
        if (E2Col(icol)%present) then
            !> The rows between one sample of a slow column and the next are
            !> error rows by construction - nine of every ten for a 1 Hz column
            !> in a 10 Hz file - and they are not a gap in the record. Subtract
            !> that spacing, so a run no longer than one of the column's own
            !> sampling intervals reports zero rather than a permanent 0.9 s
            !> gap. A column at the file's rate has no spacing to subtract and
            !> reports exactly what it did before.
            spacing = nint(Metadata%ac_freq / ColumnAcFreq(icol)) - 1
            Essentials%LGD(icol) = &
                max(0, LongestVariableGap(Set(:, icol), nrow) - spacing) &
                / Metadata%ac_freq
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
