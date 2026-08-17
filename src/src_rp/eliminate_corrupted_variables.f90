!***************************************************************************
! eliminate_corrupted_variables.f90
! ---------------------------------
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
! \brief       If a variable is missing more of its own data than its
!              instrument is allowed to, set it as not present.
! \author      Gerardo Fratini
! \note        "Its own data" rather than "the file's rows": an instrument
!              slower than the row rate cannot fill every row, and counting the
!              rows it was never going to write as missing dropped the column
!              outright - 95 % missing for a 1 Hz analyser in a 20 Hz file, on
!              a record that is complete for its own rate.
! \sa          column_sampling.f90
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine EliminateCorruptedVariables(LocSet, nrow, ncol, skip_period, logout)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    real(kind = dbl), intent(in) :: LocSet(nrow, ncol)
    logical, intent(in)  :: logout
    logical, intent(out) :: skip_period
    !> local variables
    integer :: i
    integer :: expected
    integer :: full
    real(kind = dbl) :: freq
    real(kind = dbl), external :: ColumnAcFreq
    real(kind = dbl), external :: ColumnMaxLack


    if (logout) write(*,'(a)', advance = 'no') '  Verifying time series integrity..'
    if (logout) write(ulog,'(a)', advance = 'no') '  Verifying time series integrity..'

    do i = 1, ncol
        freq = ColumnAcFreq(i)
        !> What this column should have produced in the rows at hand, and in a
        !> whole averaging period. At the file's own rate expected is nrow and
        !> full is MaxPeriodNumRecords, so this is arithmetically the test it
        !> replaces: counting what is missing from what was expected is the
        !> same as counting error rows, until the two rates differ.
        expected = nint(dble(nrow) * freq / Metadata%ac_freq)
        full = nint(RPsetup%avrg_len * 60d0 * freq)
        if (expected - count(LocSet(:, i) /= error) &
            > full * ColumnMaxLack(i)/1d2) E2Col(i) = NullCol
    end do

    skip_period = .false.
    if ((.not. E2Col(u)%present) .or. &
        (.not. E2Col(v)%present) .or. &
        (.not. E2Col(w)%present)) skip_period = .true.

    if (logout) write(*,'(a)') ' Done.'
    if (logout) write(ulog,'(a)') ' Done.'
end subroutine EliminateCorruptedVariables
