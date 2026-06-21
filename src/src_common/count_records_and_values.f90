!***************************************************************************
! count_records_and_values.f90
! ----------------------------
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
! \brief       Count available records or valid values in a dataset
! \author      Jonathan Muller, ETH Zurich
!
!***************************************************************************
integer function CountRecordsAndValues(Set, nrow, ncol, var1, var2)
    use m_common_global_var
    implicit none
    integer, intent(in) :: nrow, ncol
    integer, optional, intent(in) :: var1, var2
    real(kind = dbl), intent(in) :: Set(nrow, ncol)

    if (.not. present(var1)) then
        CountRecordsAndValues = count(any(Set(:,1:ncol) /= error, dim=2))
    else if (.not. present(var2)) then
        CountRecordsAndValues = count(Set(:, var1) /= error)
    else
        CountRecordsAndValues = &
            count(Set(:, var1) /= error .and. Set(:, var2) /= error)
    end if
end function CountRecordsAndValues
