!***************************************************************************
! blank_missing_values.f90
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
! \brief       Reduce everything a file can say for "no reading" to the error
!              code, in the units the file is written in.
! \author      Jonathan Muller, ETH Zurich
! \note        Four things arrive meaning the same thing, and they used to be
!              handled in four different places or not at all:
!
!              * an EMPTY field. List-directed input treats it as a null value
!                and leaves the item untouched, and the row buffer starts at the
!                error code, so this already worked - and keeps working, because
!                nothing here disturbs a value that is already error.
!              * an UNPARSEABLE token - NA, N/A, missing. The field parser maps
!                it to the error code for that field alone; before that, an
!                unparseable field cost the whole row.
!              * NaN and infinity, which are the dangerous ones: gfortran reads
!                them as perfectly good reals, so they pass every `/= error`
!                test downstream and poison the mean of anything they enter.
!              * a numeric fill - -9999 by convention, or whatever the logger
!                was configured to write, which the column can declare.
!
!              Called before any unit conversion or gain/offset, because a fill
!              value is only recognisable while it still looks like itself: a
!              -9999 in a nmol mol-1 column had already become -9.999 by the
!              time CleanUpE2Set's -300 floor saw it.
! \sa          define_all_var_set.f90, import_ascii.f90
!***************************************************************************
subroutine BlankMissingValues(Vec, N, err_value)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    real(kind = sgl), intent(inout) :: Vec(N)
    !> What this column says is missing, on top of the values that always are.
    !> `error` here means "the column declares nothing", which is every
    !> metadata file written before the key existed.
    real(kind = dbl), intent(in) :: err_value

    !> Not finite: NaN and +/-infinity. Tested first and separately, because a
    !> NaN compares equal to nothing at all - including itself - so an equality
    !> test would never catch it.
    where (Vec(1:N) /= Vec(1:N)) Vec(1:N) = sngl(error)
    where (abs(Vec(1:N)) > huge(Vec(1)) / 2e0) Vec(1:N) = sngl(error)

    !> The conventional fill needs no work: -9999 IS the engine's error code, so
    !> a file already writing it arrives correct and every `/= error` test
    !> downstream already honours it.

    !> This column's own, if it states one.
    if (err_value /= error) then
        where (Vec(1:N) == sngl(err_value)) Vec(1:N) = sngl(error)
    end if
end subroutine BlankMissingValues
