!***************************************************************************
! setlocale_override.f90  (Windows)
! ---------------------------------
! Copyright © 2026, ETH Zurich, Jonathan Muller
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
! \brief       Replaces the C runtime's setlocale, for libgfortran's sake.
!
!              Without newlocale, which Windows does not have, libgfortran
!              brackets every formatted READ and WRITE with
!
!                  old_locale = setlocale(LC_NUMERIC, NULL)
!                  setlocale(LC_NUMERIC, "C")
!                  ... the transfer ...
!                  setlocale(LC_NUMERIC, old_locale)
!
!              and keeps the pointer, not a copy. Under the UCRT that string
!              belongs to the locale data the "C" call replaces and frees, so
!              the restore reads freed memory. It nearly always still holds
!              "C"; once in a while the heap has reused it, and setlocale
!              faults inside ucrtbase.dll. That is the segfault
!              pwb_table_check showed about one run in eight: MAIN__,
!              StorePwbTimelagCache, LocatePwbRow, PeriodMinutes's internal
!              READ, finalize_transfer, setlocale. GNU Fortran 15.2.0,
!              MinGW-w64 UCRT; a C program doing the same three calls sees
!              the saved string overwritten.
!
!              This program never leaves the "C" locale - nothing in it calls
!              setlocale, and a C program starts in "C" - so the answer to
!              every query is "C", from storage that is never freed, a request
!              for "C" is a no-op, and a request for any other locale is
!              refused, as setlocale refuses a locale it cannot set. That is
!              everything libgfortran asks, and it no longer allocates and
!              frees locale data around each transfer either.
!
! \note        It only reaches code linked into the executable. A libgfortran
!              loaded as libgfortran-5.dll imports the runtime's own setlocale
!              and is not affected by this file.
!
!              Windows only: on POSIX libgfortran has newlocale and never
!              calls setlocale, and replacing libc's would be another matter.
! \author      Jonathan Muller
!***************************************************************************
function EddyFlowSetlocale(category, locale) bind(C, name = 'setlocale') result(current)
    use, intrinsic :: iso_c_binding, only: c_int, c_ptr, c_char, c_null_char, &
        c_null_ptr, c_loc, c_associated, c_f_pointer
    implicit none
    integer(c_int), value :: category
    type(c_ptr), value :: locale
    type(c_ptr) :: current
    !> The one locale there is. Saved, so its address outlives every call.
    character(kind = c_char), target, save :: c_name(2) = ['C', c_null_char]
    character(kind = c_char), pointer :: wanted(:)

    current = c_null_ptr
    !> Categories are small non-negative integers; anything else is an error
    !> to the real setlocale too.
    if (category < 0) return

    !> A query.
    if (.not. c_associated(locale)) then
        current = c_loc(c_name)
        return
    end if

    !> A request: "C" is honoured by doing nothing. Two characters are read
    !> only after the first is known not to end the string.
    call c_f_pointer(locale, wanted, [1])
    if (wanted(1) /= 'C') return
    call c_f_pointer(locale, wanted, [2])
    if (wanted(2) /= c_null_char) return
    current = c_loc(c_name)
end function EddyFlowSetlocale
