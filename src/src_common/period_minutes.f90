!***************************************************************************
! period_minutes.f90
! -------------------
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
! \brief       Minutes since 1970-01-01, from 'yyyy-mm-dd' and 'HH:MM'.
!
! \details     A real time axis, so that interpolation across a gap in the
!              raw files spans the time the gap actually lasted rather than
!              the number of rows it happens to occupy.
!
!              Moved here from pwb_timelag_handle.f90 (RP-only) so the FCC
!              side (pfd_handle.f90's whole-run post-flux despiking pass)
!              can share the same conversion rather than keep a second copy
!              that could drift from it.
!
! \author      Jonathan Muller, ETH Zurich
!***************************************************************************
integer(8) function PeriodMinutes(date, time)
    character(*), intent(in) :: date, time
    integer :: y, m, d, hh, mm, ios
    integer(8) :: era, yoe, doy, doe, days

    PeriodMinutes = 0
    read(date(1:4), '(i4)', iostat=ios) y
    if (ios /= 0) return
    read(date(6:7), '(i2)', iostat=ios) m
    if (ios /= 0) return
    read(date(9:10), '(i2)', iostat=ios) d
    if (ios /= 0) return
    read(time(1:2), '(i2)', iostat=ios) hh
    if (ios /= 0) return
    read(time(4:5), '(i2)', iostat=ios) mm
    if (ios /= 0) return

    !> Days from civil (Howard Hinnant): exact for any proleptic Gregorian
    !> date, and integer throughout.
    if (m <= 2) y = y - 1
    era = int(y, 8) / 400_8
    if (int(y, 8) < 0_8 .and. mod(int(y, 8), 400_8) /= 0_8) era = era - 1_8
    yoe = int(y, 8) - era * 400_8
    if (m > 2) then
        doy = (153_8 * int(m - 3, 8) + 2_8) / 5_8 + int(d, 8) - 1_8
    else
        doy = (153_8 * int(m + 9, 8) + 2_8) / 5_8 + int(d, 8) - 1_8
    end if
    doe = yoe * 365_8 + yoe / 4_8 - yoe / 100_8 + doy
    days = era * 146097_8 + doe - 719468_8
    PeriodMinutes = days * 1440_8 + int(hh, 8) * 60_8 + int(mm, 8)
end function PeriodMinutes
