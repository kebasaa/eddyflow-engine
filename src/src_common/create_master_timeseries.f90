!***************************************************************************
! create_master_timeseries.f90
! ----------------------------
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
! \brief       Create timestamp array based on start/end timestamp and step
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine CreateTimeSeries(StartTimestamp, EndTimestamp, &
    Step, RawTimeSeries, nrow, printout)
    use m_common_global_var
    implicit none
    !> In/out variables
    integer, intent(in) :: nrow
    logical, intent(in) :: printout
    type(DateType), intent(in) :: StartTimestamp
    type(DateType), intent(inout) :: EndTimestamp
    type(DateType), intent(in) :: Step
    type(DateType), intent(out) :: RawTimeSeries(nrow)
    !> in/out variables
    integer :: cnt


    if (printout) write(*, '(a)', advance = 'no') ' Creating master time series..'

    !> create master timestamps array
    RawTimeSeries(1) = StartTimestamp
    cnt = 1
    do
        if (RawTimeSeries(cnt) >= EndTimestamp) exit
        cnt = cnt + 1
        if (cnt > nrow) exit
        RawTimeSeries(cnt) = RawTimeSeries(cnt - 1) + Step
    end do

    if (printout) write(*, '(a)') ' Done.'
end subroutine CreateTimeSeries

!***************************************************************************
!
! \brief       Calculate number of periods in time series
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
integer function NumOfPeriods(StartTimestamp, EndTimestamp, Step)
    use m_common_global_var
    implicit none
    !> In/out variables
    type(DateType), intent(in) :: StartTimestamp
    type(DateType), intent(inout) :: EndTimestamp
    type(DateType), intent(in) :: Step
    !> Local variables
    integer :: cnt
    type (DateType) :: Timestamp


    Timestamp = StartTimestamp
    cnt = 1
    do
        if (Timestamp >= EndTimestamp) exit
        cnt = cnt + 1
        Timestamp = Timestamp + Step
    end do
    NumOfPeriods = cnt - 1
end function NumOfPeriods
