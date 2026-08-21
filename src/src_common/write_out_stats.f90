!***************************************************************************
! write_out_stats.f90
! -------------------
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
! \brief       Write statistics of sensitive variables on output
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine WriteOutStats(unt, LocStats, string, N)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: unt
    type(StatsType), intent(in) :: LocStats
    integer, intent(in) :: N
    character(*), intent(in) :: string
    !local variables
    integer :: j = 0
    integer :: i = 0
    integer :: k = 0
    integer :: slots(E2NumVar)
    integer :: nslots = 0
    integer :: gasslots(GHGNumVar)
    integer :: ngas = 0
    character(LongOutstringLen) :: dataline
    character(DatumLen) :: datum = ''
    include '../src_common/interfaces.inc'


    !> The layout both this writer and the header in InitOutFiles_rp walk.
    !> They used to agree only because `u, pe` happened to enumerate the
    !> twelve names the header lists, which stopped being true when E2NumVar
    !> grew from 14 to 102.
    call StatsLayoutSlots(slots, nslots)
    ngas = 0
    do i = 1, nslots
        if (slots(i) >= firstGas .and. slots(i) <= lastGas) then
            ngas = ngas + 1
            gasslots(ngas) = slots(i)
        end if
    end do

    call clearstr(dataline)
    !> add file info
    call AddDatum(dataline, string(1:len_trim(string)), separator)
    call WriteDatumInt(N, datum, EddyFlowProj%err_label)
    call AddDatum(dataline, datum, separator)
    !> add mean values
    do i = 1, nslots
        j = slots(i)
        if (E2Col(j)%present) then
            call WriteDatumFloat(LocStats%Mean(j), datum, EddyFlowProj%err_label)
            call AddDatum(dataline, datum, separator)
        else
            call AddDatum(dataline, EddyFlowProj%err_label(1:len_trim(EddyFlowProj%err_label)), separator)
        end if
    end do
    !> add wind direction
    call WriteDatumFloat(LocStats%wind_dir, datum, EddyFlowProj%err_label)
    call AddDatum(dataline, datum, separator)
    !> add variances
    do i = 1, nslots
        j = slots(i)
        if (E2Col(j)%present) then
            call WriteDatumFloat(LocStats%Cov(j, j), datum, EddyFlowProj%err_label)
            call AddDatum(dataline, datum, separator)
        else
            call AddDatum(dataline, EddyFlowProj%err_label(1:len_trim(EddyFlowProj%err_label)), separator)
        end if
    end do

    !> add u-covariances
    if (E2Col(u)%present .and. E2Col(v)%present) then
        call WriteDatumFloat(LocStats%Cov(u, v), datum, EddyFlowProj%err_label)
        call AddDatum(dataline, datum, separator)
    else
        call AddDatum(dataline, EddyFlowProj%err_label(1:len_trim(EddyFlowProj%err_label)), separator)
    end if
    if (E2Col(u)%present .and. E2Col(w)%present) then
        call WriteDatumFloat(LocStats%Cov(u, w), datum, EddyFlowProj%err_label)
        call AddDatum(dataline, datum, separator)
    else
        call AddDatum(dataline, EddyFlowProj%err_label(1:len_trim(EddyFlowProj%err_label)), separator)
    end if
    if (E2Col(u)%present .and. E2Col(ts)%present) then
        call WriteDatumFloat(LocStats%Cov(u, ts), datum, EddyFlowProj%err_label)
        call AddDatum(dataline, datum, separator)
    else
        call AddDatum(dataline, EddyFlowProj%err_label(1:len_trim(EddyFlowProj%err_label)), separator)
    end if
    !> One covariance per configured gas. These were four arms
    !> naming co2/h2o/ch4/gas4; the header had the same four, so
    !> both stop at the fourth gas or neither does.
    do k = 1, ngas
        if (E2Col(u)%present .and. E2Col(gasslots(k))%present) then
            call WriteDatumFloat(LocStats%Cov(u, gasslots(k)), datum, EddyFlowProj%err_label)
            call AddDatum(dataline, datum, separator)
        else
            call AddDatum(dataline, EddyFlowProj%err_label(1:len_trim(EddyFlowProj%err_label)), separator)
        end if
    end do

    !> add remaining v-covariances
    if (E2Col(v)%present .and. E2Col(w)%present) then
        call WriteDatumFloat(LocStats%Cov(v, w), datum, EddyFlowProj%err_label)
        call AddDatum(dataline, datum, separator)
    else
        call AddDatum(dataline, EddyFlowProj%err_label(1:len_trim(EddyFlowProj%err_label)), separator)
    end if
    if (E2Col(v)%present .and. E2Col(ts)%present) then
        call WriteDatumFloat(LocStats%Cov(v, ts), datum, EddyFlowProj%err_label)
        call AddDatum(dataline, datum, separator)
    else
        call AddDatum(dataline, EddyFlowProj%err_label(1:len_trim(EddyFlowProj%err_label)), separator)
    end if
    do k = 1, ngas
        if (E2Col(v)%present .and. E2Col(gasslots(k))%present) then
            call WriteDatumFloat(LocStats%Cov(v, gasslots(k)), datum, EddyFlowProj%err_label)
            call AddDatum(dataline, datum, separator)
        else
            call AddDatum(dataline, EddyFlowProj%err_label(1:len_trim(EddyFlowProj%err_label)), separator)
        end if
    end do

    !> add remaining w-covariances
    if (E2Col(w)%present .and. E2Col(ts)%present) then
        call WriteDatumFloat(LocStats%Cov(w, ts), datum, EddyFlowProj%err_label)
        call AddDatum(dataline, datum, separator)
    else
        call AddDatum(dataline, EddyFlowProj%err_label(1:len_trim(EddyFlowProj%err_label)), separator)
    end if
    do k = 1, ngas
        if (E2Col(w)%present .and. E2Col(gasslots(k))%present) then
            call WriteDatumFloat(LocStats%Cov(w, gasslots(k)), datum, EddyFlowProj%err_label)
            call AddDatum(dataline, datum, separator)
        else
            call AddDatum(dataline, EddyFlowProj%err_label(1:len_trim(EddyFlowProj%err_label)), separator)
        end if
    end do
    if (E2Col(w)%present .and. E2Col(tc)%present) then
        call WriteDatumFloat(LocStats%Cov(w, tc), datum, EddyFlowProj%err_label)
        call AddDatum(dataline, datum, separator)
    else
        call AddDatum(dataline, EddyFlowProj%err_label(1:len_trim(EddyFlowProj%err_label)), separator)
    end if
    if (E2Col(w)%present .and. E2Col(pi)%present) then
        call WriteDatumFloat(LocStats%Cov(w, pi), datum, EddyFlowProj%err_label)
        call AddDatum(dataline, datum, separator)
    else
        call AddDatum(dataline, EddyFlowProj%err_label(1:len_trim(EddyFlowProj%err_label)), separator)
    end if
    if (E2Col(w)%present .and. E2Col(te)%present) then
        call WriteDatumFloat(LocStats%Cov(w, te), datum, EddyFlowProj%err_label)
        call AddDatum(dataline, datum, separator)
    else
        call AddDatum(dataline, EddyFlowProj%err_label(1:len_trim(EddyFlowProj%err_label)), separator)
    end if
    if (E2Col(w)%present .and. E2Col(pe)%present) then
        call WriteDatumFloat(LocStats%Cov(w, pe), datum, EddyFlowProj%err_label)
        call AddDatum(dataline, datum, separator)
    else
        call AddDatum(dataline, EddyFlowProj%err_label(1:len_trim(EddyFlowProj%err_label)), separator)
    end if

    !> add standard deviations
    do i = 1, nslots
        j = slots(i)
        if (E2Col(j)%present) then
            call WriteDatumFloat(LocStats%StDev(j), datum, EddyFlowProj%err_label)
            call AddDatum(dataline, datum, separator)
        else
            call AddDatum(dataline, EddyFlowProj%err_label(1:len_trim(EddyFlowProj%err_label)), separator)
        end if
    end do
    !> add skewness and kurtosis
    do i = 1, nslots
        j = slots(i)
        if (E2Col(j)%present) then
            call WriteDatumFloat(LocStats%Skw(j), datum, EddyFlowProj%err_label)
            call AddDatum(dataline, datum, separator)
        else
            call AddDatum(dataline, EddyFlowProj%err_label(1:len_trim(EddyFlowProj%err_label)), separator)
        end if
    end do
    do i = 1, nslots
        j = slots(i)
        if (E2Col(j)%present) then
            call WriteDatumFloat(LocStats%Kur(j), datum, EddyFlowProj%err_label)
            call AddDatum(dataline, datum, separator)
        else
            call AddDatum(dataline, EddyFlowProj%err_label(1:len_trim(EddyFlowProj%err_label)), separator)
        end if
    end do
    write(unt, '(a)') dataline(1:len_trim(dataline) - 1)
end subroutine WriteOutStats
