!***************************************************************************
! test_spike_detection_mauder_13.f90
! ----------------------------------
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
! \brief       Detects and replace spikes by linear interpolation \n
!              Hard-flags file for too many spikes
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine TestSpikeDetectionMauder13(Set, N, printout)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    logical, intent(in) :: printout
    real(kind = dbl), intent(inout) :: Set(N, E2NumVar)
    !> local variables
    integer :: max_pass = 10
    integer :: passes = 0
    integer :: i, ii, lN
    integer :: j
    integer :: k
    integer :: cnt = 0
    integer :: prev_good
    integer :: run_start
    integer :: run_last
    integer :: nspikes(E2NumVar)
    integer :: nspikes_sng(E2NumVar)
    integer :: tot_spikes(E2NumVar)
    integer :: tot_spikes_sng(E2NumVar)
    integer :: hflags(u:lastGas)
    real(kind = dbl) :: MAD
    real(kind = dbl) :: medx
    real(kind = dbl) :: zlim
    real(kind = dbl) :: m
    real(kind = dbl) :: q
    logical :: again
    logical :: IsSpike(N, E2NumVar)
    logical :: new_spike
    real(kind = dbl), allocatable :: tmpx(:)


    if (printout) write(*, '(a)', advance = 'no') '   Spike detection/removal test..'
    if (printout) write(ulog, '(a)', advance = 'no') '   Spike detection/removal test..'
    zlim = 7d0
    passes = 0
    nspikes = 0
    nspikes_sng = 0
    tot_spikes = 0
    tot_spikes_sng = 0
    IsSpike = .false.
100 continue
    passes = passes + 1
    do j = u, pe
        if (E2Col(j)%present) then
            cnt = 0
            !> Define array tmpx made of Set(:, j) stripped of potential error values.
            lN = N - count(Set(:, j) == error)
            allocate(tmpx(lN))
            ii = 0
            do i = 1, N
                if (Set(i, j) == error) cycle
                ii = ii + 1
                tmpx(ii) = Set(i, j)
            end do

            !> Calculate median (medx) and MAD
            call median(tmpx, lN, medx)
            tmpx = dabs(tmpx - medx)
            call median(tmpx, lN, MAD)
            deallocate(tmpx)

            !> Spike detection (and replacement if needed)
            !>
            !> In samples, as in the Vickers and Mahrt variant. The error rows
            !> were skipped without breaking the run, so `cnt` counted samples
            !> while the flagging and the interpolation addressed the `cnt`
            !> rows before `i`: with an error row inside or just before the
            !> run, the line was drawn from the error code, and a spike next to
            !> one filtered sample became about (x - 9999) / 2 - a value the
            !> later passes no longer looked at, carried into rotation and
            !> every covariance. A column slower than the file has error rows
            !> between all its samples, so there it happened at every spike.
            !>
            !> The run's rows are tracked alongside, the line is drawn from the
            !> last good sample before it, and only the sample rows are
            !> written. With no good sample before the run there is nothing to
            !> draw from, and the spike is flagged but left as it is.
            prev_good = 0
            if (Set(1, j) /= error) prev_good = 1
            run_start = 0
            run_last = 0
            if (RPsetup%filter_sr) then
                do i = 2, N-1
                    if (Set(i, j) == error) cycle

                    if(Set(i, j) > medx + (zlim * MAD / 0.6745d0) .or. &
                        Set(i, j) < medx - (zlim * MAD / 0.6745d0)) then
                        if (cnt == 0) run_start = i
                        run_last = i
                        cnt = cnt + 1
                    else
!                        if ((cnt /= 0) .and. (cnt <= sr%num_spk)) then
                        if (cnt /= 0) then
                            !> check whether it was a spike already,
                            !> if not increment the number of spikes found.
                            !> The old loop exited after its first pass, so
                            !> it only ever looked at the row before i.
                            new_spike = .not. IsSpike(run_last, j)
                            if (new_spike) then
                                nspikes(j) = nspikes(j) + 1
                                nspikes_sng(j) = nspikes_sng(j) + cnt
                            end if
                            !> update spike flags, on the samples only
                            do k = run_start, i - 1
                                if (Set(k, j) /= error) IsSpike(k, j) = .true.
                            enddo
                            !> replace with linear interpolation, from the
                            !> last good sample behind the run
                            if (prev_good /= 0) then
                                m = (Set(i, j) - Set(prev_good, j)) / dble(i - prev_good)
                                q = Set(prev_good, j)
                                do k = run_start, i - 1
                                    if (Set(k, j) /= error) &
                                        Set(k, j) = m * dble(k - prev_good) + q
                                end do
                            end if
                            cnt = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! In Mauder's method there is no concept of max num of outliers to define a spike
! so these two lines need eliminating.
!                        else if (cnt > sr%num_spk) then
!                            cnt = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                        end if
                        prev_good = i
                    end if
                end do
            else
                !> Detection only, counted the same way, so a flag-only run
                !> reports the spikes a despiking run would remove.
                do i = 2, N-1
                    if (Set(i, j) == error) cycle

                    if(Set(i, j) > medx + (zlim * MAD / 0.6745d0) .or. &
                        Set(i, j) < medx - (zlim * MAD / 0.6745d0)) then
                        if (cnt == 0) run_start = i
                        run_last = i
                        cnt = cnt + 1
                    else
                        if (cnt /= 0) then
                            !> check whether it was a spike already, if not increment the
                            !> number of spikes found
                            new_spike = .not. IsSpike(run_last, j)
                            if (new_spike) then
                                nspikes(j) = nspikes(j) + 1
                                nspikes_sng(j) = nspikes_sng(j) + cnt
                            end if
                            !> update spike flags, on the samples only
                            do k = run_start, i - 1
                                if (Set(k, j) /= error) IsSpike(k, j) = .true.
                            enddo
                            cnt = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! In Mauder's method there is no concept of max num of outliers to define a spike
! so these two lines need eliminating.
!                        else if (cnt > sr%num_spk) then
!                            cnt = 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        end if
                        prev_good = i
                    end if
                end do
            end if
        end if
    end do

    !> accumulates the number of spikes
    tot_spikes = tot_spikes + nspikes
    tot_spikes_sng = tot_spikes_sng + nspikes_sng

    !> if spikes have been detected (and removed) does the whole process again
    again = .false.
    do j = u, Pe
        if (E2Col(j)%present) then
            if (nspikes(j) /= 0) again = .true.
            nspikes(j) = 0
            nspikes_sng(j) = 0
        end if
    end do
    if (again .and. (passes < max_pass)) then
        goto 100
    end if

    !> hflags the variable if nspikes is larger than a prescribed threshold.
    !>
    !> Every gas slot, for the same reason as the Vickers 97 variant: an absent
    !> slot keeps its 9, so four gases pack the string they always did, while a
    !> fifth now gets a real outcome instead of "test not performed".
    hflags = 9
    do j = u, lastGas
        if (E2Col(j)%present) then
            if(100.d0 * (dble(tot_spikes(j)) / dble(N)) >= sr%hf_lim) then
                hflags(j) = 1
            else
                hflags(j) = 0
            end if
        end if
    end do

    !> Pack one digit per variable into the flag string
    call PackFlagString(hflags(u:lastGas), GHGNumVar, CharHF%sr)

    !> Write on output variable
    if (.not. RPsetup%filter_sr) tot_spikes_sng(u:pe) = 0
    where (E2Col(u:pe)%present) 
        Essentials%e2spikes(u:pe) = tot_spikes(u:pe)
        Essentials%m_despiking(u:pe) = tot_spikes_sng(u:pe)
    elsewhere
        Essentials%e2spikes(u:pe) = ierror
        Essentials%m_despiking(u:pe) = ierror
    endwhere
    
    if (printout) write(*,'(a)') ' Done.'
    if (printout) write(ulog,'(a)') ' Done.'
end subroutine TestSpikeDetectionMauder13
