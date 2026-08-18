!*******************************************************************************
! test_spike_detection_vickers_97.f90
! -----------------------------------
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
!*******************************************************************************
!
! \brief       Detects  and count spikes, and replace by \n
!              linear interpolation if requested. \n
!              Hard-flags file for too many spikes
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!*******************************************************************************
subroutine TestSpikeDetectionVickers97(Set, N, printout)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    logical, intent(in) :: printout
    real(kind = dbl), intent(inout) :: Set(N, E2NumVar)
    !> local variables
    integer :: win_len
    integer, parameter :: step = 100 !window advancement in samples
    integer :: max_pass = 10
    real(kind = dbl), parameter :: lim_step = 0.1d0 !increase of inliers range
    integer :: i = 0
    integer :: j = 0
    integer :: k = 0
    integer :: nn = 0
    integer :: imin = 0
    integer :: imax = 0
    integer :: cnt = 0
    !> The open run's row extent. `cnt` counts spiky samples; these locate them,
    !> because the flagging and the interpolation address rows and a column
    !> slower than the file has error rows in between. `prev_good` is the last
    !> good sample before the run, which the interpolation starts from.
    integer :: run_start = 0
    integer :: run_last = 0
    integer :: prev_good = 0
    !> How many samples the column really carries, for the spike percentage
    integer :: nsample = 0
    integer :: passes = 0
    integer :: wdw_num = 0
    integer :: wdw = 0
    integer :: nspikes(E2NumVar)
    integer :: nspikes_sng(E2NumVar)
    integer :: tot_spikes(E2NumVar)
    integer :: tot_spikes_sng(E2NumVar)
    integer :: hflags(u:lastGas)
    real(kind = dbl) :: Mean(E2NumVar) = 0.d0
    real(kind = dbl) :: StDev(E2NumVar) = 0.d0
    real(kind = dbl) :: LocMean(N, E2NumVar)
    real(kind = dbl) :: LocStDev(N, E2NumVar)
    logical :: IsSpike(N, E2NumVar)
    real(kind = dbl) :: m = 0.d0
    real(kind = dbl) :: q = 0.d0
    real(kind = dbl) :: adv_lim(E2NumVar)
    real(kind = dbl), allocatable :: XX(:, :)
    logical :: again = .false.
    logical :: new_spike
    real(kind = dbl), external :: ColumnAcFreq


    if (printout) write(*, '(a)', advance = 'no') '   Spike detection/removal test..'
    if (printout) write(ulog, '(a)', advance = 'no') '   Spike detection/removal test..'

    if (.not. RPsetup%filter_sr) max_pass = 1

    win_len = RPsetup%avrg_len / 6
    if (win_len == 0) win_len = 1
    nn = idint(dble(win_len) * Metadata%ac_freq * 60.d0) !> win length in samples
    wdw_num = idint(dble(N - nn) / 1d2) + 1  !> number of wins for current file

    !> initializations
    allocate(XX(nn, E2NumVar))
    XX = error

    LocMean = 0d0
    LocStDev = 0d0
    passes = 0
    nspikes = 0
    nspikes_sng = 0
    tot_spikes = 0
    tot_spikes_sng = 0

    !> Set different threshold for different variables.
    !> Specifically, w and every configured gas have their own thresholds
    adv_lim(u:pe) = sr%lim_u
    adv_lim(w)    = sr%lim_w
    adv_lim(firstGas:lastGas) = sr%lim_gas(firstGas:lastGas)

    IsSpike = .false.
    !> main cycle, looping over the moving window
100 continue
    passes = passes + 1
    do wdw = 1, wdw_num
        !> pick up the dataset from Set for current win
        do i = 1, nn
            where(E2Col(:)%present)
                XX(i, :) = Set(i + step * (wdw - 1), :)
            end where
        end do
        !> define min and max of central points
        imin = nn/2 - step/2 + step * (wdw - 1)
        imax = nn/2 + step/2 -1 + step * (wdw - 1)

        !> Window mean values
        call AverageNoError(XX, size(XX, 1), size(XX, 2), Mean, error)

        !> Window standard deviations
        call StDevNoError(XX, size(XX, 1), size(XX, 2), StDev, error)

        !> stick window values only to central elements
        do i = imin, imax
            where(E2Col(:)%present)
                LocMean(i, :)  = Mean(:)
                LocStDev(i, :) = StDev(:)
            end where
        end do
        !> special case first window
        if (wdw == 1) then
            do i = 1, imin
                where(E2Col(:)%present)
                    LocMean(i, :)  = Mean(:)
                    LocStDev(i, :) = StDev(:)
                end where
            end do
        end if
        !> special case last window
        if (wdw == wdw_num) then
            do i = imax, N
                where(E2Col(:)%present)
                    LocMean(i, :)  = Mean(:)
                    LocStDev(i, :) = StDev(:)
                end where
            end do
        end if
    end do

    !> spikes detection and removal (if requested) in the whole file
    do j = u, pe
        if (E2Col(j)%present) then
            cnt = 0
            !> special case first record in the file
            if (Set(1, j) /= error .and. Set(1, j) > &
                LocMean(1, j) + adv_lim(j) * LocStDev(1, j)) then
                Set(1, j) = LocMean(1, j) + adv_lim(j) * LocStDev(1, j)
                    if (.not. IsSpike(1, j)) then
                        nspikes(j) = nspikes(j) + 1
                        nspikes_sng(j) = nspikes_sng(j) + 1
                        IsSpike(1, j) = .true.
                    end if
            end if
            if (Set(1, j) /= error .and. Set(1, j) < &
                LocMean(1, j) - adv_lim(j) * LocStDev(1, j)) then
                Set(1, j) = LocMean(1, j) - adv_lim(j) * LocStDev(1, j)
                    if (.not. IsSpike(1, j)) then
                        nspikes(j) = nspikes(j) + 1
                        nspikes_sng(j) = nspikes_sng(j) + 1
                        IsSpike(1, j) = .true.
                    end if
            end if
            !> Following lines
            !>
            !> `cnt` counts spiky SAMPLES. It used to count rows, which is the
            !> same number only for a column sampled at the file's rate: a
            !> slower column carries the error code on the rows between its
            !> samples, and those rows extended the run
            !> (`if (cnt /= 0) cnt = cnt + 1`) until it passed sr%num_spk and
            !> was discarded as a physical excursion. A single-sample spike in a
            !> 1 Hz analyser on a 10 Hz file counted as ten, so with the default
            !> allowance of three the column was never despiked at all - and it
            !> reported no spikes rather than "not assessed".
            !>
            !> An error row now neither counts towards the run nor breaks it:
            !> consecutive samples are what the setting is about. The run's rows
            !> are tracked alongside, because the flagging and the interpolation
            !> address rows.
            prev_good = 0
            if (Set(1, j) /= error) prev_good = 1
            run_start = 0
            run_last = 0
            if (RPsetup%filter_sr) then
                do i = 2, N
                    if (Set(i, j) == error) cycle
                    if (Set(i, j) > &
                        LocMean(i, j) + adv_lim(j) * LocStDev(i, j) .or. &
                        Set(i, j) < &
                        LocMean(i, j) - adv_lim(j) * LocStDev(i, j)) then
                        if (cnt == 0) run_start = i
                        run_last = i
                        cnt = cnt + 1
                    else
                        if ((cnt /= 0) .and. (cnt <= sr%num_spk)) then
                            !> check whether it was a spike already, if
                            !> not increment the number of spikes found
                            new_spike = .not. IsSpike(run_last, j)
                            if (new_spike) then
                                nspikes(j) = nspikes(j) + 1
                                nspikes_sng(j) = nspikes_sng(j) + cnt
                            end if
                            !> update spike flags, on the samples only
                            do k = run_start, i - 1
                                if (Set(k, j) /= error) IsSpike(k, j) = .true.
                            enddo

                            !> replace with linear interpolation if requested
                            !> and if possible, i.e. if there is a good sample
                            !> behind the run to draw the line from.
                            !>
                            !> The line is drawn in rows, not in samples, so it
                            !> is the same line as before for a column at the
                            !> file rate and the right one for a slower column.
                            !> It writes the sample rows only: filling the error
                            !> rows in between would manufacture samples the
                            !> instrument never reported, and those rows are
                            !> what tells the rest of the run how much of the
                            !> column is real.
                            if (prev_good /= 0) then
                                m = (Set(i, j) - Set(prev_good, j)) &
                                    / dble(i - prev_good)
                                q = Set(prev_good, j)
                                do k = run_start, i - 1
                                    if (Set(k, j) /= error) &
                                        Set(k, j) = m * dble(k - prev_good) + q
                                end do
                            end if
                            cnt = 0
                        else if (cnt > sr%num_spk) then
                            cnt = 0
                        end if
                        prev_good = i
                    end if
                end do
            else
                !> Detection only, no removal - but counted the same way, so
                !> that a flag-only run reports the same spikes a despiking run
                !> would remove. See the comment above on samples versus rows.
                do i = 2, N
                    if (Set(i, j) == error) cycle
                    if (Set(i, j) > &
                            LocMean(i, j) + adv_lim(j) * LocStDev(i, j) .or. &
                            Set(i, j) < &
                            LocMean(i, j) - adv_lim(j) * LocStDev(i, j)) then
                            if (cnt == 0) run_start = i
                            run_last = i
                            cnt = cnt + 1
                    else
                        if ((cnt /= 0) .and. (cnt <= sr%num_spk)) then
                            !> check whether it was a spike already,
                            !> if not increment the number of spikes found
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
                        else if (cnt > sr%num_spk) then
                            cnt = 0
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
        adv_lim(:) = adv_lim(:) + lim_step
        goto 100
    end if
    deallocate(XX)

    !> hflags the variable if nspikes is larger than a prescribed threshold.
    !>
    !> Every gas slot, not just the historical four. A slot that is absent
    !> keeps the 9 this is initialised to, so a project with four gases packs
    !> exactly the string it did before; but bounding the loop at gas4 meant a
    !> fifth gas's spike outcome was never computed, and the FLUXNET and full
    !> outputs reported it as "test not performed" however it had gone.
    !>
    !> The percentage is of the column's own samples, not of the file's rows.
    !> Dividing by N diluted a slower column's spike count by its stride, so a
    !> 1 Hz analyser on a 10 Hz file needed ten times the spikes to reach
    !> sr%hf_lim as the sonic beside it did.
    hflags = 9
    do j = u, lastGas
        if (E2Col(j)%present) then
            nsample = nint(dble(N) * ColumnAcFreq(j) / Metadata%ac_freq)
            if (nsample <= 0) nsample = N
            if(100.d0 * (dble(tot_spikes(j)) / dble(nsample)) >= sr%hf_lim) then
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

end subroutine TestSpikeDetectionVickers97
