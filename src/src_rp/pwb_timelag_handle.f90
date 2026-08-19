!***************************************************************************
! pwb_timelag_handle.f90
! ----------------------
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
! \brief       Native pre-whitening block-bootstrap time lag detector.
! \author      Jonathan Muller
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
module m_pwb_timelag
    use m_rp_global_var
    implicit none
    private
    public :: PwbDetectGas, ResetPwbDiagnostics, ReportPwbDiagnostics, InitPwbResult, GasLabel
    public :: CountPwbDiagnostic, SameAnalyser
    public :: InitPwbTimelagCache, ReadPwbTimelagCache, WritePwbTimelagCache
    public :: LookupPwbTimelagCache, StorePwbTimelagCache, SetPwbPeriodTimestamp
    public :: PostProcessPwbTimelagCache
    public :: ResetPwbAggregateSummary, AddPwbTimelagSummaryDataset, ResolvePwbAggregateSummary

    integer :: pwb_attempts(E2NumVar) = 0
    integer :: pwb_successes(E2NumVar) = 0
    integer :: pwb_carryforwards(E2NumVar) = 0
    integer :: pwb_fallbacks(E2NumVar) = 0
    integer :: pwb_fallback_maxcov(E2NumVar) = 0
    integer :: pwb_fallback_nominal(E2NumVar) = 0
    integer :: pwb_fallback_other(E2NumVar) = 0
    integer :: pwb_instrument_shared(E2NumVar) = 0
    integer :: pwb_outside_window(E2NumVar) = 0
    logical :: pwb_bounds_warned(E2NumVar) = .false.
    logical :: pwb_block_warned(E2NumVar) = .false.

    !> Scratch shared by the four combinations of one gas, reused across gases
    !> and periods. These were allocated and freed inside the bootstrap loop -
    !> n_bootstrap * 4 allocation pairs per gas per averaging period, for
    !> arrays whose size never changes within a period. The engine has no
    !> threading, so one set of module-level buffers is enough.
    real(kind = dbl), allocatable :: sc_xb(:), sc_yb(:), sc_xc(:), sc_yc(:)
    real(kind = dbl), allocatable :: sc_ccf(:), sc_smooth(:)
    real(kind = dbl), allocatable :: sc_mean_ccf(:), sc_mean_smooth(:)
    real(kind = dbl), allocatable :: sc_hdi(:)
    integer, allocatable :: sc_boot(:)
    integer :: sc_n = 0, sc_lo = 0, sc_hi = 0, sc_nboot = 0

contains

subroutine ResetPwbDiagnostics()
    pwb_attempts = 0
    pwb_successes = 0
    pwb_carryforwards = 0
    pwb_fallbacks = 0
    pwb_fallback_maxcov = 0
    pwb_fallback_nominal = 0
    pwb_fallback_other = 0
    pwb_instrument_shared = 0
    pwb_outside_window = 0
    pwb_bounds_warned = .false.
    pwb_block_warned = .false.
end subroutine ResetPwbDiagnostics

subroutine ResetPwbAggregateSummary()
    PwbSummaryDonorCount = 0
    PwbSummarySource = 0
    PwbSummaryEvidence = 0
end subroutine ResetPwbAggregateSummary

subroutine AddPwbTimelagSummaryDataset(TimelagOpt, nrow, n)
    integer, intent(in) :: nrow, n
    type(TimeLagOptType), intent(inout) :: TimelagOpt(nrow)
    integer :: gas, origin, wsl
    integer, external :: PrimaryWaterOutSlot

    TimelagOpt(n)%tlag = error
    TimelagOpt(n)%RH = error
    do gas = firstGas, lastGas
        if (.not. E2Col(gas)%present) cycle
        if (trim(PWBResult(gas)%reliability_class) == 'S1_optimal' .or. &
            trim(PWBResult(gas)%reliability_class) == 'S2_optimal') then
            TimelagOpt(n)%tlag(gas) = Essentials%used_timelag(gas)
        else
            origin = PWBResult(gas)%origin_gas
            if (origin >= firstGas .and. origin <= lastGas .and. origin /= gas) &
                PwbSummaryDonorCount(gas, origin) = PwbSummaryDonorCount(gas, origin) + 1
        end if
    end do
    !> RH travels with the water record's own time-lag, so it is gated on the
    !> site's water rather than on slot six.
    wsl = PrimaryWaterOutSlot()
    if (E2Col(wsl)%present .and. TimelagOpt(n)%tlag(wsl) /= error &
        .and. Stats%RH >= 0d0 .and. Stats%RH <= 100d0) then
        TimelagOpt(n)%RH = Stats%RH
    else
        TimelagOpt(n)%tlag(wsl) = error
    end if
end subroutine AddPwbTimelagSummaryDataset

subroutine ResolvePwbAggregateSummary(actn)
    integer, intent(inout) :: actn(E2NumVar)
    integer :: gas, donor, best_count, candidate
    logical, external :: GasSlotIsWater

    PwbSummarySource = 0
    PwbSummaryEvidence = 0
    do gas = firstGas, lastGas
        if (E2Col(gas)%present .and. actn(gas) > 0) PwbSummarySource(gas) = gas
    end do
    do gas = firstGas, lastGas
        if (.not. E2Col(gas)%present .or. PwbSummarySource(gas) /= 0) cycle
        donor = 0
        best_count = 0
        do candidate = firstGas, lastGas
            !> Water is never a donor: its lag is RH-dependent in a way the
            !> trace gases' are not, so borrowing from it is worse than not
            !> borrowing. That is a property of the species, and it was
            !> written as `candidate == h2o` - the historical sixth slot. On a
            !> project whose water sits elsewhere that excluded whichever gas
            !> occupied slot six from ever donating, and let the real
            !> hygrometer donate to everything.
            if (candidate == gas) cycle
            if (GasSlotIsWater(candidate)) cycle
            if (PwbSummarySource(candidate) == 0) cycle
            if (PwbSummaryDonorCount(gas, candidate) > best_count) then
                donor = candidate
                best_count = PwbSummaryDonorCount(gas, candidate)
            end if
        end do
        if (donor > 0) then
            toPasGas(gas) = toPasGas(donor)
            actn(gas) = actn(donor)
            PwbSummarySource(gas) = donor
            PwbSummaryEvidence(gas) = best_count
        end if
    end do
end subroutine ResolvePwbAggregateSummary

subroutine InitPwbTimelagCache()
    if (allocated(PwbTimelagCache)) deallocate(PwbTimelagCache)
    PwbTimelagCacheN = 0
    PwbCacheLoaded = .false.
    PwbCacheDirty = .false.
end subroutine InitPwbTimelagCache

subroutine SetPwbPeriodTimestamp(date, time)
    character(*), intent(in) :: date, time

    PwbPeriodDate = date
    PwbPeriodTime = time
end subroutine SetPwbPeriodTimestamp

logical function ValidPwbPeriodTimestamp(date, time)
    character(*), intent(in) :: date, time

    ValidPwbPeriodTimestamp = len_trim(date) > 0 .and. len_trim(time) > 0 &
        .and. index(date, achar(0)) == 0 .and. index(time, achar(0)) == 0
end function ValidPwbPeriodTimestamp

!***************************************************************************
!> The exact inverse of GasLabel, over every slot.
!>
!> No literal cases: one would shadow the record-derived match and send a
!> label back to a slot that no longer holds that species.
!***************************************************************************
integer function GasIndexFromLabel(label)
    character(*), intent(in) :: label
    character(64) :: tags(GHGNumVar)
    integer :: gas

    GasIndexFromLabel = 0
    call SpectralVarTags(tags)
    do gas = firstGas, lastGas
        if (len_trim(tags(gas)) > 0 .and. trim(tags(gas)) == trim(label)) then
            GasIndexFromLabel = gas
            return
        end if
    end do
end function GasIndexFromLabel

subroutine StorePwbTimelagCache(gas, actual_lag, used_lag, row_lag, default_used, res)
    integer, intent(in) :: gas, row_lag
    real(kind = dbl), intent(in) :: actual_lag, used_lag
    logical, intent(in) :: default_used
    type(PWBResultType), intent(in) :: res
    integer :: i

    if (.not. ValidPwbPeriodTimestamp(PwbPeriodDate, PwbPeriodTime)) then
        call LogSay(' Fatal> PWB time-lag file cannot use an empty or invalid period timestamp.')
        error stop 'Invalid PWB period timestamp.'
    end if
    do i = 1, PwbTimelagCacheN
        if (PwbTimelagCache(i)%date == PwbPeriodDate .and. PwbTimelagCache(i)%time == PwbPeriodTime &
            .and. PwbTimelagCache(i)%gas == gas) then
            PwbTimelagCache(i)%actual_lag = actual_lag
            PwbTimelagCache(i)%used_lag = used_lag
            PwbTimelagCache(i)%row_lag = row_lag
            PwbTimelagCache(i)%default_used = default_used
            PwbTimelagCache(i)%result = res
            PwbCacheDirty = .true.
            return
        end if
    end do

    call StorePwbTimelagCacheAt(PwbPeriodDate, PwbPeriodTime, gas, actual_lag, &
        used_lag, row_lag, default_used, res)
    PwbCacheDirty = .true.
end subroutine StorePwbTimelagCache

subroutine LookupPwbTimelagCache(gas, found, actual_lag, used_lag, row_lag, default_used, res)
    integer, intent(in) :: gas
    logical, intent(out) :: found, default_used
    real(kind = dbl), intent(out) :: actual_lag, used_lag
    integer, intent(out) :: row_lag
    type(PWBResultType), intent(out) :: res
    integer :: i

    found = .false.
    actual_lag = error
    used_lag = error
    row_lag = 0
    default_used = .false.
    call InitPwbResult(res)
    if (.not. PwbCacheLoaded) return
    if (.not. ValidPwbPeriodTimestamp(PwbPeriodDate, PwbPeriodTime)) then
        call LogSay(' Fatal> PWB time-lag file cannot use an empty or invalid period timestamp.')
        error stop 'Invalid PWB period timestamp.'
    end if
    do i = 1, PwbTimelagCacheN
        if (PwbTimelagCache(i)%date == PwbPeriodDate .and. PwbTimelagCache(i)%time == PwbPeriodTime &
            .and. PwbTimelagCache(i)%gas == gas) then
            found = .true.
            actual_lag = PwbTimelagCache(i)%actual_lag
            used_lag = PwbTimelagCache(i)%used_lag
            row_lag = PwbTimelagCache(i)%row_lag
            default_used = PwbTimelagCache(i)%default_used
            res = PwbTimelagCache(i)%result
            return
        end if
    end do
end subroutine LookupPwbTimelagCache

!***************************************************************************
!> Fingerprint the settings a cached time-lag depends on.
!>
!> Anything not in here is a setting a user can change without the file
!> noticing, so the next run reuses a lag computed under the old value.
!>
!> Wide enough for every gas the engine can hold. At character(256) a project
!> past roughly the eighth gas ran out of room, and Fortran truncates a
!> character assignment silently - so two different lag-window configurations
!> could produce the same string and a stale entry would be accepted as
!> current.
!***************************************************************************
character(2048) function PwbCacheFingerprint()
    integer :: gas
    character(64) :: extra

    !> One loop over the records, rather than a fixed block for the four
    !> legacy gases with the rest appended: the fingerprint's shape used to
    !> depend on how many gases a project had relative to four.
    write(PwbCacheFingerprint, &
        '(a,i0,a,f8.4,a,f8.4,a,f8.4,a,f8.4,a,f8.4,a,i0,a,i0)') &
        'n=', PWBSetup%n_bootstrap, '_block=', PWBSetup%block_length_s, &
        '_valid=', PWBSetup%min_valid_frac, &
        '_hdi=', PWBSetup%hdi_thresh_s, '_dev=', PWBSetup%dev_thresh_s, &
        '_prefilter=', PWBSetup%hdi_prefilter_s, &
        '_smooth=', PWBSetup%smoothing_width, '_seed=', PWBSetup%random_seed

    do gas = firstGas, lastGas
        if (gas - firstGas + 1 > min(EddyFlowProj%gas_num, MaxNumGases)) exit
        if (.not. PWBSetup%lag_bounds_provided(gas)) cycle
        write(extra, '(a,i0,a,f10.4,a,f10.4)') ':g', gas - firstGas + 1, &
            '=', PWBSetup%min_lag(gas), ':', PWBSetup%max_lag(gas)
        PwbCacheFingerprint = &
            trim(PwbCacheFingerprint) // trim(adjustl(extra))
    end do
end function PwbCacheFingerprint

subroutine ReadPwbTimelagCache(path, recognized, valid)
    character(*), intent(in) :: path
    logical, intent(out) :: recognized, valid
    integer :: u, ios, gas, row_lag, period_seconds, origin_gas, raw_row_lag
    character(1024) :: line
    character(2048) :: fingerprint
    character(10) :: date
    character(5) :: time
    character(24) :: reliability, fallback, fill_method
    character(32) :: donor, gas_label, origin_label
    character(2) :: combo
    real(kind = dbl) :: actual_lag, used_lag, selected_lag, hdi_low, hdi_high, hdi_range
    real(kind = dbl) :: unrestricted_lag, raw_cov, eff_min, eff_max, eff_block
    logical :: default_used, edge_pinned, outside, clamped, prefiltered
    type(PWBResultType) :: res

    recognized = .false.
    valid = .false.
    call InitPwbTimelagCache()
    open(newunit=u, file=path, status='old', action='read', iostat=ios, encoding='utf-8')
    if (ios /= 0) return
    read(u, '(a)', iostat=ios) line
    if (ios /= 0) then
        close(u)
        return
    end if
    !> Versions 1 and 2 are not read. Both carried a pre_wpl/post_wpl stage
    !> column for a choice that no longer exists, and both predate the two
    !> retired speed settings - so their fingerprint could not match this
    !> build in any case. Refusing them outright says so plainly, instead of
    !> letting the fingerprint say it obscurely.
    if (trim(line) /= 'PWB_TIMELAG_CACHE_VERSION=3') then
        close(u)
        return
    end if
    recognized = .true.
    read(u, '(a)', iostat=ios) line
    if (ios /= 0 .or. index(line, 'fingerprint=') /= 1) then
        close(u)
        return
    end if
    fingerprint = line(13:len_trim(line))
    if (trim(fingerprint) /= trim(PwbCacheFingerprint())) then
        call LogSay(' Fatal> PWB time-lag file was written under different PWB settings.')
        call LogSay('        Regenerate it, or restore the settings it was written with.')
        close(u)
        return
    end if
    read(u, '(a)', iostat=ios) line
    if (ios /= 0 .or. index(line, 'project_id=') /= 1) then
        close(u)
        return
    end if
    read(u, '(a)', iostat=ios) line
    if (ios /= 0 .or. index(line, 'period_seconds=') /= 1) then
        close(u)
        return
    end if
    read(line(16:len_trim(line)), *, iostat=ios) period_seconds
    if (ios /= 0 .or. period_seconds /= RPsetup%avrg_len) then
        call LogSay(' Fatal> PWB time-lag file averaging-period duration is incompatible with this project.')
        close(u)
        return
    end if
    read(u, '(a)', iostat=ios) line
    if (ios /= 0 .or. trim(line) /= 'data') then
        close(u)
        return
    end if
    read(u, '(a)', iostat=ios) line
    if (ios /= 0 .or. index(line, 'date,time,gas,') /= 1) then
        close(u)
        return
    end if
    do
        read(u, '(a)', iostat=ios) line
        if (ios /= 0) exit
        if (len_trim(line) == 0) cycle
        call InitPwbResult(res)
        read(line, *, iostat=ios) date, time, gas_label, &
            selected_lag, raw_row_lag, hdi_low, hdi_high, hdi_range, combo, edge_pinned, &
            unrestricted_lag, outside, raw_cov, &
            actual_lag, used_lag, row_lag, default_used, &
            reliability, fill_method, donor, origin_label, fallback, &
            eff_min, eff_max, eff_block, clamped, prefiltered
        if (ios /= 0 .or. .not. ValidPwbPeriodTimestamp(date, time)) then
            close(u)
            return
        end if
        !> The gas is named, so a file whose records were reordered between
        !> two runs still restores each lag onto the species it was measured
        !> for rather than onto whatever now occupies that slot.
        gas = GasIndexFromLabel(gas_label)
        if (gas == 0) then
            call LogSay(' Fatal> PWB time-lag file names a gas this project does not measure: ' &
                // trim(gas_label))
            close(u)
            return
        end if
        origin_gas = GasIndexFromLabel(origin_label)
        res%selected_lag = selected_lag
        res%row_lag = raw_row_lag
        res%hdi_low = hdi_low
        res%hdi_high = hdi_high
        res%hdi_range = hdi_range
        res%best_combination = combo
        res%edge_pinned = edge_pinned
        res%unrestricted_peak_lag = unrestricted_lag
        res%peak_outside_window = outside
        res%raw_covariance = raw_cov
        res%reliability_class = reliability
        res%fill_method = fill_method
        res%donor_gas = donor
        res%origin_gas = origin_gas
        res%fallback_source = fallback
        res%effective_min_lag = eff_min
        res%effective_max_lag = eff_max
        res%effective_block_length_s = eff_block
        res%block_length_clamped = clamped
        res%hdi_prefiltered = prefiltered
        res%applied_lag = used_lag
        res%applied_row_lag = row_lag
        res%fallback_used = default_used .or. trim(reliability) == 'fallback'
        call StorePwbTimelagCacheAt(date, time, gas, actual_lag, used_lag, &
            row_lag, default_used, res)
    end do
    close(u)
    PwbCacheLoaded = .true.
    PwbCacheDirty = .false.
    valid = PwbTimelagCacheN > 0
end subroutine ReadPwbTimelagCache

subroutine StorePwbTimelagCacheAt(date, time, gas, actual_lag, used_lag, row_lag, default_used, res)
    character(*), intent(in) :: date, time
    integer, intent(in) :: gas, row_lag
    real(kind = dbl), intent(in) :: actual_lag, used_lag
    logical, intent(in) :: default_used
    type(PWBResultType), intent(in) :: res
    type(PWBTimelagCacheEntryType), allocatable :: tmp(:)

    allocate(tmp(PwbTimelagCacheN + 1))
    if (PwbTimelagCacheN > 0) tmp(1:PwbTimelagCacheN) = PwbTimelagCache(1:PwbTimelagCacheN)
    tmp(PwbTimelagCacheN + 1)%date = date
    tmp(PwbTimelagCacheN + 1)%time = time
    tmp(PwbTimelagCacheN + 1)%gas = gas
    tmp(PwbTimelagCacheN + 1)%actual_lag = actual_lag
    tmp(PwbTimelagCacheN + 1)%used_lag = used_lag
    tmp(PwbTimelagCacheN + 1)%row_lag = row_lag
    tmp(PwbTimelagCacheN + 1)%default_used = default_used
    tmp(PwbTimelagCacheN + 1)%result = res
    call move_alloc(tmp, PwbTimelagCache)
    PwbTimelagCacheN = PwbTimelagCacheN + 1
end subroutine StorePwbTimelagCacheAt

subroutine WritePwbTimelagCache()
    integer :: u, ios, i
    character(PathLen) :: path

    if (.not. PwbCacheDirty .or. PwbTimelagCacheN == 0 .or. Dir%main_out == 'none') return
    path = Dir%main_out(1:len_trim(Dir%main_out)) // EddyFlowProj%id(1:len_trim(EddyFlowProj%id)) &
        // PwbTimelag_FilePadding // Timestamp_FilePadding // CsvExt
    open(newunit=u, file=path, status='replace', iostat=ios, encoding='utf-8')
    if (ios /= 0) return
    write(u, '(a)') 'PWB_TIMELAG_CACHE_VERSION=3'
    write(u, '(a)') 'fingerprint=' // trim(PwbCacheFingerprint())
    write(u, '(a)') 'project_id=' // trim(EddyFlowProj%id)
    write(u, '(a,i0)') 'period_seconds=', RPsetup%avrg_len
    write(u, '(a)') 'data'
    write(u, '(a)') 'date,time,gas,' &
        // 'raw_lag_s,raw_row_lag,hdi_low_s,hdi_high_s,hdi_range_s,best_combination,edge_pinned,' &
        // 'unrestricted_peak_lag_s,peak_outside_window,raw_covariance,' &
        // 'actual_lag_s,used_lag_s,used_row_lag,default_used,' &
        // 'reliability_class,fill_method,donor_gas,origin_gas,fallback_source,' &
        // 'effective_min_lag_s,effective_max_lag_s,effective_block_length_s,block_length_clamped,' &
        // 'hdi_prefiltered'
    do i = 1, PwbTimelagCacheN
        write(u, '(a,",",a,",",a,",",f12.6,",",i0,",",f12.6,",",f12.6,",",f12.6,",",a,",",l1,' &
            // '",",f12.6,",",l1,",",f14.6,",",f12.6,",",f12.6,",",i0,",",l1,' &
            // '",",a,",",a,",",a,",",a,",",a,' &
            // '",",f12.6,",",f12.6,",",f12.6,",",l1,",",l1)') &
            trim(PwbTimelagCache(i)%date), trim(PwbTimelagCache(i)%time), &
            trim(GasLabel(PwbTimelagCache(i)%gas)), &
            PwbTimelagCache(i)%result%selected_lag, PwbTimelagCache(i)%result%row_lag, &
            PwbTimelagCache(i)%result%hdi_low, PwbTimelagCache(i)%result%hdi_high, &
            PwbTimelagCache(i)%result%hdi_range, &
            trim(PwbTimelagCache(i)%result%best_combination), PwbTimelagCache(i)%result%edge_pinned, &
            PwbTimelagCache(i)%result%unrestricted_peak_lag, &
            PwbTimelagCache(i)%result%peak_outside_window, &
            PwbTimelagCache(i)%result%raw_covariance, &
            PwbTimelagCache(i)%actual_lag, PwbTimelagCache(i)%used_lag, &
            PwbTimelagCache(i)%row_lag, PwbTimelagCache(i)%default_used, &
            trim(PwbTimelagCache(i)%result%reliability_class), &
            trim(PwbTimelagCache(i)%result%fill_method), &
            trim(PwbTimelagCache(i)%result%donor_gas), &
            trim(GasLabel(PwbTimelagCache(i)%result%origin_gas)), &
            trim(PwbTimelagCache(i)%result%fallback_source), &
            PwbTimelagCache(i)%result%effective_min_lag, PwbTimelagCache(i)%result%effective_max_lag, &
            PwbTimelagCache(i)%result%effective_block_length_s, &
            PwbTimelagCache(i)%result%block_length_clamped, &
            PwbTimelagCache(i)%result%hdi_prefiltered
    end do
    close(u)
    PwbTimelagCache_Path = path
    PwbCacheLoaded = .true.
    PwbCacheDirty = .false.
    write(*, '(a)') ' PWB half-hourly time-lag table written to: ' // trim(path)
    write(ulog, '(a)') ' PWB half-hourly time-lag table written to: ' // trim(path)
end subroutine WritePwbTimelagCache

!***************************************************************************
!> Decide every period's time-lag with the whole record in hand.
!>
!> The classifier in timelag_handle runs while periods stream past, so it can
!> only ever look backwards: a period before the first reliable detection has
!> nothing to carry forward and falls back to covariance maximisation, and a
!> gap in the middle carries the last good lag forward for as long as the gap
!> lasts, however long that is. Neither is a limitation of the method - it is
!> a limitation of not having read the rest of the run yet.
!>
!> The pre-generation pass has read the rest of the run. This walks the
!> finished table per gas in time order and settles each period from the raw
!> evidence already stored, in the order the reference post-processing uses
!> (apply_hdi_prefilter, then apply_pwbopt, then fill_tlag_gaps):
!>
!>   prefilter -> S1 -> S2 -> same-analyser share -> interpolate -> back-fill
!>   -> carry-forward past the last detection -> per-gas median -> whatever
!>   the pass itself settled on
!>
!> fill_method records which arm settled each row, so the file says how every
!> half-hour was decided and not only what it decided.
!***************************************************************************
subroutine PostProcessPwbTimelagCache()
    integer :: i, j, k, g, gas, n, nsel
    integer, allocatable :: ord(:), idx(:)
    integer(8), allocatable :: tmin(:)
    real(kind = dbl), allocatable :: lag(:), fallback_lag(:), sorted(:)
    logical, allocatable :: settled(:)
    real(kind = dbl) :: median_lag, t0, t1, span, previous
    integer :: prev, nxt, shared
    logical, external :: GasSlotIsWater
    logical :: any_shared

    if (PwbTimelagCacheN <= 0) return

    !> Chronological order over the whole table. The pass appends in period
    !> order already, but the interpolation below reads a time axis off these
    !> rows and must not depend on that being true.
    allocate(ord(PwbTimelagCacheN), tmin(PwbTimelagCacheN))
    do i = 1, PwbTimelagCacheN
        ord(i) = i
        tmin(i) = PeriodMinutes(PwbTimelagCache(i)%date, PwbTimelagCache(i)%time)
    end do
    call SortByMinutes(ord, tmin, PwbTimelagCacheN)

    !> Step 1, over every row: the HDI pre-filter. A detection wider than the
    !> pre-filter is discarded before classification, so temporal continuity
    !> (S2 below) cannot rescue it just because it happens to land near the
    !> previous lag. Zero disables it, which is what the interface writes at
    !> its "Disabled" spin position.
    !>
    !> This setting has existed in the interface, the project file and the
    !> fingerprint since PWB was added, and until now nothing read it.
    do i = 1, PwbTimelagCacheN
        PwbTimelagCache(i)%result%hdi_prefiltered = .false.
        if (PWBSetup%hdi_prefilter_s > 0d0 &
            .and. PwbTimelagCache(i)%result%hdi_range /= error &
            .and. PwbTimelagCache(i)%result%hdi_range > PWBSetup%hdi_prefilter_s) &
            PwbTimelagCache(i)%result%hdi_prefiltered = .true.
    end do

    allocate(idx(PwbTimelagCacheN), lag(PwbTimelagCacheN), &
        fallback_lag(PwbTimelagCacheN), settled(PwbTimelagCacheN), sorted(PwbTimelagCacheN))

    !> Steps 2 and 3: S1 accepts a narrow HDI outright; S2 accepts a wider one
    !> that stays close to the last accepted lag.
    do gas = firstGas, lastGas
        n = 0
        do k = 1, PwbTimelagCacheN
            i = ord(k)
            if (PwbTimelagCache(i)%gas /= gas) cycle
            n = n + 1
            idx(n) = i
            lag(n) = error
            settled(n) = .false.
        end do
        if (n == 0) cycle

        previous = error
        do j = 1, n
            i = idx(j)
            if (PwbTimelagCache(i)%result%hdi_prefiltered) cycle
            if (PwbTimelagCache(i)%result%selected_lag == error) cycle
            if (PwbTimelagCache(i)%result%edge_pinned) cycle
            if (PwbTimelagCache(i)%result%hdi_range == error) cycle
            if (PwbTimelagCache(i)%result%hdi_range < PWBSetup%hdi_thresh_s) then
                lag(j) = PwbTimelagCache(i)%result%selected_lag
                settled(j) = .true.
            elseif (previous /= error) then
                if (abs(PwbTimelagCache(i)%result%selected_lag - previous) &
                    <= PWBSetup%dev_thresh_s) then
                    lag(j) = PwbTimelagCache(i)%result%selected_lag
                    settled(j) = .true.
                end if
            end if
            if (settled(j)) then
                previous = lag(j)
                if (PwbTimelagCache(i)%result%hdi_range < PWBSetup%hdi_thresh_s) then
                    PwbTimelagCache(i)%result%reliability_class = 'S1_optimal'
                else
                    PwbTimelagCache(i)%result%reliability_class = 'S2_optimal'
                end if
                PwbTimelagCache(i)%used_lag = lag(j)
                PwbTimelagCache(i)%result%fill_method = 'native'
                PwbTimelagCache(i)%result%fallback_source = 'native'
                PwbTimelagCache(i)%result%fallback_used = .false.
                PwbTimelagCache(i)%result%origin_gas = gas
                PwbTimelagCache(i)%result%donor_gas = 'none'
            end if
        end do
        do j = 1, n
            if (.not. settled(j)) PwbTimelagCache(idx(j))%result%reliability_class = 'pending'
        end do
    end do

    !> Step 4: a gas with no detection of its own in a period takes the lag of
    !> another gas measured by the same analyser in that same period. They
    !> share a tube, so they share a delay.
    !>
    !> Two things this asks that the streaming version did not. It matches on
    !> instrument identity rather than on the model string, because two
    !> analysers of the same model are two tubes and sharing a lag between
    !> them is wrong - the same distinction timelag_handle already draws for
    !> the water covariance. And it never accepts water as the donor: water's
    !> lag is RH-dependent in a way the trace gases' is not, which is exactly
    !> why ResolvePwbAggregateSummary above refuses it. The per-period rule
    !> allowed it, so the two halves of one rule disagreed.
    !>
    !> The outer loop repeats because a gas that has just adopted a lag does
    !> not itself become a donor - only S1 and S2 rows donate - so one pass is
    !> enough in fact; the loop simply makes that explicit and terminates on
    !> the first pass that shares nothing.
    any_shared = .true.
    do while (any_shared)
        any_shared = .false.
        do i = 1, PwbTimelagCacheN
            gas = PwbTimelagCache(i)%gas
            if (gas < firstGas .or. gas > lastGas) cycle
            if (trim(PwbTimelagCache(i)%result%reliability_class) /= 'pending') cycle
            shared = 0
            do j = 1, PwbTimelagCacheN
                if (j == i) cycle
                if (PwbTimelagCache(j)%date /= PwbTimelagCache(i)%date) cycle
                if (PwbTimelagCache(j)%time /= PwbTimelagCache(i)%time) cycle
                g = PwbTimelagCache(j)%gas
                if (g < firstGas .or. g > lastGas) cycle
                if (GasSlotIsWater(g)) cycle
                if (.not. SameAnalyser(gas, g)) cycle
                if (trim(PwbTimelagCache(j)%result%reliability_class) /= 'S1_optimal' &
                    .and. trim(PwbTimelagCache(j)%result%reliability_class) /= 'S2_optimal') cycle
                shared = j
                exit
            end do
            if (shared == 0) cycle
            PwbTimelagCache(i)%used_lag = PwbTimelagCache(shared)%used_lag
            PwbTimelagCache(i)%result%reliability_class = 'S4_instrument_shared'
            PwbTimelagCache(i)%result%fill_method = 'instrument_shared'
            PwbTimelagCache(i)%result%fallback_source = 'instrument_shared'
            PwbTimelagCache(i)%result%fallback_used = .false.
            PwbTimelagCache(i)%result%donor_gas = GasLabel(PwbTimelagCache(shared)%gas)
            PwbTimelagCache(i)%result%origin_gas = PwbTimelagCache(shared)%gas
            any_shared = .true.
        end do
    end do

    !> Steps 5 to 8, per gas: fill whatever is still open.
    do gas = firstGas, lastGas
        n = 0
        do k = 1, PwbTimelagCacheN
            i = ord(k)
            if (PwbTimelagCache(i)%gas /= gas) cycle
            n = n + 1
            idx(n) = i
            !> Whatever the streaming pass settled on stands as this period's
            !> last resort - covariance maximisation where it fell back, the
            !> detection itself otherwise. Read before it is overwritten.
            fallback_lag(n) = PwbTimelagCache(i)%used_lag
            settled(n) = trim(PwbTimelagCache(i)%result%reliability_class) == 'S1_optimal' &
                .or. trim(PwbTimelagCache(i)%result%reliability_class) == 'S2_optimal' &
                .or. trim(PwbTimelagCache(i)%result%reliability_class) == 'S4_instrument_shared'
            if (settled(n)) then
                lag(n) = PwbTimelagCache(i)%used_lag
            else
                lag(n) = error
            end if
        end do
        if (n == 0) cycle

        !> Step 5: an interior gap is interpolated between the reliable lags
        !> that bracket it, on the real time axis rather than on row number -
        !> a run with missing files has periods that are not equally spaced.
        !> Carrying a lag forward holds it flat across a drifting pump;
        !> interpolating tracks the drift.
        do j = 1, n
            if (settled(j)) cycle
            prev = 0
            do k = j - 1, 1, -1
                if (settled(k)) then
                    prev = k
                    exit
                end if
            end do
            nxt = 0
            do k = j + 1, n
                if (settled(k)) then
                    nxt = k
                    exit
                end if
            end do
            if (prev == 0 .or. nxt == 0) cycle
            i = idx(j)
            t0 = dble(tmin(idx(prev)))
            t1 = dble(tmin(idx(nxt)))
            span = t1 - t0
            if (span <= 0d0) then
                PwbTimelagCache(i)%used_lag = lag(prev)
            else
                PwbTimelagCache(i)%used_lag = lag(prev) &
                    + (lag(nxt) - lag(prev)) * (dble(tmin(i)) - t0) / span
            end if
            PwbTimelagCache(i)%result%reliability_class = 'S3_interpolated'
            PwbTimelagCache(i)%result%fill_method = 'interpolated'
            PwbTimelagCache(i)%result%fallback_source = 'interpolated'
            PwbTimelagCache(i)%result%origin_gas = gas
        end do

        nxt = 0
        do j = 1, n
            if (settled(j)) then
                nxt = j
                exit
            end if
        end do
        if (nxt > 0) then
            !> Step 6: periods before the first reliable detection have
            !> nothing behind them, so they take the first one ahead. This is
            !> the case a streaming classifier cannot serve at all, and it is
            !> why a PWB run used to open on covariance maximisation.
            do j = 1, nxt - 1
                i = idx(j)
                if (trim(PwbTimelagCache(i)%result%reliability_class) /= 'pending') cycle
                PwbTimelagCache(i)%used_lag = lag(nxt)
                PwbTimelagCache(i)%result%reliability_class = 'S3_backfilled'
                PwbTimelagCache(i)%result%fill_method = 'backfilled'
                PwbTimelagCache(i)%result%fallback_source = 'backfilled'
                PwbTimelagCache(i)%result%origin_gas = gas
            end do
            !> And anything after the last reliable detection carries it
            !> forward; there is nothing ahead of it to interpolate towards.
            prev = 0
            do j = n, 1, -1
                if (settled(j)) then
                    prev = j
                    exit
                end if
            end do
            do j = prev + 1, n
                i = idx(j)
                if (trim(PwbTimelagCache(i)%result%reliability_class) /= 'pending') cycle
                PwbTimelagCache(i)%used_lag = lag(prev)
                PwbTimelagCache(i)%result%reliability_class = 'S3_carryforward'
                PwbTimelagCache(i)%result%fill_method = 'carryforward'
                PwbTimelagCache(i)%result%fallback_source = 'S3_carryforward'
                PwbTimelagCache(i)%result%origin_gas = gas
            end do
        end if

        !> Step 7: with no reliable detection anywhere, the median of what was
        !> detected at all is still better than the nominal window.
        nsel = 0
        do j = 1, n
            i = idx(j)
            if (PwbTimelagCache(i)%result%hdi_prefiltered) cycle
            if (PwbTimelagCache(i)%result%selected_lag == error) cycle
            if (PwbTimelagCache(i)%result%edge_pinned) cycle
            nsel = nsel + 1
            sorted(nsel) = PwbTimelagCache(i)%result%selected_lag
        end do
        median_lag = error
        if (nsel > 0) median_lag = Median(sorted, nsel)
        do j = 1, n
            i = idx(j)
            if (trim(PwbTimelagCache(i)%result%reliability_class) /= 'pending') cycle
            if (median_lag /= error) then
                PwbTimelagCache(i)%used_lag = median_lag
                PwbTimelagCache(i)%result%reliability_class = 'S3_median'
                PwbTimelagCache(i)%result%fill_method = 'median'
                PwbTimelagCache(i)%result%fallback_source = 'median'
                PwbTimelagCache(i)%result%origin_gas = gas
            else
                !> Step 8: nothing was detected for this gas in the whole run.
                !> Whatever the pass settled on for this period stands.
                PwbTimelagCache(i)%used_lag = fallback_lag(j)
                PwbTimelagCache(i)%result%reliability_class = 'fallback'
                PwbTimelagCache(i)%result%fill_method = 'maxcov_default'
                PwbTimelagCache(i)%result%fallback_source = 'maxcov_default'
                PwbTimelagCache(i)%result%fallback_used = .true.
            end if
        end do
    end do

    !> Row lag, applied lag and the "default used" flag follow the settled
    !> time-lag, so the production pass reads one consistent set of numbers.
    do i = 1, PwbTimelagCacheN
        if (PwbTimelagCache(i)%used_lag == error) cycle
        PwbTimelagCache(i)%row_lag = nint(PwbTimelagCache(i)%used_lag * Metadata%ac_freq)
        PwbTimelagCache(i)%actual_lag = PwbTimelagCache(i)%result%selected_lag
        if (PwbTimelagCache(i)%actual_lag == error) &
            PwbTimelagCache(i)%actual_lag = PwbTimelagCache(i)%used_lag
        PwbTimelagCache(i)%default_used = &
            trim(PwbTimelagCache(i)%result%fill_method) == 'maxcov_default'
        PwbTimelagCache(i)%result%applied_lag = PwbTimelagCache(i)%used_lag
        PwbTimelagCache(i)%result%applied_row_lag = PwbTimelagCache(i)%row_lag
    end do
    PwbCacheDirty = .true.

    !> The tallies the run log prints describe the settled table, not the
    !> streaming guesses that produced it.
    call ResetPwbDiagnostics()
    do i = 1, PwbTimelagCacheN
        call CountPwbDiagnostic(PwbTimelagCache(i)%gas, PwbTimelagCache(i)%result)
    end do

    deallocate(ord, tmin, idx, lag, fallback_lag, settled, sorted)
end subroutine PostProcessPwbTimelagCache

!***************************************************************************
!> Two gas slots measured by the same physical analyser.
!>
!> The instrument name is the identity where a project states one. The model
!> string is a fallback for projects that do not, and cannot tell two
!> analysers of the same model apart - which is why matching on it alone let
!> two LI-7200s at one site share a time lag across two different tubes.
!***************************************************************************
logical function SameAnalyser(a, b)
    integer, intent(in) :: a, b

    if (len_trim(E2Col(a)%instr_name) > 0 .and. len_trim(E2Col(b)%instr_name) > 0) then
        SameAnalyser = E2Col(a)%instr_name == E2Col(b)%instr_name
    else
        SameAnalyser = E2Col(a)%instr%model == E2Col(b)%instr%model
    end if
end function SameAnalyser

!***************************************************************************
!> Minutes since 1970-01-01, from 'yyyy-mm-dd' and 'HH:MM'.
!>
!> A real time axis, so that interpolation across a gap in the raw files
!> spans the time the gap actually lasted rather than the number of rows it
!> happens to occupy.
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

subroutine SortByMinutes(ord, tmin, n)
    integer, intent(in) :: n
    integer, intent(inout) :: ord(n)
    integer(8), intent(in) :: tmin(n)
    integer :: i, j, key

    do i = 2, n
        key = ord(i)
        j = i - 1
        do while (j >= 1)
            if (tmin(ord(j)) <= tmin(key)) exit
            ord(j+1) = ord(j)
            j = j - 1
        end do
        ord(j+1) = key
    end do
end subroutine SortByMinutes

real(kind = dbl) function Median(x, n)
    integer, intent(in) :: n
    real(kind = dbl), intent(inout) :: x(n)
    integer :: i, j
    real(kind = dbl) :: tmp

    do i = 2, n
        tmp = x(i)
        j = i - 1
        do while (j >= 1)
            if (x(j) <= tmp) exit
            x(j+1) = x(j)
            j = j - 1
        end do
        x(j+1) = tmp
    end do
    if (mod(n, 2) == 1) then
        Median = x((n + 1) / 2)
    else
        Median = 0.5d0 * (x(n/2) + x(n/2 + 1))
    end if
end function Median

subroutine PwbDetectGas(Set, nrow, ncol, gas, LocResult, success)
    use m_index_parameters
    use m_log
    implicit none
    integer, intent(in) :: nrow, ncol, gas
    real(kind = dbl), intent(in) :: Set(nrow, ncol)
    type(PWBResultType), intent(out) :: LocResult
    logical, intent(out) :: success

    integer :: min_rl, max_rl, lag, margin, eval_lo, eval_hi
    integer :: nvalid_w, nvalid_t, nvalid_s
    integer :: p_scalar, p_w, p_t
    real(kind = dbl) :: min_valid
    real(kind = dbl), allocatable :: ww(:), tt(:), ss(:)
    real(kind = dbl), allocatable :: phi_s(:), phi_w(:), phi_t(:)
    real(kind = dbl), allocatable :: s_fs(:), w_fs(:), t_fs(:)
    real(kind = dbl), allocatable :: s_fw(:), w_fw(:)
    real(kind = dbl), allocatable :: s_ft(:), t_ft(:)
    real(kind = dbl), allocatable :: raw_ccov(:)
    type(PWBResultType) :: candidate(4)
    character(2) :: combo(4)
    logical :: ok(4)

    call InitPwbResult(LocResult)
    success = .false.
    if (gas < firstGas .or. gas > lastGas) return
    if (.not. E2Col(gas)%present .or. .not. E2Col(ts)%present) then
        LocResult%fallback_used = .true.
        return
    end if

    min_rl = nint(PWBSetup%min_lag(gas) * Metadata%ac_freq)
    max_rl = nint(PWBSetup%max_lag(gas) * Metadata%ac_freq)
    LocResult%effective_min_lag = dble(min_rl) / Metadata%ac_freq
    LocResult%effective_max_lag = dble(max_rl) / Metadata%ac_freq
    if (.not. pwb_bounds_warned(gas) .and. E2Col(gas)%instr%path_type == 'closed' &
        .and. PWBSetup%min_lag(gas) < 0d0 .and. PWBSetup%max_lag(gas) > 0d0 &
        .and. PWBSetup%max_lag(gas) - PWBSetup%min_lag(gas) > 20d0) then
        write(*, '(a,a,a,f8.2,a,f8.2,a)') '  WARNING: broad symmetric PWB lag window for ', &
            trim(GasLabel(gas)), ' [', PWBSetup%min_lag(gas), ', ', &
            PWBSetup%max_lag(gas), '] s on a closed-path gas; consider physical positive bounds.'
        write(ulog, '(a,a,a,f8.2,a,f8.2,a)') '  WARNING: broad symmetric PWB lag window for ', &
            trim(GasLabel(gas)), ' [', PWBSetup%min_lag(gas), ', ', &
            PWBSetup%max_lag(gas), '] s on a closed-path gas; consider physical positive bounds.'
        pwb_bounds_warned(gas) = .true.
    end if
    if (min_rl >= max_rl) then
        LocResult%fallback_used = .true.
        return
    end if

    !> The cross-correlation is evaluated a guard band beyond the declared
    !> window, not over the mirrored symmetric range the window sits inside.
    !>
    !> Half the smoothing width is what the centred rolling mean consumes, so
    !> every lag inside the search range keeps a genuine smoothed value rather
    !> than one carried in from the edge. The two seconds on top of that are
    !> there because the window is a guess: a peak just outside it used to be
    !> indistinguishable from an ordinary edge-pinned failure, and is now
    !> reported. The applied lag is still taken from inside the declared
    !> window - the guard band informs, it does not widen the search.
    margin = max(max(1, PWBSetup%smoothing_width) / 2, nint(2d0 * Metadata%ac_freq))
    eval_lo = max(min_rl - margin, -(nrow - 2))
    eval_hi = min(max_rl + margin, nrow - 2)
    if (eval_lo >= eval_hi) then
        LocResult%fallback_used = .true.
        return
    end if

    allocate(ww(nrow), tt(nrow), ss(nrow))
    ww = Set(:, w)
    tt = Set(:, ts)
    ss = Set(:, gas)

    nvalid_w = count(ww /= error)
    nvalid_t = count(tt /= error)
    nvalid_s = count(ss /= error)
    min_valid = max(0d0, min(1d0, PWBSetup%min_valid_frac)) * dble(nrow)
    if (dble(nvalid_w) < min_valid .or. dble(nvalid_t) < min_valid .or. dble(nvalid_s) < min_valid) then
        LocResult%fallback_used = .true.
        deallocate(ww, tt, ss)
        return
    end if

    call FillMissingLinear(ww, nrow)
    call FillMissingLinear(tt, nrow)
    call FillMissingLinear(ss, nrow)

    if (.not. IsStationary(ww, nrow) .or. .not. IsStationary(tt, nrow) .or. .not. IsStationary(ss, nrow)) then
        call DifferenceSeries(ww, nrow)
        call DifferenceSeries(tt, nrow)
        call DifferenceSeries(ss, nrow)
    end if

    call FitArAic(ss, nrow, phi_s, p_scalar)
    call FitArAic(ww, nrow, phi_w, p_w)
    call FitArAic(tt, nrow, phi_t, p_t)

    allocate(s_fs(nrow), w_fs(nrow), t_fs(nrow))
    allocate(s_fw(nrow), w_fw(nrow), s_ft(nrow), t_ft(nrow))
    call ApplyArFilter(ss, nrow, phi_s, p_scalar, s_fs)
    call ApplyArFilter(ww, nrow, phi_s, p_scalar, w_fs)
    call ApplyArFilter(tt, nrow, phi_s, p_scalar, t_fs)
    call ApplyArFilter(ss, nrow, phi_w, p_w, s_fw)
    call ApplyArFilter(ww, nrow, phi_w, p_w, w_fw)
    call ApplyArFilter(ss, nrow, phi_t, p_t, s_ft)
    call ApplyArFilter(tt, nrow, phi_t, p_t, t_ft)

    call EnsurePwbScratch(nrow, eval_lo, eval_hi, max(1, PWBSetup%n_bootstrap))

    combo = (/'cw', 'wc', 'ct', 'tc'/)
    call RunPwbCombination(w_fs, s_fs, nrow, min_rl, max_rl, eval_lo, eval_hi, &
        gas, combo(1), candidate(1), ok(1))
    call RunPwbCombination(w_fw, s_fw, nrow, min_rl, max_rl, eval_lo, eval_hi, &
        gas, combo(2), candidate(2), ok(2))
    call RunPwbCombination(t_fs, s_fs, nrow, min_rl, max_rl, eval_lo, eval_hi, &
        gas, combo(3), candidate(3), ok(3))
    call RunPwbCombination(t_ft, s_ft, nrow, min_rl, max_rl, eval_lo, eval_hi, &
        gas, combo(4), candidate(4), ok(4))

    call SelectBestCandidate(candidate, ok, LocResult, success)
    if (LocResult%peak_outside_window) pwb_outside_window(gas) = pwb_outside_window(gas) + 1
    if (success) then
        allocate(raw_ccov(min_rl:max_rl))
        call ComputeCcovWindow(ww, ss, nrow, min_rl, max_rl, raw_ccov)
        lag = LocResult%row_lag
        if (lag >= min_rl .and. lag <= max_rl) LocResult%raw_covariance = raw_ccov(lag)
        deallocate(raw_ccov)
    end if

    deallocate(ww, tt, ss)
    if (allocated(phi_s)) deallocate(phi_s)
    if (allocated(phi_w)) deallocate(phi_w)
    if (allocated(phi_t)) deallocate(phi_t)
    deallocate(s_fs, w_fs, t_fs, s_fw, w_fw, s_ft, t_ft)
end subroutine PwbDetectGas

subroutine InitPwbResult(res)
    type(PWBResultType), intent(out) :: res
    res%selected_lag = error
    res%row_lag = 0
    res%applied_lag = error
    res%applied_row_lag = 0
    res%hdi_low = error
    res%hdi_high = error
    res%hdi_range = error
    res%reliability_class = 'failed'
    res%best_combination = '--'
    res%fallback_source = 'none'
    res%fill_method = 'none'
    res%donor_gas = 'none'
    res%origin_gas = 0
    res%edge_pinned = .false.
    res%fallback_used = .false.
    res%block_length_clamped = .false.
    res%effective_min_lag = error
    res%effective_max_lag = error
    res%effective_block_length_s = error
    res%raw_covariance = error
    res%ccf_at_mode = 0d0
    res%unrestricted_peak_lag = error
    res%peak_outside_window = .false.
    res%hdi_prefiltered = .false.
end subroutine InitPwbResult

subroutine FillMissingLinear(x, n)
    integer, intent(in) :: n
    real(kind = dbl), intent(inout) :: x(n)
    integer :: i, j, k
    real(kind = dbl) :: x0, x1

    j = 0
    do i = 1, n
        if (x(i) /= error) then
            j = i
            exit
        end if
    end do
    if (j == 0) return
    if (j > 1) x(1:j-1) = x(j)

    i = j + 1
    do while (i <= n)
        if (x(i) /= error) then
            i = i + 1
        else
            j = i - 1
            k = i
            do
                if (k > n) exit
                if (x(k) /= error) exit
                k = k + 1
            end do
            if (k > n) then
                x(i:n) = x(j)
                exit
            end if
            x0 = x(j)
            x1 = x(k)
            do i = j + 1, k - 1
                x(i) = x0 + (x1 - x0) * dble(i - j) / dble(k - j)
            end do
            i = k + 1
        end if
    end do
end subroutine FillMissingLinear

logical function IsStationary(x, n)
    integer, intent(in) :: n
    real(kind = dbl), intent(in) :: x(n)
    integer :: i
    real(kind = dbl) :: meanx, sse, cum, rho
    real(kind = dbl), parameter :: cv_1pct = 0.00537748023783321d0

    meanx = sum(x) / dble(n)
    sse = 0d0
    cum = 0d0
    rho = 0d0
    do i = 1, n
        sse = sse + (x(i) - meanx)**2
    end do
    if (sse <= 0d0) then
        IsStationary = .true.
        return
    end if
    do i = 1, n
        cum = cum + x(i) - meanx
        rho = rho + cum**2
    end do
    rho = rho / (dble(n)**2 * sse)
    IsStationary = rho < cv_1pct
end function IsStationary

subroutine DifferenceSeries(x, n)
    integer, intent(in) :: n
    real(kind = dbl), intent(inout) :: x(n)
    integer :: i
    do i = n, 2, -1
        x(i) = x(i) - x(i-1)
    end do
    x(1) = 0d0
end subroutine DifferenceSeries

!***************************************************************************
!> AR(p) by AIC, with p searched up to floor(100*log10(N)).
!>
!> There was a cap on that search, offered as a speed option. It saved very
!> little - the autocovariance pass here is a percent or so of a PWB period,
!> against the bootstrap cross-correlations - and it works against the whole
!> point of pre-whitening: an order below the optimum leaves autocorrelation
!> in the residuals, which is exactly what broadens the CCF peak the method
!> exists to sharpen. The reference says so in as many words.
!***************************************************************************
subroutine FitArAic(x, n, phi_best, p_best)
    integer, intent(in) :: n
    real(kind = dbl), intent(in) :: x(n)
    real(kind = dbl), allocatable, intent(out) :: phi_best(:)
    integer, intent(out) :: p_best
    integer :: max_lag, p, i
    real(kind = dbl) :: meanx, sigma2, kappa, best_aic, aic
    real(kind = dbl), allocatable :: acf(:), phi(:), phi_old(:)

    max_lag = min(int(floor(100d0 * log10(dble(max(n, 2))))), n - 1)
    if (max_lag < 1) then
        allocate(phi_best(0))
        p_best = 0
        return
    end if
    allocate(acf(0:max_lag))
    meanx = sum(x) / dble(n)
    do p = 0, max_lag
        acf(p) = 0d0
        do i = 1, n - p
            acf(p) = acf(p) + (x(i) - meanx) * (x(i+p) - meanx)
        end do
        acf(p) = acf(p) / dble(n)
    end do
    if (acf(0) <= 0d0) then
        allocate(phi_best(0))
        p_best = 0
        deallocate(acf)
        return
    end if

    allocate(phi(1:max_lag), phi_old(1:max_lag), phi_best(0))
    p_best = 0
    best_aic = dble(n) * log(acf(0))
    sigma2 = acf(0)
    phi = 0d0
    do p = 1, max_lag
        if (p == 1) then
            kappa = acf(1) / sigma2
        else
            kappa = (acf(p) - dot_product(phi(1:p-1), acf(p-1:1:-1))) / sigma2
        end if
        phi_old = phi
        if (p > 1) phi(1:p-1) = phi_old(1:p-1) - kappa * phi_old(p-1:1:-1)
        phi(p) = kappa
        sigma2 = sigma2 * (1d0 - kappa**2)
        if (sigma2 <= 0d0) exit
        aic = dble(n) * log(sigma2) + 2d0 * dble(p)
        if (aic < best_aic) then
            best_aic = aic
            p_best = p
            if (allocated(phi_best)) deallocate(phi_best)
            allocate(phi_best(p_best))
            phi_best = phi(1:p_best)
        end if
    end do
    deallocate(acf, phi, phi_old)
end subroutine FitArAic

subroutine ApplyArFilter(x, n, phi, p, y)
    integer, intent(in) :: n, p
    real(kind = dbl), intent(in) :: x(n)
    real(kind = dbl), intent(in) :: phi(:)
    real(kind = dbl), intent(out) :: y(n)
    integer :: i, j
    real(kind = dbl) :: meanx

    meanx = sum(x) / dble(n)
    y = x - meanx
    if (p <= 0) return
    do i = n, 1, -1
        y(i) = x(i) - meanx
        do j = 1, min(p, i - 1)
            y(i) = y(i) - phi(j) * (x(i-j) - meanx)
        end do
        if (i <= p) y(i) = 0d0
    end do
end subroutine ApplyArFilter

!> Size the shared scratch to this gas's period, reallocating only on change.
subroutine EnsurePwbScratch(n, lo, hi, nboot)
    integer, intent(in) :: n, lo, hi, nboot

    if (allocated(sc_xb) .and. sc_n == n .and. sc_lo == lo &
        .and. sc_hi == hi .and. sc_nboot == nboot) return
    if (allocated(sc_xb)) deallocate(sc_xb, sc_yb, sc_xc, sc_yc)
    if (allocated(sc_ccf)) deallocate(sc_ccf, sc_smooth, sc_mean_ccf, sc_mean_smooth)
    if (allocated(sc_boot)) deallocate(sc_boot, sc_hdi)
    allocate(sc_xb(n), sc_yb(n), sc_xc(n), sc_yc(n))
    allocate(sc_ccf(lo:hi), sc_smooth(lo:hi), sc_mean_ccf(lo:hi), sc_mean_smooth(lo:hi))
    allocate(sc_boot(nboot), sc_hdi(nboot))
    sc_n = n
    sc_lo = lo
    sc_hi = hi
    sc_nboot = nboot
end subroutine EnsurePwbScratch

subroutine RunPwbCombination(x, y, n, min_rl, max_rl, eval_lo, eval_hi, gas, combo, res, ok)
    integer, intent(in) :: n, min_rl, max_rl, eval_lo, eval_hi, gas
    real(kind = dbl), intent(in) :: x(n), y(n)
    character(2), intent(in) :: combo
    type(PWBResultType), intent(out) :: res
    logical, intent(out) :: ok
    integer :: b, i, pos, block_len, nblocks, start
    integer :: requested_block_len, widest
    integer :: nboot, lag, best_idx, unrestricted_idx
    integer(8) :: state

    call InitPwbResult(res)
    res%best_combination = combo
    ok = .false.
    nboot = max(1, PWBSetup%n_bootstrap)
    widest = max(abs(min_rl), abs(max_rl))
    requested_block_len = nint(PWBSetup%block_length_s * Metadata%ac_freq)
    if (requested_block_len <= 0) requested_block_len = max(1, 2 * widest)
    block_len = requested_block_len
    res%block_length_clamped = .false.
    if (requested_block_len < 2 * widest .and. .not. pwb_block_warned(gas)) then
        write(*, '(a,a,a,f8.2,a,f8.2,a)') '  WARNING: PWB block length for ', &
            trim(GasLabel(gas)), ' (', &
            dble(requested_block_len) / Metadata%ac_freq, ' s) is shorter than 2*lag_max (', &
            dble(2 * widest) / Metadata%ac_freq, ' s).'
        write(ulog, '(a,a,a,f8.2,a,f8.2,a)') '  WARNING: PWB block length for ', &
            trim(GasLabel(gas)), ' (', &
            dble(requested_block_len) / Metadata%ac_freq, ' s) is shorter than 2*lag_max (', &
            dble(2 * widest) / Metadata%ac_freq, ' s).'
        pwb_block_warned(gas) = .true.
    end if
    block_len = min(max(1, block_len), n)
    res%effective_block_length_s = dble(block_len) / Metadata%ac_freq
    nblocks = (n + block_len - 1) / block_len

    res%effective_min_lag = dble(min_rl) / Metadata%ac_freq
    res%effective_max_lag = dble(max_rl) / Metadata%ac_freq
    sc_mean_ccf = 0d0
    state = PwbStreamSeed(gas, combo)

    do b = 1, nboot
        pos = 1
        do i = 1, nblocks
            start = 1 + RandBelow(state, max(1, n - block_len + 1))
            call CopyBlock(x, y, n, start, block_len, sc_xb, sc_yb, pos)
            if (pos > n) exit
        end do
        call ComputeCcfWindow(sc_xb, sc_yb, n, eval_lo, eval_hi, sc_ccf)
        call SmoothAndFill(sc_ccf, eval_lo, eval_hi, max(1, PWBSetup%smoothing_width), sc_smooth)
        best_idx = ArgmaxAbs(sc_smooth(min_rl:max_rl), min_rl, max_rl)
        sc_boot(b) = best_idx
        sc_mean_ccf = sc_mean_ccf + sc_ccf
    end do
    sc_mean_ccf = sc_mean_ccf / dble(nboot)
    call SmoothAndFill(sc_mean_ccf, eval_lo, eval_hi, max(1, PWBSetup%smoothing_width), sc_mean_smooth)

    lag = MapLagEstimate(sc_boot, nboot)
    do i = 1, nboot
        sc_hdi(i) = dble(sc_boot(i)) / Metadata%ac_freq
    end do
    call Hdi95(sc_hdi, nboot, res%hdi_low, res%hdi_high)
    res%row_lag = lag
    res%selected_lag = dble(lag) / Metadata%ac_freq
    res%hdi_range = res%hdi_high - res%hdi_low
    res%edge_pinned = lag == min_rl .or. lag == max_rl
    res%ccf_at_mode = abs(sc_mean_smooth(lag))

    !> What the guard band saw. The applied lag is the restricted one above;
    !> this only reports whether the declared window was where the signal is.
    unrestricted_idx = ArgmaxAbs(sc_mean_smooth, eval_lo, eval_hi)
    res%unrestricted_peak_lag = dble(unrestricted_idx) / Metadata%ac_freq
    res%peak_outside_window = unrestricted_idx < min_rl .or. unrestricted_idx > max_rl

    res%reliability_class = 'detected'
    ok = .not. res%edge_pinned
end subroutine RunPwbCombination

integer function MapLagEstimate(samples, n)
    integer, intent(in) :: n
    integer, intent(in) :: samples(n)
    integer :: i, grid, lo, hi, best
    real(kind = dbl) :: mean_s, var_s, sd_s, bw, dens, best_dens, z

    lo = minval(samples)
    hi = maxval(samples)
    if (lo == hi) then
        MapLagEstimate = lo
        return
    end if

    mean_s = sum(dble(samples)) / dble(n)
    var_s = 0d0
    do i = 1, n
        var_s = var_s + (dble(samples(i)) - mean_s)**2
    end do
    sd_s = sqrt(max(0d0, var_s / max(1d0, dble(n - 1))))
    bw = max(1d0, 1.06d0 * sd_s * dble(n)**(-0.2d0))

    best = lo
    best_dens = -1d0
    do grid = lo, hi
        dens = 0d0
        do i = 1, n
            z = (dble(grid) - dble(samples(i))) / bw
            dens = dens + exp(-0.5d0 * z * z)
        end do
        if (dens > best_dens) then
            best_dens = dens
            best = grid
        end if
    end do
    MapLagEstimate = best
end function MapLagEstimate

subroutine CopyBlock(x, y, n, start, block_len, xb, yb, pos)
    integer, intent(in) :: n, start, block_len
    integer, intent(inout) :: pos
    real(kind = dbl), intent(in) :: x(n), y(n)
    real(kind = dbl), intent(inout) :: xb(n), yb(n)
    integer :: j, src
    do j = 0, block_len - 1
        if (pos > n) exit
        src = min(n, start + j)
        xb(pos) = x(src)
        yb(pos) = y(src)
        pos = pos + 1
    end do
end subroutine CopyBlock

!***************************************************************************
!> The bootstrap stream for one gas and one pre-whitening combination in one
!> averaging period.
!>
!> The period is in the seed. It was not, so every half-hour drew the very
!> same sequence of block start positions for a given gas: the resampling was
!> one fixed template applied to the whole run rather than an independent
!> draw per period, and the spread it produced - which is the HDI, which is
!> what S1 is decided on - was not an independent sample of anything. Mixing
!> the timestamp in keeps a period reproducible, which the stored table
!> depends on, while making periods independent of one another.
!***************************************************************************
integer(8) function PwbStreamSeed(gas, combo)
    integer, intent(in) :: gas
    character(2), intent(in) :: combo
    integer :: i

    PwbStreamSeed = int(PWBSetup%random_seed, 8)
    call MixIn(PwbStreamSeed, int(gas, 8))
    do i = 1, len(combo)
        call MixIn(PwbStreamSeed, int(iachar(combo(i:i)), 8))
    end do
    do i = 1, len_trim(PwbPeriodDate)
        call MixIn(PwbStreamSeed, int(iachar(PwbPeriodDate(i:i)), 8))
    end do
    do i = 1, len_trim(PwbPeriodTime)
        call MixIn(PwbStreamSeed, int(iachar(PwbPeriodTime(i:i)), 8))
    end do
    if (PwbStreamSeed == 0_8) PwbStreamSeed = 88172645463325252_8
end function PwbStreamSeed

!> One xorshift round over an absorbed value. Shifts and exclusive-or only,
!> so there is no multiplication to overflow.
subroutine MixIn(state, val)
    integer(8), intent(inout) :: state
    integer(8), intent(in) :: val

    state = ieor(state, val)
    state = ieor(state, ishft(state, 13))
    state = ieor(state, ishft(state, -7))
    state = ieor(state, ishft(state, 17))
end subroutine MixIn

!***************************************************************************
!> A uniform integer in [0, upper).
!>
!> This was a linear congruential generator carrying glibc's multiplier and
!> increment - which are chosen for modulus 2^31 - applied with modulus
!> 2^31-1, and then reduced with mod() straight off the low bits. The low
!> bits of any LCG are its worst, and a bootstrap replicate here draws about
!> ninety block starts in a row from them; correlated starts make replicates
!> resemble one another, which narrows the HDI and hands out S1 to detections
!> that have not earned it.
!>
!> xorshift64 has no multiplication, so nothing can overflow, and the value
!> is taken from the top bits. ISHFT is a logical shift, so the shifted
!> result is non-negative whatever the sign of the state.
!***************************************************************************
integer function RandBelow(state, upper)
    integer(8), intent(inout) :: state
    integer, intent(in) :: upper

    state = ieor(state, ishft(state, 13))
    state = ieor(state, ishft(state, -7))
    state = ieor(state, ishft(state, 17))
    RandBelow = int(mod(ishft(state, -33), int(max(1, upper), 8)), 4)
end function RandBelow

!***************************************************************************
!> Normalised cross-correlation over [min_rl, max_rl].
!>
!> There was an option to skip the normalisation. It saved two passes over
!> the series against one per lag - well under a percent - and the divisor is
!> constant across lags, so within one combination it changed nothing. Across
!> combinations it changed everything: SelectBestCandidate compares the four
!> at their modes, and the four pair the scalar against vertical wind and
!> against sonic temperature, whose covariances carry different physical
!> units. Unnormalised, the winner was decided by that unit scale.
!***************************************************************************
subroutine ComputeCcfWindow(x, y, n, min_rl, max_rl, ccf)
    integer,  intent(in) :: n, min_rl, max_rl
    real(kind = dbl), intent(in)  :: x(n), y(n)
    real(kind = dbl), intent(out) :: ccf(min_rl:max_rl)
    integer :: lag, i, nn
    real(kind = dbl) :: mx, my, vx, vy, denom, cov

    mx = sum(x) / dble(n)
    my = sum(y) / dble(n)
    sc_xc = x - mx
    sc_yc = y - my
    vx = sum(sc_xc * sc_xc)
    vy = sum(sc_yc * sc_yc)
    denom = sqrt(vx * vy)
    if (denom <= 0d0) then
        ccf = 0d0
        return
    end if

    do lag = min_rl, max_rl
        nn = n - abs(lag)
        if (nn <= 1) then
            ccf(lag) = 0d0
            cycle
        end if
        cov = 0d0
        if (lag >= 0) then
            do i = 1, nn
                cov = cov + sc_xc(i) * sc_yc(i + lag)
            end do
        else
            do i = 1, nn
                cov = cov + sc_xc(i - lag) * sc_yc(i)
            end do
        end if
        ccf(lag) = cov / denom
    end do
end subroutine ComputeCcfWindow

subroutine ComputeCcovWindow(x, y, n, min_rl, max_rl, ccov)
    integer, intent(in) :: n, min_rl, max_rl
    real(kind = dbl), intent(in) :: x(n), y(n)
    real(kind = dbl), intent(out) :: ccov(min_rl:max_rl)
    integer :: lag, i, nn, s1, s2
    real(kind = dbl) :: mx, my
    do lag = min_rl, max_rl
        nn = n - abs(lag)
        if (nn <= 1) then
            ccov(lag) = error
            cycle
        end if
        mx = 0d0; my = 0d0
        do i = 1, nn
            if (lag >= 0) then
                s1 = i; s2 = i + lag
            else
                s1 = i - lag; s2 = i
            end if
            mx = mx + x(s1)
            my = my + y(s2)
        end do
        mx = mx / dble(nn); my = my / dble(nn)
        ccov(lag) = 0d0
        do i = 1, nn
            if (lag >= 0) then
                s1 = i; s2 = i + lag
            else
                s1 = i - lag; s2 = i
            end if
            ccov(lag) = ccov(lag) + (x(s1) - mx) * (y(s2) - my)
        end do
        ccov(lag) = ccov(lag) / dble(nn)
    end do
end subroutine ComputeCcovWindow

!***************************************************************************
!> Centred rolling mean, with the edges carried in from the nearest computed
!> value (the reference's na.locf).
!>
!> The divisor is the number of terms actually summed. It was the requested
!> width, while the loop sums 2*(width/2)+1 of them - so an even width, which
!> the interface allowed, scaled every smoothed value by (width+1)/width. A
!> uniform scale does not move an argmax, but it is still not the mean it
!> claims to be, and it fed the cross-combination comparison.
!***************************************************************************
subroutine SmoothAndFill(x, min_rl, max_rl, width, y)
    integer, intent(in) :: min_rl, max_rl, width
    real(kind = dbl), intent(in) :: x(min_rl:max_rl)
    real(kind = dbl), intent(out) :: y(min_rl:max_rl)
    integer :: i, j, half, first_valid, last_valid
    half = width / 2
    first_valid = min_rl + half
    last_valid = max_rl - half
    do i = first_valid, last_valid
        y(i) = 0d0
        do j = i - half, i + half
            y(i) = y(i) + x(j)
        end do
        y(i) = y(i) / dble(2 * half + 1)
    end do
    if (first_valid <= last_valid) then
        y(min_rl:first_valid - 1) = y(first_valid)
        y(last_valid + 1:max_rl) = y(last_valid)
    else
        y = x
    end if
end subroutine SmoothAndFill

integer function ArgmaxAbs(x, min_rl, max_rl)
    integer, intent(in) :: min_rl, max_rl
    real(kind = dbl), intent(in) :: x(min_rl:max_rl)
    integer :: i
    ArgmaxAbs = min_rl
    do i = min_rl, max_rl
        if (abs(x(i)) > abs(x(ArgmaxAbs))) ArgmaxAbs = i
    end do
end function ArgmaxAbs

subroutine Hdi95(x, n, lo, hi)
    integer, intent(in) :: n
    real(kind = dbl), intent(inout) :: x(n)
    real(kind = dbl), intent(out) :: lo, hi
    integer :: i, j, m, best
    real(kind = dbl) :: tmp, width
    do i = 2, n
        tmp = x(i)
        j = i - 1
        do while (j >= 1)
            if (x(j) <= tmp) exit
            x(j+1) = x(j)
            j = j - 1
        end do
        x(j+1) = tmp
    end do
    m = max(1, int(floor(0.95d0 * dble(n))))
    if (m >= n) then
        lo = x(1); hi = x(n); return
    end if
    best = 1
    width = x(1+m) - x(1)
    do i = 2, n - m
        if (x(i+m) - x(i) < width) then
            best = i
            width = x(i+m) - x(i)
        end if
    end do
    lo = x(best)
    hi = x(best + m)
end subroutine Hdi95

subroutine SelectBestCandidate(candidate, ok, res, success)
    type(PWBResultType), intent(in) :: candidate(4)
    logical, intent(in) :: ok(4)
    type(PWBResultType), intent(out) :: res
    logical, intent(out) :: success
    integer :: i, best

    call InitPwbResult(res)
    success = .false.

    !> Highest |mean smoothed CCF| at the mode lag, as in the reference.
    !>
    !> Deliberate deviation: a candidate whose mode did not land on the window
    !> edge is preferred before magnitude is consulted at all. The reference
    !> picks on magnitude alone and only then asks whether the winner is
    !> edge-pinned, which throws the period away when an unpinned candidate
    !> was available. Where no candidate is unpinned the two agree.
    best = 0
    do i = 1, 4
        if (ok(i)) then
            if (best == 0) then
                best = i
            elseif (candidate(i)%ccf_at_mode > candidate(best)%ccf_at_mode) then
                best = i
            end if
        end if
    end do
    if (best == 0) then
        do i = 1, 4
            if (best == 0) then
                best = i
            elseif (candidate(i)%ccf_at_mode > candidate(best)%ccf_at_mode) then
                best = i
            end if
        end do
    end if
    if (best > 0) then
        res = candidate(best)
        success = .not. res%edge_pinned
    end if
end subroutine SelectBestCandidate

subroutine CountPwbDiagnostic(gas, res)
    integer, intent(in) :: gas
    type(PWBResultType), intent(in) :: res

    if (gas < firstGas .or. gas > lastGas) return
    pwb_attempts(gas) = pwb_attempts(gas) + 1
    if (res%fallback_used .or. trim(res%reliability_class) == 'fallback') then
        pwb_fallbacks(gas) = pwb_fallbacks(gas) + 1
        select case(trim(res%fallback_source))
            case('maxcov_default')
                pwb_fallback_maxcov(gas) = pwb_fallback_maxcov(gas) + 1
            case('nominal/default')
                pwb_fallback_nominal(gas) = pwb_fallback_nominal(gas) + 1
            case default
                pwb_fallback_other(gas) = pwb_fallback_other(gas) + 1
        end select
    elseif (trim(res%reliability_class) == 'S4_instrument_shared') then
        pwb_instrument_shared(gas) = pwb_instrument_shared(gas) + 1
    elseif (index(res%reliability_class, 'S3_') == 1) then
        !> Every gap-filled arm - interpolated, back-filled, carried forward,
        !> median - counts here, which is what "not detected in this period"
        !> means to a reader of the summary.
        pwb_carryforwards(gas) = pwb_carryforwards(gas) + 1
    else
        pwb_successes(gas) = pwb_successes(gas) + 1
    end if
end subroutine CountPwbDiagnostic

subroutine ReportPwbDiagnostics()
    integer :: gas
    integer :: total_attempts, total_successes, total_carryforwards, total_fallbacks
    integer :: total_instrument_shared
    integer :: total_fallback_maxcov, total_fallback_nominal, total_fallback_other

    total_attempts = sum(pwb_attempts(firstGas:lastGas))
    if (total_attempts == 0) return

    total_successes = sum(pwb_successes(firstGas:lastGas))
    total_instrument_shared = sum(pwb_instrument_shared(firstGas:lastGas))
    total_carryforwards = sum(pwb_carryforwards(firstGas:lastGas))
    total_fallbacks = sum(pwb_fallbacks(firstGas:lastGas))
    total_fallback_maxcov = sum(pwb_fallback_maxcov(firstGas:lastGas))
    total_fallback_nominal = sum(pwb_fallback_nominal(firstGas:lastGas))
    total_fallback_other = sum(pwb_fallback_other(firstGas:lastGas))

    write(*, '(a)')
    write(ulog, '(a)')
    call LogSay(' PWB time-lag detection summary:')
    do gas = firstGas, lastGas
        if (pwb_attempts(gas) > 0) then
            write(*, '(a, a, a, i0, a, i0, a, i0, a, i0, a, i0, a, i0, a, i0, a, i0, a)') &
                '  ', trim(GasLabel(gas)), &
                ': attempts=', pwb_attempts(gas), &
                ', S1/S2=', pwb_successes(gas), &
                ', S4_shared=', pwb_instrument_shared(gas), &
                ', S3=', pwb_carryforwards(gas), &
                ', fallback=', pwb_fallbacks(gas), &
                ' (maxcov/default=', pwb_fallback_maxcov(gas), &
                ', nominal/default=', pwb_fallback_nominal(gas), &
                ', other=', pwb_fallback_other(gas), ')'
            write(ulog, '(a, a, a, i0, a, i0, a, i0, a, i0, a, i0, a, i0, a, i0, a, i0, a)') &
                '  ', trim(GasLabel(gas)), &
                ': attempts=', pwb_attempts(gas), &
                ', S1/S2=', pwb_successes(gas), &
                ', S4_shared=', pwb_instrument_shared(gas), &
                ', S3=', pwb_carryforwards(gas), &
                ', fallback=', pwb_fallbacks(gas), &
                ' (maxcov/default=', pwb_fallback_maxcov(gas), &
                ', nominal/default=', pwb_fallback_nominal(gas), &
                ', other=', pwb_fallback_other(gas), ')'
        end if
    end do

    !> The guard band's report. A gas whose peak keeps landing outside the
    !> window it was given has a window problem, not a detection problem, and
    !> nothing else in the run could tell the two apart.
    do gas = firstGas, lastGas
        if (pwb_attempts(gas) > 0 .and. pwb_outside_window(gas) * 4 > pwb_attempts(gas)) then
            write(*, '(a,a,a,i0,a,i0,a)') '  WARNING: PWB peak for ', trim(GasLabel(gas)), &
                ' fell outside the configured lag window in ', pwb_outside_window(gas), &
                ' of ', pwb_attempts(gas), ' periods; review that window.'
            write(ulog, '(a,a,a,i0,a,i0,a)') '  WARNING: PWB peak for ', trim(GasLabel(gas)), &
                ' fell outside the configured lag window in ', pwb_outside_window(gas), &
                ' of ', pwb_attempts(gas), ' periods; review that window.'
        end if
    end do

    if (total_successes == 0 .and. total_instrument_shared == 0 &
        .and. total_carryforwards == 0 .and. total_fallbacks > 0) then
        write(*, '(a, i0, a, i0, a, i0, a)') '  WARNING: all PWB detections fell back: maxcov/default=', &
            total_fallback_maxcov, ', nominal/default=', total_fallback_nominal, &
            ', other=', total_fallback_other, '.'
        write(ulog, '(a, i0, a, i0, a, i0, a)') '  WARNING: all PWB detections fell back: maxcov/default=', &
            total_fallback_maxcov, ', nominal/default=', total_fallback_nominal, &
            ', other=', total_fallback_other, '.'
        call LogSay('  Review the PWB half-hourly time-lag file before interpreting method 5 as native PWB.')
    end if
end subroutine ReportPwbDiagnostics

character(32) function GasLabel(gas)
    integer, intent(in) :: gas
    character(32) :: tags(GHGNumVar)

    !> Every slot named from its record. These strings go into the PWB
    !> half-hourly time-lag file and are read back out of it, so the writer
    !> and the reader must agree. They did not: the literal cases pinned slots
    !> 5-8 to co2/h2o/ch4/gas4, so a project with water at slot 9 wrote 'h2o'
    !> for slot 9 and read it back as slot 6 - the stored lag was applied to
    !> the wrong gas, and the right one fell back to its nominal window.
    GasLabel = 'unknown'
    if (gas < firstGas .or. gas > lastGas) return
    call SpectralVarTags(tags)
    if (len_trim(tags(gas)) > 0) GasLabel = tags(gas)
end function GasLabel

!> TimelagOptGasLabel stood here: a second naming for the same slots, which
!> spelled the first four co2/h2o/ch4/4th_gas by *position* and deferred to
!> GasLabel only past the fourth. It was kept so an existing optimisation file
!> would still match, and it cost the file its meaning: on a project whose
!> records are ordered COS, CO2, H2O it headed the COS block 'co2', the CO2
!> block 'h2o', and wrote a 'ch4' block for a project that measures no methane.
!>
!> Writer and reader shared it, so a single project still round-tripped and no
!> flux moved - but the name means "slot five" while every reader takes it for
!> a species, and reordering records between two runs therefore restored a
!> cached window onto the wrong gas. That is the defect GasLabel above already
!> carries the note for; having fixed it for the table and not for the summary
!> left the two files disagreeing about what a gas is called.

end module m_pwb_timelag
