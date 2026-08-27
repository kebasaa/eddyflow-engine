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
    use m_pwb_core
    implicit none
    private
    public :: PwbDetectGas, ResetPwbDiagnostics, ReportPwbDiagnostics, InitPwbResult, GasLabel
    public :: CountPwbDiagnostic, SameAnalyser
    public :: InitPwbTimelagCache, ReadPwbTimelagCache, WritePwbTimelagCache
    public :: LookupPwbTimelagCache, StorePwbTimelagCache, SetPwbPeriodTimestamp
    public :: PostProcessPwbTimelagCache
    public :: ResetPwbAggregateSummary, AddPwbTimelagSummaryDataset, ResolvePwbAggregateSummary
    public :: RecordPwbTimelagOptPeriod, RebuildPwbTimelagOptFromCache

    integer :: pwb_attempts(E2NumVar) = 0
    integer :: pwb_successes(E2NumVar) = 0
    integer :: pwb_carryforwards(E2NumVar) = 0
    integer :: pwb_fallbacks(E2NumVar) = 0
    integer :: pwb_fallback_maxcov(E2NumVar) = 0
    integer :: pwb_fallback_nominal(E2NumVar) = 0
    integer :: pwb_fallback_other(E2NumVar) = 0
    integer :: pwb_instrument_shared(E2NumVar) = 0
    integer :: pwb_outside_window(E2NumVar) = 0
    !> Anything the classification ladder below does not name. There should
    !> never be any; it exists so that a class added later is loud rather
    !> than quietly counted as a detection, which is what happened to
    !> S4_instrument_filled.
    integer :: pwb_unclassified(E2NumVar) = 0
    logical :: pwb_bounds_warned(E2NumVar) = .false.

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
    pwb_unclassified = 0
    pwb_bounds_warned = .false.
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

!***************************************************************************
!> \brief Record what a period contributes that the settled table cannot say.
!>
!> The cache-generation counterpart of AddPwbTimelagSummaryDataset, which gated
!> membership of the aggregate dataset on the STREAMING classification. That
!> classification is a guess made having read only the periods before this one,
!> and PostProcessPwbTimelagCache overrules it once the whole run has been
!> read - so the dataset deciding the aggregate windows, and the RH classes
!> with them, was built from guesses the table had already discarded.
!>
!> Everything the table can say is left to RebuildPwbTimelagOptFromCache. Only
!> the humidity cannot: it is a property of this period's own statistics and is
!> nowhere in the cache. So it is recorded here, ungated, and the rebuild
!> decides whether it stands.
!***************************************************************************
subroutine RecordPwbTimelagOptPeriod(TimelagOpt, nrow, n)
    integer, intent(in) :: nrow, n
    type(TimeLagOptType), intent(inout) :: TimelagOpt(nrow)

    if (.not. allocated(PwbOptDate)) allocate(PwbOptDate(nrow), PwbOptTime(nrow))
    if (n < 1 .or. n > nrow) return

    !> No lags yet: the rebuild fills them from the settled table.
    TimelagOpt(n)%tlag = error
    !> Provisional and ungated - the raw humidity, not yet tested against
    !> whether water settled. Only the table knows that, so the rebuild
    !> applies the gate.
    if (Stats%RH >= 0d0 .and. Stats%RH <= 100d0) then
        TimelagOpt(n)%RH = Stats%RH
    else
        TimelagOpt(n)%RH = error
    end if
    PwbOptDate(n) = PwbPeriodDate
    PwbOptTime(n) = PwbPeriodTime
end subroutine RecordPwbTimelagOptPeriod

!***************************************************************************
!> \brief Fill the aggregate time-lag dataset from the settled table.
!>
!> A period contributes its own lag for a gas only where the FINISHED table
!> settled that gas from its own evidence - S1 or S2. Borrowed, carried,
!> interpolated and back-filled lags are excluded, exactly as the streaming
!> version excluded them. What changes is which rows those are, and that the
!> lag taken is the settled one rather than whatever happened to be applied
!> while the walk was still going.
!>
!> Rows and periods are both in period order, so this walks a cursor rather
!> than scanning: a season is some twenty thousand rows against four thousand
!> periods, and the nested form of that is quadratic for nothing.
!***************************************************************************
subroutine RebuildPwbTimelagOptFromCache(TimelagOpt, nrow, n)
    integer, intent(in) :: nrow, n
    type(TimeLagOptType), intent(inout) :: TimelagOpt(nrow)
    integer :: i, k, cursor, gas, wsl
    integer, external :: PrimaryWaterOutSlot

    if (n < 1 .or. .not. allocated(PwbOptDate)) return

    do k = 1, n
        TimelagOpt(k)%tlag = error
    end do

    cursor = 1
    do i = 1, PwbTimelagCacheN
        gas = PwbTimelagCache(i)%gas
        if (gas < firstGas .or. gas > lastGas) cycle
        if (.not. E2Col(gas)%present) cycle
        if (trim(PwbTimelagCache(i)%result%reliability_class) /= 'S1_optimal' &
            .and. trim(PwbTimelagCache(i)%result%reliability_class) /= 'S2_optimal') cycle
        if (PwbTimelagCache(i)%used_lag == error) cycle

        k = PeriodSlot(i, cursor, n)
        if (k == 0) cycle
        cursor = k
        TimelagOpt(k)%tlag(gas) = PwbTimelagCache(i)%used_lag
    end do

    !> RH travels with the water record's own time-lag, so it is gated on the
    !> site's water rather than on slot six - and on whether the TABLE settled
    !> that water, which is the part the streaming version could not know.
    wsl = PrimaryWaterOutSlot()
    do k = 1, n
        if (E2Col(wsl)%present .and. TimelagOpt(k)%tlag(wsl) /= error &
            .and. TimelagOpt(k)%RH /= error) cycle
        TimelagOpt(k)%RH = error
        TimelagOpt(k)%tlag(wsl) = error
    end do
end subroutine RebuildPwbTimelagOptFromCache

!***************************************************************************
!> \brief Which period cache row i belongs to, searching forward from cursor.
!>
!> Zero when no period claims it, which drops the row rather than guessing.
!***************************************************************************
integer function PeriodSlot(i, cursor, n)
    integer, intent(in) :: i, cursor, n
    integer :: k

    do k = cursor, n
        if (PwbOptDate(k) == PwbTimelagCache(i)%date .and. &
            PwbOptTime(k) == PwbTimelagCache(i)%time) then
            PeriodSlot = k
            return
        end if
    end do
    !> Behind the cursor, which an ordered table cannot put it - but a wrong
    !> answer here silently drops a period out of the aggregate, so it is
    !> worth the second look rather than the assumption.
    do k = 1, min(cursor - 1, n)
        if (PwbOptDate(k) == PwbTimelagCache(i)%date .and. &
            PwbOptTime(k) == PwbTimelagCache(i)%time) then
            PeriodSlot = k
            return
        end if
    end do
    PeriodSlot = 0
end function PeriodSlot

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
            !>
            !> And never across analysers. This summary had no instrument test
            !> of any kind, so on a multi-analyser site it could hand one
            !> tube's optimised window to a gas measured down another - the
            !> per-period rule two functions down refuses exactly that, and
            !> the aggregate file feeds the same windows back into a later
            !> run through SetTimelags.
            if (candidate == gas) cycle
            if (GasSlotIsWater(candidate)) cycle
            if (.not. SameAnalyser(gas, candidate)) cycle
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
        '(a,i0,a,f8.4,a,f8.4,a,f8.4,a,f8.4,a,f8.4,a,i0,a,i0,a,f8.3)') &
        'n=', PWBSetup%n_bootstrap, '_block=', PWBSetup%block_length_s, &
        '_valid=', PWBSetup%min_valid_frac, &
        '_hdi=', PWBSetup%hdi_thresh_s, '_dev=', PWBSetup%dev_thresh_s, &
        '_prefilter=', PWBSetup%hdi_prefilter_s, &
        '_smooth=', PWBSetup%smoothing_width, '_seed=', PWBSetup%random_seed, &
        '_carry=', PWBSetup%max_carry_h

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
    logical :: default_used, edge_pinned, outside, clamped, prefiltered, differenced
    real(kind = dbl) :: carry_hours, tlag_pw, corr_pw, cv_99
    real(kind = dbl) :: maxcov_lag
    integer :: ar_s, ar_w, ar_t
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
    !> Only the current version is read. Every earlier one predates a
    !> setting that is now in the fingerprint, so none of them could have
    !> matched this build in any case; refusing them by version says so
    !> plainly rather than letting the fingerprint say it obscurely.
    if (trim(line) /= 'PWB_TIMELAG_CACHE_VERSION=5') then
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
            eff_min, eff_max, eff_block, clamped, prefiltered, carry_hours, &
            tlag_pw, corr_pw, cv_99, differenced, ar_s, ar_w, ar_t, &
            maxcov_lag
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
        res%carry_hours = carry_hours
        res%tlag_pw = tlag_pw
        res%corr_pw = corr_pw
        res%cv_99 = cv_99
        res%differenced = differenced
        res%ar_order_scalar = ar_s
        res%ar_order_w = ar_w
        res%ar_order_t = ar_t
        res%maxcov_lag = maxcov_lag
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
    write(u, '(a)') 'PWB_TIMELAG_CACHE_VERSION=5'
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
        // 'hdi_prefiltered,carry_hours,' &
        // 'tlag_pw_s,corr_pw,cv_99,differenced,ar_order_scalar,ar_order_w,ar_order_t,maxcov_lag_s'
    do i = 1, PwbTimelagCacheN
        write(u, '(a,",",a,",",a,",",f12.6,",",i0,",",f12.6,",",f12.6,",",f12.6,",",a,",",l1,' &
            // '",",f12.6,",",l1,",",f14.6,",",f12.6,",",f12.6,",",i0,",",l1,' &
            // '",",a,",",a,",",a,",",a,",",a,' &
            // '",",f12.6,",",f12.6,",",f12.6,",",l1,",",l1,",",f12.6,' &
            // '",",f12.6,",",f12.8,",",f12.8,",",l1,",",i0,",",i0,",",i0,",",f12.6)') &
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
            PwbTimelagCache(i)%result%hdi_prefiltered, &
            PwbTimelagCache(i)%result%carry_hours, &
            PwbTimelagCache(i)%result%tlag_pw, PwbTimelagCache(i)%result%corr_pw, &
            PwbTimelagCache(i)%result%cv_99, PwbTimelagCache(i)%result%differenced, &
            PwbTimelagCache(i)%result%ar_order_scalar, PwbTimelagCache(i)%result%ar_order_w, &
            PwbTimelagCache(i)%result%ar_order_t, &
            PwbTimelagCache(i)%result%maxcov_lag
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
    real(kind = dbl) :: median_lag, t0, t1, span, previous, dist, limit
    real(kind = dbl) :: previous_minutes
    logical :: stale
    integer :: prev, nxt, shared, origin
    logical, external :: GasSlotIsWater

    if (PwbTimelagCacheN <= 0) return

    !> Chronological order over the whole table. The pass appends in period
    !> order already, but everything below reads a real time axis off these
    !> rows and must not depend on that being true.
    allocate(ord(PwbTimelagCacheN), tmin(PwbTimelagCacheN))
    do i = 1, PwbTimelagCacheN
        ord(i) = i
        tmin(i) = PeriodMinutes(PwbTimelagCache(i)%date, PwbTimelagCache(i)%time)
    end do
    call SortByMinutes(ord, tmin, PwbTimelagCacheN)

    !> How far a lag may travel to a period that detected none, in minutes.
    !> Zero disables the limit, which is the paper's unbounded carry.
    limit = PWBSetup%max_carry_h * 60d0

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

    !> Step 1b: every field describing HOW a period was settled starts blank,
    !> because from here on this routine owns all of them.
    !>
    !> The arms below each set what they decide, but not one of them sets
    !> everything: interpolate, carry-forward, back-fill and median all leave
    !> donor_gas alone, and median leaves carry_hours too. A field an arm does
    !> not set keeps what the STREAMING pass left in the row - a value the
    !> table has just finished overruling.
    !>
    !> That is invisible in a serial run, because the leftover is at least the
    !> same leftover every time. It is not invisible across processes: a worker
    !> starting cold leaves different ones, and this is what a parallel run
    !> disagreed with a serial run about - three interpolated rows whose
    !> donor_gas read h2o and co2 in one and co2 and none in the other, with
    !> the settled lag, the class and origin_gas all identical. Both were
    !> meaningless; they were merely differently meaningless.
    !>
    !> Blanking them here rather than patching the four arms is deliberate: an
    !> arm added later inherits the defined value instead of a stale one.
    do i = 1, PwbTimelagCacheN
        PwbTimelagCache(i)%result%donor_gas = 'none'
        PwbTimelagCache(i)%result%carry_hours = 0d0
        PwbTimelagCache(i)%result%fallback_source = 'none'
        PwbTimelagCache(i)%result%fallback_used = .false.
    end do

    allocate(idx(PwbTimelagCacheN), lag(PwbTimelagCacheN), &
        fallback_lag(PwbTimelagCacheN), settled(PwbTimelagCacheN), &
        sorted(PwbTimelagCacheN))

    !> Steps 2 and 3: S1 accepts a narrow HDI outright; S2 accepts a wider one
    !> that stays close to the last accepted lag. An S2 acceptance updates the
    !> reference, so a run of S2 periods can drift away from the S1 that
    !> anchored it - the paper's section 2.3 is ambiguous and its S3 wording
    !> implies this reading, which is dyco's too.
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
        previous_minutes = 0d0
        do j = 1, n
            i = idx(j)
            if (PwbTimelagCache(i)%result%hdi_prefiltered) cycle
            if (PwbTimelagCache(i)%result%selected_lag == error) cycle
            if (PwbTimelagCache(i)%result%edge_pinned) cycle
            if (PwbTimelagCache(i)%result%hdi_range == error) cycle
            !> An anchor too old to hand out is too old to argue from. S2
            !> accepts a vague detection for sitting close to the last
            !> accepted lag, so letting it lean on one the carry limit has
            !> already expired would walk the series along on evidence that no
            !> longer counts - and would quietly reintroduce the unbounded
            !> carry through the classifier, while the fills below expired
            !> correctly.
            stale = limit > 0d0 .and. previous /= error &
                .and. (dble(tmin(i)) - previous_minutes) > limit
            if (PwbTimelagCache(i)%result%hdi_range < PWBSetup%hdi_thresh_s) then
                lag(j) = PwbTimelagCache(i)%result%selected_lag
                settled(j) = .true.
            elseif (previous /= error .and. .not. stale) then
                if (abs(PwbTimelagCache(i)%result%selected_lag - previous) &
                    <= PWBSetup%dev_thresh_s) then
                    lag(j) = PwbTimelagCache(i)%result%selected_lag
                    settled(j) = .true.
                end if
            end if
            if (settled(j)) then
                previous = lag(j)
                previous_minutes = dble(tmin(i))
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
                PwbTimelagCache(i)%result%carry_hours = 0d0
            end if
        end do
        do j = 1, n
            if (.not. settled(j)) PwbTimelagCache(idx(j))%result%reliability_class = 'pending'
        end do

        !> Step 4: the gas's OWN lag, in its three forms, best first.
        !>
        !> All three come before any borrowing, and that is the order dyco
        !> argues for: two gases down one tube still have different delays - a
        !> systematic 0.35 s between CH4 and N2O is ordinary - so taking a
        !> neighbour's lag trades a stale number for a biased one. What
        !> decides when staleness has become the larger error is max_carry_h,
        !> which bounds every one of these three and past which the donor
        !> below takes over. This engine used to share from the analyser
        !> first, which had that argument backwards.
        !>
        !> Interpolation leads because it is the only form that follows a
        !> drifting pump rather than holding its lag flat across the drift.
        do j = 1, n
            if (settled(j)) cycle
            i = idx(j)
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
            t0 = dble(tmin(idx(prev)))
            t1 = dble(tmin(idx(nxt)))
            !> Bounded by the distance to the NEARER anchor: an interpolation
            !> spanning more than the limit in both directions is no better
            !> informed than a carry that would have expired.
            dist = min(dble(tmin(i)) - t0, t1 - dble(tmin(i)))
            if (limit > 0d0 .and. dist > limit) cycle
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
            PwbTimelagCache(i)%result%fallback_used = .false.
            PwbTimelagCache(i)%result%carry_hours = dist / 60d0
        end do

        !> Carried forward, then filled backward. Both bounded by the same
        !> limit, which they have to be: with detections either side of a long
        !> unusable stretch, an unbounded backward fill would cover from the
        !> later detection precisely the span the forward carry was just
        !> forbidden to cross, and the limit would achieve nothing.
        do j = 1, n
            i = idx(j)
            if (trim(PwbTimelagCache(i)%result%reliability_class) /= 'pending') cycle
            prev = 0
            do k = j - 1, 1, -1
                if (settled(k)) then
                    prev = k
                    exit
                end if
            end do
            if (prev == 0) cycle
            dist = dble(tmin(i)) - dble(tmin(idx(prev)))
            if (limit > 0d0 .and. dist > limit) cycle
            PwbTimelagCache(i)%used_lag = lag(prev)
            PwbTimelagCache(i)%result%reliability_class = 'S3_carryforward'
            PwbTimelagCache(i)%result%fill_method = 'carryforward'
            PwbTimelagCache(i)%result%fallback_source = 'S3_carryforward'
            PwbTimelagCache(i)%result%origin_gas = gas
            PwbTimelagCache(i)%result%fallback_used = .false.
            PwbTimelagCache(i)%result%carry_hours = dist / 60d0
        end do

        do j = 1, n
            i = idx(j)
            if (trim(PwbTimelagCache(i)%result%reliability_class) /= 'pending') cycle
            nxt = 0
            do k = j + 1, n
                if (settled(k)) then
                    nxt = k
                    exit
                end if
            end do
            if (nxt == 0) cycle
            dist = dble(tmin(idx(nxt))) - dble(tmin(i))
            if (limit > 0d0 .and. dist > limit) cycle
            PwbTimelagCache(i)%used_lag = lag(nxt)
            PwbTimelagCache(i)%result%reliability_class = 'S3_backfilled'
            PwbTimelagCache(i)%result%fill_method = 'backfilled'
            PwbTimelagCache(i)%result%fallback_source = 'backfilled'
            PwbTimelagCache(i)%result%origin_gas = gas
            PwbTimelagCache(i)%result%fallback_used = .false.
            PwbTimelagCache(i)%result%carry_hours = dist / 60d0
        end do

        !> A period the gas could not reach at all, with a limit in force.
        do j = 1, n
            i = idx(j)
            if (trim(PwbTimelagCache(i)%result%reliability_class) == 'pending' &
                .and. limit > 0d0) &
                PwbTimelagCache(i)%result%reliability_class = 'S3_expired'
        end do
    end do

    !> Step 5: a period no form of the gas's own lag could reach takes the lag
    !> of another gas measured by the same analyser in that same period. They
    !> share a tube, so they share a delay - approximately, which is why this
    !> comes after all three of the gas's own forms rather than before them.
    !>
    !> The donor is decided by analyser identity, not by the model string: two
    !> analysers of the same model are two tubes, and matching on the model
    !> made them one - the same distinction timelag_handle already draws for
    !> the water covariance. And never from water, whose lag is
    !> humidity-dependent in a way the trace gases' is not, which is exactly
    !> why ResolvePwbAggregateSummary refuses it.
    do i = 1, PwbTimelagCacheN
        gas = PwbTimelagCache(i)%gas
        if (gas < firstGas .or. gas > lastGas) cycle
        if (trim(PwbTimelagCache(i)%result%reliability_class) /= 'pending' &
            .and. trim(PwbTimelagCache(i)%result%reliability_class) /= 'S3_expired') cycle
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
        PwbTimelagCache(i)%result%carry_hours = 0d0
    end do

    !> Whatever the streaming pass settled on stands as the last resort of
    !> step 8 - covariance maximisation where it fell back, the detection
    !> itself otherwise. Captured here, by cache row, because three passes now
    !> run between this and the arm that reads it, and any of them may
    !> overwrite used_lag.
    do i = 1, PwbTimelagCacheN
        fallback_lag(i) = PwbTimelagCache(i)%used_lag
    end do

    !> Step 6, per gas: the median of what the gas itself detected.
    do gas = firstGas, lastGas
        n = 0
        do k = 1, PwbTimelagCacheN
            i = ord(k)
            if (PwbTimelagCache(i)%gas /= gas) cycle
            n = n + 1
            idx(n) = i
        end do
        if (n == 0) cycle

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
        if (nsel > 0) median_lag = MedianOf(sorted, nsel)

        if (median_lag == error) cycle
        do j = 1, n
            i = idx(j)
            if (trim(PwbTimelagCache(i)%result%reliability_class) /= 'pending' &
                .and. trim(PwbTimelagCache(i)%result%reliability_class) /= 'S3_expired') cycle
            !> Every one of the lags being averaged here is one the rule has
            !> just rejected, so this is a last resort and is labelled as one.
            PwbTimelagCache(i)%used_lag = median_lag
            PwbTimelagCache(i)%result%reliability_class = 'S3_median'
            PwbTimelagCache(i)%result%fill_method = 'median'
            PwbTimelagCache(i)%result%fallback_source = 'median'
            PwbTimelagCache(i)%result%origin_gas = gas
            PwbTimelagCache(i)%result%fallback_used = .false.
        end do
    end do

    !> Step 7: the tube-mate again, this time from a donor that was itself
    !> filled rather than detected.
    !>
    !> Step 5 above only borrows from a donor the rule trusted outright, which
    !> is the right first answer. But a gas whose every detection was rejected
    !> - carbonyl sulfide, whose HDI routinely spans the whole search window -
    !> reaches step 6 with nothing to take a median of, and then had only its
    !> own rejected covariance maximisation left. That put COS at 10.6 s on a
    !> tube whose delay is 16.2 s, while the CO2 beside it in the same tube
    !> carried a perfectly good interpolated 16.5 s. Preferring "my own
    !> rejected number" to "my tube-mate's filled one" is the wrong way round.
    !>
    !> Same donor rule as step 5 - one analyser, never water, same period -
    !> and a distinct label, because a donor that was itself filled is weaker
    !> evidence than one the rule trusted. Donors filled by THIS pass are
    !> excluded, so a borrowed lag cannot be borrowed onward.
    do i = 1, PwbTimelagCacheN
        gas = PwbTimelagCache(i)%gas
        if (gas < firstGas .or. gas > lastGas) cycle
        if (trim(PwbTimelagCache(i)%result%reliability_class) /= 'pending' &
            .and. trim(PwbTimelagCache(i)%result%reliability_class) /= 'S3_expired') cycle
        shared = 0
        do j = 1, PwbTimelagCacheN
            if (j == i) cycle
            if (PwbTimelagCache(j)%date /= PwbTimelagCache(i)%date) cycle
            if (PwbTimelagCache(j)%time /= PwbTimelagCache(i)%time) cycle
            g = PwbTimelagCache(j)%gas
            if (g < firstGas .or. g > lastGas) cycle
            if (GasSlotIsWater(g)) cycle
            if (.not. SameAnalyser(gas, g)) cycle
            if (trim(PwbTimelagCache(j)%result%reliability_class) == 'pending') cycle
            if (trim(PwbTimelagCache(j)%result%reliability_class) == 'S3_expired') cycle
            if (trim(PwbTimelagCache(j)%result%reliability_class) == 'fallback') cycle
            if (trim(PwbTimelagCache(j)%result%reliability_class) == 'S4_instrument_filled') cycle
            shared = j
            exit
        end do
        if (shared == 0) cycle
        PwbTimelagCache(i)%used_lag = PwbTimelagCache(shared)%used_lag
        PwbTimelagCache(i)%result%reliability_class = 'S4_instrument_filled'
        PwbTimelagCache(i)%result%fill_method = 'instrument_filled'
        PwbTimelagCache(i)%result%fallback_source = 'instrument_filled'
        PwbTimelagCache(i)%result%fallback_used = .false.
        PwbTimelagCache(i)%result%donor_gas = GasLabel(PwbTimelagCache(shared)%gas)
        PwbTimelagCache(i)%result%origin_gas = PwbTimelagCache(shared)%gas
        PwbTimelagCache(i)%result%carry_hours = 0d0
    end do

    !> Step 8: nothing reached this period, so it falls back to its own
    !> covariance maximum - which is what the label has always claimed.
    !>
    !> It used to hand back whatever the STREAMING pass had settled on. For a
    !> gas the rule rejected everywhere those are the same thing, the streaming
    !> pass having had no previous lag to carry either. They part company when
    !> the HDI pre-filter - which runs here and not there - discards every
    !> detection a gas had: the streaming pass had settled and carried, so what
    !> came back was a carried lag wearing the maxcov_default label, and which
    !> lag it was depended on where the pass began.
    do i = 1, PwbTimelagCacheN
        if (PwbTimelagCache(i)%gas < firstGas .or. PwbTimelagCache(i)%gas > lastGas) cycle
        if (trim(PwbTimelagCache(i)%result%reliability_class) /= 'pending' &
            .and. trim(PwbTimelagCache(i)%result%reliability_class) /= 'S3_expired') cycle
        !> maxcov_lag is this period's own covariance maximum, taken at
        !> detection time. fallback_lag - whatever the streaming pass had
        !> settled on - is kept only for a row carrying no maxcov of its
        !> own, which a table written by this build cannot produce.
        if (PwbTimelagCache(i)%result%maxcov_lag /= error) then
            PwbTimelagCache(i)%used_lag = PwbTimelagCache(i)%result%maxcov_lag
        else
            PwbTimelagCache(i)%used_lag = fallback_lag(i)
        end if
        PwbTimelagCache(i)%result%reliability_class = 'fallback'
        PwbTimelagCache(i)%result%fill_method = 'maxcov_default'
        PwbTimelagCache(i)%result%fallback_source = 'maxcov_default'
        PwbTimelagCache(i)%result%fallback_used = .true.
    end do

    !> fallback_used means "no evidence reached this period", and only the
    !> terminal arm above leaves it set. It is raised at detection time for a
    !> period that produced no usable lag of its own, so every arm that later
    !> finds one - interpolated, carried, filled backward, shared, median -
    !> has to lower it again, or the run summary reports a period that WAS
    !> settled from evidence as a fallback.
    !>
    !> Row lag, applied lag and the "default used" flag follow the settled
    !> time-lag, so the production pass reads one consistent set of numbers.
    !>
    !> applied_lag is the lag measured back off the record shift, not the one
    !> that was asked for. Interpolation makes fractional lags ordinary and
    !> the data can only move by whole records, so the two now differ by up to
    !> half a sample as a matter of course - and it is the shifted one that
    !> reached the fluxes.
    do i = 1, PwbTimelagCacheN
        if (PwbTimelagCache(i)%used_lag == error) cycle
        PwbTimelagCache(i)%row_lag = nint(PwbTimelagCache(i)%used_lag * Metadata%ac_freq)
        PwbTimelagCache(i)%actual_lag = PwbTimelagCache(i)%result%selected_lag
        if (PwbTimelagCache(i)%actual_lag == error) &
            PwbTimelagCache(i)%actual_lag = PwbTimelagCache(i)%used_lag
        PwbTimelagCache(i)%default_used = &
            trim(PwbTimelagCache(i)%result%fill_method) == 'maxcov_default'
        PwbTimelagCache(i)%result%applied_lag = &
            dble(PwbTimelagCache(i)%row_lag) / Metadata%ac_freq
        PwbTimelagCache(i)%result%applied_row_lag = PwbTimelagCache(i)%row_lag
    end do
    PwbCacheDirty = .true.

    !> The tallies the run log prints describe the settled table, not the
    !> streaming guesses that produced it.
    call ResetPwbDiagnostics()
    do i = 1, PwbTimelagCacheN
        call CountPwbDiagnostic(PwbTimelagCache(i)%gas, PwbTimelagCache(i)%result)
    end do

    !> And so does the donor tally the aggregate summary picks each gas's
    !> lender from. It was accumulated by AddPwbTimelagSummaryDataset as the
    !> pre-pass went, off the STREAMING classification - so it carried the
    !> streaming classifier's state into the summary, and
    !> ResolvePwbAggregateSummary reads it to choose a donor by count.
    !>
    !> Two things were wrong with that. The counts described guesses the
    !> table above has since overruled, so the pre-pass summary and the one
    !> production writes at the end of the same run could disagree - the
    !> latter is already settled-table-derived, because a cache hit fills
    !> PWBResult from the cache. And the streaming chain depends on where a
    !> pass started, which is what stops the pre-pass being split across
    !> processes: a worker beginning cold classifies its first periods
    !> differently and tallies a different lender.
    !>
    !> Counted here instead, off the same rows the tallies above use. The
    !> arm is the mirror of the streaming one: a row the settled table did
    !> not settle from its own evidence borrowed, and origin_gas names who
    !> from.
    PwbSummaryDonorCount = 0
    do i = 1, PwbTimelagCacheN
        gas = PwbTimelagCache(i)%gas
        if (gas < firstGas .or. gas > lastGas) cycle
        if (.not. E2Col(gas)%present) cycle
        if (trim(PwbTimelagCache(i)%result%reliability_class) == 'S1_optimal' .or. &
            trim(PwbTimelagCache(i)%result%reliability_class) == 'S2_optimal') cycle
        origin = PwbTimelagCache(i)%result%origin_gas
        if (origin >= firstGas .and. origin <= lastGas .and. origin /= gas) &
            PwbSummaryDonorCount(gas, origin) = &
                PwbSummaryDonorCount(gas, origin) + 1
    end do

    deallocate(ord, tmin, idx, lag, fallback_lag, settled, sorted)
end subroutine PostProcessPwbTimelagCache

!***************************************************************************
!> Two gas slots measured by the same physical analyser - proven, not assumed.
!>
!> A time lag is a property of one tube. Borrowing one across analysers is
!> not an approximation, it is a different measurement, so this answers false
!> unless both records NAME an instrument and the names are the same.
!>
!> In particular there is no fall-back to the model string. Two LI-7200s at
!> one site share a model and nothing else; matching on it handed one
!> analyser's lag to the other's gases. And a record that names no instrument
!> cannot be shown to share anything, so it neither donates nor receives -
!> absence of evidence is not identity. A project that wants sharing states
!> instr_<K>_name, and ReportPwbDiagnostics says so when it does not.
!>
!> The model is compared as well as the name: equal names with unequal models
!> is a malformed project, not a shared analyser.
!***************************************************************************
logical function SameAnalyser(a, b)
    integer, intent(in) :: a, b

    SameAnalyser = .false.
    if (a < 1 .or. b < 1) return
    if (len_trim(E2Col(a)%instr_name) == 0) return
    if (len_trim(E2Col(b)%instr_name) == 0) return
    if (E2Col(a)%instr_name /= E2Col(b)%instr_name) return
    if (E2Col(a)%instr%model /= E2Col(b)%instr%model) return
    SameAnalyser = .true.
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

subroutine PwbDetectGas(Set, nrow, ncol, gas, LocResult, success)
    use m_index_parameters
    use m_log
    implicit none
    integer, intent(in) :: nrow, ncol, gas
    real(kind = dbl), intent(in) :: Set(nrow, ncol)
    type(PWBResultType), intent(out) :: LocResult
    logical, intent(out) :: success

    integer :: min_rl, max_rl, lag, margin, trail, eval_lo, eval_hi
    integer :: nvalid_w, nvalid_t, nvalid_s
    real(kind = dbl) :: min_valid
    real(kind = dbl), allocatable :: ww(:), tt(:), ss(:)
    real(kind = dbl), allocatable :: s_fs(:), w_fs(:), t_fs(:)
    real(kind = dbl), allocatable :: s_fw(:), w_fw(:)
    real(kind = dbl), allocatable :: s_ft(:), t_ft(:)
    real(kind = dbl), allocatable :: raw_ccov(:)
    type(PwbPreWhitenType) :: pw
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
    !> The first term is what the centred rolling mean consumes, so every lag
    !> inside the search range keeps a genuine smoothed value rather than one
    !> carried in from the edge. SmoothAndFill's window reaches `lead` below a
    !> position and `trail` above it, and trail = width/2 is the larger of the
    !> two whatever the parity, so covering trail covers both.
    !>
    !> The two seconds on top of that are there because the window is a guess:
    !> a peak just outside it used to be indistinguishable from an ordinary
    !> edge-pinned failure, and is now reported. The applied lag is still taken
    !> from inside the declared window - the guard band informs, it does not
    !> widen the search.
    !>
    !> One sample of headroom below nrow, because the differencing branch
    !> leaves n_eff = nrow - 1 and the evaluated range has to stay inside it.
    trail = max(1, PWBSetup%smoothing_width) / 2
    margin = max(trail, nint(2d0 * Metadata%ac_freq))
    eval_lo = max(min_rl - margin, -(nrow - 3))
    eval_hi = min(max_rl + margin, nrow - 3)
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

    call EnsurePwbScratch(nrow, eval_lo, eval_hi, max(1, PWBSetup%n_bootstrap))

    allocate(s_fs(nrow), w_fs(nrow), t_fs(nrow))
    allocate(s_fw(nrow), w_fw(nrow), s_ft(nrow), t_ft(nrow))
    allocate(raw_ccov(min_rl:max_rl))

    !> Stationarity, the AR fits, the filtered series, the raw
    !> cross-covariance and the full-data pre-whitened CCF, in one place with
    !> no engine state - which is what lets the reference test drive exactly
    !> this arithmetic against RFlux's frozen output.
    call PwbPreWhiten(ss, ww, tt, nrow, min_rl, max_rl, error, pw, &
        s_fs, w_fs, t_fs, s_fw, w_fw, s_ft, t_ft, raw_ccov, sc_xc, sc_yc)

    LocResult%differenced = pw%differenced
    LocResult%ar_order_scalar = pw%p_scalar
    LocResult%ar_order_w = pw%p_w
    LocResult%ar_order_t = pw%p_t
    LocResult%tlag_pw = dble(pw%tlag_pw_rl) / Metadata%ac_freq
    LocResult%corr_pw = pw%corr_pw
    LocResult%cv_99 = BartlettCv99(pw%n_eff)

    !> The bootstrap runs over n_eff, not nrow: the differencing branch
    !> returns one sample fewer, and resampling the stale tail would feed it a
    !> value no difference produced.
    combo = (/'cw', 'wc', 'ct', 'tc'/)
    call RunPwbCombination(w_fs, s_fs, pw%n_eff, min_rl, max_rl, eval_lo, eval_hi, &
        gas, combo(1), candidate(1), ok(1))
    call RunPwbCombination(w_fw, s_fw, pw%n_eff, min_rl, max_rl, eval_lo, eval_hi, &
        gas, combo(2), candidate(2), ok(2))
    call RunPwbCombination(t_fs, s_fs, pw%n_eff, min_rl, max_rl, eval_lo, eval_hi, &
        gas, combo(3), candidate(3), ok(3))
    call RunPwbCombination(t_ft, s_ft, pw%n_eff, min_rl, max_rl, eval_lo, eval_hi, &
        gas, combo(4), candidate(4), ok(4))

    call SelectBestCandidate(candidate, ok, LocResult, success)
    !> SelectBestCandidate copies a candidate wholesale, so the per-period
    !> diagnostics above have to be restored onto the winner.
    LocResult%differenced = pw%differenced
    LocResult%ar_order_scalar = pw%p_scalar
    LocResult%ar_order_w = pw%p_w
    LocResult%ar_order_t = pw%p_t
    LocResult%tlag_pw = dble(pw%tlag_pw_rl) / Metadata%ac_freq
    LocResult%corr_pw = pw%corr_pw
    LocResult%cv_99 = BartlettCv99(pw%n_eff)
    if (LocResult%peak_outside_window) pwb_outside_window(gas) = pwb_outside_window(gas) + 1
    lag = LocResult%row_lag
    if (success .and. lag >= min_rl .and. lag <= max_rl) LocResult%raw_covariance = raw_ccov(lag)

    deallocate(ww, tt, ss, raw_ccov)
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
    res%carry_hours = 0d0
    res%maxcov_lag = error
    res%tlag_pw = error
    res%corr_pw = error
    res%cv_99 = error
    res%differenced = .false.
    res%ar_order_scalar = 0
    res%ar_order_w = 0
    res%ar_order_t = 0
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

    !> The block is derived per gas, not taken as given.
    !>
    !> A resampling block shorter than the lag range cannot contain the lag
    !> structure the bootstrap exists to preserve, so R couples the two as
    !> l = 2*LAG.MAX and dyco floors that coupling at the configured value:
    !> block = max(configured, 2*widest bound). This used to warn and then use
    !> the short block anyway, which on a project with a 25 s window and the
    !> default 20 s block meant every gas resampled under-blocked while the
    !> log said so and nothing acted on it.
    block_len = max(requested_block_len, 2 * widest)
    res%block_length_clamped = block_len > requested_block_len
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
        call ComputeCcfWindow(sc_xb, sc_yb, n, eval_lo, eval_hi, sc_ccf, sc_xc, sc_yc)
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
    elseif (trim(res%reliability_class) == 'S4_instrument_shared' .or. &
        trim(res%reliability_class) == 'S4_instrument_filled') then
        !> Both arms are the same answer to a reader: this gas did not
        !> detect a lag, and took one measured down the same tube. They
        !> differ in when the settled table reached for the neighbour,
        !> which is a distinction for the half-hourly file and not for a
        !> summary line - the S3 arms are lumped here for the same reason.
        pwb_instrument_shared(gas) = pwb_instrument_shared(gas) + 1
    elseif (index(res%reliability_class, 'S3_') == 1) then
        !> Every gap-filled arm - interpolated, back-filled, carried forward,
        !> median - counts here, which is what "not detected in this period"
        !> means to a reader of the summary.
        pwb_carryforwards(gas) = pwb_carryforwards(gas) + 1
    elseif (trim(res%reliability_class) == 'S1_optimal' .or. &
        trim(res%reliability_class) == 'S2_optimal') then
        !> Named, not inferred from what is left over. As a catch-all this
        !> branch reported S4_instrument_filled as a detection - six of
        !> them for cos on base_pwb_cache, where the settled table holds
        !> none at all - because that class reaches here and nothing above
        !> claimed it. A column headed S1/S2 has to mean S1 or S2.
        pwb_successes(gas) = pwb_successes(gas) + 1
    else
        pwb_unclassified(gas) = pwb_unclassified(gas) + 1
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
                ', S4_borrowed=', pwb_instrument_shared(gas), &
                ', S3=', pwb_carryforwards(gas), &
                ', fallback=', pwb_fallbacks(gas), &
                ' (maxcov/default=', pwb_fallback_maxcov(gas), &
                ', nominal/default=', pwb_fallback_nominal(gas), &
                ', other=', pwb_fallback_other(gas), ')'
            write(ulog, '(a, a, a, i0, a, i0, a, i0, a, i0, a, i0, a, i0, a, i0, a, i0, a)') &
                '  ', trim(GasLabel(gas)), &
                ': attempts=', pwb_attempts(gas), &
                ', S1/S2=', pwb_successes(gas), &
                ', S4_borrowed=', pwb_instrument_shared(gas), &
                ', S3=', pwb_carryforwards(gas), &
                ', fallback=', pwb_fallbacks(gas), &
                ' (maxcov/default=', pwb_fallback_maxcov(gas), &
                ', nominal/default=', pwb_fallback_nominal(gas), &
                ', other=', pwb_fallback_other(gas), ')'
        end if
    end do

    !> A class the ladder does not name would otherwise be invisible: it
    !> would simply be missing from a line whose columns are meant to add
    !> up to attempts. Nothing produces one today.
    do gas = firstGas, lastGas
        if (pwb_unclassified(gas) > 0) then
            write(*, '(a,i0,a,a)') '  NOTE: ', pwb_unclassified(gas), &
                ' period(s) of ' // trim(GasLabel(gas)) // &
                ' carry a reliability class this summary does not know.'
            write(ulog, '(a,i0,a,a)') '  NOTE: ', pwb_unclassified(gas), &
                ' period(s) of ' // trim(GasLabel(gas)) // &
                ' carry a reliability class this summary does not know.'
        end if
    end do

    !> A gas whose record does not name its instrument can neither donate a
    !> time lag nor receive one, because nothing proves it shares a tube with
    !> anything. That is deliberate - a lag belongs to one tube, and the model
    !> string cannot tell two analysers of the same model apart - but a user
    !> whose gas fell through to the median instead of borrowing deserves to
    !> know it was the missing name that decided it.
    do gas = firstGas, lastGas
        if (pwb_attempts(gas) > 0 .and. len_trim(E2Col(gas)%instr_name) == 0) then
            write(*, '(a,a,a)') '  NOTE: ', trim(GasLabel(gas)), &
                ' names no instrument, so it neither donates a time lag to' &
                // ' nor takes one from another gas.'
            write(ulog, '(a,a,a)') '  NOTE: ', trim(GasLabel(gas)), &
                ' names no instrument, so it neither donates a time lag to' &
                // ' nor takes one from another gas.'
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
