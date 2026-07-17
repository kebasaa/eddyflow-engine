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
    public :: PwbDetectGas, ResetPwbDiagnostics, ReportPwbDiagnostics, InitPwbResult, WritePwbDiagnostic, GasLabel
    public :: InitPwbTimelagCache, ReadPwbTimelagCache, WritePwbTimelagCache
    public :: LookupPwbTimelagCache, StorePwbTimelagCache

    logical :: pwb_diag_header_written = .false.
    integer :: pwb_attempts(E2NumVar) = 0
    integer :: pwb_successes(E2NumVar) = 0
    integer :: pwb_carryforwards(E2NumVar) = 0
    integer :: pwb_fallbacks(E2NumVar) = 0
    integer :: pwb_fallback_maxcov(E2NumVar) = 0
    integer :: pwb_fallback_nominal(E2NumVar) = 0
    integer :: pwb_fallback_other(E2NumVar) = 0
    integer :: pwb_instrument_shared(E2NumVar) = 0
    logical :: pwb_bounds_warned(E2NumVar) = .false.
    logical :: pwb_block_warned(E2NumVar) = .false.

contains

subroutine ResetPwbDiagnostics()
    pwb_diag_header_written = .false.
    pwb_attempts = 0
    pwb_successes = 0
    pwb_carryforwards = 0
    pwb_fallbacks = 0
    pwb_fallback_maxcov = 0
    pwb_fallback_nominal = 0
    pwb_fallback_other = 0
    pwb_instrument_shared = 0
    pwb_bounds_warned = .false.
    pwb_block_warned = .false.
end subroutine ResetPwbDiagnostics

subroutine InitPwbTimelagCache()
    if (allocated(PwbTimelagCache)) deallocate(PwbTimelagCache)
    PwbTimelagCacheN = 0
    PwbCacheLoaded = .false.
    PwbCacheDirty = .false.
end subroutine InitPwbTimelagCache

subroutine StorePwbTimelagCache(gas, stage, actual_lag, used_lag, row_lag, default_used, res)
    integer, intent(in) :: gas, row_lag
    character(*), intent(in) :: stage
    real(kind = dbl), intent(in) :: actual_lag, used_lag
    logical, intent(in) :: default_used
    type(PWBResultType), intent(in) :: res
    type(PWBTimelagCacheEntryType), allocatable :: tmp(:)
    integer :: i

    do i = 1, PwbTimelagCacheN
        if (PwbTimelagCache(i)%date == Stats%date .and. PwbTimelagCache(i)%time == Stats%time &
            .and. PwbTimelagCache(i)%gas == gas .and. trim(PwbTimelagCache(i)%stage) == trim(stage)) then
            PwbTimelagCache(i)%actual_lag = actual_lag
            PwbTimelagCache(i)%used_lag = used_lag
            PwbTimelagCache(i)%row_lag = row_lag
            PwbTimelagCache(i)%default_used = default_used
            PwbTimelagCache(i)%result = res
            PwbCacheDirty = .true.
            return
        end if
    end do

    allocate(tmp(PwbTimelagCacheN + 1))
    if (PwbTimelagCacheN > 0) tmp(1:PwbTimelagCacheN) = PwbTimelagCache(1:PwbTimelagCacheN)
    tmp(PwbTimelagCacheN + 1)%date = Stats%date
    tmp(PwbTimelagCacheN + 1)%time = Stats%time
    tmp(PwbTimelagCacheN + 1)%gas = gas
    tmp(PwbTimelagCacheN + 1)%stage = stage
    tmp(PwbTimelagCacheN + 1)%actual_lag = actual_lag
    tmp(PwbTimelagCacheN + 1)%used_lag = used_lag
    tmp(PwbTimelagCacheN + 1)%row_lag = row_lag
    tmp(PwbTimelagCacheN + 1)%default_used = default_used
    tmp(PwbTimelagCacheN + 1)%result = res
    call move_alloc(tmp, PwbTimelagCache)
    PwbTimelagCacheN = PwbTimelagCacheN + 1
    PwbCacheDirty = .true.
end subroutine StorePwbTimelagCache

subroutine LookupPwbTimelagCache(gas, stage, found, actual_lag, used_lag, row_lag, default_used, res)
    integer, intent(in) :: gas
    character(*), intent(in) :: stage
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
    do i = 1, PwbTimelagCacheN
        if (PwbTimelagCache(i)%date == Stats%date .and. PwbTimelagCache(i)%time == Stats%time &
            .and. PwbTimelagCache(i)%gas == gas .and. trim(PwbTimelagCache(i)%stage) == trim(stage)) then
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

character(256) function PwbCacheFingerprint()
    write(PwbCacheFingerprint, '(a,4(l1,":"),a,8(f10.4,":"),a,i0,a,f8.4,a,f8.4,a,f8.4,a,f8.4,a,f8.4,a,i0,a,i0,a,l1,a,i0,a,l1)') &
        'provided=', PWBSetup%lag_bounds_provided(co2), PWBSetup%lag_bounds_provided(h2o), &
        PWBSetup%lag_bounds_provided(ch4), PWBSetup%lag_bounds_provided(gas4), 'bounds=', &
        merge(PWBSetup%min_lag(co2), 0d0, PWBSetup%lag_bounds_provided(co2)), &
        merge(PWBSetup%max_lag(co2), 0d0, PWBSetup%lag_bounds_provided(co2)), &
        merge(PWBSetup%min_lag(h2o), 0d0, PWBSetup%lag_bounds_provided(h2o)), &
        merge(PWBSetup%max_lag(h2o), 0d0, PWBSetup%lag_bounds_provided(h2o)), &
        merge(PWBSetup%min_lag(ch4), 0d0, PWBSetup%lag_bounds_provided(ch4)), &
        merge(PWBSetup%max_lag(ch4), 0d0, PWBSetup%lag_bounds_provided(ch4)), &
        merge(PWBSetup%min_lag(gas4), 0d0, PWBSetup%lag_bounds_provided(gas4)), &
        merge(PWBSetup%max_lag(gas4), 0d0, PWBSetup%lag_bounds_provided(gas4)), &
        'n=', PWBSetup%n_bootstrap, '_block=', PWBSetup%block_length_s, '_valid=', PWBSetup%min_valid_frac, &
        '_hdi=', PWBSetup%hdi_thresh_s, '_dev=', PWBSetup%dev_thresh_s, '_prefilter=', PWBSetup%hdi_prefilter_s, &
        '_smooth=', PWBSetup%smoothing_width, '_seed=', PWBSetup%random_seed, '_approx=', PWBSetup%approx_ccf, &
        '_ar=', PWBSetup%max_ar_order, '_pre=', PWBSetup%detect_prewpl
end function PwbCacheFingerprint

subroutine ReadPwbTimelagCache(path, recognized, valid)
    character(*), intent(in) :: path
    logical, intent(out) :: recognized, valid
    integer :: u, ios, gas, row_lag, period_seconds
    character(1024) :: line
    character(256) :: fingerprint
    character(10) :: date
    character(5) :: time
    character(8) :: stage
    character(24) :: reliability, fallback
    character(8) :: donor
    real(kind = dbl) :: actual_lag, used_lag, selected_lag, hdi_low, hdi_high, hdi_range
    logical :: default_used
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
    if (index(line, 'PWB_TIMELAG_CACHE_VERSION=') == 1) then
        recognized = .true.
    elseif (trim(line) /= 'PWB_TIMELAG_CACHE_VERSION=1') then
        close(u)
        return
    end if
    if (trim(line) /= 'PWB_TIMELAG_CACHE_VERSION=1') then
        close(u)
        return
    end if
    read(u, '(a)', iostat=ios) line
    if (ios /= 0 .or. index(line, 'fingerprint=') /= 1) then
        close(u)
        return
    end if
    fingerprint = line(13:len_trim(line))
    if (trim(fingerprint) /= trim(PwbCacheFingerprint())) then
        write(*, '(a)') ' Fatal> PWB time-lag cache is incompatible with the selected PWB settings.'
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
        write(*, '(a)') ' Fatal> PWB time-lag cache averaging-period duration is incompatible with this project.'
        close(u)
        return
    end if
    read(u, '(a)', iostat=ios) line
    if (ios /= 0 .or. trim(line) /= 'data') then
        close(u)
        return
    end if
    read(u, '(a)', iostat=ios) line
    if (ios /= 0 .or. index(line, 'date,time,gas,stage,') /= 1) then
        close(u)
        return
    end if
    do
        read(u, '(a)', iostat=ios) line
        if (ios /= 0) exit
        call InitPwbResult(res)
        read(line, *, iostat=ios) date, time, gas, stage, actual_lag, used_lag, row_lag, default_used, &
            reliability, fallback, donor, selected_lag, hdi_low, hdi_high, hdi_range
        if (ios /= 0) then
            close(u)
            return
        end if
        res%reliability_class = reliability
        res%fallback_source = fallback
        res%donor_gas = donor
        res%selected_lag = selected_lag
        res%hdi_low = hdi_low
        res%hdi_high = hdi_high
        res%hdi_range = hdi_range
        res%applied_lag = used_lag
        res%applied_row_lag = row_lag
        res%fallback_used = default_used .or. trim(reliability) == 'fallback'
        call StorePwbTimelagCacheAt(date, time, gas, stage, actual_lag, used_lag, row_lag, default_used, res)
    end do
    close(u)
    PwbCacheLoaded = .true.
    PwbCacheDirty = .false.
    valid = PwbTimelagCacheN > 0
end subroutine ReadPwbTimelagCache

subroutine StorePwbTimelagCacheAt(date, time, gas, stage, actual_lag, used_lag, row_lag, default_used, res)
    character(*), intent(in) :: date, time, stage
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
    tmp(PwbTimelagCacheN + 1)%stage = stage
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
        // '_pwb_timelag_cache' // Timestamp_FilePadding // CsvExt
    open(newunit=u, file=path, status='replace', iostat=ios, encoding='utf-8')
    if (ios /= 0) return
    write(u, '(a)') 'PWB_TIMELAG_CACHE_VERSION=1'
    write(u, '(a)') 'fingerprint=' // trim(PwbCacheFingerprint())
    write(u, '(a)') 'project_id=' // trim(EddyFlowProj%id)
    write(u, '(a,i0)') 'period_seconds=', RPsetup%avrg_len
    write(u, '(a)') 'data'
    write(u, '(a)') 'date,time,gas,stage,actual_lag_s,used_lag_s,row_lag,default_used,reliability_class,' &
        // 'fallback_source,donor_gas,selected_lag_s,hdi_low_s,hdi_high_s,hdi_range_s'
    do i = 1, PwbTimelagCacheN
        write(u, '(a,",",a,",",i0,",",a,",",f12.6,",",f12.6,",",i0,",",l1,",",a,",",a,' &
            // '",",a,",",f12.6,",",f12.6,",",f12.6,",",f12.6)') &
            trim(PwbTimelagCache(i)%date), trim(PwbTimelagCache(i)%time), PwbTimelagCache(i)%gas, &
            trim(PwbTimelagCache(i)%stage), PwbTimelagCache(i)%actual_lag, PwbTimelagCache(i)%used_lag, &
            PwbTimelagCache(i)%row_lag, PwbTimelagCache(i)%default_used, &
            trim(PwbTimelagCache(i)%result%reliability_class), trim(PwbTimelagCache(i)%result%fallback_source), &
            trim(PwbTimelagCache(i)%result%donor_gas), PwbTimelagCache(i)%result%selected_lag, &
            PwbTimelagCache(i)%result%hdi_low, PwbTimelagCache(i)%result%hdi_high, PwbTimelagCache(i)%result%hdi_range
    end do
    close(u)
    PwbTimelagCache_Path = path
    PwbCacheLoaded = .true.
    PwbCacheDirty = .false.
    write(*, '(a)') ' PWB per-period time-lag cache written to: ' // trim(path)
end subroutine WritePwbTimelagCache

subroutine PwbDetectGas(Set, nrow, ncol, gas, LocResult, success)
    implicit none
    integer, intent(in) :: nrow, ncol, gas
    real(kind = dbl), intent(in) :: Set(nrow, ncol)
    type(PWBResultType), intent(out) :: LocResult
    logical, intent(out) :: success

    integer :: min_rl, max_rl, lag
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
    if (gas < co2 .or. gas > gas4) return
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
        pwb_bounds_warned(gas) = .true.
    end if
    if (min_rl >= max_rl) then
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

    combo = (/'cw', 'wc', 'ct', 'tc'/)
    call RunPwbCombination(w_fs, s_fs, nrow, min_rl, max_rl, gas, combo(1), candidate(1), ok(1))
    call RunPwbCombination(w_fw, s_fw, nrow, min_rl, max_rl, gas, combo(2), candidate(2), ok(2))
    call RunPwbCombination(t_fs, s_fs, nrow, min_rl, max_rl, gas, combo(3), candidate(3), ok(3))
    call RunPwbCombination(t_ft, s_ft, nrow, min_rl, max_rl, gas, combo(4), candidate(4), ok(4))

    call SelectBestCandidate(candidate, ok, LocResult, success)
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
    res%donor_gas = ''
    res%edge_pinned = .false.
    res%fallback_used = .false.
    res%block_length_clamped = .false.
    res%effective_min_lag = error
    res%effective_max_lag = error
    res%effective_block_length_s = error
    res%raw_covariance = error
    res%ccf_at_mode = 0d0
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

subroutine FitArAic(x, n, phi_best, p_best)
    integer, intent(in) :: n
    real(kind = dbl), intent(in) :: x(n)
    real(kind = dbl), allocatable, intent(out) :: phi_best(:)
    integer, intent(out) :: p_best
    integer :: max_lag, p, i
    real(kind = dbl) :: meanx, sigma2, kappa, best_aic, aic
    real(kind = dbl), allocatable :: acf(:), phi(:), phi_old(:)

    max_lag = min(int(floor(100d0 * log10(dble(max(n, 2))))), n - 1)
    if (PWBSetup%max_ar_order > 0) max_lag = min(max_lag, PWBSetup%max_ar_order)
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

subroutine RunPwbCombination(x, y, n, min_rl, max_rl, gas, combo, res, ok)
    integer, intent(in) :: n, min_rl, max_rl, gas
    real(kind = dbl), intent(in) :: x(n), y(n)
    character(2), intent(in) :: combo
    type(PWBResultType), intent(out) :: res
    logical, intent(out) :: ok
    integer :: b, i, pos, block_len, nblocks, start, state
    integer :: requested_block_len
    integer :: nboot, lag, best_idx, full_min_rl, full_max_rl
    integer, allocatable :: boot_lags(:), counts(:)
    real(kind = dbl), allocatable :: xb(:), yb(:), ccf(:), smooth(:), hdi_samples(:)
    real(kind = dbl), allocatable :: mean_smooth(:), mean_ccf(:)

    call InitPwbResult(res)
    res%best_combination = combo
    ok = .false.
    nboot = max(1, PWBSetup%n_bootstrap)
    full_max_rl = max(abs(min_rl), abs(max_rl))
    full_min_rl = -full_max_rl
    requested_block_len = nint(PWBSetup%block_length_s * Metadata%ac_freq)
    if (requested_block_len <= 0) requested_block_len = max(1, 2 * full_max_rl)
    block_len = requested_block_len
    res%block_length_clamped = .false.
    if (requested_block_len < 2 * full_max_rl .and. .not. pwb_block_warned(gas)) then
        write(*, '(a,a,a,f8.2,a,f8.2,a)') '  WARNING: PWB block length for ', &
            trim(GasLabel(gas)), ' (', &
            dble(requested_block_len) / Metadata%ac_freq, ' s) is shorter than 2*lag_max (', &
            dble(2 * full_max_rl) / Metadata%ac_freq, ' s).'
        pwb_block_warned(gas) = .true.
    end if
    block_len = min(max(1, block_len), n)
    res%effective_block_length_s = dble(block_len) / Metadata%ac_freq
    nblocks = (n + block_len - 1) / block_len

    res%effective_min_lag = dble(min_rl) / Metadata%ac_freq
    res%effective_max_lag = dble(max_rl) / Metadata%ac_freq
    allocate(boot_lags(nboot), xb(n), yb(n), ccf(full_min_rl:full_max_rl), smooth(full_min_rl:full_max_rl))
    allocate(counts(min_rl:max_rl), mean_smooth(full_min_rl:full_max_rl), mean_ccf(full_min_rl:full_max_rl))
    counts = 0
    mean_ccf = 0d0
    state = PWBSetup%random_seed + 7919 * max(1, gas)
    do i = 1, len_trim(combo)
        state = state + 104729 * i * iachar(combo(i:i))
    end do
    state = max(1, state)

    do b = 1, nboot
        pos = 1
        do i = 1, nblocks
            start = 1 + LcgRandInt(state, max(1, n - block_len + 1))
            call CopyBlock(x, y, n, start, block_len, xb, yb, pos)
            if (pos > n) exit
        end do
        call ComputeCcfWindow(xb, yb, n, full_min_rl, full_max_rl, ccf, PWBSetup%approx_ccf)
        call SmoothAndFill(ccf, full_min_rl, full_max_rl, max(1, PWBSetup%smoothing_width), smooth)
        best_idx = ArgmaxAbs(smooth(min_rl:max_rl), min_rl, max_rl)
        boot_lags(b) = best_idx
        counts(best_idx) = counts(best_idx) + 1
        mean_ccf = mean_ccf + ccf
    end do
    mean_ccf = mean_ccf / dble(nboot)
    call SmoothAndFill(mean_ccf, full_min_rl, full_max_rl, max(1, PWBSetup%smoothing_width), mean_smooth)

    lag = MapLagEstimate(boot_lags, nboot)
    allocate(hdi_samples(nboot))
    do i = 1, nboot
        hdi_samples(i) = dble(boot_lags(i)) / Metadata%ac_freq
    end do
    call Hdi95(hdi_samples, nboot, res%hdi_low, res%hdi_high)
    res%row_lag = lag
    res%selected_lag = dble(lag) / Metadata%ac_freq
    res%hdi_range = res%hdi_high - res%hdi_low
    res%edge_pinned = lag == min_rl .or. lag == max_rl
    res%ccf_at_mode = abs(mean_smooth(lag))
    res%reliability_class = 'detected'
    ok = .not. res%edge_pinned
    deallocate(boot_lags, xb, yb, ccf, smooth, counts, hdi_samples, mean_smooth, mean_ccf)
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

integer function LcgRandInt(state, upper)
    integer, intent(inout) :: state
    integer, intent(in) :: upper
    integer(8) :: s64
    s64 = int(state, 8)
    s64 = mod(1103515245_8 * s64 + 12345_8, 2147483647_8)
    state = int(s64, 4)
    if (state < 0) state = -state
    LcgRandInt = mod(state, upper)
end function LcgRandInt

subroutine ComputeCcfWindow(x, y, n, min_rl, max_rl, ccf, approx)
    integer,  intent(in) :: n, min_rl, max_rl
    real(kind = dbl), intent(in)  :: x(n), y(n)
    real(kind = dbl), intent(out) :: ccf(min_rl:max_rl)
    logical,  intent(in) :: approx
    integer :: lag, i, nn
    real(kind = dbl) :: mx, my, vx, vy, denom, cov
    real(kind = dbl), allocatable :: xc(:), yc(:)
    !> Single-pass CCF using the computational formula (König-Huygens identity).
    !> When approx=.true., variance normalisation is skipped — valid when only
    !> the argmax is needed and N >> lag range (variance varies <1% across lags).
    allocate(xc(n), yc(n))
    mx = sum(x) / dble(n)
    my = sum(y) / dble(n)
    xc = x - mx
    yc = y - my
    vx = sum(xc * xc)
    vy = sum(yc * yc)
    denom = sqrt(vx * vy)
    if (denom <= 0d0) then
        ccf = 0d0
        deallocate(xc, yc)
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
                cov = cov + xc(i) * yc(i + lag)
            end do
        else
            do i = 1, nn
                cov = cov + xc(i - lag) * yc(i)
            end do
        end if
        if (approx) then
            ccf(lag) = cov
        else
            ccf(lag) = cov / denom
        end if
    end do
    deallocate(xc, yc)
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
        y(i) = y(i) / dble(width)
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

    !> Select by highest |mean_smooth_ccf| at mode lag (matching R/Python)
    !> First try non-edge-pinned candidates
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
    !> If all edge-pinned, pick the one with highest CCF anyway
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

subroutine WritePwbDiagnostic(gas, res)
    integer, intent(in) :: gas
    type(PWBResultType), intent(in) :: res
    integer :: u, ios
    character(PathLen) :: path

    call CountPwbDiagnostic(gas, res)
    if (Dir%main_out == 'none') return
    path = Dir%main_out(1:len_trim(Dir%main_out)) &
        // EddyFlowProj%id(1:len_trim(EddyFlowProj%id)) &
        // PwbTimelagDiag_FilePadding // Timestamp_FilePadding // CsvExt
    open(newunit = u, file = path, status = 'unknown', position = 'append', iostat = ios, encoding = 'utf-8')
    if (ios /= 0) return
    if (.not. pwb_diag_header_written) then
        write(u, '(a)') 'date,time,gas,raw_selected_lag_s,raw_row_lag,applied_lag_s,applied_row_lag,hdi_low_s,' &
            // 'hdi_high_s,hdi_range_s,reliability_class,best_combination,' &
            // 'edge_pinned,fallback_used,fallback_source,donor_gas,effective_min_lag_s,effective_max_lag_s,' &
            // 'effective_block_length_s,block_length_clamped,raw_covariance'
        pwb_diag_header_written = .true.
    end if
    write(u, '(a,",",a,",",a,",",f10.4,",",i8,",",f10.4,",",i8,' &
        // '",",f10.4,",",f10.4,",",f10.4,",",a,",",a,",",l1,",",l1,",",a,",",a,' &
        // '",",f10.4,",",f10.4,",",f10.4,",",l1,",",f14.6)') &
        trim(Stats%date), trim(Stats%time), trim(GasLabel(gas)), res%selected_lag, res%row_lag, &
        res%applied_lag, res%applied_row_lag, &
        res%hdi_low, res%hdi_high, res%hdi_range, trim(res%reliability_class), &
        trim(res%best_combination), res%edge_pinned, res%fallback_used, &
        trim(res%fallback_source), trim(res%donor_gas), res%effective_min_lag, res%effective_max_lag, &
        res%effective_block_length_s, res%block_length_clamped, &
        res%raw_covariance
    close(u)
end subroutine WritePwbDiagnostic

subroutine CountPwbDiagnostic(gas, res)
    integer, intent(in) :: gas
    type(PWBResultType), intent(in) :: res

    if (gas < co2 .or. gas > gas4) return
    pwb_attempts(gas) = pwb_attempts(gas) + 1
    if (res%fallback_used) then
        pwb_fallbacks(gas) = pwb_fallbacks(gas) + 1
        select case(trim(res%fallback_source))
            case('maxcov_default')
                pwb_fallback_maxcov(gas) = pwb_fallback_maxcov(gas) + 1
            case('nominal/default')
                pwb_fallback_nominal(gas) = pwb_fallback_nominal(gas) + 1
            case default
                pwb_fallback_other(gas) = pwb_fallback_other(gas) + 1
        end select
    elseif (trim(res%reliability_class) == 'S3_carryforward') then
        pwb_carryforwards(gas) = pwb_carryforwards(gas) + 1
    elseif (trim(res%reliability_class) == 'S4_instrument_shared') then
        pwb_instrument_shared(gas) = pwb_instrument_shared(gas) + 1
    else
        pwb_successes(gas) = pwb_successes(gas) + 1
    end if
end subroutine CountPwbDiagnostic

subroutine ReportPwbDiagnostics()
    integer :: gas, u, ios
    integer :: total_attempts, total_successes, total_carryforwards, total_fallbacks
    integer :: total_instrument_shared
    integer :: total_fallback_maxcov, total_fallback_nominal, total_fallback_other
    character(PathLen) :: path

    total_attempts = sum(pwb_attempts(co2:gas4))
    if (total_attempts == 0) return

    total_successes = sum(pwb_successes(co2:gas4))
    total_instrument_shared = sum(pwb_instrument_shared(co2:gas4))
    total_carryforwards = sum(pwb_carryforwards(co2:gas4))
    total_fallbacks = sum(pwb_fallbacks(co2:gas4))
    total_fallback_maxcov = sum(pwb_fallback_maxcov(co2:gas4))
    total_fallback_nominal = sum(pwb_fallback_nominal(co2:gas4))
    total_fallback_other = sum(pwb_fallback_other(co2:gas4))

    write(*, '(a)')
    write(*, '(a)') ' PWB time-lag detection summary:'
    do gas = co2, gas4
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
        end if
    end do
    if (total_successes == 0 .and. total_instrument_shared == 0 &
        .and. total_carryforwards == 0 .and. total_fallbacks > 0) then
        write(*, '(a, i0, a, i0, a, i0, a)') '  WARNING: all PWB detections fell back: maxcov/default=', &
            total_fallback_maxcov, ', nominal/default=', total_fallback_nominal, &
            ', other=', total_fallback_other, '.'
        write(*, '(a)') '  Review the PWB diagnostics file before interpreting method 5 as native PWB.'
    end if

    if (Dir%main_out == 'none') return
    path = Dir%main_out(1:len_trim(Dir%main_out)) &
        // EddyFlowProj%id(1:len_trim(EddyFlowProj%id)) &
        // PwbSummary_FilePadding // Timestamp_FilePadding // CsvExt
    open(newunit = u, file = path, status = 'replace', iostat = ios, encoding = 'utf-8')
    if (ios /= 0) return
    write(u, '(a)') 'gas,attempts,S1_S2_optimal,S4_instrument_shared,S3_carryforward,' &
        // 'fallback,maxcov_default,nominal_default,other_fallback'
    do gas = co2, gas4
        if (pwb_attempts(gas) > 0) then
            write(u, '(a,",",i0,",",i0,",",i0,",",i0,",",i0,",",i0,",",i0,",",i0)') trim(GasLabel(gas)), &
                pwb_attempts(gas), pwb_successes(gas), pwb_instrument_shared(gas), &
                pwb_carryforwards(gas), pwb_fallbacks(gas), pwb_fallback_maxcov(gas), &
                pwb_fallback_nominal(gas), pwb_fallback_other(gas)
        end if
    end do
    write(u, '(a,",",i0,",",i0,",",i0,",",i0,",",i0,",",i0,",",i0,",",i0)') 'all', total_attempts, &
        total_successes, total_instrument_shared, total_carryforwards, total_fallbacks, &
        total_fallback_maxcov, total_fallback_nominal, total_fallback_other
    close(u)
end subroutine ReportPwbDiagnostics

character(8) function GasLabel(gas)
    integer, intent(in) :: gas
    select case(gas)
        case(co2)
            GasLabel = 'co2'
        case(h2o)
            GasLabel = 'h2o'
        case(ch4)
            GasLabel = 'ch4'
        case(gas4)
            GasLabel = 'gas4'
        case default
            GasLabel = 'unknown'
    end select
end function GasLabel

end module m_pwb_timelag
