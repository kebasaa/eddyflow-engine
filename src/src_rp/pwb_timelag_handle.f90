!***************************************************************************
! pwb_timelag_handle.f90
! ----------------------
! Native implementation of the Vitale et al. (2024) pre-whitening
! block-bootstrap time lag detector, ported from the trusted Python reference
! in eddyflow-build-script/literature/pre-whitening block-bootstrap PWB.
!***************************************************************************
module m_pwb_timelag
    use m_rp_global_var
    implicit none
    private
    public :: PwbDetectGas, ResetPwbDiagnostics

    logical :: pwb_diag_header_written = .false.

contains

subroutine ResetPwbDiagnostics()
    pwb_diag_header_written = .false.
end subroutine ResetPwbDiagnostics

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
        call WritePwbDiagnostic(gas, LocResult)
        return
    end if

    min_rl = nint(PWBSetup%min_lag(gas) * Metadata%ac_freq)
    max_rl = nint(PWBSetup%max_lag(gas) * Metadata%ac_freq)
    if (min_rl >= max_rl) then
        LocResult%fallback_used = .true.
        call WritePwbDiagnostic(gas, LocResult)
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
        call WritePwbDiagnostic(gas, LocResult)
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
    call RunPwbCombination(w_fs, s_fs, nrow, min_rl, max_rl, combo(1), candidate(1), ok(1))
    call RunPwbCombination(w_fw, s_fw, nrow, min_rl, max_rl, combo(2), candidate(2), ok(2))
    call RunPwbCombination(t_fs, s_fs, nrow, min_rl, max_rl, combo(3), candidate(3), ok(3))
    call RunPwbCombination(t_ft, s_ft, nrow, min_rl, max_rl, combo(4), candidate(4), ok(4))

    call SelectBestCandidate(candidate, ok, min_rl, max_rl, LocResult, success)
    if (success) then
        allocate(raw_ccov(min_rl:max_rl))
        call ComputeCcovWindow(ww, ss, nrow, min_rl, max_rl, raw_ccov)
        lag = LocResult%row_lag
        if (lag >= min_rl .and. lag <= max_rl) LocResult%raw_covariance = raw_ccov(lag)
        deallocate(raw_ccov)
    else
        LocResult%fallback_used = .true.
    end if

    call WritePwbDiagnostic(gas, LocResult)

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
    res%hdi_low = error
    res%hdi_high = error
    res%hdi_range = error
    res%reliability_class = 'failed'
    res%best_combination = '--'
    res%edge_pinned = .false.
    res%fallback_used = .false.
    res%raw_covariance = error
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

subroutine RunPwbCombination(x, y, n, min_rl, max_rl, combo, res, ok)
    integer, intent(in) :: n, min_rl, max_rl
    real(kind = dbl), intent(in) :: x(n), y(n)
    character(2), intent(in) :: combo
    type(PWBResultType), intent(out) :: res
    logical, intent(out) :: ok
    integer :: b, i, pos, block_len, nblocks, start, state
    integer :: nboot, lag, best_idx
    integer, allocatable :: boot_lags(:), counts(:)
    real(kind = dbl), allocatable :: xb(:), yb(:), ccf(:), smooth(:), hdi_samples(:)

    call InitPwbResult(res)
    res%best_combination = combo
    ok = .false.
    nboot = max(1, PWBSetup%n_bootstrap)
    block_len = nint(PWBSetup%block_length_s * Metadata%ac_freq)
    if (block_len <= 0) block_len = max(1, 2 * max(abs(min_rl), abs(max_rl)))
    block_len = min(max(1, block_len), n)
    nblocks = (n + block_len - 1) / block_len

    allocate(boot_lags(nboot), xb(n), yb(n), ccf(min_rl:max_rl), smooth(min_rl:max_rl))
    allocate(counts(min_rl:max_rl))
    counts = 0
    state = max(1, PWBSetup%random_seed + 7919 * max(1, gas4 - co2 + 1) + iachar(combo(1:1)))

    do b = 1, nboot
        pos = 1
        do i = 1, nblocks
            start = 1 + LcgRandInt(state, max(1, n - block_len + 1))
            call CopyBlock(x, y, n, start, block_len, xb, yb, pos)
            if (pos > n) exit
        end do
        call ComputeCcfWindow(xb, yb, n, min_rl, max_rl, ccf)
        call SmoothAndFill(ccf, min_rl, max_rl, max(1, PWBSetup%smoothing_width), smooth)
        best_idx = ArgmaxAbs(smooth, min_rl, max_rl)
        boot_lags(b) = best_idx
        counts(best_idx) = counts(best_idx) + 1
    end do

    lag = min_rl
    do i = min_rl, max_rl
        if (counts(i) > counts(lag)) lag = i
    end do
    allocate(hdi_samples(nboot))
    do i = 1, nboot
        hdi_samples(i) = dble(boot_lags(i)) / Metadata%ac_freq
    end do
    call Hdi95(hdi_samples, nboot, res%hdi_low, res%hdi_high)
    res%row_lag = lag
    res%selected_lag = dble(lag) / Metadata%ac_freq
    res%hdi_range = res%hdi_high - res%hdi_low
    res%edge_pinned = lag == min_rl .or. lag == max_rl
    if (res%edge_pinned) then
        res%reliability_class = 'edge_pinned'
        ok = .false.
    elseif (PWBSetup%hdi_prefilter_s > 0d0 .and. res%hdi_range > PWBSetup%hdi_prefilter_s) then
        res%reliability_class = 'prefiltered'
        ok = .false.
    elseif (res%hdi_range < PWBSetup%hdi_thresh_s) then
        res%reliability_class = 'S1_optimal'
        ok = .true.
    elseif (res%hdi_range < PWBSetup%hdi_thresh_s + PWBSetup%dev_thresh_s) then
        res%reliability_class = 'S2_optimal'
        ok = .true.
    else
        res%reliability_class = 'uncertain'
        ok = .false.
    end if
    deallocate(boot_lags, xb, yb, ccf, smooth, counts, hdi_samples)
end subroutine RunPwbCombination

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
    state = mod(1103515245 * state + 12345, 2147483647)
    if (state < 0) state = -state
    LcgRandInt = mod(state, upper)
end function LcgRandInt

subroutine ComputeCcfWindow(x, y, n, min_rl, max_rl, ccf)
    integer, intent(in) :: n, min_rl, max_rl
    real(kind = dbl), intent(in) :: x(n), y(n)
    real(kind = dbl), intent(out) :: ccf(min_rl:max_rl)
    integer :: lag, i, nn, s1, s2
    real(kind = dbl) :: mx, my, vx, vy, cov
    do lag = min_rl, max_rl
        mx = 0d0; my = 0d0; vx = 0d0; vy = 0d0; cov = 0d0; nn = n - abs(lag)
        if (nn <= 1) then
            ccf(lag) = 0d0
            cycle
        end if
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
        do i = 1, nn
            if (lag >= 0) then
                s1 = i; s2 = i + lag
            else
                s1 = i - lag; s2 = i
            end if
            cov = cov + (x(s1) - mx) * (y(s2) - my)
            vx = vx + (x(s1) - mx)**2
            vy = vy + (y(s2) - my)**2
        end do
        if (vx <= 0d0 .or. vy <= 0d0) then
            ccf(lag) = 0d0
        else
            ccf(lag) = cov / sqrt(vx * vy)
        end if
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

subroutine SmoothAndFill(x, min_rl, max_rl, width, y)
    integer, intent(in) :: min_rl, max_rl, width
    real(kind = dbl), intent(in) :: x(min_rl:max_rl)
    real(kind = dbl), intent(out) :: y(min_rl:max_rl)
    integer :: i, j, lo, hi, n
    do i = min_rl, max_rl
        lo = max(min_rl, i - width / 2)
        hi = min(max_rl, i + width / 2)
        y(i) = 0d0
        n = 0
        do j = lo, hi
            y(i) = y(i) + x(j)
            n = n + 1
        end do
        if (n > 0) y(i) = y(i) / dble(n)
    end do
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

subroutine SelectBestCandidate(candidate, ok, min_rl, max_rl, res, success)
    type(PWBResultType), intent(in) :: candidate(4)
    logical, intent(in) :: ok(4)
    integer, intent(in) :: min_rl, max_rl
    type(PWBResultType), intent(out) :: res
    logical, intent(out) :: success
    integer :: i, best

    call InitPwbResult(res)
    success = .false.
    best = 0
    do i = 1, 4
        if (ok(i)) then
            if (best == 0) then
                best = i
            elseif (candidate(i)%hdi_range < candidate(best)%hdi_range) then
                best = i
            end if
        end if
    end do
    if (best == 0) then
        do i = 1, 4
            if (best == 0) then
                best = i
            elseif (candidate(i)%hdi_range /= error .and. candidate(i)%hdi_range < candidate(best)%hdi_range) then
                best = i
            end if
        end do
        if (best > 0) res = candidate(best)
        res%fallback_used = .true.
        return
    end if
    res = candidate(best)
    success = res%row_lag > min_rl .and. res%row_lag < max_rl
end subroutine SelectBestCandidate

subroutine WritePwbDiagnostic(gas, res)
    integer, intent(in) :: gas
    type(PWBResultType), intent(in) :: res
    integer :: u, ios
    character(PathLen) :: path

    if (Dir%main_out == 'none') return
    path = Dir%main_out(1:len_trim(Dir%main_out)) // 'eddyflow_pwb_timelag_diagnostics.csv'
    open(newunit = u, file = path, status = 'unknown', position = 'append', iostat = ios, encoding = 'utf-8')
    if (ios /= 0) return
    if (.not. pwb_diag_header_written) then
        write(u, '(a)') 'date,time,gas,selected_lag_s,row_lag,hdi_low_s,' &
            // 'hdi_high_s,hdi_range_s,reliability_class,best_combination,' &
            // 'edge_pinned,fallback_used,raw_covariance'
        pwb_diag_header_written = .true.
    end if
    write(u, '(a,",",a,",",a,",",f10.4,",",i8,",",f10.4,",",f10.4,",",f10.4,",",a,",",a,",",l1,",",l1,",",f14.6)') &
        trim(Stats%date), trim(Stats%time), trim(GasLabel(gas)), res%selected_lag, res%row_lag, &
        res%hdi_low, res%hdi_high, res%hdi_range, trim(res%reliability_class), &
        trim(res%best_combination), res%edge_pinned, res%fallback_used, res%raw_covariance
    close(u)
end subroutine WritePwbDiagnostic

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
