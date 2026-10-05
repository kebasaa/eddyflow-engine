!***************************************************************************
! pwb_compare_main.f90
! --------------------
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
! EddyFlow® is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
!***************************************************************************
!
! \brief       Runs the engine's whole full-rate PWB chain for one period on
!              block starts it is given, for tests/pwb_dyco/compare.py.
!
!              pwb_reference checks the deterministic half against RFlux's
!              frozen numbers. This goes the rest of the way: the four
!              block bootstraps, the smoothing, the per-replicate peaks, the
!              mode, the HDI and the choice of combination - the same core
!              routines PwbDetectGas calls, in the same order, with the
!              resamples supplied rather than drawn. compare.py hands dyco
!              the very same resamples, so everything can be compared exactly.
!
!              Input (text): one header line
!                  n hz min_rl max_rl swidth block_req nblocks nboot
!              then n lines "w tsonic scalar" (-9999 = missing), then for
!              each combination cw, wc, ct, tc and each replicate one line of
!              nblocks 1-based block starts.
!
!              Output: key=value lines on stdout.
!
!              Linked against m_pwb_core and m_numeric_kinds only, like
!              pwb_reference: it is what proves this chain needs no engine
!              state, and it stops linking the moment that stops being true.
! \author      Jonathan Muller
!***************************************************************************
program pwb_compare_main
    use m_numeric_kinds
    use m_pwb_core
    implicit none

    character(2), parameter :: combo(4) = (/'cw', 'wc', 'ct', 'tc'/)
    real(kind = dbl), parameter :: missing = -9999d0
    character(1024) :: path
    integer :: u, ios, n, i, c, b, k
    integer :: min_rl, max_rl, swidth, block_req, nblocks, nboot
    integer :: trail, margin, eval_lo, eval_hi, widest, block_len, nb_eng
    integer :: n_eff, lag(4), unrestricted(4), best
    real(kind = dbl) :: hz
    real(kind = dbl) :: hdi_lo(4), hdi_hi(4), ccf_at_mode(4)
    logical :: ok(4)
    real(kind = dbl), allocatable :: ww(:), tt(:), ss(:)
    real(kind = dbl), allocatable :: s_fs(:), w_fs(:), t_fs(:)
    real(kind = dbl), allocatable :: s_fw(:), w_fw(:), s_ft(:), t_ft(:)
    real(kind = dbl), allocatable :: raw_ccov(:), xc(:), yc(:), xb(:), yb(:)
    real(kind = dbl), allocatable :: ccf(:), smooth(:), mean_ccf(:), hdi_buf(:)
    real(kind = dbl), allocatable :: mean_smooth(:, :)
    integer, allocatable :: starts(:, :, :), boot(:, :)
    type(PwbPreWhitenType) :: pw

    if (command_argument_count() < 1) then
        write(*, '(a)') 'usage: pwb_compare <input>'
        stop 2
    end if
    call get_command_argument(1, path)
    open(newunit = u, file = trim(path), status = 'old', action = 'read', iostat = ios)
    if (ios /= 0) then
        write(*, '(a)') 'error: cannot open ' // trim(path)
        stop 3
    end if
    read(u, *) n, hz, min_rl, max_rl, swidth, block_req, nblocks, nboot
    allocate(ww(n), tt(n), ss(n))
    do i = 1, n
        read(u, *) ww(i), tt(i), ss(i)
    end do
    allocate(starts(nblocks, nboot, 4))
    do c = 1, 4
        do b = 1, nboot
            read(u, *) (starts(k, b, c), k = 1, nblocks)
        end do
    end do
    close(u)

    !> PwbDetectGas, step by step.
    call FillMissingLinear(ww, n, missing)
    call FillMissingLinear(tt, n, missing)
    call FillMissingLinear(ss, n, missing)

    trail = max(1, swidth) / 2
    margin = max(trail, nint(2d0 * hz))
    eval_lo = max(min_rl - margin, -(n - 3))
    eval_hi = min(max_rl + margin, n - 3)

    allocate(s_fs(n), w_fs(n), t_fs(n), s_fw(n), w_fw(n), s_ft(n), t_ft(n))
    allocate(raw_ccov(min_rl:max_rl), xc(n), yc(n), xb(n), yb(n))
    call PwbPreWhiten(ss, ww, tt, n, min_rl, max_rl, missing, pw, &
        s_fs, w_fs, t_fs, s_fw, w_fw, s_ft, t_ft, raw_ccov, xc, yc)
    n_eff = pw%n_eff

    !> RunPwbCombination's block length, from n_eff as there.
    widest = max(abs(min_rl), abs(max_rl))
    block_len = max(block_req, 2 * widest)
    block_len = min(max(1, block_len), n_eff)
    nb_eng = (n_eff + block_len - 1) / block_len
    if (nb_eng /= nblocks) then
        write(*, '(a,i0,a,i0)') 'error: engine wants nblocks=', nb_eng, ' got ', nblocks
        stop 4
    end if
    if (maxval(starts) > n_eff - block_len + 1 .or. minval(starts) < 1) then
        write(*, '(a)') 'error: a block start lies outside [1, n_eff - block_len + 1]'
        stop 5
    end if

    allocate(ccf(eval_lo:eval_hi), smooth(eval_lo:eval_hi), mean_ccf(eval_lo:eval_hi))
    allocate(mean_smooth(eval_lo:eval_hi, 4), boot(nboot, 4), hdi_buf(nboot))
    do c = 1, 4
        select case (c)
            case (1)
                call RunOne(w_fs, s_fs)
            case (2)
                call RunOne(w_fw, s_fw)
            case (3)
                call RunOne(t_fs, s_fs)
            case (4)
                call RunOne(t_ft, s_ft)
        end select
    end do
    best = PwbBestCombination(ccf_at_mode, ok, 4)

    write(*, '(a,i0)')      'n=', n
    write(*, '(a,i0)')      'n_eff=', n_eff
    write(*, '(a,l1)')      'differenced=', pw%differenced
    write(*, '(a,i0)')      'ar_order_scalar=', pw%p_scalar
    write(*, '(a,i0)')      'ar_order_w=', pw%p_w
    write(*, '(a,i0)')      'ar_order_t=', pw%p_t
    write(*, '(a,es24.15)') 'phi1_scalar=', pw%phi1_scalar
    write(*, '(a,es24.15)') 'phi1_w=', pw%phi1_w
    write(*, '(a,es24.15)') 'phi1_t=', pw%phi1_t
    write(*, '(a,i0)')      'tlag_pw=', pw%tlag_pw_rl
    write(*, '(a,es24.15)') 'corr_pw=', pw%corr_pw
    write(*, '(a,i0)')      'eval_lo=', eval_lo
    write(*, '(a,i0)')      'eval_hi=', eval_hi
    write(*, '(a,i0)')      'block_len=', block_len
    do c = 1, 4
        write(*, '(a,a,a,i0)')      'mode_', combo(c), '=', lag(c)
        write(*, '(a,a,a,es24.15)') 'hdi_lo_', combo(c), '=', hdi_lo(c)
        write(*, '(a,a,a,es24.15)') 'hdi_hi_', combo(c), '=', hdi_hi(c)
        write(*, '(a,a,a,es24.15)') 'ccf_at_mode_', combo(c), '=', ccf_at_mode(c)
        write(*, '(a,a,a,l1)')      'ok_', combo(c), '=', ok(c)
        write(*, '(a,a,a)', advance = 'no') 'lags_', combo(c), '='
        do b = 1, nboot
            write(*, '(i0,a)', advance = 'no') boot(b, c), ' '
        end do
        write(*, '(a)')
        write(*, '(a,a,a)', advance = 'no') 'mean_smooth_', combo(c), '='
        do k = eval_lo, eval_hi
            write(*, '(es24.15,a)', advance = 'no') mean_smooth(k, c), ' '
        end do
        write(*, '(a)')
    end do
    write(*, '(a,a)') 'best=', combo(best)

contains

    subroutine RunOne(x, y)
        real(kind = dbl), intent(in) :: x(n), y(n)

        call PwbBootstrapCombination(x(1:n_eff), y(1:n_eff), n_eff, min_rl, max_rl, &
            eval_lo, eval_hi, swidth, block_len, nblocks, nboot, starts(:, :, c), &
            boot(:, c), mean_smooth(:, c), xb(1:n_eff), yb(1:n_eff), xc(1:n_eff), &
            yc(1:n_eff), ccf, smooth, mean_ccf)
        call PwbSummariseBootstrap(boot(:, c), nboot, hz, eval_lo, eval_hi, &
            mean_smooth(:, c), hdi_buf, lag(c), hdi_lo(c), hdi_hi(c), &
            ccf_at_mode(c), unrestricted(c))
        ok(c) = lag(c) /= min_rl .and. lag(c) /= max_rl
    end subroutine RunOne
end program pwb_compare_main
