!***************************************************************************
! ccf_bitexact_main.f90
! ---------------------
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
! \brief       Proves the blocked cross-correlation is bit for bit the one it
!              replaced.
!
!              ComputeCcfWindow is the inner loop of the PWB pre-pass and was
!              rewritten to sum a block of lags at once. The claim that made
!              that acceptable is that every lag keeps its own accumulator and
!              receives the same products in the same order, so no result can
!              move by even one bit. This holds the old routine, frozen exactly
!              as it was, and compares the two on the raw bit pattern of every
!              element: thousands of random shapes, every window edge for small
!              n, a zero-variance series, and the production shape.
!
!              "Close" is not a pass. A routine that differs in the last bit
!              moves a lag only rarely, and an error that is rare is the kind
!              no fixture ever finds.
!
! \note        Built by `mingw32-make ccfcheck`, deliberately not in `all`, and
!              linked against m_pwb_core.o itself - so what is checked is the
!              shipped routine, compiled with the shipped flags.
!              static_checks/test_pwb_ccf_blocked_static.py runs it.
! \author      Jonathan Muller
!***************************************************************************
module ccf_frozen_reference
    use m_numeric_kinds
    implicit none
    private
    public :: RefCcf
contains
    !> ComputeCcfWindow as it stood before the blocked rewrite. Do not touch:
    !> it is the definition of the right answer.
    subroutine RefCcf(x, y, n, min_rl, max_rl, ccf, xc, yc)
        integer, intent(in) :: n, min_rl, max_rl
        real(kind = dbl), intent(in) :: x(n), y(n)
        real(kind = dbl), intent(out) :: ccf(min_rl:max_rl)
        real(kind = dbl), intent(inout) :: xc(n), yc(n)
        integer :: lag, i, nn
        real(kind = dbl) :: mx, my, vx, vy, denom, cov

        mx = sum(x) / dble(n)
        my = sum(y) / dble(n)
        xc = x - mx
        yc = y - my
        vx = sum(xc * xc)
        vy = sum(yc * yc)
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
                    cov = cov + xc(i) * yc(i + lag)
                end do
            else
                do i = 1, nn
                    cov = cov + xc(i - lag) * yc(i)
                end do
            end if
            ccf(lag) = cov / denom
        end do
    end subroutine RefCcf
end module ccf_frozen_reference

program ccf_bitexact_main
    use, intrinsic :: iso_fortran_env, only: int64
    use m_numeric_kinds
    use m_pwb_core, only: ComputeCcfWindow
    use ccf_frozen_reference
    implicit none
    integer :: trial, n, lo, hi, i, cases, mismatched
    real(kind = dbl) :: r, t_ref, t_new
    real(kind = dbl), allocatable :: x(:), y(:), xc(:), yc(:), a(:), b(:)
    integer(int64) :: c0, c1, rate
    integer, allocatable :: seed(:)

    call random_seed(size = i)
    allocate(seed(i))
    seed = 20261001
    call random_seed(put = seed)
    cases = 0
    mismatched = 0

    !> Windows wholly negative, wholly positive and straddling zero; reaching
    !> past +-(n-1), where a lag has fewer than two records; block sizes that
    !> are and are not multiples of the lag block.
    do trial = 1, 4000
        call random_number(r)
        n = 2 + int(r * 1500)
        call random_number(r)
        lo = -n - 3 + int(r * (2 * n + 6))
        call random_number(r)
        hi = lo + int(r * (n + 40))
        call check(n, lo, hi, .false.)
    end do
    !> Small n, every window: the edges, exhaustively.
    do n = 2, 12
        do lo = -n - 2, n + 2
            do hi = lo, n + 2
                call check(n, lo, hi, .false.)
            end do
        end do
    end do
    !> Zero variance: every element must come back 0.
    call check(500, -20, 40, .true.)
    !> The production shape: 20 Hz x 30 min, a 0-25 s window, 2 s guard band.
    do trial = 1, 20
        call check(36000, -40, 540, .false.)
    end do

    write(*, '(a,i0,a,i0)') ' cases: ', cases, '   mismatched elements: ', mismatched

    !> Speed at the production shape: 4 combinations x 99 replicates.
    n = 36000
    lo = -40
    hi = 540
    call fresh(n, lo, hi, .false.)
    call system_clock(c0, rate)
    do i = 1, 396
        call RefCcf(x, y, n, lo, hi, a, xc, yc)
    end do
    call system_clock(c1)
    t_ref = dble(c1 - c0) / dble(rate)
    call system_clock(c0, rate)
    do i = 1, 396
        call ComputeCcfWindow(x, y, n, lo, hi, b, xc, yc)
    end do
    call system_clock(c1)
    t_new = dble(c1 - c0) / dble(rate)
    write(*, '(a,f8.3,a,f8.3,a,f6.2,a)') ' 396 calls at the production shape: old ', &
        t_ref, ' s, new ', t_new, ' s (', t_ref / max(t_new, 1d-9), 'x)'

    if (mismatched == 0) then
        write(*, '(a)') ' BIT-IDENTICAL'
    else
        write(*, '(a)') ' MISMATCH'
        error stop 1
    end if

contains

    subroutine fresh(nn, l, h, constant)
        integer, intent(in) :: nn, l, h
        logical, intent(in) :: constant
        if (allocated(x)) deallocate(x, y, xc, yc, a, b)
        allocate(x(nn), y(nn), xc(nn), yc(nn), a(l:h), b(l:h))
        call random_number(x)
        call random_number(y)
        x = (x - 0.5d0) * 7.3d0
        y = y * 1d-3 + 412.7d0
        if (constant) x = 3.25d0
    end subroutine fresh

    subroutine check(nn, l, h, constant)
        integer, intent(in) :: nn, l, h
        logical, intent(in) :: constant
        integer :: k, bad
        call fresh(nn, l, h, constant)
        a = -1d300
        b = 1d300
        call RefCcf(x, y, nn, l, h, a, xc, yc)
        call ComputeCcfWindow(x, y, nn, l, h, b, xc, yc)
        bad = 0
        do k = l, h
            if (transfer(a(k), 0_int64) /= transfer(b(k), 0_int64)) bad = bad + 1
        end do
        if (bad > 0 .and. mismatched == 0) &
            write(*, '(a,i0,a,i0,a,i0,a,i0)') ' first mismatch: n=', nn, &
                ' window ', l, ':', h, ', elements ', bad
        mismatched = mismatched + bad
        cases = cases + 1
    end subroutine check

end program ccf_bitexact_main
