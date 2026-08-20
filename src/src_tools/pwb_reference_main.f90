!***************************************************************************
! pwb_reference_main.f90
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
! \brief       Run the deterministic PWB chain over one CSV and print what it
!              found, so it can be compared with RFlux's output for the same
!              input.
!
! \details     Not part of either shipped executable. It exists so that
!              static_checks/test_pwb_reference_static.py can drive exactly
!              the arithmetic the engine uses against the frozen R values
!              dyco pins its own implementation to.
!
!              It links m_pwb_core and m_numeric_kinds and nothing else. That
!              is the check: if a routine in the chain ever reaches for
!              PWBSetup, E2Col or the log, this program stops linking, and
!              the numbers stop being reproducible outside a full run.
!
!              Usage:  pwb_reference <csv> <hz> <lag_max_s>
!                      pwb_reference smooth <width> <n>
!
!              The CSV has a header line and three columns, scalar,tsonic,w.
!
!              The second mode smooths a fixed synthetic series and prints it.
!              It exists because no frozen R value touches the smoother: the
!              pinned quantities come from the UNSMOOTHED cross-correlation and
!              from the raw cross-covariance, which is how a wrong window for
!              even widths survived the reference test.
!
! \author      Jonathan Muller
! \note
! \sa          dyco-main/tests/test_pwb_reference.py
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
program pwb_reference_main
    use m_numeric_kinds
    use m_pwb_core
    implicit none

    character(1024) :: path, arg
    integer :: u, ios, n, i, hz, lag_rl
    real(kind = dbl) :: lag_max_s
    real(kind = dbl), allocatable :: ss(:), tt(:), ww(:)
    real(kind = dbl), allocatable :: s_fs(:), w_fs(:), t_fs(:)
    real(kind = dbl), allocatable :: s_fw(:), w_fw(:), s_ft(:), t_ft(:)
    real(kind = dbl), allocatable :: ccov(:), xc(:), yc(:)
    type(PwbPreWhitenType) :: res
    !> Matches the engine's own missing-value code, so the one branch that
    !> writes it behaves identically here.
    real(kind = dbl), parameter :: missing = -9999d0

    if (command_argument_count() < 3) then
        write(*, '(a)') 'usage: pwb_reference <csv> <hz> <lag_max_s>'
        write(*, '(a)') '       pwb_reference smooth <width> <n>'
        stop 2
    end if
    call get_command_argument(1, path)
    if (trim(path) == 'smooth') then
        !> Plain `stop`, not `stop 0`: gfortran prints "STOP 0" for the latter,
        !> which lands in the harness's captured stderr for no reason.
        call SmoothMode()
        stop
    end if
    call get_command_argument(2, arg)
    read(arg, *) hz
    call get_command_argument(3, arg)
    read(arg, *) lag_max_s
    lag_rl = nint(lag_max_s * dble(hz))

    !> Two passes: count, then read. The fixtures are small and this keeps the
    !> program free of a growth policy that could differ from the engine's.
    open(newunit=u, file=trim(path), status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write(*, '(a)') 'error: cannot open ' // trim(path)
        stop 3
    end if
    read(u, '(a)', iostat=ios) arg
    n = 0
    do
        read(u, '(a)', iostat=ios) arg
        if (ios /= 0) exit
        if (len_trim(arg) == 0) cycle
        n = n + 1
    end do
    rewind(u)
    allocate(ss(n), tt(n), ww(n))
    read(u, '(a)', iostat=ios) arg
    i = 0
    do
        read(u, '(a)', iostat=ios) arg
        if (ios /= 0) exit
        if (len_trim(arg) == 0) cycle
        i = i + 1
        call ParseTriplet(arg, ss(i), tt(i), ww(i))
    end do
    close(u)

    allocate(s_fs(n), w_fs(n), t_fs(n), s_fw(n), w_fw(n), s_ft(n), t_ft(n))
    allocate(ccov(-lag_rl:lag_rl), xc(n), yc(n))

    call PwbPreWhiten(ss, ww, tt, n, -lag_rl, lag_rl, missing, res, &
        s_fs, w_fs, t_fs, s_fw, w_fw, s_ft, t_ft, ccov, xc, yc)

    write(*, '(a,i0)')       'n=', n
    write(*, '(a,l1)')       'differenced=', res%differenced
    write(*, '(a,i0)')       'ar_order_scalar=', res%p_scalar
    write(*, '(a,i0)')       'ar_order_w=', res%p_w
    write(*, '(a,i0)')       'ar_order_t=', res%p_t
    write(*, '(a,es24.15)')  'phi1_scalar=', res%phi1_scalar
    write(*, '(a,es24.15)')  'phi1_w=', res%phi1_w
    write(*, '(a,es24.15)')  'phi1_t=', res%phi1_t
    write(*, '(a,i0)')       'pww=', res%tlag_pw_rl
    write(*, '(a,es24.15)')  'cor_pww=', res%corr_pw
    write(*, '(a,i0)')       'mcw=', res%ccov_rl
    write(*, '(a,es24.15)')  'cov_mcw=', res%cov_mcw

    deallocate(ss, tt, ww, s_fs, w_fs, t_fs, s_fw, w_fw, s_ft, t_ft, ccov, xc, yc)

contains

!***************************************************************************
!> Smooth a fixed synthetic series and print it, one value per line.
!>
!> A ramp with a single spike on it. The ramp makes an off-by-one in the
!> window bounds visible as a constant offset, and the spike makes the
!> lead/trail asymmetry of an even width visible as a shifted plateau - which
!> a symmetric window of the same nominal width would place differently.
!***************************************************************************
subroutine SmoothMode()
    character(1024) :: a2, a3
    integer :: width, nn, i
    real(kind = dbl), allocatable :: v(:), sm(:)

    call get_command_argument(2, a2)
    read(a2, *) width
    call get_command_argument(3, a3)
    read(a3, *) nn
    if (width < 1 .or. nn < 1) then
        write(*, '(a)') 'error: width and n must both be positive'
        stop 2
    end if

    allocate(v(nn), sm(nn))
    do i = 1, nn
        v(i) = dble(i)
    end do
    v(max(1, nn / 2)) = v(max(1, nn / 2)) + 100d0

    call SmoothAndFill(v, 1, nn, width, sm)
    write(*, '(a,i0)') 'width=', width
    write(*, '(a,i0)') 'n=', nn
    do i = 1, nn
        write(*, '(a,i0,a,es24.15)') 'y[', i, ']=', sm(i)
    end do
    deallocate(v, sm)
end subroutine SmoothMode

!> One comma-separated line of three reals. List-directed input would do it,
!> but only where the separator is a comma AND the locale agrees; splitting on
!> the comma explicitly keeps the fixture format the fixture's business.
subroutine ParseTriplet(line, a, b, c)
    character(*), intent(in) :: line
    real(kind = dbl), intent(out) :: a, b, c
    integer :: p1, p2

    p1 = index(line, ',')
    p2 = index(line(p1+1:), ',') + p1
    read(line(1:p1-1), *) a
    read(line(p1+1:p2-1), *) b
    read(line(p2+1:), *) c
end subroutine ParseTriplet

end program pwb_reference_main
