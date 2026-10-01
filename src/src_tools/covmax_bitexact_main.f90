!***************************************************************************
! covmax_bitexact_main.f90
! ------------------------
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
! \brief       Proves CrossCovarianceSeries is bit for bit the per-lag path
!              CovMax used to take.
!
!              That path, frozen here exactly as it was: for each lag, copy the
!              shifted pair into an array, copy that into another, hand it to
!              CovarianceMatrixNoError and read Cov(1,2) from the 2x2 matrix it
!              returns. The replacement reads the shifted elements in place and
!              computes that one entry. Every element is compared on its raw
!              bits, over error-code densities from none to all - the QCLS
!              gases of the Yatir project are 95 % error code - and windows
!              straddling zero and wider than the series.
!
! \note        Built by `mingw32-make covmaxcheck`, not in `all`, linked
!              against the shipped m_covmax_core.o.
!              static_checks/test_covmax_series_static.py runs it.
! \author      Jonathan Muller
!***************************************************************************
module covmax_frozen_reference
    use m_numeric_kinds
    implicit none
    private
    public :: FrozenSeries
contains

    !> CovMax's loop over lags, non-detrending path, as it stood. Do not touch.
    subroutine FrozenSeries(Col1, Col2, nrow, lagmin, lagmax, err, out)
        integer, intent(in) :: nrow, lagmin, lagmax
        real(kind = dbl), intent(in) :: Col1(nrow), Col2(nrow), err
        real(kind = dbl), intent(out) :: out(lagmin:lagmax)
        integer :: i, ii, N2
        real(kind = dbl), allocatable :: ShSet(:, :)
        real(kind = dbl), allocatable :: ShPrimes(:, :)
        real(kind = dbl) :: CovMat(2, 2)

        do i = lagmin, lagmax
            N2 = nrow - abs(i)
            allocate(ShSet(N2, 2))
            allocate(ShPrimes(N2, 2))
            do ii = 1, N2
                if (i < 0) then
                    ShSet(ii, 1) = Col1(ii - i)
                    ShSet(ii, 2) = Col2(ii)
                else
                    ShSet(ii, 1) = Col1(ii)
                    ShSet(ii, 2) = Col2(ii + i)
                end if
            end do
            ShPrimes = ShSet
            call FrozenCovMatrix(ShPrimes, size(ShPrimes, 1), size(ShPrimes, 2), CovMat, err)
            out(i) = CovMat(1, 2)
            deallocate(ShSet)
            deallocate(ShPrimes)
        end do
    end subroutine FrozenSeries

    !> CovarianceMatrixNoError, verbatim but for the module it takes dbl from.
    subroutine FrozenCovMatrix(Set, nrow, ncol, Cov, err_float)
        use m_numeric_kinds
        implicit none
        !> in/out variables
        integer, intent(in) :: nrow, ncol
        real(kind = dbl), intent(in) :: err_float
        real(kind = dbl), intent(in) :: Set(nrow, ncol)
        real(kind = dbl), intent(out) :: Cov(ncol, ncol)
        !> local variables
        integer :: i = 0
        integer :: j = 0
        integer :: k = 0
        integer :: Nact = 0
        real(kind = dbl) :: sumi
        real(kind = dbl) :: sumj

        do i = 1, ncol
            do j = 1, ncol
                sumi = 0d0
                sumj = 0d0
                Cov(i, j) = 0d0
                Nact = 0
                do k = 1, nrow
                    if (Set(k, i) /= err_float .and. Set(k, j) /= err_float) then
                        Nact = Nact + 1
                        Cov(i, j) = Cov(i, j) + Set(k, i) * Set(k, j)
                        sumi = sumi + Set(k, i)
                        sumj = sumj + Set(k, j)
                    end if
                end do
                if (Nact /= 0) then
                    sumi = sumi / dble(Nact)
                    sumj = sumj / dble(Nact)
                    Cov(i, j) = Cov(i, j) / dble(Nact)
                    Cov(i, j) = Cov(i, j) - sumi * sumj
                else
                    Cov(i, j) = err_float
                end if
            end do
        end do
    end subroutine FrozenCovMatrix

end module covmax_frozen_reference

program covmax_bitexact_main
    use, intrinsic :: iso_fortran_env, only: int64
    use m_numeric_kinds
    use m_covmax_core, only: CrossCovarianceSeries
    use covmax_frozen_reference
    implicit none
    real(kind = dbl), parameter :: err = -9999d0
    real(kind = dbl), parameter :: densities(6) = [0d0, 0.01d0, 0.3d0, 0.95d0, 0.999d0, 1d0]
    integer :: trial, n, lo, hi, i, d, cases, mismatched
    real(kind = dbl) :: r, t_ref, t_new
    real(kind = dbl), allocatable :: x(:), y(:), a(:), b(:)
    integer(int64) :: c0, c1, rate
    integer, allocatable :: seed(:)

    call random_seed(size = i)
    allocate(seed(i))
    seed = 20261001
    call random_seed(put = seed)
    cases = 0
    mismatched = 0

    do d = 1, size(densities)
        do trial = 1, 600
            call random_number(r)
            n = 1 + int(r * 900)
            call random_number(r)
            lo = -n - 3 + int(r * (2 * n + 6))
            call random_number(r)
            hi = lo + int(r * (n + 40))
            call check(n, lo, hi, densities(d))
        end do
        !> Every window for small n.
        do n = 1, 9
            do lo = -n - 2, n + 2
                do hi = lo, n + 2
                    call check(n, lo, hi, densities(d))
                end do
            end do
        end do
    end do
    !> The production shape at Yatir: 36 000 rows at 20 Hz, a 0-25 s window, a
    !> gas present on one row in twenty.
    do trial = 1, 6
        call check(36000, 0, 500, 0.95d0)
        call check(36000, 0, 500, 0d0)
    end do

    write(*, '(a,i0,a,i0)') ' cases: ', cases, '   mismatched elements: ', mismatched

    n = 36000
    lo = 0
    hi = 500
    call fresh(n, lo, hi, 0.95d0)
    call system_clock(c0, rate)
    call FrozenSeries(x, y, n, lo, hi, err, a)
    call system_clock(c1)
    t_ref = dble(c1 - c0) / dble(rate)
    call system_clock(c0, rate)
    call CrossCovarianceSeries(x, y, n, lo, hi, err, b)
    call system_clock(c1)
    t_new = dble(c1 - c0) / dble(rate)
    write(*, '(a,f8.4,a,f8.4,a,f7.2,a)') ' one gas, one period at Yatir shape: old ', &
        t_ref, ' s, new ', t_new, ' s (', t_ref / max(t_new, 1d-9), 'x)'

    if (mismatched == 0) then
        write(*, '(a)') ' BIT-IDENTICAL'
    else
        write(*, '(a)') ' MISMATCH'
        error stop 1
    end if

contains

    subroutine fresh(nn, l, h, density)
        integer, intent(in) :: nn, l, h
        real(kind = dbl), intent(in) :: density
        real(kind = dbl), allocatable :: m(:)
        if (allocated(x)) deallocate(x, y, a, b)
        allocate(x(nn), y(nn), a(l:h), b(l:h), m(nn))
        call random_number(x)
        call random_number(y)
        x = (x - 0.5d0) * 2.1d0
        y = y * 3d-2 + 412.7d0
        call random_number(m)
        where (m < density) y = err
        call random_number(m)
        where (m < density / 10d0) x = err
    end subroutine fresh

    subroutine check(nn, l, h, density)
        integer, intent(in) :: nn, l, h
        real(kind = dbl), intent(in) :: density
        integer :: k, bad
        call fresh(nn, l, h, density)
        a = -1d300
        b = 1d300
        call FrozenSeries(x, y, nn, l, h, err, a)
        call CrossCovarianceSeries(x, y, nn, l, h, err, b)
        bad = 0
        do k = l, h
            if (transfer(a(k), 0_int64) /= transfer(b(k), 0_int64)) bad = bad + 1
        end do
        if (bad > 0 .and. mismatched == 0) &
            write(*, '(a,i0,a,i0,a,i0,a,f6.3,a,i0)') ' first mismatch: n=', nn, &
                ' window ', l, ':', h, ' density ', density, ', elements ', bad
        mismatched = mismatched + bad
        cases = cases + 1
    end subroutine check

end program covmax_bitexact_main
