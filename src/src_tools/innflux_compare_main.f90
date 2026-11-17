!***************************************************************************
! innflux_compare_main.f90
! ------------------------
! Copyright (C) 2026, ETH Zurich, Jonathan Muller
!
! This file is part of EddyFlow.
!
! EddyFlow (TM) is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! EddyFlow (TM) is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with EddyFlow (TM).  If not, see <http://www.gnu.org/licenses/>.
!
!***************************************************************************
!
! \brief       Comparison driver for tests/innflux/compare.py: runs the
!              engine's innFLUX-derived routines (Striednig et al. 2020,
!              innFLUX v1.2.0) on cases the Python side writes, so the port
!              of the MATLAB code and the engine can be compared number for
!              number.
!
!              Reads a case file named on the command line and writes one
!              result line per case to standard output. Case kinds:
!
!              peak    lo hi slo shi smooth_len sign_mode debaseline
!                      then cov(lo..hi), one value per line.
!                      -> "peak <lag> <found>"
!              series  n slo shi smooth_len sign_mode debaseline
!                      then n lines "w c". The covariance is computed by
!                      CrossCovarianceSeries over the window the peak search
!                      needs, exactly as CovMax does.
!                      -> "series <lag> <found> <cov at lag>"
!
!              Lags in rows, EddyFlow's sign: positive = scalar delayed.
!              Error code -9999. Built by the CMake target innfluxcmp.
! \author      Jonathan Muller
!***************************************************************************
program innflux_compare
    use m_numeric_kinds
    use m_covmax_core, only: CrossCovarianceSeries
    use m_lagpeak_core, only: FindCovPeak, PeakPadding
    implicit none
    real(kind = dbl), parameter :: err = -9999d0
    character(len = 512) :: path
    character(len = 16) :: kind
    integer :: u
    integer :: ios
    integer :: lo
    integer :: hi
    integer :: slo
    integer :: shi
    integer :: smooth_len
    integer :: sign_mode
    integer :: idbl
    integer :: n
    integer :: k
    integer :: lag
    integer :: pad
    logical :: found
    real(kind = dbl), allocatable :: cov(:)
    real(kind = dbl), allocatable :: w(:)
    real(kind = dbl), allocatable :: c(:)

    call get_command_argument(1, path)
    if (len_trim(path) == 0) then
        write(*, '(a)') 'usage: innflux_compare <cases.txt>'
        stop 2
    end if
    open(newunit = u, file = trim(path), status = 'old', action = 'read', iostat = ios)
    if (ios /= 0) then
        write(*, '(a)') 'cannot open ' // trim(path)
        stop 2
    end if

    do
        read(u, *, iostat = ios) kind
        if (ios /= 0) exit
        select case (trim(kind))
            case ('peak')
                read(u, *) lo, hi, slo, shi, smooth_len, sign_mode, idbl
                allocate(cov(lo:hi))
                do k = lo, hi
                    read(u, *) cov(k)
                end do
                call FindCovPeak(cov, lo, hi, slo, shi, smooth_len, sign_mode, &
                    idbl /= 0, err, lag, found)
                write(*, '(a, i0, 1x, l1)') 'peak ', lag, found
                deallocate(cov)
            case ('series')
                read(u, *) n, slo, shi, smooth_len, sign_mode, idbl
                allocate(w(n), c(n))
                do k = 1, n
                    read(u, *) w(k), c(k)
                end do
                pad = PeakPadding(smooth_len, sign_mode)
                lo = max(slo - pad, -(n - 1))
                hi = min(shi + pad, n - 1)
                allocate(cov(lo:hi))
                call CrossCovarianceSeries(w, c, n, lo, hi, err, cov)
                call FindCovPeak(cov, lo, hi, slo, shi, smooth_len, sign_mode, &
                    idbl /= 0, err, lag, found)
                write(*, '(a, i0, 1x, l1, 1x, es24.16)') 'series ', lag, found, cov(lag)
                deallocate(w, c, cov)
            case default
                write(*, '(a)') 'unknown case kind ' // trim(kind)
                stop 2
        end select
    end do
    close(u)
end program innflux_compare
