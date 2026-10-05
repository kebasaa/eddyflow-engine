!***************************************************************************
! m_covmax_core.f90
! -----------------
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
! \brief       The lagged cross-covariance covariance maximisation searches,
!              computed without building each shifted pair.
!
!              CovMax is what time-lag optimisation runs for every gas of every
!              period, and what PWB runs for its terminal fallback. For each lag
!              it used to allocate two arrays, copy the shifted pair into one
!              and that into the other, and hand the pair to
!              CovarianceMatrixNoError - which computes the full 2x2 matrix,
!              four passes over the data, of which only Cov(1,2) was read.
!
!              This computes Cov(1,2) alone, reading the shifted elements where
!              they already are. It performs exactly the operations the matrix
!              routine performs for that entry: the same pairwise test against
!              the error code, the same products, sums and count accumulated in
!              the same order, and the same C/N - (si/N)(sj/N) at the end. So
!              every value is bitwise the one it replaces - which
!              src_tools/covmax_bitexact_main.f90 checks against a frozen copy
!              of the old path, sparse columns included.
!
! \note        Pure, and a function of its arguments only, so it can be checked
!              outside the engine like m_pwb_core.
! \author      Jonathan Muller
! \sa          timelag_handle.f90, stats_operator_no_error.f90
!***************************************************************************
module m_covmax_core
    use m_numeric_kinds
    implicit none
    private

    public :: CrossCovarianceSeries

contains

    !***************************************************************************
    !> \brief Cov(1,2) of the pair (col1, col2) shifted by every lag in
    !>        [lagmin, lagmax], counting only pairs where neither is err.
    !>
    !> A negative lag pairs col1(k - lag) with col2(k); a non-negative one
    !> pairs col1(k) with col2(k + lag); k runs from 1 to n - |lag|. A lag with
    !> no valid pair - including one as wide as the series - gives err.
    !***************************************************************************
    pure subroutine CrossCovarianceSeries(col1, col2, n, lagmin, lagmax, err, cov)
        integer, intent(in) :: n
        integer, intent(in) :: lagmin
        integer, intent(in) :: lagmax
        real(kind = dbl), intent(in) :: col1(n)
        real(kind = dbl), intent(in) :: col2(n)
        real(kind = dbl), intent(in) :: err
        real(kind = dbl), intent(out) :: cov(lagmin:lagmax)
        integer :: lag
        integer :: k
        integer :: n2
        integer :: nact
        real(kind = dbl) :: a
        real(kind = dbl) :: b
        real(kind = dbl) :: c
        real(kind = dbl) :: si
        real(kind = dbl) :: sj

        do lag = lagmin, lagmax
            n2 = n - abs(lag)
            c = 0d0
            si = 0d0
            sj = 0d0
            nact = 0
            if (lag < 0) then
                do k = 1, n2
                    a = col1(k - lag)
                    b = col2(k)
                    if (a /= err .and. b /= err) then
                        nact = nact + 1
                        c = c + a * b
                        si = si + a
                        sj = sj + b
                    end if
                end do
            else
                do k = 1, n2
                    a = col1(k)
                    b = col2(k + lag)
                    if (a /= err .and. b /= err) then
                        nact = nact + 1
                        c = c + a * b
                        si = si + a
                        sj = sj + b
                    end if
                end do
            end if
            !> The matrix routine's normalisation, step for step.
            if (nact /= 0) then
                si = si / dble(nact)
                sj = sj / dble(nact)
                c = c / dble(nact)
                c = c - si * sj
            else
                c = err
            end if
            cov(lag) = c
        end do
    end subroutine CrossCovarianceSeries

end module m_covmax_core
