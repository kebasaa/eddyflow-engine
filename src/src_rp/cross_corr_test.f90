!***************************************************************************
! cross_corr_test.f90
! -------------------
! Copyright © 2026-    , ETH Zurich, Jonathan Muller
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
! \brief       CCF ratio test with/without repeated values
! \author      Jonathan Muller, ETH Zurich
!
!***************************************************************************
! Reference: Vitale, D. et al. (2020). Biogeosciences Discussion.
subroutine CrossCorrTest(Set, nrow, ncol)
    use m_rp_global_var
    implicit none
    integer, intent(in) :: nrow, ncol
    real(kind = dbl), intent(in) :: Set(nrow, ncol)
    integer, parameter :: lagmin = -25
    integer, parameter :: lagmax = 25
    integer :: icol, irec
    real(kind = dbl) :: dedup_set(nrow, ncol)
    real(kind = dbl) :: raw_ccf(lagmin:lagmax), dup_ccf(lagmin:lagmax)
    real(kind = dbl) :: cov_val, sig_raw(1), sig_dup(1)
    real(kind = dbl), external :: LaggedCovarianceNoError

    write(*, '(a)', advance = 'no') &
        '  Evaluating R2 on CCFs with and without repeated values..'

    dedup_set = Set
    do icol = u, gas4
        if (OutVarPresent(icol)) then
            do irec = 2, nrow
                if (dabs(Set(irec, icol) - Set(irec-1, icol)) < 1d-8) then
                    dedup_set(irec, icol) = error
                end if
            end do
        end if
    end do

    do icol = ts, gas4
        if (OutVarPresent(icol)) then
            call CrossCorrelation(Set(:, w), Set(:, icol), nrow, lagmin, lagmax, raw_ccf)
            call CrossCorrelation(dedup_set(:, w), dedup_set(:, icol), &
                                  nrow, lagmin, lagmax, dup_ccf)
            cov_val = LaggedCovarianceNoError( &
                raw_ccf(lagmin), dup_ccf(lagmin), lagmax-lagmin+1, 0, error)
            call StDevNoError(raw_ccf(lagmin), lagmax-lagmin+1, 1, sig_raw, error)
            call StDevNoError(dup_ccf(lagmin), lagmax-lagmin+1, 1, sig_dup, error)
        end if
    end do

    write(*, *) ' Done.'
end subroutine CrossCorrTest
