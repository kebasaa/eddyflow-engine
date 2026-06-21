!***************************************************************************
! fisher.f90
! ----------
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
! \brief       Correlation-matrix difference test (Fisher)
! \author      Jonathan Muller, ETH Zurich
!
!***************************************************************************
subroutine fisher(Set, nrow, ncol)
    use m_rp_global_var
    implicit none
    integer, intent(in) :: nrow, ncol
    real(kind = dbl), intent(in) :: Set(nrow, ncol)
    integer :: ci, cj, rec
    real(kind = dbl) :: full_corr(ncol, ncol)
    real(kind = dbl) :: nodup_set(nrow, ncol)
    real(kind = dbl) :: nodup_corr(ncol, ncol)

    write(*, '(a)', advance = 'no') &
        '  Computing correlation-matrix difference with/without repeated values..'

    nodup_set = Set
    do ci = u, gas4
        do rec = 2, nrow
            if (dabs(Set(rec, ci) - Set(rec-1, ci)) < 1d-8) &
                nodup_set(rec, ci) = error
        end do
    end do

    call CorrelationMatrixNoError(Set, nrow, ncol, full_corr, error)
    call CorrelationMatrixNoError(nodup_set, nrow, ncol, nodup_corr, error)

    do ci = u, gas4
        do cj = u, gas4
            if (full_corr(ci,cj) /= error .and. nodup_corr(ci,cj) /= error) then
                Essentials%CorrDiff(ci, cj) = dabs(full_corr(ci,cj) - nodup_corr(ci,cj))
            else
                Essentials%CorrDiff(ci, cj) = error
            end if
        end do
    end do

    write(*, *) ' Done.'
end subroutine fisher
