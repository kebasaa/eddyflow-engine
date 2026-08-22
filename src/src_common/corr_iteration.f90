!***************************************************************************
! corr_iteration.f90
! ------------------
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
! \brief       How far the gas fluxes moved between two passes of the
!              iterative correction, as a percentage.
! \author      Jonathan Muller
! \note
!              The convergence measure for the iterative correction, and the
!              only thing the early exit tests. Shared by both applications
!              because the loop exists in both and they must agree on what
!              "converged" means.
!
!              THE WORST GAS, NOT EACH ONE. Every gas is corrected at the
!              same z/L, so what the loop actually converges is the
!              stability, and the per-gas changes move together. Reporting
!              the largest is the useful summary: it is the one that decides
!              whether the period converged, and a column per gas would be
!              near-copies of it.
!
!              EddyUH reports the same quantity per variable, as covsvar,
!              and tests none of it - its loop runs a fixed number of passes
!              (EddyUH.m:722-903). Here it is both reported and, when a
!              tolerance is stated, acted on.
!
!              RELATIVE TO THE PREVIOUS PASS. A flux that is error-coded in
!              either pass is skipped rather than counted as an infinite
!              change, and one that was exactly zero before is skipped too -
!              there is no relative change from zero, and treating it as
!              total would make a period that has nothing to converge look
!              like the worst one in the run.
!***************************************************************************
real(kind = dbl) function WorstRelativeChange(before, after)
    use m_common_global_var
    implicit none
    !> in/out variables
    real(kind = dbl), intent(in) :: before(GHGNumVar)
    real(kind = dbl), intent(in) :: after(GHGNumVar)
    !> local variables
    integer :: gas
    real(kind = dbl) :: dev

    WorstRelativeChange = error
    do gas = firstGas, lastGas
        if (.not. E2Col(gas)%present) cycle
        if (before(gas) == error .or. after(gas) == error) cycle
        if (before(gas) == 0d0) cycle
        dev = dabs((after(gas) - before(gas)) / before(gas)) * 100d0
        if (WorstRelativeChange == error) then
            WorstRelativeChange = dev
        else if (dev > WorstRelativeChange) then
            WorstRelativeChange = dev
        end if
    end do
end function WorstRelativeChange
