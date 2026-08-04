!***************************************************************************
! fix_timelag_opt_dataset.f90
! ---------------------------
! Copyright © 2007-2011, Eco2s team, Gerardo Fratini
! Copyright © 2011-2026, LI-COR Biosciences, Gerardo Fratini
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
! \brief       Eliminate error codes for easier following processing
!              Needs to create a new RH column for each gas
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine FixTimelagOptDataset(TimelagOpt, nrow, toSet, ton, actn, tlncol)
    use m_rp_global_var
    implicit none
    !> In/out variables
    integer, intent(in) :: nrow
    integer, intent(in) :: ton
    integer, intent(in) :: tlncol
    type(TimeLagOptType), intent(in):: TimelagOpt(nrow)
    type(TimeLagDatasetType), intent(out):: toSet(ton)
    integer, intent(out) :: actn(tlncol)
    !> Local variables
    integer :: i
    integer :: gas
    logical, external :: GasSlotIsWater

    !> Every configured gas, and water identified by its record.
    !>
    !> Bounded at the fourth slot, actn stayed zero for every gas past it, and
    !> actn is what OptimizeTimelags reads as its sample count and what
    !> WriteOutTimelagOptimization gates its rows on - so the whole time-lag
    !> optimisation was inert for those gases while AddToTimelagOptDataset,
    !> already widened, went on feeding them in.
    !>
    !> Water is the one classed by relative humidity, so it is skipped when RH
    !> is missing and carries RH into the dataset when it is not. Asked as
    !> `gas == h2o` that named record two; a project whose water sits elsewhere
    !> both dropped the real hygrometer's RH and applied the RH treatment to
    !> whatever gas held slot six.
    toSet = TimelagDatasetType(0d0, 0d0)
    actn = 0
    do i = 1, ton
        do gas = firstGas, lastGas
            if (gas - firstGas + 1 > min(EddyFlowProj%gas_num, MaxNumGases)) exit
            if(GasSlotIsWater(gas) .and. TimelagOpt(i)%RH == error) cycle
            if(TimelagOpt(i)%tlag(gas) /= error) then
                actn(gas) = actn(gas) + 1
                toSet(actn(gas))%tlag(gas) = TimelagOpt(i)%tlag(gas)
                if (GasSlotIsWater(gas)) toSet(actn(gas))%RH = TimelagOpt(i)%RH
            end if
        end do
    end do
end subroutine FixTimelagOptDataset

