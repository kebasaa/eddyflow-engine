!***************************************************************************
! filter_dataset_for_physical_thresholds.f90
! ------------------------------------------
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
! \brief       Absolute-limits filter for gas species (vectorised)
! \author      Jonathan Muller, ETH Zurich
!
!***************************************************************************
subroutine FilterDatasetForPhysicalThresholds(Set, N, M, FilterWhat)
    use m_rp_global_var
    integer, intent(in) :: N, M
    logical, intent(in) :: FilterWhat(M)
    real(kind = dbl), intent(inout) :: Set(N, M)
    integer :: gas
    real(kind = dbl) :: dens_scale
    logical, external :: GasSlotIsWater

    !> One pass per configured gas. This was four hand-copied blocks naming
    !> co2/h2o/ch4/gas4, differing in exactly one thing - the molar-density
    !> scale, which is 1d3 for every species but water. That is a property of
    !> the species and not of the slot, the same correction Fluxes0 needed,
    !> so it is asked of the record.
    !>
    !> The guard matters more than the bound. This filter and
    !> test_absolute_limits consult the *same* al%gas_min/gas_max pair, and
    !> only the four historical slots get them from fixed project keys; past
    !> those they come from the per-gas records, and a project naming a gas
    !> without them leaves the pair at whatever memory held. Testing against
    !> that rejects every value and replaces the entire series with the error
    !> code - the incident this repository already recorded once. The test
    !> declines in exactly this case, so the filter must decline on the same
    !> condition or the two disagree about the same gas.
    do gas = firstGas, lastGas
        if (.not. E2Col(gas)%present .or. .not. FilterWhat(gas)) cycle
        if (al%gas_max(gas) <= al%gas_min(gas)) cycle

        if (E2Col(gas)%measure_type == 'molar_density') then
            !> Water is held on the mmol basis, every other species on umol.
            if (GasSlotIsWater(gas)) then
                dens_scale = Ambient%Va
            else
                dens_scale = Ambient%Va * 1d3
            end if
            where (Set(:,gas) /= error .and. &
                   (Set(:,gas)*dens_scale < al%gas_min(gas) .or. &
                    Set(:,gas)*dens_scale > al%gas_max(gas)))
                Set(:,gas) = error
            end where
        else
            where (Set(:,gas) /= error .and. &
                   (Set(:,gas) < al%gas_min(gas) .or. Set(:,gas) > al%gas_max(gas)))
                Set(:,gas) = error
            end where
        end if
    end do
end subroutine FilterDatasetForPhysicalThresholds
