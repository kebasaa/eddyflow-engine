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

    if (E2Col(co2)%present .and. FilterWhat(co2)) then
        if (E2Col(co2)%measure_type == 'molar_density') then
            where (Set(:,co2) /= error .and. &
                   (Set(:,co2)*Ambient%Va*1d3 < al%co2_min .or. &
                    Set(:,co2)*Ambient%Va*1d3 > al%co2_max))
                Set(:,co2) = error
            end where
        else
            where (Set(:,co2) /= error .and. &
                   (Set(:,co2) < al%co2_min .or. Set(:,co2) > al%co2_max))
                Set(:,co2) = error
            end where
        end if
    end if

    if (E2Col(h2o)%present .and. FilterWhat(h2o)) then
        if (E2Col(h2o)%measure_type == 'molar_density') then
            where (Set(:,h2o) /= error .and. &
                   (Set(:,h2o)*Ambient%Va < al%h2o_min .or. &
                    Set(:,h2o)*Ambient%Va > al%h2o_max))
                Set(:,h2o) = error
            end where
        else
            where (Set(:,h2o) /= error .and. &
                   (Set(:,h2o) < al%h2o_min .or. Set(:,h2o) > al%h2o_max))
                Set(:,h2o) = error
            end where
        end if
    end if

    if (E2Col(ch4)%present .and. FilterWhat(ch4)) then
        if (E2Col(ch4)%measure_type == 'molar_density') then
            where (Set(:,ch4) /= error .and. &
                   (Set(:,ch4)*Ambient%Va*1d3 < al%ch4_min .or. &
                    Set(:,ch4)*Ambient%Va*1d3 > al%ch4_max))
                Set(:,ch4) = error
            end where
        else
            where (Set(:,ch4) /= error .and. &
                   (Set(:,ch4) < al%ch4_min .or. Set(:,ch4) > al%ch4_max))
                Set(:,ch4) = error
            end where
        end if
    end if

    if (E2Col(gas4)%present .and. FilterWhat(gas4)) then
        if (E2Col(gas4)%measure_type == 'molar_density') then
            where (Set(:,gas4) /= error .and. &
                   (Set(:,gas4)*Ambient%Va*1d3 < al%gas4_min .or. &
                    Set(:,gas4)*Ambient%Va*1d3 > al%gas4_max))
                Set(:,gas4) = error
            end where
        else
            where (Set(:,gas4) /= error .and. &
                   (Set(:,gas4) < al%gas4_min .or. Set(:,gas4) > al%gas4_max))
                Set(:,gas4) = error
            end where
        end if
    end if
end subroutine FilterDatasetForPhysicalThresholds
