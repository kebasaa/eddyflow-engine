!***************************************************************************
! adjust_timelag_opt_settings.f90
! -------------------------------
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
! \brief       Adjust time-lag opt settings if user did not set or set improperly
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine AdjustTimelagOptSettings()
    use m_rp_global_var
    implicit none
    !> local variables
    integer :: gas
    real(kind = dbl) :: nominal
    real(kind = dbl) :: mult(GHGNumVar)
    real(kind = dbl) :: tube_time(GHGNumVar)
    real(kind = dbl) :: tube_volume(GHGNumVar)
    real(kind = dbl) :: cell_time(GHGNumVar)
    real(kind = dbl) :: cell_volume(GHGNumVar)
    real(kind = dbl) :: safety
    logical, external :: GasSlotIsWater


    !> Initialization to zero of all timelags
    cell_volume(:) = 0d0
    E2Col(:)%min_tl = 0d0
    E2Col(:)%max_tl = 0d0

    !> Initialize multiplier
    !>
    !> Water is the active gas - it adsorbs on the tube wall, so its lag can run
    !> far past the transit time and the search window has to be widened for it.
    !> That is a property of the species, so it is asked of the record. Keyed on
    !> slot six it widened whatever record two held and left a hygrometer
    !> declared anywhere else with the passive window, roughly a third as wide,
    !> so the true water lag could fall outside the range ever searched - and a
    !> lag that is never searched is silently reported as the default one.
    mult(:) = 2d0     !< For passive gases
    do gas = firstGas, lastGas
        if (GasSlotIsWater(gas)) mult(gas) = 10d0  !< For active gases
    end do
    safety = 0.3d0    !< Safety margin for min/max setting

    !> Transit time in cell and sampling lines of closed path instruments
    where (E2Col(firstGas:lastGas)%instr%path_type == 'closed')
        tube_volume(firstGas:lastGas) = &
            (p * (E2Col(firstGas:lastGas)%instr%tube_d / 2d0)**2 * &
            E2Col(firstGas:lastGas)%instr%tube_l)
        tube_time(firstGas:lastGas) =  tube_volume(firstGas:lastGas) &
            / E2Col(firstGas:lastGas)%instr%tube_f

        cell_volume(firstGas:lastGas) = &
            (p * (E2Col(firstGas:lastGas)%instr%hpath_length / 2d0)**2 * &
                                E2Col(firstGas:lastGas)%instr%vpath_length)
        cell_time(firstGas:lastGas) = cell_volume(firstGas:lastGas) &
            / E2Col(firstGas:lastGas)%instr%tube_f
    elsewhere
        tube_time(firstGas:lastGas) = 0d0
        cell_time(firstGas:lastGas) = 0d0
    end where

    !> If user didn't set min and max time-lags, does so by using tube properties for closed path
    !> and distances for open path
    do gas = firstGas, lastGas
        if (E2Col(gas)%present) then
            if (TOSetup%min_lag(gas) < TlagDeriveThreshold &
                .or. TOSetup%max_lag(gas) < TlagDeriveThreshold) then
                if (E2Col(gas)%instr%path_type == 'closed') then
                    !> Closed path
                    nominal = tube_time(gas) + cell_time(gas)
                    E2Col(gas)%min_tl = max(0d0, nominal - 2d0)
                    E2Col(gas)%max_tl = min(nominal + mult(gas) * nominal, &
                        RPsetup%avrg_len * 60d0) + safety
                else
                    !> Open path
                    E2Col(gas)%min_tl = - dsqrt(E2Col(gas)%instr%hsep**2 &
                        + E2Col(gas)%instr%vsep**2) * 2d0 - safety
                    E2Col(gas)%max_tl = + dsqrt(E2Col(gas)%instr%hsep**2 &
                        + E2Col(gas)%instr%vsep**2) * 2d0 + safety
                end if
            else
                E2Col(gas)%min_tl = TOSetup%min_lag(gas)
                E2Col(gas)%max_tl = TOSetup%max_lag(gas)
            end if
        end if
    end do
end subroutine AdjustTimelagOptSettings
