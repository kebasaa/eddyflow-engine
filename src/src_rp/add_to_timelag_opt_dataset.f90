!***************************************************************************
! add_to_timelag_opt_dataset.f90
! ------------------------------
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
! \brief       Store calculated time-lags and other variables used for the
!              time-lag optimization, if all conditions are met
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine AddToTimelagOptDataset(TimelagOpt, nrow, n)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow
    integer, intent(in) :: n
    type(TimeLagOptType), intent(inout):: TimelagOpt(nrow)


    !> Passive gases
    if (E2Col(co2)%present &
        .and. dabs(Flux0%co2) > TOSetup%co2_min_flux &
        .and. Essentials%used_timelag(co2) /= E2Col(co2)%max_tl &
        .and. Essentials%used_timelag(co2) /= E2Col(co2)%min_tl) then
            TimelagOpt(n)%tlag(co2) = Essentials%used_timelag(co2)
    else
        TimelagOpt(n)%tlag(co2) = error
    end if

    if (E2Col(ch4)%present &
        .and. Flux0%ch4 > TOSetup%ch4_min_flux &
        .and. Essentials%used_timelag(ch4) /= E2Col(ch4)%max_tl &
        .and. Essentials%used_timelag(ch4) /= E2Col(ch4)%min_tl) then
        TimelagOpt(n)%tlag(ch4) = Essentials%used_timelag(ch4)
    else
        TimelagOpt(n)%tlag(ch4) = error
    end if

    if (E2Col(gas4)%present &
        .and. Flux0%gas4 > TOSetup%gas4_min_flux &
        .and. Essentials%used_timelag(gas4) /= E2Col(gas4)%max_tl &
        .and. Essentials%used_timelag(gas4) /= E2Col(gas4)%min_tl) then
        TimelagOpt(n)%tlag(gas4) = Essentials%used_timelag(gas4)
    else
        TimelagOpt(n)%tlag(gas4) = error
    end if

    !> Water vapor and RH
    if (E2Col(h2o)%present) then
        if (Flux0%LE > TOSetup%le_min_flux &
            .and. Essentials%used_timelag(h2o) /= E2Col(h2o)%max_tl &
            .and. Essentials%used_timelag(h2o) /= E2Col(h2o)%min_tl) then
            TimelagOpt(n)%tlag(h2o) = Essentials%used_timelag(h2o)
        else
            TimelagOpt(n)%tlag(h2o) = error
        end if
        if (Stats%RH >= 0d0 .and. Stats%RH <= 100d0) then
            TimelagOpt(n)%RH = Stats%RH
        else
            TimelagOpt(n)%RH = error
        end if
    else
        TimelagOpt(n)%tlag(h2o) = error
        TimelagOpt(n)%RH = error
    end if
end subroutine AddToTimelagOptDataset
