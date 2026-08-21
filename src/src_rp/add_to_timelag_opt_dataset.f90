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
    !> local variables
    integer :: gas
    integer :: wsl
    real(kind = dbl) :: flux_test
    character(32) :: label
    include '../src_common/interfaces.inc'

    wsl = PrimaryWaterOutSlot()

    !> Passive gases. Three unrolled arms before, naming co2, ch4 and the
    !> fourth slot, so a fifth gas never contributed to the time-lag
    !> optimisation dataset however strong its flux.
    !>
    !> NOTE: the CO2 arm compared dabs(flux) against the threshold and the
    !> other two compared the signed flux. That is preserved rather than
    !> unified - CO2 is the one gas routinely measured with a negative flux,
    !> so on the others a signed comparison and a magnitude comparison agree
    !> in practice, but making them agree by fiat would move numbers. Worth
    !> deciding deliberately.
    do gas = firstGas, lastGas
        if (gas == wsl) cycle
        TimelagOpt(n)%tlag(gas) = error
        if (.not. E2Col(gas)%present) cycle

        label = GasOutputLabel(gas)
        call lowercase(label)
        if (trim(adjustl(label)) == 'co2') then
            flux_test = dabs(Flux0%gas(gas))
        else
            flux_test = Flux0%gas(gas)
        end if

        if (flux_test > TOSetup%gas_min_flux(gas) &
            .and. Essentials%used_timelag(gas) /= E2Col(gas)%max_tl &
            .and. Essentials%used_timelag(gas) /= E2Col(gas)%min_tl) &
            TimelagOpt(n)%tlag(gas) = Essentials%used_timelag(gas)
    end do

    !> Water vapor and RH. Gated on the site's latent heat flux rather than on
    !> its own flux, which is why it keeps an arm of its own.
    if (E2Col(wsl)%present) then
        if (Flux0%LE > TOSetup%le_min_flux &
            .and. Essentials%used_timelag(wsl) /= E2Col(wsl)%max_tl &
            .and. Essentials%used_timelag(wsl) /= E2Col(wsl)%min_tl) then
            TimelagOpt(n)%tlag(wsl) = Essentials%used_timelag(wsl)
        else
            TimelagOpt(n)%tlag(wsl) = error
        end if
        if (Stats%RH >= 0d0 .and. Stats%RH <= 100d0) then
            TimelagOpt(n)%RH = Stats%RH
        else
            TimelagOpt(n)%RH = error
        end if
    else
        TimelagOpt(n)%tlag(wsl) = error
        TimelagOpt(n)%RH = error
    end if
end subroutine AddToTimelagOptDataset
