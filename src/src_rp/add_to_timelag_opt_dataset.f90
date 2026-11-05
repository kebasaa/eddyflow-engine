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
    logical :: fluxOk
    real(kind = dbl) :: le
    real(kind = dbl) :: rh
    real(kind = dbl), external :: PeriodWaterRH
    include '../src_common/interfaces.inc'

    !> Zero with no hygrometer, so no gas is taken for water: the fallback
    !> slot six holds another gas then, which this kept out of the gases'
    !> dataset and filed as water.
    wsl = PrimaryWaterSlot()
    TimelagOpt(n)%RH = error

    !> What the preliminary time-lag assessment's table shows: every lag the
    !> detector found, its window, and the two verdicts the gates below
    !> reach, so a pile-up at a window edge is shown rather than dropped.
    TimelagOpt(n)%RawLag = error
    TimelagOpt(n)%WinMin = error
    TimelagOpt(n)%WinMax = error
    TimelagOpt(n)%Flux = error
    TimelagOpt(n)%AtEdge = .false.
    TimelagOpt(n)%FluxOk = .false.

    !> Passive gases. Three unrolled arms before, naming co2, ch4 and the
    !> fourth slot, so a fifth gas never contributed to the time-lag
    !> optimisation dataset however strong its flux.
    !>
    !> Every gas is gated on the magnitude of its flux. The threshold is
    !> there because a weak covariance gives a poorly defined correlation
    !> peak, and the sign does not change how well it is defined. Only a
    !> column labelled exactly co2 used to be compared by magnitude; any
    !> other gas was compared signed, so an uptake gas such as COS fed the
    !> optimisation nothing but its occasional emission half-hours. A
    !> missing flux is excluded outright - dabs(error) clears any threshold.
    !> Water keeps the signed latent heat test below.
    do gas = firstGas, lastGas
        if (gas == wsl) cycle
        if (GasSlotIsWater(gas)) cycle
        TimelagOpt(n)%tlag(gas) = error
        if (.not. E2Col(gas)%present) cycle

        fluxOk = Flux0%gas(gas) /= error &
            .and. dabs(Flux0%gas(gas)) > TOSetup%gas_min_flux(gas)

        !> Nested: OffWindowEdge is a function, and in one .and. it need
        !> not be evaluated.
        if (fluxOk) then
            if (OffWindowEdge(gas)) TimelagOpt(n)%tlag(gas) = Essentials%used_timelag(gas)
        end if
        call NoteRawLag(gas, Flux0%gas(gas), fluxOk)
    end do

    !> Water vapour and RH, every hygrometer. Each is gated on its own latent
    !> heat flux - the threshold is a latent heat flux, not a molar one,
    !> which is why water keeps an arm of its own - and carries the RH its
    !> lag is classed by. A second hygrometer used to go through the passive
    !> arm above, gated on its molar flux against a threshold meant for other
    !> species, and took the designated hygrometer's RH.
    !>
    !> The designated one is gated on Flux0%LE itself, so its determinations
    !> are exactly what they were; the others' latent heat is formed the way
    !> Fluxes0_rp forms that one.
    do gas = firstGas, lastGas
        if (gas /= wsl) then
            if (.not. GasSlotIsWater(gas)) cycle
        end if
        TimelagOpt(n)%tlag(gas) = error
        if (.not. E2Col(gas)%present) cycle
        if (gas == wsl) then
            le = Flux0%LE
        elseif (Flux0%gas(gas) /= error .and. Ambient%lambda > 0d0) then
            le = Flux0%gas(gas) * Ambient%lambda * MW_H2O * 1d-3
        else
            le = error
        end if
        if (le /= error .and. le > TOSetup%le_min_flux) then
            if (OffWindowEdge(gas)) TimelagOpt(n)%tlag(gas) = Essentials%used_timelag(gas)
        end if
        call NoteRawLag(gas, le, le /= error .and. le > TOSetup%le_min_flux)
        rh = PeriodWaterRH(gas)
        if (rh >= 0d0 .and. rh <= 100d0) TimelagOpt(n)%RH(gas) = rh
    end do

contains

    !> Whether a gas's used lag is clear of its search window's edges, to
    !> within half a sample. A lag is a whole number of samples and a window
    !> need not be - most derived ones, nominal - 2 s or a separation-based
    !> bound, are not - so the exact comparison this replaces almost never
    !> recognised a lag pinned at the edge, and those entered the dataset
    !> the optimised lag and window are fitted to.
    logical function OffWindowEdge(gas)
        integer, intent(in) :: gas
        real(kind = dbl) :: lag
        real(kind = dbl) :: halfSample

        OffWindowEdge = .false.
        lag = Essentials%used_timelag(gas)
        if (lag == error) return
        halfSample = 0d0
        if (Metadata%ac_freq > 0d0) halfSample = 0.5d0 / Metadata%ac_freq
        OffWindowEdge = dabs(lag - E2Col(gas)%min_tl) > halfSample &
            .and. dabs(lag - E2Col(gas)%max_tl) > halfSample
    end function OffWindowEdge

    !> The raw lag of one gas, with its window and verdicts. The edge test
    !> allows half a sample, as OffWindowEdge does.
    subroutine NoteRawLag(gas, flux, fluxOk)
        integer, intent(in) :: gas
        real(kind = dbl), intent(in) :: flux
        logical, intent(in) :: fluxOk
        real(kind = dbl) :: lag
        real(kind = dbl) :: halfSample

        lag = Essentials%actual_timelag(gas)
        if (lag == error) return
        halfSample = 0d0
        if (Metadata%ac_freq > 0d0) halfSample = 0.5d0 / Metadata%ac_freq
        TimelagOpt(n)%RawLag(gas) = lag
        TimelagOpt(n)%WinMin(gas) = E2Col(gas)%min_tl
        TimelagOpt(n)%WinMax(gas) = E2Col(gas)%max_tl
        TimelagOpt(n)%Flux(gas) = flux
        TimelagOpt(n)%FluxOk(gas) = fluxOk
        TimelagOpt(n)%AtEdge(gas) = &
            dabs(lag - E2Col(gas)%min_tl) <= halfSample &
            .or. dabs(lag - E2Col(gas)%max_tl) <= halfSample
    end subroutine NoteRawLag
end subroutine AddToTimelagOptDataset
