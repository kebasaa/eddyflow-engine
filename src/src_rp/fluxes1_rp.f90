!***************************************************************************
! fluxes1_rp.f90
! --------------
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
! \brief       Calculates fluxes at Level 1. Mainly spectral corrections
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo        Merge with corresponding FCC sub into a common one
!***************************************************************************
subroutine Fluxes1_rp()
    use m_rp_global_var
    implicit none
    real(kind = dbl)  :: Cox
    integer :: msl
    integer :: wsl
    integer :: wsx
    include '../src_common/interfaces_1.inc'

    call LogSayNoAdv('  Calculating fluxes Level 1..')

    Flux1 = errFlux

    !> Water's own slot. The oxygen correction below is a krypton /
    !> Lyman-alpha hygrometer correction - genuinely about water - and the
    !> E/ET/LE terms are the primary water's. Both resolve the slot rather
    !> than assuming record two holds water.
    wsl = PrimaryWaterSlot()
    !> An always-in-bounds stand-in for wsl inside guard expressions.
    !>
    !> Fortran does not mandate short-circuit `.and.`, so
    !> `wsl >= firstGas .and. X(wsl)` still evaluates X(wsl) - and with no
    !> hygrometer wsl is 0, which is out of bounds. The wsl >= firstGas test
    !> still decides the outcome; wsx only keeps the subscript legal while it
    !> is being decided.
    wsx = max(wsl, firstGas)

    !> First, apply oxygen correction to Krypton and Lyman-alpha hygrometers,
    !> according to van Dijk et al. (2003, JAOT, eq. 13b).
    !>
    !> Over every hygrometer, each with its own ko/kw. This is genuinely a
    !> water-vapour correction - oxygen absorbs in the same band the
    !> instrument uses to measure water - so the H2O assumption stays; only
    !> the slot is resolved. Applied to the primary water alone, a second
    !> krypton on the same site was never corrected at all.
    do msl = firstGas, lastGas
        if (.not. GasSlotIsWater(msl)) cycle
        if (.not. E2Col(msl)%present) cycle
        select case (E2Col(msl)%Instr%model(1:len_trim(E2Col(msl)%Instr%model) - 2))
            case('open_path_krypton','closed_path_krypton', &
                    'open_path_lyman','closed_path_lyman')
                !> That hygrometer's own extinction coefficients. Absent, the
                !> correction is not performed for it - never performed with
                !> another instrument's numbers.
                if (E2Col(msl)%Instr%ko /= error .and. E2Col(msl)%Instr%kw /= 0d0 &
                    .and. Ambient%Ta > 0d0 .and. Ambient%Bowen /= error &
                    .and. Ambient%lambda > 0) then
                    Cox = 1d0 + 0.23d0 * E2Col(msl)%Instr%ko / E2Col(msl)%Instr%kw &
                        * Ambient%Bowen * Ambient%lambda / Ambient%Ta
                    Stats%Cov(w, msl) = Cox * Stats%Cov(w, msl)
                    Stats%Cov(msl, msl) = Cox**2 * Stats%Cov(msl, msl)
                    !> Alternative formulation by T.W. Horst
                    !> http://www.eol.ucar.edu/instrumentation/sounding&
                    !> &/isfs/isff-support-center/how-tos/&
                    !> $corrections-to-sensible-and-latent-heat-flux-measurements
                    !Stats%Cov(w, msl) = Stats%Cov(w, msl) / (1 - 8d0 * 0.23d0 &
                    !* E2Col(msl)%Instr%ko / E2Col(msl)%Instr%kw * Ambient%bowen)
                end if
        end select
    end do

    !> Sensible heat flux, H in [W m-2]
    Flux1%H = Flux0%H

    !> Internal sensible heat flux, Hint in [W m-2]
    !> Pass-through: the whole gas block, not the four slots those names pick
    !> out. Which four they picked out depended on where water sat, so two
    !> projects differing only in record order disagreed about H_CELL.
    Flux1%Hi_gas(firstGas:lastGas) = Flux0%Hi_gas(firstGas:lastGas)

    !> Level 1 all gases.
    !>
    !> Closed path: Level 1 is Level 0 unchanged. Open path: apply the
    !> bandpass correction factor. One loop over the configured gases,
    !> replacing four near-identical blocks. BPCF%of is indexed by the w_*
    !> covariance labels, which carry the same numbering as the gas slots.
    do msl = firstGas, lastGas
        if (E2Col(msl)%Instr%path_type /= 'closed' .and. BPCF%of(msl) /= error) then
            Flux1%gas(msl) = Flux0%gas(msl) * BPCF%of(msl)
        else
            Flux1%gas(msl) = Flux0%gas(msl)
        end if
        if (Flux0%gas(msl) == error) Flux1%gas(msl) = error
    end do

    !> The water flux carries evapotranspiration and latent heat with it.
    !> Those are scalars - one per project, from the primary H2O slot - so
    !> they are corrected here rather than inside the loop.
    if (wsl >= firstGas .and. BPCF%of(wsx) /= error) then
    if (E2Col(wsl)%Instr%path_type /= 'closed') then
        Flux1%E   = Flux0%E   * BPCF%of(wsl)
        Flux1%ET  = Flux0%ET  * BPCF%of(wsl)
        Flux1%LE  = Flux0%LE  * BPCF%of(wsl)
    else
        Flux1%E   = Flux0%E
        Flux1%ET  = Flux0%ET
        Flux1%LE  = Flux0%LE
    end if
    else
        Flux1%E   = Flux0%E
        Flux1%ET  = Flux0%ET
        Flux1%LE  = Flux0%LE
    end if
    if (wsl >= firstGas) then
    if (Flux0%gas(wsl) == error) then
        Flux1%E   = error
        Flux1%ET  = error
        Flux1%LE  = error
    end if
    end if

    !> Level 1 evapotranspiration fluxes with H2O covariances
    !> at time-lags of other scalars. Do nothing, no spectral correction needed
    Flux1%E_gas(firstGas:lastGas) = Flux0%E_gas(firstGas:lastGas)

    !> Momentum flux [kg m-1 s-2] and friction velocity [m s-1]
    if (BPCF%of(w_u) /= error) then
        Flux1%tau = Flux0%tau * BPCF%of(w_u)
        if (Ambient%us /= error) &
            Ambient%us = Ambient%us * dsqrt(BPCF%of(w_u))
    else
        Flux1%tau = Flux0%tau
    end if
    if (Flux0%tau == error) Flux1%tau = error
    Flux1%ustar = Ambient%us

    call LogSay(' Done.')
end subroutine Fluxes1_rp
