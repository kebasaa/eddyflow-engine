!***************************************************************************
! fluxes1.f90
! -----------
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
! \todo
!***************************************************************************
subroutine Fluxes1(lEx)
    use m_fx_global_var
    implicit none
    integer :: msl
    !> In/out variables
    type(ExType), intent(inout) :: lEx
    !> local variables
    real(kind = dbl)  :: Cox
    integer :: gas
    integer :: wsl
    include '../src_common/interfaces_1.inc'

    Flux1 = errFlux

    !> First, apply oxygen correction to Krypton and Lyman-alpha hygrometers,
    !> according to van Dijk et al. (2003, JAOT, eq. 13b)
    !>
    !> Per hygrometer, not per slot. This asked lEx%instr(ih2o) - the water
    !> role of a five-wide instrument numbering that ran alongside the gas
    !> slots - and corrected lEx%cov_w(h2o). A site with two hygrometers had
    !> only one of them corrected, and a site whose water is not record two
    !> had the correction applied to whatever gas held slot six.
    do gas = firstGas, lastGas
        if (.not. GasSlotIsWater(gas)) cycle
        select case (lEx%gas_instr(gas)%model(1:max(1, &
            len_trim(lEx%gas_instr(gas)%model) - 2)))
            case('open_path_krypton','closed_path_krypton', &
                    'open_path_lyman','closed_path_lyman')
                if (lEx%gas_instr(gas)%ko /= error &
                    .and. lEx%gas_instr(gas)%kw /= 0d0 &
                    .and. lEx%Ta > 0d0 .and. lEx%Bowen /= error &
                    .and. lEx%lambda > 0d0) then
                    Cox = 1d0 + 0.23d0 * lEx%gas_instr(gas)%ko &
                        / lEx%gas_instr(gas)%kw &
                        * lEx%Bowen * lEx%lambda / lEx%Ta
                    lEx%cov_w(gas) = Cox * lEx%cov_w(gas)
                    lEx%var(gas) = Cox**2 * lEx%var(gas)
                    !> Alternative formulation by T.W. Horst
                    !> http://www.eol.ucar.edu/instrumentation/&
                    !> &sounding/isfs/isff-support-center/how-tos/&
                    !> &corrections-to-sensible-and-latent-heat-flux-measurements
                    !lEx%cov_w(gas) = lEx%cov_w(gas) / (1 - 8d0 * 0.23d0 &
                    !* lEx%gas_instr(gas)%ko / lEx%gas_instr(gas)%kw * lEx%Bowen)
                endif
        end select
    end do

    !> Sensible heat flux, H in [W m-2]
    Flux1%H = lEx%Flux0%H

    !> Internal sensible heat flux, Hint in [W m-2]
    !> Pass-through: the whole gas block. This one still named the literal
    !> h2o slot, so on a project whose water is not record two it copied
    !> whatever gas sat in slot six and left the real water behind.
    Flux1%Hi_gas(firstGas:lastGas) = lEx%Flux0%Hi_gas(firstGas:lastGas)

    !> Level 1 all gases.
    !>
    !> Closed path: Level 1 is Level 0 unchanged. Open path: apply the
    !> bandpass correction. One loop over the configured gases, replacing four
    !> near-identical blocks; gas_instr is indexed by gas slot, so this reaches
    !> past the four the instrument-role index could address.
    do msl = firstGas, lastGas
        if (lEx%gas_instr(msl)%path_type /= 'closed' .and. BPCF%of(msl) /= error) then
            Flux1%gas(msl) = lEx%Flux0%gas(msl) * BPCF%of(msl)
        else
            Flux1%gas(msl) = lEx%Flux0%gas(msl)
        end if
        if (lEx%Flux0%gas(msl) == error) Flux1%gas(msl) = error
    end do

    !> Evapotranspiration and latent heat travel with the water flux and are
    !> scalars, taken from the primary H2O slot. The ex file carries LE, not
    !> E, so E is derived here before the correction is applied.
    !>
    !> The slot is resolved, not assumed. This said gas_instr(h2o), BPCF%of(w_h2o)
    !> and Flux0%gas(h2o) - slot six - so a project whose water is not record two
    !> corrected the latent heat with another species' transfer function, and
    !> invalidated it on that species' flux rather than on water's. The comment
    !> above already said "the primary H2O slot"; the code did not.
    !>
    !> With no water configured PrimaryWaterSlot returns 0, and the fluxes pass
    !> through uncorrected - which is what the path_type test evaluated to when
    !> the slot held nothing.
    wsl = PrimaryWaterSlot()
    !> Guarded: with no hygrometer LE is `error`, and dividing it by a lambda
    !> that is itself `error` produced a finite nonsense number that then flowed
    !> into Flux1%E.
    if (lEx%Flux0%LE /= error .and. lEx%lambda /= error .and. lEx%lambda /= 0d0) then
        lEx%Flux0%E = lEx%Flux0%LE / lEx%lambda
    else
        lEx%Flux0%E = error
    end if
    if (wsl >= firstGas) then
        if (lEx%gas_instr(wsl)%path_type /= 'closed' .and. BPCF%of(wsl) /= error) then
            Flux1%E   = lEx%Flux0%E   * BPCF%of(wsl)
            Flux1%ET  = lEx%Flux0%ET  * BPCF%of(wsl)
            Flux1%LE  = lEx%Flux0%LE  * BPCF%of(wsl)
        else
            Flux1%E   = lEx%Flux0%E
            Flux1%ET  = lEx%Flux0%ET
            Flux1%LE  = lEx%Flux0%LE
        end if
        if (lEx%Flux0%gas(wsl) == error) then
            lEx%Flux0%E = error
            Flux1%E     = error
            Flux1%ET    = error
            Flux1%LE    = error
        end if
    else
        Flux1%E   = lEx%Flux0%E
        Flux1%ET  = lEx%Flux0%ET
        Flux1%LE  = lEx%Flux0%LE
    end if


    !> Level 1 evapotranspiration fluxes with H2O covariances at time-lags
    !> of other scalars. Do nothing, no spectral correction needed
    Flux1%E_gas(firstGas:lastGas) = lEx%Flux0%E_gas(firstGas:lastGas)

    !> Momentum flux [kg m-1 s-2] and friction velocity [m s-1]
    if (BPCF%of(w_u) /= error) then
        Flux1%tau = lEx%Flux0%tau * BPCF%of(w_u)
        Flux1%ustar = lEx%Flux0%ustar * dsqrt(BPCF%of(w_u))
    else
        Flux1%tau = lEx%Flux0%tau
        Flux1%ustar = lEx%Flux0%ustar
    end if
    if (lEx%Flux0%tau == error) Flux1%tau = error
    if (lEx%Flux0%ustar == error) Flux1%ustar = error
end subroutine Fluxes1
