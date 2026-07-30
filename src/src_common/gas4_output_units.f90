!***************************************************************************
! gas4_output_units.f90
! ---------------------
! Copyright © 2026, ETH Zurich, Jonathan Muller
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
! \brief       Defines gas4 full-output scales and labels from metadata units.
! \author      Jonathan Muller
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
!***************************************************************************
!
! \brief       Mole-basis scale of a gas's FLUXNET columns, by species.
! \author      Jonathan Muller
! \note        The FLUXNET row carries no units line, and the FP-In format
!              fixes the basis per species: CO2 in umol mol-1, H2O in
!              mmol mol-1, every other species in nmol mol-1. Internally the
!              engine holds every trace gas in umol mol-1 and water in
!              mmol mol-1 (ConvertTraceGasUnits), so these factors take the
!              internal value to the column's basis and nothing else.
!
!              This replaces a gain hard-coded for whatever gases occupied the
!              ch4 and gas4 slots. That named a position rather than a
!              species, so the same gas was reported in nmol mol-1 from slot 7
!              and in umol mol-1 from slot 9. The species comes from the gas
!              record, which is the authority everywhere else in this effort.
!
!              Deliberately NOT derived from the input unit: the raw signal was
!              already normalised on input, so scaling by the input factor
!              again would count it twice. The full output is the place where
!              the project's own unit is honoured, because it carries a units
!              row to declare it.
!***************************************************************************
real(kind = dbl) function FluxnetGasScale(gas_slot)
    use m_common_global_var
    implicit none
    integer, intent(in) :: gas_slot
    character(32) :: species
    integer :: rec

    !> Anything not named by a record is a trace gas by default.
    FluxnetGasScale = 1d3
    rec = gas_slot - firstGas + 1
    if (rec < 1 .or. rec > min(EddyFlowProj%gas_num, MaxNumGases)) return

    species = EddyFlowProj%gas(rec)%var
    call uppercase(species)
    select case (trim(adjustl(species)))
        case ('CO2', 'H2O'); FluxnetGasScale = 1d0
        case default;        FluxnetGasScale = 1d3
    end select
end function FluxnetGasScale

!***************************************************************************
!
! \brief       Scale of a gas's FLUXNET vertical-advection column, by species.
! \author      Jonathan Muller
! \note        Advection is w times the molar density, so it starts from the
!              mmol m-3 basis rather than the mole basis: one further factor of
!              1d3 for every species except water, whose target already is the
!              mmol basis. Reproduces the 1d3 / none / 1d6 three-way switch
!              this replaces, without naming a slot.
!***************************************************************************
real(kind = dbl) function FluxnetGasAdvScale(gas_slot)
    use m_common_global_var
    implicit none
    integer, intent(in) :: gas_slot
    character(32) :: species
    integer :: rec
    real(kind = dbl), external :: FluxnetGasScale

    rec = gas_slot - firstGas + 1
    if (rec >= 1 .and. rec <= min(EddyFlowProj%gas_num, MaxNumGases)) then
        species = EddyFlowProj%gas(rec)%var
        call uppercase(species)
        if (trim(adjustl(species)) == 'H2O') then
            FluxnetGasAdvScale = 1d0
            return
        end if
    end if
    FluxnetGasAdvScale = FluxnetGasScale(gas_slot) * 1d3
end function FluxnetGasAdvScale

subroutine Gas4FullOutputUnits(unit_in, flux_scale, dens_scale, &
    flux_label, conc_label, mixr_label, dens_label)
    use m_common_global_var
    implicit none
    !> in/out variables
    character(*), intent(in) :: unit_in
    real(kind = dbl), intent(out) :: flux_scale
    real(kind = dbl), intent(out) :: dens_scale
    character(*), intent(out) :: flux_label
    character(*), intent(out) :: conc_label
    character(*), intent(out) :: mixr_label
    character(*), intent(out) :: dens_label

    select case (trim(adjustl(unit_in)))
        case ('ppb', 'nmol_mol', 'nmol/mol')
            flux_scale = 1d3
            dens_scale = 1d6
            flux_label = '[nmol+1s-1m-2]'
            conc_label = '[nmol+1mol_a-1]'
            mixr_label = '[nmol+1mol_d-1]'
            dens_label = '[nmol+1m-3]'
        case ('pmol_mol', 'pmol/mol')
            flux_scale = 1d6
            dens_scale = 1d9
            flux_label = '[pmol+1s-1m-2]'
            conc_label = '[pmol+1mol_a-1]'
            mixr_label = '[pmol+1mol_d-1]'
            dens_label = '[pmol+1m-3]'
        case default
            flux_scale = 1d0
            dens_scale = 1d0
            flux_label = '[' // char(194) // char(181) // 'mol+1s-1m-2]'
            conc_label = '[' // char(194) // char(181) // 'mol+1mol_a-1]'
            mixr_label = '[' // char(194) // char(181) // 'mol+1mol_d-1]'
            dens_label = '[mmol+1m-3]'
    end select
end subroutine Gas4FullOutputUnits
