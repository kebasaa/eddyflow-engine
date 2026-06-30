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
            flux_label = '[' // char(181) // 'mol+1s-1m-2]'
            conc_label = '[' // char(181) // 'mol+1mol_a-1]'
            mixr_label = '[' // char(181) // 'mol+1mol_d-1]'
            dens_label = '[mmol+1m-3]'
    end select
end subroutine Gas4FullOutputUnits
