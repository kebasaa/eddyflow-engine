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
    use m_typedef, only: ProcessingRowCount, IsH2OProcessingRow, EnsureExProcessingRows
    implicit none
    !> In/out variables
    type(ExType), intent(inout) :: lEx
    !> local variables
    integer :: i, nrows
    logical :: mirrored_co2, mirrored_h2o, mirrored_ch4, mirrored_other
    real(kind = dbl) :: row_bpcf

    Flux1 = errFlux
    call EnsureExProcessingRows(lEx, EddyFlowProj%processing)
    nrows = ProcessingRowCount(lEx%processing)

    !> Sensible heat flux, H in [W m-2]
    Flux1%H = lEx%Flux0%H

    mirrored_co2 = .false.
    mirrored_h2o = .false.
    mirrored_ch4 = .false.
    mirrored_other = .false.
    do i = 1, nrows
        row_bpcf = lEx%Flux0%gas(i)%bpcf
        if (row_bpcf == error) row_bpcf = 1d0
        Flux1%gas(i) = lEx%Flux0%gas(i)
        if (lEx%Flux0%gas(i)%flux0 /= error) then
            Flux1%gas(i)%flux1 = lEx%Flux0%gas(i)%flux0 * row_bpcf
        else
            Flux1%gas(i)%flux1 = error
        end if
        if (IsH2OProcessingRow(lEx%processing%rows(i))) then
            if (lEx%Flux0%gas(i)%evap0 == error .and. lEx%Flux0%gas(i)%le0 /= error &
                .and. lEx%lambda > 0d0) &
                Flux1%gas(i)%evap1 = lEx%Flux0%gas(i)%le0 / lEx%lambda * row_bpcf
            if (lEx%Flux0%gas(i)%evap0 /= error) Flux1%gas(i)%evap1 = lEx%Flux0%gas(i)%evap0 * row_bpcf
            if (lEx%Flux0%gas(i)%et0 /= error) Flux1%gas(i)%et1 = lEx%Flux0%gas(i)%et0 * row_bpcf
            if (lEx%Flux0%gas(i)%le0 /= error) Flux1%gas(i)%le1 = lEx%Flux0%gas(i)%le0 * row_bpcf
        end if
        select case (trim(lEx%processing%rows(i)%gas_name))
            case ('co2')
                if (.not. mirrored_co2) then
                    Flux1%co2 = Flux1%gas(i)%flux1
                    Flux1%Hi_co2 = Flux1%gas(i)%internal_heat1
                    mirrored_co2 = .true.
                end if
            case ('h2o')
                if (.not. mirrored_h2o) then
                    Flux1%h2o = Flux1%gas(i)%flux1
                    Flux1%E = Flux1%gas(i)%evap1
                    Flux1%ET = Flux1%gas(i)%et1
                    Flux1%LE = Flux1%gas(i)%le1
                    Flux1%Hi_h2o = Flux1%gas(i)%internal_heat1
                    mirrored_h2o = .true.
                end if
            case ('ch4')
                if (.not. mirrored_ch4) then
                    Flux1%ch4 = Flux1%gas(i)%flux1
                    Flux1%Hi_ch4 = Flux1%gas(i)%internal_heat1
                    mirrored_ch4 = .true.
                end if
            case default
                if (.not. mirrored_other) then
                    Flux1%gas4 = Flux1%gas(i)%flux1
                    Flux1%Hi_gas4 = Flux1%gas(i)%internal_heat1
                    mirrored_other = .true.
                end if
        end select
    end do

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
