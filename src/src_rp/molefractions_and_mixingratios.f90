!***************************************************************************
! molefractions_and_mixingratios.f90
! ----------------------------------
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
! \brief       Calculate average mole fractions, mixing ratios and
!              molar density
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo        change name (include densities)
!***************************************************************************
subroutine MoleFractionsAndMixingRatios()
    use m_rp_global_var
    implicit none
    !> local variables
    integer :: gas
    integer :: msl
    real(kind = dbl) :: LocVa(GHGNumVar)
    !> The water this gas is diluted by, as a mixing ratio and a mole
    !> fraction, from whichever source its reference names.
    real(kind = dbl) :: waterR
    real(kind = dbl) :: waterChi
    logical :: haveWater
    logical, external :: GasSlotIsWater

    !> Initialization
    do gas = firstGas, lastGas
        Stats%d(gas) = error
        Stats%chi(gas) = error
        Stats%r(gas) = error
        !> Set correct air molar volume, depending
        !> on instrument path type (closed or open)
        if (E2Col(gas)%present &
            .and. E2Col(gas)%instr%path_type == 'closed') then
            LocVa(gas) = E2Col(gas)%Va
        else
            LocVa(gas) = Ambient%Va
        end if
    end do

    !> Pass 1: every water record. Water converts without a dilution term -
    !> it is the dilutant - so it has to be done before the gases that name
    !> it. This was a single block for the h2o slot, which is record two and
    !> holds water only by convention: a project that declared its water
    !> elsewhere had that block compute a non-water gas's chi/r/d, and then
    !> diluted every real gas by it.
    do gas = firstGas, lastGas
        if (.not. GasSlotIsWater(gas)) cycle
        select case (E2Col(gas)%measure_type)
            case ('mixing_ratio')
                Stats%r(gas)   = Stats%Mean(gas)
                Stats%chi(gas) = Stats%Mean(gas) / (1d0 + Stats%Mean(gas) / 1d3)
                if (LocVa(gas) > 0d0) then
                    Stats%d(gas) = Stats%chi(gas) / LocVa(gas)
                else
                    Stats%d(gas) = error
                end if
            case ('mole_fraction')
                Stats%chi(gas) = Stats%Mean(gas)
                Stats%r(gas)   = Stats%Mean(gas) / (1.d0 - Stats%Mean(gas) / 1d3)
                if (LocVa(gas) > 0d0) then
                    Stats%d(gas) = Stats%chi(gas) / LocVa(gas)
                else
                    Stats%d(gas) = error
                end if
            case ('molar_density')
                Stats%d(gas) = Stats%Mean(gas)
                if (LocVa(gas) > 0d0) then
                    Stats%chi(gas) = Stats%Mean(gas) * LocVa(gas)
                    Stats%r(gas)   = Stats%chi(gas) &
                        / (1.d0 - Stats%chi(gas) / 1d3)
                else
                    Stats%chi(gas) = error
                    Stats%r(gas) = error
                end if
            case default
                Stats%d(gas) = error
                Stats%r(gas) = error
                Stats%chi(gas) = error
        end select
    end do

    !> Pass 2: every other gas, diluted by the water *it* names. A gas whose
    !> moisture reference does not resolve is left undiluted rather than
    !> diluted by an arbitrary slot - the same fallback the single-water code
    !> used when Stats%r(h2o) was in error.
    do gas = firstGas, lastGas
        if (GasSlotIsWater(gas)) cycle
        !> The water this gas names, resolved once: the biomet humidity where
        !> that is what it names, the hygrometer's own statistics otherwise.
        !> `haveWater` false means nothing resolved, and the gas is then left
        !> undiluted rather than diluted by an arbitrary slot - the fallback
        !> the single-water code used when Stats%r(h2o) was in error.
        msl = E2Col(gas)%moist_ref
        haveWater = .false.
        waterR   = error
        waterChi = error
        if (msl == biometMoistRef) then
            haveWater = .true.
            waterR   = Ambient%r_biomet
            waterChi = Ambient%chi_biomet
        elseif (msl >= firstGas .and. msl <= lastGas .and. msl /= gas) then
            haveWater = .true.
            waterR   = Stats%r(msl)
            waterChi = Stats%chi(msl)
        end if

        select case (E2Col(gas)%measure_type)
            case('mixing_ratio')
                Stats%r(gas)   = Stats%Mean(gas)
                if (haveWater .and. waterR /= error) then
                    Stats%chi(gas) = Stats%Mean(gas) &
                        / (1.d0 + waterR * 1d-3)
                else
                    Stats%chi(gas) = Stats%r(gas)
                end if
                if (LocVa(gas) > 0d0) then
                    Stats%d(gas) = Stats%chi(gas) / LocVa(gas) * 1d-3
                else
                    Stats%d(gas) = error
                end if
            case('mole_fraction')
                Stats%chi(gas) = Stats%Mean(gas)
                if (haveWater .and. waterChi /= error) then
                    Stats%r(gas) = Stats%chi(gas) &
                        / (1.d0 - waterChi * 1d-3)
                else
                    Stats%r(gas) = Stats%chi(gas)
                end if
                if (LocVa(gas) > 0d0) then
                    Stats%d(gas) = Stats%chi(gas) / LocVa(gas) * 1d-3
                else
                    Stats%d(gas) = error
                end if
            case('molar_density')
                Stats%d(gas) = Stats%Mean(gas)
                if (LocVa(gas) > 0d0) then
                    Stats%chi(gas) = Stats%Mean(gas) * LocVa(gas) * 1d3
                    if (haveWater .and. waterChi /= error) then
                        Stats%r(gas) = Stats%chi(gas) &
                            / (1.d0 - waterChi * 1d-3)
                    else
                        Stats%r(gas) = Stats%chi(gas)
                    end if
                else
                    Stats%chi(gas) = error
                    Stats%r(gas) = error
                end if
        end select
    end do
end subroutine MoleFractionsAndMixingRatios
