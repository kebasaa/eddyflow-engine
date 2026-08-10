!***************************************************************************
! air_and_cell_parameters.f90
! ---------------------------
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
! \brief       Calculate ambient and cell average T, P and molar volume
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine AirAndCellParameters()
    use m_rp_global_var
    implicit none
    !> local variables
    integer :: gas
    integer :: cellBase

    !> Air temperature/pressure estimates
    !> Last true condition determines which temperature is used
    Stats%T = Stats%Mean(ts)
    if(Stats%Mean(te)  > 220d0 .and. Stats%Mean(te) < 340d0) Stats%T = Stats%Mean(te)
    if(biomet%val(bTa) > 220d0 .and. biomet%val(bTa) < 340d0) Stats%T = biomet%val(bTa)

    !> Last true condition determines which pressure is used
    Stats%Pr = Metadata%bar_press
    if(Stats%Mean(pe)  > 40000 .and. Stats%Mean(pe)  < 110000) Stats%Pr = Stats%Mean(pe)
    if(biomet%val(bPa) > 40000 .and. biomet%val(bPa) < 110000) Stats%Pr = biomet%val(bPa)

    !> Ambient air molar volume [m+3 mol-1] and air mass density [kg m-3]
    if (Stats%Pr > 0d0 .and. Stats%T /= error) then
        Ambient%Va = Ru * Stats%T / Stats%Pr
    else
        Ambient%Va = error
    end if

    !> Cell Temperature, if applicable \n
    if (Stats%Mean(tc) /= error) then
        Ambient%Tcell = Stats%Mean(tc)
    else
        Ambient%Tcell = Stats%T
    end if

    !> Cell pressure, if applicable \n
    if(Stats%Mean(pi) /= error) then
        Ambient%Pcell = Stats%Mean(pi)
    elseif (Stats%Pr /= error) then
        Ambient%Pcell = Stats%Pr
    else
        Ambient%Pcell = Metadata%bar_press
    end if

    !> Using cell temperature, for each gas column related to a closed-path analyser,
    !> determine cell air molar volume. If Tcell == Stats%T (that is, if no internal temperature
    !> is provided) and Pcell == Stats%Pr, then cell air molar volume is equal to ambient air molar volume
    !> Cell conditions of the analyser that measured each gas.
    !>
    !> cell_ref is the base slot of that instrument's cell block, resolved in
    !> DefineE2Set; a gas whose analyser has no cell record points at the first
    !> block, which is where a single-analyser project's data has always been.
    !> Falls back to the scalars above whenever a block carries no reading, so
    !> a project with one cell record behaves exactly as it did.
    E2Col%Va = error
    Ambient%Tcell_at = error
    Ambient%Pcell_at = error
    do gas = firstGas, lastGas
        if (.not. E2Col(gas)%present) cycle

        !> Zero when no cell record belongs to this gas's analyser. It used to
        !> be coerced to firstCell here, which handed the first analyser's cell
        !> to a gas measured in a different one; DefineE2Set now resolves the
        !> field to zero in that case and every reader declines together.
        !> The site scalars stand in, as they do for a gas whose block carries
        !> no reading.
        cellBase = E2Col(gas)%cell_ref
        if (cellBase < firstCell .or. cellBase > lastCell) then
            Ambient%Tcell_at(gas) = Ambient%Tcell
            Ambient%Pcell_at(gas) = Ambient%Pcell
            if (E2Col(gas)%instr%path_type == 'closed') then
                if (Ambient%Pcell_at(gas) > 0d0 &
                    .and. Ambient%Tcell_at(gas) /= error) then
                    E2Col(gas)%Va = Ru * Ambient%Tcell_at(gas) / Ambient%Pcell_at(gas)
                else
                    E2Col(gas)%Va = error
                end if
            end if
            cycle
        end if

        if (Stats%Mean(cellBase) /= error) then
            Ambient%Tcell_at(gas) = Stats%Mean(cellBase)
        else
            Ambient%Tcell_at(gas) = Ambient%Tcell
        end if

        if (Stats%Mean(cellBase + 3) /= error) then
            Ambient%Pcell_at(gas) = Stats%Mean(cellBase + 3)
        else
            Ambient%Pcell_at(gas) = Ambient%Pcell
        end if

        if (E2Col(gas)%instr%path_type == 'closed') then
            if (Ambient%Pcell_at(gas) > 0d0 &
                .and. Ambient%Tcell_at(gas) /= error) then
                E2Col(gas)%Va = Ru * Ambient%Tcell_at(gas) / Ambient%Pcell_at(gas)
            else
                E2Col(gas)%Va = error
            end if
        end if
    end do
end subroutine AirAndCellParameters
