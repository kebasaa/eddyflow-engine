!***************************************************************************
! point_by_point_to_mixing_ratio.f90
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
! \brief       Convert mole fractions (moles per mole wet air) or \n
!              molar density to mixing ratios (moles per mole dry air)
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine PointByPointToMixingRatio(Set, nrow, ncol, printout)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    logical, intent(in) :: printout
    real(kind = dbl), intent(inout) :: Set(nrow, ncol)
    !> Local variables
    real(kind = dbl), allocatable :: H2Omf(:, :)
    real(kind = dbl) :: Va(nrow)
    integer :: watSlot(MaxNumGases)
    integer :: nwat
    integer :: gas, msl, k, iw
    !> Base slot of the cell block belonging to the analyser that measured the
    !> gas being converted: cellBase + 0..3 is cell_t, int_t_1, int_t_2, int_p.
    integer :: cellBase
    logical :: cellVaAvailable
    logical, external :: GasSlotIsWater

    !> Every water measurement the project describes, in slot order. This
    !> used to be the single h2o slot, which is record two and holds water
    !> only by convention; with one hygrometer per analyser there is more
    !> than one, and a gas has to be diluted by the water its *own* analyser
    !> measured.
    nwat = 0
    do gas = firstGas, lastGas
        if (gas - firstGas + 1 > min(EddyFlowProj%gas_num, MaxNumGases)) exit
        if (.not. GasSlotIsWater(gas)) cycle
        if (.not. E2Col(gas)%present) cycle
        if (E2Col(gas)%instr%path_type /= 'closed') cycle
        nwat = nwat + 1
        watSlot(nwat) = gas
    end do
    if (nwat == 0) return

    !> Indexed by water, not by slot: at 20 Hz over half an hour a
    !> (nrow, GHGNumVar) array would be some 157 MB.
    allocate(H2Omf(nrow, nwat))

    if (printout) write(*, '(a)', advance = 'no') &
        '  WPL step: converting into mixing ratios wherever possible..'

    !> Pass 1: each water slot's mole fraction, taken from the *raw* series,
    !> and then that slot's own conversion. Water must be converted before
    !> any gas that references it, and H2Omf must be taken before the
    !> conversion - which is why this is a separate pass rather than one
    !> loop with a water special case.
    do iw = 1, nwat
        msl = watSlot(iw)

        !> Air molar volume in the cell of the analyser that measured *this*
        !> water, from that analyser's own cell block.
        !>
        !> This read tc and pi - the first block - and gated them on whether
        !> that block's instrument matched this hygrometer's, which was a way
        !> of saying "only if these cell columns belong to this analyser". It
        !> gave the right answer for instrument one and no answer for any
        !> other, so a second hygrometer reporting molar density was left
        !> unconvertible. cell_ref states the block directly, and the gas loop
        !> below now reads it the same way.
        !>
        !> Zero when no cell record belongs to this hygrometer's analyser, in
        !> which case there is no cell molar volume to be had - the first block
        !> is a different instrument's cell, not this one's.
        !> Nested, not one .and. chain: Fortran does not promise to stop
        !> evaluating a compound condition once it is decided, so an
        !> unresolved cell_ref of zero reaches E2Col(0) and the run dies in
        !> the WPL step. The same note stands over the diagnostics loop in
        !> DefineE2Set, for the same reason.
        cellBase = E2Col(msl)%cell_ref
        Va(:) = error
        if (cellBase >= firstCell .and. cellBase <= lastCell) then
            if (E2Col(cellBase)%present .and. E2Col(cellBase + 3)%present) then
                where (Set(:, cellBase + 3) > 0d0 .and. Set(:, cellBase) > 0d0)
                    Va(:) = Ru * Set(:, cellBase) / Set(:, cellBase + 3)
                elsewhere
                    Va(:) = error
                end where
            end if
        end if
        cellVaAvailable = .not. all(Va(:) == error)

        select case (E2Col(msl)%measure_type)
            case ('mixing_ratio')
                where(Set(:, msl) /=  error)
                    H2Omf(:, iw) = Set(:, msl) / (1d0 + Set(:, msl) / 1d3)
                elsewhere
                    H2Omf(:, iw) = error
                end where
            case ('mole_fraction')
                H2Omf(:, iw) = Set(:, msl)
            case ('molar_density')
                where (Va(:) /= error .and. Set(:, msl) /= error)
                    H2Omf(:, iw) = Set(:, msl) * Va(:)
                elsewhere
                    H2Omf(:, iw) = error
                end where
            case default
                H2Omf(:, iw) = error
        end select

        !> Water's own conversion
        if (E2Col(msl)%measure_type == 'mole_fraction') then
            E2Col(msl)%measure_type = 'mixing_ratio'
            where(H2Omf(:, iw) /= error .and. Set(:, msl) /= error)
                Set(:, msl) = Set(:, msl) / (1.d0 - H2Omf(:, iw) * 1d-3)
            elsewhere
                Set(:, msl) = error
            endwhere
        elseif (E2Col(msl)%measure_type == 'molar_density' &
            .and. cellVaAvailable) then
            E2Col(msl)%measure_type = 'mixing_ratio'
            where(Va(:) /= error .and. H2Omf(:, iw) /= error)
                Set(:, msl) = Set(:, msl) * Va(:) &
                    / (1.d0 - H2Omf(:, iw) * 1d-3)
            elsewhere
                Set(:, msl) = error
            endwhere
        end if
    end do

    !> Pass 2: every other gas, diluted by the water it names. Water slots
    !> are skipped - pass 1 already converted them, and converting twice
    !> would apply the correction to an already-corrected series.
    do gas = firstGas, lastGas
        if (gas - firstGas + 1 > min(EddyFlowProj%gas_num, MaxNumGases)) exit
        if (.not. E2Col(gas)%present) cycle
        if (GasSlotIsWater(gas)) cycle

        !> The gas's own moisture reference, resolved in DefineE2Set.
        msl = E2Col(gas)%moist_ref
        if (msl < firstGas .or. msl > lastGas) cycle
        iw = 0
        do k = 1, nwat
            if (watSlot(k) == msl) iw = k
        end do
        if (iw == 0) cycle

        !> Whatever water the record names, including one on another analyser.
        !>
        !> A same-analyser test stood here and declined the pairing outright,
        !> so a gas whose own instrument carries no hygrometer was not diluted
        !> at all - while MoistTerms went on taking its sigma and rho_w from
        !> the very water this refused. Declining is not the neutral choice it
        !> looks like: it leaves the gas a wet mole fraction corrected with
        !> terms derived from a water it is held not to share.
        !>
        !> Borrowing across analysers is a compromise, and the instantaneous
        !> humidity in another cell is the weaker half of it.
        !> ExceptionHandler(106) says so once per gas, and declaring an H2O on
        !> the gas's own analyser makes ResolveGasRef prefer it.

        !> Cell conditions of the analyser that measured THIS gas.
        !>
        !> Asked of tc and pi - the first cell block - and gated on the
        !> *water's* model. Both were the same thing while one global cell
        !> served every analyser; with per-instrument blocks the first is
        !> whichever analyser happens to hold cell record one, and the second
        !> asks the wrong instrument entirely. A gas was converted with
        !> another analyser's cell temperature and pressure, which is exactly
        !> the mixing this routine exists to avoid. cell_ref is resolved in
        !> DefineE2Set and is what AirAndCellParameters reads.
        !> Zero when no cell record belongs to this gas's analyser.
        !> Nested, not one .and. chain: Fortran does not promise to stop
        !> evaluating a compound condition once it is decided, so an
        !> unresolved cell_ref of zero reaches E2Col(0) and the run dies in
        !> the WPL step. The same note stands over the diagnostics loop in
        !> DefineE2Set, for the same reason.
        cellBase = E2Col(gas)%cell_ref
        Va(:) = error
        if (cellBase >= firstCell .and. cellBase <= lastCell) then
            if (E2Col(cellBase)%present .and. E2Col(cellBase + 3)%present) then
                where (Set(:, cellBase + 3) > 0d0 .and. Set(:, cellBase) > 0d0)
                    Va(:) = Ru * Set(:, cellBase) / Set(:, cellBase + 3)
                elsewhere
                    Va(:) = error
                end where
            end if
        end if
        cellVaAvailable = .not. all(Va(:) == error)

        if (E2Col(gas)%measure_type == 'mole_fraction') then
            E2Col(gas)%measure_type = 'mixing_ratio'
            where(H2Omf(:, iw) /= error .and. Set(:, gas) /= error)
                Set(:, gas) = Set(:, gas) / (1.d0 - H2Omf(:, iw) * 1d-3)
            elsewhere
                Set(:, gas) = error
            endwhere
        elseif (E2Col(gas)%measure_type == 'molar_density' &
            .and. cellVaAvailable) then
            E2Col(gas)%measure_type = 'mixing_ratio'
            where(Va(:) /= error .and. H2Omf(:, iw) /= error &
                .and. Set(:, gas) /= error)
                Set(:, gas) = Set(:, gas) * Va(:) * 1d3 &
                    / (1.d0 - H2Omf(:, iw) * 1d-3)
            elsewhere
                Set(:, gas) = error
            endwhere
        end if
    end do

    deallocate(H2Omf)
    if (printout) write(*,'(a)') ' Done.'
end subroutine PointByPointToMixingRatio
