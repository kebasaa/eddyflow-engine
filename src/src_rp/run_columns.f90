!***************************************************************************
! run_columns.f90
! ---------------
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
! \brief       Which gases a run had, over the whole run.
!
!              What the end of a pass writes - the time-lag optimisation and
!              PWB summaries, the windows the optimiser fits, the donors, the
!              storage cleaning - goes by E2Col, which says what the LAST period
!              had: a gas that ran out of data in the final half-hour was left
!              out of the summary as if it had never been measured, although
!              every other period's lags for it were in the dataset.
!
!              So each processed period notes the gases it had, keeping each
!              gas's column description from the last period that had it, and
!              before the end of the pass, a gas the last period lacks is put
!              back as that period left it. A gas the last period has is not
!              touched: its description is the last period's, as before.
!
!              A pass split across workers merges the notes piece by piece in
!              time order, a later piece's overriding an earlier's, so the
!              result is a single pass's. A worker notes only its own piece -
!              not the periods it repeats before it.
!
! \author      Jonathan Muller
! \sa          eddyflow-rp_main.f90, prepass_parallel.f90,
!              production_parallel.f90
!***************************************************************************
module m_run_columns
    use m_rp_global_var
    implicit none
    private

    public :: ResetRunColumns, NoteRunColumns, UseRunColumns
    public :: WriteRunColumns, MergeRunColumns

    !> Each gas slot's column description from the last period that had it.
    type(ColType) :: RunCol(E2NumVar)
    logical :: RunColSeen(E2NumVar) = .false.

contains

    subroutine ResetRunColumns()
        RunColSeen = .false.
    end subroutine ResetRunColumns

    !> A period has been processed: note the gases it had.
    subroutine NoteRunColumns()
        integer :: j

        do j = firstGas, lastGas
            if (.not. E2Col(j)%present) cycle
            RunCol(j) = E2Col(j)
            RunColSeen(j) = .true.
        end do
    end subroutine NoteRunColumns

    !> Before the end of a pass: the gases the last period lacked, as the
    !> last period that had them left them.
    subroutine UseRunColumns()
        integer :: j

        do j = firstGas, lastGas
            if (E2Col(j)%present .or. .not. RunColSeen(j)) cycle
            E2Col(j) = RunCol(j)
        end do
    end subroutine UseRunColumns

    !> A worker's notes, onto its dump.
    subroutine WriteRunColumns(u)
        integer, intent(in) :: u

        write(u) RunColSeen
        write(u) RunCol
    end subroutine WriteRunColumns

    !> The parent: a later piece's notes, read from its dump, over its own.
    subroutine MergeRunColumns(u)
        integer, intent(in) :: u
        integer :: j
        logical :: seen(E2NumVar)
        type(ColType) :: cols(E2NumVar)

        read(u) seen
        read(u) cols
        do j = 1, E2NumVar
            if (.not. seen(j)) cycle
            RunCol(j) = cols(j)
            RunColSeen(j) = .true.
        end do
    end subroutine MergeRunColumns

end module m_run_columns
