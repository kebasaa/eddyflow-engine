!***************************************************************************
! fix_dataset_for_spectra.f90
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
! \brief       replace gaps with linear interpolation of neighboring data
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine FixDatasetForSpectra(Set, nrow, ncol, nrow2)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    integer, intent(out) :: nrow2
    real(kind = dbl) :: Set(nrow, ncol)
    !> local variables
    integer :: i
    integer :: j
    integer :: tnrow
    integer :: expected
    integer :: stride
    integer :: offset
    !> How many of the column's own intervals to look at before deciding where
    !> its samples sit, and the widest stride worth tallying - a rate ratio
    !> past this is not a sub-sampled instrument, it is a misdeclared one.
    integer, parameter :: PhaseIntervals = 20
    integer, parameter :: MaxPhaseStride = 1000
    integer :: tally(0:MaxPhaseStride - 1)
    real(kind = dbl), external :: ColumnAcFreq


    !> If more than 30% of the data is missing, don't compute spectra
    !> because linear interpolation probably too severly affect spectral shape
    !> This filter is totally arbitrary, only based on anecdotal evidence
    !>
    !> A third of what the COLUMN should have produced, not a third of the
    !> rows. An instrument slower than the row rate cannot fill them - a 1 Hz
    !> column in a 10 Hz file is nine tenths error rows - so measured against
    !> the rows it was always over the threshold and its spectra were never
    !> computed, whatever the data. At the file's own rate expected is nrow and
    !> this is the test it replaces, exactly.
    do j = 1, ncol
        expected = nint(dble(nrow) * ColumnAcFreq(j) / Metadata%ac_freq)
        if (expected - count(Set(1:nrow, j) /= error) > expected / 3) &
            SpecCol(j)%present = .false.
    end do

    !> Where a slower column's real samples sit, recorded HERE because the
    !> interpolation below is about to remove the only evidence of it. A
    !> column at the file's own rate has every row and so has phase zero.
    !>
    !> The commonest offset over many intervals, not the first one found: a
    !> column whose very first sample happens to be missing - an ordinary gap
    !> at the start of a period - would otherwise report phase zero and have
    !> every one of its rebuilt samples read an interpolated blend instead of a
    !> measurement. A real sampling pattern gives one offset an overwhelming
    !> majority; a column with no clear winner is not on a regular grid, and
    !> zero is as good an answer as any for it.
    SpecPhase(1:ncol) = 0
    do j = 1, ncol
        stride = nint(Metadata%ac_freq / ColumnAcFreq(j))
        if (stride <= 1 .or. stride > MaxPhaseStride) cycle
        tally(0:stride - 1) = 0
        do i = 1, min(nrow, PhaseIntervals * stride)
            if (Set(i, j) /= error) then
                offset = mod(i - 1, stride)
                tally(offset) = tally(offset) + 1
            end if
        end do
        if (any(tally(0:stride - 1) > 0)) &
            SpecPhase(j) = maxloc(tally(0:stride - 1), dim = 1) - 1
    end do

    !> nrow2 is the smallest nrow of all columns
    nrow2 = nrow
    do j = 1, GHGNumVar
        if (SpecCol(j)%present) then
            call ReplaceGapWithLinearInterpolation(Set(1:nrow, j), size(Set, 1), tnrow, error)
            if (tnrow < nrow2) nrow2 = tnrow
        end if
    end do
end subroutine FixDatasetForSpectra
