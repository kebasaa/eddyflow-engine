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
    integer :: nslot
    integer :: nreal
    !> The widest stride worth following - a rate ratio past this is not a
    !> sub-sampled instrument, it is a misdeclared one.
    integer, parameter :: MaxSampleStride = 1000
    integer, allocatable :: rows(:)
    logical, allocatable :: isreal(:)
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

    !> Where a slower column's real samples are, recorded HERE because the
    !> interpolation below is about to remove the only evidence of them.
    !>
    !> The rows themselves, not one phase per column. A phase was the
    !> commonest offset in the first twenty intervals, which describes a
    !> column on a fixed grid exactly - and on such a column the rows are
    !> exactly the phase rows, so nothing it computed moves. But the Yatir
    !> laser drifts through every row position in a half-hour, and there no
    !> phase is right for more than a stretch: most rebuilt samples read an
    !> interpolated blend of two real ones. A missed sample gets a slot at its
    !> place in the sequence; its value is the interpolated one, as before.
    allocate(rows(nrow), isreal(nrow))
    do j = 1, ncol
        SpecSamples(j)%n = 0
        if (allocated(SpecSamples(j)%r)) deallocate(SpecSamples(j)%r)
        stride = nint(Metadata%ac_freq / ColumnAcFreq(j))
        if (stride <= 1 .or. stride > MaxSampleStride) cycle
        call SlowColumnSampleRows(Set(1:nrow, j), nrow, stride, error, &
            rows, isreal, nslot, nreal)
        if (nslot < 2) cycle
        allocate(SpecSamples(j)%r(nslot))
        SpecSamples(j)%r = rows(1:nslot)
        SpecSamples(j)%n = nslot
    end do
    deallocate(rows, isreal)

    !> The shift ReplaceGapWithLinearInterpolation is about to apply to each
    !> column: it removes leading error rows by moving the column up.
    SpecLead = 0
    do j = 1, GHGNumVar
        if (.not. SpecCol(j)%present) cycle
        do i = 1, nrow
            if (Set(i, j) /= error) exit
        end do
        if (i <= nrow) SpecLead(j) = i - 1
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
