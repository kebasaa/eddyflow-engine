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
subroutine FixDatasetForSpectra(Set, nrow, ncol, lo, nrow2)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    !> The window every (co)spectrum is computed on: rows lo .. lo + nrow2 - 1
    !> of Set, the same rows for every column.
    integer, intent(out) :: lo
    integer, intent(out) :: nrow2
    real(kind = dbl) :: Set(nrow, ncol)
    !> local variables
    integer :: j
    integer :: hi
    integer :: first
    integer :: last
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

    !> One window for every column: from the first row on which all of them
    !> have a value to the last. Each column used to be interpolated over all
    !> its rows, and ReplaceGapWithLinearInterpolation drops a column's
    !> LEADING error rows by moving the rest of it up - so any column that
    !> opened with a gap was paired with w that many rows out of step in
    !> every cospectrum, and in the covariance it was normalised by. A
    !> negative lag opens a column with exactly such a gap, so for it the
    !> move undid the lag compensation; data missing or filtered at the start
    !> of a period did the same by however many rows it was. On CH-LAE the
    !> MIRO's normalising covariance matched its flux in a fifth to a half of
    !> the half-hours, the sonic's and LI analyser's in three quarters or more.
    !>
    !> Each column is interpolated over its own span, from its first value to
    !> its last, where it opens and closes on a value - so nothing moves, only
    !> interior gaps are filled. Its own span and not the common window: a
    !> slower column has a value only every few rows, and the window can open
    !> between two of them, which would hand the interpolation a leading gap
    !> to shift away again. The lag-shift error rows at either end lie outside
    !> the window, which is what the trim by RowLags in the caller used to
    !> approximate.
    !> j <= lastGas: past the project's last gas the slots are cells now.
    lo = 1
    hi = nrow
    do j = 1, lastGas
        if (.not. SpecCol(j)%present) cycle
        do first = 1, nrow
            if (Set(first, j) /= error) exit
        end do
        if (first > nrow) then
            SpecCol(j)%present = .false.
            cycle
        end if
        do last = nrow, first, -1
            if (Set(last, j) /= error) exit
        end do
        lo = max(lo, first)
        hi = min(hi, last)
        call ReplaceGapWithLinearInterpolation(Set(first:last, j), &
            last - first + 1, tnrow, error)
    end do
    !> Nothing moves now, so a slower column's samples are read where
    !> SpecSamples recorded them; SlowColumnSpectra still subtracts this.
    SpecLead = 0
    nrow2 = max(0, hi - lo + 1)
end subroutine FixDatasetForSpectra
