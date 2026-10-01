!***************************************************************************
! column_sampling.f90
! -------------------
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
! \brief       How fast a column is sampled, and how much of its own data may
!              be missing.
! \author      Jonathan Muller, ETH Zurich
! \note        The raw file has one row rate, but the instruments writing into
!              those rows need not share it: a 1 Hz analyser in a 20 Hz file
!              writes one row in twenty and leaves the other nineteen at the
!              error code. Measured against the row grid that column is 95 %
!              missing, and every completeness test dropped it - the data was
!              complete for its own rate and was discarded anyway.
!
!              These two answer for a column rather than for the file, so a
!              test can ask what THIS column should have produced. Both fall
!              back to the file-wide setting, which is what leaves a
!              single-rate site behaving exactly as before.
!
!              RP-side because the allowance is a RawProcess project setting;
!              a src_common file is compiled into FCC too, where RPsetup does
!              not exist.
! \sa          eliminate_corrupted_variables.f90
!***************************************************************************

!***************************************************************************
!
! \brief       The rate a column is sampled at [Hz]: its instrument's, or the
!              file's when the instrument states none.
! \note        The anemometer needs no special case - it is simply an
!              instrument whose rate is the file rate. A column with no
!              matched instrument carries NullInstrument, whose ac_freq is the
!              error sentinel, and so takes the file rate too.
!
!              Capped at the file's rate, because the rows are all there is: an
!              instrument declared faster than the file it was written into
!              would be expected to deliver more samples than the period has
!              rows, and every one of its columns would be dropped for missing
!              data that could not have been recorded. The setting says how
!              much SLOWER an instrument is.
!
!***************************************************************************
real(kind = dbl) function ColumnAcFreq(icol)
    use m_rp_global_var
    implicit none
    integer, intent(in) :: icol

    ColumnAcFreq = Metadata%ac_freq
    if (icol < 1 .or. icol > size(E2Col)) return
    if (E2Col(icol)%instr%ac_freq > 0d0) &
        ColumnAcFreq = min(E2Col(icol)%instr%ac_freq, Metadata%ac_freq)
end function ColumnAcFreq

!***************************************************************************
!
! \brief       The share of its OWN expected samples a column may be missing,
!              in percent.
! \note        Falls back to RPsetup%max_lack for a column whose instrument the
!              project says nothing about, and for one with no instrument at
!              all (slot 0). That fallback is what makes the project-wide
!              setting "the anemometer's" without naming the anemometer
!              anywhere.
!
!***************************************************************************
real(kind = dbl) function ColumnMaxLack(icol)
    use m_rp_global_var
    implicit none
    integer, intent(in) :: icol
    integer :: slot

    ColumnMaxLack = RPsetup%max_lack
    if (icol < 1 .or. icol > size(E2Col)) return

    slot = E2Col(icol)%instr%slot
    if (slot < 1 .or. slot > MaxNumInstruments) return
    if (RPsetup%instr_max_lack(slot) >= 0d0) &
        ColumnMaxLack = RPsetup%instr_max_lack(slot)
end function ColumnMaxLack

!***************************************************************************
!
! \brief       Where a slower column's real samples are, in the file's rows,
!              with a slot for each sample its instrument missed.
! \note        A slower instrument's samples do not sit at a fixed phase. The
!              Yatir laser runs at about 0.987 Hz against the 20 Hz sonic
!              clock, so its spacing is mostly 20 rows, often 21, and its
!              samples pass through every row position within a half-hour.
!              Anything that reads such a column at a fixed stride reads empty
!              rows. This walks the column instead and returns the rows its
!              values are actually on.
!
!              A gap wider than one and a half intervals is a missed sample,
!              and gets as many slots as whole intervals fit in it, spread
!              evenly; its value is left to the caller (usually the error
!              code, then filled). Nothing is inserted before the first sample
!              or after the last.
!
!              rows(1:ns) are the positions, isreal(1:ns) says which hold a
!              value, nreal counts them. Both arrays must hold nrow entries,
!              which is always enough: a slot is only ever inserted inside a
!              gap at least two rows wide.
! \sa          pwb_timelag_handle.f90 (PwbDetectSlowGas), spectral_analysis.f90
!
!***************************************************************************
subroutine SlowColumnSampleRows(col, nrow, stride, missing, rows, isreal, ns, nreal)
    use m_numeric_kinds
    implicit none
    integer, intent(in) :: nrow
    integer, intent(in) :: stride
    real(kind = dbl), intent(in) :: col(nrow)
    real(kind = dbl), intent(in) :: missing
    integer, intent(out) :: rows(nrow)
    logical, intent(out) :: isreal(nrow)
    integer, intent(out) :: ns
    integer, intent(out) :: nreal
    integer :: i
    integer :: j
    integer :: prev
    integer :: gap
    integer :: nmiss

    ns = 0
    nreal = 0
    prev = 0
    do i = 1, nrow
        if (col(i) == missing) cycle
        if (prev > 0) then
            gap = i - prev
            if (2 * gap > 3 * stride) then
                nmiss = nint(dble(gap) / dble(stride)) - 1
                do j = 1, nmiss
                    ns = ns + 1
                    rows(ns) = prev + nint(dble(j) * dble(gap) / dble(nmiss + 1))
                    isreal(ns) = .false.
                end do
            end if
        end if
        ns = ns + 1
        rows(ns) = i
        isreal(ns) = .true.
        nreal = nreal + 1
        prev = i
    end do
end subroutine SlowColumnSampleRows
