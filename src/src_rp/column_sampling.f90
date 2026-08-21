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
