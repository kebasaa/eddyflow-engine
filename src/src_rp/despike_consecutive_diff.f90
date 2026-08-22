!***************************************************************************
! despike_consecutive_diff.f90
! ----------------------------
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
! \brief       Replace a sample that steps too far from the one before it.
! \author      Jonathan Muller
! \note
!              WHAT THIS IS. A rate-of-change limit, not a statistical
!              outlier test. If |x(k) - x(k-1)| exceeds a stated limit, x(k)
!              is replaced by x(k-1) and counted as a spike. Nothing here is
!              scaled by a standard deviation and nothing here iterates,
!              which is the whole difference between this and the two methods
!              beside it - Vickers & Mahrt walks a window widening a sigma
!              multiplier until no more spikes are found, and Mauder does the
!              same on a median absolute deviation. This looks only at
!              neighbours, in one pass, against a number the user states in
!              the variable's own units.
!
!              WHY IT IS HERE. EddyUH's spi_method 1
!              (EC_Software_Common/EddyUH_despike.m:66-76), and the method
!              the CH-LAE carbonyl sulfide project actually ran. Without it,
!              a run configured to match that project is despiked by a
!              different rule than the one that produced the fluxes it is
!              being compared against.
!
!              LIMITS. Zero means "not stated" and such a column is left
!              alone - the same shape as EddyUH's dlim, where a NaN entry
!              means the variable is not despiked. Which columns were skipped
!              is reported, because a method that silently does nothing when
!              nobody filled in its numbers is worse than one that refuses.
!
!              MISSING VALUES. EddyUH replaces a NaN by its predecessor and
!              COUNTS IT AS A SPIKE, before any method runs
!              (EddyUH_despike.m:39-60). This does not: a gap is a gap, the
!              engine has its own machinery for it, and folding gaps into the
!              spike count would make the spike percentage - which drives the
!              hard flag - a measure of two different things at once. An
!              error-coded sample is therefore skipped, and the comparison
!              resumes from the next valid pair.
! \sa          TestSpikeDetectionVickers97, TestSpikeDetectionMauder13
!***************************************************************************
subroutine DespikeConsecutiveDiff(Set, N, printout)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    logical, intent(in) :: printout
    real(kind = dbl), intent(inout) :: Set(N, E2NumVar)
    !> local variables
    integer :: i
    integer :: j
    integer :: prev
    integer :: nsample
    integer :: nspikes(E2NumVar)
    integer :: nreplaced(E2NumVar)
    integer :: hflags(u:lastGas)
    real(kind = dbl) :: step_lim(E2NumVar)
    real(kind = dbl), external :: ColumnAcFreq
    character(LongOutstringLen) :: skipped


    if (printout) write(*, '(a)') '   Spike detection/removal test, &
        &method: consecutive difference..'
    if (printout) write(ulog, '(a)') '   Spike detection/removal test, &
        &method: consecutive difference..'

    !> One limit per column, in that column's own units. The sonic components
    !> are named separately rather than sharing sr%lim_u's grouping: these are
    !> absolute limits, and a kelvin is not a metre per second.
    step_lim = 0d0
    step_lim(u)  = sr%step_u
    step_lim(v)  = sr%step_v
    step_lim(w)  = sr%step_w
    step_lim(ts) = sr%step_ts
    step_lim(firstGas:lastGas) = sr%step_gas(firstGas:lastGas)

    nspikes = 0
    nreplaced = 0
    skipped = ''
    do j = u, lastGas
        if (.not. E2Col(j)%present) cycle
        if (step_lim(j) <= 0d0) then
            !> Named, not just counted. The user has to know which column was
            !> passed over, or a half-filled settings page looks like a clean
            !> run.
            skipped = trim(skipped) // ' ' // trim(E2Col(j)%label)
            cycle
        end if

        !> prev is the last VALID sample, not simply i-1, so a gap does not
        !> manufacture a step out of the values either side of it.
        prev = 0
        do i = 1, N
            if (Set(i, j) == error) cycle
            if (prev /= 0) then
                if (dabs(Set(i, j) - Set(prev, j)) > step_lim(j)) then
                    nspikes(j) = nspikes(j) + 1
                    nreplaced(j) = nreplaced(j) + 1
                    !> Counted either way, replaced only when asked - the
                    !> same division the other two methods make, and the
                    !> tail below zeroes the replacement count exactly as
                    !> they do rather than inventing a second spelling.
                    if (RPsetup%filter_sr) Set(i, j) = Set(prev, j)
                end if
            end if
            prev = i
        end do
    end do

    if (len_trim(skipped) > 0) call LogSay('  No step limit stated for:' &
        // trim(skipped) // ' - left undespiked.')

    !> Hard flag on the same rule the other two methods use: the percentage is
    !> of the COLUMN's own samples, not of the file's rows, so a slower
    !> analyser is not diluted by its stride.
    hflags = 9
    do j = u, lastGas
        if (.not. E2Col(j)%present) cycle
        nsample = nint(dble(N) * ColumnAcFreq(j) / Metadata%ac_freq)
        if (nsample <= 0) nsample = N
        if (100d0 * (dble(nspikes(j)) / dble(nsample)) >= sr%hf_lim) then
            hflags(j) = 1
        else
            hflags(j) = 0
        end if
    end do

    !> Pack one digit per variable into the flag string
    call PackFlagString(hflags(u:lastGas), GHGNumVar, CharHF%sr)

    !> Write on output variable. Term for term what the other two methods do
    !> at their own tails - test_despike_consecutive_diff_static.py asserts
    !> the three agree, because a fourth spelling of "how a spike count
    !> reaches the output" is how one of them would quietly stop reporting.
    if (.not. RPsetup%filter_sr) nreplaced(u:pe) = 0
    where (E2Col(u:pe)%present)
        Essentials%e2spikes(u:pe) = nspikes(u:pe)
        Essentials%m_despiking(u:pe) = nreplaced(u:pe)
    elsewhere
        Essentials%e2spikes(u:pe) = ierror
        Essentials%m_despiking(u:pe) = ierror
    endwhere

    if (printout) write(*,'(a)') ' Done.'
    if (printout) write(ulog,'(a)') ' Done.'
end subroutine DespikeConsecutiveDiff
