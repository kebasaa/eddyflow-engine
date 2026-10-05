!***************************************************************************
! pwb_stream.f90
! --------------
! Copyright © 2026, ETH Zurich, Jonathan Muller
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
! \brief       The PWB streaming classifier, in two halves: what a period's
!              own data say, and what the stream makes of it.
!
!              The production pass classifies each gas's PWB detection against
!              the last lag it settled on - S1, S2, an analyser neighbour's, a
!              carried one, or the covariance maximum. That memory reaches back
!              without limit, which is what kept this pass from being split.
!
!              It is split this way instead. PwbGatherEvidence does everything
!              that depends on the period alone - the cached row, the
!              detection, the covariance maximum - and none of it reads the
!              stream. PwbClassifyPeriod is the stream: it takes that evidence
!              and the per-gas memory, and settles the lags, exactly as the
!              single loop in TimeLagHandle did. A run that is not split calls
!              one after the other, as TimeLagHandle always did in one go.
!
!              A split run gathers evidence in its workers, in parallel - that
!              is the expensive half - and the parent feeds it through the
!              classifier in time order, which costs nothing. What comes out
!              for each period, a verdict, goes back to the worker that
!              computes that period's fluxes.
!
! \author      Jonathan Muller
! \sa          timelag_handle.f90, production_parallel.f90
!***************************************************************************
module m_pwb_stream
    use m_rp_global_var
    use m_pwb_timelag, only: InitPwbResult, GasLabel, SameAnalyser, &
        LookupPwbTimelagCache, StorePwbTimelagCache, CountPwbDiagnostic, &
        PwbDetectGas, SetPwbPeriodTimestamp
    implicit none
    private
    public :: PwbEvidenceType, PwbVerdictType
    public :: PwbGatherEvidence, PwbClassifyPeriod
    public :: PwbTakeVerdict, PwbApplyVerdict, PwbReplayEvidence
    public :: PwbEvidenceOnly, PwbLastEvidence
    public :: WritePwbEvidence, ReadPwbEvidence, WritePwbVerdict, ReadPwbVerdict

    !> Everything the classifier needs from one period, and nothing it
    !> remembers from the ones before.
    type :: PwbEvidenceType
        integer :: pcount = 0
        character(10) :: date = ''
        character(5) :: time = ''
        real(kind = dbl) :: ac_freq = 0d0
        logical :: present(E2NumVar) = .false.
        real(kind = dbl) :: def_tl(E2NumVar) = 0d0
        integer :: def_rl(E2NumVar) = 0
        !> A row of the time-lag table for this period, if there was one.
        logical :: cache_found(E2NumVar) = .false.
        real(kind = dbl) :: cache_actual(E2NumVar) = 0d0
        real(kind = dbl) :: cache_used(E2NumVar) = 0d0
        integer :: cache_row(E2NumVar) = 0
        logical :: cache_default(E2NumVar) = .false.
        !> The table's result where it had a row, the detection's otherwise.
        type(PWBResultType) :: res(E2NumVar)
        logical :: success(E2NumVar) = .false.
        !> The covariance maximum with the default on the window's edges.
        real(kind = dbl) :: mc_actual(E2NumVar) = 0d0
        real(kind = dbl) :: mc_used(E2NumVar) = 0d0
        integer :: mc_row(E2NumVar) = 0
        logical :: mc_default(E2NumVar) = .false.
        !> Which gases share an analyser, as E2Col says this period.
        logical :: same(E2NumVar, E2NumVar) = .false.
    end type PwbEvidenceType

    !> What the classifier settled for one period: everything the detection
    !> call leaves behind for the rest of the period to use. The lag arrays
    !> whole, since a gas a period lacks keeps whatever an earlier period left
    !> there; the results only for the gases present, since every other slot
    !> the call touches it resets - a parent of a year's run keeps one of these
    !> per half-hour, and the whole PWBResult array is some 25 kB.
    type :: PwbVerdictType
        integer :: pcount = 0
        integer :: row_lags(E2NumVar) = 0
        real(kind = dbl) :: act_tlag(E2NumVar) = 0d0
        real(kind = dbl) :: tlag(E2NumVar) = 0d0
        logical :: def_used(E2NumVar) = .false.
        integer :: own_row_lags(E2NumVar) = 0
        integer :: ngas = 0
        integer, allocatable :: gas(:)
        type(PWBResultType), allocatable :: res(:)
    end type PwbVerdictType

    !> Set by a worker that only gathers evidence: TimeLagHandle then leaves
    !> the period's evidence in PwbLastEvidence and classifies nothing.
    logical :: PwbEvidenceOnly = .false.
    type(PwbEvidenceType) :: PwbLastEvidence

contains

    !***************************************************************************
    !> \brief This period's evidence, for every present gas.
    !***************************************************************************
    subroutine PwbGatherEvidence(Set, nrow, ncol, def_rl, min_rl, max_rl, ev)
        integer, intent(in) :: nrow
        integer, intent(in) :: ncol
        real(kind = dbl), intent(in) :: Set(nrow, ncol)
        integer, intent(in) :: def_rl(ncol)
        integer, intent(in) :: min_rl(ncol)
        integer, intent(in) :: max_rl(ncol)
        type(PwbEvidenceType), intent(out) :: ev
        integer :: j
        integer :: k
        logical :: found
        external :: ApplyCovMaxDefaultFallback

        ev%date = PwbPeriodDate
        ev%time = PwbPeriodTime
        ev%ac_freq = Metadata%ac_freq
        do j = 1, E2NumVar
            call InitPwbResult(ev%res(j))
        end do
        ev%present(ts:pe) = E2Col(ts:pe)%present
        ev%def_tl(ts:pe) = E2Col(ts:pe)%def_tl
        ev%def_rl(ts:pe) = def_rl(ts:pe)

        do j = firstGas, lastGas
            if (.not. E2Col(j)%present) cycle
            call LookupPwbTimelagCache(j, found, ev%cache_actual(j), &
                ev%cache_used(j), ev%cache_row(j), ev%cache_default(j), ev%res(j))
            ev%cache_found(j) = found
            if (found) cycle
            call PwbDetectGas(Set, nrow, ncol, j, ev%res(j), ev%success(j))
            !> This period's covariance maximum, taken whether or not anything
            !> here needs it - see TimeLagHandle for why it must be the
            !> period's own. Pass 3 below falls back to the same call, with
            !> the same arguments, so this one value serves both.
            call ApplyCovMaxDefaultFallback(Set, nrow, ncol, j, .true., &
                def_rl(j), min_rl(j), max_rl(j), &
                ev%mc_actual(j), ev%mc_used(j), ev%mc_row(j), ev%mc_default(j))
            ev%res(j)%maxcov_lag = ev%mc_used(j)
        end do

        do j = firstGas, lastGas
            do k = firstGas, lastGas
                ev%same(j, k) = SameAnalyser(j, k)
            end do
        end do
    end subroutine PwbGatherEvidence

    !***************************************************************************
    !> \brief Settle this period's lags from its evidence and the stream's
    !>        memory, and move the memory on.
    !>
    !> Writes what the detection call always wrote: RowLags, PWBResult, the
    !> own-evidence row lags, the per-gas memory, ActTLag/TLag/DefTlagUsed - and
    !> adds the period to the diagnostics and the time-lag table.
    !***************************************************************************
    subroutine PwbClassifyPeriod(ev, ActTLag, TLag, DefTlagUsed)
        type(PwbEvidenceType), intent(in) :: ev
        real(kind = dbl), intent(inout) :: ActTLag(E2NumVar)
        real(kind = dbl), intent(inout) :: TLag(E2NumVar)
        logical, intent(inout) :: DefTlagUsed(E2NumVar)
        integer :: j
        integer :: k
        logical :: cache_hit(E2NumVar)
        logical :: pwb_success
        type(PWBResultType) :: lPwbResult
        logical, external :: GasSlotIsWater

        DefTlagUsed = .false.
        cache_hit = .false.
        do j = ts, pe
            call InitPwbResult(PWBResult(j))
        end do

        !> Pass 1: S1/S2 classification of each gas's own evidence
        do j = firstGas, lastGas
            if (.not. ev%present(j)) cycle
            if (ev%cache_found(j)) then
                lPwbResult = ev%res(j)
                cache_hit(j) = .true.
                PWBResult(j) = lPwbResult
                RowLags(j) = ev%cache_row(j)
                TLag(j) = ev%cache_used(j)
                ActTLag(j) = ev%cache_actual(j)
                DefTlagUsed(j) = ev%cache_default(j)
                pwb_raw_OwnRowLags(j) = ev%cache_row(j)
                if (trim(lPwbResult%reliability_class) == 'S1_optimal' .or. &
                    trim(lPwbResult%reliability_class) == 'S2_optimal' .or. &
                    trim(lPwbResult%reliability_class) == 'S4_instrument_shared') then
                    pwb_last_optimal_lag(j) = ev%cache_used(j)
                    pwb_last_optimal_origin(j) = lPwbResult%origin_gas
                    pwb_has_previous(j) = .true.
                end if
                cycle
            end if
            lPwbResult = ev%res(j)
            pwb_success = ev%success(j)

            !> The lag this period's own evidence gives, whatever the
            !> classifier below makes of it: the detection where it succeeded
            !> off the window edge, else the covariance maximum.
            if (pwb_success .and. .not. lPwbResult%edge_pinned) then
                pwb_raw_OwnRowLags(j) = lPwbResult%row_lag
            else
                pwb_raw_OwnRowLags(j) = ev%mc_row(j)
            end if

            if (pwb_success .and. .not. lPwbResult%edge_pinned) then
                if (lPwbResult%hdi_range < PWBSetup%hdi_thresh_s) then
                    lPwbResult%reliability_class = 'S1_optimal'
                    lPwbResult%fill_method = 'native'
                    RowLags(j) = lPwbResult%row_lag
                    TLag(j) = lPwbResult%selected_lag
                    ActTLag(j) = lPwbResult%selected_lag
                    DefTlagUsed(j) = .false.
                    pwb_last_optimal_lag(j) = lPwbResult%selected_lag
                    pwb_last_optimal_origin(j) = j
                    lPwbResult%origin_gas = j
                    pwb_has_previous(j) = .true.
                elseif (pwb_has_previous(j) .and. &
                    abs(lPwbResult%selected_lag - pwb_last_optimal_lag(j)) &
                    <= PWBSetup%dev_thresh_s) then
                    lPwbResult%reliability_class = 'S2_optimal'
                    lPwbResult%fill_method = 'native'
                    RowLags(j) = lPwbResult%row_lag
                    TLag(j) = lPwbResult%selected_lag
                    ActTLag(j) = lPwbResult%selected_lag
                    DefTlagUsed(j) = .false.
                    pwb_last_optimal_lag(j) = lPwbResult%selected_lag
                    pwb_last_optimal_origin(j) = j
                    lPwbResult%origin_gas = j
                    pwb_has_previous(j) = .true.
                else
                    lPwbResult%reliability_class = 'pending'
                end if
            else
                lPwbResult%reliability_class = 'pending'
            end if
            if (lPwbResult%applied_lag == error .and. &
                trim(lPwbResult%reliability_class) /= 'pending') then
                lPwbResult%applied_lag = TLag(j)
                lPwbResult%applied_row_lag = RowLags(j)
            end if
            PWBResult(j) = lPwbResult
        end do

        !> Pass 2: a gas with no detection of its own takes its analyser
        !> neighbour's - never water's. See TimeLagHandle.
        do j = firstGas, lastGas
            if (.not. ev%present(j)) cycle
            if (trim(PWBResult(j)%reliability_class) /= 'pending') cycle
            do k = firstGas, lastGas
                if (k == j) cycle
                if (.not. ev%present(k)) cycle
                if (GasSlotIsWater(k)) cycle
                if (.not. ev%same(j, k)) cycle
                if (trim(PWBResult(k)%reliability_class) /= 'S1_optimal' &
                    .and. trim(PWBResult(k)%reliability_class) /= 'S2_optimal') cycle
                PWBResult(j)%reliability_class = 'S4_instrument_shared'
                PWBResult(j)%fallback_used = .false.
                PWBResult(j)%fill_method = 'instrument_shared'
                PWBResult(j)%fallback_source = 'instrument_shared'
                PWBResult(j)%donor_gas = GasLabel(k)
                PWBResult(j)%origin_gas = merge(PWBResult(k)%origin_gas, k, &
                    PWBResult(k)%origin_gas > 0)
                PWBResult(j)%applied_lag = TLag(k)
                PWBResult(j)%applied_row_lag = RowLags(k)
                TLag(j) = TLag(k)
                RowLags(j) = RowLags(k)
                ActTLag(j) = ActTLag(k)
                DefTlagUsed(j) = .false.
                pwb_last_optimal_lag(j) = TLag(k)
                pwb_last_optimal_origin(j) = PWBResult(j)%origin_gas
                pwb_has_previous(j) = .true.
                exit
            end do
        end do

        !> Pass 3: S3 carry-forward or maxcov/default fallback for the rest
        do j = firstGas, lastGas
            if (.not. ev%present(j)) cycle
            if (trim(PWBResult(j)%reliability_class) /= 'pending') cycle
            if (pwb_has_previous(j)) then
                PWBResult(j)%reliability_class = 'S3_carryforward'
                PWBResult(j)%fill_method = 'carryforward'
                PWBResult(j)%fallback_source = 'S3_carryforward'
                PWBResult(j)%origin_gas = pwb_last_optimal_origin(j)
                if (PWBResult(j)%origin_gas > 0) &
                    PWBResult(j)%donor_gas = GasLabel(PWBResult(j)%origin_gas)
                TLag(j) = pwb_last_optimal_lag(j)
                if (PWBResult(j)%selected_lag /= error) then
                    ActTLag(j) = PWBResult(j)%selected_lag
                else
                    ActTLag(j) = pwb_last_optimal_lag(j)
                end if
                RowLags(j) = nint(pwb_last_optimal_lag(j) * ev%ac_freq)
                DefTlagUsed(j) = .false.
            else
                ActTLag(j) = ev%mc_actual(j)
                TLag(j) = ev%mc_used(j)
                RowLags(j) = ev%mc_row(j)
                DefTlagUsed(j) = ev%mc_default(j)
                PWBResult(j)%reliability_class = 'fallback'
                PWBResult(j)%fill_method = 'maxcov_default'
                PWBResult(j)%fallback_used = .true.
            end if
            if (PWBResult(j)%applied_lag == error) then
                PWBResult(j)%applied_lag = TLag(j)
                PWBResult(j)%applied_row_lag = RowLags(j)
            end if
        end do

        !> Finalize: fallback_source labels, diagnostics, the table
        do j = firstGas, lastGas
            if (.not. ev%present(j)) cycle
            if (PWBResult(j)%fallback_used .and. trim(PWBResult(j)%fallback_source) == 'none') &
                PWBResult(j)%fallback_source = 'maxcov_default'
            if (.not. PWBResult(j)%fallback_used .and. trim(PWBResult(j)%fallback_source) == 'none') &
                PWBResult(j)%fallback_source = 'native'
            if (trim(PWBResult(j)%fill_method) == 'none') PWBResult(j)%fill_method = 'native'
            if (.not. cache_hit(j)) then
                !> Counted here for the live path. A pre-generation run
                !> recounts from the settled table afterwards, so these
                !> streaming guesses never reach the summary.
                call CountPwbDiagnostic(j, PWBResult(j))
                call StorePwbTimelagCache(j, ActTLag(j), TLag(j), &
                    RowLags(j), DefTlagUsed(j), PWBResult(j))
            end if
        end do

        !> Non-gas scalars (ts, etc.) keep their nominal lags
        do j = ts, pe
            if (j >= firstGas .and. j <= lastGas) cycle
            if (ev%present(j)) then
                RowLags(j) = ev%def_rl(j)
                TLag(j) = ev%def_tl(j)
                ActTLag(j) = ev%def_tl(j)
                DefTlagUsed(j) = .true.
            else
                RowLags(j) = 0
                TLag(j) = 0d0
                ActTLag(j) = 0d0
            end if
        end do
    end subroutine PwbClassifyPeriod

    !***************************************************************************
    !> \brief Parent of a split run: classify a worker's evidence, in time
    !>        order, exactly as the detection call of a single pass would have.
    !>
    !> Each period's detection call starts from RowLags zeroed - the main loop
    !> does that before it - and leaves its lags in the pwb_raw arrays, which
    !> carry from one call to the next for gases a period does not have; the
    !> verdict takes all of it as it stands.
    !***************************************************************************
    subroutine PwbReplayEvidence(ev, v)
        type(PwbEvidenceType), intent(in) :: ev
        type(PwbVerdictType), intent(out) :: v

        call SetPwbPeriodTimestamp(ev%date, ev%time)
        RowLags = 0
        call PwbClassifyPeriod(ev, pwb_raw_ActTLag, pwb_raw_TLag, pwb_raw_DefTlagUsed)
        pwb_raw_Result = PWBResult
        call PwbTakeVerdict(ev%pcount, ev%present, v)
    end subroutine PwbReplayEvidence

    !> What the detection call just left behind, as a verdict for period p;
    !> present is which slots that call saw present.
    subroutine PwbTakeVerdict(p, present, v)
        integer, intent(in) :: p
        logical, intent(in) :: present(E2NumVar)
        type(PwbVerdictType), intent(out) :: v
        integer :: j

        v%pcount = p
        v%row_lags = RowLags
        v%act_tlag = pwb_raw_ActTLag
        v%tlag = pwb_raw_TLag
        v%def_used = pwb_raw_DefTlagUsed
        v%own_row_lags = pwb_raw_OwnRowLags
        v%ngas = count(present(firstGas:lastGas))
        allocate(v%gas(v%ngas), v%res(v%ngas))
        v%ngas = 0
        do j = firstGas, lastGas
            if (.not. present(j)) cycle
            v%ngas = v%ngas + 1
            v%gas(v%ngas) = j
            v%res(v%ngas) = PWBResult(j)
        end do
    end subroutine PwbTakeVerdict

    !> A worker of a split run: take the period's verdict in place of the
    !> detection call. PWBResult as that call leaves it: every slot from ts to
    !> pe reset, and the gases present classified.
    subroutine PwbApplyVerdict(v)
        type(PwbVerdictType), intent(in) :: v
        integer :: i
        integer :: j

        RowLags = v%row_lags
        pwb_raw_ActTLag = v%act_tlag
        pwb_raw_TLag = v%tlag
        pwb_raw_DefTlagUsed = v%def_used
        pwb_raw_OwnRowLags = v%own_row_lags
        do j = ts, pe
            call InitPwbResult(PWBResult(j))
        end do
        do i = 1, v%ngas
            PWBResult(v%gas(i)) = v%res(i)
        end do
    end subroutine PwbApplyVerdict

    !***************************************************************************
    !> \brief One period's evidence to a stream file, the present gases only.
    !***************************************************************************
    subroutine WritePwbEvidence(u, ev)
        integer, intent(in) :: u
        type(PwbEvidenceType), intent(in) :: ev
        integer :: j
        integer :: k
        integer :: n

        write(u) ev%pcount, ev%date, ev%time, ev%ac_freq
        write(u) ev%present, ev%def_tl, ev%def_rl
        n = count(ev%present(firstGas:lastGas))
        write(u) n
        do j = firstGas, lastGas
            if (.not. ev%present(j)) cycle
            write(u) j, ev%cache_found(j), ev%cache_actual(j), ev%cache_used(j), &
                ev%cache_row(j), ev%cache_default(j), ev%res(j), ev%success(j), &
                ev%mc_actual(j), ev%mc_used(j), ev%mc_row(j), ev%mc_default(j)
            do k = firstGas, lastGas
                if (ev%present(k)) write(u) ev%same(j, k)
            end do
        end do
    end subroutine WritePwbEvidence

    subroutine ReadPwbEvidence(u, ev)
        integer, intent(in) :: u
        type(PwbEvidenceType), intent(out) :: ev
        integer :: i
        integer :: j
        integer :: k
        integer :: n

        do j = 1, E2NumVar
            call InitPwbResult(ev%res(j))
        end do
        read(u) ev%pcount, ev%date, ev%time, ev%ac_freq
        read(u) ev%present, ev%def_tl, ev%def_rl
        read(u) n
        do i = 1, n
            read(u) j, ev%cache_found(j), ev%cache_actual(j), ev%cache_used(j), &
                ev%cache_row(j), ev%cache_default(j), ev%res(j), ev%success(j), &
                ev%mc_actual(j), ev%mc_used(j), ev%mc_row(j), ev%mc_default(j)
            do k = firstGas, lastGas
                if (ev%present(k)) read(u) ev%same(j, k)
            end do
        end do
    end subroutine ReadPwbEvidence

    subroutine WritePwbVerdict(u, v)
        integer, intent(in) :: u
        type(PwbVerdictType), intent(in) :: v

        write(u) v%pcount, v%row_lags, v%act_tlag, v%tlag, v%def_used, &
            v%own_row_lags, v%ngas
        if (v%ngas > 0) write(u) v%gas, v%res
    end subroutine WritePwbVerdict

    subroutine ReadPwbVerdict(u, v)
        integer, intent(in) :: u
        type(PwbVerdictType), intent(out) :: v

        read(u) v%pcount, v%row_lags, v%act_tlag, v%tlag, v%def_used, &
            v%own_row_lags, v%ngas
        allocate(v%gas(v%ngas), v%res(v%ngas))
        if (v%ngas > 0) read(u) v%gas, v%res
    end subroutine ReadPwbVerdict

end module m_pwb_stream
