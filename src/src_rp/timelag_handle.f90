!***************************************************************************
! timelag_handle.f90
! ------------------
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
! \brief       Calculates time lags (in terms of data rows) for all scalars \n
!              not measured by the anemometer. Also calculates covariances \n
!              of H2O and Cell T with time-lags of other scalars (from the \n
!              same instrument) for proper WPL of closed path systems.
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine TimeLagHandle(TlagMeth, Set, nrow, ncol, ActTLag, TLag, &
    DefTlagUsed, InTimelagOpt)
    use m_rp_global_var
    use m_pwb_timelag
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    character(*), intent(in) :: TlagMeth
    logical, intent(in) :: InTimelagOpt
    logical, intent(out) :: DefTlagUsed(ncol)
    real(kind = dbl), intent(out) :: ActTLag(ncol)
    real(kind = dbl), intent(out) :: TLag(ncol)
    real(kind = dbl), intent(inout) :: Set(nrow, ncol)
    !> local variables
    integer :: i = 0
    integer :: j = 0
    integer :: k = 0
    logical :: skip_apply
    logical :: cache_found
    logical :: cache_default_used
    logical :: cache_hit(E2NumVar)
    integer :: def_rl(ncol)
    integer :: min_rl(ncol)
    integer :: max_rl(ncol)
    real(kind = dbl) :: ColW(nrow)
    real(kind = dbl) :: ColH2O(nrow)
    real(kind = dbl) :: ColTC(nrow)
    real(kind = dbl) :: TmpSet(nrow, ncol)
    type(PWBResultType) :: lPwbResult
    logical :: pwb_success
    character(8) :: cache_stage
    real(kind = dbl) :: cache_actual_lag
    real(kind = dbl) :: cache_used_lag
    integer :: cache_row_lag
    !> The hygrometer that corrects the gas being handled, from that gas's
    !> own record. Was the site's water for every gas.
    integer :: msl
    !> Base slot of the cell block belonging to the analyser that measured the
    !> gas being handled: cellBase + 0..3 is cell_t, int_t_1, int_t_2, int_p.
    integer :: cellBase
    !> Row lag the water covariance is taken at: the gas's own where the two
    !> share an analyser, the hygrometer's where they do not.
    integer :: lagRow
    include '../src_common/interfaces_1.inc'

    skip_apply = pwb_detect_only_mode
    pwb_detect_only_mode = .false.

    if  (.not. InTimelagOpt .and. .not. skip_apply) write(*, '(a)', advance = 'no') &
        '  Compensating time-lags..'

    !> for E2Set scalars, initialise auxiliary vars to zero
    def_rl(:) = 0
    min_rl(:) = 0
    max_rl(:) = 0
    !> Define "row-lags" for scalars, using time-lags
    !> retrieved from metadata file
    where (E2Col(ts:pe)%present)
        def_rl(ts:pe) = nint(E2Col(ts:pe)%def_tl * Metadata%ac_freq)
        min_rl(ts:pe) = nint(E2Col(ts:pe)%min_tl * Metadata%ac_freq)
        max_rl(ts:pe) = nint(E2Col(ts:pe)%max_tl * Metadata%ac_freq)
    end where

    DefTlagUsed = .false.
    cache_hit = .false.
    do j = ts, pe
        call InitPwbResult(PWBResult(j))
    end do

    !> calculate actual time-lags according to the chosen method
    select case(TlagMeth)
        case ('constant')
            !> constant timelags are set equal to default values (user selected)
            RowLags(ts:pe) = def_rl(ts:pe)
            TLag(ts:pe)    = E2Col(ts:pe)%def_tl
            ActTLag(ts:pe) = E2Col(ts:pe)%def_tl
            DefTlagUsed(ts:pe) = .true.
        case ('maxcov', 'maxcov&default')
            !> covariance maximization method, with or without default
            do j = ts, pe
                !> Only for present variables,
                !> with both min and max "row lags" /= 0
                if (E2Col(j)%present &
                    .and. (min_rl(j) /= 0 .or. max_rl(j) /= 0)) then
                    call ApplyCovMaxDefaultFallback(Set, nrow, ncol, j, &
                        TlagMeth == 'maxcov&default', def_rl(j), &
                        min_rl(j), max_rl(j), ActTLag(j), TLag(j), &
                        RowLags(j), DefTlagUsed(j))
                else
                    RowLags(j) = 0
                    TLag(j) = 0d0
                    ActTLag(j) = 0d0
               end if
            end do
        case ('pwb')
            if (pwb_raw_detection_done .and. .not. skip_apply) then
                ActTLag = pwb_raw_ActTLag
                TLag = pwb_raw_TLag
                DefTlagUsed = pwb_raw_DefTlagUsed
                PWBResult = pwb_raw_Result
                pwb_raw_detection_done = .false.
            else
            if (skip_apply) then
                cache_stage = 'pre_wpl'
            else
                cache_stage = 'post_wpl'
            end if
            !> Pass 1: Run PWB detection and S1/S2 classification for all gases
            do j = firstGas, lastGas
                if (.not. E2Col(j)%present) cycle
                call LookupPwbTimelagCache(j, cache_stage, cache_found, cache_actual_lag, &
                    cache_used_lag, cache_row_lag, cache_default_used, lPwbResult)
                if (cache_found) then
                    cache_hit(j) = .true.
                    PWBResult(j) = lPwbResult
                    RowLags(j) = cache_row_lag
                    TLag(j) = cache_used_lag
                    ActTLag(j) = cache_actual_lag
                    DefTlagUsed(j) = cache_default_used
                    if (trim(lPwbResult%reliability_class) == 'S1_optimal' .or. &
                        trim(lPwbResult%reliability_class) == 'S2_optimal' .or. &
                        trim(lPwbResult%reliability_class) == 'S4_instrument_shared') then
                        pwb_last_optimal_lag(j) = cache_used_lag
                        pwb_last_optimal_origin(j) = lPwbResult%origin_gas
                        pwb_has_previous(j) = .true.
                    end if
                    cycle
                end if
                call PwbDetectGas(Set, nrow, ncol, j, lPwbResult, pwb_success)

                if (pwb_success .and. .not. lPwbResult%edge_pinned) then
                    if (lPwbResult%hdi_range < PWBSetup%hdi_thresh_s) then
                        lPwbResult%reliability_class = 'S1_optimal'
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

            !> Pass 2: Same-instrument lag sharing for gases that didn't get S1/S2
            do j = firstGas, lastGas
                if (.not. E2Col(j)%present) cycle
                if (trim(PWBResult(j)%reliability_class) /= 'pending') cycle
                do k = firstGas, lastGas
                    if (k == j) cycle
                    if (.not. E2Col(k)%present) cycle
                    if (E2Col(k)%instr%model /= E2Col(j)%instr%model) cycle
                    if (trim(PWBResult(k)%reliability_class) /= 'S1_optimal' &
                        .and. trim(PWBResult(k)%reliability_class) /= 'S2_optimal') cycle
                    PWBResult(j)%reliability_class = 'S4_instrument_shared'
                    PWBResult(j)%fallback_used = .false.
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

            !> Pass 3: S3 carry-forward or maxcov/default fallback for remaining gases
            do j = firstGas, lastGas
                if (.not. E2Col(j)%present) cycle
                if (trim(PWBResult(j)%reliability_class) /= 'pending') cycle
                if (pwb_has_previous(j)) then
                    PWBResult(j)%reliability_class = 'S3_carryforward'
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
                    RowLags(j) = nint(pwb_last_optimal_lag(j) * Metadata%ac_freq)
                    DefTlagUsed(j) = .false.
                else
                    call ApplyCovMaxDefaultFallback(Set, nrow, ncol, j, &
                        .true., def_rl(j), min_rl(j), max_rl(j), &
                        ActTLag(j), TLag(j), RowLags(j), DefTlagUsed(j))
                    PWBResult(j)%reliability_class = 'fallback'
                    PWBResult(j)%fallback_used = .true.
                end if
                if (PWBResult(j)%applied_lag == error) then
                    PWBResult(j)%applied_lag = TLag(j)
                    PWBResult(j)%applied_row_lag = RowLags(j)
                end if
            end do

            !> Finalize: set fallback_source labels and write diagnostics
            do j = firstGas, lastGas
                if (.not. E2Col(j)%present) cycle
                if (PWBResult(j)%fallback_used .and. trim(PWBResult(j)%fallback_source) == 'none') &
                    PWBResult(j)%fallback_source = 'maxcov_default'
                if (.not. PWBResult(j)%fallback_used .and. trim(PWBResult(j)%fallback_source) == 'none') &
                    PWBResult(j)%fallback_source = 'native'
                if (.not. cache_hit(j)) then
                    call WritePwbDiagnostic(j, PWBResult(j))
                    call StorePwbTimelagCache(j, cache_stage, ActTLag(j), TLag(j), &
                        RowLags(j), DefTlagUsed(j), PWBResult(j))
                end if
            end do

            !> Handle non-gas scalars (ts, etc.)
            do j = ts, pe
                if (j >= firstGas .and. j <= lastGas) cycle
                if (E2Col(j)%present) then
                    RowLags(j) = def_rl(j)
                    TLag(j) = E2Col(j)%def_tl
                    ActTLag(j) = E2Col(j)%def_tl
                    DefTlagUsed(j) = .true.
                else
                    RowLags(j) = 0
                    TLag(j) = 0d0
                    ActTLag(j) = 0d0
                end if
            end do
            end if  !> pwb_raw_detection_done bypass
        case ('none')
            !> not compensating for timelags
            RowLags(ts:pe) = 0
            TLag(ts:pe) = 0d0
    end select

    if (.not. skip_apply .and. .not. InTimelagOpt) then
        !> For closed path instruments, calculate H2O covariances
        !> for time-lags of other scalars from the same instrument
        !>
        !> One pass per configured gas, replacing three unrolled arms that
        !> named co2, ch4 and the fourth slot. A gas is skipped against its
        !> own hygrometer by identity rather than by position: that
        !> covariance is the water flux itself, not a cross term.
        !> Each gas against the water that corrects IT, not against the site's.
        !>
        !> This asked the primary hygrometer three questions on every gas's
        !> behalf: is the primary present, is the primary closed-path, and is
        !> this gas on the same analyser as the primary. A gas whose
        !> moist_ref is a second hygrometer therefore got no covariance at
        !> all, and so no water-flux term in its WPL correction - while its
        !> sigma and rho_w came from that second hygrometer. The two halves of
        !> one term disagreed about which water they meant.
        !>
        !> The same-model test was doing what moist_ref does properly:
        !> ResolveGasRef already prefers a hygrometer on the gas's own
        !> instrument. Asking the reference directly covers the case the model
        !> string cannot - two analysers of the same model - and drops the
        !> primary-only restriction with it.
        Stats%h2ocov_tl = error
        ColW(1:nrow) = Set(1:nrow, w)
        do j = firstGas, lastGas
            if (.not. E2Col(j)%present) cycle
            msl = E2Col(j)%moist_ref
            !> Out of slot range covers the biomet reference as well as an
            !> unresolved one, and declining is right for both: this is a
            !> covariance of w with a high-frequency water signal, and a
            !> half-hourly RH sensor has none. Fluxes0 then finds
            !> h2ocov_tl at error and leaves E_gas unset, which is the
            !> honest answer rather than a missing term to explain later.
            if (msl < firstGas .or. msl > lastGas) cycle
            if (j == msl) cycle
            if (.not. E2Col(msl)%present) cycle
            !> The covariance is only meaningful where the water is sampled
            !> through the same cell the gas is, which is what closed-path
            !> means here.
            if (E2Col(msl)%instr%path_type /= 'closed') cycle
            !> At whichever lag the water itself has.
            !>
            !> A gas and a hygrometer sharing an analyser share a tube, so the
            !> gas's own lag is the water's too. Down a different tube it is
            !> not, and this used to decline the pairing on that ground -
            !> which left the gas with no water-flux term at all while
            !> MoistTerms went on taking its sigma and rho_w from that same
            !> hygrometer. Two halves of one term disagreeing about which
            !> water they meant, which is the fault the paragraph above
            !> records for a second hygrometer, returning by another route.
            !>
            !> The quantity wanted is that hygrometer's own water flux, so it
            !> is evaluated at that hygrometer's lag. Borrowing across
            !> analysers is a compromise and ExceptionHandler(106) says so;
            !> declaring an H2O on the gas's own analyser makes ResolveGasRef
            !> prefer it and the two lags become one again.
            if (E2Col(j)%instr%model == E2Col(msl)%instr%model) then
                lagRow = RowLags(j)
            else
                lagRow = RowLags(msl)
            end if
            if (lagRow <= 0) cycle
            ColH2O(1:nrow) = Set(1:nrow, msl)
            call CovarianceW(ColW, ColH2O, size(ColW), &
                lagRow, Stats%h2ocov_tl(j))
        end do

        !> Cell temperature covariance, each gas against the cell it was
        !> measured in, at its own time lag.
        !>
        !> This read `tc` - firstCell, the *first* cell block - and admitted a
        !> gas only if its analyser matched that block's. One global cell made
        !> the two the same thing. With per-instrument blocks the first is
        !> whichever analyser happens to hold cell record one, so every gas on
        !> any other analyser failed the test, got no covariance, and lost the
        !> cell-temperature term of its WPL correction entirely - reported as
        !> H_CELL = -9999 while the analyser that owned record one had it.
        !>
        !> On CH-LAE the MIRO owns record one, so the LI-7200's CO2 and H2O
        !> were the ones going without. v7.2.5 does the same, its single cell
        !> record belonging to the MIRO; this is not a regression, it is the
        !> defect becoming expressible now that cells are per-instrument.
        !>
        !> cell_ref is resolved in DefineE2Set and is what AirAndCellParameters
        !> and PointByPointToMixingRatio already read. A project with one cell
        !> block falls back to firstCell and behaves exactly as before.
        Stats%tc_cov_tl = error
        ColW(1:nrow) = Set(1:nrow, w)
        do j = firstGas, lastGas
            if (.not. E2Col(j)%present) cycle
            !> Zero when no cell record belongs to this gas's analyser, in
            !> which case there is no cell temperature to correlate against -
            !> the first block belongs to a different instrument.
            cellBase = E2Col(j)%cell_ref
            if (cellBase < firstCell .or. cellBase > lastCell) cycle
            if (.not. E2Col(cellBase)%present) cycle
            if (RowLags(j) <= 0) cycle
            !> Inside the loop: the column is this gas's cell now, not one
            !> array filled once for every gas.
            ColTC(1:nrow) = Set(1:nrow, cellBase)
            call CovarianceW(ColW, ColTC, size(ColTC), &
                RowLags(j), Stats%tc_cov_tl(j))
        end do
    end if

    if (.not. skip_apply) then
        !> Align data according to relevant time-lags,
        !> filling remaining with error code.
        do j = u, pe
            if (E2Col(j)%present) then
                if (RowLags(j) >= 0) then
                    !> For positive lags
                    do i = 1, nrow - RowLags(j)
                        TmpSet(i, j) = Set(i + RowLags(j), j)
                    end do
                    do i = nrow - Rowlags(j) + 1, nrow
                        TmpSet(i, j) = error
                    end do
                else
                    !> For negative lags
                    do i = 1, abs(RowLags(j))
                        TmpSet(i, j) = error
                    end do
                    do i = abs(RowLags(j)) + 1, nrow
                        TmpSet(i, j) = Set(i + RowLags(j), j)
                    end do
                end if
            else
                TmpSet(1:nrow, j) = error
            end if
        end do
        Set = TmpSet
        if  (.not. InTimelagOpt) write(*,'(a)') ' Done.'
    end if
end subroutine TimeLagHandle

subroutine ApplyCovMaxDefaultFallback(Set, nrow, ncol, gas, use_default_on_edge, &
    def_rl, min_rl, max_rl, actual_tlag, used_tlag, row_lag, def_tlag_used)
    use m_rp_global_var
    implicit none
    integer, intent(in) :: nrow, ncol, gas
    integer, intent(in) :: def_rl, min_rl, max_rl
    logical, intent(in) :: use_default_on_edge
    real(kind = dbl), intent(in) :: Set(nrow, ncol)
    real(kind = dbl), intent(out) :: actual_tlag
    real(kind = dbl), intent(out) :: used_tlag
    integer, intent(out) :: row_lag
    logical, intent(out) :: def_tlag_used
    real(kind = dbl) :: FirstCol(nrow)
    real(kind = dbl) :: SecondCol(nrow)

    FirstCol(:)  = Set(:, RPSetup%covmax_var)
    SecondCol(:) = Set(:, gas)
    call CovMax(min_rl, max_rl, FirstCol, SecondCol, size(FirstCol), &
        actual_tlag, row_lag)
    used_tlag = actual_tlag
    def_tlag_used = .false.
    if (use_default_on_edge .and. ((row_lag == min_rl) .or. (row_lag == max_rl))) then
        def_tlag_used = .true.
        used_tlag = dble(def_rl) / Metadata%ac_freq
        row_lag = def_rl
    end if
end subroutine ApplyCovMaxDefaultFallback

!*******************************************************************************
!
! \brief       Performs covariance analysis for determining the "optimal" \n
!              time lag, the one that maximizes the covariance.
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!*******************************************************************************
subroutine CovMax(lagmin, lagmax, Col1, Col2, nrow, TLag, RLag)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow
    integer, intent(in) :: lagmin
    integer, intent(in) :: lagmax
    real(kind = dbl), intent(in) :: Col1(nrow)
    real(kind = dbl), intent(in) :: Col2(nrow)
    integer, intent(out) :: RLag
    real(kind = dbl), intent(out) :: TLag
    !> local variables
    integer :: i = 0
    integer :: ii = 0
    integer :: N2
    real(kind = dbl), allocatable :: ShSet(:, :)
    real(kind = dbl), allocatable :: ShPrimes(:, :)
    real(kind = dbl) :: CovMat(2,2)
    real(kind = dbl) :: Cov
    real(kind = dbl) :: MaxCov

    Cov = 0.d0
    MaxCov = 0.d0
    TLag = 0.d0
    do i = lagmin, lagmax
        N2 = nrow - abs(i)
        allocate(ShSet(N2, 2))
        allocate(ShPrimes(N2, 2))

        !> Align the two timeseries at the current time-lag 
        do ii = 1, N2
            if (i < 0) then
                ShSet(ii, 1) = Col1(ii - i)
                ShSet(ii, 2) = Col2(ii)
            else
                ShSet(ii, 1) = Col1(ii)
                ShSet(ii, 2) = Col2(ii + i)
            end if
        end do


        !> Linear detrending
        ! call VariableLinearDetrending(ShSet(:, 1), ShPrimes(:, 1), N2)
        ! call VariableLinearDetrending(ShSet(:, 2), ShPrimes(:, 2), N2)
        if (RPSetup%covmax_stocdet) then
            !> Stochastic detrending
            call VariableStochasticDetrending(ShSet(:, 1), ShPrimes(:, 1), N2)
            call VariableStochasticDetrending(ShSet(:, 2), ShPrimes(:, 2), N2)
        else
            !> Block average
            ShPrimes = ShSet
        end if

        call CovarianceMatrixNoError(ShPrimes, size(ShPrimes, 1), size(ShPrimes, 2), CovMat, error)
        Cov = CovMat(1, 2)

        !> Max cov and actual time lag
        if (abs(Cov) > MaxCov) then
            MaxCov = abs(Cov)
            TLag = dble(i) / Metadata%ac_freq
            RLag = i
        end if
        deallocate(ShSet)
        deallocate(ShPrimes)
    end do
end subroutine CovMax


!***************************************************************************
!
! \brief       Calculate covariance between two arrays using an imposed  \n
!              time-lag.
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine CovarianceW(col1, col2, nrow, lag, cov)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow
    integer, intent(in) :: lag
    real(kind = dbl), intent(in) :: col1(nrow)
    real(kind = dbl), intent(in) :: col2(nrow)
    real(kind = dbl), intent(out) :: cov
    !> local variables
    integer :: i
    integer :: N2
    real(kind = dbl) ::sum1
    real(kind = dbl) ::sum2

    sum1 = 0d0
    sum2 = 0d0
    Cov = 0d0
    N2 = 0
    do i = 1, nrow - lag
        if (col1(i) /= error .and. col2(i+lag) /= error) then
            N2 = N2 + 1
            Cov = Cov + col1(i) * col2(i+lag)
            sum1 = sum1 + col1(i)
            sum2 = sum2 + col2(i+lag)
        end if
    end do

    if (N2 /= 0) then
        sum1 = sum1 / dble(N2)
        sum2 = sum2 / dble(N2)
        cov = cov / dble(N2)
        cov = cov - sum1 * sum2
    else
        cov = error
    end if
end subroutine CovarianceW

!***************************************************************************
!
! \brief       Stochastic Detrending
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine VariableStochasticDetrending(Var, Primes, N)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    real(kind = dbl), intent(in) :: Var(N)
    real(kind = dbl), intent(out) :: Primes(N)
    !> local variables
    integer :: i

    Primes(1) = error
    do i = 2, N
        if (Var(i) /= error .and. Var(i-1) /= error) then
            Primes(i) = Var(i) - Var(i-1)
        else 
            Primes(i) = error
        end if
    end do
end subroutine VariableStochasticDetrending

!***************************************************************************
!
! \brief       Linear detrending of one time series
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine VariableLinearDetrending(Var, Primes, N)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    real(kind = dbl), intent(in) :: Var(N)
    real(kind = dbl), intent(out) :: Primes(N)
    !> Local variables
    real(kind = dbl) :: Trend(N)

    call CalculateTrend(Var, Trend, N)
    call Detrend(Var, Trend, Primes, N)

end subroutine VariableLinearDetrending

!***************************************************************************
!
! \brief       Remove trend from time series
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine Detrend(Var, Trend, Primes, N)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    real(kind = dbl), intent(in) :: Var(N)
    real(kind = dbl), intent(in) :: Trend(N)
    real(kind = dbl), intent(out) :: Primes(N)


    Primes = error
    where (Var /= error .and. Trend /= error)
        Primes = Var - Trend
    end where
end subroutine Detrend

!***************************************************************************
!
! \brief       Calculate linear trend in time series
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine CalculateTrend(Var, Trend, N)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: N
    real(kind = dbl), intent(in) :: Var(N)
    real(kind = dbl), intent(out) :: Trend(N)
    !> local variables
    integer :: i
    integer :: nn
    integer :: mm
    real(kind = dbl) :: sumx1
    real(kind = dbl) :: sumx2
    real(kind = dbl) :: mean = 0d0
    real(kind = dbl) :: sumtime
    real(kind = dbl) :: sumtime2
    real(kind = dbl) :: b


    !> Linear regression
    sumx1 = 0d0
    sumx2 = 0d0
    sumtime = 0d0
    sumtime2 = 0d0
    nn = 0
    do i = 1, N
        if (Var(i) /= error) then
            nn = nn + 1
            sumx1 = sumx1 + (Var(i) * (dble(nn - 1)))
            sumx2 = sumx2 + Var(i)
            sumtime = sumtime + (dble(nn - 1))
            sumtime2 = sumtime2 + (dble(nn - 1))**2
        end if
    end do
    if (nn /= 0) then
        mean = sumx2 / dble(nn)
    end if

    !> Trend
    mm = 0
    b = (sumx1 - (sumx2 * sumtime) / dble(nn)) / (sumtime2 - (sumtime * sumtime) / dble(nn))
    do i = 1, N
        mm = mm + 1
        if (Var(i) /= error) then
            Trend(i) = mean + b * (dble(mm - 1) - sumtime / dble(nn))
        else
            Trend(i) = error
        end if
    end do        
end subroutine CalculateTrend
