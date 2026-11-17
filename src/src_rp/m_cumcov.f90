!***************************************************************************
! m_cumcov.f90
! ------------
! Copyright (C) 2026, ETH Zurich, Jonathan Muller
!
! This file is part of EddyFlow.
!
! EddyFlow (TM) is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! EddyFlow (TM) is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with EddyFlow (TM).  If not, see <http://www.gnu.org/licenses/>.
!
!***************************************************************************
!
! \brief       Time lag from cumulated cross-covariances, tlag_meth 6, after
!              innFLUX (Striednig et al. 2020, Atmos. Meas. Tech. 13,
!              1447-1465, Sect. 2.3.2; innFLUX_step2.m), who credit Park et
!              al. (2013) with cumulatively adding covariance functions; and
!              its optional clock-drift correction after innFLUX v1.2.0's
!              driftfix tool (2023, unpublished; M. Graus, M. Striednig).
!
!              A weak flux - most trace gases a PTR-MS measures - has no peak
!              a single half hour can find; picking the largest |cov| in each
!              period then picks noise, and biases the flux upwards. Summed
!              over many periods the noise averages away and the transport
!              lag, which a closed-path system holds nearly constant, stands
!              out. innFLUX sums |cov| over every period that passes four
!              filters - u*, the completeness of the sonic and of the gas, and
!              the steady-state test - inside user-chosen date segments, and
!              takes one lag per segment from the summed function, smoothed
!              with a Hamming window, as its largest value.
!
!              The time-lag pre-pass computes each period's cross-covariance
!              over the gas's search window (plus what the smoothing and the
!              drift need) and writes it, with the period's filter values, to
!              a scratch file per process. Nothing is summed there, so a
!              split pre-pass gives the same records as a serial one; the
!              parent reads them all back, orders them by period, filters and
!              sums. The production pass then applies, period by period, the
!              lag of the segment the period falls in, plus the drift.
!
!              Clock drift: the per-period lags (innFLUX's lagtime1) of the
!              gases marked as drift references are smoothed with a robust
!              LOWESS, split where the lag jumps (m_driftfix_core), and the
!              result is applied to every gas on the same analyser: each
!              period's function is shifted by the drift before summing, and
!              the drift is added back to the lag found (innFLUX's lagtime3).
!
! \author      Jonathan Muller
! \sa          m_lagpeak_core.f90, m_driftfix_core.f90, timelag_handle.f90
!***************************************************************************
module m_cumcov
    use m_rp_global_var
    use m_covmax_core, only: CrossCovarianceSeries
    use m_lagpeak_core, only: FindCovPeak, PeakPadding, FokenRelDev, &
        PeakPositive, PeakDerivative
    use m_driftfix_core, only: SegmentedDrift, DaysFromCivil, InterpLinear
    implicit none
    private

    public :: ParseCumCovSegments
    public :: CumCovBeginPrepass
    public :: CumCovGatherPeriod
    public :: CumCovFlushPeriod
    public :: CumCovEndSlice
    public :: CumCovResolve
    public :: CumCovPeriodLags
    public :: WriteCcfOutput

    !> Set by the reader when tlag_meth is 6.
    logical, public :: CumCovGenerate = .false.
    !> The lag of every variable for the period being processed, in seconds;
    !> error where the method has none and the nominal lag is used.
    real(kind = dbl), public :: CumCovLagNow(E2NumVar)

    !> innFLUX's steady-state test uses five sub-periods.
    integer, parameter :: SstSubPeriods = 5

    !> What one gas of one period contributes.
    type :: CumCovRec
        integer :: ord = 0
        character(10) :: date = ''
        character(5) :: time = ''
        integer :: gas = 0
        real(kind = dbl) :: ustar = -9999d0
        real(kind = dbl) :: csonic = 0d0
        real(kind = dbl) :: cgas = 0d0
        real(kind = dbl) :: sst = -9999d0
        integer :: lag1 = 0
        integer :: wmin = 0
        integer :: wmax = 0
        integer :: lo = 0
        integer :: hi = -1
        real(kind = sgl), allocatable :: ccf(:)
    end type CumCovRec

    integer :: uScratch = 0
    !> newunit= hands out negative numbers, so openness is its own flag.
    logical :: scratchOpen = .false.
    integer :: nPend = 0
    type(CumCovRec) :: pend(MaxNumGases)

    !> The resolved table the production pass reads.
    logical :: tableReady = .false.
    integer :: nPer = 0
    character(10), allocatable :: perDate(:)
    character(5), allocatable :: perTime(:)
    real(kind = dbl), allocatable :: perDay(:)
    real(kind = dbl), allocatable :: perLag(:, :)
    integer :: nSeg = 1
    real(kind = dbl) :: segBound(MaxCumSegments)
    real(kind = dbl), allocatable :: segLag(:, :)
    real(kind = dbl), allocatable :: perDrift(:, :)

contains

    !***************************************************************************
    !> \brief Segment boundaries from tlag_cum_segments: yyyy-mm-dd or
    !>        yyyy-mm-ddTHH:MM, separated by commas, blanks or quotes. Sorted.
    !***************************************************************************
    subroutine ParseCumCovSegments(text, bound, nmax, n)
        character(*), intent(in) :: text
        integer, intent(in) :: nmax
        real(kind = dbl), intent(out) :: bound(nmax)
        integer, intent(out) :: n
        character(len(text)) :: s
        character(32) :: tok
        integer :: i
        integer :: j
        integer :: k
        integer :: y
        integer :: mo
        integer :: d
        integer :: hh
        integer :: mi
        integer :: ios
        real(kind = dbl) :: t

        bound = error
        n = 0
        s = text
        do i = 1, len(s)
            if (s(i:i) == ',' .or. s(i:i) == '"' .or. s(i:i) == ';') s(i:i) = ' '
        end do
        i = 1
        do while (i <= len_trim(s))
            if (s(i:i) == ' ') then
                i = i + 1
                cycle
            end if
            j = index(s(i:), ' ')
            if (j == 0) then
                j = len_trim(s) + 1
            else
                j = i + j - 1
            end if
            tok = s(i:j - 1)
            i = j
            hh = 0
            mi = 0
            read(tok(1:4), '(i4)', iostat = ios) y
            if (ios /= 0) cycle
            read(tok(6:7), '(i2)', iostat = ios) mo
            if (ios /= 0) cycle
            read(tok(9:10), '(i2)', iostat = ios) d
            if (ios /= 0) cycle
            if (len_trim(tok) >= 16) then
                read(tok(12:13), '(i2)', iostat = ios) hh
                if (ios /= 0) hh = 0
                read(tok(15:16), '(i2)', iostat = ios) mi
                if (ios /= 0) mi = 0
            end if
            if (n >= nmax) exit
            n = n + 1
            bound(n) = DaysFromCivil(y, mo, d, hh, mi)
        end do
        !> Ascending, so a segment is the count of boundaries at or before t.
        do i = 2, n
            t = bound(i)
            k = i - 1
            do while (k >= 1)
                if (bound(k) <= t) exit
                bound(k + 1) = bound(k)
                k = k - 1
            end do
            bound(k + 1) = t
        end do
    end subroutine ParseCumCovSegments

    !> Days since 1970 of the START of the period whose end stamp is
    !> yyyy-mm-dd / HH:MM: innFLUX places segments and drift on the start time
    !> of each interval, and a period ending on a boundary belongs before it.
    function StampDays(date, time) result(days)
        character(*), intent(in) :: date
        character(*), intent(in) :: time
        real(kind = dbl) :: days
        integer :: y, mo, d, hh, mi, ios

        days = error
        read(date(1:4), '(i4)', iostat = ios) y
        if (ios /= 0) return
        read(date(6:7), '(i2)', iostat = ios) mo
        if (ios /= 0) return
        read(date(9:10), '(i2)', iostat = ios) d
        if (ios /= 0) return
        read(time(1:2), '(i2)', iostat = ios) hh
        if (ios /= 0) hh = 0
        read(time(4:5), '(i2)', iostat = ios) mi
        if (ios /= 0) mi = 0
        days = DaysFromCivil(y, mo, d, hh, mi) - dble(RPsetup%avrg_len) / 1440d0
    end function StampDays

    !> The segment a time (days) falls in, 1 .. nb + 1.
    pure integer function SegmentOf(t, bound, nb)
        real(kind = dbl), intent(in) :: t
        integer, intent(in) :: nb
        real(kind = dbl), intent(in) :: bound(:)
        integer :: k

        SegmentOf = 1
        do k = 1, nb
            if (t >= bound(k)) SegmentOf = k + 1
        end do
    end function SegmentOf

    !> This process's scratch file: piece 1 is the parent's.
    function ScratchPath(k) result(path)
        integer, intent(in) :: k
        character(PathLen) :: path
        character(8) :: tag

        write(tag, '(i2.2)') k
        path = trim(SharedDir()) // 'cumcov_b' // trim(tag) // '.bin'
    end function ScratchPath

    function TablePath() result(path)
        character(PathLen) :: path
        path = trim(SharedDir()) // 'cumcov_lags.bin'
    end function TablePath

    !> The parent's temporary folder, which every process of the run can see.
    !> A worker's own TmpDir is private to it; the parent hands it its own as
    !> --batch-tmp. Without this a worker's records never reached the parent,
    !> and a production worker never found the resolved lags.
    function SharedDir() result(dir)
        character(PathLen) :: dir

        if (BatchIndex > 0 .and. len_trim(BatchTmpDir) > 0) then
            dir = trim(BatchTmpDir) // slash
        else
            dir = TmpDir
        end if
    end function SharedDir

    !***************************************************************************
    !> \brief Open this process's scratch file at the start of the pre-pass.
    !***************************************************************************
    subroutine CumCovBeginPrepass()
        integer :: ios

        nPend = 0
        open(newunit = uScratch, file = trim(ScratchPath(max(BatchIndex, 1))), &
            form = 'unformatted', access = 'stream', status = 'replace', iostat = ios)
        scratchOpen = ios == 0
        if (.not. scratchOpen) then
            call LogSay(' Warning: the cumulated-covariance scratch file could not ' &
                // 'be opened; this piece contributes no periods.')
        end if
    end subroutine CumCovBeginPrepass

    !***************************************************************************
    !> \brief One period's cross-covariances and filter values, held until the
    !>        period's u* is known (CumCovFlushPeriod).
    !>
    !> Set is the series as covariance maximisation sees it, before the lag
    !> is applied; the gas windows are those AdjustTimelagOptSettings set.
    !***************************************************************************
    subroutine CumCovGatherPeriod(Set, nrow, ncol, ord, date, time)
        integer, intent(in) :: nrow
        integer, intent(in) :: ncol
        real(kind = dbl), intent(in) :: Set(nrow, ncol)
        integer, intent(in) :: ord
        character(*), intent(in) :: date
        character(*), intent(in) :: time
        real(kind = dbl), external :: ColumnAcFreq
        real(kind = dbl), allocatable :: cov(:)
        real(kind = dbl) :: expected
        integer :: gas
        integer :: pad
        integer :: wmin
        integer :: wmax
        integer :: lo
        integer :: hi
        integer :: lag1
        integer :: nw
        logical :: found

        nPend = 0
        if (.not. scratchOpen) return
        expected = Metadata%ac_freq * dble(RPsetup%avrg_len) * 60d0
        if (expected <= 0d0) return
        nw = count(Set(:, u) /= error .and. Set(:, v) /= error .and. Set(:, w) /= error)

        pad = PeakPadding(RPSetup%tlag_cum_filter_len, PeakDerivative)
        if (RPSetup%tlag_drift_meth == 1) &
            pad = pad + nint(RPSetup%tlag_drift_max_s * Metadata%ac_freq)

        do gas = firstGas, lastGas
            if (gas > ncol) exit
            if (.not. E2Col(gas)%present) cycle
            wmin = nint(E2Col(gas)%min_tl * Metadata%ac_freq)
            wmax = nint(E2Col(gas)%max_tl * Metadata%ac_freq)
            if (wmin == 0 .and. wmax == 0) cycle
            if (wmax < wmin) cycle
            lo = max(wmin - pad, -(nrow - 1))
            hi = min(wmax + pad, nrow - 1)
            if (hi < lo) cycle
            allocate(cov(lo:hi))
            call CrossCovarianceSeries(Set(:, w), Set(:, gas), nrow, lo, hi, error, cov)
            !> innFLUX's lagtime1: the period's own peak by its derivative rule.
            call FindCovPeak(cov, lo, hi, wmin, wmax, RPSetup%tlag_cum_filter_len, &
                PeakDerivative, .false., error, lag1, found)
            if (.not. found) lag1 = huge(1)

            nPend = nPend + 1
            pend(nPend)%ord = ord
            pend(nPend)%date = date
            pend(nPend)%time = time
            pend(nPend)%gas = gas
            pend(nPend)%csonic = dble(nw) / expected
            pend(nPend)%cgas = dble(count(Set(:, gas) /= error)) &
                / (expected * ColumnAcFreq(gas) / Metadata%ac_freq)
            if (found) then
                pend(nPend)%sst = FokenRelDev(Set(:, w), Set(:, gas), nrow, lag1, &
                    SstSubPeriods, error)
            else
                pend(nPend)%sst = error
            end if
            pend(nPend)%lag1 = lag1
            pend(nPend)%wmin = wmin
            pend(nPend)%wmax = wmax
            pend(nPend)%lo = lo
            pend(nPend)%hi = hi
            if (allocated(pend(nPend)%ccf)) deallocate(pend(nPend)%ccf)
            allocate(pend(nPend)%ccf(lo:hi))
            where (cov == error)
                pend(nPend)%ccf = real(error, sgl)
            elsewhere
                pend(nPend)%ccf = real(cov, sgl)
            end where
            deallocate(cov)
        end do
    end subroutine CumCovGatherPeriod

    !***************************************************************************
    !> \brief Write the held records once the period's u* is known. A period
    !>        skipped after CumCovGatherPeriod is simply never flushed.
    !***************************************************************************
    subroutine CumCovFlushPeriod(ustar)
        real(kind = dbl), intent(in) :: ustar
        integer :: k

        if (.not. scratchOpen) return
        do k = 1, nPend
            write(uScratch) pend(k)%ord, pend(k)%date, pend(k)%time, pend(k)%gas, &
                ustar, pend(k)%csonic, pend(k)%cgas, pend(k)%sst, pend(k)%lag1, &
                pend(k)%wmin, pend(k)%wmax, pend(k)%lo, pend(k)%hi
            write(uScratch) pend(k)%ccf
        end do
        nPend = 0
    end subroutine CumCovFlushPeriod

    subroutine CumCovEndSlice()
        if (scratchOpen) close(uScratch)
        scratchOpen = .false.
        nPend = 0
    end subroutine CumCovEndSlice

    !> Whether a record passes innFLUX's four filters.
    logical function Passes(r)
        type(CumCovRec), intent(in) :: r

        Passes = .false.
        if (r%lag1 == huge(1)) return
        if (r%ustar == error .or. r%ustar <= RPSetup%tlag_cum_min_ustar) return
        if (r%csonic <= RPSetup%tlag_cum_min_compl_sonic) return
        if (r%cgas <= RPSetup%tlag_cum_min_compl_gas) return
        if (r%sst == error .or. r%sst >= RPSetup%tlag_cum_max_sst) return
        Passes = .true.
    end function Passes

    !> Same analyser as innFLUX's "the same tracer data files": by name and model.
    logical function OnAnalyser(a, b)
        integer, intent(in) :: a
        integer, intent(in) :: b

        OnAnalyser = len_trim(E2Col(a)%instr_name) > 0 &
            .and. E2Col(a)%instr_name == E2Col(b)%instr_name &
            .and. E2Col(a)%instr%model == E2Col(b)%instr%model
    end function OnAnalyser

    !***************************************************************************
    !> \brief Read every piece's records back, sum, and resolve the lags.
    !***************************************************************************
    subroutine CumCovResolve(nChunks)
        integer, intent(in) :: nChunks
        type(CumCovRec), allocatable :: recs(:)
        type(CumCovRec), allocatable :: tmp(:)
        integer, allocatable :: order(:)
        integer :: nrec
        integer :: k
        integer :: u2
        integer :: ios
        integer :: i
        integer :: j
        integer :: gas
        integer :: g2
        integer :: s
        integer :: lag
        integer :: d
        integer :: L
        integer :: glo
        integer :: ghi
        integer :: gwmin
        integer :: gwmax
        integer :: lag2
        integer :: ppad
        integer :: maxd
        integer :: npt
        integer :: njump
        integer :: perDayCount
        integer :: nref
        logical :: found
        logical :: jumps_ok
        logical :: ok
        real(kind = dbl), allocatable :: sumAbs(:, :)
        real(kind = dbl), allocatable :: sumCov(:, :)
        integer, allocatable :: nUsed(:)
        real(kind = dbl), allocatable :: px(:)
        real(kind = dbl), allocatable :: py(:)
        real(kind = dbl), allocatable :: drift(:)
        integer, allocatable :: perOf(:)
        integer :: usum
        logical :: sumOpen
        real(kind = dbl) :: dmed
        logical :: done(E2NumVar)
        type(CumCovRec) :: r
        character(64) :: msg

        !> The files written here are plain CSV, header first, so every reader
        !> and check_columns.py takes them as they are; the method and its
        !> source are named here, in the run log, and in the documentation.
        call LogSay(' Resolving the cumulated-covariance time lags (innFLUX, ' &
            // 'Striednig et al. 2020, after Park et al. 2013):')

        !> Read every piece.
        nrec = 0
        allocate(recs(1024))
        do k = 1, max(1, nChunks)
            open(newunit = u2, file = trim(ScratchPath(k)), form = 'unformatted', &
                access = 'stream', status = 'old', iostat = ios)
            if (ios /= 0) cycle
            do
                read(u2, iostat = ios) r%ord, r%date, r%time, r%gas, r%ustar, &
                    r%csonic, r%cgas, r%sst, r%lag1, r%wmin, r%wmax, r%lo, r%hi
                if (ios /= 0) exit
                if (allocated(r%ccf)) deallocate(r%ccf)
                allocate(r%ccf(r%lo:r%hi))
                read(u2, iostat = ios) r%ccf
                if (ios /= 0) exit
                if (nrec >= size(recs)) then
                    allocate(tmp(2 * size(recs)))
                    do i = 1, nrec
                        call move_alloc(recs(i)%ccf, tmp(i)%ccf)
                        tmp(i)%ord = recs(i)%ord
                        tmp(i)%date = recs(i)%date
                        tmp(i)%time = recs(i)%time
                        tmp(i)%gas = recs(i)%gas
                        tmp(i)%ustar = recs(i)%ustar
                        tmp(i)%csonic = recs(i)%csonic
                        tmp(i)%cgas = recs(i)%cgas
                        tmp(i)%sst = recs(i)%sst
                        tmp(i)%lag1 = recs(i)%lag1
                        tmp(i)%wmin = recs(i)%wmin
                        tmp(i)%wmax = recs(i)%wmax
                        tmp(i)%lo = recs(i)%lo
                        tmp(i)%hi = recs(i)%hi
                    end do
                    call move_alloc(tmp, recs)
                end if
                nrec = nrec + 1
                recs(nrec) = r
            end do
            close(u2, status = 'delete')
        end do

        write(msg, '(i0, a, i0, a)') nrec, ' period-gas records from ', &
            max(1, nChunks), ' piece(s)'
        call LogSay('  ' // trim(msg))

        !> Period order, then gas order: the order a serial walk wrote them,
        !> so the sums are the same whatever the split.
        allocate(order(nrec))
        order = [(i, i = 1, nrec)]
        do i = 2, nrec
            j = order(i)
            k = i - 1
            do while (k >= 1)
                if (recs(order(k))%ord < recs(j)%ord) exit
                if (recs(order(k))%ord == recs(j)%ord .and. &
                    recs(order(k))%gas <= recs(j)%gas) exit
                order(k + 1) = order(k)
                k = k - 1
            end do
            order(k + 1) = j
        end do

        !> The periods, once each, in time.
        nPer = 0
        allocate(perOf(nrec))
        if (allocated(perDate)) deallocate(perDate, perTime, perDay)
        allocate(perDate(max(1, nrec)), perTime(max(1, nrec)), perDay(max(1, nrec)))
        do i = 1, nrec
            j = order(i)
            if (nPer == 0) then
                nPer = 1
            elseif (recs(j)%ord /= recs(order(i - 1))%ord) then
                nPer = nPer + 1
            end if
            perOf(j) = nPer
            perDate(nPer) = recs(j)%date
            perTime(nPer) = recs(j)%time
            perDay(nPer) = StampDays(recs(j)%date, recs(j)%time)
        end do

        nSeg = RPSetup%tlag_cum_nseg + 1
        segBound = RPSetup%tlag_cum_seg
        if (allocated(segLag)) deallocate(segLag)
        allocate(segLag(nSeg, E2NumVar))
        segLag = error
        if (allocated(perLag)) deallocate(perLag)
        allocate(perLag(max(1, nPer), E2NumVar))
        perLag = error
        if (allocated(perDrift)) deallocate(perDrift)
        allocate(perDrift(max(1, nPer), E2NumVar))
        perDrift = 0d0

        !> Clock drift, per analyser that has a reference gas.
        if (RPSetup%tlag_drift_meth == 1 .and. nPer > 0) then
            perDayCount = max(1, nint(1440d0 / dble(RPsetup%avrg_len)))
            maxd = nint(RPSetup%tlag_drift_max_s * Metadata%ac_freq)
            done = .false.
            allocate(drift(nPer))
            do gas = firstGas, lastGas
                if (.not. RPSetup%tlag_drift_ref(gas) .or. done(gas)) cycle
                !> Every reference gas on this analyser, pooled as innFLUX pools
                !> its selected tracers.
                nref = 0
                do g2 = firstGas, lastGas
                    if (RPSetup%tlag_drift_ref(g2) .and. OnAnalyser(gas, g2)) then
                        nref = nref + 1
                        done(g2) = .true.
                    end if
                end do
                npt = 0
                allocate(px(nrec), py(nrec))
                do i = 1, nrec
                    j = order(i)
                    if (.not. RPSetup%tlag_drift_ref(recs(j)%gas)) cycle
                    if (.not. OnAnalyser(gas, recs(j)%gas)) cycle
                    if (.not. Passes(recs(j))) cycle
                    npt = npt + 1
                    px(npt) = perDay(perOf(j))
                    py(npt) = dble(recs(j)%lag1) / Metadata%ac_freq
                end do
                if (npt > 0) then
                    call SegmentedDrift(px(1:npt), py(1:npt), npt, perDay(1:nPer), nPer, &
                        perDayCount * nref, &
                        max(3, nint(RPSetup%tlag_drift_span_d * perDayCount)) * nref, &
                        RPSetup%tlag_drift_jump_thr, error, drift, njump, jumps_ok)
                    !> innFLUX shifts every function by the whole smoothed lag,
                    !> so the summed peak lands near zero, and adds the drift
                    !> back afterwards. Here the drift is taken relative to its
                    !> median: the summed peak then stays where the transport
                    !> lag is, inside the gas's search window, and the result,
                    !> segment lag plus relative drift, is the same.
                    call MedianOf(drift, nPer, dmed)
                    do g2 = firstGas, lastGas
                        if (.not. OnAnalyser(gas, g2)) cycle
                        do i = 1, nPer
                            if (drift(i) /= error) perDrift(i, g2) = drift(i) - dmed
                        end do
                    end do
                    write(msg, '(i0, a, i0, a)') njump, ' jump(s), ', npt, ' lags'
                    call LogSay('  Clock drift for ' // trim(E2Col(gas)%instr_name) &
                        // ': ' // trim(msg))
                else
                    call LogSay('  Clock drift for ' // trim(E2Col(gas)%instr_name) &
                        // ': no reference lag passed the filters; not corrected.')
                end if
                deallocate(px, py)
            end do
            deallocate(drift)
            !> Bounded, because the windows were gathered with that much room.
            where (perDrift > RPSetup%tlag_drift_max_s) perDrift = RPSetup%tlag_drift_max_s
            where (perDrift < -RPSetup%tlag_drift_max_s) perDrift = -RPSetup%tlag_drift_max_s
        end if

        !> Sum and resolve, gas by gas.
        open(newunit = usum, file = trim(OutPath('cumcov_timelags')), &
            status = 'replace', iostat = ios, encoding = 'utf-8')
        sumOpen = ios == 0
        if (sumOpen) then
            write(usum, '(a)') 'gas,segment,from,to,n_periods,timelag_s'
        end if
        ppad = PeakPadding(RPSetup%tlag_cum_filter_len, PeakPositive)
        do gas = firstGas, lastGas
            !> The window of the first record of this gas fixes the grid.
            gwmin = 0
            gwmax = -1
            do i = 1, nrec
                if (recs(order(i))%gas == gas) then
                    gwmin = recs(order(i))%wmin
                    gwmax = recs(order(i))%wmax
                    exit
                end if
            end do
            if (gwmax < gwmin) cycle
            glo = gwmin - ppad
            ghi = gwmax + ppad
            allocate(sumAbs(glo:ghi, nSeg), sumCov(glo:ghi, nSeg), nUsed(nSeg))
            sumAbs = 0d0
            sumCov = 0d0
            nUsed = 0
            do i = 1, nrec
                j = order(i)
                if (recs(j)%gas /= gas) cycle
                if (recs(j)%wmin /= gwmin .or. recs(j)%wmax /= gwmax) cycle
                if (.not. Passes(recs(j))) cycle
                d = 0
                if (RPSetup%tlag_drift_meth == 1) &
                    d = nint(perDrift(perOf(j), gas) * Metadata%ac_freq)
                !> Shifted by the drift: ccf'(L) = ccf(L + d), innFLUX's
                !> circshift on its reversed axis.
                if (glo + d < recs(j)%lo .or. ghi + d > recs(j)%hi) cycle
                ok = .true.
                do L = glo, ghi
                    if (recs(j)%ccf(L + d) == real(error, sgl)) then
                        ok = .false.
                        exit
                    end if
                end do
                if (.not. ok) cycle
                s = SegmentOf(perDay(perOf(j)), segBound, nSeg - 1)
                do L = glo, ghi
                    sumAbs(L, s) = sumAbs(L, s) + abs(dble(recs(j)%ccf(L + d)))
                    sumCov(L, s) = sumCov(L, s) + dble(recs(j)%ccf(L + d))
                end do
                nUsed(s) = nUsed(s) + 1
            end do
            do s = 1, nSeg
                if (nUsed(s) == 0) cycle
                call FindCovPeak(sumAbs(:, s), glo, ghi, gwmin, gwmax, &
                    RPSetup%tlag_cum_filter_len, PeakPositive, .false., error, lag2, found)
                if (found) segLag(s, gas) = dble(lag2) / Metadata%ac_freq
            end do
            call WriteCumCovFunction(gas, glo, ghi, sumAbs, sumCov, nUsed)
            if (sumOpen) call WriteCumCovSegments(usum, gas, nUsed)
            deallocate(sumAbs, sumCov, nUsed)
        end do

        if (sumOpen) close(usum)

        !> Every period's lag: its segment's, plus its drift.
        do i = 1, nPer
            s = SegmentOf(perDay(i), segBound, nSeg - 1)
            do gas = firstGas, lastGas
                if (segLag(s, gas) == error) cycle
                perLag(i, gas) = segLag(s, gas) + perDrift(i, gas)
            end do
        end do
        if (RPSetup%tlag_drift_meth == 1) call WriteDriftTable()
        call WriteTable()
        tableReady = .true.
        call LogSay(' Done.')
    end subroutine CumCovResolve

    !> Output files are named like every other: <out>/<id>_<what>_<stamp>.csv
    function OutPath(what) result(path)
        character(*), intent(in) :: what
        character(PathLen) :: path

        path = Dir%main_out(1:len_trim(Dir%main_out)) &
            // EddyFlowProj%id(1:len_trim(EddyFlowProj%id)) &
            // '_' // trim(what) // Timestamp_FilePadding // CsvExt
    end function OutPath

    !> The gas's label, numbered like the FLUXNET columns when another gas
    !> before it has the same one: co2, then co2_2 for a second analyser's.
    function GasTag(gas) result(tag)
        integer, intent(in) :: gas
        character(32) :: tag
        character(32), external :: GasOutputLabel
        integer :: g
        integer :: k

        tag = GasOutputLabel(gas)
        if (len_trim(tag) == 0) write(tag, '(a, i0)') 'gas', gas - firstGas + 1
        k = 1
        do g = firstGas, gas - 1
            if (.not. E2Col(g)%present) cycle
            if (GasOutputLabel(g) == GasOutputLabel(gas)) k = k + 1
        end do
        if (k > 1) write(tag, '(a, a, i0)') trim(tag), '_', k
    end function GasTag

    !> The summed functions of one gas, per segment: lag (s), sum|cov|, sum cov.
    subroutine WriteCumCovFunction(gas, glo, ghi, sumAbs, sumCov, nUsed)
        integer, intent(in) :: gas
        integer, intent(in) :: glo
        integer, intent(in) :: ghi
        real(kind = dbl), intent(in) :: sumAbs(glo:ghi, nSeg)
        real(kind = dbl), intent(in) :: sumCov(glo:ghi, nSeg)
        integer, intent(in) :: nUsed(nSeg)
        integer :: uf
        integer :: ios
        integer :: L
        integer :: s
        character(32) :: c

        open(newunit = uf, file = trim(OutPath('cumcov_function_' // trim(GasTag(gas)))), &
            status = 'replace', iostat = ios, encoding = 'utf-8')
        if (ios /= 0) return
        write(uf, '(a)', advance = 'no') 'lag_s'
        do s = 1, nSeg
            write(c, '(i0)') s
            write(uf, '(a)', advance = 'no') ',sum_abs_cov_seg' // trim(c) // ',sum_cov_seg' // trim(c)
        end do
        write(uf, '(a)') ''
        write(uf, '(a)', advance = 'no') 'n_periods'
        do s = 1, nSeg
            write(uf, '(a, i0, a, i0)', advance = 'no') ',', nUsed(s), ',', nUsed(s)
        end do
        write(uf, '(a)') ''
        do L = glo, ghi
            write(uf, '(f0.4)', advance = 'no') dble(L) / Metadata%ac_freq
            do s = 1, nSeg
                write(uf, '(a, es15.7, a, es15.7)', advance = 'no') ',', sumAbs(L, s), &
                    ',', sumCov(L, s)
            end do
            write(uf, '(a)') ''
        end do
        close(uf)
    end subroutine WriteCumCovFunction

    !> One line per segment of this gas in the summary file.
    subroutine WriteCumCovSegments(us, gas, nUsed)
        integer, intent(in) :: us
        integer, intent(in) :: gas
        integer, intent(in) :: nUsed(nSeg)
        integer :: s
        character(16) :: a
        character(16) :: b

        do s = 1, nSeg
            a = 'start'
            b = 'end'
            if (s > 1) a = DayStamp(segBound(s - 1))
            if (s < nSeg) b = DayStamp(segBound(s))
            write(us, '(a, a, i0, a, a, a, a, a, i0, a, f0.4)') trim(GasTag(gas)), ',', s, &
                ',', trim(a), ',', trim(b), ',', nUsed(s), ',', segLag(s, gas)
        end do
    end subroutine WriteCumCovSegments

    !> yyyy-mm-ddTHH:MM of days since 1970 (Hinnant's civil_from_days).
    function DayStamp(days) result(stamp)
        real(kind = dbl), intent(in) :: days
        character(16) :: stamp
        integer :: z, era, doe, yoe, y, doy, mp, d, m, mins

        z = floor(days)
        mins = nint((days - dble(z)) * 1440d0)
        if (mins >= 1440) then
            z = z + 1
            mins = mins - 1440
        end if
        z = z + 719468
        if (z >= 0) then
            era = z / 146097
        else
            era = (z - 146096) / 146097
        end if
        doe = z - era * 146097
        yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365
        y = yoe + era * 400
        doy = doe - (365 * yoe + yoe / 4 - yoe / 100)
        mp = (5 * doy + 2) / 153
        d = doy - (153 * mp + 2) / 5 + 1
        if (mp < 10) then
            m = mp + 3
        else
            m = mp - 9
        end if
        if (m <= 2) y = y + 1
        write(stamp, '(i4.4, a, i2.2, a, i2.2, a, i2.2, a, i2.2)') y, '-', m, '-', d, &
            'T', mins / 60, ':', mod(mins, 60)
    end function DayStamp

    !> Median of the values of a(1:n) that are not error; 0 if none.
    subroutine MedianOf(a, n, med)
        integer, intent(in) :: n
        real(kind = dbl), intent(in) :: a(n)
        real(kind = dbl), intent(out) :: med
        real(kind = dbl) :: b(n)
        real(kind = dbl) :: t
        integer :: m
        integer :: i
        integer :: k

        m = 0
        do i = 1, n
            if (a(i) == error) cycle
            m = m + 1
            b(m) = a(i)
        end do
        med = 0d0
        if (m == 0) return
        do i = 2, m
            t = b(i)
            k = i - 1
            do while (k >= 1)
                if (b(k) <= t) exit
                b(k + 1) = b(k)
                k = k - 1
            end do
            b(k + 1) = t
        end do
        if (mod(m, 2) == 1) then
            med = b((m + 1) / 2)
        else
            med = 0.5d0 * (b(m / 2) + b(m / 2 + 1))
        end if
    end subroutine MedianOf

    subroutine WriteDriftTable()
        integer :: ud
        integer :: ios
        integer :: i
        integer :: gas

        open(newunit = ud, file = trim(OutPath('tlag_drift')), status = 'replace', &
            iostat = ios, encoding = 'utf-8')
        if (ios /= 0) return
        write(ud, '(a)', advance = 'no') 'date,time'
        do gas = firstGas, lastGas
            if (.not. E2Col(gas)%present) cycle
            write(ud, '(a)', advance = 'no') ',drift_s_' // trim(GasTag(gas))
        end do
        write(ud, '(a)') ''
        do i = 1, nPer
            write(ud, '(a)', advance = 'no') perDate(i) // ',' // perTime(i)
            do gas = firstGas, lastGas
                if (.not. E2Col(gas)%present) cycle
                write(ud, '(a, f0.4)', advance = 'no') ',', perDrift(i, gas)
            end do
            write(ud, '(a)') ''
        end do
        close(ud)
    end subroutine WriteDriftTable

    !> The table production workers read back.
    subroutine WriteTable()
        integer :: ut
        integer :: ios

        open(newunit = ut, file = trim(TablePath()), form = 'unformatted', &
            access = 'stream', status = 'replace', iostat = ios)
        if (ios /= 0) return
        write(ut) nPer, nSeg, E2NumVar
        write(ut) segBound
        write(ut) segLag
        if (nPer > 0) write(ut) perDate(1:nPer), perTime(1:nPer), perDay(1:nPer), &
            perLag(1:nPer, :), perDrift(1:nPer, :)
        close(ut)
    end subroutine WriteTable

    subroutine ReadTable()
        integer :: ut
        integer :: ios
        integer :: nv

        open(newunit = ut, file = trim(TablePath()), form = 'unformatted', &
            access = 'stream', status = 'old', iostat = ios)
        if (ios /= 0) return
        read(ut) nPer, nSeg, nv
        if (nv /= E2NumVar) then
            close(ut)
            return
        end if
        read(ut) segBound
        if (allocated(segLag)) deallocate(segLag)
        allocate(segLag(nSeg, E2NumVar))
        read(ut) segLag
        if (allocated(perDate)) deallocate(perDate, perTime, perDay, perLag, perDrift)
        allocate(perDate(max(1, nPer)), perTime(max(1, nPer)), perDay(max(1, nPer)))
        allocate(perLag(max(1, nPer), E2NumVar), perDrift(max(1, nPer), E2NumVar))
        if (nPer > 0) read(ut) perDate(1:nPer), perTime(1:nPer), perDay(1:nPer), &
            perLag(1:nPer, :), perDrift(1:nPer, :)
        close(ut)
        tableReady = .true.
    end subroutine ReadTable

    !***************************************************************************
    !> \brief The lags of the period being processed, into CumCovLagNow.
    !>
    !> A period the pre-pass saw takes its own entry. Any other takes its
    !> segment's lag plus the drift interpolated in time, and a gas with no
    !> lag in its segment is left at error, which TimeLagHandle answers with
    !> the nominal lag.
    !***************************************************************************
    subroutine CumCovPeriodLags(date, time)
        character(*), intent(in) :: date
        character(*), intent(in) :: time
        integer :: i
        integer :: gas
        integer :: s
        real(kind = dbl) :: t

        CumCovLagNow = error
        if (.not. tableReady) call ReadTable()
        if (.not. tableReady) return
        do i = 1, nPer
            if (perDate(i) == date(1:min(10, len(date))) &
                .and. perTime(i) == time(1:min(5, len(time)))) then
                CumCovLagNow(firstGas:lastGas) = perLag(i, firstGas:lastGas)
                return
            end if
        end do
        t = StampDays(date, time)
        if (t == error) return
        s = SegmentOf(t, segBound, nSeg - 1)
        do gas = firstGas, lastGas
            if (segLag(s, gas) == error) cycle
            CumCovLagNow(gas) = segLag(s, gas)
            if (nPer > 0) CumCovLagNow(gas) = CumCovLagNow(gas) &
                + InterpLinear(perDay(1:nPer), perDrift(1:nPer, gas), nPer, t)
        end do
    end subroutine CumCovPeriodLags

    !***************************************************************************
    !> \brief The period's cross-covariance function of every gas that asks
    !>        for it (gas_<i>_out_ccf), after innFLUX's step 1 covariance files
    !>        (Striednig et al. 2020). One file per period in eddyflow_ccf/:
    !>        lag in seconds (positive: the gas delayed) by gas, every
    !>        out_ccf_decim-th lag over +/- out_ccf_max_lag_s. Taken before the
    !>        time lag is applied, on the series covariance maximisation sees.
    !***************************************************************************
    subroutine WriteCcfOutput(Set, nrow, ncol, start_date, start_time)
        integer, intent(in) :: nrow
        integer, intent(in) :: ncol
        real(kind = dbl), intent(in) :: Set(nrow, ncol)
        character(*), intent(in) :: start_date
        character(*), intent(in) :: start_time
        real(kind = dbl), allocatable :: cov(:, :)
        real(kind = dbl) :: one(1)
        integer :: gas
        integer :: m
        integer :: step
        integer :: L
        integer :: k
        integer :: nlag
        integer :: uc
        integer :: ios
        integer :: ng
        integer :: sel(E2NumVar)
        character(16) :: stamp
        character(PathLen) :: path

        if (.not. any(RPSetup%out_ccf)) return
        ng = 0
        do gas = firstGas, lastGas
            if (gas > ncol) exit
            if (.not. RPSetup%out_ccf(gas) .or. .not. E2Col(gas)%present) cycle
            ng = ng + 1
            sel(ng) = gas
        end do
        if (ng == 0) return
        m = min(nint(RPSetup%out_ccf_max_lag_s * Metadata%ac_freq), nrow - 1)
        step = max(1, RPSetup%out_ccf_decim)
        nlag = 2 * (m / step) + 1
        allocate(cov(nlag, ng))
        do k = 1, ng
            do L = -(m / step), m / step
                call CrossCovarianceSeries(Set(:, w), Set(:, sel(k)), nrow, &
                    L * step, L * step, error, one)
                cov(L + m / step + 1, k) = one(1)
            end do
        end do

        stamp = start_date(1:4) // start_date(6:7) // start_date(9:10) // '-' &
            // start_time(1:2) // start_time(4:5)
        path = CcfDir(1:len_trim(CcfDir)) // trim(stamp) // '_ccf' &
            // Timestamp_FilePadding // CsvExt
        open(newunit = uc, file = trim(path), status = 'replace', iostat = ios, &
            encoding = 'utf-8')
        if (ios /= 0) then
            deallocate(cov)
            return
        end if
        write(uc, '(a)', advance = 'no') 'lag_s'
        do k = 1, ng
            write(uc, '(a)', advance = 'no') ',' // trim(GasTag(sel(k)))
        end do
        write(uc, '(a)') ''
        do L = -(m / step), m / step
            write(uc, '(f0.4)', advance = 'no') dble(L * step) / Metadata%ac_freq
            do k = 1, ng
                if (cov(L + m / step + 1, k) == error) then
                    write(uc, '(a)', advance = 'no') ',' // trim(EddyFlowProj%err_label)
                else
                    write(uc, '(a, es15.7)', advance = 'no') ',', cov(L + m / step + 1, k)
                end if
            end do
            write(uc, '(a)') ''
        end do
        close(uc)
        deallocate(cov)
    end subroutine WriteCcfOutput

end module m_cumcov
