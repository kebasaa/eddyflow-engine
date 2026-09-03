!***************************************************************************
! borrow_timelag.f90
! ------------------
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
! \brief       Take a tube-mate's time lag where a gas cannot resolve its own,
!              after Nemitz et al. (2018).
! \author      Jonathan Muller
!
! \note
! Gases drawn down one tube share a transport delay. A species whose own
! cross-covariance peak cannot be told from noise has nothing to detect, and
! whatever the maximisation returned for it is close to a coin toss over the
! search window; a species beside it in the same tube that *can* resolve its
! peak is measuring the delay they both have.
!
! Two things trigger it, both from EddyUH's EC_Software_FluxCalc/
! EddyUH_SC_Flux2.m:324: the covariance at the gas's own lag failing to clear
! tlag_borrow_snr times its detection limit - the paper's three - or the lag
! landing on an end of the search window, which is where a maximisation goes
! when there is no interior peak to find.
!
! WHICH NOISE FLOOR, AND WHICH DONOR, ARE BOTH CHOICES - and the defaults are
! not EddyUH's.
!
! EddyUH compares against unc3(:,3), Lenschow et al. (2000) instrumental noise
! propagated to a covariance as sqrt(sigma2_noise * sigma2_w / N), and borrows
! specifically from the analyser's carbon dioxide, hard-coded by variable name.
! The default here is the Wienhold et al. (1994) detection limit - the scatter
! of the very function whose peak is being tested, which is the question being
! asked - and the best-resolved tube-mate. tlag_borrow_noise and
! tlag_borrow_donor select EddyUH's instead.
!
! Read EddyUH's own comment at that line with care: it calls unc3 a detection
! limit, and it is not one. Both quantities are noise floors in covariance
! units, but the Lenschow estimate measures the ANALYSER and the Wienhold one
! measures the COVARIANCE. The first port of this took the comment at its word
! and divided by the wrong thing.
!
! Neither of EddyUH's two routines for this is reachable in the shipped
! package: EddyUH_SC_Flux2.m is called only from a commented-out block, and
! the live EddyUH_SC_Flux.m copies the carbon dioxide lag onto carbonyl
! sulfide unconditionally, with no test at all.
!
! \sa          pwb_timelag_handle.f90, whose S4_instrument_shared pass also
!              borrows across tube-mates. A different question, kept separate:
!              that one asks whether a bootstrap could not make up its mind
!              across periods, this one whether a peak stands above the noise
!              in this one.
!***************************************************************************
subroutine BorrowTimelagBelowDetectionLimit(Set, nrow, ncol, min_rl, max_rl, &
    TLag, DefTlagUsed)
    use m_rp_global_var
    use m_pwb_timelag, only: SameAnalyser
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    real(kind = dbl), intent(in) :: Set(nrow, ncol)
    integer, intent(in) :: min_rl(ncol)
    integer, intent(in) :: max_rl(ncol)
    !> Passed rather than reached through Essentials, which these alias: the
    !> caller holds them as dummy arguments, so writing the globals from here
    !> would be modifying a dummy through a second path.
    !>
    !> ActTLag is deliberately absent. It reports what each gas own
    !> maximisation found, and borrowing must not overwrite that - the two
    !> differing is the record that a lag was taken from elsewhere.
    real(kind = dbl), intent(inout) :: TLag(ncol)
    logical, intent(inout) :: DefTlagUsed(ncol)
    !> local variables
    logical, external :: GasSlotIsWater
    integer, external :: CarbonOnAnalyserOf
    real(kind = dbl), external :: LenschowFluxNoise
    integer :: gas, donor, chosen
    !> Heap, not automatic arrays. At 20 Hz over half an hour each of these
    !> is some 288 kB, and this runs from deep inside the period loop with
    !> TimeLagHandle's own frame already on the stack. Left automatic, and
    !> with a second one allocated per call inside the signal-to-noise
    !> helper, the pair overran and corrupted a neighbouring value: one gas's
    !> published time lag came out as a denormal, in periods where nothing
    !> had been borrowed for it at all.
    real(kind = dbl), allocatable :: ColW(:)
    real(kind = dbl), allocatable :: ColGas(:)
    real(kind = dbl) :: cov
    !> Signal-to-noise at each gas's own lag: the covariance there over its
    !> detection limit. Computed once, and used twice - to decide who can
    !> keep their lag, and to rank who is worth borrowing from.
    real(kind = dbl) :: snr(E2NumVar)
    real(kind = dbl) :: best
    logical :: trusted(E2NumVar)

    if (.not. RPSetup%tlag_borrow_meth) return
    !> Nothing to compare a covariance against otherwise. The interface greys
    !> the control for the same reason; this is the engine saying so for a
    !> hand-written project.
    !>
    !> Only the detection limit has to be switched on separately - it is
    !> computed elsewhere and read from Essentials. The Lenschow noise is
    !> measured here, from the series in hand, so choosing it needs nothing
    !> else enabled. That asymmetry is deliberate: it means the EddyUH-faithful
    !> combination can be selected on its own.
    if (RPSetup%tlag_borrow_noise == 'detlim' .and. &
        RPSetup%detlim_meth /= 'wienhold_94') return

    allocate(ColW(nrow))
    allocate(ColGas(nrow))
    ColW(1:nrow) = Set(1:nrow, w)

    !> Who may be borrowed from, decided before anything is borrowed.
    !>
    !> A gas that took its neighbour's lag is weaker evidence than one that
    !> found its own, and letting a borrowed lag be borrowed again would walk
    !> a single detection along a whole tube. Snapshotting the trusted set
    !> first is what stops that, and it is the rule the block-bootstrap
    !> borrowing follows too.
    snr = error
    trusted = .false.
    do gas = firstGas, lastGas
        if (.not. Eligible(gas)) cycle
        snr(gas) = SignalToNoise(gas)
        if (DefTlagUsed(gas)) cycle
        if (snr(gas) == error) cycle
        trusted(gas) = snr(gas) >= RPSetup%tlag_borrow_snr
    end do

    do gas = firstGas, lastGas
        if (.not. Eligible(gas)) cycle
        if (trusted(gas)) cycle

        !> A gas already carrying the nominal default is a candidate too: a
        !> tube-mate's measured delay beats a number typed into the metadata.
        if (.not. (DefTlagUsed(gas) .or. OnWindowEdge(gas) &
                   .or. snr(gas) == error &
                   .or. snr(gas) < RPSetup%tlag_borrow_snr)) cycle

        !> Who to borrow from. Two rules, and the default is not EddyUH's.
        !>
        !> 'best_resolved' ranks the trusted tube-mates by signal-to-noise.
        !> Taking the first in slot order - which is what a naive loop does -
        !> put whichever gas happened to sit lowest in the metadata in charge,
        !> and on a tube where nothing clears the threshold comfortably that
        !> produced pairings that read backwards: carbon dioxide, the
        !> strongest flux on the analyser, taking its lag from nitrous oxide.
        !>
        !> 'carbon_dioxide' is what EddyUH does, hard-coded by variable name
        !> at EddyUH_SC_Flux2.m:325-331 with no user switch. It is the right
        !> answer wherever carbon dioxide IS the best-resolved channel, which
        !> on a trace-gas analyser it usually is, and it has the merit of
        !> being the same choice in every period - a lag population that does
        !> not change donor halfway through the day is easier to defend. It
        !> refuses rather than falling back when that analyser measures no
        !> carbon dioxide, or when its carbon dioxide is not itself trusted.
        chosen = 0
        if (RPSetup%tlag_borrow_donor == 'carbon_dioxide') then
            donor = CarbonOnAnalyserOf(gas)
            if (donor >= firstGas .and. donor <= lastGas) then
                if (donor /= gas .and. trusted(donor)) chosen = donor
            end if
        else
            best = 0d0
            do donor = firstGas, lastGas
                if (donor == gas) cycle
                if (.not. trusted(donor)) cycle
                if (.not. SameAnalyser(gas, donor)) cycle
                if (snr(donor) > best) then
                    best = snr(donor)
                    chosen = donor
                end if
            end do
        end if
        if (chosen == 0) cycle

        RowLags(gas) = RowLags(chosen)
        TLag(gas) = TLag(chosen)
        !> ActTLag is left alone. It reports what this gas's own maximisation
        !> found, TLag what was applied, and the two differing is the record
        !> that this happened.
        !> Not its own detection, which is what this flag means everywhere
        !> else. It does not distinguish borrowed from nominal; the log line
        !> names the donor.
        DefTlagUsed(gas) = .true.

        call LogSay('  Time lag for ' // trim(E2Col(gas)%label) &
            // ' taken from ' // trim(E2Col(chosen)%label) &
            // ': its own covariance does not clear the detection limit.')
    end do

    deallocate(ColW)
    deallocate(ColGas)

contains

    !> A closed-path gas that is not water.
    !>
    !> Water is excluded on both sides. Its lag is the one every other gas's
    !> water covariance is taken at, so moving it moves those too; and it is
    !> rarely the weak signal here - at the site this was written for it is
    !> the strongest.
    logical function Eligible(g)
        integer, intent(in) :: g

        Eligible = .false.
        if (.not. E2Col(g)%present) return
        if (GasSlotIsWater(g)) return
        if (E2Col(g)%instr%path_type /= 'closed') return
        Eligible = .true.
    end function Eligible

    !> How far this gas's covariance at its own lag stands above the noise,
    !> as a multiple of its detection limit. The error code where the
    !> question cannot be asked - no limit, or no covariance to divide.
    !> ColGas and cov are the host's, allocated once, rather than a fresh
    !> automatic array on every call.
    real(kind = dbl) function SignalToNoise(g)
        integer, intent(in) :: g
        real(kind = dbl) :: floor_

        SignalToNoise = error
        !> Which noise floor. The detection limit is this engine's own -
        !> the scatter of the cross-covariance away from the peak, so it
        !> measures what the covariance itself does when there is no flux.
        !> The Lenschow noise is EddyUH's, and measures the analyser rather
        !> than the covariance; EddyUH's own comment at that line calls it a
        !> detection limit, which is how the first port of this came to
        !> divide by the wrong quantity.
        if (RPSetup%tlag_borrow_noise == 'lenschow_00') then
            floor_ = LenschowFluxNoise(Set, nrow, ncol, g)
        else
            floor_ = Essentials%detlim(g)
        end if
        if (floor_ == error) return
        if (floor_ <= 0d0) return
        !> CovarianceW takes either sign of lag, but a zero lag means the run
        !> compensates none and there is no peak to judge.
        if (RowLags(g) == 0) return
        ColGas(1:nrow) = Set(1:nrow, g)
        call CovarianceW(ColW, ColGas, nrow, RowLags(g), cov)
        if (cov == error) return
        SignalToNoise = dabs(cov) / floor_
    end function SignalToNoise

    !> Did the maximisation stop at an end of the window rather than on a peak
    !> inside it? Unreachable when covmax_debaseline is on, because the ends
    !> score zero there by construction - which is why the signal-to-noise
    !> test above is the load-bearing one of the two.
    logical function OnWindowEdge(g)
        integer, intent(in) :: g

        OnWindowEdge = .false.
        if (min_rl(g) == 0 .and. max_rl(g) == 0) return
        OnWindowEdge = RowLags(g) == min_rl(g) .or. RowLags(g) == max_rl(g)
    end function OnWindowEdge

end subroutine BorrowTimelagBelowDetectionLimit
