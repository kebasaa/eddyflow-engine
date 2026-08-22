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
! WHAT THIS THRESHOLDS AGAINST IS NOT WHAT EDDYUH USES, deliberately. EddyUH
! compares against unc3(:,3), Lenschow et al. (2000) instrumental noise
! propagated to a covariance as sqrt(sigma2_noise * sigma2_w / N). This
! compares against the Wienhold et al. (1994) detection limit, the scatter of
! the cross-covariance function away from its peak. Both are noise floors in
! covariance units and of similar size, but the second is the scatter of the
! very function whose peak is being tested, which is the question being
! asked - and it is what detlim_meth already computes, so the two settings
! compose rather than each needing their own estimator.
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
    if (RPSetup%detlim_meth /= 'wienhold_94') return

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

        !> The best-resolved tube-mate, not the first one in slot order.
        !>
        !> Taking the first put whichever gas happened to sit lowest in the
        !> metadata in charge, and on a tube where nothing clears the
        !> threshold comfortably that produced pairings that read backwards -
        !> carbon dioxide, the strongest flux on the analyser, taking its lag
        !> from nitrous oxide. Ranking by signal-to-noise says what was
        !> actually meant: borrow from whoever resolved their peak best.
        chosen = 0
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

        SignalToNoise = error
        if (Essentials%detlim(g) == error) return
        if (Essentials%detlim(g) <= 0d0) return
        !> CovarianceW takes either sign of lag, but a zero lag means the run
        !> compensates none and there is no peak to judge.
        if (RowLags(g) == 0) return
        ColGas(1:nrow) = Set(1:nrow, g)
        call CovarianceW(ColW, ColGas, nrow, RowLags(g), cov)
        if (cov == error) return
        SignalToNoise = dabs(cov) / Essentials%detlim(g)
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
