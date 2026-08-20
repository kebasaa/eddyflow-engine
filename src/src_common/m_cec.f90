!***************************************************************************
! m_cec.f90
! ---------
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
! \brief       Shared Conditional Eddy Covariance implementation following:
!              Zahn et al. (2022), Agricultural and Forest Meteorology 315, 108790.
!              https://doi.org/10.1016/j.agrformet.2021.108790
!
!              Ejections are sorted into two octants by the signs of the water
!              and carbon fluctuations they carry:
!
!                O1  w' > 0, q' > 0, c' > 0   non-stomatal - it came off the
!                                             ground, enriched in both
!                O2  w' > 0, q' > 0, c' < 0   stomatal - it passed through
!                                             stomata, which take carbon up
!
!              The conditional covariance of any scalar over those two octants
!              gives the ratio of its two components, and that ratio applied to
!              the authoritative total gives the components themselves.
!
!              The indicator functions are built from (w', q', c') and nothing
!              else, and Zahn et al. say they "remain the same" for the second
!              scalar. So the octants are a property of the CO2/water PAIR, and
!              any number of further species - carbonyl sulfide, methane - can
!              be partitioned by the same pass, at the cost of two more sums
!              each. Nothing here may make an octant depend on the scalar being
!              accumulated.
!
!              Split in three because the three parts belong in three places:
!              BuildCecPrimes and ExtractCecDescriptor run in RP, on the raw
!              high-frequency series; ApplyCecDescriptor runs wherever the final
!              corrected totals live, which is FCC when FCC follows.
! \author      Jonathan Muller
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
module m_cec
    use ieee_arithmetic
    use m_common_global_var
    implicit none
    private

    public :: ResetCecDescriptor
    public :: ResetCecFlux
    public :: BuildCecPrimes
    public :: ExtractCecDescriptor
    public :: ApplyCecDescriptor

contains

subroutine ResetCecDescriptor(descriptor)
    type(CECDescriptorType), intent(out) :: descriptor

    integer :: k

    descriptor%meth = 0
    descriptor%carbon_slot = 0
    descriptor%water_slot = 0
    descriptor%n_target = 0
    descriptor%n_valid = 0
    descriptor%n_O1 = 0
    descriptor%n_O2 = 0
    descriptor%frac_O1 = error
    descriptor%frac_O2 = error
    do k = 1, MaxNumCecTargets
        descriptor%target(k)%slot = 0
        descriptor%target(k)%f_O1 = error
        descriptor%target(k)%f_O2 = error
        descriptor%target(k)%r = error
        descriptor%target(k)%status = cec_rejected
        descriptor%target(k)%valid = .false.
    end do
end subroutine ResetCecDescriptor

subroutine ResetCecFlux(flux)
    type(CECFluxType), intent(out) :: flux

    integer :: k

    do k = 1, MaxNumCecTargets
        flux%comp(k)%nonstomatal = error
        flux%comp(k)%stomatal = error
        flux%comp(k)%total = error
        flux%comp(k)%status = cec_rejected
    end do
    flux%E_cec_ET = error
    flux%Tr_cec_ET = error
end subroutine ResetCecFlux

!***************************************************************************
!
! \brief       Screened, gap-filled, detrended fluctuations for one pairing.
! \author      Jonathan Muller
! \note        Runs on the RAW series, before the run's own detrending, and in
!              the order Zahn et al. Sec. 3.2 states: fill the short gaps, then
!              delete what the analyser diagnostics condemn, then detrend.
!
!              This used to work on the finished E2Primes instead, which got
!              the order wrong twice over. The screen ran before the gap filler
!              rather than after, so every rejected run of four samples or fewer
!              was linearly reconstructed and classified into an octant anyway -
!              the screen only ever removed runs of five and longer. And
!              detrending had already happened, so the trend that produced the
!              fluctuations had been fitted through the very samples being
!              rejected.
!
!              The working set is compact - w and the pairing's targets, six or
!              seven columns - so redoing the detrending costs almost nothing
!              and, crucially, touches neither E2Set nor E2Primes. The fluxes
!              themselves are computed from those and must not move because a
!              partition asked for a stricter screen.
!***************************************************************************
subroutine BuildCecPrimes(pair, rawset, nrow, ncol, userset, nuser, tconst, &
    setup, primes, ok)
    type(CECResolvedPairType), intent(in) :: pair
    integer, intent(in) :: nrow
    integer, intent(in) :: ncol
    real(kind = dbl), intent(in) :: rawset(nrow, ncol)
    integer, intent(in) :: nuser
    !> Assumed shape, not (nrow, max(nuser,1)): a site with no custom columns
    !> passes a zero-column array, and an explicit-shape dummy claiming one
    !> column would be describing storage that is not there.
    real(kind = dbl), intent(in) :: userset(:, :)
    integer, intent(in) :: tconst
    type(CECSetupType), intent(in) :: setup
    real(kind = dbl), intent(out) :: primes(nrow, MaxNumCecTargets + 1)
    logical, intent(out) :: ok

    integer :: k
    integer :: n
    integer :: ntarget
    integer :: slot
    integer :: sig_col
    integer :: ncompact
    integer :: slots(MaxNumCecTargets)
    logical :: sig_is_rssi
    real(kind = dbl), allocatable :: work(:, :)
    type(ColType), allocatable :: lcol(:)
    type(StatsType), allocatable :: lstats
    integer, external :: CecSignalColumnFor
    logical, external :: CecSignalIsRssi
    external :: CecTargetSlots

    primes = error
    ok = .false.
    if (nrow < 2) return

    call CecTargetSlots(pair, slots, ntarget)
    do k = 1, ntarget
        slot = slots(k)
        if (slot < firstGas .or. slot > lastGas) return
        if (slot > ncol) return
        !> A column the project declares is not necessarily a column this
        !> period has: EliminateCorruptedVariables clears the flag for a period
        !> whose values went out of range, and the detrending routines then
        !> leave that column of Primes untouched. Reading it would be reading
        !> whatever the freshly allocated array happened to contain, which
        !> ieee_is_finite will often accept.
        if (.not. E2Col(slot)%present) return
    end do

    ncompact = ntarget + 1
    allocate(work(nrow, ncompact))
    allocate(lcol(ncompact))
    allocate(lstats)

    work(:, 1) = rawset(:, w)
    do k = 1, ntarget
        work(:, k + 1) = rawset(:, slots(k))
    end do

    do k = 1, ncompact
        call InterpolateShortCecGaps(work(:, k), setup%max_gap_fill)
    end do

    if (setup%signal_strength > 0d0 .and. nuser > 0) then
        do k = 1, ntarget
            sig_col = CecSignalColumnFor(slots(k))
            if (sig_col <= 0 .or. sig_col > nuser) cycle
            if (sig_col > size(userset, 2)) cycle
            sig_is_rssi = CecSignalIsRssi(slots(k), sig_col)
            call FilterCecSignalStrength(work(:, k + 1), userset(:, sig_col), &
                setup%signal_strength, sig_is_rssi)
        end do
    end if

    !> Forced present so Fluctuations writes the column whatever the period
    !> held. A detrending routine leaves an absent column of Primes untouched,
    !> and untouched here means whatever the freshly allocated array contained;
    !> a column of error values is a period with nothing in it, which the
    !> caller can see, and uninitialised memory is not.
    lcol(1) = E2Col(w)
    lcol(1)%present = .true.
    do k = 1, ntarget
        lcol(k + 1) = E2Col(slots(k))
    end do

    !> Only the block-average branch reads these, and it must read the mean of
    !> the screened series - not the one BasicStats took before the screening.
    lstats%Mean = error
    do k = 1, ncompact
        lstats%Mean(k) = CecMean(work(:, k))
    end do

    call Fluctuations(work, primes(:, 1:ncompact), nrow, ncompact, tconst, &
        lstats, lcol)

    n = 0
    do k = 1, nrow
        if (CecValueIsValid(primes(k, 1))) n = n + 1
    end do
    ok = n >= 2

    deallocate(work, lcol, lstats)
end subroutine BuildCecPrimes

!***************************************************************************
!
! \brief       Octant statistics and component ratios for one pairing.
! \author      Jonathan Muller
! \note        `primes` is what BuildCecPrimes produced: column one is w', and
!              column 1+k is target k, water first and carbon second.
!***************************************************************************
subroutine ExtractCecDescriptor(pair, primes, nrow, stationarity_carbon, &
    stationarity_water, descriptor, setup)
    type(CECResolvedPairType), intent(in) :: pair
    integer, intent(in) :: nrow
    real(kind = dbl), intent(in) :: primes(nrow, MaxNumCecTargets + 1)
    integer, intent(in) :: stationarity_carbon
    integer, intent(in) :: stationarity_water
    type(CECDescriptorType), intent(out) :: descriptor
    type(CECSetupType), intent(in) :: setup

    integer :: i
    integer :: k
    integer :: ntarget
    integer :: slots(MaxNumCecTargets)
    integer :: iw
    integer :: ic
    real(kind = dbl) :: sum_O1(MaxNumCecTargets)
    real(kind = dbl) :: sum_O2(MaxNumCecTargets)
    real(kind = dbl) :: sigma_w
    real(kind = dbl) :: sigma_q
    real(kind = dbl) :: sigma_c
    logical :: usable
    external :: CecTargetSlots

    call ResetCecDescriptor(descriptor)
    if (nrow < 2) return

    call CecTargetSlots(pair, slots, ntarget)
    descriptor%meth = pair%meth
    descriptor%carbon_slot = pair%carbon_slot
    descriptor%water_slot = pair%water_slot
    descriptor%n_target = ntarget
    do k = 1, ntarget
        descriptor%target(k)%slot = slots(k)
    end do

    !> Columns of `primes`, not gas slots: water and carbon are the first two
    !> targets by construction, so they sit right behind w'.
    iw = 1 + cecTargetWater
    ic = 1 + cecTargetCarbon

    call CecStandardDeviations(primes(:, 1), primes(:, ic), primes(:, iw), &
        sigma_w, sigma_c, sigma_q)

    sum_O1 = 0d0
    sum_O2 = 0d0
    do i = 1, nrow
        usable = CecValueIsValid(primes(i, 1))
        do k = 1, ntarget
            if (.not. CecValueIsValid(primes(i, k + 1))) usable = .false.
        end do
        if (.not. usable) cycle

        !> A point the hole rejects still counts here. The paper normalises the
        !> sample fluxes by every sample, not by the ones that survived, and
        !> the ratio is invariant to the choice anyway - it is the octant
        !> fractions the hole is allowed to thin.
        descriptor%n_valid = descriptor%n_valid + 1
        if (primes(i, 1) > 0d0 .and. primes(i, iw) > 0d0 &
            .and. primes(i, ic) > 0d0) then
            if (.not. CecPassesHyperbolicThreshold(primes(i, 1), primes(i, iw), &
                primes(i, ic), setup%h, sigma_w, sigma_q, sigma_c)) cycle
            descriptor%n_O1 = descriptor%n_O1 + 1
            do k = 1, ntarget
                sum_O1(k) = sum_O1(k) + primes(i, 1) * primes(i, k + 1)
            end do
        else if (primes(i, 1) > 0d0 .and. primes(i, iw) > 0d0 &
            .and. primes(i, ic) < 0d0) then
            if (.not. CecPassesHyperbolicThreshold(primes(i, 1), primes(i, iw), &
                primes(i, ic), setup%h, sigma_w, sigma_q, sigma_c)) cycle
            descriptor%n_O2 = descriptor%n_O2 + 1
            do k = 1, ntarget
                sum_O2(k) = sum_O2(k) + primes(i, 1) * primes(i, k + 1)
            end do
        end if
    end do

    !> Zahn et al. retained periods with at least 90% instantaneous data.
    if (dble(descriptor%n_valid) < setup%min_valid * dble(nrow)) return
    if (setup%max_stationarity > 0d0) then
        if (stationarity_carbon == ierror .or. stationarity_water == ierror) return
        if (dble(stationarity_carbon) > setup%max_stationarity &
            .or. dble(stationarity_water) > setup%max_stationarity) return
    end if

    descriptor%frac_O1 = dble(descriptor%n_O1) / dble(descriptor%n_valid)
    descriptor%frac_O2 = dble(descriptor%n_O2) / dble(descriptor%n_valid)
    if (descriptor%frac_O1 + descriptor%frac_O2 < setup%min_o1_o2) return

    do k = 1, ntarget
        descriptor%target(k)%f_O1 = sum_O1(k) / dble(descriptor%n_valid)
        descriptor%target(k)%f_O2 = sum_O2(k) / dble(descriptor%n_valid)
        call ClassifyCecTarget(descriptor%target(k), descriptor%frac_O1, &
            descriptor%frac_O2, setup)
    end do
end subroutine ExtractCecDescriptor

!***************************************************************************
!
! \brief       Which of the five verdicts one target's sample fluxes support.
! \author      Jonathan Muller
! \note        The occupancy fallbacks are the pairing's, not the target's:
!              they ask how many points the octants hold, which is the same
!              question whichever scalar is being summed over them.
!***************************************************************************
subroutine ClassifyCecTarget(target, frac_O1, frac_O2, setup)
    type(CECTargetType), intent(inout) :: target
    real(kind = dbl), intent(in) :: frac_O1
    real(kind = dbl), intent(in) :: frac_O2
    type(CECSetupType), intent(in) :: setup

    target%r = error
    target%status = cec_rejected
    target%valid = .false.

    !> Too few points in one octant and its component is taken as negligible,
    !> the whole flux going to the other. Zahn et al. Sec. 2.4.
    if (frac_O1 < setup%min_octant) then
        target%status = cec_all_stomatal
        target%valid = .true.
    else if (frac_O2 < setup%min_octant) then
        target%status = cec_all_nonstomatal
        target%valid = .true.
    else if (target%f_O2 == 0d0) then
        !> Nothing in the stomatal octant to divide by: the whole flux is the
        !> other component, which is what a ratio of infinity would have said.
        target%status = cec_all_nonstomatal
        target%valid = .true.
    else if (target%f_O1 == 0d0) then
        !> And the mirror of it. Reached through the ratio this would be
        !> total/(1 + 1/0), which is the right answer arrived at through a
        !> division by zero.
        target%status = cec_all_stomatal
        target%valid = .true.
    else
        target%r = target%f_O1 / target%f_O2
        !> Exactly -1 is a division by zero however wide the band is, so it is
        !> rejected even when the guard below is switched off. Unreachable in
        !> practice - the two sums would have to cancel to the last bit - but
        !> an infinity in an output column is not a thing to leave to luck.
        if (target%r == -1d0) then
            target%status = cec_singular
            return
        end if
        !> Only opposite-signed components can cancel, and only their ratio can
        !> approach -1. Two sinks, or two sources, put r on the positive side
        !> where there is no singularity to fall into - which is why the water
        !> partition never suffers this, and why a soil that takes up carbonyl
        !> sulfide alongside the canopy does not either.
        if (target%r < 0d0 .and. CecIsSingular(target%r, setup%singular_band)) then
            target%status = cec_singular
        else
            target%status = cec_normal
            target%valid = .true.
        end if
    end if
end subroutine ClassifyCecTarget

!***************************************************************************
!
! \brief       Turn the ratios into fluxes, against the authoritative totals.
! \author      Jonathan Muller
! \note        `totals` is indexed like the descriptor's targets. It comes from
!              Flux3, so it carries every correction the ratios do not - the
!              density and frequency-response terms the paper says must be in
!              the total and cannot be in an instantaneous fluctuation.
!***************************************************************************
subroutine ApplyCecDescriptor(descriptor, totals, flux)
    type(CECDescriptorType), intent(in) :: descriptor
    real(kind = dbl), intent(in) :: totals(MaxNumCecTargets)
    type(CECFluxType), intent(out) :: flux

    integer :: k
    real(kind = dbl) :: total

    call ResetCecFlux(flux)
    if (descriptor%meth == 0) return

    do k = 1, descriptor%n_target
        if (.not. CecTargetIsWanted(k, descriptor%meth)) cycle

        total = totals(k)
        flux%comp(k)%total = total
        flux%comp(k)%status = descriptor%target(k)%status
        if (.not. descriptor%target(k)%valid) cycle
        if (total == error) cycle

        if (CecTotalContradictsOctants(descriptor%target(k), total)) then
            flux%comp(k)%status = cec_wrong_sign
            cycle
        end if

        select case (descriptor%target(k)%status)
            case (cec_normal)
                flux%comp(k)%nonstomatal = total / (1d0 + 1d0 / descriptor%target(k)%r)
                flux%comp(k)%stomatal = total / (1d0 + descriptor%target(k)%r)
            case (cec_all_stomatal)
                flux%comp(k)%nonstomatal = 0d0
                flux%comp(k)%stomatal = total
            case (cec_all_nonstomatal)
                flux%comp(k)%nonstomatal = total
                flux%comp(k)%stomatal = 0d0
        end select
    end do

    if (flux%comp(cecTargetWater)%nonstomatal /= error &
        .and. flux%comp(cecTargetWater)%stomatal /= error) then
        flux%E_cec_ET = flux%comp(cecTargetWater)%nonstomatal * h2o_to_ET
        flux%Tr_cec_ET = flux%comp(cecTargetWater)%stomatal * h2o_to_ET
    end if
end subroutine ApplyCecDescriptor

!> Water on 1 and 2, carbon on 1 and 3, extras whenever the pairing runs at
!> all: an extra species needs the octants, not the choice of which of the two
!> gases that define them is itself reported.
logical function CecTargetIsWanted(k, meth)
    integer, intent(in) :: k
    integer, intent(in) :: meth

    select case (k)
        case (cecTargetWater)
            CecTargetIsWanted = meth == 1 .or. meth == 2
        case (cecTargetCarbon)
            CecTargetIsWanted = meth == 1 .or. meth == 3
        case default
            CecTargetIsWanted = .true.
    end select
end function CecTargetIsWanted

!***************************************************************************
!
! \brief       Does the authoritative total disagree with the octant sums?
! \author      Jonathan Muller
! \note        Substituting Eq. 10 into Eq. 12 and cancelling gives the two
!              components as
!
!                nonstomatal = total * f_O1 / (f_O1 + f_O2)
!                stomatal    = total * f_O2 / (f_O1 + f_O2)
!
!              so each carries the sign of its own sample flux only when the
!              total points the same way as the two of them together. When it
!              does not, the arithmetic still returns numbers that sum to the
!              total - and they are a negative respiration beside a positive
!              photosynthesis, or a negative transpiration.
!
!              One condition covers every species and both regimes. Water is
!              the case where both sample fluxes are positive by construction,
!              w' and q' being positive in both octants, so a downward total -
!              dewfall, condensation on the canopy - is caught. Carbon is the
!              case where they have opposite signs and can cancel, so it is
!              caught only when they cancel the wrong way, which is a real
!              observation about a period rather than a property of the gas.
!              Carbonyl sulfide is whichever of the two its soil makes it: a
!              soil source against canopy uptake behaves like carbon, a soil
!              that takes it up alongside the canopy behaves like water.
!
!              The ratio is still reported. What it is telling you when this
!              fires is that the octants and the corrected total describe
!              different half-hours, and that is worth reading.
!***************************************************************************
logical function CecTotalContradictsOctants(target, total)
    type(CECTargetType), intent(in) :: target
    real(kind = dbl), intent(in) :: total

    real(kind = dbl) :: sampled

    CecTotalContradictsOctants = .false.
    if (total == 0d0) return
    if (.not. CecValueIsValid(target%f_O1)) return
    if (.not. CecValueIsValid(target%f_O2)) return

    sampled = target%f_O1 + target%f_O2
    if (sampled == 0d0) return

    CecTotalContradictsOctants = sampled * total < 0d0
end function CecTotalContradictsOctants

!***************************************************************************
!
! \brief       Delete the samples an analyser's own diagnostic condemns.
! \author      Jonathan Muller
! \note        The threshold is one number read in two directions. RSSI is a
!              received-signal strength, so low is dirty; the AGC of a
!              pre-5.3.0 LI-7500 is the gain the instrument had to apply to see
!              through its windows, so HIGH is dirty. Both are conventionally
!              compared against 70, and this used to test "below threshold"
!              whichever it was - so on an old open-path analyser it kept the
!              dirtiest samples and threw away the cleanest.
!***************************************************************************
subroutine FilterCecSignalStrength(values, signal_strength, threshold, is_rssi)
    real(kind = dbl), intent(inout) :: values(:)
    real(kind = dbl), intent(in) :: signal_strength(:)
    real(kind = dbl), intent(in) :: threshold
    logical, intent(in) :: is_rssi

    integer :: i
    integer :: n

    n = min(size(values), size(signal_strength))
    do i = 1, n
        if (.not. CecValueIsValid(signal_strength(i))) cycle
        if (is_rssi) then
            if (signal_strength(i) < threshold) values(i) = error
        else
            if (signal_strength(i) > threshold) values(i) = error
        end if
    end do
end subroutine FilterCecSignalStrength

!> Standard deviations of the three conditioning series over the samples all
!> three of them have. They normalise the hyperbolic hole below, so they must
!> come from the same screened, gap-filled arrays the octants are built from -
!> not from Stats%StDev, which is computed on the unscreened E2Primes.
subroutine CecStandardDeviations(w_prime, c_prime, q_prime, sigma_w, sigma_c, sigma_q)
    real(kind = dbl), intent(in) :: w_prime(:)
    real(kind = dbl), intent(in) :: c_prime(:)
    real(kind = dbl), intent(in) :: q_prime(:)
    real(kind = dbl), intent(out) :: sigma_w
    real(kind = dbl), intent(out) :: sigma_c
    real(kind = dbl), intent(out) :: sigma_q

    integer :: i
    integer :: n
    real(kind = dbl) :: sw, sc, sq
    real(kind = dbl) :: sw2, sc2, sq2

    sigma_w = 0d0
    sigma_c = 0d0
    sigma_q = 0d0
    n = 0
    sw = 0d0
    sc = 0d0
    sq = 0d0
    sw2 = 0d0
    sc2 = 0d0
    sq2 = 0d0

    do i = 1, min(size(w_prime), min(size(c_prime), size(q_prime)))
        if (.not. CecValueIsValid(w_prime(i))) cycle
        if (.not. CecValueIsValid(c_prime(i))) cycle
        if (.not. CecValueIsValid(q_prime(i))) cycle
        n = n + 1
        sw = sw + w_prime(i)
        sw2 = sw2 + w_prime(i)**2
        sc = sc + c_prime(i)
        sc2 = sc2 + c_prime(i)**2
        sq = sq + q_prime(i)
        sq2 = sq2 + q_prime(i)**2
    end do
    if (n < 2) return

    sigma_w = dsqrt(max(0d0, sw2 / dble(n) - (sw / dble(n))**2))
    sigma_c = dsqrt(max(0d0, sc2 / dble(n) - (sc / dble(n))**2))
    sigma_q = dsqrt(max(0d0, sq2 / dble(n) - (sq / dble(n))**2))
end subroutine CecStandardDeviations

real(kind = dbl) function CecMean(values)
    real(kind = dbl), intent(in) :: values(:)

    integer :: i
    integer :: n
    real(kind = dbl) :: total

    CecMean = error
    n = 0
    total = 0d0
    do i = 1, size(values)
        if (.not. CecValueIsValid(values(i))) cycle
        n = n + 1
        total = total + values(i)
    end do
    if (n > 0) CecMean = total / dble(n)
end function CecMean

!***************************************************************************
!
! \brief       The hyperbolic hole of Thomas et al. (2008), Eqs. 7a/7b.
! \author      Jonathan Muller
! \note        Leave out the events nearest the origin, whose sign is noise
!              rather than a surface signature. That is what makes the method
!              usable when the turbulence is weak: without it, instrument noise
!              around c' = 0 populates both octants about equally, and the
!              carbon ratio walks straight into the r ~ -1 singularity.
!
!              H is dimensionless and scales the instantaneous flux against
!              sigma_w*sigma_s, which is what lets one value mean the same thing
!              at every site. It used to be compared against the raw product, so
!              what it meant depended on whether the gas was carried in umol/mol
!              or mmol/m3, and no single number was usable anywhere.
!
!              The hole belongs to the PAIR, not to the scalar being
!              accumulated. Testing each scalar against its own sigma would make
!              the octants depend on which scalar was being summed, and the
!              property that lets an arbitrary species be partitioned by these
!              octants would be gone.
!***************************************************************************
logical function CecPassesHyperbolicThreshold(w_prime, q_prime, c_prime, h, &
    sigma_w, sigma_q, sigma_c)
    real(kind = dbl), intent(in) :: w_prime
    real(kind = dbl), intent(in) :: q_prime
    real(kind = dbl), intent(in) :: c_prime
    real(kind = dbl), intent(in) :: h
    real(kind = dbl), intent(in) :: sigma_w
    real(kind = dbl), intent(in) :: sigma_q
    real(kind = dbl), intent(in) :: sigma_c

    if (h <= 0d0) then
        CecPassesHyperbolicThreshold = .true.
    else
        CecPassesHyperbolicThreshold = &
            abs(w_prime * q_prime) >= h * sigma_w * sigma_q .and. &
            abs(w_prime * c_prime) >= h * sigma_w * sigma_c
    end if
end function CecPassesHyperbolicThreshold

!> Two components nearly cancelling puts 1 + r on top of zero, so the partition
!> blows up on a total that is itself near zero. Zahn et al. reject the band
!> -1.2 < r < -0.8 and say the width is dataset-dependent, hence the setting.
logical function CecIsSingular(r, band)
    real(kind = dbl), intent(in) :: r
    real(kind = dbl), intent(in) :: band

    if (band <= 0d0) then
        CecIsSingular = .false.
    else
        CecIsSingular = abs(r + 1d0) < band
    end if
end function CecIsSingular

subroutine InterpolateShortCecGaps(values, max_gap)
    real(kind = dbl), intent(inout) :: values(:)
    integer, intent(in) :: max_gap

    integer :: first_gap
    integer :: gap_length
    integer :: i
    integer :: k
    real(kind = dbl) :: increment

    i = 2
    do while (i < size(values))
        if (CecValueIsValid(values(i))) then
            i = i + 1
            cycle
        end if

        first_gap = i
        do while (i <= size(values))
            if (CecValueIsValid(values(i))) exit
            i = i + 1
        end do
        gap_length = i - first_gap

        if (i <= size(values)) then
            if (gap_length <= max_gap .and. CecValueIsValid(values(first_gap - 1)) &
                .and. CecValueIsValid(values(i))) then
                increment = (values(i) - values(first_gap - 1)) / dble(gap_length + 1)
                do k = 1, gap_length
                    values(first_gap + k - 1) = values(first_gap - 1) + increment * dble(k)
                end do
            end if
        end if
    end do
end subroutine InterpolateShortCecGaps

logical function CecValueIsValid(value)
    real(kind = dbl), intent(in) :: value

    CecValueIsValid = ieee_is_finite(value) .and. value /= error
end function CecValueIsValid

end module m_cec
