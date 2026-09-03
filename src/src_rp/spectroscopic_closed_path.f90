!***************************************************************************
! spectroscopic_closed_path.f90
! -----------------------------
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
! \brief       Remove the water-vapour spectroscopic bias from a closed-path
!              laser analyser, sample by sample.
! \author      Jonathan Muller
!
! \note
! Water vapour broadens the absorption lines a laser analyser measures, so the
! mixing ratio it reports depends on humidity beyond simple dilution. Peltola
! et al. (2014) describe the effect as a quadratic scaling of the analyser's
! sensitivity,
!
!     x_reported = x_true * (1 + a*chi_q + b*chi_q^2)
!
! with a and b per gas and per analyser, so the correction is a division by
! that polynomial. Applied here point by point, after Chen et al. (2010),
! using the water the SAME analyser read at the SAME sample: the broadening
! happens in the cell at that instant, so no time-lag alignment enters - the
! two columns are already the readings of one moment.
!
! Doing it on the series rather than on the flux earns three things. The mean
! and the covariance are both corrected, by construction and consistently.
! The full nonlinearity is kept, where a flux-level form would carry only the
! first-order term. And it reaches gases reported as a mixing ratio, which
! Level 2 returns early on - at a site whose laser reports mixing ratios,
! which is the common case for carbonyl sulfide and nitrous oxide, a
! flux-level correction placed with the density terms would do nothing at all.
!
! THE COEFFICIENTS HERE ARE SPECTROSCOPIC ONLY, and are not EddyUH's.
!
! EddyUH divides by the same polynomial - Functions_Library/dilucorr.m, with
! the water in mol/mol, identically to this routine - but folds the dilution
! into it: its a = -1, b = 0 is the pure-dilution case, giving a divisor of
! (1 - chi_q). EddyFlow corrects the density separately and more fully, so
! here the identity is a = b = 0 and the coefficient carries only the
! spectroscopic excess.
!
! Mapping one to the other, D_here = D_EddyUH / (1 - chi_q), which to second
! order in chi_q is
!
!     a_here = a_EddyUH + 1
!     b_here = b_EddyUH + a_EddyUH + 1
!
! and the pure-dilution pair maps to (0, 0) exactly, as it must. A published
! value entered unconverted would count the dilution a second time: EddyUH's
! table gives -1.39 for nitrous oxide on an Aerodyne cw-QCL, which is -0.39
! here.
!
! \sa          EC_Software_FluxCalc/EddyUH_spect_CP.m and
!              Functions_Library/dilucorr.m in EddyUH 1.7b, which are the
!              flux-level and point-by-point forms of the same correction.
!              EddyUH's "crosstalk" coefficients are these coefficients: the
!              name is historical. crosstalk_coeff_CP.m is a table of
!              published (a, b) pairs that its setup dialogs use to seed
!              set_Gan.crosstalkcoeff, which EddyUH_spect_CP.m then reads as
!              its a and b. EddyUH_ctc_CP.m, which looks like a separate
!              additive crosstalk correction, is unreachable - nothing calls
!              it, and it would fail if anything did, since it multiplies a
!              scalar water value by a two-element coefficient vector.
!              EddyUH's flux-level form folds the dilution in with the
!              spectroscopy, because its own WPL is applied first and the
!              formula partly undoes it; substituting one into the other,
!              the dilution cancels and what is left is the division above.
!              EddyFlow's WPL is a different and fuller formulation, with
!              cell temperature and pressure terms EddyUH's closed-path WPL
!              does not carry, so only the spectroscopic part is taken here
!              and the density work is left where it already is.
!***************************************************************************
subroutine SpectroscopicClosedPath(Set, nrow, ncol, printout)
    !> RP, not common: the correction acts on the raw series, which only RP
    !> holds, and FCC never redoes it. RPSetup lives here too.
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow, ncol
    logical, intent(in) :: printout
    real(kind = dbl), intent(inout) :: Set(nrow, ncol)
    !> local variables
    logical, external :: GasSlotIsWater
    real(kind = dbl), allocatable :: ChiQ(:, :)
    integer :: watSlot(MaxNumGases)
    integer :: nwat
    integer :: gas, msl, iw, iwat
    integer :: ncorrected
    real(kind = dbl) :: coef_a, coef_b
    logical :: anyCoefficients

    if (RPSetup%spectro_meth /= 'chen_10') return

    !> Nothing to do unless some column actually declares a coefficient. A
    !> project may switch the method on and leave the metadata at zero, which
    !> is the identity - saying so here keeps the pass off the data entirely.
    anyCoefficients = .false.
    do gas = firstGas, lastGas
        if (.not. E2Col(gas)%present) cycle
        if (E2Col(gas)%instr%path_type /= 'closed') cycle
        if (E2Col(gas)%spectro_a /= 0d0 .or. E2Col(gas)%spectro_b /= 0d0) &
            anyCoefficients = .true.
    end do
    if (.not. anyCoefficients) return

    !> The hygrometers whose readings can serve as chi_q. Closed path only:
    !> the broadening is a property of the cell the gas passed through, and an
    !> open-path hygrometer metres away is measuring different air.
    nwat = 0
    do gas = firstGas, lastGas
        if (.not. GasSlotIsWater(gas)) cycle
        if (.not. E2Col(gas)%present) cycle
        if (E2Col(gas)%instr%path_type /= 'closed') cycle
        nwat = nwat + 1
        watSlot(nwat) = gas
    end do
    if (nwat == 0) return

    if (printout) write(*, '(a)', advance = 'no') &
        '  Removing the spectroscopic effect of water vapour..'
    if (printout) write(ulog, '(a)', advance = 'no') &
        '  Removing the spectroscopic effect of water vapour..'

    !> Pass one: every hygrometer's mole fraction in mol/mol, taken from the
    !> series as measured. A separate pass because a hygrometer may itself be
    !> corrected below, and every gas has to be divided by the water that was
    !> read, not by a water that has already moved.
    allocate(ChiQ(nrow, nwat))
    do iw = 1, nwat
        call WaterMoleFraction(Set, nrow, ncol, watSlot(iw), ChiQ(:, iw))
    end do

    !> Pass two: divide each declared column by its own analyser's water.
    ncorrected = 0
    do gas = firstGas, lastGas
        if (.not. E2Col(gas)%present) cycle
        if (E2Col(gas)%instr%path_type /= 'closed') cycle
        coef_a = E2Col(gas)%spectro_a
        coef_b = E2Col(gas)%spectro_b
        if (coef_a == 0d0 .and. coef_b == 0d0) cycle

        !> Water against itself is self-broadening. Real, and not part of what
        !> Peltola et al. publish - they derive the effect of water on another
        !> gas's lines. EddyUH corrects it with a flux-level coefficient whose
        !> derivation its own comment says is unpublished. Offered here in the
        !> same point-by-point form as everything else, and only when asked
        !> for, so that a run cannot acquire it by accident.
        if (GasSlotIsWater(gas) .and. .not. RPSetup%spectro_water) cycle

        !> Which hygrometer corrects this gas. moist_ref already prefers one
        !> on the gas's own analyser, which is the only kind that can be
        !> right here; a reference on another instrument read different air.
        msl = E2Col(gas)%moist_ref
        if (GasSlotIsWater(gas)) msl = gas
        iwat = 0
        do iw = 1, nwat
            if (watSlot(iw) == msl) iwat = iw
        end do
        if (iwat == 0) cycle
        if (E2Col(gas)%instr%model /= E2Col(msl)%instr%model) cycle

        call DivideBySensitivity(Set(:, gas), ChiQ(:, iwat), nrow, coef_a, coef_b)
        ncorrected = ncorrected + 1
    end do

    deallocate(ChiQ)
    if (printout) write(*, '(a, i0, a)') ' Done (', ncorrected, ' column(s)).'
    if (printout) write(ulog, '(a, i0, a)') ' Done (', ncorrected, ' column(s)).'

contains

    !***********************************************************************
    !> This hygrometer's reading as a mole fraction in mol/mol.
    !>
    !> Shapes follow PointByPointToMixingRatio, which does the same
    !> conversion a stage later and for a different purpose. Molar density is
    !> declined rather than converted: it would need the cell molar volume,
    !> and a hygrometer reporting density on an analyser whose cell block is
    !> absent would silently produce a chi_q of nothing. A column the
    !> conversion cannot reach leaves error, and every sample it touches is
    !> left uncorrected.
    !***********************************************************************
    subroutine WaterMoleFraction(LocSet, lnrow, lncol, slot, chi)
        implicit none
        integer, intent(in) :: lnrow, lncol, slot
        real(kind = dbl), intent(in) :: LocSet(lnrow, lncol)
        real(kind = dbl), intent(out) :: chi(lnrow)

        chi(:) = error
        select case (E2Col(slot)%measure_type)
            case ('mole_fraction')
                !> mmol/mol on this grid, as everywhere else in the raw set.
                where (LocSet(:, slot) /= error) chi(:) = LocSet(:, slot) * 1d-3
            case ('mixing_ratio')
                !> Dry-air mixing ratio to mole fraction, then to mol/mol.
                where (LocSet(:, slot) /= error) &
                    chi(:) = (LocSet(:, slot) / (1d0 + LocSet(:, slot) / 1d3)) * 1d-3
        end select
    end subroutine WaterMoleFraction

    !***********************************************************************
    !> Divide a column by 1 + a*chi + b*chi^2, sample by sample.
    !>
    !> A sample whose water is missing is left alone rather than voided. The
    !> correction is a fraction of a per cent at ordinary humidity, so
    !> carrying a handful of uncorrected samples costs less than discarding
    !> them - and discarding them would put gaps into a series the despiking
    !> and the spectra have already been told is continuous.
    !>
    !> A sensitivity outside a sane band is refused rather than applied.
    !> Physically it sits within a few per cent of one; anything below
    !> min_sensitivity would scale the reading by more than tenfold, which no
    !> instrument characterisation means and which a single bad water sample
    !> can produce. Guarding only against zero would still admit a divisor of
    !> 1d-15 and turn a sub-percent bias into an infinity that then travels
    !> into the statistics and the spectra.
    !***********************************************************************
    subroutine DivideBySensitivity(col, chi, lnrow, a_coef, b_coef)
        implicit none
        integer, intent(in) :: lnrow
        real(kind = dbl), intent(inout) :: col(lnrow)
        real(kind = dbl), intent(in) :: chi(lnrow)
        real(kind = dbl), intent(in) :: a_coef, b_coef
        !> Heap, not an automatic array: at 20 Hz over half an hour this is
        !> some 288 kB, and it is called once per corrected column from deep
        !> inside the period loop.
        real(kind = dbl), allocatable :: d(:)
        real(kind = dbl), parameter :: min_sensitivity = 0.1d0

        allocate(d(lnrow))
        d(:) = error
        where (chi(:) /= error) &
            d(:) = 1d0 + a_coef * chi(:) + b_coef * chi(:) * chi(:)
        where (col(:) /= error .and. chi(:) /= error &
               .and. d(:) >= min_sensitivity) &
            col(:) = col(:) / d(:)
        deallocate(d)
    end subroutine DivideBySensitivity

end subroutine SpectroscopicClosedPath
