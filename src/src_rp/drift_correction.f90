!***************************************************************************
! drift_correction.f90
! --------------------
! Copyright © 2011-2026, LI-COR Biosciences, Gerardo Fratini
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
! \brief       Corrects concentration biases due to instrumental drifts and
!              based on calibration-check data
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine DriftCorrection(Set, nrow, ncol, locCol, ncol2, nCalibEvents, InitialTimestamp)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow
    integer, intent(in) :: ncol
    integer, intent(in) :: ncol2
    integer, intent(in) :: nCalibEvents
    type(ColType), intent(in) :: locCol(ncol2)
    type(DateType), intent(in) :: InitialTimestamp
    real(kind = dbl), intent(inout) :: Set(nrow, ncol)
    !> Local variables
    integer :: i
    integer :: gas
    integer :: wsl
    integer :: msl
    integer :: nPeriods
    integer, external :: NumOfPeriods
    real(kind = dbl) :: lDrift(GHGNumVar)
    real(kind = dbl) :: MeanAbs(1)
    real(kind = dbl) :: TempFact
    real(kind = dbl) :: abs_scale
    real(kind = dbl) :: broadening
    !> Mole fraction of the water this gas names, from whichever source.
    real(kind = dbl) :: chi_moist
    character(32) :: label
    include '../src_common/interfaces.inc'


    !> Calculate best guesses of cell and air temperature and pressure
    call AirAndCellParameters()

    !> Calculate drift for the current period, depending on the chosen method
    !>
    !> Over the whole gas block rather than co2:h2o. The arrays behind these -
    !> Calib%offset, %ri, %rf and refCounts - have always been GHGNumVar wide;
    !> only the slices were two long.
    lDrift(firstGas:lastGas) = error
    select case (trim(adjustl(DriftCorr%method)))
        case ('linear')
            do i = 1, nCalibEvents
                if (Calib(i-1)%ts <= InitialTimestamp &
                    .and. InitialTimestamp <= Calib(i)%ts) then
                    nPeriods = &
                        NumOfPeriods(Calib(i-1)%ts, InitialTimestamp, DateStep)
                    lDrift(firstGas:lastGas) = &
                        Calib(i-1)%offset(firstGas:lastGas) &
                        + (Calib(i)%offset(firstGas:lastGas) &
                           - Calib(i-1)%offset(firstGas:lastGas)) &
                        / Calib(i)%numPeriods * nPeriods
                    exit
                end if
            end do

        case ('signal_strength')
            !> In case of signal strength proxy, calculate
            !> drift based on signal strength
            !> Temperature dependency (LI-7200 manual REv5, Eq. 3-32)

            if (Ambient%Tcell > 0d0) then
                TempFact = 0.6d0 + 0.4d0 / (1d0 &
                    + DriftCorr%b * dexp(DriftCorr%c * (Ambient%Tcell - 273.15d0)))
            else
                TempFact = 1d0
            end if

            !> Detect relevant drift data (all data for periods between
            !> t1 and t2 are stored in Calib(t1))
            do i = 1, nCalibEvents
                if (Calib(i-1)%ts <= InitialTimestamp &
                    .and. InitialTimestamp <= Calib(i)%ts) then
                    where (refCounts(firstGas:lastGas) /= error &
                        .and. Calib(i)%ri(firstGas:lastGas) /= error &
                        .and. Calib(i)%rf(firstGas:lastGas) /= error &
                        .and. Calib(i)%offset(firstGas:lastGas) /= error)
                        lDrift(firstGas:lastGas) = &
                            (refCounts(firstGas:lastGas) &
                             - Calib(i)%ri(firstGas:lastGas) * TempFact) &
                            / (Calib(i)%rf(firstGas:lastGas) &
                               - Calib(i)%ri(firstGas:lastGas) * TempFact) &
                            * Calib(i)%offset(firstGas:lastGas)
                    elsewhere
                        lDrift(firstGas:lastGas) = error
                    end where
                    exit
                end if
            end do
    end select


    !> This call only to calculate chi_h2o, needed for equivalent pressure
    call MoleFractionsAndMixingRatios()

    !> If chi could not be calculated, set it to zero, which in this
    !> context means not accounting for broadening effects. From the site's
    !> water record, not the sixth slot, which is water by convention only.
    !>
    !> PrimaryWaterSlot, not the fallback variant: the fallback answers
    !> histGas2 when a project has no water, and this line would then write a
    !> zero into whatever real trace gas holds that slot. With no hygrometer
    !> there is simply nothing to zero, and broadening stays 1 below.
    wsl = PrimaryWaterSlot()
    if (wsl >= firstGas) then
        if (Stats%chi(wsl) == error) Stats%chi(wsl) = 0d0
    end if

    !> Convert to density/press, correct the absorptance, convert back.
    !>
    !> Two unrolled channels before, co2 and h2o, spelled out across eight
    !> blocks. They differ in exactly two things, and both are properties of
    !> the species rather than of the slot:
    !>
    !>   - the mmol basis, which is water's alone (GasSlotIsWater);
    !>   - the equivalent-pressure factor P_ec = P*[1 + 0.15 chi_h2o], which
    !>     is the water-vapour band-broadening of the *CO2* band (LI-7200
    !>     manual Rev 5). It is not a general property of a trace gas, so it
    !>     is applied where the record says CO2 and nowhere else. The h2o
    !>     channel never had it either.
    !>
    !> A gas whose inverse polynomial is `error` has no calibration curve and
    !> is left untouched - not zeroed, and not run through a polynomial of
    !> error codes, which is what an unguarded loop would do.
    do gas = firstGas, lastGas
        if (gas > ncol) exit
        if (DriftCorr%inv_cal(0, gas) == error) cycle
        if (DriftCorr%dir_cal(0, gas) == error) cycle

        abs_scale = 1d0
        if (GasSlotIsWater(gas)) abs_scale = 1d3

        label = GasOutputLabel(gas)
        call lowercase(label)
        !> Pressure broadening of the CO2 band by water vapour, from the
        !> humidity in the same cell as this CO2 - not the site's. With two
        !> CO2 analysers the site's water is the wrong sample for one of them,
        !> and the effect is a property of the gas mixture the detector sees.
        broadening = 1d0
        msl = E2Col(gas)%moist_ref
        if (trim(adjustl(label)) == 'co2') then
            !> Whichever source this CO2 names, hygrometer or biomet. The
            !> biomet is not the cell's own sample, but it is a measurement
            !> of the air being drawn into it, and it is what the user asked
            !> for by naming it.
            chi_moist = error
            if (msl == biometMoistRef) then
                chi_moist = Ambient%chi_biomet
            elseif (msl >= firstGas .and. msl <= lastGas) then
                chi_moist = Stats%chi(msl)
            end if
            if (chi_moist /= error) &
                broadening = 1d0 + 0.15d0 * chi_moist * 1d-3
        end if

        if (locCol(gas)%measure_type /= 'molar_density') then
            where (Set(:, gas) /= error)
                Set(:, gas) = Set(:, gas) / Ru / Ambient%Tcell &
                    * abs_scale / broadening
            end where
        else
            where (Set(:, gas) /= error)
                Set(:, gas) = Set(:, gas) / (Ambient%Pcell / 1d3 * broadening)
            end where
        end if

        !> Convert densities/press to absorbances/press
        call PolyVal(DriftCorr%inv_cal(0:6, gas), 6, Set(:, gas), &
            size(Set, 1), Set(:, gas))

        !> Calculate the mean absorptance of this column. The whole-Set call
        !> this replaces averaged every column to use one of them; inside a
        !> loop over sixty-four slots that is sixty-four passes over the raw
        !> data per averaging period, for the same number.
        call AverageNoError(Set(:, gas:gas), size(Set, 1), 1, MeanAbs, error)

        !> Remove absorptance/press drift if detected for current period
        !> Note: it can be demonstrated (see Fratini et al. 2014, BG,
        !> Eqs. 10-11), that:
        !> abs_theor = (mean_abs_meas - bias) + d_abs_meas / (1 - bias*P[kPa])
        !> where all terms are intended as normalized by pressure.
        if (lDrift(gas) /= error) then
            where (Set(:, gas) /= error .and. MeanAbs(1) /= error)
                Set(:, gas) = (MeanAbs(1) - lDrift(gas)) + &
                    (Set(:, gas) - MeanAbs(1)) &
                    / (1d0 - lDrift(gas) * Ambient%Pcell * 1d-3)
            end where
        end if

        !> Convert absorptances/press back to density/press
        call PolyVal(DriftCorr%dir_cal(0:6, gas), 6, Set(:, gas), &
            size(Set, 1), Set(:, gas))

        !> Convert density/press back to concentration or density
        if (locCol(gas)%measure_type /= 'molar_density') then
            where (Set(:, gas) /= error)
                Set(:, gas) = Set(:, gas) * Ru * Ambient%Tcell &
                    / abs_scale * broadening
            end where
        else
            where (Set(:, gas) /= error)
                Set(:, gas) = Set(:, gas) * (Ambient%Pcell / 1d3 * broadening)
            end where
        end if
    end do

end subroutine DriftCorrection

!***************************************************************************
!
! \brief       Calculate mean reference counts from raw data if available
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine ReferenceCounts(Set, nrow, ncol)
    use m_rp_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow
    integer, intent(in) :: ncol
    real(kind = dbl), intent(in) :: Set(nrow, ncol)
    !> local variables
    integer :: i
    integer :: slot
    real(kind = dbl) :: meanSet(ncol)
    integer, external :: GasSlotFromDynMDTag


    !> Calculate mean values from raw data
    call AverageNoError(Set, nrow,  ncol, meanSet, error)

    !> Extract mean counts from average, if available.
    !>
    !> Through the same resolver the dynamic metadata file's `<gas>_ref`
    !> columns go through, so the two cannot disagree about which gas
    !> `n2o_ref` names. Spelled out here as two literal cases, they would not
    !> have had to: the signal-strength method reads its reference counts from
    !> the raw data and its offsets from the metadata file, and a mismatch
    !> would silently pair one gas's counts with another's calibration.
    refCounts = error
    do i = 1, NumCol
        slot = GasSlotFromDynMDTag(Col(i)%var, '_ref')
        if (slot > 0) refCounts(slot) = meanSet(i)
    end do
end subroutine ReferenceCounts
