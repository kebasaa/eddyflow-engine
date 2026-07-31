!***************************************************************************
! bpcf_aux_subs.f90
! -----------------
! Copyright © 2007-2011, Eco2s team, Gerardo Fratini
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
! \brief      Contains auxiliary subroutines for BPCF calculations.
! \author     Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************

!***************************************************************************
! \brief       Set all band-pass TF components to given value.
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine SetTransferFunctionsToValue(BPTF, nfreq, val)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nfreq
    type(BPTFType), intent(out) :: BPTF(nfreq)
    real(kind = dbl), intent(in) :: val
    !> local variables
    integer :: var

    !> Every variable, not the first eight. `BPTF` is `intent(out)`, so a slot
    !> this loop skips is left undefined - and SpectralCorrectionFactors then
    !> finds no usable band-pass value and reports no correction factor. That
    !> is why a gas past the fourth got none under *every* method, analytic
    !> ones included, with nothing to do with the cospectra file.
    do var = u, lastGas
        BPTF(1:nfreq)%HP(var)  = val
        BPTF(1:nfreq)%EXP(var) = val
        BPTF(1:nfreq)%LP(var)  = LPTFType(val, val, val, val, val, val, &
            val, val, val, val, val)
    end do
end subroutine SetTransferFunctionsToValue

!***************************************************************************
! \brief       Calculates total transfer function (band-pass transfer function)
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine BandPassTransferFunction(BPTF, var1, var2, varout, nfreq)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: var1
    integer, intent(in) :: var2
    integer, intent(in) :: varout
    integer, intent(in) :: nfreq
    type(BPTFType), intent(inout) :: BPTF(nfreq)

    BPTF%BP(varout) = BPTF%LP(var1)%dirga     * BPTF%LP(var2)%dirga  &
                    * BPTF%LP(var1)%dsonic    * BPTF%LP(var2)%dsonic &
                    * BPTF%LP(var1)%wirga     * BPTF%LP(var2)%wirga  &
                    * BPTF%LP(var1)%wsonic    * BPTF%LP(var2)%wsonic &
                    * BPTF%LP(var1)%sver      * BPTF%LP(var2)%sver   &
                    * BPTF%LP(var1)%shor      * BPTF%LP(var2)%shor   &
                    * BPTF%LP(var1)%t         * BPTF%LP(var2)%t      &
                    * BPTF%LP(var1)%ba_sonic  * BPTF%LP(var2)%ba_sonic &
                    * BPTF%LP(var1)%ba_irga   * BPTF%LP(var2)%ba_irga  &
                    * BPTF%LP(var1)%zoh_sonic * BPTF%LP(var2)%zoh_sonic &
                    * BPTF%HP(var1)           * BPTF%HP(var2) &
                    * BPTF%EXP(var1)          * BPTF%EXP(var2)
end subroutine BandPassTransferFunction

!***************************************************************************
! \brief       Calculate spectral correction factors for all fluxes
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine SpectralCorrectionFactors(Cosp, var, nf, nfreq, BPTF)
    use m_common_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nfreq
    integer, intent(in) :: var
    real(kind = dbl), intent(in) :: nf(nfreq)
    real(kind = dbl), intent(in) :: Cosp(nfreq)
    type(BPTFType), intent(in) :: BPTF(nfreq)
    !> local variables
    integer :: k = 0! \file        src/bpcf_aux_subs.f90

    integer :: err_cnt = 0
    real(kind = dbl) :: IntCO
    real(kind = dbl) :: IntTFCO
    real(kind = dbl) :: nf_min, nf_max
    real(kind = dbl) :: df


    !> If cospectrum is made up of only error codes \n
    !> set correction factors to error as well
    err_cnt = 0
    do k = 1, nfreq
        if (Cosp(k) == error) err_cnt = err_cnt + 1
    end do
    if (err_cnt == nfreq) then
        BPCF%of(var) = error
        return
    end if

    !> Artificial frequency range, large enough to accommodate all cases
    nf_min = 1d0/5000d0
    nf_max = 100d0

    !> Integrals of cospectrum and filtered cospectrum
    IntCO = 0d0
    IntTFCO = 0d0
    do k = 1, nfreq - 1
        if (nf(k) > nf_min .and. nf(k + 1) < nf_max .and. &
            Cosp(k) /= error .and. BPTF(k)%BP(var) /= error) then
            df = nf(k + 1) - nf(k)
            IntCO = IntCO + Cosp(k) * df
            IntTFCO = IntTFCO + BPTF(k)%BP(var) * Cosp(k) * df
        end if
    end do
    if (IntTFCO /= 0d0) then
        BPCF%of(var) = IntCO / IntTFCO
    else
        BPCF%of(var) = error
    end if
end subroutine SpectralCorrectionFactors

!***************************************************************************
! \brief       Retrieve transfer function parameters, based on current month or RH
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine RetrieveLPTFpars(lEx, tf_shape, LocSetup)
    use m_common_global_var
    implicit none
    integer :: gas
    !> Optional input arguments
    type(ExType), optional, intent(in) :: lEx
    character(*), optional, intent(in) :: tf_shape
    type(FCCsetupType), optional, intent(in) :: LocSetup
    !> local variables
    integer :: RH
    integer :: month
    real(kind = dbl)  :: A
    real(kind = dbl)  :: B
    real(kind = dbl)  :: C
    real(kind = dbl)  :: lRH
    logical, external :: GasSlotIsWater
    integer, external :: PrimaryWaterSlot
    integer :: wsl

    f_c(firstGas:lastGas) = error
    f_2(firstGas:lastGas) = error

    select case (tf_shape)
        case('iir')
            !> calculate H2O cut-off frequency from current RH, using
            !> exponential fit parameters
            !> Every hygrometer, not the one in the fixed slot. The three
            !> exponential coefficients are a single set for the whole project
            !> - a documented limit of the data model, recorded in
            !> fit_rh_to_cutoff - so a second hygrometer reuses the primary's
            !> RH dependence rather than having none at all. Before this, only
            !> slot 6 was given a cut-off and every other water record fell
            !> through to the analytic correction without saying so.
            A = RegPar(dum, dum)%e1
            B = RegPar(dum, dum)%e2
            C = RegPar(dum, dum)%e3
            lRH = lEx%RH * 1d-2
            do gas = firstGas, lastGas
                if (.not. GasSlotIsWater(gas)) cycle
                if (.not. lEx%var_present(gas)) cycle
                f_c(gas) = dexp(A * lRH**2 + B * lRH + C)
            end do
            !> select relevant tranfer function parameters
            !> according to the month, for CO2, CH4, GAS4
            call char2int(lEx%end_date(6:7), month, 2)
            !> Every configured gas but the water slot, whose cutoff comes
            !> from the RH class above. Guarded on the class index: a gas the
            !> spectral assessment never classified has none, and RegPar would
            !> be indexed out of range rather than merely give a wrong number.
            do gas = firstGas, lastGas
                if (GasSlotIsWater(gas)) cycle
                if (.not. lEx%var_present(gas)) cycle
                if (LocSetup%SA%class(gas, month) < 1 .or. &
                    LocSetup%SA%class(gas, month) > MaxGasClasses) cycle
                f_c(gas) = RegPar(gas, LocSetup%SA%class(gas, month))%fc
            end do

        case('sigma')
            !> select relevant tranfer function parameters
            !> according to the RH-class, for H2O
            !> Same rule as the iir arm: every hygrometer, each taking the
            !> primary water's RH-class parameter, since that is the only
            !> water slot the assessment fits.
            wsl = PrimaryWaterSlot()
            if (wsl >= firstGas) then
                do RH = RH10, RH90
                    if(lEx%RH > dfloat(RH)*10d0 - 5d0 &
                       .and. lEx%RH < dfloat(RH)*10d0 + 5d0) then
                        do gas = firstGas, lastGas
                            if (.not. GasSlotIsWater(gas)) cycle
                            if (.not. lEx%var_present(gas)) cycle
                            f_2(gas) = RegPar(wsl, RH)%f2
                        end do
                        exit
                    end if
                end do
            end if
            call char2int(lEx%end_date(6:7), month, 2)
            !> Every configured gas but the water slot, whose cutoff comes
            !> from the RH class above. Guarded on the class index: a gas the
            !> spectral assessment never classified has none, and RegPar would
            !> be indexed out of range rather than merely give a wrong number.
            do gas = firstGas, lastGas
                if (GasSlotIsWater(gas)) cycle
                if (.not. lEx%var_present(gas)) cycle
                if (LocSetup%SA%class(gas, month) < 1 .or. &
                    LocSetup%SA%class(gas, month) > MaxGasClasses) cycle
                f_2(gas) = RegPar(gas, LocSetup%SA%class(gas, month))%f2
            end do
    end select
end subroutine RetrieveLPTFpars

!***************************************************************************
! \brief       Calculates spectral correction factors based on the procedure \n
!              described in Horst, 1997, Boundary-Layer Meteorology 82: 219�233, 1997. \n
!              Based on  analytical cospectra and \n
!              analytical transfer functions, parameterized in-situ.
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine CorrectionFactorsHorst97(lBPCF, lEx)
    use m_common_global_var
    implicit none
    integer :: gas
    !> in/out variables
    type(SpectralType), intent(inout) :: lBPCF
    !> Optional input arguments
    type(ExType), optional, intent(in) :: lEx
    !> local variables
    real(kind = dbl) :: Nm
    real(kind = dbl) :: alpha
    real(kind = dbl) :: zeta
    real(kind = dbl) :: t_c

    !> Define peak frequency and alpha exponent, basing on stability, \n
    !> see Horst (1997, BLM), eq. 11. Note that his eq. 5 is equal to eq. 3 \n
    !> in Ibrom et al (2007, AFM), with 2*pi*t_c = 1/fc. Thus, the fc derived with \n
    !> the procedure described in Ibrom et al. 2007 can be used to calculate \n
    !> tau_c, and hence to derive the BPCF with Horst's method.
    if (lEx%Flux0%zL <= 0d0) then
        Nm = 0.085d0
        alpha = 7d0 / 8d0
    else
        Nm = 2d0 - 1.915d0 / (1d0 + 0.5d0 * lEx%Flux0%zL)
        alpha = 1d0
    end if
    zeta = lEx%instr(sonic)%height - lEx%disp_height

    !> Correction factor per configured gas. f_c must be a real cutoff: an
    !> unclassified gas leaves it at the error sentinel, and dividing by that
    !> would produce a plausible-looking factor from a value that means
    !> "not available".
    do gas = firstGas, lastGas
        if (.not. lEx%var_present(gas)) cycle
        if (f_c(gas) == error .or. f_c(gas) <= 0d0) cycle
        t_c = 1d0 / (2d0 * p * f_c(gas))
        lBPCF%of(gas) = (1d0 + 2d0 * p * lEx%WS * t_c * Nm / zeta )**alpha
    end do
end subroutine CorrectionFactorsHorst97

!***************************************************************************
! \brief       Calculate correction factors, according to Ibrom et al (2007) \n
!              Agr. For. Meteorol., 147, 149-156. \n
!              Integral method. Derive cut-off frequencies from ensemble averages \n
!              of all available spectra, sorted for RH (H2) and months (other gases), eq. 6 \n
!              Derive low-pass correction factors from degraded temperature time-series, eq. 9. \n
!              analytical transfer functions, parameterized in-situ.
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine CorrectionFactorsIbrom07(gas, lBPCF, lEx)
    use m_common_global_var
    implicit none
    !> in/out variables
    !> One gas slot per call. The four do_co2/do_h2o/do_ch4/do_gas4 logicals
    !> this replaces could only ever name the historical four, and every caller
    !> but one already set exactly one of them.
    integer, intent(in) :: gas
    type(SpectralType), intent(inout) :: lBPCF
    !> Optional input arguments
    type(ExType), optional, intent(in) :: lEx

    if (gas < firstGas .or. gas > lastGas) return
    if (.not. lEx%var_present(gas)) return
    !> An unclassified gas leaves f_c at the error sentinel; adding that to the
    !> denominator gives a plausible-looking factor from a missing value.
    if (f_c(gas) == error) return

    if (lEx%Flux0%zL >= 0d0) then
        lBPCF%of(gas) = StPar(1) * lEx%WS / (StPar(2) + f_c(gas)) + 1d0
    else
        lBPCF%of(gas) = UnPar(1) * lEx%WS / (UnPar(2) + f_c(gas)) + 1d0
    end if
end subroutine CorrectionFactorsIbrom07

!***************************************************************************
! \brief       Build up IIR and SIGMA transfer functions, \n
!              with file specific cut-offs
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine ExperimentalLPTF(shape, nf, N, BPTF)
    use m_common_global_var
    implicit none
    integer :: gas
    !> in/out variables
    integer, intent(in) :: N
    real(kind = dbl), intent(in) :: nf(N)
    character(*), intent(in) :: shape
    type(BPTFType), intent(out) :: BPTF(N)

    !> experimental transfer function, Fratini et al. 2012, Eq. 1 and 3
    if (shape == 'iir') then
        do gas = firstGas, lastGas
            if (f_c(gas) /= error) then
                BPTF(1:N)%EXP(gas) = 1d0 / (1d0 + (nf(1:N) / f_c(gas))**2)
            else
                BPTF(:)%EXP(gas) = 1d0
            end if
        end do

    !> experimental transfer function, see Aubinet et al. (2001, AFM)
    elseif (shape == 'sigma') then
        do gas = firstGas, lastGas
            if (f_2(gas) /= error) then
                BPTF(1:N)%EXP(gas) = dexp(-0.346574d0 * (nf(1:N) / f_2(gas))**2)
            else
                BPTF(:)%EXP(gas) = 1d0
            end if
        end do
    end if
end subroutine ExperimentalLPTF
