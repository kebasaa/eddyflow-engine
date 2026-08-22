!***************************************************************************
! bpcf_Horst_97.f90
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
! \brief       Calculates spectral correction factors based on the procedure \n
!              described in Horst, 1997, Boundary-Layer Meteorology 82: 219-233, 1997. \n
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
subroutine BPCF_Horst97(measuring_height, displ_height, loc_var_present, wind_speed, zL, &
    ac_frequency, avrg_length, detrending_time_constant, detrending_method, lEx, LocSetup)
    use m_common_global_var
    implicit none
    !> in/out variables
    real(kind = dbl), intent(in) :: measuring_height
    real(kind = dbl), intent(in) :: displ_height
    logical, intent(in) :: loc_var_present(GHGNumVar)
    real(kind = dbl), intent(in) :: wind_speed
    real(kind = dbl), intent(in) :: zL
    real(kind = dbl), intent(in) :: ac_frequency
    integer, intent(in) :: avrg_length
    integer, intent(in) :: detrending_time_constant
    character(2), intent(in) :: detrending_method
    !> Optional input arguments
    type(ExType), optional, intent(in) :: lEx
    type(FCCsetupType), optional, intent(in) :: LocSetup
    !> local variables
    integer :: gas
    integer, parameter :: nseconds = 7200
    integer, parameter :: nfreq = 500
    real(kind = dbl) :: kf(nfreq)
    real(kind = dbl) :: nf(nfreq)
    type(SpectralType) :: Cospectrum(nfreq)
    type(BPTFType) :: BPTF(nfreq)
    type(SpectralType) :: BPCF_Horst
    !> local variables
    integer :: i

    !> add analytic high-pass transfer functions, if requested
    if (EddyFlowProj%lf_meth == 'analytic') then
        !> Log natural frequencies in an artificial freq range
        !> f_min = 1 / 2h --> f_max = 10 Hz
        !> Natural frequency array
        nf(1) = 1d0 / nseconds
        do i = 2, nfreq
            nf(i) = dexp(dlog(nf(1)) + (dlog(10d0)-dlog(nf(1)))/dfloat(nfreq) * dfloat(i))
        end do

        !> Initialize all transfer functions to 1
        call SetTransferFunctionsToValue(BPTF, nfreq, 1d0)

        !> Calculate analytic high-pass transfer functions
        call AnalyticHighPassTransferFunction(nf, size(nf), w, ac_frequency, avrg_length, &
            detrending_method, detrending_time_constant, BPTF)

        do gas = firstGas, lastGas
            if (loc_var_present(gas)) call AnalyticHighPassTransferFunction(nf, size(nf), gas, ac_frequency, avrg_length, &
                detrending_method, detrending_time_constant, BPTF)
        end do

        !> normalized frequency vector, kf
        kf(:) = nf(:) * dabs((measuring_height - displ_height) / wind_speed)

        !> analytical cospectra after Moncrieff et al. (1997, JH)
        call CospectralModel(nf, kf, Cospectrum, zL, nfreq)

        !> combined tf (only high-pass analytic)
        !> One call per configured gas. The w_* covariance labels carry the
        !> same numbering as the gas slots, so the slot indexes both.
        do gas = firstGas, lastGas
            if (loc_var_present(gas)) call BandPassTransferFunction(BPTF, w, gas, gas, nfreq)
        end do

        !> calculate correction factors
        do gas = firstGas, lastGas
            if (loc_var_present(gas)) &
                call SpectralCorrectionFactors(Cospectrum%of(gas), gas, nf, nfreq, BPTF)
        end do
    end if

    !> file-specific cut-off frequencies
    call RetrieveLPTFpars(lEx, 'iir', LocSetup)

    !> calculate correction factors
    BPCF_Horst%of(firstGas:lastGas) = 1d0
    call CorrectionFactorsHorst97(BPCF_Horst, lEx)

    !> Combine low freq (analytic) + high freq (in situ) correction factors
    where (BPCF%of(firstGas:lastGas) /= error .and. BPCF_Horst%of(firstGas:lastGas) /= error)
        BPCF%of(firstGas:lastGas) = BPCF%of(firstGas:lastGas) * BPCF_Horst%of(firstGas:lastGas)
    elsewhere (BPCF_Horst%of(firstGas:lastGas) /= error)
        BPCF%of(firstGas:lastGas) = BPCF_Horst%of(firstGas:lastGas)
    end where
end subroutine BPCF_Horst97
