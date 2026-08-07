!***************************************************************************
! bpcf_fratini_12.f90
! -------------------
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
! \brief       Calculate spectral correction factors according to \n
!              Fratini et al. 2012 (AFM)
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine BPCF_Fratini12(loc_var_present, LocInstr, wind_speed, t_air, ac_frequency, avrg_length, &
        detrending_time_constant, detrending_method, nfull, nfreq, LocFileList, lEx, LocSetup)
    use m_common_global_var
    implicit none
    logical :: low_flux
    integer :: gas
    !> in/out variables
    logical, intent(in) :: loc_var_present(GHGNumVar)
    type(InstrumentType), intent(in) :: LocInstr(GHGNumVar)
    real(kind = dbl), intent(in) :: wind_speed
    real(kind = dbl), intent(in) :: t_air
    real(kind = dbl), intent(in) :: ac_frequency
    integer, intent(in) :: avrg_length
    integer, intent(in) :: detrending_time_constant
    character(2), intent(in) :: detrending_method
    integer, intent(in) :: nfull
    integer, intent(in) :: nfreq
    !> Optional input arguments
    type(ExType), optional, intent(in) :: lEx
    type(FileListType), optional, intent(in) :: LocFileList(nfull)
    type(FCCsetupType), optional, intent(in) :: LocSetup
    !> local variables
    integer :: i
    integer :: indx
    logical :: wanted(GHGNumVar)
    logical :: skip
    type(SpectraSetType) :: fullCospectra(nfreq)
    type(BPTFType), allocatable   :: BPTF(:)
    type(BPTFType), allocatable   :: hBPTF(:)
    real(kind = dbl), allocatable :: nf(:)
    real(kind = dbl) :: min_bpcf_f12(GHGNumVar)
    real(kind = dbl) :: max_bpcf_f12(GHGNumVar)
    type(DateType) :: Timestamp
    logical, external :: GasSlotIsWater

    !> Plausibility band for the correction factor the direct method returns.
    !>
    !> This was four `data` statements naming co2/h2o/ch4/gas4, so every slot
    !> past the fourth held whatever the saved storage did - zero in practice.
    !> That does not make the test permissive, it inverts it: `BPCF >= 0` is
    !> always true, so a gas past the fourth was pushed onto the Ibrom 2007
    !> fallback in every period and never kept its direct factor.
    !>
    !> The band is a property of the species, not of the slot: water's upper
    !> bound is four times the others because its tube attenuation is that
    !> much larger, so a correction factor that would be absurd for a trace
    !> gas is ordinary for water.
    min_bpcf_f12(:) = 0.8d0
    max_bpcf_f12(:) = 5d0
    do i = firstGas, lastGas
        if (GasSlotIsWater(i)) max_bpcf_f12(i) = 20d0
    end do

    !> Detect name of file to be read
    indx = nint(error)
    do i = 1, nfull
        call DateTimeToDateType(lEx%end_date, lEx%end_time, Timestamp)
        if (LocFileList(i)%timestamp == Timestamp) then
            indx = i
            exit
        end if
    enddo

    !> Set cospectra to be retrieved
    wanted(w_u:w_w) = .false.
    wanted(w_ts) = .true.
!    wanted(firstGas:lastGas) = loc_var_present(co2:w_gas4)
    wanted(firstGas:lastGas) = .false.

    !> Read full co-spectrum of H from file
    if (indx /= nint(error)) &
        call ImportFullCospectra(LocFileList(indx), fullCospectra, nfreq, wanted, skip)

    if (.not. skip) then
        if (.not. allocated(nf)) allocate(nf(nfreq))
        if (.not. allocated(BPTF)) allocate(BPTF(nfreq))
        if (.not. allocated(hBPTF)) allocate(hBPTF(nfreq))
        nf(1:nfreq) = fullCospectra(1:nfreq)%fn

        !> Initialize all transfer functions to 1
        call SetTransferFunctionsToValue(BPTF,  nfreq, 1d0)
        call SetTransferFunctionsToValue(hBPTF, nfreq, 1d0)

        !> File-specific cut-off frequencies
        call RetrieveLPTFpars(lEx, 'iir', LocSetup)

        !> In-situ low-pass transfer function
        call ExperimentalLPTF('iir', nf, nfreq, BPTF)

        !> Combined TF (actually only low-pass insitu)
        do gas = firstGas, lastGas
            if (loc_var_present(gas)) &
                call BandPassTransferFunction(BPTF, w, gas, gas, nfreq)
        end do

        !> Calculate TF to apply to cospectrum of H before using it as a model:
        !> it's applied as H_theor(k) = H_meas(k) / hBPTF%BP(k)
        if (LocSetup%SA%add_sonic_lptf) then
            !> Add analytic components to the experimental transfer function, for \n
            !> sonic path averaging and finite dynamic response
            call AnalyticLowPassTransferFunction(nf, size(nf),  w, LocInstr, &
                loc_var_present, wind_speed, t_air, hBPTF)
            call AnalyticLowPassTransferFunction(nf, size(nf), ts, LocInstr, &
                loc_var_present, wind_speed, t_air, hBPTF)
            !> reset to 1 analytic transfer functions that are substituted by in-situ ones
            do i = 1, nfreq
                hBPTF(i)%LP%t      = 1d0
                hBPTF(i)%LP%wirga  = 1d0
                hBPTF(i)%LP%sver   = 1d0
                hBPTF(i)%LP%shor   = 1d0
                hBPTF(i)%LP%m      = 1d0
                hBPTF(i)%LP%dirga  = 1d0
            end do

            !> Add analytic components to the experimental transfer function, for \n
            !> sonic data filtering in the LI-7550 (BA and ZOH)
            call LI7550_AnalogSignalsTransferFunctions(nf, size(nf), u, &
                ac_frequency, loc_var_present, hBPTF)
            call LI7550_AnalogSignalsTransferFunctions(nf, size(nf), w, &
                ac_frequency, loc_var_present, hBPTF)
            call LI7550_AnalogSignalsTransferFunctions(nf, size(nf), ts, &
                ac_frequency, loc_var_present, hBPTF)
            !> reset to 1 BA and ZOH low-pass transfer functions if the case
            if (.not. EddyFlowProj%hf_correct_ghg_ba) then
                do i = 1, nfreq
                    hBPTF(i)%LP%ba_sonic = 1d0
                    hBPTF(i)%LP%ba_irga = 1d0  !< Redundant, but does not harm
                end do
            end if
            if (.not. EddyFlowProj%hf_correct_ghg_zoh) then
                do i = 1, nfreq
                    hBPTF(i)%LP%zoh_sonic = 1d0
                end do
            end if
        end if

        !> Add analytic high-pass transfer function, if requested
        if (EddyFlowProj%lf_meth == 'analytic') then
            call AnalyticHighPassTransferFunction(nf, size(nf), w, ac_frequency, avrg_length, &
                detrending_method, detrending_time_constant, hBPTF)
            call AnalyticHighPassTransferFunction(nf, size(nf), ts, ac_frequency, avrg_length, &
                detrending_method, detrending_time_constant, hBPTF)
        end if

        !> Combine high-pass TF and sonic-related TF
        if (loc_var_present(ts))  call BandPassTransferFunction(hBPTF, w, ts,  w_ts,  nfreq)

        !> Apply to measured H co-spectrum to reconstruct a "more theoretical" model co-spectrum
        where (fullCospectra(1:nfreq)%of(w_ts) /= error .and. hBPTF(1:nfreq)%BP(w_ts) /= error)
            fullCospectra(1:nfreq)%of(w_ts) = fullCospectra(1:nfreq)%of(w_ts) / hBPTF(1:nfreq)%BP(w_ts)
        end where

        !> Calculate correction factors after Fratini et al. (2012, AFM)
        !>
        !> The model cospectrum is the *measured* w/T one, for every gas: that
        !> substitution is the method. Only the slot being corrected varies,
        !> which is what the loop iterates.
        !>
        !> Indexing the first argument by `gas` instead reads a slot `wanted`
        !> above deliberately excludes from the import, so it is all error -
        !> SpectralCorrectionFactors returns error, the plausibility band below
        !> rejects it, and every gas falls to Ibrom 2007 in every period. The
        !> direct method then never runs, silently.
        do gas = firstGas, lastGas
            if (loc_var_present(gas)) &
                call SpectralCorrectionFactors(fullCospectra%of(w_ts), gas, nf, nfreq, BPTF)
        end do

        !> Calculate correction factors after revision of Laubach and Fratini, unpublished
!        if(loc_var_present(co2)) &
!            call SpectralCorrectionFactorsLaubach(fullCospectra%of(w_co2),  co2,  nf, nfreq, BPTF)
!        if(loc_var_present(h2o)) &
!            call SpectralCorrectionFactorsLaubach(fullCospectra%of(w_h2o),  h2o,  nf, nfreq, BPTF)
!        if(loc_var_present(ch4)) &
!            call SpectralCorrectionFactorsLaubach(fullCospectra%of(w_ch4),  ch4,  nf, nfreq, BPTF)
!        if(loc_var_present(gas4)) &
!            call SpectralCorrectionFactorsLaubach(fullCospectra%of(w_gas4), gas4,  nf, nfreq, BPTF)

        !> Recalculate spectral correction factors following the
        !> approach of Ibrom et al. 2007 (or Fratini et al. 2012, Eq. 4) in the following cases:
        !> 1) Fluxes too low (either sensible heat or concerned gas)
        !> 2) Unrealistic correction factors calculated from direct method
        !> One test per configured gas. Water keeps its own thresholds - the
        !> latent-heat flux and its minimum, not a gas flux - which is the same
        !> carve-out water has everywhere else in this work.
        do gas = firstGas, lastGas
            if (.not. loc_var_present(gas)) cycle
            if (GasSlotIsWater(gas)) then
                !> This hygrometer's own latent heat flux. Screening
                !> every water slot on the site's meant a second
                !> hygrometer was judged on the primary's.
                if (lEx%Flux0%gas(gas) /= error &
                    .and. lEx%lambda /= error) then
                    low_flux = dabs(lEx%Flux0%H) < LocSetup%SA%min_un_H &
                        .or. dabs(lEx%Flux0%gas(gas) * lEx%lambda &
                            * MW_H2O * 1d-3) < LocSetup%SA%min_un_LE
                else
                    low_flux = dabs(lEx%Flux0%H) < LocSetup%SA%min_un_H &
                        .or. dabs(lEx%Flux0%LE) < LocSetup%SA%min_un_LE
                end if
            else
                low_flux = dabs(lEx%Flux0%H) < LocSetup%SA%min_un_H &
                    .or. dabs(lEx%Flux0%gas(gas)) < LocSetup%SA%min_un_gas(gas)
            end if
            if (low_flux .or. BPCF%of(gas) <= min_bpcf_f12(gas) &
                .or. BPCF%of(gas) >= max_bpcf_f12(gas)) &
                call CorrectionFactorsIbrom07(gas, BPCF, lEx)
        end do

        if(allocated(nf)) deallocate(nf)
        if(allocated(BPTF)) deallocate(BPTF)
        if(allocated(hBPTF)) deallocate(hBPTF)
    else
        BPCF%of(:) = 1d0
    end if
end subroutine BPCF_Fratini12

!subroutine SpectralCorrectionFactorsLaubach(Cosp, var, nf, nfreq, BPTF)
!    use m_common_global_var
!    implicit none
!    !> in/out variables
!    integer, intent(in) :: nfreq
!    integer, intent(in) :: var
!    real(kind = dbl), intent(in) :: nf(nfreq)
!    real(kind = dbl), intent(in) :: Cosp(nfreq)
!    type(BPTFType), intent(in) :: BPTF(nfreq)
!    !> local variables
!    integer :: k = 0! \file        src/bpcf_aux_subs.f90
!
!    integer :: err_cnt = 0
!    real(kind = dbl) :: IntCO
!    real(kind = dbl) :: IntTFCO
!    real(kind = dbl) :: nf_min, nf_max
!    real(kind = dbl) :: df
!
!
!    !> If cospectrum is made up of only error codes \n
!    !> set correction factors to error as well
!    err_cnt = 0
!    do k = 1, nfreq
!        if (Cosp(k) == error) err_cnt = err_cnt + 1
!    end do
!    if (err_cnt == nfreq) then
!        BPCF%of(var) = error
!        return
!    end if
!
!    !> Artificial frequency range, large enough to accomodate all cases
!    nf_min = 1d0/5000d0
!    nf_max = 100d0
!
!    !> Integrals of cospectrum and filtered cospectrum
!    IntCO = 0d0
!    IntTFCO = 0d0
!    do k = 1, nfreq - 1
!        if (nf(k) > nf_min .and. nf(k + 1) < nf_max .and. &
!            Cosp(k) /= error .and. BPTF(k)%BP(var) /= error .and. BPTF(k)%BP(var) /= 0d0) then
!            df = nf(k + 1) - nf(k)
!            IntCO = IntCO + Cosp(k) * df
!            IntTFCO = IntTFCO + Cosp(k) / BPTF(k)%BP(var) * df
!        end if
!    end do
!    if (IntTFCO /= 0d0) then
!        BPCF%of(var) = IntTFCO / IntCO
!    else
!        BPCF%of(var) = error
!    end if
!end subroutine SpectralCorrectionFactorsLaubach
