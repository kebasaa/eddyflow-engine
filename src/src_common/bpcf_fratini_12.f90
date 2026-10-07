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
    !> The first full-cospectra file's length. No longer used for sizing -
    !> each period is sized from its own file below - but kept, because the
    !> explicit interface of this routine names it.
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
    !> This period's own full cospectrum, as many rows as its own file has.
    !> It was sized from the first file of the run, so a shorter period - a
    !> slower acquisition rate, or just a gap - left the tail of the array
    !> unread and integrated whatever memory held, and a longer one was cut
    !> at the first file's Nyquist frequency.
    type(SpectraSetType), allocatable :: fullCospectra(:)
    integer :: nrow
    type(BPTFType), allocatable   :: BPTF(:)
    type(BPTFType), allocatable   :: hBPTF(:)
    real(kind = dbl), allocatable :: nf(:)
    real(kind = dbl) :: min_bpcf_f12(GHGNumVar)
    real(kind = dbl) :: max_bpcf_f12(GHGNumVar)
    type(DateType) :: Timestamp
    logical, external :: GasSlotIsWater
    !> Which file holds each period's full cospectra: the list's indices in
    !> timestamp order, built once per list - see FirstFileAt.
    integer, allocatable, save :: Order(:)
    type(DateType), allocatable, save :: SortedTs(:)
    integer, save :: nIndexed = -1
    character(PathLen), save :: IndexedFirst = ''
    character(PathLen), save :: IndexedLast = ''

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

    !> Detect name of file to be read: the first in the list whose timestamp
    !> is this period's end
    call DateTimeToDateType(lEx%end_date, lEx%end_time, Timestamp)
    indx = FirstFileAt(Timestamp)

    !> Set cospectra to be retrieved
    wanted(w_u:w_w) = .false.
    wanted(w_ts) = .true.
!    wanted(firstGas:lastGas) = loc_var_present(co2:w_gas4)
    wanted(firstGas:lastGas) = .false.

    !> Read full co-spectrum of H from file, sized from the file itself.
    !> skip starts true: with no file for this period it was tested unset.
    skip = .true.
    nrow = 0
    if (indx /= nint(error)) then
        call FullCospectraLength(LocFileList(indx)%path, nrow)
        if (nrow > 1) then
            allocate(fullCospectra(nrow))
            call ImportFullCospectra(LocFileList(indx), fullCospectra, nrow, wanted, skip)
        end if
    end if

    if (.not. skip) then
        if (.not. allocated(nf)) allocate(nf(nrow))
        if (.not. allocated(BPTF)) allocate(BPTF(nrow))
        if (.not. allocated(hBPTF)) allocate(hBPTF(nrow))
        nf(1:nrow) = fullCospectra(1:nrow)%fn

        !> Initialize all transfer functions to 1
        call SetTransferFunctionsToValue(BPTF,  nrow, 1d0)
        call SetTransferFunctionsToValue(hBPTF, nrow, 1d0)

        !> File-specific cut-off frequencies
        call RetrieveLPTFpars(lEx, 'iir', LocSetup)

        !> In-situ low-pass transfer function
        call ExperimentalLPTF('iir', nf, nrow, BPTF)

        !> Combined TF (actually only low-pass insitu)
        do gas = firstGas, lastGas
            if (loc_var_present(gas)) &
                call BandPassTransferFunction(BPTF, w, gas, gas, nrow)
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
            do i = 1, nrow
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
                do i = 1, nrow
                    hBPTF(i)%LP%ba_sonic = 1d0
                    hBPTF(i)%LP%ba_irga = 1d0  !< Redundant, but does not harm
                end do
            end if
            if (.not. EddyFlowProj%hf_correct_ghg_zoh) then
                do i = 1, nrow
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
        if (loc_var_present(ts))  call BandPassTransferFunction(hBPTF, w, ts,  w_ts,  nrow)

        !> Apply to measured H co-spectrum to reconstruct a "more theoretical" model co-spectrum
        where (fullCospectra(1:nrow)%of(w_ts) /= error .and. hBPTF(1:nrow)%BP(w_ts) /= error)
            fullCospectra(1:nrow)%of(w_ts) = fullCospectra(1:nrow)%of(w_ts) / hBPTF(1:nrow)%BP(w_ts)
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
                call SpectralCorrectionFactors(fullCospectra%of(w_ts), gas, nf, nrow, BPTF)
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

contains

    !***************************************************************************
    !> \brief The lowest list index whose timestamp equals ts, or the error
    !>        code if none does.
    !>
    !> What a scan of the list from the top finds, which is what this was:
    !> every period walked the whole list, converting the period's own date
    !> again at each step - a year of half-hours against a year of files,
    !> some 300 million steps. The list is fixed for the run, so it is sorted
    !> once and each period is a binary search.
    !>
    !> Sorted by the five fields EqualDates compares, in turn, and ties by list
    !> index - so equal timestamps sit together, lowest index first, and two
    !> are equal here exactly when the scan's == would say so.
    !***************************************************************************
    integer function FirstFileAt(ts)
        type(DateType), intent(in) :: ts
        integer :: lo
        integer :: hi
        integer :: mid

        if (.not. IndexCurrent()) call BuildIndex()
        FirstFileAt = nint(error)
        lo = 1
        hi = nIndexed + 1
        do while (lo < hi)
            mid = (lo + hi) / 2
            if (Before(SortedTs(mid), ts)) then
                lo = mid + 1
            else
                hi = mid
            end if
        end do
        if (lo <= nIndexed) then
            if (SortedTs(lo) == ts) FirstFileAt = Order(lo)
        end if
    end function FirstFileAt

    !> The index was built for this list: same length, same first and last file.
    logical function IndexCurrent()
        IndexCurrent = nIndexed == nfull
        if (.not. IndexCurrent .or. nfull == 0) return
        IndexCurrent = IndexedFirst == LocFileList(1)%path &
            .and. IndexedLast == LocFileList(nfull)%path
    end function IndexCurrent

    !> Bottom-up merge sort of the list's indices, by timestamp then index.
    subroutine BuildIndex()
        integer, allocatable :: tmp(:)
        integer :: width
        integer :: left
        integer :: mid_
        integer :: right
        integer :: a
        integer :: b
        integer :: k

        if (allocated(Order)) deallocate(Order)
        if (allocated(SortedTs)) deallocate(SortedTs)
        allocate(Order(max(1, nfull)), SortedTs(max(1, nfull)), tmp(max(1, nfull)))
        do k = 1, nfull
            Order(k) = k
        end do
        width = 1
        do while (width < nfull)
            left = 1
            do while (left <= nfull)
                mid_ = min(left + width - 1, nfull)
                right = min(left + 2 * width - 1, nfull)
                a = left
                b = mid_ + 1
                k = left
                do while (a <= mid_ .and. b <= right)
                    if (Precedes(Order(b), Order(a))) then
                        tmp(k) = Order(b)
                        b = b + 1
                    else
                        tmp(k) = Order(a)
                        a = a + 1
                    end if
                    k = k + 1
                end do
                do while (a <= mid_)
                    tmp(k) = Order(a)
                    a = a + 1
                    k = k + 1
                end do
                do while (b <= right)
                    tmp(k) = Order(b)
                    b = b + 1
                    k = k + 1
                end do
                left = left + 2 * width
            end do
            Order(1:nfull) = tmp(1:nfull)
            width = 2 * width
        end do
        do k = 1, nfull
            SortedTs(k) = LocFileList(Order(k))%timestamp
        end do
        nIndexed = nfull
        if (nfull > 0) then
            IndexedFirst = LocFileList(1)%path
            IndexedLast = LocFileList(nfull)%path
        end if
        deallocate(tmp)
    end subroutine BuildIndex

    !> List entry i before entry j: earlier timestamp, or the same and lower index.
    logical function Precedes(i, j)
        integer, intent(in) :: i
        integer, intent(in) :: j

        if (Before(LocFileList(i)%timestamp, LocFileList(j)%timestamp)) then
            Precedes = .true.
        elseif (Before(LocFileList(j)%timestamp, LocFileList(i)%timestamp)) then
            Precedes = .false.
        else
            Precedes = i < j
        end if
    end function Precedes

    !> d1 before d2 comparing year, month, day, hour, minute in turn.
    logical function Before(d1, d2)
        type(DateType), intent(in) :: d1
        type(DateType), intent(in) :: d2

        if (d1%Year /= d2%Year) then
            Before = d1%Year < d2%Year
        elseif (d1%Month /= d2%Month) then
            Before = d1%Month < d2%Month
        elseif (d1%Day /= d2%Day) then
            Before = d1%Day < d2%Day
        elseif (d1%Hour /= d2%Hour) then
            Before = d1%Hour < d2%Hour
        else
            Before = d1%Minute < d2%Minute
        end if
    end function Before
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
