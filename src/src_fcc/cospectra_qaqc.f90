!*******************************************************************************
! cospectra_qaqc.f90
! ------------------
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
!*******************************************************************************
!
! \brief       Set (co)spectra to error if user-provided quality criteria \n
!              are not met
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!*******************************************************************************
subroutine CospectraQAQC(BinSpec, BinCosp, nrow, lEx, &
    BinCospForStable, BinCospForUnstable, skip_spectra, skip_cospectra)
    use m_fx_global_var
    implicit none
    !> in/out variables
    integer, intent(in) :: nrow
    type(ExType), intent(in) :: lEx
    type(SpectraSetType), intent(inout) :: BinSpec(nrow)
    type(SpectraSetType), intent(inout) :: BinCosp(nrow)
    type(SpectraSetType) :: BinCospForStable(nrow)
    type(SpectraSetType) :: BinCospForUnstable(nrow)
    logical, intent(out) :: skip_spectra
    logical, intent(out) :: skip_cospectra
    !> Local variables
    integer :: i
    !> One digit per variable, taken from the transposed VM97 strings. These
    !> were character(9), which covered u,v,w,ts and exactly four gases, so a
    !> fifth gas could never be flagged here and its cospectra were kept
    !> regardless of what the tests said.
    character(FlagStrLen) :: hf_sr
    character(FlagStrLen) :: hf_do
    character(FlagStrLen) :: hf_sk, sf_sk
    character(FlagStrLen) :: hf_ds, sf_ds
    integer :: STFlg(GHGNumVar)
    integer :: DTFlg(GHGNumVar)
    integer :: qc_tau, qc_H, qc_co2, qc_h2o, qc_ch4, qc_gas4
    integer :: month
    integer :: gas
    integer :: sort
    real(kind = dbl) :: flux
    real(kind = dbl) :: gas_flux
    real(kind = dbl) :: lo, hi
    logical :: all_low
    logical, external :: GasSlotIsWater
    logical :: usable_wt
    logical :: vm_ok(GHGNumVar)
    logical :: foken_ok(GHGNumVar)
    logical :: wind_vm_bad


    !> Initialization
    skip_spectra   = .false.
    skip_cospectra = .false.
    BinCospForStable = BinCosp
    BinCospForUnstable = BinCosp
    usable_wt = any(BinSpec%of(w) /= error .and. BinSpec%of(ts) /= error)
    vm_ok = .true.
    foken_ok = .true.
    if (usable_wt) SADiagUsableWT = SADiagUsableWT + 1

    if (lEx%ustar < FCCsetup%SA%min_un_ustar .or. &
        lEx%ustar > FCCsetup%SA%max_ustar) SADiagRejectedUstar = SADiagRejectedUstar + 1


    !> Razor blade on spectra and co-spectra for unstable case, \n
    !> based on corresponding fluxes.
    !> Unstable case
    if (dabs(lEx%Flux0%H) > FCCsetup%SA%min_un_H &
        .and. dabs(lEx%Flux0%H) < FCCsetup%SA%max_H) then

        !> One test per configured gas, replacing four hand-written blocks
        !> that named co2/ch4/gas4 and water. Water is tested on its latent
        !> heat flux and its own thresholds - the carve-out water has
        !> everywhere else in this work - and every other gas on its own.
        !>
        !> A threshold still at the sentinel means the project never set one:
        !> that test is skipped rather than applied at zero, which would read
        !> as "accept everything" for a minimum and "reject everything" for a
        !> maximum. The gas then contributes no term to the all-fluxes-low
        !> conjunction below either.
        all_low = .true.
        do gas = firstGas, lastGas
            if (gas - firstGas + 1 > &
                min(EddyFlowProj%gas_num, MaxNumGases)) exit
            if (.not. fcc_var_present(gas)) cycle
            if (GasSlotIsWater(gas)) then
                gas_flux = dabs(lEx%Flux0%LE)
                lo = FCCsetup%SA%min_un_LE
                hi = FCCsetup%SA%max_LE
            else
                gas_flux = dabs(lEx%Flux0%gas(gas))
                lo = FCCsetup%SA%min_un_gas(gas)
                hi = FCCsetup%SA%max_gas(gas)
            end if
            if (gas_flux == dabs(error)) cycle
            if ((lo /= error .and. gas_flux < lo) .or. &
                (hi /= error .and. gas_flux > hi)) then
                SADiagRejectedFlux(gas) = SADiagRejectedFlux(gas) + 1
                BinSpec%of(gas) = error
                BinCospForUnstable%of(gas) = error
            end if
            if (lo /= error .and. gas_flux >= lo) all_low = .false.
        end do
        if (all_low) then
            BinSpec = ErrSpec
            BinCospForUnstable = ErrSpec
            skip_spectra = .true.
        end if
    else
        BinSpec = ErrSpec
        BinCospForUnstable = ErrSpec
        skip_spectra = .true.
    end if

    !> Filter co-spectra for u*
    if (lEx%ustar < FCCsetup%SA%min_un_ustar &
        .or. lEx%ustar > FCCsetup%SA%max_ustar) then
        BinCospForUnstable = ErrSpec
        BinSpec = ErrSpec
    end if

    !> Stable case
    if (dabs(lEx%Flux0%H) > FCCsetup%SA%min_st_H &
        .and. dabs(lEx%Flux0%H) < FCCsetup%SA%max_H) then

        !> Same shape as the unstable case above, on the stable thresholds.
        all_low = .true.
        do gas = firstGas, lastGas
            if (gas - firstGas + 1 > &
                min(EddyFlowProj%gas_num, MaxNumGases)) exit
            if (.not. fcc_var_present(gas)) cycle
            if (GasSlotIsWater(gas)) then
                gas_flux = dabs(lEx%Flux0%LE)
                lo = FCCsetup%SA%min_st_LE
                hi = FCCsetup%SA%max_LE
            else
                gas_flux = dabs(lEx%Flux0%gas(gas))
                lo = FCCsetup%SA%min_st_gas(gas)
                hi = FCCsetup%SA%max_gas(gas)
            end if
            if (gas_flux == dabs(error)) cycle
            if ((lo /= error .and. gas_flux < lo) .or. &
                (hi /= error .and. gas_flux > hi)) &
                BinCospForStable%of(gas) = error
            if (lo /= error .and. gas_flux >= lo) all_low = .false.
        end do
        if (all_low) then
            BinCospForStable = ErrSpec
            skip_cospectra = .true.
        end if
    else
        BinCospForStable = ErrSpec
        skip_cospectra = .true.
    end if

    !> Filter co-spectra for u*
    if (lEx%ustar < FCCsetup%SA%min_st_ustar &
        .or. lEx%ustar > FCCsetup%SA%max_ustar) then
        BinCospForStable = ErrSpec
        skip_cospectra = .true.
    end if

    !> Filter based on results of Vickers and Mahrt (1997) quality tests
    !> if requested
    if (FCCsetup%SA%filter_cosp_by_vm_flags) then
        !> Position 1 of each transposed string is the filler digit, so the
        !> variable digits start at 2 and run to the end.
        hf_sr(1:GHGNumVar) = lEx%vm_flags(1)(2:FlagStrLen)
        hf_do(1:GHGNumVar) = lEx%vm_flags(3)(2:FlagStrLen)

        hf_sk(1:GHGNumVar) = lEx%vm_flags(5)(2:FlagStrLen)
        sf_sk(1:GHGNumVar) = lEx%vm_flags(6)(2:FlagStrLen)

        hf_ds(1:GHGNumVar) = lEx%vm_flags(7)(2:FlagStrLen)
        sf_ds(1:GHGNumVar) = lEx%vm_flags(8)(2:FlagStrLen)

        !> If vertical wind speed is flagged, all cospectra are eliminated
        wind_vm_bad = hf_sr(w:w) == '1' .or. hf_do(w:w) == '1' &
            .or. hf_sk(w:w) == '1' .or. hf_ds(w:w) == '1'
        if (wind_vm_bad) then
            BinCospForUnstable = ErrSpec
        end if

        !> Elimination of individual (co)spectra based on the flags on
        !> the relevant variable
        do i = u, lastGas
            if (hf_sr(i:i) == '1' .or. hf_do(i:i) == '1' &
                .or. hf_sk(i:i) == '1' .or. hf_ds(i:i) == '1') then
                if (i >= firstGas) SADiagRejectedVM(i) = SADiagRejectedVM(i) + 1
                BinSpec%of(i) = error
                BinCospForUnstable%of(i) = error
            end if
        end do
        do i = firstGas, lastGas
            vm_ok(i) = .not. wind_vm_bad .and. .not. (hf_sr(i:i) == '1' &
                .or. hf_do(i:i) == '1' .or. hf_sk(i:i) == '1' .or. hf_ds(i:i) == '1')
        end do
    end if

    !> Filter based on results of Foken quality tests if requested.
    !> Regardless of user's choice on how to flag fluxes, here the 0/1/2 scheme
    !> of Mauder and Foken 2004 is used
    if (FCCsetup%SA%foken_lim >= 0) then
        !> Partial flags
        !> Stationarity flags
        call PartialFlagLF(nint(lEx%F_SS(co2)), STFlg(w_co2))
        call PartialFlagLF(nint(lEx%F_SS(h2o)), STFlg(w_h2o))
        call PartialFlagLF(nint(lEx%F_SS(ch4)), STFlg(w_ch4))
        call PartialFlagLF(nint(lEx%F_SS(gas4)), STFlg(w_gas4))
        call PartialFlagLF(nint(lEx%H_SS),  STFlg(w_ts))
        call PartialFlagLF(nint(lEx%TAU_SS),   STFlg(w_u))
        !> Developed turbulence flags
        call PartialFlagLF(nint(lEx%U_ITC), DTFlg(u))
        call PartialFlagLF(nint(lEx%W_ITC), DTFlg(w))
        call PartialFlagLF(nint(lEx%TS_ITC), DTFlg(ts))
        DTFlg(u)  = max(DTFlg(u),  DTFlg(w))

        !> Composite flags
        call GTK2Flag(STFlg(w_u),   DTFlg(u), qc_tau)
        call GTK2Flag(STFlg(w_ts),  DTFlg(w), qc_H)
        call GTK2Flag(STFlg(w_co2), DTFlg(w), qc_co2)
        call GTK2Flag(STFlg(w_h2o), DTFlg(w), qc_h2o)
        call GTK2Flag(STFlg(w_ch4), DTFlg(w), qc_ch4)
        call GTK2Flag(STFlg(w_gas4), DTFlg(w), qc_gas4)

        !> Actual (co)spectra elimination
        if (qc_H < FCCsetup%SA%foken_lim &
            .and. qc_tau < FCCsetup%SA%foken_lim) then
            if (qc_h2o >= FCCsetup%SA%foken_lim) then
                SADiagRejectedFoken(h2o) = SADiagRejectedFoken(h2o) + 1
                BinSpec%of(h2o) = error
                BinCospForUnstable%of(h2o) = error
            end if
            if (qc_co2 >= FCCsetup%SA%foken_lim)  then
                SADiagRejectedFoken(co2) = SADiagRejectedFoken(co2) + 1
                BinSpec%of(co2) = error
                BinCospForUnstable%of(co2) = error
            end if
            if (qc_ch4 >= FCCsetup%SA%foken_lim)  then
                SADiagRejectedFoken(ch4) = SADiagRejectedFoken(ch4) + 1
                BinSpec%of(ch4) = error
                BinCospForUnstable%of(ch4) = error
            end if
            if (qc_gas4 >= FCCsetup%SA%foken_lim) then
                SADiagRejectedFoken(gas4) = SADiagRejectedFoken(gas4) + 1
                BinSpec%of(gas4) = error
                BinCospForUnstable%of(gas4) = error
            end if
            if (qc_h2o >= FCCsetup%SA%foken_lim &
                .and. qc_co2 >= FCCsetup%SA%foken_lim &
                .and. qc_ch4 >= FCCsetup%SA%foken_lim &
                .and. qc_gas4 >= FCCsetup%SA%foken_lim) then
                BinSpec = ErrSpec
                BinCospForUnstable = ErrSpec
                skip_spectra = .true.
            end if
        else
            BinSpec = ErrSpec
            BinCospForUnstable = ErrSpec
            skip_spectra = .true.
        end if
        foken_ok(h2o) = qc_H < FCCsetup%SA%foken_lim .and. qc_tau < FCCsetup%SA%foken_lim &
            .and. qc_h2o < FCCsetup%SA%foken_lim
        foken_ok(co2) = qc_H < FCCsetup%SA%foken_lim .and. qc_tau < FCCsetup%SA%foken_lim &
            .and. qc_co2 < FCCsetup%SA%foken_lim
        foken_ok(ch4) = qc_H < FCCsetup%SA%foken_lim .and. qc_tau < FCCsetup%SA%foken_lim &
            .and. qc_ch4 < FCCsetup%SA%foken_lim
        foken_ok(gas4) = qc_H < FCCsetup%SA%foken_lim .and. qc_tau < FCCsetup%SA%foken_lim &
            .and. qc_gas4 < FCCsetup%SA%foken_lim
    end if

    !> Keep flux candidates that passed every non-flux quality requirement.
    !> These support informational, data-driven threshold suggestions only.
    call char2int(lEx%end_date(6:7), month, 2)
    do i = firstGas, lastGas
        if (.not. lEx%var_present(i) .or. .not. usable_wt) cycle
        if (.not. vm_ok(i) .or. .not. foken_ok(i)) cycle
        sort = 0
        if (i == h2o) then
            if (lEx%RH > 5d0 .and. lEx%RH < 95d0) sort = nint(lEx%RH / 10d0)
            flux = dabs(lEx%Flux0%LE)
        else
            if (month >= JAN .and. month <= DEC) sort = FCCsetup%SA%class(i, month)
            select case (i)
                case (co2)
                    flux = dabs(lEx%Flux0%gas(co2))
                case (ch4)
                    flux = dabs(lEx%Flux0%gas(ch4))
                case default
                    flux = dabs(lEx%Flux0%gas(gas4))
            end select
        end if
        if (flux == dabs(error)) cycle
        if (sort == 0) cycle
        if (dabs(lEx%Flux0%H) > FCCsetup%SA%min_un_H .and. &
            dabs(lEx%Flux0%H) < FCCsetup%SA%max_H .and. &
            lEx%ustar >= FCCsetup%SA%min_un_ustar .and. &
            lEx%ustar <= FCCsetup%SA%max_ustar) &
            call RecordSpectralAssessmentFluxCandidate(i, SADiagUnstable, flux, sort)
        if (dabs(lEx%Flux0%H) > FCCsetup%SA%min_st_H .and. &
            dabs(lEx%Flux0%H) < FCCsetup%SA%max_H .and. &
            lEx%ustar >= FCCsetup%SA%min_st_ustar .and. &
            lEx%ustar <= FCCsetup%SA%max_ustar) &
            call RecordSpectralAssessmentFluxCandidate(i, SADiagStable, flux, sort)
    end do

    do i = firstGas, lastGas
        if (any(BinSpec%of(i) /= error)) SADiagAccepted(i) = SADiagAccepted(i) + 1
    end do

    !> For time sorted cospectra use milder filtering
    BinCosp = BinCospForStable

end subroutine CospectraQAQC
