!***************************************************************************
! write_out_fluxnet_fcc.f90
! -------------------------
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
! \brief       Write results on output files
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine WriteOutFluxnetFcc(lEx)
    use m_fx_global_var
    implicit none
    !> in/out variables
    Type(ExType), intent(in) :: lEx
    character(16000) :: csv_row

    !> local variables
    integer :: var
    integer :: i
    integer :: j
    integer :: gas
    !> The site's water record. LE and ET are one per site and follow it;
    !> they used to read the sixth slot, which is water by convention only.
    integer :: wsl
    integer :: vi
    !> Gases the fixed part of the row carries columns for. Mirrors what
    !> InitFluxnetFile_rp sized its header loops from.
    integer :: n_layout_gas
    !> One field per variable: a filler digit plus the 8 test outcomes. Nine
    !> characters counts tests, not gases, so this width is unaffected by the
    !> gas capacity - unlike lEx%vm_flags, which is the transpose.
    character(9) :: vm97flags(GHGNumVar)
    include '../src_common/interfaces_1.inc'


    !> Gases the fixed part of the row carries columns for.
    n_layout_gas = min(EddyFlowProj%gas_num, MaxNumGases)

    call clearstr(csv_row)
    !> Timestamp
    !> Start/end imestamps
    call AddDatum(csv_row, trim(adjustl(lEx%start_timestamp)), separator)
    call AddDatum(csv_row, trim(adjustl(lEx%end_timestamp)), separator)
    call AddFloatDatumToDataline(lEx%DOY_start, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%DOY_end, csv_row, EddyFlowProj%err_label)

    !> Filename
    call AddCharDatumToDataline(lEx%fname, csv_row, EddyFlowProj%err_label)

    !> Potential radiation and daytime
    call AddFloatDatumToDataline(lEx%RP, csv_row, EddyFlowProj%err_label)
    call AddIntDatumToDataline(lEx%nighttime_int, csv_row, EddyFlowProj%err_label)

    !> Number of records
    call AddIntDatumToDataline(lEx%nr_theor, csv_row, EddyFlowProj%err_label)        
    call AddIntDatumToDataline(lEx%nr_files, csv_row, EddyFlowProj%err_label)
    call AddIntDatumToDataline(lEx%nr_after_custom_flags, csv_row, EddyFlowProj%err_label)
    call AddIntDatumToDataline(lEx%nr_after_wdf, csv_row, EddyFlowProj%err_label)
    call AddIntDatumToDataline(lEx%nr(u), csv_row, EddyFlowProj%err_label)
    do var = ts, ts + n_layout_gas
        call AddIntDatumToDataline(lEx%nr(var), csv_row, EddyFlowProj%err_label)
    end do
    call AddIntDatumToDataline(lEx%nr_w(u), csv_row, EddyFlowProj%err_label)
    do var = ts, ts + n_layout_gas
        call AddIntDatumToDataline(lEx%nr_w(var), csv_row, EddyFlowProj%err_label)
    end do

    !> Final fluxes
    call AddFloatDatumToDataline(Flux3%Tau, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Flux3%H, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Flux3%LE, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Flux3%ET, csv_row, EddyFlowProj%err_label)
    do gas = firstGas, ts + n_layout_gas
        call AddFloatDatumToDataline(Flux3%gas(gas), csv_row, EddyFlowProj%err_label, &
            gain=FluxnetGasScale(gas), offset=0d0)
    end do

    !> Random uncertainties
    call AddFloatDatumToDataline(lEx%rand_uncer(u), csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%rand_uncer(ts), csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%rand_uncer_LE, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%rand_uncer_ET, csv_row, EddyFlowProj%err_label)
    do gas = firstGas, ts + n_layout_gas
        call AddFloatDatumToDataline(lEx%rand_uncer(gas), csv_row, EddyFlowProj%err_label, &
            gain=FluxnetGasScale(gas), offset=0d0)
    end do

    !> Storage fluxes
    call AddFloatDatumToDataline(lEx%Stor%H, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%Stor%LE, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%Stor%ET, csv_row, EddyFlowProj%err_label)
    do gas = firstGas, ts + n_layout_gas
        call AddFloatDatumToDataline(lEx%Stor%of(gas), csv_row, EddyFlowProj%err_label)
    end do

    !> Advection fluxes
    do gas = firstGas, ts + n_layout_gas
        if (lEx%rot_w /= error .and. lEx%d(gas) >= 0d0) then
            if (lEx%rot_w /= error .and. lEx%d(gas) /= error) then
                call AddFloatDatumToDataline(lEx%rot_w * lEx%d(gas), &
                    csv_row, EddyFlowProj%err_label, &
                    gain=FluxnetGasAdvScale(gas), offset=0d0)
            else
                call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
            end if
        else
            call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        end if
    end do

    !> Turbulence and micromet
    !> Unrotated and rotated wind components
    call AddFloatDatumToDataline(lEx%unrot_u, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%unrot_v, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%unrot_w, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%rot_u, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%rot_v, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%rot_w, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%WS, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%MWS, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%WD, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%WD_SIGMA, csv_row, EddyFlowProj%err_label)

    !> Turbulence
    call AddFloatDatumToDataline(Flux3%ustar, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%TKE, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%L, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%zL, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%bowen, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%Tstar, csv_row, EddyFlowProj%err_label)

    !> Thermodynamics
    !> Temperature, pressure, RH, VPD, e, es, etc.
    call AddFloatDatumToDataline(lEx%Ts, csv_row, EddyFlowProj%err_label, gain=1d0, offset=-273.15d0)
    call AddFloatDatumToDataline(lEx%Ta, csv_row, EddyFlowProj%err_label, gain=1d0, offset=-273.15d0)
    call AddFloatDatumToDataline(lEx%Pa, csv_row, EddyFlowProj%err_label, gain=1d-3, offset=0d0)
    call AddFloatDatumToDataline(lEx%RH, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%Va, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%RHO%a, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%RhoCp, csv_row, EddyFlowProj%err_label)
    if (lEx%RHO%a > 0) then
        call AddFloatDatumToDataline(lEx%RhoCp / lEx%RHO%a, csv_row, EddyFlowProj%err_label)
    else
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
    end if

    !> Water
    call AddFloatDatumToDataline(lEx%RHO%w, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%e, csv_row, EddyFlowProj%err_label, gain=1d-2, offset=0d0)
    call AddFloatDatumToDataline(lEx%es, csv_row, EddyFlowProj%err_label, gain=1d-2, offset=0d0)
    call AddFloatDatumToDataline(lEx%Q, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%VPD, csv_row, EddyFlowProj%err_label, gain=1d-2, offset=0d0)
    call AddFloatDatumToDataline(lEx%Tdew, csv_row, EddyFlowProj%err_label, gain=1d0, offset=-273.15d0)
    !> Dry air
    call AddFloatDatumToDataline(lEx%Pd, csv_row, EddyFlowProj%err_label, gain=1d-3, offset=0d0)
    call AddFloatDatumToDataline(lEx%RHO%d, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%Vd, csv_row, EddyFlowProj%err_label)
    !> Specific heat of evaporation
    call AddFloatDatumToDataline(lEx%lambda, csv_row, EddyFlowProj%err_label)
    !> Wet to dry air density ratio
    call AddFloatDatumToDataline(lEx%sigma, csv_row, EddyFlowProj%err_label)

    !> Gas concentrations/densities
    do gas = firstGas, ts + n_layout_gas
        call AddIntDatumToDataline(lEx%measure_type_int(gas), csv_row, EddyFlowProj%err_label)
        call AddFloatDatumToDataline(lEx%d(gas), csv_row, EddyFlowProj%err_label)
        call AddFloatDatumToDataline(lEx%r(gas), csv_row, EddyFlowProj%err_label, &
            gain=FluxnetGasScale(gas), offset=0d0)
        call AddFloatDatumToDataline(lEx%chi(gas), csv_row, EddyFlowProj%err_label, &
            gain=FluxnetGasScale(gas), offset=0d0)
    end do

    !> Time lags
    do gas = firstGas, ts + n_layout_gas
        call AddFloatDatumToDataline(lEx%act_tlag(gas), csv_row, EddyFlowProj%err_label)
        call AddFloatDatumToDataline(lEx%used_tlag(gas), csv_row, EddyFlowProj%err_label)
        call AddFloatDatumToDataline(lEx%nom_tlag(gas), csv_row, EddyFlowProj%err_label)
        call AddFloatDatumToDataline(lEx%min_tlag(gas), csv_row, EddyFlowProj%err_label)
        call AddFloatDatumToDataline(lEx%max_tlag(gas), csv_row, EddyFlowProj%err_label)
    end do

    !> PWB lag source, one per gas. RP writes these and its header declares
    !> them; FCC used to omit them, which shifted every later column left of
    !> the name above it in the published file.
    do gas = firstGas, ts + n_layout_gas
        call AddFloatDatumToDataline(lEx%pwb_source(gas), csv_row, EddyFlowProj%err_label)
    end do

    !> Stats
    do var = u, ts + n_layout_gas
        if (var == ts) then
            call AddFloatDatumToDataline(lEx%stats%median(var), csv_row, &
                EddyFlowProj%err_label, gain=1d0, offset=-273.15d0)
        else
            call AddFloatDatumToDataline(lEx%stats%median(var), csv_row, EddyFlowProj%err_label)
        end if
    end do
    do var = u, ts + n_layout_gas
        if (var == ts) then
            call AddFloatDatumToDataline(lEx%stats%Q1(var), csv_row, &
                EddyFlowProj%err_label, gain=1d0, offset=-273.15d0)
        else
            call AddFloatDatumToDataline(lEx%stats%Q1(var), csv_row, EddyFlowProj%err_label)
        end if
    end do
    do var = u, ts + n_layout_gas
        if (var == ts) then
            call AddFloatDatumToDataline(lEx%stats%Q3(var), csv_row, &
                EddyFlowProj%err_label, gain=1d0, offset=-273.15d0)
        else
            call AddFloatDatumToDataline(lEx%stats%Q3(var), csv_row, EddyFlowProj%err_label)
        end if
    end do
    do var = u, ts + n_layout_gas
        call AddFloatDatumToDataline(sqrt(lEx%stats%Cov(var, var)), csv_row, EddyFlowProj%err_label)
    end do
    do var = u, ts + n_layout_gas
        call AddFloatDatumToDataline(lEx%stats%Skw(var), csv_row, EddyFlowProj%err_label)
    end do
    do var = u, ts + n_layout_gas
        call AddFloatDatumToDataline(lEx%stats%Kur(var), csv_row, EddyFlowProj%err_label)
    end do
    call AddFloatDatumToDataline(lEx%stats%Cov(w, u), csv_row, EddyFlowProj%err_label)
    do var = ts, ts + n_layout_gas
        call AddFloatDatumToDataline(lEx%stats%Cov(w, var), csv_row, EddyFlowProj%err_label)
    end do
    !> Upper triangle over the configured gases, in the header's pair order.
    do gas = firstGas, ts + n_layout_gas - 1
        do var = gas + 1, ts + n_layout_gas
            call AddFloatDatumToDataline(lEx%stats%Cov(gas, var), csv_row, EddyFlowProj%err_label)
        end do
    end do

    !> Footprint
    call AddFloatDatumToDataline(Foot%peak, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Foot%offset, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Foot%x10, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Foot%x30, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Foot%x50, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Foot%x70, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Foot%x80, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Foot%x90, csv_row, EddyFlowProj%err_label)

    !> Fluxes Level 0 (uncorrected)
    call AddFloatDatumToDataline(lEx%Flux0%ustar, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%Flux0%L, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%Flux0%zL, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%Flux0%Tau, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%Flux0%H, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%Flux0%LE, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%Flux0%ET, csv_row, EddyFlowProj%err_label)
    do gas = firstGas, ts + n_layout_gas
        call AddFloatDatumToDataline(lEx%Flux0%gas(gas), csv_row, EddyFlowProj%err_label, &
            gain=FluxnetGasScale(gas), offset=0d0)
    end do
    !> Fluxes Level 1 
    call AddFloatDatumToDataline(Flux1%Tau, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Flux1%H, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Flux1%LE, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Flux1%ET, csv_row, EddyFlowProj%err_label)
    do gas = firstGas, ts + n_layout_gas
        call AddFloatDatumToDataline(Flux1%gas(gas), csv_row, EddyFlowProj%err_label, &
            gain=FluxnetGasScale(gas), offset=0d0)
    end do
    !> Fluxes Level 2
    call AddFloatDatumToDataline(Flux2%Tau, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Flux2%H, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Flux2%LE, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(Flux2%ET, csv_row, EddyFlowProj%err_label)
    do gas = firstGas, ts + n_layout_gas
        call AddFloatDatumToDataline(Flux2%gas(gas), csv_row, EddyFlowProj%err_label, &
            gain=FluxnetGasScale(gas), offset=0d0)
    end do

    !> Cell values
    call AddFloatDatumToDataline(lEx%Tcell, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%Pcell, csv_row, EddyFlowProj%err_label)
    do gas = firstGas, ts + n_layout_gas
        call AddFloatDatumToDataline(lEx%Vcell(gas), csv_row, EddyFlowProj%err_label)
    end do
    !> Echoed unchanged: these are already SI, unlike the two scalars above.
    do gas = firstGas, ts + n_layout_gas
        call AddFloatDatumToDataline(lEx%Tcell_at(gas), csv_row, EddyFlowProj%err_label)
    end do
    do gas = firstGas, ts + n_layout_gas
        call AddFloatDatumToDataline(lEx%Pcell_at(gas), csv_row, EddyFlowProj%err_label)
    end do
    do gas = firstGas, ts + n_layout_gas
        call AddFloatDatumToDataline(lEx%cov_w_pcell(gas), csv_row, EddyFlowProj%err_label)
    end do
    do gas = firstGas, ts + n_layout_gas
        if (GasSlotIsWater(gas)) cycle
        call AddFloatDatumToDataline(lEx%Flux0%E_gas(gas), csv_row, EddyFlowProj%err_label)
    end do
    do gas = firstGas, ts + n_layout_gas
        call AddFloatDatumToDataline(lEx%Flux0%Hi_gas(gas), csv_row, EddyFlowProj%err_label)
    end do

    !> Burba terms
    if (lEx%Burba%h_bot + lEx%Burba%h_top + lEx%Burba%h_spar /= 0.0) then
        call AddFloatDatumToDataline(lEx%Burba%h_bot, csv_row, EddyFlowProj%err_label)
        call AddFloatDatumToDataline(lEx%Burba%h_top, csv_row, EddyFlowProj%err_label)
        call AddFloatDatumToDataline(lEx%Burba%h_spar, csv_row, EddyFlowProj%err_label)
    else
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
    end if
    !> LI-7700 multipliers
    call AddFloatDatumToDataline(lEx%Mul7700%A, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%Mul7700%B, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%Mul7700%C, csv_row, EddyFlowProj%err_label)

    !> Spectral correction factors. The two water entries are LE and ET,
    !> which take the site's water record - not slot six, which holds water
    !> only when record two happens to.
    wsl = PrimaryWaterOutSlot()
    call AddFloatDatumToDataline(BPCF%of(w_u), csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(BPCF%of(w_ts), csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(BPCF%of(wsl), csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(BPCF%of(wsl), csv_row, EddyFlowProj%err_label)
    do gas = firstGas, ts + n_layout_gas
        call AddFloatDatumToDataline(BPCF%of(gas), csv_row, EddyFlowProj%err_label)
    end do

    !> Degraded covariances
    call AddFloatDatumToDataline(lEx%degT%cov, csv_row, EddyFlowProj%err_label)
    do i = 1, 9
        call AddFloatDatumToDataline(lEx%degT%dcov(i), csv_row, EddyFlowProj%err_label)
    end do
    do var = u, ts + n_layout_gas
        call AddIntDatumToDataline(lEx%spikes(var), csv_row, EddyFlowProj%err_label)
    end do

    !> Write first string from Chunks
    !> M_CUSTOM_FLAGS thru VM97_NSW_RNS
    call AddDatum(csv_row, trim(fluxnetChunks%s(1)), separator)

    !> VM97 flags, here organized per variable instead of per test.
    !>
    !> Re-derived from the transposed form, so this loop has to reproduce
    !> exactly the fields RP wrote: u,v,w,ts then one per configured gas. FCC
    !> has no access to RP's layout lists, but it reads the same project, and
    !> SelectFluxnetGasSlots assigns slot firstGas+k-1 over this same range.
    n_layout_gas = min(EddyFlowProj%gas_num, MaxNumGases)
    if (lEx%vm_flags(1) == '-9999') then
        do j = 1, 4 + n_layout_gas
            call AddCharDatumToDataline(EddyFlowProj%err_label, csv_row, EddyFlowProj%err_label)
        end do
    else
        !> u,v,w,ts then the configured gases - one contiguous range, because
        !> firstGas is ts + 1 and the gas slots run from there.
        do var = u, ts + n_layout_gas
            vi = var + 1
            vm97flags(var)(1 : 1) = '8'
            vm97flags(var)(2 : 2) = lEx%vm_flags(1)(vi:vi)
            vm97flags(var)(3 : 3) = lEx%vm_flags(2)(vi:vi)
            vm97flags(var)(4 : 4) = lEx%vm_flags(3)(vi:vi)
            vm97flags(var)(5 : 5) = lEx%vm_flags(4)(vi:vi)
            vm97flags(var)(6 : 6) = lEx%vm_flags(5)(vi:vi)
            vm97flags(var)(7 : 7) = lEx%vm_flags(6)(vi:vi)
            vm97flags(var)(8 : 8) = lEx%vm_flags(7)(vi:vi)
            vm97flags(var)(9 : 9) = lEx%vm_flags(8)(vi:vi)
            call AddCharDatumToDataline(trim(vm97flags(var)), csv_row, EddyFlowProj%err_label)
        end do
    end if

    !> Uncomment to reintroduce flags for last 3 tests
    !> The third argument is the error label, not the separator.
    !>
    !> These four passed `separator`, so whenever WriteDatumChar decided a cell
    !> was missing it substituted a comma - which is a field boundary, not a
    !> value, so the row gained a column and every field after it shifted. It
    !> stayed hidden because WriteDatumChar's "missing" test is the literal
    !> '899999999', and the time-lag cell only reaches that exact string at
    !> eight gases: one filler digit plus one per gas, all nines for a test
    !> that was not performed. Four and five gases produced a shorter cell and
    !> never tripped it.
    call AddCharDatumToDataline(lEx%vm_tlag_hf, csv_row, EddyFlowProj%err_label)
    call AddCharDatumToDataline(lEx%vm_tlag_sf, csv_row, EddyFlowProj%err_label)
    call AddCharDatumToDataline(lEx%vm_aoa_hf, csv_row, EddyFlowProj%err_label)
    call AddCharDatumToDataline(lEx%vm_nshw_hf, csv_row, EddyFlowProj%err_label)

    !> Write second string from Chunks
    call AddDatum(csv_row, fluxnetChunks%s(2), separator)

    !> Foken's QC details
    call AddFloatDatumToDataline(lEx%TAU_SS, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%H_SS, csv_row, EddyFlowProj%err_label)
    do gas = firstGas, ts + n_layout_gas
        call AddFloatDatumToDataline(lEx%F_SS(gas), csv_row, EddyFlowProj%err_label)
    end do
    call AddFloatDatumToDataline(lEx%U_ITC, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%W_ITC, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%TS_ITC, csv_row, EddyFlowProj%err_label)

    !> Write second string from Chunks
    call AddDatum(csv_row, fluxnetChunks%s(3), separator)

    !> Foken's final flags. Re-derived here rather than read back, so this loop
    !> has to reproduce exactly what RP emitted: momentum, sensible heat and
    !> the two water fluxes, then one per configured gas.
    call AddIntDatumToDataline(QCFlag%tau, csv_row, EddyFlowProj%err_label)
    call AddIntDatumToDataline(QCFlag%H, csv_row, EddyFlowProj%err_label)
    call AddIntDatumToDataline(QCFlag%gas(wsl), csv_row, EddyFlowProj%err_label)
    call AddIntDatumToDataline(QCFlag%gas(wsl), csv_row, EddyFlowProj%err_label)
    do gas = firstGas, ts + n_layout_gas
        call AddIntDatumToDataline(QCFlag%gas(gas), csv_row, EddyFlowProj%err_label)
    end do

    !> LI-COR's IRGAs diagnostics breakdown
    do i = 1, 29
        call AddFloatDatumToDataline(lEx%licor_flags(i), csv_row, EddyFlowProj%err_label)
    end do

    !> AGC/RSSI
    call AddFloatDatumToDataline(lEx%agc72, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%agc75, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%rssi77, csv_row, EddyFlowProj%err_label)

    !> Write third string from Chunks
    !> WBOOST_APPLIED thru AXES_ROTATION_METHOD
    call AddDatum(csv_row, fluxnetChunks%s(4), separator)

    !> Rotation angles
    call AddFloatDatumToDataline(lEx%yaw, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%pitch, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%roll, csv_row, EddyFlowProj%err_label)

    !> Detrending method and time constant
    call AddIntDatumToDataline(lEx%det_meth_int, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%det_timec, csv_row, EddyFlowProj%err_label)

    !> Write forth string from Chunks
    !> TIMELAG_DETECTION_METHOD thru FOOTPRINT_MODEL
    call AddDatum(csv_row, fluxnetChunks%s(5), separator)

    select case(trim(adjustl(foot_model_used)))
    case('none')
        call AddDatum(csv_row, trim(adjustl(EddyFlowProj%err_label)), separator)
    case('kljun_04')
        call AddIntDatumToDataline(0, csv_row, EddyFlowProj%err_label)
    case('kormann_meixner_01')
        call AddIntDatumToDataline(1, csv_row, EddyFlowProj%err_label)
    case('hsieh_00')
        call AddIntDatumToDataline(2, csv_row, EddyFlowProj%err_label)
    end select

    !> Metadata
    call AddIntDatumToDataline(lEx%logger_swver%major, csv_row, EddyFlowProj%err_label)
    call AddIntDatumToDataline(lEx%logger_swver%minor, csv_row, EddyFlowProj%err_label)
    call AddIntDatumToDataline(lEx%logger_swver%revision, csv_row, EddyFlowProj%err_label)
    !>> Site info
    call AddFloatDatumToDataline(lEx%lat, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%lon, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%alt, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%canopy_height, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%disp_height, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%rough_length, csv_row, EddyFlowProj%err_label)
    !>> Acquisition setup
    call AddFloatDatumToDataline(lEx%file_length, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%ac_freq, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%avrg_length, csv_row, EddyFlowProj%err_label)
    !>> Master sonic height and north offset
    call AddDatum(csv_row, trim(lEx%instr(sonic)%firm), separator)
    call AddDatum(csv_row, trim(lEx%instr(sonic)%model), separator)
    call AddFloatDatumToDataline(lEx%instr(sonic)%height, csv_row, EddyFlowProj%err_label)
    call AddDatum(csv_row, lEx%instr(sonic)%wformat, separator)
    call AddDatum(csv_row, lEx%instr(sonic)%wref, separator)
    call AddFloatDatumToDataline(lEx%instr(sonic)%north_offset, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%instr(sonic)%hpath_length, csv_row, EddyFlowProj%err_label, gain=1d2, offset=0d0)
    call AddFloatDatumToDataline(lEx%instr(sonic)%vpath_length, csv_row, EddyFlowProj%err_label, gain=1d2, offset=0d0)
    call AddFloatDatumToDataline(lEx%instr(sonic)%tau, csv_row, EddyFlowProj%err_label)

    !>> Gas analysers, one block per configured gas.
    !>
    !> Re-emitted from the slot-indexed array rather than by instrument role,
    !> which only ever addressed four. ReadExRecord fills that array for every
    !> gas - the historical four by mirroring the role-indexed one - so this
    !> reproduces the block RP wrote, at whatever width the project has.
    do gas = firstGas, firstGas + min(EddyFlowProj%gas_num, MaxNumGases) - 1
        if (gas > lastGas) exit
        call AddDatum(csv_row, trim(lEx%gas_instr(gas)%firm), separator)
        call AddDatum(csv_row, trim(lEx%gas_instr(gas)%model), separator)
        call AddFloatDatumToDataline(lEx%gas_instr(gas)%nsep, csv_row, EddyFlowProj%err_label, gain=1d2, offset=0d0)
        call AddFloatDatumToDataline(lEx%gas_instr(gas)%esep, csv_row, EddyFlowProj%err_label, gain=1d2, offset=0d0)
        call AddFloatDatumToDataline(lEx%gas_instr(gas)%vsep, csv_row, EddyFlowProj%err_label, gain=1d2, offset=0d0)
        call AddFloatDatumToDataline(lEx%gas_instr(gas)%tube_l, csv_row, EddyFlowProj%err_label, gain=1d2, offset=0d0)
        call AddFloatDatumToDataline(lEx%gas_instr(gas)%tube_d, csv_row, EddyFlowProj%err_label, gain=1d3, offset=0d0)
        call AddFloatDatumToDataline(lEx%gas_instr(gas)%tube_f, csv_row, EddyFlowProj%err_label, gain=6d4, offset=0d0)
        if (GasSlotIsWater(gas)) then
            call AddFloatDatumToDataline(lEx%gas_instr(gas)%kw, csv_row, EddyFlowProj%err_label)
            call AddFloatDatumToDataline(lEx%gas_instr(gas)%ko, csv_row, EddyFlowProj%err_label)
        end if
        call AddFloatDatumToDataline(lEx%gas_instr(gas)%hpath_length, csv_row, EddyFlowProj%err_label, gain=1d2, offset=0d0)
        call AddFloatDatumToDataline(lEx%gas_instr(gas)%vpath_length, csv_row, EddyFlowProj%err_label, gain=1d2, offset=0d0)
        call AddFloatDatumToDataline(lEx%gas_instr(gas)%tau, csv_row, EddyFlowProj%err_label)
    end do

    !> Custom variables

    call AddIntDatumToDataline(lEx%ncustom, csv_row, EddyFlowProj%err_label)
    if (lEx%ncustom > 0) then
        do i = 1, lEx%ncustom
            call AddFloatDatumToDataline(lEx%user_var(i), csv_row, EddyFlowProj%err_label)
        end do
    end if

    !> Per-gas water vapour terms, echoed back in the same position and order
    !> as RP wrote them (see WriteOutFluxnet and InitFluxnetFile_rp). FCC
    !> rewrites the whole row, so omitting these would leave the header
    !> describing columns the data no longer has.
    call AddIntDatumToDataline(lEx%n_gas_moist, csv_row, EddyFlowProj%err_label)
    do i = 1, lEx%n_gas_moist
        gas = lEx%gas_moist_slot(i)
        call AddIntDatumToDataline(gas, csv_row, EddyFlowProj%err_label)
        call AddFloatDatumToDataline(lEx%rhow_at(gas), csv_row, EddyFlowProj%err_label)
        call AddFloatDatumToDataline(lEx%sigma_at(gas), csv_row, EddyFlowProj%err_label)
    end do

    !> Analyser of each gas past the four historical slots, replayed exactly as
    !> read (SI units, no conversion) and in the same order.
    call AddIntDatumToDataline(lEx%n_gas_instr, csv_row, EddyFlowProj%err_label)
    do i = 1, lEx%n_gas_instr
        gas = lEx%gas_instr_slot(i)
        call AddIntDatumToDataline(gas, csv_row, EddyFlowProj%err_label)
        call AddCharDatumToDataline(lEx%gas_instr(gas)%firm, csv_row, EddyFlowProj%err_label)
        call AddCharDatumToDataline(lEx%gas_instr(gas)%model, csv_row, EddyFlowProj%err_label)
        call AddFloatDatumToDataline(lEx%gas_instr(gas)%nsep, csv_row, EddyFlowProj%err_label)
        call AddFloatDatumToDataline(lEx%gas_instr(gas)%esep, csv_row, EddyFlowProj%err_label)
        call AddFloatDatumToDataline(lEx%gas_instr(gas)%vsep, csv_row, EddyFlowProj%err_label)
        call AddFloatDatumToDataline(lEx%gas_instr(gas)%tube_l, csv_row, EddyFlowProj%err_label)
        call AddFloatDatumToDataline(lEx%gas_instr(gas)%tube_d, csv_row, EddyFlowProj%err_label)
        call AddFloatDatumToDataline(lEx%gas_instr(gas)%tube_f, csv_row, EddyFlowProj%err_label)
        call AddFloatDatumToDataline(lEx%gas_instr(gas)%hpath_length, csv_row, EddyFlowProj%err_label)
        call AddFloatDatumToDataline(lEx%gas_instr(gas)%vpath_length, csv_row, EddyFlowProj%err_label)
        call AddFloatDatumToDataline(lEx%gas_instr(gas)%tau, csv_row, EddyFlowProj%err_label)
        call AddFloatDatumToDataline(lEx%gas_instr(gas)%kw, csv_row, EddyFlowProj%err_label)
        call AddFloatDatumToDataline(lEx%gas_instr(gas)%ko, csv_row, EddyFlowProj%err_label)
    end do


    !> Write sisxth string from Chunks
    !> Biomet data
    call AddDatum(csv_row, fluxnetChunks%s(6), separator)

    !> Preserve RP's high-frequency CEC descriptor in FCC output.
    call AddFloatDatumToDataline(lEx%cec%r_ET, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%cec%r_Fc, csv_row, EddyFlowProj%err_label)
    call AddIntDatumToDataline(lEx%cec%n_valid, csv_row, EddyFlowProj%err_label)
    call AddIntDatumToDataline(lEx%cec%n_O1, csv_row, EddyFlowProj%err_label)
    call AddIntDatumToDataline(lEx%cec%n_O2, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%cec%frac_O1, csv_row, EddyFlowProj%err_label)
    call AddFloatDatumToDataline(lEx%cec%frac_O2, csv_row, EddyFlowProj%err_label)
    call AddIntDatumToDataline(merge(1, 0, lEx%cec%h2o_valid), csv_row, EddyFlowProj%err_label)
    call AddIntDatumToDataline(merge(1, 0, lEx%cec%co2_valid), csv_row, EddyFlowProj%err_label)
    call AddIntDatumToDataline(lEx%cec%h2o_status, csv_row, EddyFlowProj%err_label)
    call AddIntDatumToDataline(lEx%cec%co2_status, csv_row, EddyFlowProj%err_label)

    !> Replace error codes with user-defined error code
    csv_row = replace2(csv_row, ',-9999,', ',' // trim(EddyFlowProj%err_label) // ',')
    csv_row = replace2(csv_row, ',NaN,',   ',' // trim(EddyFlowProj%err_label) // ',')
    csv_row = replace2(csv_row, ',+Inf,', ',' // trim(EddyFlowProj%err_label) // ',')
    csv_row = replace2(csv_row, ',-Inf,', ',' // trim(EddyFlowProj%err_label) // ',')
    csv_row = replace2(csv_row, ',Inf,', ',' // trim(EddyFlowProj%err_label) // ',')
    csv_row = replace2(csv_row, ',+Infinity,', ',' // trim(EddyFlowProj%err_label) // ',')
    csv_row = replace2(csv_row, ',-Infinity,', ',' // trim(EddyFlowProj%err_label) // ',')
    csv_row = replace2(csv_row, ',Infinity,', ',' // trim(EddyFlowProj%err_label) // ',')

    write(uflxnt, '(a)') csv_row(1:len_trim(csv_row) - 1)

end subroutine WriteOutFluxnetFcc
