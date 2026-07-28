!***************************************************************************
! read_ex_record.f90
! ------------------
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
! \brief       Read one record of essentials file. Based on the requested
!              record number, either reads following record (rec_num < 0)
!              or open the file and look for the actual rec_num
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine ReadExRecord(FilePath, unt, rec_num, lEx, ValidRecord, EndOfFileReached)
    use m_common_global_var
    !> In/out variables
    character(*), intent(in) :: FilePath
    integer, intent(in) :: rec_num
    logical, intent(out) :: ValidRecord
    logical, intent(out) :: EndOfFileReached
    type (ExType), intent(out) :: lEx
    integer, intent(inout) :: unt
    !> Local variables
    integer :: flag
    integer :: gas
    integer :: open_status
    integer :: read_status
    integer :: i
    integer :: var
    integer :: ix
    integer :: cec_h2o_valid
    integer :: cec_co2_valid
    integer :: cec_start
    integer :: remaining_fields
    integer :: n_gas_moist
    real(kind = dbl) :: moist_rhow
    real(kind = dbl) :: moist_sigma
    integer :: n_gas_instr
    integer :: n_gas_extra
    real(kind = dbl) :: xn, xnw, xmt, xd, xr, xchi, xf0, xf3, xf1, xf2, xscf
    real(kind = dbl) :: xru, xstor, xtla, xtlu, xtln, xtlmin, xtlmax, xpwb
    real(kind = dbl) :: xmed, xq1, xq3, xsig, xskw, xkur, xwcov, xvcell, xspk
    character(32) :: instr_firm
    character(32) :: instr_model
    real(kind = dbl) :: instr_nsep, instr_esep, instr_vsep
    real(kind = dbl) :: instr_tube_l, instr_tube_d, instr_tube_f
    real(kind = dbl) :: instr_hpath, instr_vpath, instr_tau
    real(kind = dbl) :: instr_kw, instr_ko
    character(9) :: vm97flags(GHGNumVar)
    character(16000) :: dataline
    character(16000) :: cec_line
    real(kind = dbl) :: aux(32)
    !> Field counts of the ex-file layout.
    !>
    !> This file is navigated positionally: after each list-directed read the
    !> consumed fields are discarded by skipping exactly that many commas. The
    !> counts used to be written as bare literals (263, 23, 12, 38, ...), which
    !> hid the fact that almost every one of them is a function of how many gas
    !> slots the layout carries. Getting one wrong does not raise an error - it
    !> shifts every subsequent field by one and the record is silently misread.
    !>
    !> Expressed here as sums of the groups they cover, in the same order as
    !> the reads below, so the layout is stated once and widening the gas
    !> capacity moves the offsets with it. WriteOutFluxnet and
    !> InitFluxnetFile must emit these same groups in this same order.
    integer, parameter :: nExGas = gas4 - co2 + 1     !< gas slots in the layout
    integer, parameter :: nExVar = gas4 - u + 1       !< u,v,w,ts + gas slots
    integer, parameter :: nExScal = gas4 - ts + 1     !< ts + gas slots

    !> Main record: everything from DOY_START through the spike counts.
    integer, parameter :: nMainFields = &
          4                         &  !< DOY_start, DOY_end, fname, RP
        + 2                         &  !< nighttime_int, nr_theor
        + 3                         &  !< nr_files, nr_after_custom_flags, nr_after_wdf
        + 2 * (1 + nExScal)         &  !< nr, nr_w
        + 8                         &  !< skipped: final fluxes
        + 2                         &  !< rand_uncer u, ts
        + 2 + nExGas                &  !< rand_uncer LE/ET, per gas
        + 3 + nExGas                &  !< storage H/LE/ET, per gas
        + 4                         &  !< skipped: advection fluxes
        + 6                         &  !< unrotated and rotated u,v,w
        + 10                        &  !< WS .. Tstar
        + 8                         &  !< Ts .. Cp
        + 6                         &  !< RHO%w .. Tdew
        + 5                         &  !< Pd .. sigma
        + 4 * nExGas                &  !< measure_type, d, r, chi per gas
        + 5 * nExGas                &  !< timelag quintuplet per gas
        + 4                         &  !< skipped: PWB timelag source
        + 3 * nExVar                &  !< median, Q1, Q3
        + 3 * nExVar                &  !< variance, skewness, kurtosis
        + 1 + nExScal               &  !< Cov(w,u), Cov(w, ts:gas4)
        + (nExGas * (nExGas - 1)) / 2 & !< gas-gas covariances
        + 8                         &  !< skipped: footprint
        + 3                         &  !< Flux0 ustar, L, zL
        + 4                         &  !< Flux0 Tau, H, LE, ET
        + nExGas                    &  !< Flux0 per gas
        + 16                        &  !< skipped: fluxes level 1 and 2
        + 2 + nExGas                &  !< Tcell, Pcell, Vcell per gas
        + (nExGas - 1)              &  !< cell E per gas except h2o
        + nExGas                    &  !< cell Hi per gas
        + 3                         &  !< Burba terms
        + 3                         &  !< LI-7700 multipliers
        + 8                         &  !< skipped: spectral correction factors
        + 1 + 9                     &  !< degraded T covariance and its 9 lags
        + nExVar                       !< spike counts

    integer, parameter :: nNrexFields  = 3 + 4 + 4 + 3 * nExGas   !< NREX chunk
    integer, parameter :: nVmFields    = nExVar + 4               !< VM97 flags
    integer, parameter :: nLgdFields   = 3 * nExVar + (4 + nExGas) + (2 + nExGas)
    integer, parameter :: nSsItcFields = (2 + nExGas) + 3         !< SS + ITC
    integer, parameter :: nSsTestFields = (2 + nExGas) + 3        !< SS/ITC tests
    integer, parameter :: nLicorFields = (4 + nExGas) + 29        !< SSITC + IRGA flags
    integer, parameter :: nAgcFields   = 3                        !< AGC/RSSI
    integer, parameter :: nWboostFields = 3                       !< WBOOST .. AXES_ROT
    integer, parameter :: nRotFields   = 5                        !< angles + detrending
    integer, parameter :: nTlagMethFields = 4                     !< TLAG .. SPEC_CORR
    integer, parameter :: nMetaFields  = 13 + 9 + nExGas * 11 + 2 + 1
    integer, parameter :: nCecFields   = 11                       !< CEC descriptor
    !> Per-gas moisture block: a count, then this many fields per gas
    !> (slot, rhow, sigma).
    integer, parameter :: nGasMoistFields = 3
    !> Per-gas analyser block: a count, then this many fields per gas
    !> (slot, firm, model, nsep, esep, vsep, tube_l, tube_d, tube_f,
    !>  hpath, vpath, tau, kw, ko).
    integer, parameter :: nGasInstrFields = 14
    !> Per-gas family block for gases past the four historical slots: a count,
    !> then this many fields per gas.
    integer, parameter :: nGasExtraFields = 29
    include 'interfaces_1.inc'

    ! integer, external :: strCharIndex

    !> If rec_num > 0,open file and moves to the requested record
    if (rec_num > 0) then
        open(udf, file = trim(adjustl(FilePath)), status = 'old', iostat = open_status)
        if (open_status /= 0) call ExceptionHandler(60)
        unt = udf
        !> Skip header and all records until the requested one
        do i = 1, rec_num
            read(unt, *)
        end do
    end if

    !> Read data line
    ValidRecord = .true.
    EndOfFileReached = .false.
    read(unt, '(a)', iostat = read_status) dataline

    !> Controls on what was read
    if (read_status > 0) then
        ValidRecord = .false.
        if (rec_num > 0) close(unt)
        return
    end if
    if (read_status < 0) then
        EndOfFileReached = .true.
        if (rec_num > 0) close(unt)
        return
    end if

    !> Replace error code with -9999
    dataline = replace2(dataline, trim(EddyFlowProj%err_label), '-9999')

    !> Read timestamps and eliminate them from dataline
    lEx%start_timestamp = dataline(1:12)
    dataline = dataline(14: len_trim(dataline))
    lEx%end_timestamp = dataline(1:12)
    dataline = dataline(14: len_trim(dataline))
    lEx%end_date = lEx%end_timestamp(1:4) // '-' // lEx%end_timestamp(5:6) // '-' // lEx%end_timestamp(7:8) 
    lEx%end_time = lEx%end_timestamp(9:10) // ':' // lEx%end_timestamp(11:12)  

    !> Extract some data
    read(dataline, *, iostat = read_status) lEx%DOY_start, lEx%DOY_end, lEx%fname, lEx%RP, &
        lEx%nighttime_int, lEx%nr_theor, &
        lEx%nr_files, lEx%nr_after_custom_flags, lEx%nr_after_wdf, &
        lEx%nr(u), lEx%nr(ts:gas4), lEx%nr_w(u), lEx%nr_w(ts:gas4), &
        aux(1:8), & !< Skip final fluxes
        lEx%rand_uncer(u), lEx%rand_uncer(ts), &
        lEx%rand_uncer_LE, lEx%rand_uncer_ET, lEx%rand_uncer(co2:gas4), &
        lEx%Stor%H, lEx%Stor%LE, lEx%Stor%ET, lEx%Stor%of(co2:gas4), &
        aux(1:4), & !< Skip advection fluxes
        lEx%unrot_u, lEx%unrot_v, lEx%unrot_w, lEx%rot_u, lEx%rot_v, lEx%rot_w, &
        lEx%WS, lEx%MWS, lEx%WD, lEx%WD_SIGMA, lEx%ustar, lEx%TKE, lEx%L, lEx%zL, lEx%Bowen, lEx%Tstar, &
        lEx%Ts, lEx%Ta, lEx%Pa, lEx%RH, lEx%Va, lEx%RHO%a, lEx%RhoCp, lEx%Cp, &
        lEx%RHO%w, lEx%e, lEx%es, lEx%Q, lEx%VPD, lEx%Tdew, &
        lEx%Pd, lEx%RHO%d, lEx%Vd, lEx%lambda, lEx%sigma, &
        lEx%measure_type_int(co2), lEx%d(co2), lEx%r(co2), lEx%chi(co2), &
        lEx%measure_type_int(h2o), lEx%d(h2o), lEx%r(h2o), lEx%chi(h2o), &
        lEx%measure_type_int(ch4), lEx%d(ch4), lEx%r(ch4), lEx%chi(ch4), &
        lEx%measure_type_int(gas4), lEx%d(gas4), lEx%r(gas4), lEx%chi(gas4), &
        lEx%act_tlag(co2), lEx%used_tlag(co2), lEx%nom_tlag(co2), lEx%min_tlag(co2), lEx%max_tlag(co2), &
        lEx%act_tlag(h2o), lEx%used_tlag(h2o), lEx%nom_tlag(h2o), lEx%min_tlag(h2o), lEx%max_tlag(h2o),&
        lEx%act_tlag(ch4), lEx%used_tlag(ch4), lEx%nom_tlag(ch4), lEx%min_tlag(ch4), lEx%max_tlag(ch4),&
        lEx%act_tlag(gas4), lEx%used_tlag(gas4), lEx%nom_tlag(gas4), lEx%min_tlag(gas4), lEx%max_tlag(gas4), &
        lEx%pwb_source(co2:gas4), &
        lEx%stats%median(u:gas4), lEx%stats%Q1(u:gas4), lEx%stats%Q3(u:gas4), &
        (lEx%stats%Cov(var, var), var=u, gas4), lEx%stats%Skw(u:gas4), lEx%stats%Kur(u:gas4), &
        lEx%stats%Cov(w, u), lEx%stats%Cov(w, ts:gas4), lEx%stats%Cov(co2, h2o:gas4), &
        lEx%stats%Cov(h2o, ch4:gas4), lEx%stats%Cov(ch4, gas4), &
        aux(1:8), & !< Skip footprint
        lEx%Flux0%ustar, lEx%Flux0%L, lEx%Flux0%zL, &
        lEx%Flux0%Tau, lEx%Flux0%H, lEx%Flux0%LE, lEx%Flux0%ET, &
        lEx%Flux0%gas(co2), lEx%Flux0%gas(h2o), lEx%Flux0%gas(ch4), lEx%Flux0%gas(gas4), &
        aux(1:16), & !< Skip fluxes level 1 and 2
        lEx%Tcell, lEx%Pcell, lEx%Vcell(co2:gas4), &
        lEx%Flux0%E_gas(co2), lEx%Flux0%E_gas(ch4), lEx%Flux0%E_gas(gas4), &
        lEx%Flux0%Hi_gas(co2), lEx%Flux0%Hi_gas(h2o), lEx%Flux0%Hi_gas(ch4), lEx%Flux0%Hi_gas(gas4), &
        lEx%Burba%h_bot, lEx%Burba%h_top, lEx%Burba%h_spar, &
        lEx%Mul7700%A, lEx%Mul7700%B, lEx%Mul7700%C, &
        aux(1:8), & !< Skip SCFs
        lEx%degT%cov, lEx%degT%dcov(1:9), &
        lEx%spikes(u:gas4)
    if (read_status /= 0) then
        call InvalidateRecord()
        return
    end if
    if (lEx%fname == 'not_enough_data') then
        call InvalidateRecord()
        return
    end if
    ix = strCharIndex(dataline, ',', nMainFields)
    if (ix <= 0) then
        call InvalidateRecord()
        return
    end if
    dataline = dataline(ix+1: len_trim(dataline))


    !> Copy NREX chunk
    ix = strCharIndex(dataline, ',', nNrexFields)
    if (ix <= 0) then
        call InvalidateRecord()
        return
    end if
    fluxnetChunks%s(1) = dataline(1: ix-1)
    dataline = dataline(ix+1: len_trim(dataline))

    !> Read out VM flags and Foken QC details.
    !> Bound by gas4, not GHGNumVar: WriteOutFluxnet emits one field per
    !> variable over "do var = u, gas4". The two were numerically equal while
    !> GHGNumVar was 8, but they are not once the gas capacity grows, and a
    !> list-directed read of too many items silently swallows the fields that
    !> follow.
    read(dataline, *, iostat = read_status) vm97flags(u:gas4), &
        lEx%vm_tlag_hf, lEx%vm_tlag_sf, lEx%vm_aoa_hf, lEx%vm_nshw_hf
    if (read_status /= 0) then
        call InvalidateRecord()
        return
    end if
    ix = strCharIndex(dataline, ',', nVmFields)
    if (ix <= 0) then
        call InvalidateRecord()
        return
    end if
    dataline = dataline(ix+1: len_trim(dataline))

    !> Rearrage VM flags per test, instead of per variable
    if (vm97flags(u) == '-9999') then
        lEx%vm_flags = '-9999'
    else
        do flag = 1, 8
            lEx%vm_flags(flag)(1:1) = '8'
            lEx%vm_flags(flag)(2:2) = vm97flags(u)(flag + 1: flag + 1)
            lEx%vm_flags(flag)(3:3) = vm97flags(v)(flag + 1: flag + 1)
            lEx%vm_flags(flag)(4:4) = vm97flags(w)(flag + 1: flag + 1)
            lEx%vm_flags(flag)(5:5) = vm97flags(ts)(flag + 1: flag + 1)
            do gas = co2, gas4
                if (vm97flags(gas)(1:1) == '8') then
                    lEx%vm_flags(flag)(gas + 1 : gas + 1) = vm97flags(gas)(flag + 1: flag + 1)
                else
                    lEx%vm_flags(flag)(gas + 1 : gas + 1) = '9'
                end if
            end do
        end do
    end if

    !> Copy LGD/KID/ZCD/CORRDIFF/NSR chunk
    ix = strCharIndex(dataline, ',', nLgdFields)
    if (ix <= 0) then
        call InvalidateRecord()
        return
    end if
    fluxnetChunks%s(2) = dataline(1: ix-1)
    dataline = dataline(ix+1: len_trim(dataline))

    read(dataline, *, iostat = read_status) &
        lEx%TAU_SS, lEx%H_SS, lEx%FC_SS, lEx%FH2O_SS, &
        lEx%FCH4_SS, lEx%FGS4_SS, lEx%U_ITC, lEx%W_ITC, lEx%TS_ITC
    if (read_status /= 0) then
        call InvalidateRecord()
        return
    end if
    ix = strCharIndex(dataline, ',', nSsItcFields)
    if (ix <= 0) then
        call InvalidateRecord()
        return
    end if
    dataline = dataline(ix+1: len_trim(dataline))

    !> Copy another piece
    ix = strCharIndex(dataline, ',', nSsTestFields)
    if (ix <= 0) then
        call InvalidateRecord()
        return
    end if
    fluxnetChunks%s(3) = dataline(1: ix-1)
    dataline = dataline(ix+1: len_trim(dataline))

    !> Read licor IRGA flags
    read(dataline, *, iostat = read_status) aux(1:8), lEx%licor_flags(1:29)
    if (read_status /= 0) then
        call InvalidateRecord()
        return
    end if
    ix = strCharIndex(dataline, ',', nLicorFields)
    if (ix <= 0) then
        call InvalidateRecord()
        return
    end if
    dataline = dataline(ix+1: len_trim(dataline))

    !> Read AGC/RSSI
    read(dataline, *, iostat = read_status) lEx%agc72,lEx%agc75,lEx%rssi77
    if (read_status /= 0) then
        call InvalidateRecord()
        return
    end if
    ix = strCharIndex(dataline, ',', nAgcFields)
    if (ix <= 0) then
        call InvalidateRecord()
        return
    end if
    dataline = dataline(ix+1: len_trim(dataline))

    !> Copy WBOOST_APPLIED thru AXES_ROTATION_METHOD
    ix = strCharIndex(dataline, ',', nWboostFields)
    if (ix <= 0) then
        call InvalidateRecord()
        return
    end if
    fluxnetChunks%s(4) = dataline(1: ix-1)
    dataline = dataline(ix+1: len_trim(dataline))

    !> Read rotation angles and detrending method/time constant
    read(dataline, *, iostat = read_status) &
        lEx%yaw, lEx%pitch, lEx%roll, lEx%det_meth_int, lEx%det_timec
    if (read_status /= 0) then
        call InvalidateRecord()
        return
    end if
    ix = strCharIndex(dataline, ',', nRotFields)
    if (ix <= 0) then
        call InvalidateRecord()
        return
    end if
    dataline = dataline(ix+1: len_trim(dataline))

    !> Copy TIMELAG_DETECTION_METHOD thru SPECTRAL_CORRECTION_METHOD
    ix = strCharIndex(dataline, ',', nTlagMethFields)
    if (ix <= 0) then
        call InvalidateRecord()
        return
    end if
    fluxnetChunks%s(5) = dataline(1: ix-1)
    dataline = dataline(ix+1: len_trim(dataline))

    !> Read out metadata
    read(dataline, *, iostat = read_status) aux(1), &
        lEx%logger_swver%major,lEx%logger_swver%minor,lEx%logger_swver%revision, &
        lEx%lat, lEx%lon, lEx%alt, &
        lEx%canopy_height, lEx%disp_height, lEx%rough_length, &
        lEx%file_length, lEx%ac_freq, lEx%avrg_length, &
        lEx%instr(sonic)%firm, lEx%instr(sonic)%model, lEx%instr(sonic)%height, &
        lEx%instr(sonic)%wformat, lEx%instr(sonic)%wref, lEx%instr(sonic)%north_offset, &
        lEx%instr(sonic)%hpath_length, lEx%instr(sonic)%vpath_length, lEx%instr(sonic)%tau, &
        lEx%instr(ico2)%firm, lEx%instr(ico2)%model, lEx%instr(ico2)%nsep, lEx%instr(ico2)%esep, &
        lEx%instr(ico2)%vsep, lEx%instr(ico2)%tube_l, lEx%instr(ico2)%tube_d, &
        lEx%instr(ico2)%tube_f, &
        lEx%instr(ico2)%hpath_length, lEx%instr(ico2)%vpath_length, lEx%instr(ico2)%tau, &
        lEx%instr(ih2o)%firm, lEx%instr(ih2o)%model, lEx%instr(ih2o)%nsep, lEx%instr(ih2o)%esep, &
        lEx%instr(ih2o)%vsep, lEx%instr(ih2o)%tube_l, lEx%instr(ih2o)%tube_d, &
        lEx%instr(ih2o)%tube_f, lEx%instr(ih2o)%kw, lEx%instr(ih2o)%ko, &
        lEx%instr(ih2o)%hpath_length, lEx%instr(ih2o)%vpath_length, lEx%instr(ih2o)%tau, &
        lEx%instr(ich4)%firm, lEx%instr(ich4)%model, lEx%instr(ich4)%nsep, lEx%instr(ich4)%esep, &
        lEx%instr(ich4)%vsep, lEx%instr(ich4)%tube_l, lEx%instr(ich4)%tube_d, &
        lEx%instr(ich4)%tube_f, &
        lEx%instr(ich4)%hpath_length, lEx%instr(ich4)%vpath_length, lEx%instr(ich4)%tau, &
        lEx%instr(igas4)%firm, lEx%instr(igas4)%model, lEx%instr(igas4)%nsep, lEx%instr(igas4)%esep, &
        lEx%instr(igas4)%vsep, lEx%instr(igas4)%tube_l, lEx%instr(igas4)%tube_d, &
        lEx%instr(igas4)%tube_f, &
        lEx%instr(igas4)%hpath_length, lEx%instr(igas4)%vpath_length, lEx%instr(igas4)%tau, &
        lEx%ncustom
    if (read_status /= 0) then
        call InvalidateRecord()
        return
    end if
    ix = strCharIndex(dataline, ',', nMetaFields)
    if (ix <= 0) then
        call InvalidateRecord()
        return
    end if
    dataline = dataline(ix+1: len_trim(dataline))

    !> Read custom variables
    if (lEx%ncustom > 0) then
        do i = 1, lEx%ncustom
            read(dataline, *, iostat = read_status) lEx%user_var(i)
            if (read_status /= 0) then
                call InvalidateRecord()
                return
            end if
            ix = strCharIndex(dataline, ',', 1)
            if (ix <= 0) then
                call InvalidateRecord()
                return
            end if
            dataline = dataline(ix+1: len_trim(dataline))
        end do
    end if

    !> Read the per-gas water vapour terms written by WriteOutFluxnet.
    !> Self-describing: a count, then that many (slot, rhow, sigma) triples.
    !> The slot is explicit because configured gases need not be contiguous.
    !> Slots not named here keep `error`, which makes the consumers fall back
    !> to the single global sigma/RHO%w.
    lEx%rhow_at = error
    lEx%sigma_at = error
    lEx%n_gas_moist = 0
    lEx%gas_moist_slot = 0
    n_gas_moist = 0
    if (len_trim(dataline) > 0) then
        read(dataline, *, iostat = read_status) n_gas_moist
        if (read_status /= 0) then
            call InvalidateRecord()
            return
        end if
        ix = strCharIndex(dataline, ',', 1)
        if (ix <= 0) then
            call InvalidateRecord()
            return
        end if
        dataline = dataline(ix+1: len_trim(dataline))

        if (n_gas_moist < 0 .or. n_gas_moist > GHGNumVar) then
            call InvalidateRecord()
            return
        end if
        do i = 1, n_gas_moist
            read(dataline, *, iostat = read_status) gas, moist_rhow, moist_sigma
            if (read_status /= 0) then
                call InvalidateRecord()
                return
            end if
            if (gas >= firstGas .and. gas <= lastGas) then
                lEx%rhow_at(gas) = moist_rhow
                lEx%sigma_at(gas) = moist_sigma
                lEx%n_gas_moist = lEx%n_gas_moist + 1
                lEx%gas_moist_slot(lEx%n_gas_moist) = gas
            end if
            ix = strCharIndex(dataline, ',', nGasMoistFields)
            if (ix <= 0) then
                call InvalidateRecord()
                return
            end if
            dataline = dataline(ix+1: len_trim(dataline))
        end do
    end if

    !> Read the analyser of each gas past the four historical slots.
    !> Self-describing like the block above: a count, then that many records
    !> led by the gas slot. Values are SI as written; unlike the GA_* columns
    !> they are not unit-converted here.
    lEx%n_gas_instr = 0
    lEx%gas_instr_slot = 0
    n_gas_instr = 0
    if (len_trim(dataline) > 0) then
        read(dataline, *, iostat = read_status) n_gas_instr
        if (read_status /= 0) then
            call InvalidateRecord()
            return
        end if
        ix = strCharIndex(dataline, ',', 1)
        if (ix <= 0) then
            call InvalidateRecord()
            return
        end if
        dataline = dataline(ix+1: len_trim(dataline))

        if (n_gas_instr < 0 .or. n_gas_instr > GHGNumVar) then
            call InvalidateRecord()
            return
        end if
        do i = 1, n_gas_instr
            read(dataline, *, iostat = read_status) gas, &
                instr_firm, instr_model, instr_nsep, instr_esep, instr_vsep, &
                instr_tube_l, instr_tube_d, instr_tube_f, &
                instr_hpath, instr_vpath, instr_tau, instr_kw, instr_ko
            if (read_status /= 0) then
                call InvalidateRecord()
                return
            end if
            if (gas >= firstGas .and. gas <= lastGas) then
                lEx%gas_instr(gas)%firm = instr_firm
                lEx%gas_instr(gas)%model = instr_model
                lEx%gas_instr(gas)%nsep = instr_nsep
                lEx%gas_instr(gas)%esep = instr_esep
                lEx%gas_instr(gas)%vsep = instr_vsep
                lEx%gas_instr(gas)%tube_l = instr_tube_l
                lEx%gas_instr(gas)%tube_d = instr_tube_d
                lEx%gas_instr(gas)%tube_f = instr_tube_f
                lEx%gas_instr(gas)%hpath_length = instr_hpath
                lEx%gas_instr(gas)%vpath_length = instr_vpath
                lEx%gas_instr(gas)%tau = instr_tau
                lEx%gas_instr(gas)%kw = instr_kw
                lEx%gas_instr(gas)%ko = instr_ko
                lEx%gas_instr(gas)%category = 'irga'
                select case (IrgaPathTypeFromModel(lEx%gas_instr(gas)%model))
                    case ('open')
                        lEx%gas_instr(gas)%path_type = 'open'
                    case default
                        lEx%gas_instr(gas)%path_type = 'closed'
                end select
                if (instr_nsep /= error .and. instr_esep /= error) then
                    lEx%gas_instr(gas)%hsep = dsqrt(instr_nsep**2 + instr_esep**2)
                elseif (instr_nsep /= error) then
                    lEx%gas_instr(gas)%hsep = instr_nsep
                elseif (instr_esep /= error) then
                    lEx%gas_instr(gas)%hsep = instr_esep
                end if
                lEx%n_gas_instr = lEx%n_gas_instr + 1
                lEx%gas_instr_slot(lEx%n_gas_instr) = gas
            end if
            ix = strCharIndex(dataline, ',', nGasInstrFields)
            if (ix <= 0) then
                call InvalidateRecord()
                return
            end if
            dataline = dataline(ix+1: len_trim(dataline))
        end do
    end if

    !> Read the per-gas families of the gases past the four historical slots
    !> straight into the slot-indexed arrays the first four already use, so
    !> downstream code addresses every gas the same way.
    n_gas_extra = 0
    if (len_trim(dataline) > 0) then
        read(dataline, *, iostat = read_status) n_gas_extra
        if (read_status /= 0) then
            call InvalidateRecord()
            return
        end if
        ix = strCharIndex(dataline, ',', 1)
        if (ix <= 0) then
            call InvalidateRecord()
            return
        end if
        dataline = dataline(ix+1: len_trim(dataline))
        if (n_gas_extra < 0 .or. n_gas_extra > GHGNumVar) then
            call InvalidateRecord()
            return
        end if
        do i = 1, n_gas_extra
            read(dataline, *, iostat = read_status) gas, &
                xn, xnw, xmt, xd, xr, xchi, xf0, xf3, xf1, xf2, xscf, xru, xstor, &
                xtla, xtlu, xtln, xtlmin, xtlmax, xpwb, &
                xmed, xq1, xq3, xsig, xskw, xkur, xwcov, xvcell, xspk
            if (read_status /= 0) then
                call InvalidateRecord()
                return
            end if
            if (gas >= firstGas .and. gas <= lastGas) then
                lEx%nr(gas) = nint(xn)
                lEx%nr_w(gas) = nint(xnw)
                lEx%measure_type_int(gas) = nint(xmt)
                lEx%d(gas) = xd
                lEx%r(gas) = xr
                lEx%chi(gas) = xchi
                lEx%Flux0%gas(gas) = xf0
                lEx%rand_uncer(gas) = xru
                lEx%Stor%of(gas) = xstor
                lEx%act_tlag(gas) = xtla
                lEx%used_tlag(gas) = xtlu
                lEx%nom_tlag(gas) = xtln
                lEx%min_tlag(gas) = xtlmin
                lEx%max_tlag(gas) = xtlmax
                lEx%pwb_source(gas) = xpwb
                lEx%stats%median(gas) = xmed
                lEx%stats%Q1(gas) = xq1
                lEx%stats%Q3(gas) = xq3
                lEx%stats%Cov(gas, gas) = xsig
                lEx%stats%Skw(gas) = xskw
                lEx%stats%Kur(gas) = xkur
                lEx%stats%Cov(w, gas) = xwcov
                lEx%Vcell(gas) = xvcell
                lEx%spikes(gas) = nint(xspk)
            end if
            ix = strCharIndex(dataline, ',', nGasExtraFields)
            if (ix <= 0) then
                call InvalidateRecord()
                return
            end if
            dataline = dataline(ix+1: len_trim(dataline))
        end do
    end if

    !> Split and read the CEC descriptor appended after variable biomet data.
    lEx%cec%r_ET = error
    lEx%cec%r_Fc = error
    lEx%cec%frac_O1 = error
    lEx%cec%frac_O2 = error
    lEx%cec%n_valid = 0
    lEx%cec%n_O1 = 0
    lEx%cec%n_O2 = 0
    lEx%cec%h2o_status = cec_rejected
    lEx%cec%co2_status = cec_rejected
    lEx%cec%h2o_valid = .false.
    lEx%cec%co2_valid = .false.
    if (len_trim(dataline) > 0) then
        remaining_fields = count([(dataline(i:i) == ',', &
            i = 1, len_trim(dataline))]) + 1
        cec_line = ''
        if (remaining_fields == nCecFields) then
            cec_line = dataline(1:len_trim(dataline))
            dataline = ''
        elseif (remaining_fields > nCecFields) then
            cec_start = strCharIndex(dataline, ',', remaining_fields - nCecFields)
            if (cec_start > 0) then
                cec_line = dataline(cec_start + 1:len_trim(dataline))
                dataline = dataline(1:cec_start - 1)
            end if
        end if
        if (len_trim(cec_line) > 0) then
            read(cec_line, *, iostat = read_status) lEx%cec%r_ET, lEx%cec%r_Fc, &
                lEx%cec%n_valid, lEx%cec%n_O1, lEx%cec%n_O2, &
                lEx%cec%frac_O1, lEx%cec%frac_O2, cec_h2o_valid, cec_co2_valid, &
                lEx%cec%h2o_status, lEx%cec%co2_status
            if (read_status == 0) then
                lEx%cec%h2o_valid = cec_h2o_valid == 1
                lEx%cec%co2_valid = cec_co2_valid == 1
            end if
        end if
    end if

    !> Put remaining into last chunk
    fluxnetChunks%s(6) = dataline(1: len_trim(dataline))

    !> Complete essentials information based on retrieved ones
    call CompleteEssentials(lEx)

    !> Close file only if it wasn't open on entrance
    if (rec_num > 0) close(unt)
contains

subroutine InvalidateRecord()
    ValidRecord = .false.
    if (rec_num > 0) close(unt)
end subroutine InvalidateRecord

end subroutine ReadExRecord

!***************************************************************************
!
! \brief       Complete essentials information, based on those retrieved \n
!              from the file be useful to other programs
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine CompleteEssentials(lEx)
    use m_common_global_var
    use m_typedef, only: IrgaPathTypeFromModel
    implicit none
    !> in/out variables
    type(ExType), intent(inout) :: lEx
    !> local variables
    integer :: igas
    integer :: gas
    integer :: var

    if (lEx%fname == 'not_enough_data') then
        lEx%not_enough_data = .True.
    else
        lEx%not_enough_data = .False.
    end if

    lEx%var_present = .false.
    if (lEx%WS /= error) lEx%var_present(u:w) = .true.
    if (lEx%Ts /= error) lEx%var_present(ts)  = .true.
    !> Over every gas slot, not just the fixed four: a gas past the fourth
    !> that is left absent here is later gated out of the flux correction.
    do gas = firstGas, lastGas
        if (lEx%Flux0%gas(gas) /= error) lEx%var_present(gas) = .true.
    end do

    !> Units adjustments
    if (lEx%Flux0%gas(ch4) /= error) lEx%Flux0%gas(ch4) = lEx%Flux0%gas(ch4) * 1d-3
    if (lEx%Flux0%gas(gas4) /= error) lEx%Flux0%gas(gas4) = lEx%Flux0%gas(gas4) * 1d-3
    if (lEx%rand_uncer(ch4) /= error) lEx%rand_uncer(ch4) = lEx%rand_uncer(ch4) * 1d-3
    if (lEx%rand_uncer(gas4) /= error) lEx%rand_uncer(gas4) = lEx%rand_uncer(gas4) * 1d-3
    if (lEx%Ts /= error) lEx%Ts = lEx%Ts + 273.15d0
    if (lEx%Ta /= error) lEx%Ta = lEx%Ta + 273.15d0
    if (lEx%Tdew /= error) lEx%Tdew = lEx%Tdew + 273.15d0
    if (lEx%Pa /= error) lEx%Pa = lEx%Pa * 1d3
    if (lEx%Pd /= error) lEx%Pd = lEx%Pd * 1d3
    if (lEx%e /= error) lEx%e = lEx%e * 1d2
    if (lEx%es /= error) lEx%es = lEx%es * 1d2
    if (lEx%Pa /= error) lEx%VPD = lEx%VPD * 1d2
    if (lEx%r(ch4) /= error) lEx%r(ch4) = lEx%r(ch4) * 1d-3
    if (lEx%chi(ch4) /= error) lEx%chi(ch4) = lEx%chi(ch4) * 1d-3
    if (lEx%r(gas4) /= error) lEx%r(gas4) = lEx%r(gas4) * 1d-3
    if (lEx%chi(gas4) /= error) lEx%chi(gas4) = lEx%chi(gas4) * 1d-3
    if (lEx%stats%median(ts) /= error) lEx%stats%median(ts) = lEx%stats%median(ts) + 273.15d0
    if (lEx%stats%Q1(ts) /= error) lEx%stats%Q1(ts) = lEx%stats%Q1(ts) + 273.15d0
    if (lEx%stats%Q3(ts) /= error) lEx%stats%Q3(ts) = lEx%stats%Q3(ts) + 273.15d0
    if (lEx%instr(sonic)%hpath_length /= error) lEx%instr(sonic)%hpath_length = lEx%instr(sonic)%hpath_length * 1d-2
    if (lEx%instr(sonic)%vpath_length /= error) lEx%instr(sonic)%vpath_length = lEx%instr(sonic)%vpath_length * 1d-2
    if (lEx%instr(sonic)%nsep /= error) lEx%instr(sonic)%nsep = lEx%instr(sonic)%nsep * 1d-2
    if (lEx%instr(sonic)%esep /= error) lEx%instr(sonic)%esep = lEx%instr(sonic)%esep * 1d-2
    if (lEx%Burba%h_bot == error) lEx%Burba%h_bot = 0d0
    if (lEx%Burba%h_top == error) lEx%Burba%h_top = 0d0
    if (lEx%Burba%h_spar == error) lEx%Burba%h_spar = 0d0    

    !> Variances were actually read as standard deviations
    do var = u, gas4 
        if (lEx%var_present(var)) &
            lEx%stats%Cov(var, var) = lEx%stats%Cov(var, var)**2
    end do

    lEx%instr(sonic)%category = 'sonic'
    lEx%instr(ico2:igas4)%category = 'irga'
    !> Determine whether igas analysers are open or closed path
    do igas = ico2, igas4
        select case (IrgaPathTypeFromModel(lEx%instr(igas)%model))
            case ('open')
                lEx%instr(igas)%path_type = 'open'
            case default
                lEx%instr(igas)%path_type = 'closed'
                if (lEx%instr(igas)%tube_d /= error) &
                    lEx%instr(igas)%tube_d = lEx%instr(igas)%tube_d * 1d-3
                if (lEx%instr(igas)%tube_l /= error) &
                    lEx%instr(igas)%tube_l = lEx%instr(igas)%tube_l * 1d-2
                if (lEx%instr(igas)%tube_f /= error) &
                    lEx%instr(igas)%tube_f = lEx%instr(igas)%tube_f / 6d4
        end select
        if (lEx%instr(igas)%vsep /= error) lEx%instr(igas)%vsep = lEx%instr(igas)%vsep * 1d-2
        if (lEx%instr(igas)%nsep /= error) lEx%instr(igas)%nsep = lEx%instr(igas)%nsep * 1d-2
        if (lEx%instr(igas)%esep /= error) lEx%instr(igas)%esep = lEx%instr(igas)%esep * 1d-2
        if (lEx%instr(igas)%hpath_length /= error) lEx%instr(igas)%hpath_length = lEx%instr(igas)%hpath_length * 1d-2
        if (lEx%instr(igas)%vpath_length /= error) lEx%instr(igas)%vpath_length = lEx%instr(igas)%vpath_length * 1d-2

        if (lEx%instr(igas)%nsep /= error .and. lEx%instr(igas)%esep /= error) then
            lEx%instr(igas)%hsep = dsqrt(lEx%instr(igas)%nsep**2 + lEx%instr(igas)%esep**2)
        elseif (lEx%instr(igas)%nsep /= error) then
            lEx%instr(igas)%hsep = lEx%instr(igas)%nsep
        elseif (lEx%instr(igas)%esep /= error) then
            lEx%instr(igas)%hsep = lEx%instr(igas)%esep
        end if
    end do

    !> Mirror the four historical analysers into the slot-indexed array, after
    !> the unit conversions above so both views agree. Gases past gas4 were
    !> filled from their own block earlier, already in SI. With this, every
    !> configured gas is addressable the same way and flux code no longer has
    !> to know that instrument role and gas slot are different numberings.
    do igas = ico2, igas4
        lEx%gas_instr(igas - ico2 + firstGas) = lEx%instr(igas)
    end do

    !> Understand software version (AGC (or RSSI) value is negative)
    !> LI-7200
    if (lEx%agc72 < 0 .and. lEx%agc72 /= error) then
        lEx%agc72 =  - lEx%agc72
    else
        co2_new_sw_ver = .true.
    end if
    !> LI-7500
    if (lEx%agc75 < 0 .and. lEx%agc75 /= error) then
        lEx%agc75 =  - lEx%agc75
    else
        co2_new_sw_ver = .true.
    end if

    !> Detrending method from integers to strings
    select case(lEx%det_meth_int)
        case(0)
            lEx%det_meth = 'ba'
        case(1)
            lEx%det_meth = 'ld'
        case(2)
            lEx%det_meth = 'rm'
        case(3)
            lEx%det_meth = 'ew'
    end select

    !> Measurement type from integers to strings. Runs over every gas slot:
    !> the extra-gas block above fills measure_type_int past the fourth, and
    !> the flux code selects on the string form.
    do gas = firstGas, lastGas
        select case(lEx%measure_type_int(gas))
            case(0)
                lEx%measure_type(gas) = 'mixing_ratio'
            case(1)
                lEx%measure_type(gas) = 'mole_fraction'
            case(2)
                lEx%measure_type(gas) = 'molar_density'
        end select
    end do

    !> Daytime
    lEx%daytime = lEx%nighttime_int == 0

    !> Legacy values to be later replaced with newer (left-hand sides) *********
    lEx%file_records = lEx%nr_files
    lEx%used_records = lEx%nr_after_wdf
    lEx%tlag = lEx%used_tlag
    lEx%def_tlag = lEx%used_tlag == lEx%nom_tlag
    do var = u, gas4
        lEx%var(var) = lEx%stats%Cov(var, var)
    end do
    lEx%cov_w(u) = lEx%stats%cov(w, u)
    lEx%cov_w(ts:gas4) = lEx%stats%cov(w, ts:gas4)
end subroutine CompleteEssentials
