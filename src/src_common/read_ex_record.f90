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
    integer :: igas
    integer :: open_status
    integer :: read_status
    integer :: i
    integer :: var
    integer :: ix
    !> Gases the fixed part of the row carries a column set for: every gas the
    !> project configures, present or not, which is what InitFluxnetFile_rp
    !> sizes its loops from. Clamped, because a project naming more gases than
    !> the build supports would otherwise index past the layout arrays.
    integer :: n_layout_gas
    integer :: cec_h2o_valid
    integer :: cec_co2_valid
    integer :: cec_start
    integer :: remaining_fields
    integer :: n_gas_moist
    real(kind = dbl) :: moist_rhow
    real(kind = dbl) :: moist_sigma
    integer :: n_gas_instr
    character(32) :: instr_firm
    character(32) :: instr_model
    real(kind = dbl) :: instr_nsep, instr_esep, instr_vsep
    real(kind = dbl) :: instr_tube_l, instr_tube_d, instr_tube_f
    real(kind = dbl) :: instr_hpath, instr_vpath, instr_tau
    real(kind = dbl) :: instr_kw, instr_ko
    character(9) :: vm97flags(GHGNumVar)
    character(16000) :: dataline
    character(16000) :: cec_line
    !> Scratch for fields that are read only to be discarded. Must hold the
    !> largest such group, which is the level 1 and 2 fluxes: two flux families
    !> of four plus one per gas each. At 32 a project with more than 28 gases
    !> would have written past the end of it.
    real(kind = dbl) :: aux(2 * (4 + MaxNumGases))
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
    !> The three group widths the main record is built from. Runtime now:
    !> they used to be gas4-relative constants, which is what pinned the whole
    !> record at four gases.
    integer :: nExGas      !< gas slots in the layout
    integer :: nExWater    !< how many of them are hygrometers
    integer :: nExVar      !< u,v,w,ts + gas slots
    integer :: nExScal     !< ts + gas slots
    !> Last configured gas slot, and the implied-do indices of the main read.
    integer :: lastCfg
    integer :: mgas, mgi, mgj
    !> Cell water flux arrives for every gas but the water slot, so it cannot
    !> be an implied-do; read into this and scatter afterwards. Separate from
    !> aux, which the skip groups in the same read are using.
    real(kind = dbl) :: e_gas_buf(MaxNumGases)

    !> Main record: everything from DOY_START through the spike counts.
    !>
    !> Runtime, like every other per-gas count: thirteen gas families live
    !> inside the single read below, so the moment one of them follows the
    !> configured gas count they all must, and so must this.
    integer :: nMainFields

    !> NREX chunk. Three per-gas runs, so its width is a runtime
    !> quantity; the chunk itself is copied verbatim and never parsed.
    integer :: nNrexFields
    !> VM97 flags: u,v,w,ts, one field per configured gas, then the four
    !> per-test scalars. Runtime for the same reason as the NREX chunk - but
    !> unlike it, these fields are parsed rather than copied, so FCC's re-emit
    !> had to widen with them.
    integer :: nVmFields
    !> LGD/KID/ZCD, then the correlation differences and the Mahrt 98
    !> nonstationarity ratios. Three variable-shaped families and two
    !> flux-shaped ones. Copied verbatim into fluxnetChunks%s(2), so only its
    !> width follows the gas count.
    integer :: nLgdFields
    !> Foken statistics: SS per flux, ITC on u/w/ts. Parsed into lEx%F_SS and
    !> re-emitted by FCC, so this one is not a chunk copy.
    integer :: nSsItcFields
    !> Partial Foken flags: SS per flux, ITC on u/w/ts. Copied into
    !> fluxnetChunks%s(3), so only its width follows the gas count.
    integer :: nSsTestFields
    !> The final Foken flags, which are discarded here, followed by the 29
    !> LI-COR IRGA flags. Only the leading term is per gas: the IRGA flags are
    !> keyed to three instrument models, not to gas slots.
    integer :: nLicorFields
    integer, parameter :: nAgcFields   = 3                        !< AGC/RSSI
    integer, parameter :: nWboostFields = 3                       !< WBOOST .. AXES_ROT
    integer, parameter :: nRotFields   = 5                        !< angles + detrending
    integer, parameter :: nTlagMethFields = 4                     !< TLAG .. SPEC_CORR
    !> Metadata block. The analyser part is per configured gas and therefore
    !> runtime-sized, so only the fixed prefix and the per-gas width are
    !> constants; the old single nMetaFields covered all three at once and
    !> could only ever describe a four-gas layout.
    integer, parameter :: nMetaFixedFields = 13 + 9  !< ident, geometry, sonic
    integer, parameter :: nMetaGasFields = 11        !< +2 for the water slot
    integer :: n_meta_gas
    integer, parameter :: nCecFields   = 11                       !< CEC descriptor
    !> Per-gas moisture block: a count, then this many fields per gas
    !> (slot, rhow, sigma).
    integer, parameter :: nGasMoistFields = 3
    !> Per-gas analyser block: a count, then this many fields per gas
    !> (slot, firm, model, nsep, esep, vsep, tube_l, tube_d, tube_f,
    !>  hpath, vpath, tau, kw, ko).
    integer, parameter :: nGasInstrFields = 14
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

    !> Widths of the main record, from the gas count the project configures.
    n_layout_gas = min(EddyFlowProj%gas_num, MaxNumGases)
    lastCfg = ts + n_layout_gas
    nExGas  = n_layout_gas
    !> How many of the configured gases are water. Two families are sized by
    !> it: the in-cell water flux, which is written for every gas but a
    !> hygrometer, and the krypton pair, which only a hygrometer carries. Both
    !> used to assume exactly one, so a second hygrometer moved the record
    !> width and FCC could no longer read what RP had written.
    nExWater = 0
    do mgas = firstGas, lastCfg
        if (GasSlotIsWater(mgas)) nExWater = nExWater + 1
    end do
    nExVar  = 4 + n_layout_gas
    nExScal = 1 + n_layout_gas
    nMainFields = &
          4                         &  !< DOY_start, DOY_end, fname, RP
        + 2                         &  !< nighttime_int, nr_theor
        + 3                         &  !< nr_files, nr_after_custom_flags, nr_after_wdf
        + 2 * (1 + nExScal)         &  !< nr, nr_w
        + 4 + nExGas                &  !< skipped: final fluxes
        + 2                         &  !< rand_uncer u, ts
        + 2 + nExGas                &  !< rand_uncer LE/ET, per gas
        + 3 + nExGas                &  !< storage H/LE/ET, per gas
        + nExGas                    &  !< skipped: advection fluxes
        + 6                         &  !< unrotated and rotated u,v,w
        + 10                        &  !< WS .. Tstar
        + 8                         &  !< Ts .. Cp
        + 6                         &  !< RHO%w .. Tdew
        + 5                         &  !< Pd .. sigma
        + 4 * nExGas                &  !< measure_type, d, r, chi per gas
        + 5 * nExGas                &  !< timelag quintuplet per gas
        + nExGas                    &  !< PWB timelag source, per gas
        + 3 * nExVar                &  !< median, Q1, Q3
        + 3 * nExVar                &  !< variance, skewness, kurtosis
        + 1 + nExScal               &  !< Cov(w,u), Cov(w, ts:gas4)
        + (nExGas * (nExGas - 1)) / 2 & !< gas-gas covariances
        + 8                         &  !< skipped: footprint
        + 3                         &  !< Flux0 ustar, L, zL
        + 4                         &  !< Flux0 Tau, H, LE, ET
        + nExGas                    &  !< Flux0 per gas
        + 2 * (4 + nExGas)          &  !< skipped: fluxes level 1 and 2
        + 2 + nExGas                &  !< Tcell, Pcell, Vcell per gas
        + 3 * nExGas                &  !< per-gas cell T, cell P, w/cell-P cov
        + (nExGas - nExWater)       &  !< cell E per gas except hygrometers
        + nExGas                    &  !< cell Hi per gas
        + 3                         &  !< Burba terms
        + 3                         &  !< LI-7700 multipliers
        + 4 + nExGas                &  !< skipped: spectral correction factors
        + 1 + 9                     &  !< degraded T covariance and its 9 lags
        + nExVar                       !< spike counts

    !> Extract some data
    read(dataline, *, iostat = read_status) lEx%DOY_start, lEx%DOY_end, lEx%fname, lEx%RP, &
        lEx%nighttime_int, lEx%nr_theor, &
        lEx%nr_files, lEx%nr_after_custom_flags, lEx%nr_after_wdf, &
        lEx%nr(u), lEx%nr(ts:lastCfg), lEx%nr_w(u), lEx%nr_w(ts:lastCfg), &
        aux(1 : 4 + n_layout_gas), & !< Skip final fluxes
        lEx%rand_uncer(u), lEx%rand_uncer(ts), &
        lEx%rand_uncer_LE, lEx%rand_uncer_ET, lEx%rand_uncer(firstGas:lastCfg), &
        lEx%Stor%H, lEx%Stor%LE, lEx%Stor%ET, lEx%Stor%of(firstGas:lastCfg), &
        aux(1 : n_layout_gas), & !< Skip advection fluxes
        lEx%unrot_u, lEx%unrot_v, lEx%unrot_w, lEx%rot_u, lEx%rot_v, lEx%rot_w, &
        lEx%WS, lEx%MWS, lEx%WD, lEx%WD_SIGMA, lEx%ustar, lEx%TKE, lEx%L, lEx%zL, lEx%Bowen, lEx%Tstar, &
        lEx%Ts, lEx%Ta, lEx%Pa, lEx%RH, lEx%Va, lEx%RHO%a, lEx%RhoCp, lEx%Cp, &
        lEx%RHO%w, lEx%e, lEx%es, lEx%Q, lEx%VPD, lEx%Tdew, &
        lEx%Pd, lEx%RHO%d, lEx%Vd, lEx%lambda, lEx%sigma, &
        (lEx%measure_type_int(mgas), lEx%d(mgas), lEx%r(mgas), lEx%chi(mgas), &
            mgas = firstGas, lastCfg), &
        (lEx%act_tlag(mgas), lEx%used_tlag(mgas), lEx%nom_tlag(mgas), &
            lEx%min_tlag(mgas), lEx%max_tlag(mgas), mgas = firstGas, lastCfg), &
        lEx%pwb_source(firstGas:lastCfg), &
        lEx%stats%median(u:lastCfg), lEx%stats%Q1(u:lastCfg), lEx%stats%Q3(u:lastCfg), &
        (lEx%stats%Cov(var, var), var=u, lastCfg), &
        lEx%stats%Skw(u:lastCfg), lEx%stats%Kur(u:lastCfg), &
        lEx%stats%Cov(w, u), lEx%stats%Cov(w, ts:lastCfg), &
        ((lEx%stats%Cov(mgi, mgj), mgj = mgi + 1, lastCfg), &
            mgi = firstGas, lastCfg - 1), &
        aux(1:8), & !< Skip footprint
        lEx%Flux0%ustar, lEx%Flux0%L, lEx%Flux0%zL, &
        lEx%Flux0%Tau, lEx%Flux0%H, lEx%Flux0%LE, lEx%Flux0%ET, &
        lEx%Flux0%gas(firstGas:lastCfg), &
        aux(1 : 2 * (4 + n_layout_gas)), & !< Skip fluxes level 1 and 2
        lEx%Tcell, lEx%Pcell, lEx%Vcell(firstGas:lastCfg), &
        lEx%Tcell_at(firstGas:lastCfg), lEx%Pcell_at(firstGas:lastCfg), &
        lEx%cov_w_pcell(firstGas:lastCfg), &
        !> One field per gas that is *not* a hygrometer, which is what
        !> WriteOutFluxnet emits - it cycles on GasSlotIsWater, so it skips
        !> every water record, not one. Reading `n_layout_gas - 1` assumed
        !> exactly one: on a project with two hygrometers the buffer swallowed
        !> the first field of the block after it and every value from there to
        !> the end of the record came back one field out of step. nMainFields
        !> above already counts this correctly as nExGas - nExWater; only the
        !> read did not.
        e_gas_buf(1 : max(nExGas - nExWater, 0)), &
        lEx%Flux0%Hi_gas(firstGas:lastCfg), &
        lEx%Burba%h_bot, lEx%Burba%h_top, lEx%Burba%h_spar, &
        lEx%Mul7700%A, lEx%Mul7700%B, lEx%Mul7700%C, &
        aux(1 : 4 + n_layout_gas), & !< Skip SCFs
        lEx%degT%cov, lEx%degT%dcov(1:9), &
        lEx%spikes(u:lastCfg)
    if (read_status /= 0) then
        call InvalidateRecord()
        return
    end if
    if (lEx%fname == 'not_enough_data') then
        call InvalidateRecord()
        return
    end if

    !> Scatter the cell water flux back onto its slots. The row carries one
    !> field per gas except the water slot, which no implied-do can express,
    !> so it was read into a buffer above.
    lEx%Flux0%E_gas(firstGas:lastCfg) = error
    mgi = 0
    do mgas = firstGas, lastCfg
        if (GasSlotIsWater(mgas)) cycle
        mgi = mgi + 1
        lEx%Flux0%E_gas(mgas) = e_gas_buf(mgi)
    end do

    ix = strCharIndex(dataline, ',', nMainFields)
    if (ix <= 0) then
        call InvalidateRecord()
        return
    end if
    dataline = dataline(ix+1: len_trim(dataline))


    !> Copy NREX chunk
    nNrexFields = 3 + 4 + 4 + 3 * min(EddyFlowProj%gas_num, MaxNumGases)
    ix = strCharIndex(dataline, ',', nNrexFields)
    if (ix <= 0) then
        call InvalidateRecord()
        return
    end if
    fluxnetChunks%s(1) = dataline(1: ix-1)
    dataline = dataline(ix+1: len_trim(dataline))

    !> Read out VM flags and Foken QC details.
    !> Bound by what WriteOutFluxnet emits: u,v,w,ts then one field per
    !> configured gas. Reading a fixed u:gas4 would consume four gas fields
    !> whatever the project holds - too few for five gases, and for three it
    !> would silently swallow the field that follows.
    n_layout_gas = min(EddyFlowProj%gas_num, MaxNumGases)
    read(dataline, *, iostat = read_status) vm97flags(u:ts), &
        vm97flags(firstGas : firstGas + n_layout_gas - 1), &
        lEx%vm_tlag_hf, lEx%vm_tlag_sf, lEx%vm_aoa_hf, lEx%vm_nshw_hf
    if (read_status /= 0) then
        call InvalidateRecord()
        return
    end if
    nVmFields = 4 + n_layout_gas + 4
    ix = strCharIndex(dataline, ',', nVmFields)
    if (ix <= 0) then
        call InvalidateRecord()
        return
    end if
    dataline = dataline(ix+1: len_trim(dataline))

    !> Rearrage VM flags per test, instead of per variable.
    !>
    !> Transposed: the array index becomes the test and the string position the
    !> variable. Positions are filled to FlagStrLen with '9' - "not performed"
    !> - so a slot the project does not configure reads as absent rather than
    !> as whatever the previous record left there.
    if (vm97flags(u) == '-9999') then
        lEx%vm_flags = '-9999'
    else
        do flag = 1, 8
            lEx%vm_flags(flag) = repeat('9', FlagStrLen)
            lEx%vm_flags(flag)(1:1) = '8'
            lEx%vm_flags(flag)(2:2) = vm97flags(u)(flag + 1: flag + 1)
            lEx%vm_flags(flag)(3:3) = vm97flags(v)(flag + 1: flag + 1)
            lEx%vm_flags(flag)(4:4) = vm97flags(w)(flag + 1: flag + 1)
            lEx%vm_flags(flag)(5:5) = vm97flags(ts)(flag + 1: flag + 1)
            do gas = firstGas, firstGas + n_layout_gas - 1
                if (gas > lastGas) exit
                if (vm97flags(gas)(1:1) == '8') then
                    lEx%vm_flags(flag)(gas + 1 : gas + 1) = vm97flags(gas)(flag + 1: flag + 1)
                else
                    lEx%vm_flags(flag)(gas + 1 : gas + 1) = '9'
                end if
            end do
        end do
    end if

    !> Copy LGD/KID/ZCD/CORRDIFF/NSR chunk
    nLgdFields = 3 * (4 + n_layout_gas) + (4 + n_layout_gas) + (2 + n_layout_gas)
    ix = strCharIndex(dataline, ',', nLgdFields)
    if (ix <= 0) then
        call InvalidateRecord()
        return
    end if
    fluxnetChunks%s(2) = dataline(1: ix-1)
    dataline = dataline(ix+1: len_trim(dataline))

    read(dataline, *, iostat = read_status) &
        lEx%TAU_SS, lEx%H_SS, &
        lEx%F_SS(firstGas : firstGas + n_layout_gas - 1), &
        lEx%U_ITC, lEx%W_ITC, lEx%TS_ITC
    if (read_status /= 0) then
        call InvalidateRecord()
        return
    end if
    nSsItcFields = (2 + n_layout_gas) + 3
    ix = strCharIndex(dataline, ',', nSsItcFields)
    if (ix <= 0) then
        call InvalidateRecord()
        return
    end if
    dataline = dataline(ix+1: len_trim(dataline))

    !> Copy the partial Foken flags
    nSsTestFields = (2 + n_layout_gas) + 3
    ix = strCharIndex(dataline, ',', nSsTestFields)
    if (ix <= 0) then
        call InvalidateRecord()
        return
    end if
    fluxnetChunks%s(3) = dataline(1: ix-1)
    dataline = dataline(ix+1: len_trim(dataline))

    !> Read licor IRGA flags. The final Foken flags come first and are thrown
    !> away - FCC recomputes them - but they still have to be consumed, and
    !> there are now four of them plus one per configured gas.
    read(dataline, *, iostat = read_status) aux(1 : 4 + n_layout_gas), &
        lEx%licor_flags(1:29)
    if (read_status /= 0) then
        call InvalidateRecord()
        return
    end if
    nLicorFields = (4 + n_layout_gas) + 29
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

    !> Read out metadata, in three parts.
    !>
    !> The analyser blocks used to sit in this one list, four of them named
    !> after fixed instrument roles, so a fifth gas had nowhere to be read
    !> from. They are now a loop over the configured gases, matching the
    !> header, which means the field count is a runtime quantity and the
    !> single read had to be split around it.
    read(dataline, *, iostat = read_status) aux(1), &
        lEx%logger_swver%major,lEx%logger_swver%minor,lEx%logger_swver%revision, &
        lEx%lat, lEx%lon, lEx%alt, &
        lEx%canopy_height, lEx%disp_height, lEx%rough_length, &
        lEx%file_length, lEx%ac_freq, lEx%avrg_length, &
        lEx%instr(sonic)%firm, lEx%instr(sonic)%model, lEx%instr(sonic)%height, &
        lEx%instr(sonic)%wformat, lEx%instr(sonic)%wref, lEx%instr(sonic)%north_offset, &
        lEx%instr(sonic)%hpath_length, lEx%instr(sonic)%vpath_length, lEx%instr(sonic)%tau
    if (read_status /= 0) then
        call InvalidateRecord()
        return
    end if
    ix = strCharIndex(dataline, ',', nMetaFixedFields)
    if (ix <= 0) then
        call InvalidateRecord()
        return
    end if
    dataline = dataline(ix+1: len_trim(dataline))

    !> One analyser block per configured gas, in slot order, exactly as
    !> InitFluxnetFile_rp writes them.
    !>
    !> All into gas_instr, converted to SI here and not later.
    !>
    !> The first four used to go through lEx%instr under an instrument-role
    !> index, be converted by a separate loop far below, and be mirrored
    !> across; gases past the fourth were converted inline here by a duplicate
    !> of that arithmetic. Both paths are this one now.
    !>
    !> Converting HERE and not after is what makes it correct: the
    !> self-describing analyser block further down carries the same fields
    !> already in SI and overwrites whatever it lists. The old mirror happened
    !> to run after that block and clobbered slots 5-8 back to converted
    !> values; a conversion pass placed there instead would scale that block's
    !> SI values a second time.
    do i = 1, min(EddyFlowProj%gas_num, MaxNumGases)
        gas = firstGas + i - 1
        if (gas > lastGas) exit
        if (GasSlotIsWater(gas)) then
            read(dataline, *, iostat = read_status) &
                instr_firm, instr_model, instr_nsep, instr_esep, instr_vsep, &
                instr_tube_l, instr_tube_d, instr_tube_f, instr_kw, instr_ko, &
                instr_hpath, instr_vpath, instr_tau
            n_meta_gas = nMetaGasFields + 2
        else
            read(dataline, *, iostat = read_status) &
                instr_firm, instr_model, instr_nsep, instr_esep, instr_vsep, &
                instr_tube_l, instr_tube_d, instr_tube_f, &
                instr_hpath, instr_vpath, instr_tau
            n_meta_gas = nMetaGasFields
        end if
        if (read_status /= 0) then
            call InvalidateRecord()
            return
        end if
        lEx%gas_instr(gas)%firm = instr_firm
        lEx%gas_instr(gas)%model = instr_model
        lEx%gas_instr(gas)%category = 'irga'
        lEx%gas_instr(gas)%tau = instr_tau
        if (GasSlotIsWater(gas)) then
            lEx%gas_instr(gas)%kw = instr_kw
            lEx%gas_instr(gas)%ko = instr_ko
        end if
        !> merge, not a guarded assignment: a gas with no analyser carries the
        !> error code in every field, and that has to reach gas_instr. Guarding
        !> the assignment leaves the field at its initialised zero instead, and
        !> a zero separation reads as a real measurement.
        lEx%gas_instr(gas)%nsep = &
            merge(instr_nsep * 1d-2, instr_nsep, instr_nsep /= error)
        lEx%gas_instr(gas)%esep = &
            merge(instr_esep * 1d-2, instr_esep, instr_esep /= error)
        lEx%gas_instr(gas)%vsep = &
            merge(instr_vsep * 1d-2, instr_vsep, instr_vsep /= error)
        lEx%gas_instr(gas)%hpath_length = &
            merge(instr_hpath * 1d-2, instr_hpath, instr_hpath /= error)
        lEx%gas_instr(gas)%vpath_length = &
            merge(instr_vpath * 1d-2, instr_vpath, instr_vpath /= error)
        select case (IrgaPathTypeFromModel(lEx%gas_instr(gas)%model))
            case ('open')
                lEx%gas_instr(gas)%path_type = 'open'
                lEx%gas_instr(gas)%tube_l = instr_tube_l
                lEx%gas_instr(gas)%tube_d = instr_tube_d
                lEx%gas_instr(gas)%tube_f = instr_tube_f
            case default
                lEx%gas_instr(gas)%path_type = 'closed'
                lEx%gas_instr(gas)%tube_l = &
                    merge(instr_tube_l * 1d-2, instr_tube_l, instr_tube_l /= error)
                lEx%gas_instr(gas)%tube_d = &
                    merge(instr_tube_d * 1d-3, instr_tube_d, instr_tube_d /= error)
                lEx%gas_instr(gas)%tube_f = &
                    merge(instr_tube_f / 6d4, instr_tube_f, instr_tube_f /= error)
        end select
        if (instr_nsep /= error .and. instr_esep /= error) then
            lEx%gas_instr(gas)%hsep = dsqrt(lEx%gas_instr(gas)%nsep**2 &
                + lEx%gas_instr(gas)%esep**2)
        elseif (instr_nsep /= error) then
            lEx%gas_instr(gas)%hsep = lEx%gas_instr(gas)%nsep
        elseif (instr_esep /= error) then
            lEx%gas_instr(gas)%hsep = lEx%gas_instr(gas)%esep
        end if
        ix = strCharIndex(dataline, ',', n_meta_gas)
        if (ix <= 0) then
            call InvalidateRecord()
            return
        end if
        dataline = dataline(ix+1: len_trim(dataline))
    end do

    read(dataline, *, iostat = read_status) lEx%ncustom
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
            !> Past the fourth slot only. The GA_* block above owns co2..gas4
            !> and is the authority for them, exactly as InitFluxnetFile_rp
            !> describes when it writes these columns "for the remaining
            !> gases". The old code enforced that by accident - it mirrored
            !> the converted GA_* values over the top of whatever this block
            !> had left - so with the mirror gone the precedence has to be
            !> stated. Letting this block win for the first four would also
            !> mix units: it is SI, the GA_* columns are metadata units.
            if (gas > gas4 .and. gas >= firstGas .and. gas <= lastGas) then
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
    !> The inverse of WriteOutFluxnet's per-species column scale. This
    !> subroutine has no interfaces include of its own, unlike ReadExRecord.
    real(kind = dbl), external :: FluxnetGasScale

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
    !>
    !> Bounded by the declared gas count, not by lastGas. Flux0 is not reset
    !> between records, so a slot the project never declared holds 0 rather
    !> than the error sentinel and would test as present - a phantom gas that
    !> every var_present-gated loop then processes. It surfaced as the whole
    !> spectral correction falling back to Moncrieff, because the phantom has
    !> no class and so no fitted cutoff to look up.
    do gas = firstGas, lastGas
        if (gas - firstGas + 1 > min(EddyFlowProj%gas_num, MaxNumGases)) exit
        if (lEx%Flux0%gas(gas) /= error) lEx%var_present(gas) = .true.
    end do

    !> Units adjustments.
    !>
    !> The gas quantities are the exact inverse of what WriteOutFluxnet
    !> applied, so they must come from the same function. This was a literal
    !> 1d-3 on the ch4 and gas4 slots, which is only the right inverse while
    !> those slots hold trace gases: a project with CO2 in slot four was
    !> divided by a thousand it had never been multiplied by. CO2 and H2O
    !> scale by 1, so covering every slot changes nothing for them.
    do gas = co2, ts + min(EddyFlowProj%gas_num, MaxNumGases)
        if (gas > lastGas) exit
        if (lEx%Flux0%gas(gas) /= error) &
            lEx%Flux0%gas(gas) = lEx%Flux0%gas(gas) / FluxnetGasScale(gas)
        if (lEx%rand_uncer(gas) /= error) &
            lEx%rand_uncer(gas) = lEx%rand_uncer(gas) / FluxnetGasScale(gas)
        if (lEx%r(gas) /= error) lEx%r(gas) = lEx%r(gas) / FluxnetGasScale(gas)
        if (lEx%chi(gas) /= error) lEx%chi(gas) = lEx%chi(gas) / FluxnetGasScale(gas)
    end do
    if (lEx%Ts /= error) lEx%Ts = lEx%Ts + 273.15d0
    if (lEx%Ta /= error) lEx%Ta = lEx%Ta + 273.15d0
    if (lEx%Tdew /= error) lEx%Tdew = lEx%Tdew + 273.15d0
    if (lEx%Pa /= error) lEx%Pa = lEx%Pa * 1d3
    if (lEx%Pd /= error) lEx%Pd = lEx%Pd * 1d3
    if (lEx%e /= error) lEx%e = lEx%e * 1d2
    if (lEx%es /= error) lEx%es = lEx%es * 1d2
    if (lEx%Pa /= error) lEx%VPD = lEx%VPD * 1d2
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

    !> Variances were actually read as standard deviations.
    !>
    !> Every configured gas: the read above fills stats%Cov over u..lastCfg,
    !> so squaring only as far as the fourth gas left every slot past it
    !> holding a standard deviation in a slot the rest of the code reads as a
    !> variance. It never surfaced because the copy below stopped at the same
    !> place and those columns came out as zero instead - the two bounds have
    !> to move together, or the zero becomes a plausible wrong number.
    do var = u, lastGas
        if (lEx%var_present(var)) &
            lEx%stats%Cov(var, var) = lEx%stats%Cov(var, var)**2
    end do

    lEx%instr(sonic)%category = 'sonic'

    !> Analyser units are converted at the read, per gas, and the
    !> self-describing block that follows it carries SI already - so there is
    !> nothing left to convert here. This used to be a loop over the
    !> instrument-role indices ico2..igas4 plus a mirror into gas_instr, whose
    !> only real job was to run *after* that block and undo it for the first
    !> four slots.

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
    !> Variances and w-covariances for every configured gas, not four.
    !>
    !> write_out_full_fcc already loops firstGas..lastGas over both, and the
    !> ex file already carries them - stats%Cov is read over u..lastCfg. Only
    !> this copy stopped at the fourth gas, so <gas>_var and w/<gas>_cov came
    !> out as exactly 0.00000 for every gas past it. That is a claim about the
    !> data, not a missing value: a variance of zero says the series was
    !> constant.
    !>
    !> Found by moving water from slot 6 to slot 9 between two fixtures and
    !> watching h2o_var and n2o_var trade places - the zero followed the slot,
    !> not the species.
    lEx%var(firstGas:lastGas) = error
    lEx%cov_w(firstGas:lastGas) = error
    do var = u, ts
        lEx%var(var) = lEx%stats%Cov(var, var)
    end do
    do var = firstGas, lastGas
        if (.not. lEx%var_present(var)) cycle
        lEx%var(var) = lEx%stats%Cov(var, var)
        lEx%cov_w(var) = lEx%stats%cov(w, var)
    end do
    lEx%cov_w(u) = lEx%stats%cov(w, u)
    lEx%cov_w(ts) = lEx%stats%cov(w, ts)
end subroutine CompleteEssentials
