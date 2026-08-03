!***************************************************************************
! read_ini_rp.f90
! ---------------
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
! \brief       Reads file "processing.eddypro", placed in the predifined program \n
!              folder "prog_folder\ini\" and stores relevant variables
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine ReadIniRP(key)
    use m_rp_global_var
    implicit none
    !> in/out variables
    character(*), intent(in) :: key
    !> local variables
    logical :: IniFileNotFound


    write(*,'(a)') ' Reading EddyFlow project file: ' &
                     // PrjPath(1:len_trim(PrjPath)) // '..'

    !> parse processing.eddypro file and store [Project] variables,
    !> common to all programs
    call ParseIniFile(PrjPath, 'Project', EPPrjNTags, EPPrjCTags,&
        size(EPPrjNTags), size(EPPrjCTags), EPPrjNTagFound, EPPrjCTagFound, &
        IniFileNotFound)

    if (IniFileNotFound) call ExceptionHandler(21)
    call WriteProcessingProjectVariables()

    !> parse processing.eddypro file and store all numeric and character tags
    call ParseIniFile(PrjPath, key, SNTags, SCTags, size(SNTags), size(SCTags),&
        SNTagFound, SCTagFound, IniFileNotFound)

    if (IniFileNotFound) call ExceptionHandler(21)
    !> selects only tags needed in this software, and store
    !> them in relevant variables
    call WriteVariablesRP()

    write(*,'(a)')   ' Done.'
end subroutine ReadIniRP

!*******************************************************************************
!
! \brief       Looks in "SNTags" and "SCTags" and retrieve variables used for \n
!              express processing.
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!*******************************************************************************
subroutine WriteVariablesRP()
    use m_rp_global_var
    implicit none
    !> local variables
    integer :: i
    integer :: j
    integer :: init_an_flags
    integer :: leap_an_flags
    integer :: init_an_wsect
    integer :: leap_an_wsect
    integer :: init_prof_z
    integer :: hlen
    integer :: gasslot
    character(64) :: raw_gas_tag(GHGNumVar)
    logical :: proceed

    !> Initializations
    Dir%main_in      = 'none'
    AuxFile%pf       = 'none'

    !> Flags for elimination of individual data points
    leap_an_flags = 3
    init_an_flags = 73 - leap_an_flags
    NumRawFlags = 0
    RPsetup%filter_by_raw_flags = .false.
    do i = 1, MaxNumRawFlags
        if (SNTagFound(init_an_flags + i*leap_an_flags) .and. &
            nint(SNTags(init_an_flags + i*leap_an_flags)%value) > 0) then
            NumRawFlags = NumRawFlags + 1
            RawFlag(NumRawFlags)%col = &
                nint(SNTags(init_an_flags + i*leap_an_flags)%value)
            RawFlag(NumRawFlags)%threshold = &
                dble(SNTags(init_an_flags + i*leap_an_flags + 1)%value)
            RawFlag(NumRawFlags)%upper = .true.
            if (SNTags(init_an_flags + i*leap_an_flags + 2)%value > 0) &
                RawFlag(NumRawFlags)%upper = .false.
        end if
    end do
    if (NumRawFlags > 0) RPsetup%filter_by_raw_flags = .true.

    !> In/out directories
    Dir%main_in = SCTags(1)%value(1:len_trim(SCTags(1)%value))
    if (len_trim(Dir%main_in) == 0) Dir%main_in = 'none'

    !> Wind direction filter option and corresponding sectors
    RPSetup%apply_wdf = SCTags(99)%value(1:1) == '1'
    if (RPSetup%apply_wdf) then
        leap_an_wsect = 2
        init_an_wsect = 373 - leap_an_wsect
        RPSetup%wdf_num_secs = 0
        do i = 1, MaxNumWdfSectors
            if (SNTagFound(init_an_wsect + i*leap_an_wsect)) then
                RPSetup%wdf_num_secs = RPSetup%wdf_num_secs + 1
                RPSetup%wdf_start(RPSetup%wdf_num_secs) = &
                    SNTags(init_an_wsect + i*leap_an_wsect)%value
                RPSetup%wdf_end(RPSetup%wdf_num_secs) = &
                    SNTags(init_an_wsect + i*leap_an_wsect + 1)%value
            end if
        end do
    end if

    !> Everything about raw statistical tests
    Test%sr = SCTags(3)%value(1:1) == '1'
    Test%ar = SCTags(4)%value(1:1) == '1'
    Test%do = SCTags(5)%value(1:1) == '1'
    Test%al = SCTags(6)%value(1:1) == '1'
    Test%sk = SCTags(7)%value(1:1) == '1'
    Test%ds = SCTags(8)%value(1:1) == '1'
    Test%tl = SCTags(9)%value(1:1) == '1'
    Test%aa = SCTags(10)%value(1:1) == '1'
    Test%ns = SCTags(11)%value(1:1) == '1'

    !> method of spike removal
    RPSetup%despike_vickers97 = SCTags(90)%value(1:1) == '0'

    !> Spike removal test
    sr%num_spk = idint(SNTags(1)%value)
    sr%lim_u = SNTags(2)%value
    sr%hf_lim = SNTags(3)%value
    ar%lim = idint(SNTags(4)%value)
    ar%bins = idint(SNTags(5)%value)
    ar%hf_lim = SNTags(6)%value
    sr%lim_w   = SNTags(54)%value
    !> Legacy per-gas spike limits fill the four legacy slots; per-gas records
    !> override them below when the project uses the record format.
    sr%lim_gas       = SNTags(55)%value
    sr%lim_gas(co2)  = SNTags(55)%value
    sr%lim_gas(h2o)  = SNTags(56)%value
    sr%lim_gas(ch4)  = SNTags(57)%value
    sr%lim_gas(gas4) = SNTags(58)%value

    !> Per-gas records override the legacy slots. Reading them here, rather
    !> than in place of the above, keeps a legacy project working untouched:
    !> with no records the loop does nothing.
    do i = 1, min(EddyFlowProj%gas_num, MaxNumGases)
        if (SNTagFound(rpGasOriginN + (i - 1) * rpGasLeapN)) &
            sr%lim_gas(firstGas + i - 1) = &
                SNTags(rpGasOriginN + (i - 1) * rpGasLeapN)%value
    end do

    !> Dropout test
    do%extlim_dw = SNTags(7)%value
    do%hf1_lim = SNTags(8)%value
    do%hf2_lim = SNTags(9)%value

    !> Absolute limits test
    al%u_max = SNTags(10)%value
    al%w_max = SNTags(11)%value
    al%t_min = SNTags(12)%value
    al%t_max = SNTags(13)%value
    al%gas_min(co2) = SNTags(14)%value
    al%gas_max(co2) = SNTags(15)%value
    al%gas_min(h2o) = SNTags(16)%value
    al%gas_max(h2o) = SNTags(17)%value
    al%gas_min(ch4) = SNTags(59)%value
    al%gas_max(ch4) = SNTags(60)%value
    al%gas_min(gas4) = SNTags(61)%value
    al%gas_max(gas4) = SNTags(62)%value

    !> Per-gas records override the legacy slots (offsets 1 and 2 from the
    !> record origin: sr_lim, al_min, al_max, ...).
    do i = 1, min(EddyFlowProj%gas_num, MaxNumGases)
        if (SNTagFound(rpGasOriginN + (i - 1) * rpGasLeapN + 1)) &
            al%gas_min(firstGas + i - 1) = &
                SNTags(rpGasOriginN + (i - 1) * rpGasLeapN + 1)%value
        if (SNTagFound(rpGasOriginN + (i - 1) * rpGasLeapN + 2)) &
            al%gas_max(firstGas + i - 1) = &
                SNTags(rpGasOriginN + (i - 1) * rpGasLeapN + 2)%value
    end do

    !> Skewness and Kurtosis
    sk%hf_skmin = SNTags(18)%value
    sk%hf_skmax = SNTags(19)%value
    sk%sf_skmin = SNTags(20)%value
    sk%sf_skmax = SNTags(21)%value
    sk%hf_kumin = SNTags(22)%value
    sk%hf_kumax = SNTags(23)%value
    sk%sf_kumin = SNTags(24)%value
    sk%sf_kumax = SNTags(25)%value

    !> Discontinuities
    ds%hf_uv = SNTags(26)%value
    ds%hf_w = SNTags(27)%value
    ds%hf_t = SNTags(28)%value
    ds%hf_gas(co2) = SNTags(29)%value
    ds%hf_gas(h2o) = SNTags(30)%value
    ds%hf_gas(ch4) = SNTags(63)%value
    ds%hf_gas(gas4) = SNTags(64)%value
    ds%hf_var = SNTags(31)%value
    ds%sf_uv = SNTags(32)%value
    ds%sf_w = SNTags(33)%value
    ds%sf_t = SNTags(34)%value
    ds%sf_gas(co2) = SNTags(35)%value
    ds%sf_gas(h2o) = SNTags(36)%value
    ds%sf_gas(ch4) = SNTags(65)%value
    ds%sf_gas(gas4) = SNTags(66)%value

    !> Per-gas records override the legacy slots (offsets 3 and 4: sr_lim,
    !> al_min, al_max, ds_hf, ds_sf, ...).
    do i = 1, min(EddyFlowProj%gas_num, MaxNumGases)
        if (SNTagFound(rpGasOriginN + (i - 1) * rpGasLeapN + 3)) &
            ds%hf_gas(firstGas + i - 1) = &
                SNTags(rpGasOriginN + (i - 1) * rpGasLeapN + 3)%value
        if (SNTagFound(rpGasOriginN + (i - 1) * rpGasLeapN + 4)) &
            ds%sf_gas(firstGas + i - 1) = &
                SNTags(rpGasOriginN + (i - 1) * rpGasLeapN + 4)%value
    end do
    ds%sf_var = SNTags(37)%value

    !> Timelag
    tl%hf_lim = SNTags(38)%value
    tl%sf_lim = SNTags(39)%value
    tl%def_gas(co2) = SNTags(40)%value
    tl%def_gas(h2o) = SNTags(41)%value
    tl%def_gas(ch4) = SNTags(67)%value
    tl%def_gas(gas4) = SNTags(68)%value

    !> Per-gas records override the legacy slots (offset 5: sr_lim, al_min,
    !> al_max, ds_hf, ds_sf, tl_def, ...).
    do i = 1, min(EddyFlowProj%gas_num, MaxNumGases)
        if (SNTagFound(rpGasOriginN + (i - 1) * rpGasLeapN + 5)) &
            tl%def_gas(firstGas + i - 1) = &
                SNTags(rpGasOriginN + (i - 1) * rpGasLeapN + 5)%value
    end do

    !> Angle of attack
    aa%min = SNTags(42)%value
    aa%max = SNTags(43)%value
    aa%lim = SNTags(44)%value

    !> Non-statiorarity of horizontal wind
    ns%hf_lim = SNTags(45)%value

    !> select angle-of-attack calibration option
    if (nint(SNTags(290)%value) < 0) then
        RPsetup%calib_aoa = 'automatic'
    else
        select case (nint(SNTags(290)%value))
            case (1)
                RPSetup%calib_aoa = 'nakai_12'
            case (2)
                RPSetup%calib_aoa = 'nakai_06'
            case default
                RPSetup%calib_aoa = 'none'
        end select
    end if
    !> Select whether to apply the w-boost correction to WM/WMPro sonics
    RPsetup%calib_wboost = SCTags(67)%value(1:1) == '1'

    !> Cross-wind correction
    RPsetup%calib_cw = SCTags(13)%value(1:1) == '1'
    !> select whether to look for raw files in sub-folders
    RPsetup%recurse = SCTags(19)%value(1:1) == '1'
    !> select whether to output binned (co)spectra
    RPsetup%out_bin_sp = SCTags(26)%value(1:1) == '1'

    !> select whether to output binned ogives
    RPsetup%out_bin_og = SCTags(51)%value(1:1) == '1'

    !> Regardless of user selection, if ensemble averaged (co)spectra
    !> have been requested and binned spectra files are not available,
    !> need to create them.
    if ( (EddyFlowProj%out_avrg_cosp .or. EddyFlowProj%out_avrg_spec) &
        .and. .not. EddyFlowProj%binned_spec_avail) RPsetup%out_bin_sp = .true.

    !> select output file
    !> Clear both arrays first: only the slots below are assigned explicitly,
    !> and out_full_cosp(w_w) never is, yet callers such as InitOutFilesRP and
    !> SpectralAnalysis scan the whole 1..GHGNumVar range. Without this the
    !> unassigned slots are read undefined.
    RPsetup%out_full_sp   = .false.
    RPsetup%out_full_cosp = .false.

    RPsetup%out_full_sp(u)   = SCTags(27)%value(1:1) == '1'
    RPsetup%out_full_sp(v)   = SCTags(28)%value(1:1) == '1'
    RPsetup%out_full_sp(w)   = SCTags(29)%value(1:1) == '1'
    RPsetup%out_full_sp(ts)  = SCTags(30)%value(1:1) == '1'
    RPsetup%out_full_sp(co2) = SCTags(31)%value(1:1) == '1'
    RPsetup%out_full_sp(h2o) = SCTags(32)%value(1:1) == '1'
    RPsetup%out_full_sp(ch4) = SCTags(33)%value(1:1) == '1'
    RPsetup%out_full_sp(gas4) = SCTags(34)%value(1:1) == '1'

    RPsetup%out_full_cosp(w_u)   = SCTags(35)%value(1:1) == '1'
    RPsetup%out_full_cosp(w_v)   = SCTags(36)%value(1:1) == '1'
    RPsetup%out_full_cosp(w_ts)  = SCTags(37)%value(1:1) == '1'
    RPsetup%out_full_cosp(w_co2) = SCTags(38)%value(1:1) == '1'
    RPsetup%out_full_cosp(w_h2o) = SCTags(39)%value(1:1) == '1'
    RPsetup%out_full_cosp(w_ch4) = SCTags(40)%value(1:1) == '1'
    RPsetup%out_full_cosp(w_gas4) = SCTags(41)%value(1:1) == '1'

    RPsetup%out_st(1) = SCTags(42)%value(1:1) == '1'
    RPsetup%out_st(2) = SCTags(43)%value(1:1) == '1'
    RPsetup%out_st(3) = SCTags(44)%value(1:1) == '1'
    RPsetup%out_st(4) = SCTags(45)%value(1:1) == '1'
    RPsetup%out_st(5) = SCTags(46)%value(1:1) == '1'
    RPsetup%out_st(6) = SCTags(47)%value(1:1) == '1'
    RPsetup%out_st(7) = SCTags(48)%value(1:1) == '1'

    RPsetup%out_raw(1) = SCTags(68)%value(1:1) == '1'
    RPsetup%out_raw(2) = SCTags(69)%value(1:1) == '1'
    RPsetup%out_raw(3) = SCTags(70)%value(1:1) == '1'
    RPsetup%out_raw(4) = SCTags(71)%value(1:1) == '1'
    RPsetup%out_raw(5) = SCTags(72)%value(1:1) == '1'
    RPsetup%out_raw(6) = SCTags(73)%value(1:1) == '1'
    RPsetup%out_raw(7) = SCTags(74)%value(1:1) == '1'

    RPsetup%out_raw_var(u)   = SCTags(75)%value(1:1) == '1'
    RPsetup%out_raw_var(v)   = SCTags(76)%value(1:1) == '1'
    RPsetup%out_raw_var(w)   = SCTags(77)%value(1:1) == '1'
    RPsetup%out_raw_var(ts)  = SCTags(78)%value(1:1) == '1'
    RPsetup%out_raw_var(co2) = SCTags(79)%value(1:1) == '1'
    RPsetup%out_raw_var(h2o) = SCTags(80)%value(1:1) == '1'
    RPsetup%out_raw_var(ch4) = SCTags(81)%value(1:1) == '1'
    RPsetup%out_raw_var(gas4) = SCTags(82)%value(1:1) == '1'
    RPsetup%out_raw_var(te)  = SCTags(83)%value(1:1) == '1'
    RPsetup%out_raw_var(pe)  = SCTags(84)%value(1:1) == '1'

    !> Per-gas output selections. They override the four legacy flags above
    !> for the gases those can reach, and are the only way to reach a fifth.
    !> SCTagFound-guarded, so a project without records is untouched.
    !>
    !> All three arrays are indexed by variable slot: w_u..w_gas4 are aliases
    !> of u..gas4 rather than a separate numbering, and SpectralAnalysis pairs
    !> out_full_cosp(var) with Stats%cov(w, var). So gas record i, which
    !> ApplyGasRecords places in slot firstGas + i - 1, uses that index in all
    !> three. Note out_raw_var and not out_raw: the latter is seven processing
    !> stages, not a per-variable selection.
    do i = 1, MaxNumGases
        gasslot = firstGas + i - 1
        if (gasslot > lastGas) exit

        if (SCTagFound(rpGasOriginC + (i - 1) * rpGasLeapC)) &
            RPsetup%out_full_sp(gasslot) = &
                SCTags(rpGasOriginC + (i - 1) * rpGasLeapC)%value(1:1) == '1'

        if (SCTagFound(rpGasOriginC + (i - 1) * rpGasLeapC + 1)) &
            RPsetup%out_full_cosp(gasslot) = &
                SCTags(rpGasOriginC + (i - 1) * rpGasLeapC + 1)%value(1:1) == '1'

        if (SCTagFound(rpGasOriginC + (i - 1) * rpGasLeapC + 2)) &
            RPsetup%out_raw_var(gasslot) = &
                SCTags(rpGasOriginC + (i - 1) * rpGasLeapC + 2)%value(1:1) == '1'
    end do

    !> If no spectral output is selected, identify this situation for skipping
    !> completely the spectral analysis.
    !> Bounded by lastGas rather than gas4: a project whose fifth gas asked
    !> for spectra would otherwise have the whole analysis skipped.
    RPsetup%do_spectral_analysis = .false.
    if (RPsetup%out_bin_sp .or. RPsetup%out_bin_og &
        .or. any(RPsetup%out_full_sp(u:lastGas)) &
        .or. any(RPsetup%out_full_cosp(w_u:w_v)) &
        .or. any(RPsetup%out_full_cosp(w_ts:lastGas))) &
        RPsetup%do_spectral_analysis = .true.

    !> If no variable was selected for output, force out_raw to false
    !> regardless of user setting
    if (.not. any(RPsetup%out_raw_var(u:pe))) RPsetup%out_raw = .false.

    !> Raw dataset dir
    proceed = .false.
    do i = 1, 7
        if (RPsetup%out_raw(i)) then
            proceed = .true.
            exit
        end if
    end do

    !> Define header of raw dataset files
    raw_out_header = '   '
    hlen = 3
    if (RPsetup%out_raw_var(u))   then
        raw_out_header = raw_out_header(1:hlen) // 'u'
        hlen = hlen + 25
    end if
    if (RPsetup%out_raw_var(v))   then
        raw_out_header = raw_out_header(1:hlen) // 'v'
        hlen = hlen + 25
    end if
    if (RPsetup%out_raw_var(w))   then
        raw_out_header = raw_out_header(1:hlen) // 'w'
        hlen = hlen + 25
    end if
    if (RPsetup%out_raw_var(ts))  then
        raw_out_header = raw_out_header(1:hlen) // 'ts'
        hlen = hlen + 25
    end if
    !> One header name per selected gas slot, over the same range OutRawData
    !> writes columns for (it loops out_raw_var over every column). Enumerating
    !> only the historical four here would give a fifth gas a data column with
    !> no name, leaving header and rows a field apart.
    !>
    !> Names come from the project configuration rather than from E2Col, whose
    !> %var is still empty at this point: DefineE2Set has not run yet. A
    !> project with no gas records keeps the historical names exactly.
    call FullOutputGasTags(raw_gas_tag)
    do i = 1, MaxNumGases
        gasslot = firstGas + i - 1
        if (gasslot > lastGas) exit
        if (.not. RPsetup%out_raw_var(gasslot)) cycle

        !> Disambiguated, so two records of the same species do not give two
        !> columns of the same name. This used to take the record's %var
        !> verbatim, and a site measuring CO2 on two analysers wrote `co2`
        !> twice - a header a reader cannot key on. FullOutputGasTags is
        !> record-derived for exactly this reason and is safe before
        !> DefineE2Set has run; its stems carry a trailing underscore.
        if (len_trim(raw_gas_tag(gasslot)) > 1) then
            raw_out_header = raw_out_header(1:hlen) &
                // raw_gas_tag(gasslot)(1:len_trim(raw_gas_tag(gasslot)) - 1)
        else if (EddyFlowProj%gas_num > 0 .and. i <= EddyFlowProj%gas_num &
            .and. len_trim(EddyFlowProj%gas(i)%var) > 0) then
            raw_out_header = raw_out_header(1:hlen) &
                // trim(EddyFlowProj%gas(i)%var)
        else
            select case (gasslot)
                case (co2);  raw_out_header = raw_out_header(1:hlen) // 'co2'
                case (h2o);  raw_out_header = raw_out_header(1:hlen) // 'h2o'
                case (ch4);  raw_out_header = raw_out_header(1:hlen) // 'ch4'
                case (gas4); raw_out_header = raw_out_header(1:hlen) // '4th gas'
                case default
                    raw_out_header = raw_out_header(1:hlen) // 'gas'
            end select
        end if
        hlen = hlen + 25
    end do
    if (RPsetup%out_raw_var(te))  then
        raw_out_header = raw_out_header(1:hlen) // 'air_t'
        hlen = hlen + 25
    end if
    if (RPsetup%out_raw_var(pe))  then
        raw_out_header = raw_out_header(1:hlen) // 'air_p'
        hlen = hlen + 25
    end if

    !> Output QC details
    RPsetup%out_qc_details = SCTags(85)%value(1:1) == '1'

    !> select the averaging length. If zero, files are processed as they are
    RPsetup%avrg_len = nint(SNTags(50)%value)

    !> select detrending method
    select case (SCTags(14)%value(1:1))
        case ('0')
        Meth%det = 'ba'
        case ('1')
        Meth%det = 'ld'
        case ('2')
        Meth%det = 'rm'
        case ('3')
        Meth%det = 'ew'
        case default
        Meth%det = 'ba'
    end select
    if (Meth%det == 'ld' .or. Meth%det == 'rm' .or. Meth%det == 'ew') &
        RPsetup%Tconst = nint(SNTags(46)%value)

    !> select rotation method
    select case (SCTags(15)%value(1:1))
        case ('0')
        Meth%rot = 'none'
        case ('1')
        Meth%rot = 'double_rotation'
        case ('2')
        Meth%rot = 'triple_rotation'
        case ('3')
        Meth%rot = 'planar_fit'
        case ('4')
        Meth%rot = 'planar_fit_no_bias'
        case default
        Meth%rot = 'none'
    end select

    !> Planar fit extra settings
    RPsetup%pf_onthefly = .false.
    RPsetup%pf_assessment_only = SCTagFound(100) .and. &
        SCTags(100)%value(1:1) == '1'
    if (index(Meth%rot, 'planar_fit') /= 0) then
        !> Whether to perfom planar fit on the fly or use previous results file
        if (SCTags(56)%value(1:1) == '1') then
            RPsetup%pf_onthefly = .true.
        else
            AuxFile%pf = SCTags(57)%value(1:len_trim(SCTags(57)%value))
        end if
        !> Whether to subtract b0 from mean w
        RPsetup%pf_subtract_b0 = SCTags(96)%value(1:1) /= '1'
    end if
    !> Assessment-only planar fit applies only to an on-the-fly planar fit.
    RPsetup%pf_assessment_only = RPsetup%pf_assessment_only .and. &
        RPsetup%pf_onthefly

    !> select time lag handling method
    select case (SCTags(16)%value(1:1))
        case ('0')
        Meth%tlag = 'none'
        case ('1')
        Meth%tlag = 'constant'
        case ('2')
        Meth%tlag = 'maxcov&default'
        case ('3')
        Meth%tlag = 'maxcov'
        case ('4')
        Meth%tlag = 'tlag_opt'
        case ('5')
        Meth%tlag = 'pwb'
        case default
        Meth%tlag = 'none'
    end select

    !> Pre-whitening block-bootstrap (Vitale et al. 2024) defaults.
    PWBSetup%min_lag(co2) = -10d0
    PWBSetup%max_lag(co2) =  10d0
    PWBSetup%min_lag(h2o) = -10d0
    PWBSetup%max_lag(h2o) =  10d0
    PWBSetup%min_lag(ch4) = -10d0
    PWBSetup%max_lag(ch4) =  10d0
    PWBSetup%min_lag(gas4) = -10d0
    PWBSetup%max_lag(gas4) =  10d0
    PWBSetup%lag_bounds_provided = .false.
    PWBSetup%n_bootstrap = 99
    PWBSetup%block_length_s = 20d0
    PWBSetup%min_valid_frac = 0.3d0
    PWBSetup%hdi_thresh_s = 0.5d0
    PWBSetup%dev_thresh_s = 0.5d0
    PWBSetup%hdi_prefilter_s = 1.0d0
    PWBSetup%smoothing_width = 5
    PWBSetup%random_seed = 2024
    if (SNTagFound(406)) PWBSetup%min_lag(co2) = SNTags(406)%value
    if (SNTagFound(407)) PWBSetup%max_lag(co2) = SNTags(407)%value
    if (SNTagFound(406) .or. SNTagFound(407)) PWBSetup%lag_bounds_provided(co2) = .true.
    if (SNTagFound(408)) PWBSetup%min_lag(h2o) = SNTags(408)%value
    if (SNTagFound(409)) PWBSetup%max_lag(h2o) = SNTags(409)%value
    if (SNTagFound(408) .or. SNTagFound(409)) PWBSetup%lag_bounds_provided(h2o) = .true.
    if (SNTagFound(410)) PWBSetup%min_lag(ch4) = SNTags(410)%value
    if (SNTagFound(411)) PWBSetup%max_lag(ch4) = SNTags(411)%value
    if (SNTagFound(410) .or. SNTagFound(411)) PWBSetup%lag_bounds_provided(ch4) = .true.
    if (SNTagFound(412)) PWBSetup%min_lag(gas4) = SNTags(412)%value
    if (SNTagFound(413)) PWBSetup%max_lag(gas4) = SNTags(413)%value
    if (SNTagFound(412) .or. SNTagFound(413)) PWBSetup%lag_bounds_provided(gas4) = .true.
    !> Per-gas records override the legacy slots (offsets 9 and 10:
    !> pwb_min_lag, pwb_max_lag). lag_bounds_provided is set the same way the
    !> legacy pairs above set it, so a record-supplied window counts as
    !> user-provided rather than falling back to the default range.
    do i = 1, min(EddyFlowProj%gas_num, MaxNumGases)
        if (firstGas + i - 1 > lastGas) exit
        if (SNTagFound(rpGasOriginN + (i - 1) * rpGasLeapN + 9)) then
            PWBSetup%min_lag(firstGas + i - 1) = &
                SNTags(rpGasOriginN + (i - 1) * rpGasLeapN + 9)%value
            PWBSetup%lag_bounds_provided(firstGas + i - 1) = .true.
        end if
        if (SNTagFound(rpGasOriginN + (i - 1) * rpGasLeapN + 10)) then
            PWBSetup%max_lag(firstGas + i - 1) = &
                SNTags(rpGasOriginN + (i - 1) * rpGasLeapN + 10)%value
            PWBSetup%lag_bounds_provided(firstGas + i - 1) = .true.
        end if
    end do

    if (SNTagFound(414)) PWBSetup%n_bootstrap = max(1, nint(SNTags(414)%value))
    if (SNTagFound(415)) PWBSetup%block_length_s = SNTags(415)%value
    if (SNTagFound(416)) PWBSetup%min_valid_frac = SNTags(416)%value
    if (SNTagFound(417)) PWBSetup%hdi_thresh_s = SNTags(417)%value
    if (SNTagFound(418)) PWBSetup%dev_thresh_s = SNTags(418)%value
    if (SNTagFound(419)) PWBSetup%hdi_prefilter_s = SNTags(419)%value
    if (SNTagFound(420)) PWBSetup%smoothing_width = max(1, nint(SNTags(420)%value))
    if (SNTagFound(421)) PWBSetup%random_seed = max(1, nint(SNTags(421)%value))
    PWBSetup%approx_ccf   = .false.
    PWBSetup%max_ar_order = 0
    if (SNTagFound(422)) PWBSetup%approx_ccf   = nint(SNTags(422)%value) /= 0
    if (SNTagFound(423)) PWBSetup%max_ar_order = max(0, nint(SNTags(423)%value))
    PWBSetup%detect_prewpl = .false.
    if (SNTagFound(424)) PWBSetup%detect_prewpl = nint(SNTags(424)%value) /= 0

    !> Time lag optimizer extra settings
    RPsetup%to_onthefly = .false.
    PwbCacheGenerate = .false.
    PwbCacheLoaded = .false.
    PwbCacheDirty = .false.
    PwbCacheUpdateRequested = .false.
    RPsetup%tlag_assessment_only = SCTagFound(101) .and. &
        SCTags(101)%value(1:1) == '1'
    TimeLagOptSelected = .false.
    if (Meth%tlag == 'tlag_opt') then
        TimeLagOptSelected = .true.
        if (SCTags(91)%value(1:1) == '1') then
            RPsetup%to_onthefly = .true.
        else
            AuxFile%to = ''
            if (len_trim(SCTags(92)%value) > 0) &
                AuxFile%to = SCTags(92)%value(1:len_trim(SCTags(92)%value))
        end if
    elseif (Meth%tlag == 'pwb') then
        if (SCTags(91)%value(1:1) == '1') then
            !> PWB on-the-fly means generate a per-period PWB cache before
            !> production processing, not an aggregate time-lag optimization.
            RPsetup%to_onthefly = .true.
            PwbCacheGenerate = .true.
        else
            AuxFile%to = ''
            if (len_trim(SCTags(92)%value) > 0) &
                AuxFile%to = SCTags(92)%value(1:len_trim(SCTags(92)%value))
            !> An empty selected-file field is live PWB: detect during the
            !> production pass and persist the resulting cache afterwards.
            PwbCacheUpdateRequested = len_trim(AuxFile%to) > 0
        end if
        !> PWB aggregate summaries use the existing H2O RH-class layout even
        !> when no aggregate optimizer prepass was requested.
        TOSetup%h2o_nclass = max(1, nint(SNTags(207)%value))
        TOSetup%h2o_class_size = floor(100d0 / TOSetup%h2o_nclass)
        TOSetup%pg_range = SNTags(198)%value
    end if
    !> Assessment-only time-lag optimization applies only to an on-the-fly optimizer.
    RPsetup%tlag_assessment_only = RPsetup%tlag_assessment_only .and. &
        RPsetup%to_onthefly

    !>  tapering window
    select case (SCTags(17)%value(1:1))
        case ('0')
        RPsetup%tap_win = 'squared'
        case ('1')
        RPsetup%tap_win = 'bartlett'
        case ('2')
        RPsetup%tap_win = 'welch'
        case ('3')
        RPsetup%tap_win = 'hamming'
        case ('4')
        RPsetup%tap_win = 'hann'
    end select

    !> number of frequency bins
    Meth%spec%nbins = nint(SNTags(48)%value)

    RPsetup%tcell_filter_tconst = nint(SNTags(372)%value)

    !> max acceptable lack of data lines in a raw file
    RPsetup%max_lack = SNTags(49)%value

    !> read wind speed offsets
    RPsetup%offset(u) = 0d0
    RPsetup%offset(v) = 0d0
    RPsetup%offset(w) = 0d0
    RPsetup%offset(u) = dble(SNTags(51)%value)
    RPsetup%offset(v) = dble(SNTags(52)%value)
    RPsetup%offset(w) = dble(SNTags(53)%value)

    !> Planar fit settings
    if (index(Meth%rot, 'planar_fit') /= 0) then
        if (RPsetup%pf_onthefly) then
            PFSetup%subperiod     = SCTags(97)%value(1:1) == '1'
            PFSetup%start_date    = SCTags(49)%value(1:len_trim(SCTags(49)%value))
            PFSetup%end_date      = SCTags(50)%value(1:len_trim(SCTags(50)%value))
            PFSetup%start_time    = SCTags(22)%value(1:len_trim(SCTags(22)%value))
            PFSetup%end_time      = SCTags(23)%value(1:len_trim(SCTags(23)%value))
            PFSetup%min_per_sec   = nint(SNTags(70)%value)
            PFSetup%w_max         = SNTags(71)%value
            PFSetup%u_min         = SNTags(72)%value
            !> If w_max is found to be < 0.099, it means it has not been set, so it
            !> is forced to the max value, which implies no filtering for w_max.
            if(PFSetup%w_max  <= 0.099d0) PFSetup%w_max = 10d0

            !> Customization of wind sectors
            leap_an_wsect = 2
            init_an_wsect = 209 - leap_an_wsect
            PFSetup%num_sec = 0
            PFSetup%north_offset = SNTags(208)%value
            do i = 1, MaxNumWSect
                if (SNTagFound(init_an_wsect + i*leap_an_wsect) .and. &
                    SNTags(init_an_wsect + i*leap_an_wsect)%value > 0) then
                    PFSetup%num_sec = PFSetup%num_sec + 1
                    PFSetup%width(PFSetup%num_sec) = &
                        SNTags(init_an_wsect + i*leap_an_wsect)%value
                    PFSetup%wsect_exclude(PFSetup%num_sec) = &
                        nint(SNTags(init_an_wsect + i*leap_an_wsect + 1)%value) == 1

                end if
            end do
            if (PFSetup%num_sec == 0) then
                call ExceptionHandler(40)
                PFSetup%num_sec = 1
            elseif (PFSetup%num_sec == 1) then
                PFSetup%wsect_end(PFSetup%num_sec) = 360
            elseif (PFSetup%num_sec > 1) then
                !> Calculate ending angle of each sector
                do i = 1, PFSetup%num_sec
                    PFSetup%wsect_end(i) = nint(sum(PFSetup%width(1:i)))
                    if (PFSetup%wsect_end(i) < 0) &
                        PFSetup%wsect_end(i) = 360 + PFSetup%wsect_end(i)
                end do
                PFSetup%wsect_end(PFSetup%num_sec) = 360
            end if
        end if
        PFSetup%fix = 'clockwise'
        if(SCTags(88)%value(1:1) == '1') PFSetup%fix = 'counterclockwise'
        if(SCTags(88)%value(1:1) == '2') PFSetup%fix = 'double_rotation'
    end if

    !> Time lag optimizer settings
    !>
    !> The period is set unconditionally. WriteOutTimelagOptimization prints
    !> it into every summary it writes, including the PWB aggregate one, which
    !> runs with to_onthefly false - so these four fields were being written
    !> straight out of uninitialised memory. The result was NUL bytes in the
    !> two header lines, which is why the file reads as binary rather than as
    !> text. The processing period is the right answer for the PWB summary
    !> anyway: that is the span it aggregates over.
    TOSetup%start_date = EddyFlowProj%start_date
    TOSetup%end_date   = EddyFlowProj%end_date
    TOSetup%start_time = EddyFlowProj%start_time
    TOSetup%end_time   = EddyFlowProj%end_time

    if (RPsetup%to_onthefly) then
        TOSetup%h2o_nclass    = 0
        TOSetup%subperiod     = SCTags(98)%value(1:1) == '1'
        TOSetup%start_date    = SCTags(93)%value(1:len_trim(SCTags(93)%value))
        TOSetup%end_date      = SCTags(94)%value(1:len_trim(SCTags(94)%value))
        TOSetup%start_time    = SCTags(24)%value(1:len_trim(SCTags(24)%value))
        TOSetup%end_time      = SCTags(25)%value(1:len_trim(SCTags(25)%value))
        TOSetup%gas_min_flux(co2)  = SNTags(194)%value
        TOSetup%gas_min_flux(ch4)  = SNTags(195)%value
        TOSetup%gas_min_flux(gas4) = SNTags(196)%value
        TOSetup%le_min_flux   = SNTags(197)%value

        !> Per-gas records override the legacy slots (offset 6: sr_lim,
        !> al_min, al_max, ds_hf, ds_sf, tl_def, to_min_flux, ...).
        do i = 1, min(EddyFlowProj%gas_num, MaxNumGases)
            if (SNTagFound(rpGasOriginN + (i - 1) * rpGasLeapN + 6)) &
                TOSetup%gas_min_flux(firstGas + i - 1) = &
                    SNTags(rpGasOriginN + (i - 1) * rpGasLeapN + 6)%value
        end do
        TOSetup%pg_range      = SNTags(198)%value
        TOSetup%min_lag(co2)   = SNTags(199)%value
        TOSetup%max_lag(co2)   = SNTags(200)%value
        TOSetup%min_lag(h2o)   = SNTags(201)%value
        TOSetup%max_lag(h2o)   = SNTags(202)%value
        TOSetup%min_lag(ch4)   = SNTags(203)%value
        TOSetup%max_lag(ch4)   = SNTags(204)%value
        TOSetup%min_lag(gas4)  = SNTags(205)%value
        TOSetup%max_lag(gas4)  = SNTags(206)%value
        TOSetup%h2o_nclass    = nint(SNTags(207)%value)
        if (TOSetup%h2o_nclass > 1) then
            TOSetup%h2o_class_size = floor(100d0 / TOSetup%h2o_nclass)
        end if

        !> Per-gas records override the legacy slots (offsets 7 and 8:
        !> to_min_lag, to_max_lag). Guarded on the tag being present, so a
        !> project without records keeps exactly the values read above.
        do i = 1, min(EddyFlowProj%gas_num, MaxNumGases)
            if (firstGas + i - 1 > lastGas) exit
            if (SNTagFound(rpGasOriginN + (i - 1) * rpGasLeapN + 7)) &
                TOSetup%min_lag(firstGas + i - 1) = &
                    SNTags(rpGasOriginN + (i - 1) * rpGasLeapN + 7)%value
            if (SNTagFound(rpGasOriginN + (i - 1) * rpGasLeapN + 8)) &
                TOSetup%max_lag(firstGas + i - 1) = &
                    SNTags(rpGasOriginN + (i - 1) * rpGasLeapN + 8)%value
        end do
    end if

    !> Timelag by covariance maximization options
    select case (SCTags(59)%value(1:len_trim(SCTags(59)%value)))
    case('w')
        RPSetup%covmax_var = w
    case('ts')
        RPSetup%covmax_var = ts
    case default
        RPSetup%covmax_var = w
    end select
    RPSetup%covmax_stocdet = SCTags(60)%value(1:1) == '1'

    !> Biomet measurements
    select case (SCTags(61)%value(1:len_trim(SCTags(61)%value)))
        case('comma')
            bFileMetadata%separator = ','
        case('semicolon')
            bFileMetadata%separator = ';'
        case('space')
            bFileMetadata%separator = ' '
        case('tab')
            bFileMetadata%separator = char(9)
        case default
            bFileMetadata%separator = SCTags(61)%value(1:1)
    end select
    bFileMetadata%tstamp_ref = SCTags(62)%value(1: len_trim(SCTags(62)%value))
    bFileMetadata%nhead = nint(SNTags(192)%value)

    !> Wheter to filter for spikes and abolute limits
    RPsetup%filter_sr = SCTags(63)%value(1:1) == '1'
    RPsetup%filter_al = SCTags(64)%value(1:1) == '1'

    !> Burba correction params
    RPsetup%bu_corr = 'none'
    if(SCTags(65)%value == '1')  RPsetup%bu_corr = 'yes'
    if(SCTags(65)%value == '-1') RPsetup%bu_corr = 'none'
    RPsetup%bu_multi = SCTags(66)%value(1:1) == '1'

    !> Whether to use power-of-two samples for FFT
    RPsetup%power_of_two = SCTags(87)%value(1:1) == '1'

    !> Whether to refer wind direction to geographic north
    RPsetup%use_geo_north = SCTags(89)%value(1:1) == '1'

    magnetic_declination = 0d0
    if (RPsetup%use_geo_north) &
        magnetic_declination = nint(SNTags(193)%value)

    !> Biomet measurements numeric params
    bSetup%sel(bTa)   = nint(SNTags(111)%value)
    bSetup%sel(bPa)   = nint(SNTags(112)%value)
    bSetup%sel(bRH)   = nint(SNTags(113)%value)
    bSetup%sel(bPPFD) = nint(SNTags(114)%value)
    bSetup%sel(bLWin) = nint(SNTags(115)%value)
    bSetup%sel(bRg)   = nint(SNTags(116)%value)

    init_prof_z = 120
    bSetup%zT    = error
    bSetup%zCO2  = error
    bSetup%zH2O  = error
    bSetup%zCH4  = error
    bSetup%zGAS4 = error
    do i = 1, MaxProfNodes
        if (nint(SNTags(init_prof_z + i)%value) >= 0) &
            bSetup%zT(i) = nint(SNTags(init_prof_z + i)%value)
        if (nint(SNTags(init_prof_z + MaxProfNodes + i)%value) >= 0) &
            bSetup%zCO2(i) = nint(SNTags(init_prof_z + MaxProfNodes + i)%value)
        if (nint(SNTags(init_prof_z + 2 * MaxProfNodes + i)%value) >= 0) &
            bSetup%zH2O(i)  = nint(SNTags(init_prof_z + 2 * MaxProfNodes + i)%value)
        if (nint(SNTags(init_prof_z + 3 * MaxProfNodes + i)%value) >= 0) &
            bSetup%zCH4(i)  = nint(SNTags(init_prof_z + 3 * MaxProfNodes + i)%value)
        if (nint(SNTags(init_prof_z + 4 * MaxProfNodes + i)%value) >= 0) &
            bSetup%zGAS4(i) = nint(SNTags(init_prof_z + 4 * MaxProfNodes + i)%value)
    end do
    do i = 1, MaxProfNodes - 1
        if (bSetup%zT(i + 1) /= error .and. bSetup%zT(i) /= error) then
            bSetup%dz(1, i) = bSetup%zT(i + 1) - bSetup%zT(i)
        else
            bSetup%dz(1, i) = error
        end if
        bSetup%dz(2, i) = bSetup%zCO2(i + 1)  - bSetup%zCO2(i)
        bSetup%dz(3, i) = bSetup%zH2O(i + 1)  - bSetup%zH2O(i)
        bSetup%dz(4, i) = bSetup%zCH4(i + 1)  - bSetup%zCH4(i)
        bSetup%dz(5, i) = bSetup%zGAS4(i + 1) - bSetup%zGAS4(i)
    end do

    !> Parameters for Burba correction
    !> Multiple linear regressions
    BurbaPar%m(daytime, bot, 1)     =  SNTags(156)%value
    BurbaPar%m(daytime, bot, 2)     =  SNTags(157)%value
    BurbaPar%m(daytime, bot, 3)     =  SNTags(158)%value
    BurbaPar%m(daytime, bot, 4)     =  SNTags(159)%value
    BurbaPar%m(daytime, top, 1)     =  SNTags(160)%value
    BurbaPar%m(daytime, top, 2)     =  SNTags(161)%value
    BurbaPar%m(daytime, top, 3)     =  SNTags(162)%value
    BurbaPar%m(daytime, top, 4)     =  SNTags(163)%value
    BurbaPar%m(daytime, spar, 1)    =  SNTags(164)%value
    BurbaPar%m(daytime, spar, 2)    =  SNTags(165)%value
    BurbaPar%m(daytime, spar, 3)    =  SNTags(166)%value
    BurbaPar%m(daytime, spar, 4)    =  SNTags(167)%value
    BurbaPar%m(nighttime, bot, 1)   =  SNTags(168)%value
    BurbaPar%m(nighttime, bot, 2)   =  SNTags(169)%value
    BurbaPar%m(nighttime, bot, 3)   =  SNTags(170)%value
    BurbaPar%m(nighttime, bot, 4)   =  SNTags(171)%value
    BurbaPar%m(nighttime, top, 1)   =  SNTags(172)%value
    BurbaPar%m(nighttime, top, 2)   =  SNTags(173)%value
    BurbaPar%m(nighttime, top, 3)   =  SNTags(174)%value
    BurbaPar%m(nighttime, top, 4)   =  SNTags(175)%value
    BurbaPar%m(nighttime, spar, 1)  =  SNTags(176)%value
    BurbaPar%m(nighttime, spar, 2)  =  SNTags(177)%value
    BurbaPar%m(nighttime, spar, 3)  =  SNTags(178)%value
    BurbaPar%m(nighttime, spar, 4)  =  SNTags(179)%value
    !> Simple linear regressions
    BurbaPar%l(daytime, bot, 1)     =  SNTags(180)%value
    BurbaPar%l(daytime, bot, 2)     =  SNTags(181)%value
    BurbaPar%l(daytime, top, 1)     =  SNTags(182)%value
    BurbaPar%l(daytime, top, 2)     =  SNTags(183)%value
    BurbaPar%l(daytime, spar, 1)    =  SNTags(184)%value
    BurbaPar%l(daytime, spar, 2)    =  SNTags(185)%value
    BurbaPar%l(nighttime, bot, 1)   =  SNTags(186)%value
    BurbaPar%l(nighttime, bot, 2)   =  SNTags(187)%value
    BurbaPar%l(nighttime, top, 1)   =  SNTags(188)%value
    BurbaPar%l(nighttime, top, 2)   =  SNTags(189)%value
    BurbaPar%l(nighttime, spar, 1)  =  SNTags(190)%value
    BurbaPar%l(nighttime, spar, 2)  =  SNTags(191)%value

    !> Settings related to drift correction
    !> initializations
    DriftCorr%method = 'none'
    select case (nint(SNTags(300)%value))
        case(1)
            DriftCorr%method = 'linear'
        case(2)
            DriftCorr%method = 'signal_strength'
        case default
            DriftCorr%method = 'none'
    end select
    DriftCorr%dir_cal = error
    DriftCorr%inv_cal = error
    DriftCorr%b = error
    DriftCorr%c = error
    !> read values
    !>
    !> The tag table carries seven direct and seven inverse coefficients for
    !> each of the four legacy slots - drift_dir_co2_0..6 at 301, h2o at 308,
    !> ch4 at 315, gas4 at 322, and the inverse set from 329 - but only the
    !> first two of each were read. The CH4 and fourth-gas polynomials could
    !> be written into a project and were then silently discarded, so those
    !> gases were never drift-corrected however they were configured.
    !>
    !> Guarded, unlike the four slices this replaces. SearchLocalTags leaves a
    !> tag it did not find untouched, so an unguarded read copied whatever the
    !> table happened to hold over the `error` set two lines up - and `error`
    !> is what marks a gas as having no calibration polynomial and thus no
    !> drift correction.
    do i = co2, gas4
        do j = 0, 6
            if (SNTagFound(301 + (i - co2) * 7 + j)) &
                DriftCorr%dir_cal(j, i) = SNTags(301 + (i - co2) * 7 + j)%value
            if (SNTagFound(329 + (i - co2) * 7 + j)) &
                DriftCorr%inv_cal(j, i) = SNTags(329 + (i - co2) * 7 + j)%value
        end do
    end do
    DriftCorr%b = SNTags(370)%value
    DriftCorr%c = SNTags(371)%value

    !> adjust paths
    call AdjDir(Dir%main_in, slash)
    call AdjFilePath(AuxFile%pf, slash)
    call AdjFilePath(AuxFile%to, slash)
end subroutine WriteVariablesRP
