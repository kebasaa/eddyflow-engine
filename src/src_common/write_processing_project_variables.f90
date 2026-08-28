!***************************************************************************
! write_processing_project_variables.f90
! --------------------------------------
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
! \brief       Read EddyFlow configuration file, section [Project]
!              which is common to both RP and FCC
! \author      Gerardo Fratini
! \note
! \sa
! \bug
! \deprecated
! \test
! \todo
!***************************************************************************
subroutine WriteProcessingProjectVariables()
    use m_common_global_var
    implicit none
    !> local variables
    integer :: dot
    type(SwVerType) :: file_ini_ver
    include 'interfaces_1.inc'

    !> Refuse a project file written in a newer format than this engine knows.
    !> Without this the file is silently half-parsed: unknown keys are dropped
    !> and missing numeric tags are left undefined rather than defaulted.
    !> 'ini_version' is absent in very old files, which is fine - those predate
    !> the current format and are read on a best-effort basis, as before.
    if (EPPrjCTagFound(2)) then
        file_ini_ver = SwVerFromString(trim(adjustl(EPPrjCTags(2)%value)))
        if (.not. EqualSwVer(file_ini_ver, errSwVer)) then
            if (CompareSwVer(file_ini_ver, SwVerFromString(MaxSupportedIniVer))) then
                write(*, '(a)') '  Fatal error(96)> Project file format version ' &
                    // trim(adjustl(EPPrjCTags(2)%value)) &
                    // ' is newer than this engine supports (' &
                    // MaxSupportedIniVer // ').'
                write(ulog, '(a)') '  Fatal error(96)> Project file format version ' &
                    // trim(adjustl(EPPrjCTags(2)%value)) &
                    // ' is newer than this engine supports (' &
                    // MaxSupportedIniVer // ').'
                call ExceptionHandler(96)
            end if
        end if
    end if

    call ReadMeasurementRecords()

    !> Initializations
    Auxfile%metadata   = 'none'
    Auxfile%biomet     = 'none'
    Dir%biomet         = 'none'
    Dir%main_out       = 'none'
    EddyFlowProj%fname_template = 'none'
    EddyFlowProj%cec%h = 0d0
    EddyFlowProj%cec%singular_band = 0.2d0
    EddyFlowProj%cec%min_o1_o2 = 0.20d0
    EddyFlowProj%cec%min_octant = 0.05d0
    EddyFlowProj%cec%min_valid = 0.90d0
    EddyFlowProj%cec%signal_strength = 70d0
    EddyFlowProj%cec%max_stationarity = 25d0
    !> Zahn et al.'s criterion, so an untouched project reproduces the paper.
    EddyFlowProj%cec%stationarity_mode = cec_stat_flux
    !> Off, for the same reason: the paper applies no significance test.
    EddyFlowProj%cec%min_flux_sigma = 0d0
    EddyFlowProj%cec%max_gap_fill = 4

    !> Project general info
    select case (EPPrjCTags(16)%value(1:1))
        case('1')
            EddyFlowProj%run_mode =  'express'
        case('2')
            EddyFlowProj%run_mode =  'md_retrieval'
        case default
            EddyFlowProj%run_mode =  'advanced'
    end select

    EddyFlowProj%title  = trim(adjustl(EPPrjCTags(4)%value))
    EddyFlowProj%id     = trim(adjustl(EPPrjCTags(5)%value))
    if (EddyFlowProj%id(1:1) /= '_') then
        EddyFlowProj%id = 'eddyflow_' // trim(adjustl(EddyFlowProj%id))
    else
        EddyFlowProj%id = 'eddyflow' // trim(adjustl(EddyFlowProj%id))
    end if

    !>  file type
    select case (EPPrjCTags(6)%value(1:1))
        case ('0')
        EddyFlowProj%ftype = 'licor_ghg'
        EddyFlowProj%fext = 'ghg'
        case ('1')
        EddyFlowProj%ftype = 'generic_ascii'
        case ('2')
        EddyFlowProj%ftype = 'tob1'
        case ('3')
        EddyFlowProj%ftype = 'eddymeas_bin'
        case ('4')
        EddyFlowProj%ftype = 'edisol_bin'
        case ('5')
        EddyFlowProj%ftype = 'generic_bin'
        case ('6')
        EddyFlowProj%ftype = 'alteddy_bin'
    end select

    !> If file type is different from GHG, metadata
    !> retrieval mode is not feasible so forces into advanced mode
    if (EddyFlowProj%ftype /= 'licor_ghg' &
        .and. EddyFlowProj%run_mode ==  'md_retrieval') &
        EddyFlowProj%run_mode =  'advanced'

    !> File names prototype and related
    if (EddyFlowProj%run_env /= 'embedded') then
        EddyFlowProj%fname_template = &
            trim(adjustl(EPPrjCTags(7)%value))
        if (index(EddyFlowProj%fname_template, '.') /= 0) then
            !> File extensions
            dot = index(EddyFlowProj%fname_template, '.', .true.)
            EddyFlowProj%fext = &
            EddyFlowProj%fname_template(dot + 1:len_trim(EddyFlowProj%fname_template))
            !> ISO format
            EddyFlowLog%iso_format = index(EddyFlowProj%fname_template, 'mm') /= 0
        end if
    else
        EddyFlowLog%iso_format = .true.
        EddyFlowProj%fname_template = 'yyyy-mm-ddTHHMM'
    end if

    !> If file type is TOB1, check if user entered the format
    FileInterpreter%tob1_format = 'none'
    if (EddyFlowProj%ftype == 'tob1') then
        select case (trim(adjustl(EPPrjCTags(32)%value)))
            case ('1')
                FileInterpreter%tob1_format = 'IEEE4'
            case ('2')
                FileInterpreter%tob1_format = 'FP2'
        end select
    end if

    !> Whether to use alternative metadata file
    EddyFlowProj%use_extmd_file = EPPrjCTags(9)%value(1:1) == '1'
    AuxFile%metadata ='none'
    if(EddyFlowProj%use_extmd_file) &
        AuxFile%metadata = trim(adjustl(EPPrjCTags(10)%value))

    !> Whether to use dynamic metadata file
    EddyFlowProj%use_dynmd_file = EPPrjCTags(11)%value(1:1) == '1'
    if(EddyFlowProj%use_dynmd_file) &
        AuxFile%DynMD = trim(adjustl(EPPrjCTags(12)%value))

    !> Settings for binary raw files
    if (EddyFlowProj%ftype(1:len_trim(EddyFlowProj%ftype)) == 'generic_bin') then
        !> select line terminator in ASCII header of binary files
        select case(EPPrjCTags(13)%value(1:1))
            case ('0')
            Binary%ascii_head_eol = 'cr/lf'
            case ('1')
            Binary%ascii_head_eol = 'lf'
            case ('2')
            Binary%ascii_head_eol = 'cr'
        end select
        !> select binary files endianess
        Binary%little_endian = EPPrjCTags(14)%value(1:1) == '1'
        !> Select number of bytes per variable
        Binary%nbytes = nint(EPPrjNTags(1)%value)
        !> Select number of ASCII header lines
        Binary%head_nlines = nint(EPPrjNTags(2)%value)
    end if

    !> Master sonic
    EddyFlowProj%master_sonic = trim(adjustl(EPPrjCTags(15)%value))
    !> Variables to be used other than sonic ones
    EddyFlowProj%col(ts:pe) = nint(error)
    EddyFlowProj%col(ts)  = nint(EPPrjNTags(3)%value)
    !> The gas, cell and diagnostic columns are **not read from tags any
    !> more**. col_co2 .. col_diag_77 are retired: the records describe those
    !> measurements, and can say which analyser each came from and name the
    !> same species more than once, which one column per role never could.
    !>
    !> Their slots stay at nint(error), which is what DefineUsedVariables'
    !> "> 0" guard tests, so the legacy marking simply finds nothing and the
    !> record loops beside it do the work.
    !>
    !> Reading a retired tag would be worse than useless: the labels are
    !> blanked in the table, so SearchLocalTags never matches them and leaves
    !> %value untouched - the engine does not default missing tags.
    EddyFlowProj%col(te)  = nint(EPPrjNTags(12)%value)
    EddyFlowProj%col(pe)  = nint(EPPrjNTags(13)%value)
    !> Cleared explicitly: the "ts:pe" initialisation above does not reach the
    !> diagnostic slots, so dropping their assignments without this would
    !> leave them holding whatever was there - and DefineUsedVariables tests
    !> them with "> 0".
    EddyFlowProj%col(E2NumVar + diag72) = nint(error)
    EddyFlowProj%col(E2NumVar + diag75) = nint(error)
    EddyFlowProj%col(E2NumVar + diag77) = nint(error)
    EddyFlowProj%col(E2NumVar + diagAnem) = nint(error)
    !> Now that the slots are cleared, fill the ones the records name.
    call ApplyDiagnosticRecordColumns()
    EddyFlowProj%col(E2NumVar + diagStaA) = nint(EPPrjNTags(21)%value)
    EddyFlowProj%col(E2NumVar + diagStaD) = nint(EPPrjNTags(22)%value)

    !> The fourth gas's diffusivity and molecular weight used to be read here
    !> from gas_diff and gas_mw, gated on col_gas4. All three tags are retired,
    !> so the gate could never open - and a gas states its own overrides in its
    !> record now, applied per slot a few dozen lines below with a species
    !> default when it states none. That reaches every gas, not the fourth.

    !> Post-flux despiking (test_pfd), FCC-only, off by default.
    EddyFlowProj%test_pfd = EPPrjCTags(19)%value(1:1) == '1'

    !> biomet measurements info
    select case (EPPrjCTags(17)%value(1:1))
        case ('1')
        EddyFlowProj%biomet_data = 'embedded'
        case ('2')
        EddyFlowProj%biomet_data = 'ext_file'
        case ('3')
        EddyFlowProj%biomet_data = 'ext_dir'
        case default
        EddyFlowProj%biomet_data = 'none'
    end select
    !> biomet files/folders as applicable
    if (EddyFlowProj%biomet_data == 'ext_file') &
        AuxFile%biomet = trim(adjustl(EPPrjCTags(18)%value))
    if (EddyFlowProj%biomet_data == 'ext_dir') then
        Dir%biomet = trim(adjustl(EPPrjCTags(29)%value))
        if (len_trim(Dir%biomet) == 0) then
            Dir%biomet = 'none'
        else
            EddyFlowProj%biomet_tail = trim(adjustl(EPPrjCTags(30)%value))
            EddyFlowProj%biomet_recurse = EPPrjCTags(31)%value(1:1) == '1'
        end if
    end if

    !> If selected embedded biomet without GHG files (only possible via non-GUI
    !> file edit), set biomet to none.
    if (EddyFlowProj%biomet_data == 'embedded' &
        .and. EddyFlowProj%ftype /= 'licor_ghg') then
        call ExceptionHandler(93)
        EddyFlowProj%biomet_data = 'none'
    end if

    !> select whether to binned/full spectra files are available
    !> for current dataset
    EddyFlowProj%binned_spec_avail = EPPrjCTags(44)%value(1:1) == '1'
    EddyFlowProj%full_spec_avail   = EPPrjCTags(45)%value(1:1) == '1'

    !> select whether to output full output file
    EddyFlowProj%out_full = EPPrjCTags(21)%value(1:1) == '1'
    !> select whether to use fixed or dynamic output format
    EddyFlowProj%out_md = EPPrjCTags(39)%value(1:1) == '1'
    !> select whether to output average cospectra
    EddyFlowProj%out_avrg_cosp = EPPrjCTags(41)%value(1:1) == '1'
    !> select whether to output average spectra
    EddyFlowProj%out_avrg_spec = EPPrjCTags(43)%value(1:1) == '1'
    !> select whether to output biomet average values
    EddyFlowProj%out_biomet = EPPrjCTags(42)%value(1:1) == '1'

    !> Select whether to apply high-pass theoretical spectral correction.
    !> It is independent from the choice of the low-pass method
    select case (EPPrjCTags(22)%value(1:1))
        case ('0')
            !> Do not apply low-frequency spectral correction
            EddyFlowProj%lf_meth = 'none'
        case ('1')
            EddyFlowProj%lf_meth = 'analytic'
    end select

    !> Which analytic cospectrum the corrections are integrated against.
    !> A modifier on every analytic method rather than a method of its own -
    !> Moncrieff, Massman, Horst, Ibrom and Fratini all weight a transfer
    !> function by this shape, and all of them keep working whichever is
    !> chosen. Absent means the shape this program has always used, so an
    !> older project is unaffected.
    EddyFlowProj%cosp_model = 'moncrieff_97'
    if (EPPrjNTagFound(7)) then
        select case (nint(EPPrjNTags(7)%value))
            case (1)
                EddyFlowProj%cosp_model = 'kaimal_72'
            case (2)
                EddyFlowProj%cosp_model = 'sakai_01'
            case (3)
                EddyFlowProj%cosp_model = 'su_03'
            case (4)
                EddyFlowProj%cosp_model = 'moraes_08'
            case (5)
                EddyFlowProj%cosp_model = 'kristensen_97'
            case default
                !> Including 0, and including a value from some later version
                !> this one does not know. Falling back to the shape every
                !> correction was written against is the safe direction.
                EddyFlowProj%cosp_model = 'moncrieff_97'
        end select
    end if

    !> Iterative correction, after EddyUH.m:722-903.
    !>
    !> Off, and the defaults below are EddyUH's own, so switching it on
    !> reproduces what EddyUH does rather than something chosen here: four
    !> passes and no early exit. EddyUH's loop is `while indexITER <= 3` with
    !> the counter incremented at the top, so it runs four times; its only
    !> `break` tests the urban footprint's roughness length, requires
    !> indexITER > 3 - true on the last pass only - and sits in a branch that
    !> EddyUH_footprint.m:151 makes unreachable. Its covsvar output is a
    !> reported diagnostic, not a control.
    EddyFlowProj%corr_iter_meth = .false.
    EddyFlowProj%corr_iter_max = 4
    EddyFlowProj%corr_iter_tol = 0d0
    if (EPPrjNTagFound(8)) &
        EddyFlowProj%corr_iter_meth = nint(EPPrjNTags(8)%value) == 1
    if (EPPrjNTagFound(9)) &
        EddyFlowProj%corr_iter_max = nint(EPPrjNTags(9)%value)
    if (EPPrjNTagFound(10)) &
        EddyFlowProj%corr_iter_tol = EPPrjNTags(10)%value
    !> One pass is the un-iterated case and is what "off" already means, so a
    !> smaller number is a typed-in mistake rather than a setting. A negative
    !> tolerance would exit before the first comparison.
    if (EddyFlowProj%corr_iter_max < 1) EddyFlowProj%corr_iter_max = 4
    if (EddyFlowProj%corr_iter_tol < 0d0) EddyFlowProj%corr_iter_tol = 0d0

    !> Select low-pass spectral correction method.
    select case (EPPrjCTags(23)%value(1:1))
        case ('0')
            !> Do not apply spectral correction (e.g. open-path)
            EddyFlowProj%hf_meth = 'none'
            EddyFlowProj%hf_meth_in_situ = .false.
        case ('1')
            !> Correction after Moncrieff et al (1997, JH) fully analytical
            EddyFlowProj%hf_meth = 'moncrieff_97'
            EddyFlowProj%hf_meth_in_situ = .false.
        case ('2')
            !> Correction after Horst (1997, BLM), in-situ/analytical
            EddyFlowProj%hf_meth = 'horst_97'
            EddyFlowProj%hf_meth_in_situ = .true.
        case ('3')
            !> Correction after Ibrom et al (2007, AFM) fully in-situ
            EddyFlowProj%hf_meth = 'ibrom_07'
            EddyFlowProj%hf_meth_in_situ = .true.
        case ('4')
            !> Correction after Fratini et al. 2010, fully in-situ
            EddyFlowProj%hf_meth = 'fratini_12'
            EddyFlowProj%hf_meth_in_situ = .true.
        case ('5')
            !> Correction after Massman (2000, 2001), fully analytical
            EddyFlowProj%hf_meth = 'massman_00'
            EddyFlowProj%hf_meth_in_situ = .false.
        case ('6')
            !> Custom correction, in-situ/analytical
            EddyFlowProj%hf_meth = 'custom'
            EddyFlowProj%hf_meth_in_situ = .false.
        case default
            !> If not specified, set to none
            EddyFlowProj%hf_meth = 'none'
            EddyFlowProj%hf_meth_in_situ = .false.
    end select

    !> select whether to correct for LI-7550-related attenuations
    !> Relevant only for GHG files and logger software version < 7.7.0
    ! !>  Block-averaging
    ! EddyFlowProj%hf_correct_ghg_ba = EPPrjCTags(46)%value(1:1) == '1'
    ! !>  ZOH
    ! EddyFlowProj%hf_correct_ghg_zoh = EPPrjCTags(47)%value(1:1) == '1'
    ! if (EddyFlowProj%ftype /= 'licor_ghg') then
    !     EddyFlowProj%hf_correct_ghg_ba = .false.
    !     EddyFlowProj%hf_correct_ghg_zoh = .false.
    ! end if
    EddyFlowProj%hf_correct_ghg_ba = .false.
    EddyFlowProj%hf_correct_ghg_zoh = .false.

    EddyFlowProj%sonic_output_rate = nint(EPPrjNTags(19)%value)

    !> select whether to fill gaps with error codes
    EddyFlowProj%make_dataset = EPPrjCTags(24)%value(1:1) == '1'

    !> start/end date and time of period to be processed
    EddyFlowProj%subperiod = EPPrjCTags(40)%value(1:1) == '1'

    if (EddyFlowProj%subperiod) then
        EddyFlowProj%start_date = &
            trim(adjustl(EPPrjCTags(25)%value))
        EddyFlowProj%start_time = &
            trim(adjustl(EPPrjCTags(26)%value))
        EddyFlowProj%end_date = &
            trim(adjustl(EPPrjCTags(27)%value))
        EddyFlowProj%end_time = &
            trim(adjustl(EPPrjCTags(28)%value))
    end if

    if (len_trim(EddyFlowProj%start_date) == 0 &
        .or. len_trim(EddyFlowProj%start_time) == 0 &
        .or. len_trim(EddyFlowProj%end_date) == 0 &
        .or. len_trim(EddyFlowProj%end_time) == 0) &
        EddyFlowProj%subperiod = .false.

    !> select whether to apply WPL correction
    EddyFlowProj%wpl = EPPrjCTags(33)%value(1:1) /= '0'

    !> set error string
    EddyFlowProj%err_label = trim(adjustl(EPPrjCTags(36)%value))
    if (len_trim(EddyFlowProj%err_label) == 0 .or. EddyFlowProj%err_label == 'none') &
        EddyFlowProj%err_label = '-9999'

    !> select footprint method
    select case (EPPrjCTags(34)%value(1:1))
        case ('0')
        Meth%foot = 'none'
        case ('1')
        Meth%foot = 'kljun_04'
        case ('2')
        Meth%foot = 'kormann_meixner_01'
        case ('3')
        Meth%foot = 'hsieh_00'
        case default
        Meth%foot = 'kljun_04'
    end select

    !> select quality-flagging method
    select case (EPPrjCTags(38)%value(1:1))
        case ('0')
        Meth%qcflag = 'none'
        case ('1')
        Meth%qcflag = 'mauder_foken_04'
        case ('2')
        Meth%qcflag = 'foken_03'
        case ('3')
        Meth%qcflag = 'goeckede_06'
        case ('4')
        Meth%qcflag = 'vitale_20'
        case default
        Meth%qcflag = 'mauder_foken_04'
    end select

    !> Select whether to standardize biomets or not
    EddyFlowProj%fluxnet_standardize_biomet = EPPrjCTags(48)%value(1:1) == '1'
    EddyFlowProj%fluxnet_mode = EPPrjCTags(49)%value(1:1) == '1'

    !> Conditional Eddy Covariance (Zahn et al. 2022)
    !> The whole value, not its first character: reading one character put a
    !> ceiling of nine on a key that has no reason to have one.
    !>
    !> Two and three are accepted as well as one, and all three mean "on". The
    !> key was a three-way choice of which flux to partition before that choice
    !> moved onto the pairing, as cec_<i>_meth; a project written then still
    !> says 2 or 3, and refusing it would switch the partition off for those
    !> projects rather than telling anyone. The interface writes only 0 or 1.
    select case (trim(adjustl(EPPrjCTags(50)%value)))
        case ('1')
            EddyFlowProj%do_cec = 1
        case ('2')
            EddyFlowProj%do_cec = 2
        case ('3')
            EddyFlowProj%do_cec = 3
        case default
            EddyFlowProj%do_cec = 0
    end select
    if (EPPrjNTagFound(26)) &
        EddyFlowProj%cec%h = max(0d0, EPPrjNTags(26)%value)
    if (EPPrjNTagFound(4)) &
        EddyFlowProj%cec%singular_band = NormalizeCecBand( &
            EPPrjNTags(4)%value, 0.2d0)
    if (EPPrjNTagFound(27)) &
        EddyFlowProj%cec%min_o1_o2 = NormalizeCecPercent( &
            EPPrjNTags(27)%value, 0.20d0)
    if (EPPrjNTagFound(28)) &
        EddyFlowProj%cec%min_octant = NormalizeCecPercent( &
            EPPrjNTags(28)%value, 0.05d0)
    if (EPPrjNTagFound(29)) &
        EddyFlowProj%cec%min_valid = NormalizeCecPercent( &
            EPPrjNTags(29)%value, 0.90d0)
    if (EPPrjNTagFound(30)) &
        EddyFlowProj%cec%signal_strength = NormalizeCecSignalStrength( &
            EPPrjNTags(30)%value, 70d0)
    if (EPPrjNTagFound(31)) &
        EddyFlowProj%cec%max_gap_fill = NormalizeCecMaxGapFill( &
            EPPrjNTags(31)%value, 4)
    if (EPPrjNTagFound(32)) &
        EddyFlowProj%cec%max_stationarity = NormalizeCecStationarity( &
            EPPrjNTags(32)%value, 25d0)
    if (EPPrjNTagFound(6)) then
        !> Negative is meaningless, and reading it as "off" is the safe way to
        !> misread it.
        EddyFlowProj%cec%min_flux_sigma = max(0d0, EPPrjNTags(6)%value)
    end if
    if (EPPrjNTagFound(5)) then
        !> Anything but the ratio mode is the paper's, a value from some later
        !> version this one does not understand included. Falling back to the
        !> published criterion is the safe direction to be wrong in.
        if (nint(EPPrjNTags(5)%value) == cec_stat_ratio) then
            EddyFlowProj%cec%stationarity_mode = cec_stat_ratio
        else
            EddyFlowProj%cec%stationarity_mode = cec_stat_flux
        end if
    end if

    !> main output directory, only in Desktop mode
    if (EddyFlowProj%run_env /= 'embedded') then
        Dir%main_out = EPPrjCTags(35)%value
        if (len_trim(Dir%main_out) == 0) then
            write(*, *)
            write(ulog, *)
            call ExceptionHandler(36)
        end if
        call AdjDir(Dir%main_out, slash)
    end if

    !> Random error estimation settings
    select case (nint(EPPrjNTags(24)%value))
        case(1)
            RUsetup%meth = 'finkelstein_sims_01'
        case(2)
            RUsetup%meth = 'mann_lenschow_94'
        case(3)
            RUsetup%meth = 'mahrt_98'
        case(4)
            !> Billesbach (2011) random shuffle. Four rather than three
            !> because three is already Mahrt and the interface maps its
            !> menu onto these numbers; renumbering would silently change
            !> the method of every project that states one.
            RUsetup%meth = 'billesbach_11'
        case(5)
            !> Lenschow et al. (2000) instrumental noise, as Mauder et al.
            !> (2013) apply it. Distinct from mann_lenschow_94 above, which
            !> is a sampling error and a different paper - the shared name
            !> is the only thing they have in common.
            RUsetup%meth = 'lenschow_00'
        case default
            RUsetup%meth = 'none'
    end select
    if (RUsetup%meth /= 'none') then
        select case (nint(EPPrjNTags(23)%value))
            case(1)
                RUsetup%its_meth = 'cross_0'
            case(2)
                RUsetup%its_meth = 'full_integral'
            case default
                RUsetup%its_meth = 'cross_e'
        end select
        RUsetup%tlag_max = nint(EPPrjNTags(25)%value)
    end if

    !> Adjust paths
    call AdjFilePath(AuxFile%metadata, slash)
    call AdjFilePath(AuxFile%biomet, slash)
    call AdjDir(Dir%biomet, slash)
contains

!***************************************************************************
!> Read the indexed gas / cell / diagnostic records from the [Project] group.
!>
!> Slot positions come from the generated ProjectRecordOrigins parameters, so
!> this walks the groups by stride arithmetic exactly as ReadMetadataFile does
!> for instr_<K>_* and col_<N>_*, with no literal indices to rot when a slot is
!> appended.
!>
!> Every field is *TagFound-guarded: a missing numeric tag is left undefined by
!> SearchLocalTags rather than defaulted, so reading one blind yields garbage.
!>
!> These records ARE what drives processing. ApplyGasRecords fills E2Col from
!> them, SelectFluxnetGasSlots lays out the output from them, and the MW and
!> Dc tables below are filled per record. The flat col_co2/col_h2o/col_ch4/
!> col_gas4 keys they replaced are retired and blanked in the tag table.
!***************************************************************************
subroutine ReadMeasurementRecords()
    integer :: i
    integer :: b
    integer :: slot

    EddyFlowProj%gas_num  = 0
    EddyFlowProj%cell_num = 0
    EddyFlowProj%diag_num = 0
    EddyFlowProj%agc_num  = 0

    !> Whether the file states a gas count at all, which is a different
    !> question from whether that count is zero. See gas_num_stated.
    EddyFlowProj%gas_num_stated = EPPrjNTagFound(gasNumTag)

    if (EPPrjNTagFound(gasNumTag)) &
        EddyFlowProj%gas_num = nint(EPPrjNTags(gasNumTag)%value)
    if (EPPrjNTagFound(cellNumTag)) &
        EddyFlowProj%cell_num = nint(EPPrjNTags(cellNumTag)%value)
    if (EPPrjNTagFound(diagNumTag)) &
        EddyFlowProj%diag_num = nint(EPPrjNTags(diagNumTag)%value)
    !> An absent agc_num is not "no signal strength": it is a project written
    !> before these records existed, and CecSignalColumnFor falls back to
    !> matching a column named AGC or RSSI by name, which is all such a file
    !> has ever had.
    if (EPPrjNTagFound(agcNumTag)) &
        EddyFlowProj%agc_num = nint(EPPrjNTags(agcNumTag)%value)

    !> Clamp to what we can hold. The GUI enforces the same limits, so this
    !> only bites on a hand-edited file, but a silent overrun would be worse.
    EddyFlowProj%gas_num  = min(max(EddyFlowProj%gas_num,  0), MaxNumGases)
    EddyFlowProj%cell_num = min(max(EddyFlowProj%cell_num, 0), MaxNumCellCols)
    EddyFlowProj%diag_num = min(max(EddyFlowProj%diag_num, 0), MaxNumDiagCols)
    EddyFlowProj%agc_num  = min(max(EddyFlowProj%agc_num,  0), MaxNumAgcCols)

    do i = 1, MaxNumGases
        EddyFlowProj%gas(i) = GasRecordType('none', 'none', nint(error), 0, 0, &
                                            error, error, 0)
        if (i > EddyFlowProj%gas_num) cycle

        b = gasRecOriginC + (i - 1) * gasRecLeapC
        if (EPPrjCTagFound(b))     EddyFlowProj%gas(i)%var = &
            trim(adjustl(EPPrjCTags(b)%value))
        if (EPPrjCTagFound(b + 1)) EddyFlowProj%gas(i)%instr = &
            trim(adjustl(EPPrjCTags(b + 1)%value))

        b = gasRecOriginN + (i - 1) * gasRecLeapN
        if (EPPrjNTagFound(b))     EddyFlowProj%gas(i)%col   = nint(EPPrjNTags(b)%value)
        if (EPPrjNTagFound(b + 1)) EddyFlowProj%gas(i)%moist = nint(EPPrjNTags(b + 1)%value)
        if (EPPrjNTagFound(b + 2)) EddyFlowProj%gas(i)%cell  = nint(EPPrjNTags(b + 2)%value)
        if (EPPrjNTagFound(b + 3)) EddyFlowProj%gas(i)%mw    = dble(EPPrjNTags(b + 3)%value)
        if (EPPrjNTagFound(b + 4)) EddyFlowProj%gas(i)%diff  = dble(EPPrjNTags(b + 4)%value)
        if (EPPrjNTagFound(b + 5)) EddyFlowProj%gas(i)%fluxnet_default = &
            nint(EPPrjNTags(b + 5)%value)
    end do

    !> Carry each record's molecular weight and diffusivity onto its gas slot.
    !>
    !> MW and Dc are sized to E2NumVar but their `data` statements only fill
    !> co2:gas4, so every slot past the fourth gas holds whatever was in memory
    !> until something writes it. Nothing did: the records were read into
    !> EddyFlowProj%gas() and never applied, and the fourth slot got its values
    !> from the retired flat gas_mw/gas_diff tags alone. A gas with a garbage
    !> molecular weight produces a plausible-looking flux that is silently
    !> wrong, so this runs for every slot a record names.
    do i = 1, min(EddyFlowProj%gas_num, MaxNumGases)
        slot = firstGas + i - 1
        if (slot > lastGas) exit
        if (EddyFlowProj%gas(i)%col <= 0) cycle
        !> g mol-1 -> kg mol-1, cm+2 s-1 -> m+2 s-1, matching the units the
        !> interface writes and the defaults below.
        !>
        !> A record that carries neither still needs usable numbers, and those
        !> come from the species it names - not from the slot it happens to
        !> occupy. The fallback used to be gated on `slot > gas4`, which let
        !> the first four slots keep whatever the compile-time table put
        !> there: water declared at slot 9 was given N2O's molecular weight
        !> and diffusivity, and a gas declared at slot 6 was given water's.
        !> Only the two the interface can leave blank are wrong in practice,
        !> but the rule is the same for all of them.
        if (EddyFlowProj%gas(i)%mw > 0d0) then
            MW(slot) = sngl(EddyFlowProj%gas(i)%mw) * 1e-3
        else
            MW(slot) = DefaultMolecularWeight(EddyFlowProj%gas(i)%var)
        end if
        if (EddyFlowProj%gas(i)%diff > 0d0) then
            Dc(slot) = EddyFlowProj%gas(i)%diff * 1d-4
        else
            Dc(slot) = DefaultDiffusivity(EddyFlowProj%gas(i)%var)
        end if
        !> Say so when we had to guess. A record that carries neither value
        !> and names a species the tables do not know is given nitrous
        !> oxide's numbers, and a molecular weight wrong by a factor produces
        !> a flux wrong by the same factor while looking entirely ordinary.
        if (EddyFlowProj%gas(i)%mw <= 0d0 .and. EddyFlowProj%gas(i)%diff <= 0d0) then
            !> Nested rather than one .and. chain: gfortran warns that a
            !> function in a compound condition might not be evaluated.
            if (.not. HasSpeciesDefaults(EddyFlowProj%gas(i)%var)) &
                call ExceptionHandler(100)
        end if
    end do

    do i = 1, MaxNumCellCols
        EddyFlowProj%cell(i) = MeasRecordType('none', 'none', nint(error))
        if (i > EddyFlowProj%cell_num) cycle

        b = cellRecOriginC + (i - 1) * cellRecLeapC
        if (EPPrjCTagFound(b))     EddyFlowProj%cell(i)%var = &
            trim(adjustl(EPPrjCTags(b)%value))
        if (EPPrjCTagFound(b + 1)) EddyFlowProj%cell(i)%instr = &
            trim(adjustl(EPPrjCTags(b + 1)%value))

        b = cellRecOriginN + (i - 1) * cellRecLeapN
        if (EPPrjNTagFound(b)) EddyFlowProj%cell(i)%col = nint(EPPrjNTags(b)%value)
    end do

    do i = 1, MaxNumDiagCols
        EddyFlowProj%diag(i) = MeasRecordType('none', 'none', nint(error))
        if (i > EddyFlowProj%diag_num) cycle

        b = diagRecOriginC + (i - 1) * diagRecLeapC
        if (EPPrjCTagFound(b))     EddyFlowProj%diag(i)%var = &
            trim(adjustl(EPPrjCTags(b)%value))
        if (EPPrjCTagFound(b + 1)) EddyFlowProj%diag(i)%instr = &
            trim(adjustl(EPPrjCTags(b + 1)%value))

        b = diagRecOriginN + (i - 1) * diagRecLeapN
        if (EPPrjNTagFound(b)) EddyFlowProj%diag(i)%col = nint(EPPrjNTags(b)%value)
    end do

    !> Signal-strength columns.
    do i = 1, MaxNumAgcCols
        EddyFlowProj%agc(i) = MeasRecordType('none', 'none', nint(error))
        if (i > EddyFlowProj%agc_num) cycle

        b = agcRecOriginC + (i - 1) * agcRecLeapC
        if (EPPrjCTagFound(b))     EddyFlowProj%agc(i)%var = &
            trim(adjustl(EPPrjCTags(b)%value))
        if (EPPrjCTagFound(b + 1)) EddyFlowProj%agc(i)%instr = &
            trim(adjustl(EPPrjCTags(b + 1)%value))

        b = agcRecOriginN + (i - 1) * agcRecLeapN
        if (EPPrjNTagFound(b)) EddyFlowProj%agc(i)%col = nint(EPPrjNTags(b)%value)
    end do

    !> After the gas loop, because a pairing names gas records and there is no
    !> point resolving one against a list that is not read yet.
    call ReadCecRecords()

end subroutine ReadMeasurementRecords

!***************************************************************************
!
! \brief       The Conditional Eddy Covariance pairings the project states.
! \author      Jonathan Muller
! \note        A pairing is (one CO2 channel, one water channel, and any
!              further species to partition in the octants those two define).
!              Indices are into the gas records rather than raw columns, so a
!              re-ordered project keeps its pairings.
!
!              An absent cec_num leaves cec_num zero, and CecPairs then derives
!              one pairing per CO2 channel from the analyser layout. Stating
!              the pairings overrides that entirely - including stating none.
!***************************************************************************
subroutine ReadCecRecords()
    integer :: i
    integer :: k
    integer :: b
    integer :: ntok
    integer :: tok(MaxNumCecExtra)

    EddyFlowProj%cec_num = 0
    if (EPPrjNTagFound(cecNumTag)) &
        EddyFlowProj%cec_num = nint(EPPrjNTags(cecNumTag)%value)
    EddyFlowProj%cec_num = min(max(EddyFlowProj%cec_num, 0), MaxNumCecPairs)

    do i = 1, MaxNumCecPairs
        EddyFlowProj%cec_pair(i)%meth = 0
        EddyFlowProj%cec_pair(i)%carbon = 0
        EddyFlowProj%cec_pair(i)%water = 0
        EddyFlowProj%cec_pair(i)%extra = 0
        if (i > EddyFlowProj%cec_num) cycle

        !> A stated pairing is on unless it says otherwise, so a file that
        !> names a pairing and forgets the flag still gets a pairing.
        EddyFlowProj%cec_pair(i)%meth = 1

        b = cecRecOriginN + (i - 1) * cecRecLeapN
        if (EPPrjNTagFound(b))     EddyFlowProj%cec_pair(i)%meth = &
            nint(EPPrjNTags(b)%value)
        if (EPPrjNTagFound(b + 1)) EddyFlowProj%cec_pair(i)%carbon = &
            nint(EPPrjNTags(b + 1)%value)
        if (EPPrjNTagFound(b + 2)) EddyFlowProj%cec_pair(i)%water = &
            nint(EPPrjNTags(b + 2)%value)

        if (EddyFlowProj%cec_pair(i)%meth < 0 &
            .or. EddyFlowProj%cec_pair(i)%meth > 3) &
            EddyFlowProj%cec_pair(i)%meth = 1

        b = cecRecOriginC + (i - 1) * cecRecLeapC
        if (EPPrjCTagFound(b)) then
            call ParseCecExtraList(EPPrjCTags(b)%value, tok, ntok)
            do k = 1, ntok
                EddyFlowProj%cec_pair(i)%extra(k) = tok(k)
            end do
        end if
    end do
end subroutine ReadCecRecords

!***************************************************************************
!
! \brief       Split "6,7" into gas record indices.
! \author      Jonathan Muller
! \note        A list rather than a fixed run of numbered keys because most
!              pairings carry none and the few that do carry one or two, and
!              because the same shape already serves gas_<i>_sa_months.
!
!              Anything that is not a positive integer is dropped rather than
!              refused: an extra species that cannot be resolved costs the
!              project a column, not the run.
!***************************************************************************
subroutine ParseCecExtraList(text, tok, ntok)
    character(*), intent(in) :: text
    integer, intent(out) :: tok(MaxNumCecExtra)
    integer, intent(out) :: ntok

    character(len(text)) :: buf
    integer :: i
    integer :: value
    integer :: io_status

    tok = 0
    ntok = 0
    buf = text
    do i = 1, len_trim(buf)
        if (buf(i:i) == ',' .or. buf(i:i) == ';') buf(i:i) = ' '
    end do
    if (len_trim(buf) == 0) return

    do
        buf = adjustl(buf)
        if (len_trim(buf) == 0) exit
        read(buf, *, iostat = io_status) value
        if (io_status /= 0) exit
        if (value > 0 .and. ntok < MaxNumCecExtra) then
            ntok = ntok + 1
            tok(ntok) = value
        end if
        i = index(trim(buf), ' ')
        if (i <= 0) exit
        buf = buf(i:)
    end do
end subroutine ParseCecExtraList

!***************************************************************************
!
! \brief       Bridge the diagnostic records onto the internal column slots.
! \author      Jonathan Muller
! \note        Unlike the gas and cell records, which are applied directly to
!              E2Col, the diagnostics are consumed all over the engine through
!              EddyFlowProj%col(E2NumVar + diag*) - it is what sets
!              Diag7200%present, NumDiag and the per-analyser flag columns.
!              Marking the column "used" is not enough: without this the flags
!              are read but every INST_* output stays at its error value.
!              Populating the slot keeps all those consumers working
!              unchanged. Only the file *tags* are retired; this array is the
!              engine's own representation.
! \note        Must run after the diag slots are cleared to nint(error),
!              which happens well below the call that reads the records.
!***************************************************************************
subroutine ApplyDiagnosticRecordColumns()
    use m_common_global_var
    implicit none
    integer :: i

    do i = 1, min(EddyFlowProj%diag_num, MaxNumDiagCols)
        if (EddyFlowProj%diag(i)%col <= 0) cycle
        select case (trim(adjustl(EddyFlowProj%diag(i)%var)))
            case ('diag_72')
                EddyFlowProj%col(E2NumVar + diag72)   = EddyFlowProj%diag(i)%col
            case ('diag_75')
                EddyFlowProj%col(E2NumVar + diag75)   = EddyFlowProj%diag(i)%col
            case ('diag_77')
                EddyFlowProj%col(E2NumVar + diag77)   = EddyFlowProj%diag(i)%col
            case ('diag_anem')
                EddyFlowProj%col(E2NumVar + diagAnem) = EddyFlowProj%diag(i)%col
        end select
    end do
end subroutine ApplyDiagnosticRecordColumns

!> An occupancy limit is a percentage, and only a percentage. This used to
!> accept a fraction as well, deciding which was meant from the magnitude - so
!> the interface, whose spin boxes are labelled [%] and go down to 0.1, could
!> write 0.5 and have it read back as fifty percent.
real(kind = dbl) function NormalizeCecPercent(value, default_value)
    real(kind = dbl), intent(in) :: value
    real(kind = dbl), intent(in) :: default_value

    if (value >= 0d0 .and. value <= 100d0) then
        NormalizeCecPercent = value / 100d0
    else
        NormalizeCecPercent = default_value
    end if
end function NormalizeCecPercent

!> Half-width of the r_Fc ~ -1 band in which R and P nearly cancel and the
!> partition is a division by something near zero. Zahn et al. use 0.2; 0
!> switches the guard off.
real(kind = dbl) function NormalizeCecBand(value, default_value)
    real(kind = dbl), intent(in) :: value
    real(kind = dbl), intent(in) :: default_value

    if (value >= 0d0 .and. value <= 1d0) then
        NormalizeCecBand = value
    else
        NormalizeCecBand = default_value
    end if
end function NormalizeCecBand

real(kind = dbl) function NormalizeCecSignalStrength(value, default_value)
    real(kind = dbl), intent(in) :: value
    real(kind = dbl), intent(in) :: default_value

    if (value <= 0d0) then
        NormalizeCecSignalStrength = 0d0
    else if (value <= 100d0) then
        NormalizeCecSignalStrength = value
    else
        NormalizeCecSignalStrength = default_value
    end if
end function NormalizeCecSignalStrength

integer function NormalizeCecMaxGapFill(value, default_value)
    real(kind = dbl), intent(in) :: value
    integer, intent(in) :: default_value

    if (value >= 0d0) then
        NormalizeCecMaxGapFill = nint(value)
    else
        NormalizeCecMaxGapFill = default_value
    end if
end function NormalizeCecMaxGapFill

real(kind = dbl) function NormalizeCecStationarity(value, default_value)
    real(kind = dbl), intent(in) :: value
    real(kind = dbl), intent(in) :: default_value

    if (value >= 0d0) then
        NormalizeCecStationarity = value
    else
        NormalizeCecStationarity = default_value
    end if
end function NormalizeCecStationarity

!***************************************************************************
!> Molecular weight [kg mol-1] for a species a record names but does not
!> quantify. The interface writes an explicit mw for a gas it does not
!> recognise, so in practice these cover the species it does.
!>
!> The unrecognised default is N2O's, which is what the fourth slot has
!> always used, so nothing gets zero.
real(kind = sgl) function DefaultMolecularWeight(var)
    character(*), intent(in) :: var
    character(32) :: species

    species = var
    call uppercase(species)
    select case (trim(adjustl(species)))
        !> Carbon dioxide and nitrous oxide weigh nearly the same and were
        !> written to the same two decimals, so the table said they weigh
        !> *exactly* the same. On a site measuring both, that made the two
        !> indistinguishable by eye in the one place their identity is stated.
        !> Six figures on the atomic weights the rest of this table already
        !> uses - C 12.0107, N 14.0067, O 15.9994, H 1.00794 - separates them.
        case ('CO2'); DefaultMolecularWeight = 44.0095e-3
        case ('H2O'); DefaultMolecularWeight = MW_H2O
        case ('CH4'); DefaultMolecularWeight = 16.0425e-3
        case ('N2O'); DefaultMolecularWeight = 44.0128e-3
        case ('CO');  DefaultMolecularWeight = 28.0101e-3
        case ('SO2'); DefaultMolecularWeight = 64.066e-3
        case ('NH3'); DefaultMolecularWeight = 17.0305e-3
        case ('O3');  DefaultMolecularWeight = 47.9982e-3
        case ('NO2'); DefaultMolecularWeight = 46.0055e-3
        case ('NO');  DefaultMolecularWeight = 30.0061e-3
        case ('N2');  DefaultMolecularWeight = 28.0134e-3
        case ('O2');  DefaultMolecularWeight = 31.9988e-3
        case ('AR');  DefaultMolecularWeight = 39.948e-3
        case ('COS'); DefaultMolecularWeight = 60.075e-3
        case default; DefaultMolecularWeight = 44.0128e-3
    end select
end function DefaultMolecularWeight

!***************************************************************************
!> Whether the species tables above carry values for this gas.
!>
!> Callers use this to warn rather than to choose: a record naming a species
!> the engine does not know still gets numbers, but they are nitrous oxide's,
!> and a molecular weight that is wrong by a factor produces a flux that is
!> wrong by the same factor while looking entirely ordinary.
!***************************************************************************
logical function HasSpeciesDefaults(var)
    character(*), intent(in) :: var
    character(32) :: species

    species = var
    call uppercase(species)
    select case (trim(adjustl(species)))
        case ('CO2', 'H2O', 'CH4', 'N2O', 'CO', 'SO2', 'NH3', 'O3', &
              'NO2', 'NO', 'N2', 'O2', 'AR', 'COS')
            HasSpeciesDefaults = .true.
        case default
            HasSpeciesDefaults = .false.
    end select
end function HasSpeciesDefaults

!***************************************************************************
!> Molecular diffusivity in air [m+2 s-1], Massman (1998, Atm Env, Table 2),
!> for a species a record names but does not quantify. Same rule and same
!> unrecognised default as DefaultMolecularWeight.
real(kind = dbl) function DefaultDiffusivity(var)
    character(*), intent(in) :: var
    character(32) :: species

    species = var
    call uppercase(species)
    select case (trim(adjustl(species)))
        case ('CO2'); DefaultDiffusivity = 0.00001381d0
        case ('H2O'); DefaultDiffusivity = 0.00002178d0
        case ('CH4'); DefaultDiffusivity = 0.00001952d0
        case ('N2O'); DefaultDiffusivity = 0.00001436d0
        case ('CO');  DefaultDiffusivity = 0.00001807d0
        case ('SO2'); DefaultDiffusivity = 0.00001089d0
        case ('NH3'); DefaultDiffusivity = 0.00001978d0
        case ('O3');  DefaultDiffusivity = 0.00001444d0
        case ('NO2'); DefaultDiffusivity = 0.00001361d0
        case ('NO');  DefaultDiffusivity = 0.00001988d0
        case ('N2');  DefaultDiffusivity = 0.000019939d0
        case ('O2');  DefaultDiffusivity = 0.000020255d0
        case ('AR');  DefaultDiffusivity = 0.000019064d0
        case ('COS'); DefaultDiffusivity = 0.000012344d0
        case default; DefaultDiffusivity = 0.00001436d0
    end select
end function DefaultDiffusivity
end subroutine WriteProcessingProjectVariables
