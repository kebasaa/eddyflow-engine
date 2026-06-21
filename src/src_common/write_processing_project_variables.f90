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

    !> Initializations
    Auxfile%metadata   = 'none'
    Auxfile%biomet     = 'none'
    Dir%biomet         = 'none'
    Dir%main_out       = 'none'
    EddyFlowProj%fname_template = 'none'

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
    EddyFlowProj%col(co2) = nint(EPPrjNTags(4)%value)
    EddyFlowProj%col(h2o) = nint(EPPrjNTags(5)%value)
    EddyFlowProj%col(ch4) = nint(EPPrjNTags(6)%value)
    EddyFlowProj%col(gas4) = nint(EPPrjNTags(7)%value)
    EddyFlowProj%col(tc)  = nint(EPPrjNTags(8)%value)
    EddyFlowProj%col(ti1) = nint(EPPrjNTags(9)%value)
    EddyFlowProj%col(ti2) = nint(EPPrjNTags(10)%value)
    EddyFlowProj%col(pi)  = nint(EPPrjNTags(11)%value)
    EddyFlowProj%col(te)  = nint(EPPrjNTags(12)%value)
    EddyFlowProj%col(pe)  = nint(EPPrjNTags(13)%value)
    EddyFlowProj%col(E2NumVar + diag72) = nint(EPPrjNTags(14)%value)
    EddyFlowProj%col(E2NumVar + diag75) = nint(EPPrjNTags(15)%value)
    EddyFlowProj%col(E2NumVar + diag77) = nint(EPPrjNTags(16)%value)
    EddyFlowProj%col(E2NumVar + diagAnem) = nint(EPPrjNTags(20)%value)
    EddyFlowProj%col(E2NumVar + diagStaA) = nint(EPPrjNTags(21)%value)
    EddyFlowProj%col(E2NumVar + diagStaD) = nint(EPPrjNTags(22)%value)

    !> if a column was selected for gas4, read diffusivity. If diffusivity is
    !> below zero, defaults to gas4 diffusivity
    if (EddyFlowProj%col(gas4) > 0) then
        Dc(gas4) = EPPrjNTags(17)%value * 1d-4 !< takes from cm+2s-1 to m+2s-1
        if (Dc(gas4) <= 0) Dc(gas4) = 0.00001436d0  !< default for N2O from Massman (1998, J. Atm. Env) Table 2.
        MW(gas4) = sngl(EPPrjNTags(18)%value) * 1e-3 !< takes from g+1mol-1 to kg+1mol-1
        if (MW(gas4) <= 0) MW(gas4) = 44.01e-3  !< default for N2O
    end if

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
    !> select whether to use fixed or dynamic output format
    EddyFlowProj%fix_out_format = EPPrjCTags(37)%value(1:1) == '1'

    !> Select whether to apply high-pass theoretical spectral correction.
    !> It is independent from the choice of the low-pass method
    select case (EPPrjCTags(22)%value(1:1))
        case ('0')
            !> Do not apply low-frequency spectral correction
            EddyFlowProj%lf_meth = 'none'
        case ('1')
            EddyFlowProj%lf_meth = 'analytic'
    end select

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
        case default
        Meth%qcflag = 'mauder_foken_04'
    end select

    !> Select whether to standardize biomets or not
    EddyFlowProj%fluxnet_standardize_biomet = EPPrjCTags(48)%value(1:1) == '1'
    EddyFlowProj%fluxnet_mode = EPPrjCTags(49)%value(1:1) == '1'

    !> main output directory, only in Desktop mode
    if (EddyFlowProj%run_env /= 'embedded') then
        Dir%main_out = EPPrjCTags(35)%value
        if (len_trim(Dir%main_out) == 0) then
            write(*, *)
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
end subroutine WriteProcessingProjectVariables
