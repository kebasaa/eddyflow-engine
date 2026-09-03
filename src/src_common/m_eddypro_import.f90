!***************************************************************************
! m_eddypro_import.f90
! --------------------
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
! \brief       Run an EddyPro project, by importing it into EddyFlow's own
!              format once, beside the file it came from.
! \author      Jonathan Muller
! \note        An EddyPro project cannot be read as it stands. Two things stop
!              it: InitEnv discards a positional path that is not *.eddyflow,
!              and from format version 5.0.0 a project stating no gas_num is
!              refused outright - Fatal error(99) - which is the shape of
!              every EddyPro file there is.
!
!              The gap is otherwise small. Of the keys in an EddyPro project,
!              all but about ninety are read verbatim by this engine under the
!              same name. Those ninety are one structural difference: EddyPro
!              describes gases, cell measurements and diagnostics as four
!              fixed column slots, with per-species settings named after the
!              slot (sr_lim_co2, al_n2o_max, sa_fmax_gas4). EddyFlow describes
!              them as indexed records that name the analyser as well as the
!              column, and can carry the same species twice.
!
!              So this reads the EddyPro pair, writes an ordinary EddyFlow
!              pair beside it, and hands the engine the path to that. Nothing
!              downstream knows an import happened.
!
!              The source files are never written. The import runs only when
!              the imported project does not yet exist, because RP edits the
!              project it ran - EditIniFile on ex_file - and FCC reads that
!              edit back. Re-importing on the FCC run would discard it and
!              send FCC to the EddyPro project's stale ex_file.
!
!              The mapping is not invented here. It is the one the interface
!              already performs in EcProject::importEddyProProject and its
!              neighbours; the tables below are that code's, and
!              static_checks/test_eddypro_import_static.py holds the two in
!              step across the repositories.
! \sa          init_env.f90, write_processing_project_variables.f90
! \bug
! \deprecated
! \test        static_checks/test_eddypro_import_static.py
! \todo
!***************************************************************************
module m_eddypro_import
    use m_common_global_var
    implicit none
    private
    public :: ImportEddyProProject

    !> First lines. An EddyPro project declares itself; so does ours, and the
    !> interface refuses to open a project whose first line says neither.
    character(*), parameter :: EddyProTag  = ';EDDYPRO_PROCESSING'
    character(*), parameter :: EddyFlowTag = ';EDDYFLOW_PROCESSING'

    !> What the imported pair is called. Deliberately not <base>.metadata:
    !> that is the name the EddyPro metadata already carries, in the same
    !> folder, so writing there would destroy the file being imported.
    character(*), parameter :: ImportSuffix = '_ep_imported'

    !> The record format this import writes. It must not exceed
    !> MaxSupportedIniVer, or the file we have just written is refused by
    !> Fatal error(96) the moment it is read back.
    character(*), parameter :: ImportedIniVer = '5.1.0'

    !> EddyPro has four gas slots and no more.
    integer, parameter :: nEddyProGasSlots = 4

    !> How EddyPro spells those slots. The fourth is spelled two ways within
    !> one file - n2o in the parameter settings and in out_full_sp_* /
    !> out_full_cosp_w_*, gas4 in the spectral, time-lag and out_raw_*
    !> families - so every lookup tries both rather than carrying one table
    !> per section. A file hand-edited into the other spelling then still
    !> reads, which one table per section would not manage.
    character(4), parameter :: slotName(nEddyProGasSlots) = &
        ['co2 ', 'h2o ', 'ch4 ', 'n2o ']
    character(4), parameter :: slotAlt(nEddyProGasSlots) = &
        ['co2 ', 'h2o ', 'ch4 ', 'gas4']

    !> The metadata's own first line, and the two sections the added keys
    !> belong to. DlProject::saveProject writes this tag unconditionally and
    !> never ';EDDYFLOW_METADATA', which it only accepts on read.
    character(*), parameter :: GhgMetadataTag = ';GHG_METADATA'
    character(*), parameter :: mdSectInstruments = 'Instruments'
    character(*), parameter :: mdSectColumns = 'FileDescription'

    !> What EddyFlow added to the metadata format, and nothing else: diffed
    !> against the interface's EddyPro v6.2.2 base, these three families are
    !> the whole delta. Each value is the one the engine already applies when
    !> the key is absent, so stating it changes nothing - see the comments at
    !> each writer below for why.
    character(*), parameter :: acFreqDefault = '0.000'
    character(*), parameter :: integratesDefault = '0'
    character(*), parameter :: errorValueDefault = '-9999.0000'

    !> Sections the injected keys go into. ReadIniRP sweeps RawProcess* and
    !> ReadIniFCC sweeps FluxCorrection*, each as one scope, so a per-gas key
    !> written into the wrong family is swept by neither and is silently
    !> absent - which reads downstream as "the project stated nothing".
    character(*), parameter :: sectProject = 'Project'
    character(*), parameter :: sectSettings = 'RawProcess_Settings'
    character(*), parameter :: sectParams = 'RawProcess_ParameterSettings'
    character(*), parameter :: sectTimelag = 'RawProcess_TimelagOptimization_Settings'
    character(*), parameter :: sectSpectral = 'FluxCorrection_SpectralAnalysis_General'
    character(*), parameter :: sectTilt = 'RawProcess_TiltCorrection_Settings'
    character(*), parameter :: sectPwb = 'RawProcess_PWBTimelag_Settings'
    character(*), parameter :: sectBiomet = 'RawProcess_BiometMeasurements'

    !***********************************************************************
    !> What the import states that the source did not, in two kinds.
    !>
    !> KIND ONE, and all but one of these: the settings EddyFlow added to the
    !> project format, which no EddyPro file can state, each with THE VALUE
    !> THIS ENGINE ALREADY APPLIES WHEN THE KEY IS ABSENT. Stating them is what
    !> makes the imported file a complete EddyFlow project rather than one the
    !> interface has to migrate on open; taking the engine's own value is what
    !> makes stating them a no-op. Every one is *TagFound-guarded at its
    !> reader, so absent and stated-at-this-value take the same branch.
    !>
    !> KIND TWO, and tlag_meth alone: an EddyPro key restated because the
    !> engine's absent-key behaviour for it is WRONG. This is the one entry
    !> here that deliberately changes what a run does, and it is last in the
    !> tables so it does not read as one of the others.
    !>
    !> An EddyPro project almost always states tlag_meth, and when it does it
    !> is copied through untouched like any other key - the import takes over
    !> whatever method the source chose, and cannot produce the fifth value,
    !> PWB, which EddyPro has no way to express. But if the key is absent the
    !> reader has no *TagFound guard: read_ini_rp.f90 tests
    !> SCTags(16)%value(1:1), finds a blank, takes case default and runs
    !> Meth%tlag = 'none' - no time-lag compensation at all, silently. Two is
    !> covariance maximization with default, which is EddyPro's own default
    !> and the interface's, so the converted project runs what it displays.
    !>
    !> The reader is deliberately not changed with it: a hand-written
    !> .eddyflow that omits the key still gets 'none', and widening this to
    !> every native project is a far larger question than the import path.
    !>
    !> Keep this in step with the readers, which own these numbers:
    !>   cec_*                    write_processing_project_variables.f90:75-81
    !>   cosp_model               write_processing_project_variables.f90, the
    !>                            unconditional 'moncrieff_97' above the tag
    !>   automatic_spectra_config read_ini_fcc.f90:334
    !>   rot_pf_assessment_only   read_ini_rp.f90:500
    !>   tlag_assessment_only     read_ini_rp.f90:589
    !>   pwb_*                    read_ini_rp.f90:542-549
    !> static_checks/test_eddypro_import_static.py reads them back out of
    !> those files, so a default that moves there fails the check here.
    !>
    !> Grouped by section, and the grouping is load-bearing: a section the
    !> source never opened is emitted at the end, one header per run of
    !> entries. tlag_meth is the sole RawProcess_Settings entry and sits on
    !> its own at the end for that reason as well as the one above.
    !>
    !> flux_run_mode and biom_rh_override are in no tag table here and are
    !> written for the interface alone, so that a converted project opens
    !> there with nothing missing.
    !>
    !> pwb_approx_ccf, pwb_max_ar_order and pwb_detect_prewpl stood here and
    !> are retired: the first two were speed options that cost accuracy for
    !> under a percent of runtime, and the third chose a detection stage that
    !> is no longer a choice. A converted project simply does not carry them,
    !> which is what absent has always meant here.
    !>
    !> Deliberately NOT here, because for these "absent" is itself a decision
    !> and no value reproduces it:
    !>   gas_<i>_pwb_min_lag/_max_lag  presence sets lag_bounds_provided, and
    !>                                 set_timelags.f90:147 takes the window
    !>                                 from instrument geometry without it
    !>   pf_sect_<k>_*, wdf_sect_<k>_* presence CREATES a sector;
    !>                                 read_ini_rp.f90:687 counts them
    !>   gas_<i>_drift_dir_*/_inv_*    absent means the gas is not drift-
    !>                                 corrected; there is no default curve
    !>   instr_<K>_max_lack            absent means "follow max_lack", which
    !>                                 pinning would stop tracking
    !>   gas_<i>_fluxnet_default       absent means the lowest record of the
    !>                                 species is designated
    !> pf_sect_* and wdf_sect_* are EddyPro keys in any case, present since
    !> the fork, and are copied through whenever the source states them.
    !***********************************************************************
    integer, parameter :: nProjectDefaults = 25
    !> The type spec pads: these section names are parameters of their own
    !> natural lengths, and an array constructor needs them equal.
    character(40), parameter :: defaultSect(nProjectDefaults) = [ &
        character(40) :: &
        sectProject, sectProject, sectProject, sectProject, &
        sectProject, sectProject, sectProject, sectProject, &
        sectProject, sectProject, sectProject, &
        sectSpectral, sectSpectral, &
        sectTilt, &
        sectTimelag, &
        sectPwb, sectPwb, sectPwb, sectPwb, sectPwb, sectPwb, &
        sectPwb, sectPwb, &
        sectBiomet, &
        sectSettings]
    character(26), parameter :: defaultKey(nProjectDefaults) = [ &
        'cec_meth                  ', 'cec_h                     ', &
        'cec_min_o1_o2             ', 'cec_min_octant            ', &
        'cec_min_valid             ', 'cec_signal_strength       ', &
        'cec_max_gap_fill          ', 'cec_max_stationarity      ', &
        'cec_singular_band         ', 'cec_stationarity_mode     ', &
        'cosp_model                ', &
        'automatic_spectra_config  ', 'flux_run_mode             ', &
        'rot_pf_assessment_only    ', &
        'tlag_assessment_only      ', &
        'pwb_n_bootstrap           ', 'pwb_block_length_s        ', &
        'pwb_min_valid_frac        ', 'pwb_hdi_thresh_s          ', &
        'pwb_dev_thresh_s          ', 'pwb_hdi_prefilter_s       ', &
        'pwb_smoothing_width       ', 'pwb_random_seed           ', &
        'biom_rh_override          ', &
        'tlag_meth                 ']
    character(8), parameter :: defaultValue(nProjectDefaults) = [ &
        '0       ', '0.000   ', '20.0    ', '5.0     ', &
        '90.0    ', '70.0    ', '4       ', '25.0    ', &
        '0.200   ', '0       ', &
        '0       ', &
        '0       ', '0       ', &
        '0       ', &
        '0       ', &
        '99      ', '20.0    ', '0.300   ', '0.50    ', &
        '0.50    ', '1.00    ', '5       ', '2024    ', &
        '0       ', &
        '2       ']

    !> One key=value line of an INI file, with the section it was found in.
    !> ParseIniFile drops the section and these keys have to be written back
    !> into the family that reads them, so the import keeps its own reader.
    type :: IniPair
        character(iniLabelLen) :: sect  = ''
        character(iniLabelLen) :: label = ''
        character(iniValueLen) :: value = ''
    end type IniPair

    !> A gas record under construction. `slot` is the EddyPro slot it came
    !> from, kept because every per-gas setting is looked up by that, and the
    !> record's own index is not it: records are compacted, so a site without
    !> CO2 has water at record one.
    type :: ImpGasRec
        integer :: slot = 0
        integer :: col = 0
        character(32) :: var = 'none'
        character(64) :: instr = 'none'
        character(iniValueLen) :: mw = ''
        character(iniValueLen) :: diff = ''
    end type ImpGasRec

    type :: ImpMeasRec
        integer :: col = 0
        character(32) :: var = 'none'
        character(64) :: instr = 'none'
    end type ImpMeasRec

contains

    !***********************************************************************
    !> Import an EddyPro project, if that is what we were handed.
    !>
    !> Returns with prjPath untouched for anything else, so this is safe to
    !> call unconditionally on whatever path InitEnv resolved - including the
    !> default ini/processing.eddyflow, which is how a project renamed to our
    !> extension is still recognised rather than half-parsed.
    !***********************************************************************
    subroutine ImportEddyProProject(prjPath, swVer)
        implicit none
        !> in/out variables
        character(*), intent(inout) :: prjPath
        !> This engine's version, for the files the import stamps. Passed in
        !> rather than included: version_and_date.inc also declares a build
        !> date, which nothing here has any use for.
        character(*), intent(in) :: swVer
        !> local variables
        type(IniPair), allocatable :: prjPairs(:)
        type(IniPair), allocatable :: mdPairs(:)
        type(ImpGasRec) :: gas(nEddyProGasSlots)
        type(ImpMeasRec) :: cell(MaxNumCellCols)
        type(ImpMeasRec) :: diag(MaxNumDiagCols)
        character(32) :: mdVar(MaxNumCol)
        character(64) :: mdInstr(MaxNumCol)
        character(PathLen) :: outPrj
        character(PathLen) :: outMd
        character(PathLen) :: srcMd
        character(19) :: stamp
        integer :: nPrj
        integer :: nMd
        integer :: nGas
        integer :: nCell
        integer :: nDiag
        logical :: haveMd
        logical :: ok

        if (.not. IsEddyProProject(prjPath)) return

        call ImportedPaths(prjPath, outPrj, outMd)

        !> An imported project that is already there is used as it stands.
        !> RP writes ex_file into it and FCC reads that back, so a second
        !> import - which is what invoking FCC with the same .eddypro would
        !> otherwise do - would overwrite that with the EddyPro file's own
        !> stale path, and FCC would read a FLUXNET file from another run.
        inquire(file = trim(outPrj), exist = ok)
        if (ok) then
            call LogSayList(' Reusing the imported project ' // trim(outPrj) // '.')
            call LogSayList(' Delete it to re-import from ' // trim(prjPath) // '.')
            prjPath = outPrj
            return
        end if

        call LogSayList(' Importing EddyPro project ' // trim(prjPath) // '.')

        call ReadIniPairs(prjPath, prjPairs, nPrj, ok)
        if (.not. ok) call ExceptionHandler(21)

        call EddyProMetadataPath(prjPath, prjPairs, nPrj, srcMd, haveMd)
        mdVar = 'none'
        mdInstr = 'none'
        nMd = 0
        if (haveMd) then
            call ReadIniPairs(srcMd, mdPairs, nMd, ok)
            if (.not. ok) call ExceptionHandler(22)
            call MetadataColumns(mdPairs, nMd, mdVar, mdInstr)
        else
            !> use_pfile=0: the metadata travels inside each GHG archive, so
            !> there is no file to import and no column table to resolve the
            !> records against. The records still carry the column numbers the
            !> project states; what they cannot carry is the analyser, and the
            !> fourth slot falls back to the name of the slot.
            call LogSayList(' This project reads its metadata from the raw files, so no metadata')
            call LogSayList(' file was imported. Gas records name no analyser.')
        end if

        call BuildRecords(prjPairs, nPrj, mdVar, mdInstr, &
            gas, nGas, cell, nCell, diag, nDiag)

        call CurrentStamp(stamp)
        if (haveMd) call WriteImportedMetadata(srcMd, outMd, stamp, swVer, &
            mdPairs, nMd)
        call WriteImportedProject(prjPairs, nPrj, outPrj, outMd, haveMd, stamp, swVer, &
            gas, nGas, cell, nCell, diag, nDiag)

        call LogSayList(' Imported as ' // trim(outPrj) // '.')
        prjPath = outPrj
    end subroutine ImportEddyProProject

    !***********************************************************************
    !> Whether a file declares itself an EddyPro project.
    !>
    !> The first line rather than the extension, so a project renamed to our
    !> extension - which the engine would otherwise half-parse and refuse with
    !> Fatal error(99), naming a cause that is not the real one - is still
    !> recognised for what it is.
    !***********************************************************************
    logical function IsEddyProProject(path)
        implicit none
        !> in/out variables
        character(*), intent(in) :: path
        !> local variables
        integer :: uf
        integer :: io_status
        character(ShortInstringLen) :: dataline

        IsEddyProProject = .false.
        open(newunit = uf, file = trim(path), status = 'old', iostat = io_status)
        if (io_status /= 0) return
        do
            read(uf, '(a)', iostat = io_status) dataline
            if (io_status /= 0) exit
            call stripstr(dataline)
            !> A byte order mark, if an editor has put one there. It is three
            !> bytes ahead of the tag, so without this the file is not
            !> recognised - and it fails as a *native* project rather than
            !> saying so, because the tag line is a comment either way and
            !> what is then missing is gas_num. Notepad writes one by default,
            !> and so does PowerShell's Set-Content -Encoding utf8.
            if (len_trim(dataline) >= 3) then
                if (dataline(1:3) == achar(239) // achar(187) // achar(191)) &
                    dataline = dataline(4:)
            end if
            if (len_trim(dataline) == 0) cycle
            IsEddyProProject = index(dataline, EddyProTag) == 1
            exit
        end do
        close(uf)
    end function IsEddyProProject

    !***********************************************************************
    !> Where the imported pair goes: beside the source, under a name that
    !> cannot collide with it.
    !***********************************************************************
    subroutine ImportedPaths(prjPath, outPrj, outMd)
        implicit none
        !> in/out variables
        character(*), intent(in) :: prjPath
        character(*), intent(out) :: outPrj
        character(*), intent(out) :: outMd
        !> local variables
        character(PathLen) :: stem
        integer :: dot

        stem = prjPath
        dot = index(stem, '.', .true.)
        !> Only a real extension, never a dot in a directory name on the way
        !> here: C:/my.data/project would otherwise lose half its path.
        if (dot > index(stem, slash, .true.)) stem = stem(1:dot - 1)

        outPrj = trim(stem) // ImportSuffix // '.eddyflow'
        outMd = trim(stem) // ImportSuffix // '.metadata'
    end subroutine ImportedPaths

    !***********************************************************************
    !> Read every key=value line of an INI file, remembering the section.
    !>
    !> The line rules are StoreIniTags': strip blanks, skip comments and empty
    !> lines, truncate at an inline semicolon, split on the first '='. What is
    !> added is the section, which the tag-table reader has no use for and
    !> this does - see the section constants above.
    !>
    !> Deliberately not ParseIniFile. The tag tables have retired precisely
    !> the keys an EddyPro file is made of - they are blanked in
    !> prj/gen_project_tags.py - so SearchLocalTags can never match one, and a
    !> second full tag table for EddyPro would be a hundred more generated
    !> slots to keep in step with these.
    !***********************************************************************
    subroutine ReadIniPairs(path, pairs, n, ok)
        implicit none
        !> in/out variables
        character(*), intent(in) :: path
        type(IniPair), allocatable, intent(out) :: pairs(:)
        integer, intent(out) :: n
        logical, intent(out) :: ok
        !> local variables
        integer :: uf
        integer :: io_status
        integer :: separ
        integer :: com
        character(iniLabelLen) :: sect
        character(ShortInstringLen) :: dataline

        n = 0
        ok = .false.
        !> Allocatable, not automatic: MaxNLinesIni of these is about nine
        !> megabytes, which is not a thing to put on the stack.
        allocate(pairs(MaxNLinesIni))

        open(newunit = uf, file = trim(path), status = 'old', iostat = io_status)
        if (io_status /= 0) return
        ok = .true.

        sect = ''
        do
            read(uf, '(a)', iostat = io_status) dataline
            if (io_status /= 0) exit
            call stripstr(dataline)
            if (len_trim(dataline) == 0) cycle
            if (dataline(1:1) == ';') cycle
            if (dataline(1:1) == '[') then
                com = index(dataline, ']')
                if (com > 2) then
                    sect = dataline(2:com - 1)
                else
                    sect = ''
                end if
                cycle
            end if
            com = index(dataline, ';')
            if (com /= 0) then
                dataline(com:len_trim(dataline)) = ''
                call stripstr(dataline)
            end if
            separ = index(dataline, '=')
            if (separ <= 1) cycle
            if (n >= MaxNLinesIni) then
                call ExceptionHandler(98)
                exit
            end if
            n = n + 1
            pairs(n)%sect = sect
            pairs(n)%label = dataline(1:separ - 1)
            pairs(n)%value = dataline(separ + 1:len_trim(dataline))
        end do
        close(uf)
    end subroutine ReadIniPairs

    !***********************************************************************
    !> The value of a key, by exact label. Blank and .false. if absent.
    !>
    !> Exact, never a substring. The tag tables learnt that the expensive way:
    !> err_label was captured by fluxnet_err_label, so the missing-value token
    !> written to every output file came from the wrong key.
    !***********************************************************************
    logical function PairValue(pairs, n, label, value)
        implicit none
        !> in/out variables
        type(IniPair), intent(in) :: pairs(:)
        integer, intent(in) :: n
        character(*), intent(in) :: label
        character(*), intent(out) :: value
        !> local variables
        integer :: i

        value = ''
        PairValue = .false.
        do i = 1, n
            if (trim(adjustl(pairs(i)%label)) == trim(adjustl(label))) then
                value = pairs(i)%value
                PairValue = .true.
                return
            end if
        end do
    end function PairValue

    !> The value of a key naming one of the four gas slots, in either of the
    !> two spellings EddyPro uses for the fourth.
    logical function SlotValue(pairs, n, pre, post, slot, value)
        implicit none
        !> in/out variables
        type(IniPair), intent(in) :: pairs(:)
        integer, intent(in) :: n
        character(*), intent(in) :: pre
        character(*), intent(in) :: post
        integer, intent(in) :: slot
        character(*), intent(out) :: value

        SlotValue = PairValue(pairs, n, &
            pre // trim(slotName(slot)) // post, value)
        if (SlotValue) return
        SlotValue = PairValue(pairs, n, &
            pre // trim(slotAlt(slot)) // post, value)
    end function SlotValue

    !> An integer key, or `dflt` when it is absent or unreadable.
    integer function PairInt(pairs, n, label, dflt)
        implicit none
        !> in/out variables
        type(IniPair), intent(in) :: pairs(:)
        integer, intent(in) :: n
        character(*), intent(in) :: label
        integer, intent(in) :: dflt
        !> local variables
        character(iniValueLen) :: value
        integer :: io_status

        PairInt = dflt
        if (.not. PairValue(pairs, n, label, value)) return
        read(value, *, iostat = io_status) PairInt
        if (io_status /= 0) PairInt = dflt
    end function PairInt

    !***********************************************************************
    !> Where the EddyPro metadata is: the path the project states, or the
    !> file of the project's own name beside it.
    !>
    !> The fallback is not a convenience. proj_file is an absolute path
    !> written on the machine the project was made on, so an EddyPro project
    !> that has been copied anywhere at all states a path that is not there -
    !> the example this was built against names C:/ABACOS/eddypro/... and
    !> carries its metadata in the same folder as itself.
    !***********************************************************************
    subroutine EddyProMetadataPath(prjPath, pairs, n, mdPath, found)
        implicit none
        !> in/out variables
        character(*), intent(in) :: prjPath
        type(IniPair), intent(in) :: pairs(:)
        integer, intent(in) :: n
        character(*), intent(out) :: mdPath
        logical, intent(out) :: found
        !> local variables
        character(iniValueLen) :: stated
        character(PathLen) :: sibling
        integer :: dot

        found = .false.
        mdPath = ''
        stated = ''

        !> use_pfile=0 means the metadata is embedded in the raw GHG files.
        if (PairInt(pairs, n, 'use_pfile', 1) /= 1) return

        if (PairValue(pairs, n, 'proj_file', stated)) then
            if (len_trim(stated) > 0) then
                mdPath = stated
                call AdjFilePath(mdPath, slash)
                inquire(file = trim(mdPath), exist = found)
                if (found) return
            end if
        end if

        sibling = prjPath
        dot = index(sibling, '.', .true.)
        if (dot > index(sibling, slash, .true.)) sibling = sibling(1:dot - 1)
        sibling = trim(sibling) // '.metadata'
        inquire(file = trim(sibling), exist = found)
        if (found) then
            call LogSayList(' The metadata this project names is not there; using')
            call LogSayList(' ' // trim(sibling) // ' instead.')
            mdPath = sibling
            return
        end if

        call AbortOnMissingPath('proj_file', trim(stated), &
            'The metadata file this project names does not exist, and no file ' // &
            'called ' // trim(sibling) // ' sits beside the project either. ' // &
            'Put the metadata file in one of those two places, or point ' // &
            'proj_file at where it actually is.')
    end subroutine EddyProMetadataPath

    !***********************************************************************
    !> The variable and instrument of every column the metadata describes.
    !>
    !> One pass over the pairs rather than two hundred searches through them,
    !> and the only thing the import needs the metadata for: it is what tells
    !> a record which analyser measured its column, and it is the only place
    !> that knows what species sits in the fourth slot.
    !***********************************************************************
    subroutine MetadataColumns(pairs, n, mdVar, mdInstr)
        implicit none
        !> in/out variables
        type(IniPair), intent(in) :: pairs(:)
        integer, intent(in) :: n
        character(*), intent(out) :: mdVar(MaxNumCol)
        character(*), intent(out) :: mdInstr(MaxNumCol)
        !> local variables
        integer :: i
        integer :: col
        integer :: under
        integer :: io_status
        character(iniLabelLen) :: label
        character(iniLabelLen) :: suffix

        mdVar = 'none'
        mdInstr = 'none'

        do i = 1, n
            label = adjustl(pairs(i)%label)
            if (label(1:4) /= 'col_') cycle
            under = index(label(5:len_trim(label)), '_')
            if (under <= 1) cycle
            read(label(5:under + 3), *, iostat = io_status) col
            if (io_status /= 0) cycle
            if (col < 1 .or. col > MaxNumCol) cycle
            suffix = label(under + 5:len_trim(label))
            select case (trim(suffix))
                case ('variable')
                    mdVar(col) = pairs(i)%value
                    !> The engine matches a species by its lower-case name,
                    !> and the interface writes COS with capitals.
                    call lowercase(mdVar(col))
                case ('instrument')
                    mdInstr(col) = &
                        CanonicalInstrumentModel(trim(adjustl(pairs(i)%value)))
            end select
        end do
    end subroutine MetadataColumns

    !***********************************************************************
    !> Turn EddyPro's fixed column slots into EddyFlow's records.
    !>
    !> Records are compacted: a slot the project did not name gets no record
    !> at all. That is what the interface does, and it is why nothing
    !> downstream may read a record positionally.
    !***********************************************************************
    subroutine BuildRecords(pairs, n, mdVar, mdInstr, &
            gas, nGas, cell, nCell, diag, nDiag)
        implicit none
        !> in/out variables
        type(IniPair), intent(in) :: pairs(:)
        integer, intent(in) :: n
        character(*), intent(in) :: mdVar(MaxNumCol)
        character(*), intent(in) :: mdInstr(MaxNumCol)
        type(ImpGasRec), intent(out) :: gas(nEddyProGasSlots)
        integer, intent(out) :: nGas
        type(ImpMeasRec), intent(out) :: cell(MaxNumCellCols)
        integer, intent(out) :: nCell
        type(ImpMeasRec), intent(out) :: diag(MaxNumDiagCols)
        integer, intent(out) :: nDiag
        !> local variables
        integer :: slot
        integer :: rec
        integer :: col
        integer :: io_status
        character(iniValueLen) :: value
        real(kind = dbl) :: aux

        nGas = 0
        do slot = 1, nEddyProGasSlots
            if (.not. SlotValue(pairs, n, 'col_', '', slot, value)) cycle
            read(value, *, iostat = io_status) col
            if (io_status /= 0) cycle
            if (col <= 0 .or. col > MaxNumCol) cycle

            nGas = nGas + 1
            gas(nGas)%slot = slot
            gas(nGas)%col = col
            gas(nGas)%instr = ColumnInstrument(mdInstr, col)
            gas(nGas)%var = trim(slotName(slot))

            !> The species of the first three slots is the slot. The fourth
            !> takes whatever the site measured, and the project file records
            !> that nowhere - so it comes from the metadata. The example this
            !> was built against states col_n2o=18 over a column the metadata
            !> calls COS, and a record written as n2o would hand a 60.075
            !> g/mol gas nitrous oxide's molecular weight.
            if (slot == nEddyProGasSlots) then
                if (len_trim(mdVar(col)) > 0) then
                    if (trim(mdVar(col)) /= 'none' .and. &
                        trim(mdVar(col)) /= 'ignore') gas(nGas)%var = mdVar(col)
                end if
            end if
        end do

        !> gas_mw and gas_diff are single-valued in EddyPro and describe the
        !> open slot alone. Their -1 is a sentinel meaning "use the built-in
        !> constant", and it has to reach the record as an ABSENT value: the
        !> engine tests mw > 0 to decide whether to consult its own tables, so
        !> a literal -1 would be read as a stated molecular weight.
        do rec = 1, nGas
            if (gas(rec)%slot /= nEddyProGasSlots) cycle
            if (PairValue(pairs, n, 'gas_mw', value)) then
                read(value, *, iostat = io_status) aux
                if (io_status == 0 .and. aux > 0d0) gas(rec)%mw = value
            end if
            if (PairValue(pairs, n, 'gas_diff', value)) then
                read(value, *, iostat = io_status) aux
                if (io_status == 0 .and. aux > 0d0) gas(rec)%diff = value
            end if
        end do

        nCell = 0
        call AddMeasRecord(pairs, n, mdInstr, 'col_cell_t', 'cell_t', cell, nCell)
        call AddMeasRecord(pairs, n, mdInstr, 'col_int_t_1', 'int_t_1', cell, nCell)
        call AddMeasRecord(pairs, n, mdInstr, 'col_int_t_2', 'int_t_2', cell, nCell)
        call AddMeasRecord(pairs, n, mdInstr, 'col_int_p', 'int_p', cell, nCell)

        nDiag = 0
        call AddMeasRecord(pairs, n, mdInstr, 'col_diag_75', 'diag_75', diag, nDiag)
        call AddMeasRecord(pairs, n, mdInstr, 'col_diag_72', 'diag_72', diag, nDiag)
        call AddMeasRecord(pairs, n, mdInstr, 'col_diag_77', 'diag_77', diag, nDiag)
        call AddMeasRecord(pairs, n, mdInstr, 'col_diag_anem', 'diag_anem', diag, nDiag)
    end subroutine BuildRecords

    !> One cell or diagnostic record, if the project names a column for it.
    subroutine AddMeasRecord(pairs, n, mdInstr, label, var, recs, nRecs)
        implicit none
        !> in/out variables
        type(IniPair), intent(in) :: pairs(:)
        integer, intent(in) :: n
        character(*), intent(in) :: mdInstr(MaxNumCol)
        character(*), intent(in) :: label
        character(*), intent(in) :: var
        type(ImpMeasRec), intent(inout) :: recs(:)
        integer, intent(inout) :: nRecs
        !> local variables
        integer :: col

        col = PairInt(pairs, n, label, 0)
        if (col <= 0 .or. col > MaxNumCol) return
        if (nRecs >= size(recs)) return

        nRecs = nRecs + 1
        recs(nRecs)%col = col
        recs(nRecs)%var = var
        recs(nRecs)%instr = ColumnInstrument(mdInstr, col)
    end subroutine AddMeasRecord

    !> The analyser a column belongs to, as the record has to spell it.
    character(64) function ColumnInstrument(mdInstr, col)
        implicit none
        !> in/out variables
        character(*), intent(in) :: mdInstr(MaxNumCol)
        integer, intent(in) :: col

        ColumnInstrument = 'none'
        if (col < 1 .or. col > MaxNumCol) return
        if (len_trim(mdInstr(col)) == 0) return
        ColumnInstrument = mdInstr(col)
    end function ColumnInstrument

    !***********************************************************************
    !> Write the imported project.
    !>
    !> Everything the source states is copied through, in its own section and
    !> in file order, except the keys the records replace. What is added is
    !> the records themselves and the per-gas settings, appended to whichever
    !> section reads them.
    !>
    !> Nothing EddyFlow-only is synthesised - not cec_*, not the PWB block,
    !> not the planar-fit sectors. Every one of them is *TagFound-guarded with
    !> a default, so an absent key is the documented "leave it alone", and
    !> writing the engine's defaults into the file would dress them up as
    !> decisions somebody made.
    !***********************************************************************
    subroutine WriteImportedProject(pairs, n, outPrj, outMd, haveMd, stamp, swVer, &
            gas, nGas, cell, nCell, diag, nDiag)
        implicit none
        !> in/out variables
        type(IniPair), intent(in) :: pairs(:)
        integer, intent(in) :: n
        character(*), intent(in) :: outPrj
        character(*), intent(in) :: outMd
        logical, intent(in) :: haveMd
        character(*), intent(in) :: stamp
        character(*), intent(in) :: swVer
        type(ImpGasRec), intent(in) :: gas(nEddyProGasSlots)
        integer, intent(in) :: nGas
        type(ImpMeasRec), intent(in) :: cell(MaxNumCellCols)
        integer, intent(in) :: nCell
        type(ImpMeasRec), intent(in) :: diag(MaxNumDiagCols)
        integer, intent(in) :: nDiag
        !> local variables
        integer, parameter :: MaxOverrides = 8
        integer :: uf
        integer :: io_status
        integer :: i
        integer :: nOvr
        logical :: usedOvr(MaxOverrides)
        character(iniLabelLen) :: ovrLabel(MaxOverrides)
        character(iniValueLen) :: ovrValue(MaxOverrides)
        character(iniLabelLen) :: sect
        character(iniValueLen) :: value
        character(PathLen) :: slashed
        character(40) :: pending
        logical :: seenProject
        logical :: seenSect(nProjectDefaults)
        integer :: k

        !> The keys whose value the import decides rather than carries. They
        !> are rewritten where the source states them and appended to
        !> [Project] where it does not, so an EddyPro file carrying no
        !> ini_version - which is legal; very old ones have none - still comes
        !> out stamped with the format it is now written in.
        nOvr = 0
        ovrLabel = ''
        ovrValue = ''
        slashed = outPrj
        call ForceSlash(slashed, .false.)
        call AddOverride(ovrLabel, ovrValue, nOvr, 'file_name', slashed)
        call AddOverride(ovrLabel, ovrValue, nOvr, 'ini_version', ImportedIniVer)
        call AddOverride(ovrLabel, ovrValue, nOvr, 'sw_version', swVer)
        call AddOverride(ovrLabel, ovrValue, nOvr, 'last_change_date', stamp)
        if (haveMd) then
            slashed = outMd
            call ForceSlash(slashed, .false.)
            call AddOverride(ovrLabel, ovrValue, nOvr, 'proj_file', slashed)
        end if
        !> The master sonic is named by model key, and it is compared against
        !> the model the metadata reader has already canonicalised
        !> (define_used_variables.f90). Left as EddyPro spelled it, a Campbell
        !> site matches no column at all: nothing is flagged master_sonic, and
        !> the sonic axis adjustment, the angle-of-attack and w-boost
        !> overrides and the anemometric block of DefineE2Set go with it.
        if (PairValue(pairs, n, 'master_sonic', value)) then
            call AddOverride(ovrLabel, ovrValue, nOvr, 'master_sonic', &
                CanonicalInstrumentModel(trim(adjustl(value))))
        end if
        usedOvr = .false.

        open(newunit = uf, file = trim(outPrj), status = 'replace', iostat = io_status)
        if (io_status /= 0) then
            call LogSayList(' Fatal error(112)> ' // trim(outPrj))
            call ExceptionHandler(112)
            return
        end if

        write(uf, '(a)') EddyFlowTag

        sect = ''
        seenProject = .false.
        seenSect = .false.
        do i = 1, n
            if (pairs(i)%sect /= sect) then
                call FlushSection(uf, sect, gas, nGas, cell, nCell, diag, nDiag, &
                    pairs, n, ovrLabel, ovrValue, nOvr, usedOvr)
                if (len_trim(sect) > 0) write(uf, '(a)') ''
                sect = pairs(i)%sect
                write(uf, '(a)') '[' // trim(sect) // ']'
                if (trim(sect) == sectProject) seenProject = .true.
                do k = 1, nProjectDefaults
                    if (trim(defaultSect(k)) == trim(sect)) seenSect(k) = .true.
                end do
            end if
            if (IsConsumedKey(pairs(i)%label)) cycle
            if (WriteOverride(uf, pairs(i)%label, ovrLabel, ovrValue, nOvr, usedOvr)) cycle
            write(uf, '(a)') trim(adjustl(pairs(i)%label)) // '=' &
                // trim(pairs(i)%value)
        end do
        call FlushSection(uf, sect, gas, nGas, cell, nCell, diag, nDiag, &
            pairs, n, ovrLabel, ovrValue, nOvr, usedOvr)

        !> A [Project] the source never opened still has to carry the records:
        !> without gas_num the file we just wrote is refused by Fatal
        !> error(99), which is the very thing this import exists to prevent.
        if (.not. seenProject) then
            write(uf, '(a)') ''
            write(uf, '(a)') '[' // sectProject // ']'
            call FlushSection(uf, sectProject, gas, nGas, cell, nCell, diag, nDiag, &
                pairs, n, ovrLabel, ovrValue, nOvr, usedOvr)
        end if
        call WriteUnusedOverrides(uf, ovrLabel, ovrValue, nOvr, usedOvr)

        !> Sections the source never opened at all. An EddyPro file has no
        !> [RawProcess_PWBTimelag_Settings] - the whole method postdates it -
        !> so its block is created here rather than filled in above.
        pending = ''
        do k = 1, nProjectDefaults
            if (seenSect(k)) cycle
            if (trim(defaultSect(k)) /= trim(pending)) then
                write(uf, '(a)') ''
                write(uf, '(a)') '[' // trim(defaultSect(k)) // ']'
                pending = defaultSect(k)
            end if
            write(uf, '(a)') trim(defaultKey(k)) // '=' // trim(defaultValue(k))
        end do

        close(uf)
    end subroutine WriteImportedProject

    !> Append whatever a section owes, at the point that section ends.
    subroutine FlushSection(uf, sect, gas, nGas, cell, nCell, diag, nDiag, &
            pairs, n, ovrLabel, ovrValue, nOvr, usedOvr)
        implicit none
        !> in/out variables
        integer, intent(in) :: uf
        character(*), intent(in) :: sect
        type(ImpGasRec), intent(in) :: gas(nEddyProGasSlots)
        integer, intent(in) :: nGas
        type(ImpMeasRec), intent(in) :: cell(MaxNumCellCols)
        integer, intent(in) :: nCell
        type(ImpMeasRec), intent(in) :: diag(MaxNumDiagCols)
        integer, intent(in) :: nDiag
        type(IniPair), intent(in) :: pairs(:)
        integer, intent(in) :: n
        character(*), intent(in) :: ovrLabel(:)
        character(*), intent(in) :: ovrValue(:)
        integer, intent(in) :: nOvr
        logical, intent(inout) :: usedOvr(:)

        select case (trim(sect))
            case (sectProject)
                call WriteRecords(uf, gas, nGas, cell, nCell, diag, nDiag)
                call WriteUnusedOverrides(uf, ovrLabel, ovrValue, nOvr, usedOvr)
            case (sectSettings, sectParams, sectTimelag, sectSpectral)
                call WritePerGasSettings(uf, trim(sect), pairs, n, gas, nGas)
        end select
        call WriteSectionDefaults(uf, trim(sect), pairs, n)
    end subroutine FlushSection

    !***********************************************************************
    !> The EddyFlow-only settings belonging to a section, as it ends.
    !>
    !> A key the source already states is left alone: this fills gaps, it does
    !> not overrule. That matters for a file that has been round-tripped
    !> through the interface once already.
    !***********************************************************************
    subroutine WriteSectionDefaults(uf, sect, pairs, n)
        implicit none
        !> in/out variables
        integer, intent(in) :: uf
        character(*), intent(in) :: sect
        type(IniPair), intent(in) :: pairs(:)
        integer, intent(in) :: n
        !> local variables
        integer :: k
        character(iniValueLen) :: value

        if (len_trim(sect) == 0) return
        do k = 1, nProjectDefaults
            if (trim(defaultSect(k)) /= trim(sect)) cycle
            if (PairValue(pairs, n, trim(defaultKey(k)), value)) cycle
            write(uf, '(a)') trim(defaultKey(k)) // '=' // trim(defaultValue(k))
        end do
    end subroutine WriteSectionDefaults

    !> The gas, cell and diagnostic records.
    subroutine WriteRecords(uf, gas, nGas, cell, nCell, diag, nDiag)
        implicit none
        !> in/out variables
        integer, intent(in) :: uf
        type(ImpGasRec), intent(in) :: gas(nEddyProGasSlots)
        integer, intent(in) :: nGas
        type(ImpMeasRec), intent(in) :: cell(MaxNumCellCols)
        integer, intent(in) :: nCell
        type(ImpMeasRec), intent(in) :: diag(MaxNumDiagCols)
        integer, intent(in) :: nDiag
        !> local variables
        integer :: i
        character(16) :: tag

        write(uf, '(a)') 'gas_num=' // trim(IntText(nGas))
        do i = 1, nGas
            tag = 'gas_' // trim(IntText(i)) // '_'
            write(uf, '(a)') trim(tag) // 'var=' // trim(gas(i)%var)
            write(uf, '(a)') trim(tag) // 'instr=' // trim(gas(i)%instr)
            write(uf, '(a)') trim(tag) // 'col=' // trim(IntText(gas(i)%col))
            !> 0 is "resolve it for me" for both, which is all EddyPro can
            !> mean: it has no way to say that a gas is corrected with a
            !> particular hygrometer or sits in a particular cell.
            write(uf, '(a)') trim(tag) // 'moist=0'
            write(uf, '(a)') trim(tag) // 'cell=0'
            !> Written empty rather than omitted: the engine walks the record
            !> block by stride, so every record has to keep the same shape.
            write(uf, '(a)') trim(tag) // 'mw=' // trim(gas(i)%mw)
            write(uf, '(a)') trim(tag) // 'diff=' // trim(gas(i)%diff)
        end do

        write(uf, '(a)') 'cell_num=' // trim(IntText(nCell))
        do i = 1, nCell
            tag = 'cell_' // trim(IntText(i)) // '_'
            write(uf, '(a)') trim(tag) // 'var=' // trim(cell(i)%var)
            write(uf, '(a)') trim(tag) // 'instr=' // trim(cell(i)%instr)
            write(uf, '(a)') trim(tag) // 'col=' // trim(IntText(cell(i)%col))
        end do

        write(uf, '(a)') 'diag_num=' // trim(IntText(nDiag))
        do i = 1, nDiag
            tag = 'diag_' // trim(IntText(i)) // '_'
            write(uf, '(a)') trim(tag) // 'var=' // trim(diag(i)%var)
            write(uf, '(a)') trim(tag) // 'instr=' // trim(diag(i)%instr)
            write(uf, '(a)') trim(tag) // 'col=' // trim(IntText(diag(i)%col))
        end do
    end subroutine WriteRecords

    !***********************************************************************
    !> The per-gas settings a section owns, one block per record.
    !>
    !> A setting the source does not state is not written. The engine applies
    !> a record override whenever the tag is PRESENT, so writing an absent
    !> setting as 0 would state a decision the project never made - and 0 is a
    !> meaningful absolute limit, a meaningful default time lag and a
    !> meaningful minimum flux.
    !***********************************************************************
    subroutine WritePerGasSettings(uf, sect, pairs, n, gas, nGas)
        implicit none
        !> in/out variables
        integer, intent(in) :: uf
        character(*), intent(in) :: sect
        type(IniPair), intent(in) :: pairs(:)
        integer, intent(in) :: n
        type(ImpGasRec), intent(in) :: gas(nEddyProGasSlots)
        integer, intent(in) :: nGas
        !> local variables
        integer :: i
        integer :: slot
        logical :: isWater
        character(16) :: tag
        character(iniValueLen) :: months

        do i = 1, nGas
            slot = gas(i)%slot
            tag = 'gas_' // trim(IntText(i)) // '_'
            !> Water is the exception in four places. Its minimum-flux
            !> counterpart is to_le_min_flux and its spectral QA/QC thresholds
            !> are the LE triple; none of the four is a per-gas quantity, and
            !> EddyPro has no h2o member for any of them.
            isWater = trim(gas(i)%var) == 'h2o'

            select case (sect)
                case (sectParams)
                    call EmitSlot(uf, pairs, n, tag, 'sr_lim', 'sr_lim_', '', slot)
                    call EmitSlot(uf, pairs, n, tag, 'al_min', 'al_', '_min', slot)
                    call EmitSlot(uf, pairs, n, tag, 'al_max', 'al_', '_max', slot)
                    call EmitSlot(uf, pairs, n, tag, 'ds_hf', 'ds_hf_', '', slot)
                    call EmitSlot(uf, pairs, n, tag, 'ds_sf', 'ds_sf_', '', slot)
                    call EmitSlot(uf, pairs, n, tag, 'tl_def', 'tl_def_', '', slot)
                case (sectSettings)
                    call EmitSlot(uf, pairs, n, tag, 'out_full_sp', &
                        'out_full_sp_', '', slot)
                    call EmitSlot(uf, pairs, n, tag, 'out_full_cosp_w', &
                        'out_full_cosp_w_', '', slot)
                    call EmitSlot(uf, pairs, n, tag, 'out_raw', 'out_raw_', '', slot)
                case (sectTimelag)
                    call EmitSlot(uf, pairs, n, tag, 'to_min_lag', 'to_', '_min_lag', slot)
                    call EmitSlot(uf, pairs, n, tag, 'to_max_lag', 'to_', '_max_lag', slot)
                    if (.not. isWater) &
                        call EmitSlot(uf, pairs, n, tag, 'to_min_flux', &
                            'to_', '_min_flux', slot)
                case (sectSpectral)
                    call EmitSlot(uf, pairs, n, tag, 'sa_fmin', 'sa_fmin_', '', slot)
                    call EmitSlot(uf, pairs, n, tag, 'sa_fmax', 'sa_fmax_', '', slot)
                    call EmitSlot(uf, pairs, n, tag, 'sa_hfn_fmin', 'sa_hfn_', '_fmin', slot)
                    if (.not. isWater) then
                        call EmitSlot(uf, pairs, n, tag, 'sa_min_st', 'sa_min_st_', '', slot)
                        call EmitSlot(uf, pairs, n, tag, 'sa_min_un', 'sa_min_un_', '', slot)
                        call EmitSlot(uf, pairs, n, tag, 'sa_max', 'sa_max_', '', slot)
                        call MonthGrouping(pairs, n, slot, months)
                        if (len_trim(months) > 0) &
                            write(uf, '(a)') trim(tag) // 'sa_months=' // trim(months)
                    end if
            end select
        end do
    end subroutine WritePerGasSettings

    !> One per-gas setting, written only if the source states it.
    subroutine EmitSlot(uf, pairs, n, tag, suffix, pre, post, slot)
        implicit none
        !> in/out variables
        integer, intent(in) :: uf
        type(IniPair), intent(in) :: pairs(:)
        integer, intent(in) :: n
        character(*), intent(in) :: tag
        character(*), intent(in) :: suffix
        character(*), intent(in) :: pre
        character(*), intent(in) :: post
        integer, intent(in) :: slot
        !> local variables
        character(iniValueLen) :: value

        if (.not. SlotValue(pairs, n, pre, post, slot, value)) return
        write(uf, '(a)') trim(tag) // suffix // '=' // trim(value)
    end subroutine EmitSlot

    !***********************************************************************
    !> The months a gas pools before a transfer function is fitted, as the
    !> single list EddyFlow reads, out of the twelve start/stop pairs EddyPro
    !> spread it over.
    !>
    !> One group spanning the calendar is a placeholder rather than a
    !> decision, and is deliberately not carried: it is what a gas stating
    !> nothing gets anyway, and writing it would turn a default into a choice.
    !***********************************************************************
    subroutine MonthGrouping(pairs, n, slot, months)
        implicit none
        !> in/out variables
        type(IniPair), intent(in) :: pairs(:)
        integer, intent(in) :: n
        integer, intent(in) :: slot
        character(*), intent(out) :: months
        !> local variables
        integer :: k
        integer :: nGroup
        integer :: mstart
        integer :: mstop

        months = ''
        nGroup = 0
        do k = 1, MaxGasClasses
            mstart = SlotIntKey(pairs, n, 'sa_', &
                '_g' // trim(IntText(k)) // '_start', slot)
            mstop = SlotIntKey(pairs, n, 'sa_', &
                '_g' // trim(IntText(k)) // '_stop', slot)
            if (mstart < 1 .or. mstop < 1) cycle
            if (mstart > 12 .or. mstop > 12 .or. mstop < mstart) cycle
            nGroup = nGroup + 1
            if (nGroup > 1) months = trim(months) // ','
            months = trim(months) // trim(IntText(mstart)) // '-' &
                // trim(IntText(mstop))
        end do

        if (nGroup == 1 .and. trim(months) == '1-12') months = ''
    end subroutine MonthGrouping

    !> A slot-named integer key, or -1 when it is absent or unreadable.
    integer function SlotIntKey(pairs, n, pre, post, slot)
        implicit none
        !> in/out variables
        type(IniPair), intent(in) :: pairs(:)
        integer, intent(in) :: n
        character(*), intent(in) :: pre
        character(*), intent(in) :: post
        integer, intent(in) :: slot
        !> local variables
        character(iniValueLen) :: value
        integer :: io_status

        SlotIntKey = -1
        if (.not. SlotValue(pairs, n, pre, post, slot, value)) return
        read(value, *, iostat = io_status) SlotIntKey
        if (io_status /= 0) SlotIntKey = -1
    end function SlotIntKey

    !***********************************************************************
    !> Whether a key is one the records replace, and so must not be copied.
    !>
    !> The engine already ignores every one of these - their slots are blanked
    !> in the tag tables - but a file stating col_co2=7 beside a record saying
    !> something else invites a reader to believe the wrong one, and it is the
    !> state the interface removes on its own side.
    !>
    !> Matched exactly, never as a prefix: sa_max_h, sa_max_le, sa_max_ustar
    !> and the sa_min_st / sa_min_un LE, H and u* thresholds share these
    !> prefixes, are whole-run settings, and have to survive.
    !***********************************************************************
    logical function IsConsumedKey(label)
        implicit none
        !> in/out variables
        character(*), intent(in) :: label
        !> local variables
        integer :: slot
        integer :: k
        character(iniLabelLen) :: key

        key = adjustl(label)
        IsConsumedKey = .true.

        select case (trim(key))
            !> The fixed slots the records replace, and the fourth gas's
            !> single molecular weight and diffusivity pair.
            case ('col_co2', 'col_h2o', 'col_ch4', 'col_n2o', 'col_gas4', &
                  'col_cell_t', 'col_int_t_1', 'col_int_t_2', 'col_int_p', &
                  'col_diag_72', 'col_diag_75', 'col_diag_77', 'col_diag_anem', &
                  'gas_mw', 'gas_diff')
                return
            !> The fixed full-output format, retired: it promised four gas
            !> blocks whatever the project held, so a fifth gas fell out of
            !> the file, and the flag now selects nothing.
            case ('fix_out_format')
                return
        end select

        do slot = 1, nEddyProGasSlots
            if (MatchesSlot(key, slot, 'sr_lim_', '')) return
            if (MatchesSlot(key, slot, 'al_', '_min')) return
            if (MatchesSlot(key, slot, 'al_', '_max')) return
            if (MatchesSlot(key, slot, 'ds_hf_', '')) return
            if (MatchesSlot(key, slot, 'ds_sf_', '')) return
            if (MatchesSlot(key, slot, 'tl_def_', '')) return
            if (MatchesSlot(key, slot, 'to_', '_min_lag')) return
            if (MatchesSlot(key, slot, 'to_', '_max_lag')) return
            if (MatchesSlot(key, slot, 'to_', '_min_flux')) return
            if (MatchesSlot(key, slot, 'pwb_', '_min_lag')) return
            if (MatchesSlot(key, slot, 'pwb_', '_max_lag')) return
            if (MatchesSlot(key, slot, 'out_full_sp_', '')) return
            if (MatchesSlot(key, slot, 'out_full_cosp_w_', '')) return
            if (MatchesSlot(key, slot, 'out_raw_', '')) return
            if (MatchesSlot(key, slot, 'sa_fmin_', '')) return
            if (MatchesSlot(key, slot, 'sa_fmax_', '')) return
            if (MatchesSlot(key, slot, 'sa_hfn_', '_fmin')) return
            if (MatchesSlot(key, slot, 'sa_min_st_', '')) return
            if (MatchesSlot(key, slot, 'sa_min_un_', '')) return
            if (MatchesSlot(key, slot, 'sa_max_', '')) return
            do k = 1, MaxGasClasses
                if (MatchesSlot(key, slot, 'sa_', &
                    '_g' // trim(IntText(k)) // '_start')) return
                if (MatchesSlot(key, slot, 'sa_', &
                    '_g' // trim(IntText(k)) // '_stop')) return
            end do
        end do

        IsConsumedKey = .false.
    end function IsConsumedKey

    !> Whether a key is `pre`, one of a slot's two spellings, and `post`.
    logical function MatchesSlot(key, slot, pre, post)
        implicit none
        !> in/out variables
        character(*), intent(in) :: key
        integer, intent(in) :: slot
        character(*), intent(in) :: pre
        character(*), intent(in) :: post

        MatchesSlot = trim(key) == pre // trim(slotName(slot)) // post &
                 .or. trim(key) == pre // trim(slotAlt(slot)) // post
    end function MatchesSlot

    !> Record an override, i.e. a key whose value the import decides.
    subroutine AddOverride(ovrLabel, ovrValue, nOvr, label, value)
        implicit none
        !> in/out variables
        character(*), intent(inout) :: ovrLabel(:)
        character(*), intent(inout) :: ovrValue(:)
        integer, intent(inout) :: nOvr
        character(*), intent(in) :: label
        character(*), intent(in) :: value

        if (nOvr >= size(ovrLabel)) return
        nOvr = nOvr + 1
        ovrLabel(nOvr) = label
        ovrValue(nOvr) = value
    end subroutine AddOverride

    !> Write a key's overridden value in place of the source's, if it has one.
    logical function WriteOverride(uf, label, ovrLabel, ovrValue, nOvr, usedOvr)
        implicit none
        !> in/out variables
        integer, intent(in) :: uf
        character(*), intent(in) :: label
        character(*), intent(in) :: ovrLabel(:)
        character(*), intent(in) :: ovrValue(:)
        integer, intent(in) :: nOvr
        logical, intent(inout) :: usedOvr(:)
        !> local variables
        integer :: i

        WriteOverride = .false.
        do i = 1, nOvr
            if (usedOvr(i)) cycle
            if (trim(adjustl(label)) /= trim(ovrLabel(i))) cycle
            write(uf, '(a)') trim(ovrLabel(i)) // '=' // trim(ovrValue(i))
            usedOvr(i) = .true.
            WriteOverride = .true.
            return
        end do
    end function WriteOverride

    !> Every override the source never gave us a line to rewrite.
    subroutine WriteUnusedOverrides(uf, ovrLabel, ovrValue, nOvr, usedOvr)
        implicit none
        !> in/out variables
        integer, intent(in) :: uf
        character(*), intent(in) :: ovrLabel(:)
        character(*), intent(in) :: ovrValue(:)
        integer, intent(in) :: nOvr
        logical, intent(inout) :: usedOvr(:)
        !> local variables
        integer :: i

        do i = 1, nOvr
            if (usedOvr(i)) cycle
            write(uf, '(a)') trim(ovrLabel(i)) // '=' // trim(ovrValue(i))
            usedOvr(i) = .true.
        end do
    end subroutine WriteUnusedOverrides

    !***********************************************************************
    !> Write the imported metadata.
    !>
    !> A faithful line-for-line copy, comments and blank lines included, with
    !> four values rewritten. The key SETS of an EddyPro and an EddyFlow
    !> metadata file are the same, so nothing has to be added; what differs is
    !> the instrument vocabulary, and that difference is not cosmetic.
    !>
    !> instr_<K>_model is canonicalised by the metadata reader already. Its
    !> col_<N>_instrument counterpart is NOT: columns are bound to instruments
    !> by index(LocCol(i)%instr_name, Instr(j)%model) - a raw name searched for
    !> a canonicalised one - so on a Campbell site written by EddyPro,
    !> index('csat3_1', 'csi_csat3_1') is zero, every anemometric column keeps
    !> a null instrument and metadata validation rejects the file. Rewriting
    !> both here is what lets such a site run at all.
    !***********************************************************************
    subroutine WriteImportedMetadata(srcMd, outMd, stamp, swVer, pairs, n)
        implicit none
        !> in/out variables
        character(*), intent(in) :: srcMd
        character(*), intent(in) :: outMd
        character(*), intent(in) :: stamp
        character(*), intent(in) :: swVer
        !> The source, already parsed. Used to learn which instruments and
        !> columns it describes, so the added keys cover exactly those.
        type(IniPair), intent(in) :: pairs(:)
        integer, intent(in) :: n
        !> local variables
        integer :: uin
        integer :: uout
        integer :: io_status
        integer :: separ
        integer :: com
        character(PathLen) :: slashed
        character(ShortInstringLen) :: dataline
        character(iniLabelLen) :: label
        character(iniLabelLen) :: sect
        logical :: firstLine

        open(newunit = uin, file = trim(srcMd), status = 'old', iostat = io_status)
        if (io_status /= 0) call ExceptionHandler(22)
        open(newunit = uout, file = trim(outMd), status = 'replace', iostat = io_status)
        if (io_status /= 0) then
            close(uin)
            call LogSayList(' Fatal error(112)> ' // trim(outMd))
            call ExceptionHandler(112)
            return
        end if

        slashed = outMd
        call ForceSlash(slashed, .false.)

        sect = ''
        firstLine = .true.
        do
            read(uin, '(a)', iostat = io_status) dataline
            if (io_status /= 0) exit
            separ = index(dataline, '=')
            if (separ <= 1) then
                !> A section header, a comment, or a blank line. The header is
                !> what the two appended blocks are placed by: the engine reads
                !> metadata with a blank section key and so does not care, but
                !> the interface reads it through QSettings, which groups by
                !> section - a key in the wrong group is a key it cannot see.
                if (index(dataline, '[') == 1) then
                    com = index(dataline, ']')
                    if (trim(sect) == mdSectInstruments) &
                        call WriteInstrumentExtras(uout, pairs, n)
                    if (com > 2) then
                        sect = dataline(2:com - 1)
                    else
                        sect = ''
                    end if
                end if
                !> The first line declares the format. saveProject emits this
                !> one unconditionally and never the EDDYFLOW spelling, which
                !> it only accepts on read - so state it, rather than carrying
                !> through whatever the source happened to say.
                if (firstLine) then
                    write(uout, '(a)') GhgMetadataTag
                else
                    write(uout, '(a)') trim(dataline)
                end if
                firstLine = .false.
                cycle
            end if
            firstLine = .false.
            label = adjustl(dataline(1:separ - 1))
            !> Exact labels: instr_2_sw_version is not sw_version, and a blind
            !> substring rewrite would stamp an analyser's firmware with the
            !> engine's version number.
            select case (trim(label))
                case ('file_name')
                    write(uout, '(a)') 'file_name=' // trim(slashed)
                case ('sw_version')
                    write(uout, '(a)') 'sw_version=' // trim(swVer)
                case ('last_change_date')
                    write(uout, '(a)') 'last_change_date=' // trim(stamp)
                case default
                    if (IsModelKey(label)) then
                        call NoteUnknownModel(dataline(separ + 1:len_trim(dataline)))
                        write(uout, '(a)') trim(label) // '=' // &
                            trim(CanonicalInstrumentModel( &
                                trim(adjustl(dataline(separ + 1:len_trim(dataline))))))
                    elseif (.not. WriteUpgradedConversionKey(uout, label, &
                            dataline, separ, pairs, n)) then
                        write(uout, '(a)') trim(dataline)
                    end if
            end select
        end do

        !> Whatever the last section was, it ends here.
        if (trim(sect) == mdSectInstruments) call WriteInstrumentExtras(uout, pairs, n)
        if (trim(sect) == mdSectColumns) call WriteColumnExtras(uout, pairs, n)

        close(uin)
        close(uout)
    end subroutine WriteImportedMetadata

    !***********************************************************************
    !> Handle the four keys that carry a column's linear conversion, dropping
    !> the two that are retired and upgrading the two that move.
    !>
    !> Returns .true. when it has dealt with the line (written a replacement,
    !> or deliberately written nothing), .false. to let the caller copy it
    !> through unchanged - which is what happens for every column that is not
    !> zero_fullscale, i.e. every column of every file the interface has
    !> written in years.
    !>
    !> min_value/max_value go unconditionally: they exist only to state a
    !> zero_fullscale input range, and the engine no longer computes with
    !> them. Where they carried a real range, it is folded into a_value and
    !> b_value here, by the same algebra WriteEddyFlowMetadataVariables
    !> applies on read:
    !>     gain   = (b - a)/(max - min)
    !>     offset = (a*max - b*min)/(max - min)
    !>
    !> An unconvertible range (max == min) is left exactly as it was, keys
    !> and all. The import's job is to carry a file forward, not to decide
    !> it is unusable - ColumnValidation already rejects that file, and it
    !> should reject the imported copy for the same stated reason rather
    !> than for a different one this rewrite invented.
    !***********************************************************************
    logical function WriteUpgradedConversionKey(uout, label, dataline, &
            separ, pairs, n) result(handled)
        implicit none
        !> in/out variables
        integer, intent(in) :: uout
        character(*), intent(in) :: label
        character(*), intent(in) :: dataline
        integer, intent(in) :: separ
        type(IniPair), intent(in) :: pairs(:)
        integer, intent(in) :: n
        !> local variables
        character(iniLabelLen) :: lab
        character(iniLabelLen) :: stem
        character(iniLabelLen) :: key
        character(iniValueLen) :: conv
        integer :: cut
        real(kind = dbl) :: cmin, cmax, ca, cb
        character(32) :: numtxt

        handled = .false.
        lab = trim(adjustl(label))
        if (index(lab, 'col_') /= 1) return

        !> The conversion key itself: the type is what changes name, so it
        !> is rewritten here rather than copied. Only when the range is
        !> convertible - see the note above on why max == min is left whole.
        cut = len_trim(lab) - len('conversion') + 1
        if (cut > 1) then
            if (lab(cut:len_trim(lab)) == 'conversion') then
                if (trim(adjustl(dataline(separ + 1:len_trim(dataline)))) &
                        /= 'zero_fullscale') return
                stem = lab(1:cut - 1)
                cmin = PairReal(pairs, n, trim(stem) // 'min_value', 0d0)
                cmax = PairReal(pairs, n, trim(stem) // 'max_value', 0d0)
                if (cmax == cmin) return
                write(uout, '(a)') trim(lab) // '=gain_offset'
                handled = .true.
                return
            end if
        end if

        !> Split col_<N>_<key> on the last underscore of the four keys this
        !> handles. Each is <something>_value, so the split is the
        !> underscore before that word - found by taking the last one of the
        !> whole label and then the last one of what precedes it.
        cut = index(trim(lab), '_', .true.)
        if (cut <= 4) return
        if (trim(lab(cut + 1:)) /= 'value') return
        cut = index(trim(lab(1:cut - 1)), '_', .true.)
        if (cut <= 4) return

        stem = lab(1:cut)                       !> 'col_<N>_'
        key  = lab(cut + 1:len_trim(lab))       !> 'min_value' etc.

        !> The two retired keys go whatever the column's conversion says.
        !> An imported file is a new file, and new files do not carry them -
        !> the interface stopped writing them for the same reason. Dropping
        !> them from a column that was never zero_fullscale loses nothing:
        !> they had no other meaning, and every such column in every file
        !> the interface has written states them as a constant zero.
        if (trim(key) == 'min_value' .or. trim(key) == 'max_value') then
            handled = .true.
            return
        end if

        if (.not. PairValue(pairs, n, trim(stem) // 'conversion', conv)) return
        if (trim(adjustl(conv)) /= 'zero_fullscale') return

        cmin = PairReal(pairs, n, trim(stem) // 'min_value', 0d0)
        cmax = PairReal(pairs, n, trim(stem) // 'max_value', 0d0)
        ca   = PairReal(pairs, n, trim(stem) // 'a_value', 0d0)
        cb   = PairReal(pairs, n, trim(stem) // 'b_value', 0d0)
        if (cmax == cmin) return

        select case (trim(key))
            case ('a_value')
                write(numtxt, '(f0.6)') (cb - ca) / (cmax - cmin)
                write(uout, '(a)') trim(lab) // '=' // trim(adjustl(numtxt))
                handled = .true.
            case ('b_value')
                write(numtxt, '(f0.6)') &
                    (ca * cmax - cb * cmin) / (cmax - cmin)
                write(uout, '(a)') trim(lab) // '=' // trim(adjustl(numtxt))
                handled = .true.
        end select
    end function WriteUpgradedConversionKey

    !> A real key, or `dflt` when it is absent or unreadable.
    real(kind = dbl) function PairReal(pairs, n, label, dflt) result(val)
        implicit none
        !> in/out variables
        type(IniPair), intent(in) :: pairs(:)
        integer, intent(in) :: n
        character(*), intent(in) :: label
        real(kind = dbl), intent(in) :: dflt
        !> local variables
        character(iniValueLen) :: value
        integer :: io_status

        val = dflt
        if (.not. PairValue(pairs, n, label, value)) return
        read(value, *, iostat = io_status) val
        if (io_status /= 0) val = dflt
    end function PairReal

    !***********************************************************************
    !> The per-instrument sampling keys EddyFlow added, for every instrument
    !> the source describes.
    !>
    !> Both are inert as written. `ac_freq` of 0 means "the file's rate":
    !> ColumnAcFreq falls back to Metadata%ac_freq for anything <= 0, which is
    !> also where an absent key leaves it, since NullInstrument pre-sets the
    !> error code. `integrates` of 0 is instantaneous, the same .false. an
    !> absent key leaves behind.
    !>
    !> Stating them is nonetheless worth doing beyond format tidiness: read
    !> blind, an absent numeric tag carries whatever the shared ANTags array
    !> last held. That read is guarded now, but a file that says what it means
    !> does not depend on the guard.
    !***********************************************************************
    subroutine WriteInstrumentExtras(uf, pairs, n)
        implicit none
        !> in/out variables
        integer, intent(in) :: uf
        type(IniPair), intent(in) :: pairs(:)
        integer, intent(in) :: n
        !> local variables
        integer :: k
        logical :: wrote
        character(iniValueLen) :: value
        character(24) :: tag

        wrote = .false.
        do k = 1, MaxNumInstruments
            tag = 'instr_' // trim(IntText(k)) // '_'
            !> Only instruments the file actually describes, and by index
            !> rather than by count: nothing says the blocks are contiguous.
            if (.not. PairValue(pairs, n, trim(tag) // 'model', value)) cycle
            if (.not. PairValue(pairs, n, trim(tag) // 'ac_freq', value)) then
                write(uf, '(a)') trim(tag) // 'ac_freq=' // acFreqDefault
                wrote = .true.
            end if
            if (.not. PairValue(pairs, n, trim(tag) // 'integrates', value)) then
                write(uf, '(a)') trim(tag) // 'integrates=' // integratesDefault
                wrote = .true.
            end if
        end do
        !> The section that follows opens on the next line otherwise, which no
        !> other section in the file does. Purely how it reads.
        if (wrote) write(uf, '(a)') ''
    end subroutine WriteInstrumentExtras

    !***********************************************************************
    !> The per-column fill declaration EddyFlow added, for every column the
    !> source describes.
    !>
    !> Inert: -9999 IS the engine's error code, and BlankMissingValues guards
    !> its extra test on `err_value /= error`, so a column declaring the
    !> conventional fill takes the same branch as one declaring nothing.
    !> The interface writes it for every column including `ignore`, so this
    !> does too.
    !***********************************************************************
    subroutine WriteColumnExtras(uf, pairs, n)
        implicit none
        !> in/out variables
        integer, intent(in) :: uf
        type(IniPair), intent(in) :: pairs(:)
        integer, intent(in) :: n
        !> local variables
        integer :: c
        character(iniValueLen) :: value
        character(24) :: tag

        do c = 1, MaxNumCol
            tag = 'col_' // trim(IntText(c)) // '_'
            if (.not. PairValue(pairs, n, trim(tag) // 'variable', value)) cycle
            if (.not. PairValue(pairs, n, trim(tag) // 'error_value', value)) &
                write(uf, '(a)') trim(tag) // 'error_value=' // errorValueDefault
        end do
    end subroutine WriteColumnExtras

    !> Whether a metadata key names an instrument model.
    logical function IsModelKey(label)
        implicit none
        !> in/out variables
        character(*), intent(in) :: label
        !> local variables
        character(iniLabelLen) :: key
        integer :: under

        key = adjustl(label)
        IsModelKey = .false.
        if (len_trim(key) < 7) return
        if (key(1:6) == 'instr_') then
            IsModelKey = index(key, '_model') == len_trim(key) - 5
            return
        end if
        if (key(1:4) /= 'col_') return
        under = index(key(5:len_trim(key)), '_')
        if (under <= 1) return
        IsModelKey = trim(key(under + 5:len_trim(key))) == 'instrument'
    end function IsModelKey

    !***********************************************************************
    !> Say so when a metadata file names a model this engine does not know.
    !>
    !> The vocabulary here is EddyPro 6.2.2's plus this fork's additions, and
    !> an EddyPro 7 file can name a model neither knows. Left unsaid, that
    !> surfaces only as ExceptionHandler(25) once per raw file, with nothing
    !> to say which key was at fault - and for a GHG archive that is every
    !> file in the run.
    !>
    !> Held against metadata_file_validation.f90's own lists by
    !> static_checks/test_eddypro_import_static.py, so a model added there and
    !> not here cannot start warning about itself.
    !***********************************************************************
    subroutine NoteUnknownModel(model)
        implicit none
        !> in/out variables
        character(*), intent(in) :: model
        !> local variables
        character(32) :: base

        if (len_trim(model) == 0) return
        base = InstrumentModelBase(CanonicalInstrumentModel(trim(adjustl(model))))

        select case (trim(base))
            case ('hs_50', 'hs_100', 'r2', 'r3_50', 'r3_100', 'r3a_100', &
                  'wm', 'wmpro', 'usa1_standard', 'usa1_fast', &
                  'usoni3_classa_mp', 'usoni3_cage_mp', &
                  'csi_csat3', 'csi_csat3a', 'csi_csat3b', 'csi_csat3c', &
                  'csi_irgason_sonic', &
                  '81000', '81000v', '81000re', '81000vre', 'generic_sonic')
                return
            case ('li6262', 'li7000', 'li7200', 'li7200rs', 'li7500', &
                  'li7500a', 'li7500rs', 'li7500ds', 'li7700', &
                  'csi_ec150', 'csi_ec155', 'csi_irgason_irga', 'csi_tga200a', &
                  'miro_mga1_5', 'miro_mga4_6', 'miro_mga9_10', 'miro_mgai_n2o', &
                  'aerodyne_tildas', &
                  'generic_open_path', 'generic_closed_path', &
                  'open_path_krypton', 'open_path_lyman', &
                  'closed_path_krypton', 'closed_path_lyman')
                return
            !> An unassigned column names no instrument, which is not a fault.
            case ('none')
                return
        end select

        call LogSayList(' The metadata names an instrument model this engine does not know: ' &
            // trim(adjustl(model)) // '.')
        call LogSayList(' Every column assigned to it will fail metadata validation.')
    end subroutine NoteUnknownModel

    !> Now, in the form the project and metadata files write dates.
    subroutine CurrentStamp(stamp)
        implicit none
        !> in/out variables
        character(*), intent(out) :: stamp
        !> local variables
        character(32) :: timestring

        call hms_current_string(timestring)
        !> The hour comes out space-padded, which is not a timestamp.
        if (timestring(12:12) == ' ') timestring(12:12) = '0'
        stamp = timestring(1:10) // 'T' // timestring(12:19)
    end subroutine CurrentStamp

    !> An integer as its shortest text, for building keys and values.
    character(16) function IntText(i)
        implicit none
        !> in/out variables
        integer, intent(in) :: i

        write(IntText, '(i0)') i
    end function IntText

end module m_eddypro_import
