"""An EddyPro project is converted by the same table the interface uses.

The engine can now be handed a *.eddypro file. It reads it, writes an ordinary
EddyFlow pair beside it and runs that, so nothing downstream has to know the
project did not start life in this format.

Almost all of that conversion is a rename, and a rename is the kind of thing
that rots without saying so: nothing fails to compile, nothing raises, the run
simply proceeds with a setting the project did not make. The two tables worth
pinning are

  * which of the two spellings EddyPro uses for its fourth gas slot in which
    family - n2o in the parameter settings and in out_full_sp_* /
    out_full_cosp_w_*, gas4 in the spectral, time-lag and out_raw_* families -
    which is asserted against the interface's own copy, because the two
    repositories are built separately and nothing at build time connects them;

  * where an instrument model key is canonicalised. `csat3` and `csat3b` became
    `csi_csat3` and `csi_csat3b`, and a model key appears in three places: the
    metadata's instr_<K>_model, the metadata's col_<N>_instrument, and the
    project's master_sonic. Only the first is canonicalised by the reader. The
    other two are compared by index() against a model that HAS been - see
    read_metadata_file.f90 and define_used_variables.f90 - so a Campbell site
    whose import misses either one loses its anemometer entirely, and the
    failure surfaces as a metadata rejection naming no cause.

The rest of the file guards the decisions that are easy to undo by tidying:
records are compacted, an unstated setting is not written, water takes none of
the four LE/H thresholds, and the retired flat keys are dropped rather than
carried alongside the records that replace them.
"""

import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GUI = ROOT.parent / "eddyflow-gui"

IMPORT = "src/src_common/m_eddypro_import.f90"
INIT_ENV = "src/src_common/init_env.f90"
VALIDATION = "src/src_common/metadata_file_validation.f90"
TYPEDEF = "src/src_common/m_typedef.f90"
GEN_TAGS = "prj/gen_project_tags.py"

#: The interface's copy of the same tables, read rather than transcribed.
GUI_ECPROJECT = "src/ecproject.cpp"

#: Which slot spelling belongs to which family, as the interface states it.
#: The engine accepts BOTH spellings everywhere, which is a superset of this -
#: what is asserted is that it accepts at least what the interface writes.
GUI_SLOT_TABLES = {
    "paramSlot": ("co2", "h2o", "ch4", "n2o"),
    "specSlot": ("co2", "h2o", "ch4", "gas4"),
    "outSpSlot": ("co2", "h2o", "ch4", "n2o"),
    "outRawSlot": ("co2", "h2o", "ch4", "gas4"),
}


def read(rel, root=ROOT):
    return (root / rel).read_text(encoding="utf-8", errors="replace")


def code(rel):
    """The source with its comment lines removed.

    Every claim below is about what the code does, and this file is written in
    the house style - long prose comments naming the very keys under test. A
    substring search over the raw text would find `sr_lim_n2o` in a paragraph
    explaining it and pass with the code deleted.
    """
    return "\n".join(ln for ln in read(rel).splitlines()
                     if not ln.strip().startswith("!"))


def body(src, start, end):
    text = src[src.index(start):]
    return text[:text.index(end)]


class TheSlotSpellings(unittest.TestCase):
    """Both of EddyPro's spellings for the fourth slot are understood."""

    def setUp(self):
        self.src = code(IMPORT)

    def test_both_spelling_tables_exist(self):
        #> slotName carries what EddyPro writes in most families, slotAlt the
        #> other spelling, and every lookup tries the one then the other.
        self.assertIn("'co2 ', 'h2o ', 'ch4 ', 'n2o '", self.src)
        self.assertIn("'co2 ', 'h2o ', 'ch4 ', 'gas4'", self.src)

    def test_every_lookup_tries_both(self):
        #> SlotValue is the single door to a slot-named key. If a second door
        #> is ever added that consults slotName alone, half the fourth-slot
        #> settings become invisible depending on the family.
        fn = body(self.src, "logical function SlotValue(", "end function SlotValue")
        self.assertEqual(fn.count("slotName(slot)"), 1)
        self.assertEqual(fn.count("slotAlt(slot)"), 1)
        #> Only two places build a key from a slot, and both consult both
        #> tables. A third that consulted slotName alone would make half the
        #> fourth-slot settings visible or not depending on the family.
        #>
        #> BuildRecords is the exception and is deliberate: there it is a
        #> species name, not half a key, which is why it reads slotName only.
        #> `gas4` is not a species and must never be written as one.
        allowed = ("SlotValue", "MatchesSlot", "BuildRecords")
        for m in re.finditer(r"slot(?:Name|Alt)\(slot\)", self.src):
            owner = re.findall(r"(?:function|subroutine) (\w+)\(",
                               self.src[:m.start()])
            self.assertIn(owner[-1], allowed,
                          "a slot table is read outside %s" % (allowed,))
        build = body(self.src, "subroutine BuildRecords(",
                     "end subroutine BuildRecords")
        self.assertNotIn("slotAlt", build)

    @unittest.skipUnless((GUI / GUI_ECPROJECT).exists(),
                         "eddyflow-gui not checked out beside this repository")
    def test_the_interface_writes_nothing_we_cannot_read(self):
        gui = read(GUI_ECPROJECT, GUI)
        for name, expected in GUI_SLOT_TABLES.items():
            m = re.search(r"const char\* %s\[4\] = \{([^}]*)\}" % name, gui)
            self.assertIsNotNone(m, "%s is gone from the interface" % name)
            found = tuple(re.findall(r'"([^"]+)"', m.group(1)))
            self.assertEqual(found, expected,
                             "%s changed on the interface side" % name)
            #> Whatever it spells, we accept it: the fourth entry is in one of
            #> the two tables this engine tries.
            self.assertIn(found[3], ("n2o", "gas4"))


class TheModelKeySites(unittest.TestCase):
    """All three places a model key can appear are canonicalised."""

    def setUp(self):
        self.src = code(IMPORT)

    def test_the_metadata_model_and_column_instrument_are_rewritten(self):
        #> IsModelKey is what decides, and it has to answer for both shapes.
        fn = body(self.src, "logical function IsModelKey(", "end function IsModelKey")
        self.assertIn("'instr_'", fn)
        self.assertIn("'_model'", fn)
        self.assertIn("'col_'", fn)
        self.assertIn("'instrument'", fn)
        #> ...and the writer has to act on its answer.
        writer = body(self.src, "subroutine WriteImportedMetadata(",
                      "end subroutine WriteImportedMetadata")
        self.assertIn("IsModelKey(label)", writer)
        self.assertIn("CanonicalInstrumentModel(", writer)

    def test_the_projects_master_sonic_is_rewritten(self):
        writer = body(self.src, "subroutine WriteImportedProject(",
                      "end subroutine WriteImportedProject")
        self.assertIn("'master_sonic'", writer)
        m = re.search(r"'master_sonic'[^\n]*\n?[^\n]*CanonicalInstrumentModel",
                      writer)
        self.assertIsNotNone(
            m, "master_sonic is carried across without canonicalisation")

    def test_the_gas_records_name_a_canonical_analyser(self):
        #> A record's instrument comes from the metadata column table, so that
        #> table is where the rewrite has to happen for the records too.
        fn = body(self.src, "subroutine MetadataColumns(",
                  "end subroutine MetadataColumns")
        self.assertIn("CanonicalInstrumentModel(", fn)

    def test_the_mapping_itself_is_not_duplicated_here(self):
        #> There is exactly one table of retired spellings, in m_typedef, and
        #> exactly one static check owning it. A bare Campbell key appearing
        #> here would be a second copy of that mapping to keep in step.
        #>
        #> The canonical names are a different matter: NoteUnknownModel has to
        #> name them, because it answers "would metadata validation accept
        #> this?" and validation only ever sees canonical ones.
        for bare in ("'csat3'", "'csat3a'", "'csat3b'", "'ec150'"):
            self.assertNotIn(bare, self.src,
                             "%s is a retired spelling; the mapping from it "
                             "belongs to CanonicalInstrumentModel alone" % bare)


class TheKnownModelList(unittest.TestCase):
    """The importer warns about a model the validator would reject."""

    def setUp(self):
        self.fn = body(code(IMPORT), "subroutine NoteUnknownModel(",
                       "end subroutine NoteUnknownModel")

    def test_every_model_the_validator_accepts_is_known_here(self):
        #> Otherwise the import warns about a perfectly good instrument, which
        #> is worse than saying nothing: it teaches the reader to ignore it.
        validator = code(VALIDATION)
        models = set()
        for m in re.finditer(r"LocInstr%model\(1:len_trim\(LocInstr%model\)-2\)",
                             validator):
            block = validator[m.start():]
            block = block[:block.index("end select")]
            #> Only the case labels. The arms carry strings of their own -
            #> the krypton and Lyman-alpha arm tests LocCol%var against
            #> 'h2o' - and those are not instrument models.
            for label in re.finditer(r"case \(([^)]*)\)", block):
                models.update(re.findall(r"'([a-z0-9_]+)'", label.group(1)))
        models.discard("generic_sonic")
        self.assertGreater(len(models), 20, "the validator lists were not found")
        for model in sorted(models):
            self.assertIn("'%s'" % model, self.fn,
                          "%s is accepted by metadata validation but the "
                          "import would call it unknown" % model)


class TheRecords(unittest.TestCase):
    """The four fixed slots become a compacted list of records."""

    def setUp(self):
        self.src = code(IMPORT)
        self.build = body(self.src, "subroutine BuildRecords(",
                          "end subroutine BuildRecords")

    def test_a_slot_without_a_column_gets_no_record(self):
        #> Compaction is what makes a record's index mean nothing about its
        #> species, which the whole record format rests on.
        self.assertIn("if (col <= 0 .or. col > MaxNumCol) cycle", self.build)

    def test_the_fourth_slots_species_comes_from_the_metadata(self):
        #> The project file records it nowhere. Written as n2o, a COS column
        #> would be given nitrous oxide's molecular weight - the two differ by
        #> a third, and nothing in the output would look unusual.
        self.assertIn("mdVar(col)", self.build)
        self.assertIn("nEddyProGasSlots", self.build)

    def test_the_single_mw_and_diff_pair_belong_to_the_fourth_slot(self):
        self.assertIn("'gas_mw'", self.build)
        self.assertIn("'gas_diff'", self.build)
        self.assertIn("gas(rec)%slot /= nEddyProGasSlots", self.build)

    def test_an_unset_mw_reaches_the_record_as_absent(self):
        #> EddyPro's -1 means "use the built-in constant" and the engine tests
        #> mw > 0 to decide the same thing, so a literal -1 would be read as a
        #> stated molecular weight of minus one gram per mole.
        self.assertIn("aux > 0d0", self.build)

    def test_every_record_keeps_the_same_shape(self):
        #> The engine walks the record block by stride, so mw and diff are
        #> written empty rather than omitted.
        recs = body(self.src, "subroutine WriteRecords(", "end subroutine WriteRecords")
        self.assertIn("'mw=' // trim(gas(i)%mw)", recs)
        self.assertIn("'diff=' // trim(gas(i)%diff)", recs)

    def test_the_cell_and_diagnostic_slots_are_all_covered(self):
        for key in ("col_cell_t", "col_int_t_1", "col_int_t_2", "col_int_p",
                    "col_diag_75", "col_diag_72", "col_diag_77", "col_diag_anem"):
            self.assertIn("'%s'" % key, self.build)


class ThePerGasSettings(unittest.TestCase):
    """A setting is written only where the project made a decision."""

    def setUp(self):
        self.src = code(IMPORT)
        self.fn = body(self.src, "subroutine WritePerGasSettings(",
                       "end subroutine WritePerGasSettings")

    def test_water_takes_none_of_the_four_le_and_h_thresholds(self):
        #> Water's minimum-flux counterpart is to_le_min_flux and its spectral
        #> QA/QC thresholds are the LE triple. None is a per-gas quantity, and
        #> EddyPro has no h2o member for any of the four.
        self.assertIn("isWater = trim(gas(i)%var) == 'h2o'", self.fn)
        for setting in ("'to_min_flux'", "'sa_min_st'", "'sa_min_un'", "'sa_max'"):
            before = self.fn[:self.fn.index(setting)]
            self.assertIn(".not. isWater", before,
                          "%s is not behind the water carve-out" % setting)

    def test_an_absent_setting_is_not_written(self):
        #> The engine applies a record override whenever the tag is PRESENT,
        #> so a setting written as 0 states a decision the project never made -
        #> and 0 is a meaningful absolute limit, default time lag and minimum
        #> flux.
        emit = body(self.src, "subroutine EmitSlot(", "end subroutine EmitSlot")
        self.assertIn("if (.not. SlotValue(", emit)
        self.assertIn(") return", emit)

    def test_every_family_reaches_the_section_that_reads_it(self):
        #> ReadIniRP sweeps RawProcess* and ReadIniFCC sweeps FluxCorrection*,
        #> each as one scope, so a per-gas key in the wrong family is swept by
        #> neither and reads downstream as "the project stated nothing".
        families = {
            "sectParams": ("sr_lim", "al_min", "al_max", "ds_hf", "ds_sf", "tl_def"),
            "sectSettings": ("out_full_sp", "out_full_cosp_w", "out_raw"),
            "sectTimelag": ("to_min_lag", "to_max_lag", "to_min_flux"),
            "sectSpectral": ("sa_fmin", "sa_fmax", "sa_hfn_fmin",
                             "sa_min_st", "sa_min_un", "sa_max"),
        }
        arms = self.fn.split("case (")
        for section, settings in families.items():
            arm = [a for a in arms if a.startswith(section)]
            self.assertEqual(len(arm), 1, "no %s arm" % section)
            for setting in settings:
                self.assertIn("'%s'" % setting, arm[0],
                              "%s is not written into %s" % (setting, section))

    def test_a_calendar_wide_month_group_is_a_placeholder(self):
        #> It is what a gas stating nothing gets anyway, so carrying it would
        #> turn a default into a choice.
        fn = body(self.src, "subroutine MonthGrouping(", "end subroutine MonthGrouping")
        self.assertIn("nGroup == 1", fn)
        self.assertIn("'1-12'", fn)


class TheRetiredKeys(unittest.TestCase):
    """What the records replace is dropped, and nothing else is."""

    def setUp(self):
        self.fn = body(code(IMPORT), "logical function IsConsumedKey(",
                       "end function IsConsumedKey")

    def _retired(self):
        """The retired names, from the generator that owns them."""
        gen = read(GEN_TAGS)
        block = body(gen, "RETIRED_LABELS = {", "}")
        return set(re.findall(r'"([a-z0-9_]+)"', block))

    def test_every_retired_flat_key_is_consumed(self):
        for key in self._retired():
            self.assertIn("'%s'" % key, self.fn,
                          "%s is retired but would be copied through" % key)

    def test_the_whole_run_thresholds_are_not_swept_up(self):
        #> sa_max_h, sa_max_le, sa_max_ustar and the sa_min_st / sa_min_un LE,
        #> H and u* thresholds share a prefix with the per-gas families and are
        #> whole-run settings. Matching has to be exact, never a prefix.
        matcher = body(code(IMPORT), "logical function MatchesSlot(",
                       "end function MatchesSlot")
        self.assertIn("trim(key) ==", matcher)
        self.assertNotIn("index(", matcher)

    def test_the_still_live_column_keys_survive(self):
        #> col_ts, col_air_t and col_air_p are one per project, not one per
        #> instrument, and no record replaces them.
        for key in ("col_ts", "col_air_t", "col_air_p"):
            self.assertNotIn("'%s'" % key, self.fn)


class TheHandshake(unittest.TestCase):
    """The engine reaches the import, and the import writes our format."""

    def test_an_eddypro_path_is_not_discarded_on_the_command_line(self):
        env = code(INIT_ENV)
        self.assertIn(".eddypro", env)
        self.assertIn("ImportEddyProProject(PrjPath", env)

    def test_the_import_decides_by_the_files_own_first_line(self):
        #> So a project renamed to our extension is recognised too, rather than
        #> half-parsed and refused with Fatal error(99) naming a cause that is
        #> not the real one.
        src = code(IMPORT)
        self.assertIn("';EDDYPRO_PROCESSING'", src)
        fn = body(src, "logical function IsEddyProProject(",
                  "end function IsEddyProProject")
        self.assertIn("EddyProTag", fn)
        #> A byte order mark sits three bytes ahead of the tag. Without the
        #> guard the file is read as a native project instead, and what is
        #> reported is the missing gas_num rather than the real cause.
        self.assertIn("achar(239) // achar(187) // achar(191)", fn)

    def test_the_written_project_declares_our_format_and_version(self):
        src = code(IMPORT)
        self.assertIn("';EDDYFLOW_PROCESSING'", src)
        #> Not pinned to a literal: what matters is that the two agree, and
        #> re-stating the number here means every format bump edits a test that
        #> is not about the number.
        #>
        #> A file written newer than the engine reads is refused by Fatal
        #> error(96) the moment it is read back, which would make the import
        #> produce something it cannot run. The import stamps what it writes,
        #> so it moves with the format: it now emits cec_singular_band and the
        #> pairing keys, which a 5.0.0 reader does not know.
        import re as _re
        stamped = _re.search(r"ImportedIniVer = '([\d.]+)'", src)
        self.assertIsNotNone(stamped, "the import does not stamp a version")
        supported = _re.search(
            r"MaxSupportedIniVer = '([\d.]+)'",
            read("src/src_common/m_common_global_var.f90"))
        self.assertIsNotNone(supported)
        self.assertEqual(stamped.group(1), supported.group(1),
                         "the import writes a version the engine will refuse")

    def test_the_imported_pair_cannot_overwrite_its_source(self):
        #> <base>.metadata is the name the EddyPro metadata already carries in
        #> the same folder.
        src = code(IMPORT)
        self.assertIn("ImportSuffix = '_ep_imported'", src)

    def test_the_import_runs_once(self):
        #> RP edits the project it ran - EditIniFile on ex_file - and FCC reads
        #> that edit back. Importing again on the FCC run would discard it.
        fn = body(code(IMPORT), "subroutine ImportEddyProProject(",
                  "end subroutine ImportEddyProProject")
        self.assertIn("inquire(file = trim(outPrj), exist =", fn)
        self.assertIn("prjPath = outPrj", fn)

    def test_the_source_files_are_never_written(self):
        src = code(IMPORT)
        for stmt in re.findall(r"open\([^)]*status = '(\w+)'", src):
            self.assertIn(stmt, ("old", "replace"))
        #> The only 'replace' opens name the two imported paths.
        for m in re.finditer(r"open\(newunit = \w+, file = trim\((\w+)\),"
                             r" status = 'replace'", src):
            self.assertIn(m.group(1), ("outPrj", "outMd"))

    def test_the_eddyflow_only_settings_are_written(self):
        #> The import produces a complete EddyFlow project, not one the
        #> interface has to migrate on open. Each of these is written at the
        #> value this engine already applies when the key is absent, which is
        #> what makes writing them a no-op - see TheProjectDefaults below,
        #> which checks the values against the readers that own them.
        src = code(IMPORT)
        for key in ("cec_meth", "cec_h", "cec_min_o1_o2", "cec_min_octant",
                    "cec_min_valid", "cec_signal_strength", "cec_max_gap_fill",
                    "cec_max_stationarity", "automatic_spectra_config",
                    "flux_run_mode", "rot_pf_assessment_only",
                    "tlag_assessment_only", "biom_rh_override",
                    "cec_singular_band", "cosp_model",
                    "pwb_n_bootstrap", "pwb_block_length_s", "pwb_min_valid_frac",
                    "pwb_hdi_thresh_s", "pwb_dev_thresh_s", "pwb_hdi_prefilter_s",
                    "pwb_smoothing_width", "pwb_random_seed"):
            self.assertIn("'%s" % key, src, "%s is not written" % key)
        #> And the one entry that is NOT a no-op, kept separate here for the
        #> same reason it is kept separate in the table it checks.
        self.assertIn("'tlag_meth", src, "tlag_meth is not written")

    def test_a_setting_whose_absence_is_a_decision_is_left_absent(self):
        #> The other half of the rule, and the half that is easy to lose to a
        #> later tidy-up. For each of these there is NO value that reproduces
        #> what the engine does when the key is missing:
        #>
        #>   pf_sect_*, wdf_sect_*   presence CREATES a wind sector -
        #>                           read_ini_rp.f90 counts them, so writing
        #>                           twelve would replace the Alert(40)
        #>                           "forcing to 1 sector of 360" path. They
        #>                           are EddyPro keys in any case and are
        #>                           copied through when the source states them
        #>   gas_<i>_pwb_*_lag       presence sets lag_bounds_provided, and
        #>                           set_timelags.f90 takes the window from
        #>                           instrument geometry without it
        #>   gas_<i>_drift_*         absent means the gas is not drift-
        #>                           corrected; there is no default curve
        #>   instr_<K>_max_lack      absent means "follow max_lack"
        #>   gas_<i>_fluxnet_default absent means the lowest record of the
        #>                           species is designated
        src = code(IMPORT)
        for key in ("pf_sect_", "wdf_sect_", "pwb_min_lag", "pwb_max_lag",
                    "drift_dir_", "drift_inv_", "max_lack", "fluxnet_default"):
            self.assertNotIn("'%s" % key, src,
                             "%s is written, but its absence is a decision" % key)


class TheProjectDefaults(unittest.TestCase):
    """Each written default is the one its reader already applies."""

    READERS = ("src/src_rp/read_ini_rp.f90",
               "src/src_fcc/read_ini_fcc.f90",
               "src/src_common/write_processing_project_variables.f90")

    #: key -> the assignment that sets it when the tag is absent, spelled as
    #: the reader spells it. Checked against the readers rather than
    #: transcribed, so a default that moves there fails here instead of
    #: quietly leaving the import writing a value nothing applies any more.
    OWNED = {
        "cec_h": "EddyFlowProj%cec%h = 0d0",
        "cec_min_o1_o2": "EddyFlowProj%cec%min_o1_o2 = 0.20d0",
        "cec_min_octant": "EddyFlowProj%cec%min_octant = 0.05d0",
        "cec_min_valid": "EddyFlowProj%cec%min_valid = 0.90d0",
        "cec_signal_strength": "EddyFlowProj%cec%signal_strength = 70d0",
        "cec_max_stationarity": "EddyFlowProj%cec%max_stationarity = 25d0",
        "cec_max_gap_fill": "EddyFlowProj%cec%max_gap_fill = 4",
        "cec_singular_band": "EddyFlowProj%cec%singular_band = 0.2d0",
        "pwb_n_bootstrap": "PWBSetup%n_bootstrap = 99",
        "pwb_block_length_s": "PWBSetup%block_length_s = 20d0",
        "pwb_min_valid_frac": "PWBSetup%min_valid_frac = 0.3d0",
        "pwb_hdi_thresh_s": "PWBSetup%hdi_thresh_s = 0.5d0",
        "pwb_dev_thresh_s": "PWBSetup%dev_thresh_s = 0.5d0",
        "pwb_hdi_prefilter_s": "PWBSetup%hdi_prefilter_s = 1.0d0",
        "pwb_smoothing_width": "PWBSetup%smoothing_width = 5",
        "pwb_random_seed": "PWBSetup%random_seed = 2024",
        "automatic_spectra_config": "FCCsetup%SA%automatic_config = .false.",
        "cosp_model": "EddyFlowProj%cosp_model = 'moncrieff_97'",
    }

    def written(self):
        """What the import writes, parsed out of its own table."""
        src = read(IMPORT)
        keys = re.search(r"defaultKey\(nProjectDefaults\) = \[(.*?)\]",
                         src, re.S).group(1)
        vals = re.search(r"defaultValue\(nProjectDefaults\) = \[(.*?)\]",
                         src, re.S).group(1)
        k = [x.strip() for x in re.findall(r"'([^']*)'", keys)]
        v = [x.strip() for x in re.findall(r"'([^']*)'", vals)]
        self.assertEqual(len(k), len(v), "key and value tables differ in length")
        return dict(zip(k, v))

    def test_the_tables_line_up(self):
        table = self.written()
        n = int(re.search(r"nProjectDefaults = (\d+)", read(IMPORT)).group(1))
        self.assertEqual(len(table), n)

    def test_every_default_still_exists_in_its_reader(self):
        #> Not the value itself - that is a formatting question - but the
        #> assignment it was taken from. If this fails, the engine default
        #> moved and the import is now writing something nothing applies.
        joined = "\n".join(read(r) for r in self.READERS)
        joined = "\n".join(ln for ln in joined.splitlines()
                           if not ln.strip().startswith("!"))
        for key, assignment in self.OWNED.items():
            self.assertIn(assignment, joined,
                          "%s: the reader no longer sets it that way, so the "
                          "value the import writes may not be the default"
                          % key)

    def test_the_numbers_agree_with_those_assignments(self):
        #> The fractions are the trap: the engine holds 0.20 and the file
        #> states 20.0, because NormalizeCecFraction accepts either form.
        #> Assert the relation rather than the literal.
        table = self.written()
        for key, fraction in (("cec_min_o1_o2", 0.20),
                              ("cec_min_octant", 0.05),
                              ("cec_min_valid", 0.90)):
            self.assertAlmostEqual(float(table[key]) / 100.0, fraction, places=6,
                                   msg="%s is not the reader default" % key)
        for key, expected in (("cec_h", 0.0), ("cec_signal_strength", 70.0),
                              ("cec_max_stationarity", 25.0),
                              ("cec_max_gap_fill", 4),
                              ("cec_singular_band", 0.2),
                              ("pwb_n_bootstrap", 99),
                              ("pwb_block_length_s", 20.0),
                              ("pwb_min_valid_frac", 0.3),
                              ("pwb_hdi_thresh_s", 0.5),
                              ("pwb_dev_thresh_s", 0.5),
                              ("pwb_hdi_prefilter_s", 1.0),
                              ("pwb_smoothing_width", 5),
                              ("pwb_random_seed", 2024)):
            self.assertAlmostEqual(float(table[key]), float(expected), places=6,
                                   msg="%s disagrees with its reader" % key)
        #> The flags are all off, which is what every one of those readers
        #> initialises before it looks for the tag.
        #>
        #> cec_meth carries a requirement of its own beyond that: a project
        #> arriving from EddyPro has had no chance to say which CO2 channel
        #> pairs with which water, so the partition must not start itself.
        #> TheImportedProjectStartsWithoutCec below covers the rest of it.
        for key in ("cec_meth", "automatic_spectra_config", "flux_run_mode",
                    "rot_pf_assessment_only", "tlag_assessment_only",
                    "biom_rh_override"):
            self.assertEqual(table[key], "0", "%s should be off" % key)

    def test_the_one_corrective_default_is_not_claimed_to_be_a_no_op(self):
        """tlag_meth is the exception, and it has to stay visibly one.

        Every other entry states the value the engine already applies when the
        key is absent, so writing it changes nothing. This one states 2 where
        the engine would have applied 'none', deliberately - so it must not be
        in OWNED, where the checks above would try to find a reader assignment
        agreeing with it and would be asserting something untrue.
        """
        self.assertEqual(self.written()["tlag_meth"], "2")
        self.assertNotIn("tlag_meth", self.OWNED)

    def test_the_retired_pwb_keys_are_not_written(self):
        #> pwb_detect_prewpl chose whether detection saw the gas series before
        #> or after the pointwise mixing-ratio conversion. Both alternatives
        #> ran on rotated 20 Hz data and the conversion runs before time-lag
        #> compensation, so detecting after it puts cell temperature and water
        #> into the gas series at the wrong relative lag. There is no choice
        #> left to write.
        #>
        #> pwb_approx_ccf and pwb_max_ar_order were speed options that bought
        #> under a percent of runtime and cost accuracy for it.
        #>
        #> A converted project simply does not carry them, which is what
        #> absent has always meant in this table.
        table = self.written()
        for key in ("pwb_detect_prewpl", "pwb_detect_on_raw",
                    "pwb_approx_ccf", "pwb_max_ar_order"):
            self.assertNotIn(key, table)


class TheMetadataAdditions(unittest.TestCase):
    """The three key families EddyFlow added, at their inert values."""

    def setUp(self):
        self.src = code(IMPORT)

    def test_the_values_are_the_ones_the_engine_assumes(self):
        self.assertIn("acFreqDefault = '0.000'", self.src)
        self.assertIn("integratesDefault = '0'", self.src)
        self.assertIn("errorValueDefault = '-9999.0000'", self.src)

    def test_they_land_inside_the_right_sections(self):
        #> The engine reads metadata with a blank section key and does not
        #> care, but the interface reads it through QSettings, which groups by
        #> section - a key in the wrong group is one it cannot see.
        self.assertIn("mdSectInstruments = 'Instruments'", self.src)
        self.assertIn("mdSectColumns = 'FileDescription'", self.src)
        writer = body(self.src, "subroutine WriteImportedMetadata(",
                      "end subroutine WriteImportedMetadata")
        self.assertIn("WriteInstrumentExtras", writer)
        self.assertIn("WriteColumnExtras", writer)

    def test_only_indices_the_source_describes(self):
        #> Nothing guarantees the blocks are contiguous, so drive off the keys
        #> the file actually states rather than off a count.
        fn = body(self.src, "subroutine WriteInstrumentExtras(",
                  "end subroutine WriteInstrumentExtras")
        self.assertIn("'model'", fn)
        fn = body(self.src, "subroutine WriteColumnExtras(",
                  "end subroutine WriteColumnExtras")
        self.assertIn("'variable'", fn)

    def test_a_key_the_source_states_is_not_overruled(self):
        for name in ("WriteInstrumentExtras", "WriteColumnExtras"):
            fn = body(self.src, "subroutine %s(" % name,
                      "end subroutine %s" % name)
            self.assertIn(".not. PairValue(", fn)

    def test_the_format_tag_is_stated_not_carried(self):
        self.assertIn("GhgMetadataTag = ';GHG_METADATA'", self.src)

    def test_the_reader_no_longer_reads_those_two_blind(self):
        #> SearchLocalTags clears CharTags%value but not NumTags%value, so an
        #> unguarded numeric read returns whatever the shared array last held
        #> - the previous raw file's value, for a GHG archive.
        md = code("src/src_common/read_metadata_file.f90")
        for offset in (13, 14):
            self.assertIn(
                "if (ANTagFound(init_an_instr + i*leap_an_instr + %d))" % offset,
                md)


class TheFixtures(unittest.TestCase):
    """The regression pair exists and exercises what it was built for."""

    EP = Path("tests/regression/base_ep.eddypro")
    EP_MD = Path("tests/regression/base_ep.metadata")
    NATIVE = Path("tests/regression/base_ep_native.eddyflow")

    def setUp(self):
        for rel in (self.EP, self.EP_MD, self.NATIVE):
            if not (ROOT / rel).exists():
                self.skipTest("%s is not checked out" % rel)

    def test_the_eddypro_fixture_is_one(self):
        self.assertTrue(read(self.EP).startswith(";EDDYPRO_PROCESSING"))
        self.assertTrue(read(self.NATIVE).startswith(";EDDYFLOW_PROCESSING"))

    def test_the_anemometer_is_a_campbell(self):
        #> The rename bites on Campbell and nowhere else, so a Gill fixture
        #> would pass whether or not any of the three rewrites happened.
        md = read(self.EP_MD)
        self.assertIn("instr_1_model=csat3b_1", md)
        self.assertIn("col_1_instrument=csat3b_1", md)
        self.assertIn("master_sonic=csat3b_1", read(self.EP))
        #> ...and the native twin describes the same instrument, spelled the
        #> way the import must produce.
        self.assertIn("master_sonic=csi_csat3b_1", read(self.NATIVE))

    def test_it_covers_both_fourth_slot_spellings(self):
        ep = read(self.EP)
        self.assertIn("al_n2o_max=", ep)
        self.assertIn("sa_fmax_gas4=", ep)
        self.assertIn("out_full_cosp_w_n2o=", ep)
        self.assertIn("out_raw_gas4=", ep)

    def test_it_covers_compaction_and_the_fourth_slots_species(self):
        ep = read(self.EP)
        #> A named-but-columnless gas, so the record list is shorter than four.
        self.assertIn("col_ch4=0", ep)
        native = read(self.NATIVE)
        self.assertIn("gas_num=3", native)
        #> The fourth slot holds COS, which only the metadata knows.
        self.assertIn("gas_3_var=cos", native)
        self.assertIn("gas_3_mw=60.0750", native)

    def test_it_covers_a_real_month_grouping(self):
        self.assertIn("sa_co2_g2_start=7", read(self.EP))
        self.assertIn("gas_1_sa_months=1-6,7-12", read(self.NATIVE))

    def test_the_native_twin_is_a_complete_eddyflow_project(self):
        #> It is what the import must produce, so it carries the same
        #> EddyFlow-only settings - and none of the ones whose absence is a
        #> decision, or the pair would diff for a reason that is not the
        #> conversion.
        native = read(self.NATIVE)
        for key in ("cec_meth=0", "automatic_spectra_config=0",
                    "flux_run_mode=0", "rot_pf_assessment_only=0",
                    "tlag_assessment_only=0", "biom_rh_override=0",
                    "pwb_n_bootstrap=99", "pwb_smoothing_width=5"):
            self.assertIn(key, native)
        self.assertIn("[RawProcess_PWBTimelag_Settings]", native)
        for key in ("pf_sect_1_width=", "wdf_sect_1_start=",
                    "gas_1_pwb_min_lag=", "instr_1_max_lack="):
            self.assertNotIn(key, native)

    def test_the_native_twin_metadata_carries_the_added_keys(self):
        md = read(Path("tests/regression/base_ep_native.metadata"))
        self.assertTrue(md.startswith(";GHG_METADATA"))
        for key in ("instr_1_ac_freq=0.000", "instr_1_integrates=0",
                    "col_1_error_value=-9999.0000"):
            self.assertIn(key, md)



class TheImportedProjectStartsWithoutCec(unittest.TestCase):
    """A project arriving from EddyPro must not start partitioning by itself.

    Conditional Eddy Covariance needs to know which CO2 channel is paired with
    which water channel, and a project that has only ever existed as an EddyPro
    file has had no chance to say. It is written off, and the pairing list is
    deliberately not written at all - absent means "derive one pairing per CO2
    channel from the analyser layout", which is a decision for the interface to
    make once someone opens the converted project, not for the importer.
    """

    EP = Path("tests/regression/base_ep.eddypro")

    def test_the_key_is_written_off(self):
        #> Asserted against the table rather than the fixture, so it holds for
        #> every EddyPro file rather than for the one in the tree.
        src = read(IMPORT)
        keys = re.search(r"defaultKey\(nProjectDefaults\) = \[(.*?)\]",
                         src, re.S).group(1)
        vals = re.search(r"defaultValue\(nProjectDefaults\) = \[(.*?)\]",
                         src, re.S).group(1)
        table = dict(zip([x.strip() for x in re.findall(r"'([^']*)'", keys)],
                         [x.strip() for x in re.findall(r"'([^']*)'", vals)]))
        self.assertEqual(table["cec_meth"], "0")

    def test_no_eddypro_file_can_state_it_so_the_default_always_fires(self):
        """The gap-fill only fires when the source is silent, so this is the
        half that makes it unconditional in practice."""
        if not (ROOT / self.EP).exists():
            self.skipTest("%s is not checked out" % self.EP)
        lines = [ln for ln in read(self.EP).splitlines()
                 if ln.strip().startswith("cec_")]
        self.assertEqual(lines, [],
                         "the EddyPro fixture states a CEC key, so it no "
                         "longer shows that a real EddyPro project leaves the "
                         "importer's default to fire")

    def test_the_pairing_list_is_left_for_the_interface(self):
        #> Writing cec_num=0 would say "no pairings", which is a different
        #> thing from "you decide" - see the reader's own note in
        #> write_processing_project_variables.f90.
        src = code(IMPORT)
        for key in ("cec_num", "cec_1_meth", "cec_1_co2", "cec_1_h2o",
                    "cec_1_extra"):
            self.assertNotIn("'%s" % key, src,
                             "%s must not be written: an absent pairing list "
                             "means the layout decides" % key)


class TheImportedProjectKeepsEddyProsTimeLag(unittest.TestCase):
    """Whatever method the source chose is what the converted project runs.

    The importer has three ways to interfere with a key - discard it in
    IsConsumedKey, replace it through the override list, or fill it in from the
    defaults table when absent. tlag_meth is in none of the first two, so a
    stated value reaches the written file untouched, in its own section, with
    its own value.

    It is in the third, at 2, and only for the case where the source states
    nothing: the reader has no *TagFound guard for this key and would otherwise
    fall through to 'none' and run with no time-lag compensation at all.
    """

    EP = Path("tests/regression/base_ep.eddypro")
    NATIVE = Path("tests/regression/base_ep_native.eddyflow")

    @staticmethod
    def tlag_select():
        """Just the time-lag select case.

        read_ini_rp.f90 has several selects over a single character, and the
        detrending one comes first - so searching the whole file for
        `case ('2')` lands on Meth%det and asserts nothing about time lags.
        """
        reader = read("src/src_rp/read_ini_rp.f90")
        start = reader.index("select case (SCTags(16)%value(1:1))")
        return reader[start:reader.index("end select", start)]

    def test_a_stated_value_is_copied_through(self):
        src = code(IMPORT)

        consumed = src[src.index("logical function IsConsumedKey("):]
        consumed = consumed[:consumed.index("end function IsConsumedKey")]
        self.assertNotIn("tlag_meth", consumed,
                         "discarding it would drop the source's choice")

        #> The overrides are added in one run of AddOverride calls; naming
        #> tlag_meth among them would replace the source's value with ours.
        added = re.findall(r"call AddOverride\([^)]*'([a-z_]+)'", src)
        self.assertNotIn("tlag_meth", added)
        self.assertIn("ini_version", added, "the override list moved")

    def test_the_regression_pair_still_covers_that(self):
        """It only does while the fixture states the key."""
        for rel in (self.EP, self.NATIVE):
            if not (ROOT / rel).exists():
                self.skipTest("%s is not checked out" % rel)
        for rel in (self.EP, self.NATIVE):
            stated = [ln.strip() for ln in read(rel).splitlines()
                      if ln.strip().startswith("tlag_meth=")]
            self.assertEqual(stated, ["tlag_meth=2"],
                             "%s no longer states tlag_meth once, so the pair "
                             "stops proving the source's value survives" % rel)

    def test_the_importer_can_never_produce_pwb(self):
        """EddyPro has no fifth method, and the import must not invent one."""
        src = read(IMPORT)
        keys = re.search(r"defaultKey\(nProjectDefaults\) = \[(.*?)\]",
                         src, re.S).group(1)
        vals = re.search(r"defaultValue\(nProjectDefaults\) = \[(.*?)\]",
                         src, re.S).group(1)
        table = dict(zip([x.strip() for x in re.findall(r"'([^']*)'", keys)],
                         [x.strip() for x in re.findall(r"'([^']*)'", vals)]))
        self.assertNotEqual(table["tlag_meth"], "5")

        #> And five is still what PWB is, so the number above means what this
        #> test thinks it means.
        block = self.tlag_select()
        self.assertIn("Meth%tlag = 'pwb'", block[block.index("case ('5')"):])

    def test_the_absent_case_is_filled_with_eddypros_own_default(self):
        src = read(IMPORT)
        keys = re.search(r"defaultKey\(nProjectDefaults\) = \[(.*?)\]",
                         src, re.S).group(1)
        vals = re.search(r"defaultValue\(nProjectDefaults\) = \[(.*?)\]",
                         src, re.S).group(1)
        sects = re.search(r"defaultSect\(nProjectDefaults\) = \[(.*?)\]",
                          src, re.S).group(1)
        k = [x.strip() for x in re.findall(r"'([^']*)'", keys)]
        v = [x.strip() for x in re.findall(r"'([^']*)'", vals)]
        names = re.findall(r"\bsect[A-Z]\w*", sects)
        self.assertEqual(len(names), len(k), "defaultSect is the wrong length")

        at = k.index("tlag_meth")
        self.assertEqual(v[at], "2")
        self.assertEqual(names[at], "sectSettings",
                         "tlag_meth lives in [RawProcess_Settings]; filling it "
                         "into another section would leave the reader's own "
                         "section still empty")

        #> Two is covariance maximization with default, and the reader still
        #> maps it there.
        block = self.tlag_select()
        arm = block[block.index("case ('2')"):block.index("case ('3')")]
        self.assertIn("Meth%tlag = 'maxcov&default'", arm)

        #> And the hole this fills is still open: no *TagFound guard, and the
        #> fall-through is 'none'. If either changes, this default may have
        #> stopped being needed.
        self.assertNotIn("SCTagFound(16)", read("src/src_rp/read_ini_rp.f90"))
        self.assertIn("Meth%tlag = 'none'",
                      block[block.index("case default"):])

    def test_the_hole_it_fills_is_described_where_it_is_filled(self):
        """The entry breaks the table's own rule, so the comment has to say
        which entry and why - otherwise the next reader takes it for a no-op
        like the other twenty-two."""
        src = read(IMPORT)
        head = src[:src.index("integer, parameter :: nProjectDefaults")]
        self.assertIn("KIND TWO", head)
        self.assertIn("tlag_meth", head)
        self.assertIn("case default", head)




class TheTwoImportersAgree(unittest.TestCase):
    """An EddyPro project must convert the same way whichever door it comes in.

    There are two importers. The engine's, in m_eddypro_import.f90, rewrites the
    file and runs it. The interface's is `EcProject::importEddyProProject`,
    which renames the fourth-gas keys, hands the rest to `loadEcProject` - where
    an absent key takes `defaultEcProjectState` - and then saves the result as a
    native .eddyflow.

    So the two arrive at the same answer by different routes: the engine writes
    a literal into the file, the interface applies its own default and then
    writes that. They agree today, and nothing said so. This is the check that
    notices when one of them moves.
    """

    GUI_STATE = "src/ecprojectstate.h"
    GUI_PROJECT = "src/ecproject.cpp"

    def engine_default(self, key):
        src = read(IMPORT)
        keys = re.search(r"defaultKey\(nProjectDefaults\) = \[(.*?)\]",
                         src, re.S).group(1)
        vals = re.search(r"defaultValue\(nProjectDefaults\) = \[(.*?)\]",
                         src, re.S).group(1)
        table = dict(zip([x.strip() for x in re.findall(r"'([^']*)'", keys)],
                         [x.strip() for x in re.findall(r"'([^']*)'", vals)]))
        return table[key]

    def gui_default(self, field):
        state = read(self.GUI_STATE, GUI)
        m = re.search(r"int %s = (-?\d+);" % field, state)
        self.assertIsNotNone(m, "%s is gone from the interface's state" % field)
        return m.group(1)

    @unittest.skipUnless((GUI / "src/ecprojectstate.h").exists(),
                         "eddyflow-gui not checked out beside this repository")
    def test_the_time_lag_fallback_is_the_same_number(self):
        """Both routes must land an EddyPro file that states no tlag_meth on
        the same method.

        The engine writes the literal because its own reader would otherwise
        fall through to 'none'; the interface never had that problem, because
        an absent key takes its default. The two numbers have to match or the
        same project computes different fluxes depending on whether it was
        opened or run.
        """
        self.assertEqual(self.engine_default("tlag_meth"),
                         self.gui_default("tlag_meth"),
                         "the import default and the interface default for "
                         "tlag_meth have diverged")

    @unittest.skipUnless((GUI / "src/ecprojectstate.h").exists(),
                         "eddyflow-gui not checked out beside this repository")
    def test_neither_route_switches_the_partition_on(self):
        self.assertEqual(self.engine_default("cec_meth"), "0")
        self.assertEqual(self.gui_default("cec_meth"), "0")

    @unittest.skipUnless((GUI / "src/ecproject.cpp").exists(),
                         "eddyflow-gui not checked out beside this repository")
    def test_the_interface_writes_both_back_out(self):
        """Reading them is not enough: the interface saves the converted
        project as a native file, so both keys have to be written to it or the
        next reader is back to guessing."""
        gui = read(self.GUI_PROJECT, GUI)
        #> The whole call, field and all. Matching the key name alone passes
        #> on a rename to INI_SCREEN_SETTINGS_77, because the old name is a
        #> prefix of the new one - and it says nothing about which field the
        #> value came from.
        for call in (
            "project_ini.setValue(EcIni::INI_PROJECT_72, "
            "ec_project_state_.projectGeneral.cec_meth)",
            "project_ini.setValue(EcIni::INI_SCREEN_SETTINGS_7, "
            "ec_project_state_.screenSetting.tlag_meth)",
        ):
            self.assertIn(call, gui, "no longer written on save: %s" % call)

        defs = read("src/ecinidefs.h", GUI)
        self.assertIn('INI_PROJECT_72   = QStringLiteral("cec_meth")', defs)
        self.assertIn('INI_SCREEN_SETTINGS_7    = QStringLiteral("tlag_meth")', defs)

    @unittest.skipUnless((GUI / "src/ecproject.cpp").exists(),
                         "eddyflow-gui not checked out beside this repository")
    def test_the_interface_import_is_a_rename_and_a_load(self):
        """It must not grow a defaults table of its own.

        If it ever needs one, the two would have to be kept in step by hand -
        which is the failure this whole class exists to prevent. Today it
        renames the fourth-gas keys and defers everything else to loadEcProject,
        where the state defaults above apply.
        """
        gui = read(self.GUI_PROJECT, GUI)
        body = gui[gui.index("bool EcProject::importEddyProProject("):]
        body = body[:body.index("\n}\n")]
        self.assertIn("fourthGasKeyRenames()", body)
        self.assertIn("loadEcProject(", body)
        for key in ("cec_meth", "tlag_meth"):
            self.assertNotIn(key, body,
                             "the interface's import now handles %s itself, so "
                             "it no longer agrees with the engine by "
                             "construction" % key)



if __name__ == "__main__":
    unittest.main()
