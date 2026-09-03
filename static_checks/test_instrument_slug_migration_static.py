"""A model key is resolved to its current spelling before anything reads it.

The Campbell keys were renamed twice: the bare `csat3` and `csat3b` an EddyPro
site wrote took a manufacturer-name prefix for one release, and now carry
`csi_`. No other manufacturer moved - every Gill, Metek, Young and LI-COR key
is what EddyPro wrote - so this only ever bites on a Campbell site.

There it bites hard, and quietly. Every gate downstream is a `select case` over
`csi_*` names: metadata validation, the sonic coordinate adjustment, the
spectral transfer functions, the master-sonic overrides. A retired key matches
none of them, so ValidateMetadata takes its `case default`, passed(4) goes
false, and ExceptionHandler(25) skips the file. For a GHG archive - which this
program unzips itself, and which the interface cannot reach or rewrite - that
is every file in the run, reported only as a metadata problem per file.

CanonicalInstrumentModel is where that is caught, and it is caught at
ingestion so that every one of those select cases keeps its `csi_*`-only list.
The interface carries the same table for the metadata it reads; nothing at
build time connects the two repositories.

The middle spelling - the one release that prefixed these with the
manufacturer's name - is deliberately not handled here. That prefix is a
retired identifier in this program and
test_irga_path_classification_static.test_legacy_model_prefixes_are_not_active_identifiers
bans the string from the tree outright.
"""

from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]

TYPEDEF = "src/src_common/m_typedef.f90"
READ_MD = "src/src_common/read_metadata_file.f90"
DYN_MD = "src/src_rp/retrieve_dynamic_metadata.f90"
VALIDATION = "src/src_common/metadata_file_validation.f90"

#: Spellings a file may carry, and what each one means now.
BARE_MODEL_KEYS = {
    "csat3": "csi_csat3",
    "csat3a": "csi_csat3a",
    "csat3b": "csi_csat3b",
    "ec150": "csi_ec150",
}


def read(rel):
    return (ROOT / rel).read_text(encoding="utf-8", errors="replace")


def body(src, start, end):
    text = src[src.index(start):]
    return text[: text.index(end)]


class TheHelper(unittest.TestCase):
    def setUp(self):
        self.body = body(read(TYPEDEF),
                         "character(32) function CanonicalInstrumentModel(",
                         "end function CanonicalInstrumentModel")

    def test_every_bare_spelling_is_mapped(self):
        for bare, canonical in BARE_MODEL_KEYS.items():
            self.assertIn("'%s'" % bare, self.body)
            self.assertIn("'%s'" % canonical, self.body)

    def test_no_bare_key_is_invented(self):
        # csat3c, ec155 and tga200a arrived already carrying a prefix, so no
        # file can spell them without one. A bare case for them would be
        # matching something nothing writes.
        for never_bare in ("csat3c", "ec155", "tga200a"):
            self.assertNotIn("case ('%s')" % never_bare, self.body)

    def test_the_index_suffix_is_preserved(self):
        # The value is the key plus a trailing `_<n>` that pairs it with an
        # instrument, and every consumer strips it with the same two-character
        # arithmetic InstrumentModelBase uses. Dropping it here would leave
        # each of them chopping two characters off the key instead.
        self.assertIn("model(model_len - 1:model_len)", self.body)

    def test_it_splits_with_the_shared_helper(self):
        self.assertIn("InstrumentModelBase(model)", self.body)

    def test_it_returns_the_input_untouched_by_default(self):
        # Anything not Campbell must come back byte-identical, or a Gill or
        # LI-COR site would be rewritten on a rule that was never about it.
        self.assertIn("CanonicalInstrumentModel = model", self.body)
        self.assertIn("case default", self.body)


class EveryIngestionPoint(unittest.TestCase):
    """The three places a model key enters the program from a file."""

    def test_the_metadata_reader_normalises(self):
        # ReadMetadataFile is the single reader for both the standalone
        # metadata file and the copy inside a GHG archive, so these calls
        # cover both.
        #
        # TWO model keys enter here, not one: instr_<k>_model, and the
        # extended-.ghg instr_<k>_ef_model that overrides it where the first
        # states a generic stand-in. Both are file input and both must be
        # normalised - an ef_model left unnormalised would be the only model
        # in the program still wearing its raw spelling, and it is the one
        # every downstream select case then matches on.
        self.assertEqual(2, read(READ_MD).count("CanonicalInstrumentModel("))

    def test_the_dynamic_metadata_normalises_sonic_and_analysers(self):
        self.assertEqual(2, read(DYN_MD).count("CanonicalInstrumentModel("))

    def test_the_dynamic_metadata_pass_follows_the_field_loop(self):
        # DynMDGasOrder visits fields in file-column order, so the model is
        # not guaranteed to have been read at any given point inside the loop.
        src = read(DYN_MD)
        self.assertLess(src.index("select case (fld)"),
                        src.index("CanonicalInstrumentModel("))


class TheGatesDownstreamStayNarrow(unittest.TestCase):
    """Normalising at ingestion is what lets these keep one spelling each.

    If a legacy key ever reached them the fix would be to widen every list,
    in every file, forever. This asserts nobody has started doing that.
    """

    def test_validation_lists_only_current_keys(self):
        src = read(VALIDATION)
        for bare in BARE_MODEL_KEYS:
            self.assertNotIn(
                "'%s'" % bare, src,
                "%s reached validation - it should have been resolved at "
                "ingestion" % bare)


class ColumnsAreMatchedOnTheFilesOwnSpelling(unittest.TestCase):
    """Canonicalising the instrument does not canonicalise the columns.

    col_<n>_instrument names an instrument by its MODEL STRING, and nothing
    rewrites it. So the moment CanonicalInstrumentModel changes what the
    instrument calls itself - which is exactly what it is for - the columns are
    left pointing at a name the instrument list no longer holds.

    An archive carrying EddyPro's own `csat3b_1` (the only Campbell spelling
    EddyPro takes) had its instrument turned into `csi_csat3b_1` while all five
    sonic columns still said `csat3b_1`. No column bound to the sonic, and the
    run died with "exactly one selected u, v, w and one selected ts or sos are
    required" - which does not mention instruments at all.

    Instr%ep_label carries the file's own spelling so the columns can still be
    matched on it, and it must be read RAW - assigning it from %model after
    canonicalisation makes it a second copy of the canonical name and restores
    the bug exactly.
    """

    def setUp(self):
        self.src = read(READ_MD)

    def test_ep_label_is_read_raw_not_copied_from_model(self):
        self.assertNotIn("Instr(i)%ep_label = Instr(i)%model", self.src,
                         "ep_label must hold the FILE's spelling, not the "
                         "canonicalised one")
        self.assertIn("Instr(i)%ep_label = ACTags(init_ac_instr "
                      "+ i*leap_ac_instr + 2)%value", self.src)

    def test_the_column_loop_tries_ep_label_too(self):
        self.assertIn("Instr(j)%ep_label(1:len_trim(Instr(j)%ep_label))",
                      self.src)

    def test_the_ep_label_match_is_length_guarded(self):
        #> index(x, '') is 1, so an unguarded empty label would match the first
        #> column against the first instrument.
        self.assertIn("if (len_trim(Instr(j)%ep_label) > 0) then", self.src)

    def test_both_spellings_reach_the_same_offset(self):
        #> ep_label and %model must read the SAME tag - offset 2 of the
        #> instrument block - or they describe different instruments.
        #>
        #> Four rather than two: a Fortran substring read names its tag once
        #> for the value and once inside the len_trim that bounds it, so each
        #> of the two reads contributes two occurrences.
        self.assertEqual(
            4, self.src.count("ACTags(init_ac_instr + i*leap_ac_instr + 2)%value"),
            "the model tag is read once for %model and once for %ep_label")


if __name__ == "__main__":
    unittest.main()
