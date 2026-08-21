"""Buffers whose width follows the gas count are sized from the gas count.

Three files carry a block per gas - the raw-dataset header, the dynamic
metadata file and the calibration-events file - so their width grows with the
number of gases a project declares. All three were held in buffers whose size
was written as a literal, chosen when four gases was the ceiling.

`raw_out_header` is the one that proves the point. The loop filling it was
widened from the four historical slots to MaxNumGases; the `character(512)` it
writes into was not. Twenty-five characters a name against sixty-four gases is
1600, so it overflows at about seventeen. Nothing said so, because a literal
in a declaration has no connection to a loop bound in another file.

The other two were worse in kind if not in reach: `mdStringVars(256)` and
`text_vars(128)` were filled by `var_num = var_num + 1` with no `size()` test
at all, so a wide file wrote past the end of a stack array rather than
truncating. `Headerlabels(256)` did stop, but silently, so the fields past the
cut simply never existed as far as the rest of the run was concerned.

None of this is reachable with the fixtures available - the metadata here does
not describe enough gas columns - so it is pinned by construction rather than
observed. The arithmetic is not in doubt.
"""

from pathlib import Path
import re
import unittest


ROOT = Path(__file__).resolve().parents[1]


def code(rel):
    return "\n".join(ln for ln in (ROOT / rel).read_text(
        encoding="utf-8", errors="replace").splitlines()
        if not ln.lstrip().startswith("!"))


def const(name):
    src = (ROOT / "src/src_common/m_typedef.f90").read_text(
        encoding="utf-8", errors="replace")
    m = re.search(r"integer, parameter :: %s = (\d+)" % name, src)
    return int(m.group(1)) if m else None


class TheSizesAreDerived(unittest.TestCase):
    def test_the_capacity_constants_exist(self):
        src = code("src/src_common/m_typedef.f90")
        self.assertIn("integer, parameter :: MaxRowFields", src)
        self.assertIn("integer, parameter :: RawHeaderLen", src)

    def test_they_are_expressed_from_the_gas_capacity(self):
        """A literal here is the bug, whatever number it holds."""
        src = code("src/src_common/m_typedef.f90")
        row = src[src.index("MaxRowFields"):]
        row = row[:row.index("RawHeaderLen")]
        self.assertIn("MaxNumGases", row)
        self.assertIn("nDynMDGasFields", row)

        raw = src[src.index("integer, parameter :: RawHeaderLen"):]
        raw = raw[:raw.index("\n", raw.index("="))+200].split("\n")[0:3]
        self.assertIn("E2NumVar", "".join(raw),
                      "the raw header must be sized from the variable set it "
                      "names, not from a literal")

    def test_the_old_literal_could_not_hold_the_gases_alone(self):
        """The arithmetic that makes this a bug rather than a tidy-up.

        Twenty-five characters a name against MaxNumGases gases already
        exceeds the 512 the buffer used to be, before counting the
        anemometric four or the cell and ambient columns.
        """
        gases = const("MaxNumGases")
        self.assertIsNotNone(gases)
        self.assertGreater(
            25 * gases, 512,
            "if this ever stops being true the overflow story needs redoing")


class TheBuffersUseThoseSizes(unittest.TestCase):
    DECLS = {
        "src/src_rp/m_rp_global_var.f90":
            ["character(RawHeaderLen) :: raw_out_header",
             "integer :: DynamicMetadataOrder(MaxRowFields)"],
        "src/src_rp/retrieve_dynamic_metadata.f90":
            ["character(32) :: mdStringVars(MaxRowFields)",
             "character(32) :: mdCurrentStringVars(MaxRowFields)"],
        "src/src_rp/drift_retrieve_calibration_events.f90":
            ["character(32) :: text_vars(MaxRowFields)"],
        "src/src_rp/init_dynamic_medata.f90":
            ["character(64) :: Headerlabels(MaxRowFields)"],
    }

    def test_every_row_buffer_is_declared_from_a_constant(self):
        for rel, decls in self.DECLS.items():
            src = code(rel)
            for decl in decls:
                self.assertIn(decl, src,
                              "%s: %r is not sized from a capacity constant"
                              % (rel, decl))

    #: The exact declarations these replaced. Matched by name rather than by
    #: searching for the numbers, because 128 and 256 are perfectly ordinary
    #: sizes for the many arrays here that do not scale with the gas count.
    RETIRED = {
        "src/src_rp/m_rp_global_var.f90":
            ["character(512) :: raw_out_header",
             "DynamicMetadataOrder(256)"],
        "src/src_rp/retrieve_dynamic_metadata.f90":
            ["mdStringVars(256)", "mdCurrentStringVars(256)"],
        "src/src_rp/drift_retrieve_calibration_events.f90":
            ["text_vars(128)"],
        "src/src_rp/init_dynamic_medata.f90":
            ["Headerlabels(256)"],
    }

    def test_no_literal_sized_row_buffer_survives(self):
        for rel, gone in self.RETIRED.items():
            src = code(rel)
            for decl in gone:
                self.assertNotIn(
                    decl, src,
                    "%s: %r is back. Its width follows the gas count, so a "
                    "literal will fall behind the next time a loop is widened "
                    "- which is exactly what happened to raw_out_header."
                    % (rel, decl))


class EveryUnboundedFillIsGuarded(unittest.TestCase):
    """`var_num = var_num + 1` with no size test is how a row runs off the end."""

    FILLS = {
        "src/src_rp/retrieve_dynamic_metadata.f90": "mdStringVars",
        "src/src_rp/drift_retrieve_calibration_events.f90": "text_vars",
        "src/src_rp/init_dynamic_medata.f90": "Headerlabels",
    }

    def test_the_fill_tests_the_size_first(self):
        for rel, name in self.FILLS.items():
            src = code(rel)
            self.assertIn("size(%s)" % name, src,
                          "%s fills %s with no size() test" % (rel, name))

    def test_an_overrun_is_reported_not_silent(self):
        """Headerlabels used to stop without saying anything, so a truncated
        header read as a complete one."""
        for rel in self.FILLS:
            self.assertIn("ExceptionHandler(105)", code(rel),
                          "%s truncates without telling the user" % rel)


if __name__ == "__main__":
    unittest.main()
