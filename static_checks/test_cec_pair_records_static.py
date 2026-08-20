"""The pairing keys are appended, and nothing already written moves.

A Conditional Eddy Covariance pairing is one carbon channel, one water channel
and any further species partitioned in the octants those two define. The
project states them as `cec_num` plus `cec_<i>_meth` / `_co2` / `_h2o` /
`_extra`, and those keys had to go somewhere in a table that is addressed by
position.

That is the hazard this pins. `EPPrjNTags` and `EPPrjCTags` are positional -
a tag's index *is* its identity - and the per-gas, per-cell and per-diagnostic
record blocks are read by offset from an origin. Put a new key in the wrong
place and the origin moves, every one of the 1600 per-gas settings is re-emitted
at a new index, and a project written yesterday silently means something else
today. The generator refuses to write if that happens; this says the same thing
from outside it, so the refusal cannot be argued away by editing the generator.

`cec_singular_band` is the other case: a plain [Project] scalar, which cannot be
appended at all because the record blocks come after everything. It goes into
one of the blanks a retired key left, below the first origin.
"""

from pathlib import Path
import re
import subprocess
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
TAGS = (ROOT / "src/src_common/m_common_global_var.f90").read_text(
    encoding="utf-8", errors="replace")


def origin(name):
    m = re.search(r"integer, parameter :: %s\s*=\s*(\d+)" % name, TAGS)
    assert m is not None, "%s not found" % name
    return int(m.group(1))


class TheExistingOriginsDidNotMove(unittest.TestCase):
    """The whole point. These are the numbers every per-gas setting is read
    at; if one of them changes, files already on disk are misread."""

    def test_the_record_origins_are_where_they_were(self):
        self.assertEqual(origin("gasNumTag"), 33)
        self.assertEqual(origin("cellNumTag"), 34)
        self.assertEqual(origin("diagNumTag"), 35)
        self.assertEqual(origin("gasRecOriginN"), 36)
        self.assertEqual(origin("cellRecOriginN"), 420)
        self.assertEqual(origin("diagRecOriginN"), 452)
        self.assertEqual(origin("gasRecOriginC"), 51)
        self.assertEqual(origin("cellRecOriginC"), 179)
        self.assertEqual(origin("diagRecOriginC"), 243)
        self.assertEqual(origin("rpGasOriginN"), 425)
        self.assertEqual(origin("rpGasOriginC"), 102)
        self.assertEqual(origin("fccGasOriginN"), 110)
        self.assertEqual(origin("fccGasOriginC"), 28)

    def test_the_cec_block_comes_after_every_other_record_block(self):
        self.assertGreater(origin("cecNumTag"), origin("diagRecOriginN"))
        self.assertEqual(origin("cecRecOriginN"), origin("cecNumTag") + 1)
        self.assertGreater(origin("cecRecOriginC"), origin("diagRecOriginC"))

    def test_the_strides_are_named_not_counted_at_the_use_site(self):
        self.assertEqual(origin("cecRecLeapN"), 3)
        self.assertEqual(origin("cecRecLeapC"), 1)


class TheKeysAreRegistered(unittest.TestCase):
    def test_the_scalar_sits_in_a_blank_below_the_first_origin(self):
        m = re.search(r"EPPrjNTags\((\d+)\)%Label / 'cec_singular_band' /", TAGS)
        self.assertIsNotNone(m, "cec_singular_band is not in the tag table")
        self.assertLess(int(m.group(1)), origin("gasNumTag"),
                        "a [Project] scalar appended past the first record "
                        "origin lifts every per-gas setting after it")

    def test_the_pairing_keys_are_registered(self):
        self.assertIn("'cec_num'", TAGS)
        for field in ("meth", "co2", "h2o"):
            self.assertIn("'cec_1_%s'" % field, TAGS)
        self.assertIn("'cec_1_extra'", TAGS)

    def test_the_capacity_matches_the_engines_own_bound(self):
        typedef = (ROOT / "src/src_common/m_typedef.f90").read_text(
            encoding="utf-8", errors="replace")
        m = re.search(r"integer, parameter :: MaxNumCecPairs = (\d+)", typedef)
        self.assertIsNotNone(m)
        cap = int(m.group(1))
        self.assertIn("'cec_%d_meth'" % cap, TAGS)
        self.assertNotIn("'cec_%d_meth'" % (cap + 1), TAGS)


class TheGeneratorOwnsTheTable(unittest.TestCase):
    def test_the_block_is_what_the_generator_would_write(self):
        """Hand-editing the table is how it drifts from the generator, and the
        generator is what refuses to move an origin."""
        result = subprocess.run(
            [sys.executable, str(ROOT / "prj/gen_project_tags.py"), "--check"],
            cwd=ROOT, capture_output=True, text=True)
        self.assertEqual(result.returncode, 0,
                         "gen_project_tags.py --check failed:\n"
                         + result.stdout + result.stderr)


class ThePairingsAreReadBesideTheGasRecords(unittest.TestCase):
    def setUp(self):
        self.parser = (ROOT / "src/src_common/write_processing_project_variables.f90"
                       ).read_text(encoding="utf-8", errors="replace")

    def test_they_are_read_after_the_gas_loop(self):
        #> A pairing names gas records, so there is no point resolving one
        #> against a list that is not read yet.
        self.assertIn("call ReadCecRecords()", self.parser)
        self.assertLess(self.parser.index("EddyFlowProj%gas(i)%col"),
                        self.parser.index("call ReadCecRecords()"))

    def test_an_absent_count_means_the_layout_decides(self):
        self.assertIn("EddyFlowProj%cec_num = 0", self.parser)
        resolver = (ROOT / "src/src_common/gas_slot_resolution.f90").read_text(
            encoding="utf-8", errors="replace")
        self.assertIn("call AutoCecPairs(pairs, npairs)", resolver)
        #> One pairing per carbon channel, each with the water on its own
        #> analyser - the same rule the interface seeds into the file.
        self.assertIn("CecWaterOnAnalyserOf(gas)", resolver)

    def test_the_count_is_clamped_to_what_can_be_held(self):
        self.assertIn("min(max(EddyFlowProj%cec_num, 0), MaxNumCecPairs)",
                      self.parser)


if __name__ == "__main__":
    unittest.main()
