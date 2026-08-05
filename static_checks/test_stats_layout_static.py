"""The statistics files' header and rows walk one shared layout.

st1..st7 write five per-slot families - mean, var, st_dev, skw, kur - and name
their variables once in a header. The two agreed only by coincidence: the row
writer looped `u, pe`, and while E2NumVar was 14 that enumerated exactly the
twelve names the header listed.

E2NumVar is now 102 - 64 gas slots and 32 per-instrument cell slots - so the
rows became roughly seven times wider than their own header. Measured before
the fix, on both the 4-gas and the 8-gas fixture: header 88 fields, rows 528.
The files could not be read by column at all, and had not been since the
capacity change.

Nothing caught it because no fixture enabled them. `base_n_gas_st.eddyflow`
does now, and `run.sh` on it must give header == rows in every st file.

Both sides call StatsLayoutSlots. What this check defends is that they keep
doing so: a family added to one side and not the other reintroduces exactly
the same silent mismatch, and it is invisible unless something compares the
two widths.
"""

from pathlib import Path
import re
import unittest


ROOT = Path(__file__).resolve().parents[1]

WRITER = "src/src_common/write_out_stats.f90"
HEADER = "src/src_rp/init_outfiles_rp.f90"
HELPER = "src/src_common/gas_slot_resolution.f90"
TYPEDEF = "src/src_common/m_typedef.f90"
FIXTURE = "tests/regression/base_n_gas_st.eddyflow"


def code(path):
    return "\n".join(ln for ln in (ROOT / path).read_text(encoding="utf-8",
                                                          errors="replace").splitlines()
                     if not ln.lstrip().startswith("!"))


class OneLayoutServesBothSides(unittest.TestCase):
    def test_the_helper_exists_and_is_record_driven(self):
        src = code(HELPER)
        self.assertIn("subroutine StatsLayoutSlots", src)
        body = src.split("subroutine StatsLayoutSlots")[1].split("end subroutine")[0]
        self.assertIn("EddyFlowProj%gas_num", body,
                      "the gas entries must come from the records, so a gas "
                      "named without a column keeps its column of error codes")
        self.assertRegex(body, r"do gas = firstGas, lastGas")

    def test_both_sides_call_it(self):
        for path in (WRITER, HEADER):
            self.assertIn("call StatsLayoutSlots", code(path),
                          "%s must walk the shared layout, not its own idea "
                          "of which slots exist" % path)

    def test_the_row_writer_no_longer_walks_every_slot(self):
        """`do j = u, pe` is 102 slots now, not the 14 it was written for."""
        self.assertNotRegex(
            code(WRITER), r"do j = u, pe",
            "the per-slot families must iterate the layout list; `u, pe` "
            "enumerates every cell slot and outruns the header")

    def test_no_four_gas_covariance_arms_remain(self):
        """The three covariance groups had one arm per historical gas, and the
        header had the same four - so both stop at the fourth or neither."""
        src = code(WRITER)
        for token in ("(co2)", "(h2o)", "(ch4)", "(gas4)"):
            self.assertNotIn(token, src,
                             "%s names a fixed gas slot in the stats writer" % token)

    def test_the_header_is_built_once_and_shared(self):
        """It was seven byte-identical literals; seven copies is seven chances
        for one to drift away from the writer."""
        src = code(HEADER)
        self.assertIn("stats_header", src)
        self.assertEqual(
            len(re.findall(r"write\(ust\d, '\(a\)'\) stats_header", src)), 7,
            "all seven statistics files must write the same generated header")


class TheStatisticsTypeIsIndexedBySlot(unittest.TestCase):
    """StatsType carried seven named scalars beside its E2NumVar-wide arrays:
    h2ocov_tl_co2/ch4/gas4 and tc_cov_tl_co2/h2o/ch4/gas4. Those are the water
    and cell-temperature covariances taken at *another* gas's timelag, and
    they feed the internal sensible heat flux and the in-cell
    evapotranspiration. Named per position, they stopped at the fourth slot,
    so a fifth gas got neither however it was configured - and the water arm
    had to read the field named for the historical slot while writing its
    result to the resolved one."""

    def test_the_covariances_at_other_timelags_are_arrays(self):
        src = code(TYPEDEF)
        for token in ("h2ocov_tl_", "tc_cov_tl_"):
            self.assertNotIn(token, src,
                             "%s names a fixed gas slot; these are per-slot "
                             "quantities and must be indexed" % token)
        self.assertIn("h2ocov_tl(E2NumVar)", src)
        self.assertIn("tc_cov_tl(E2NumVar)", src)

    def test_the_water_column_they_are_taken_against_is_resolved(self):
        """The covariance is against the water that corrects THIS gas.

        Reading E2Col(h2o) took whatever species record two held - on
        base_h2o_late that was N2O, and the in-cell water flux came out three
        orders of magnitude wrong. Resolving the site's water fixed that, and
        was still one answer for every gas: a gas on a second analyser was
        matched against the primary hygrometer's model, failed, and got no
        covariance at all. It reads the gas's own moist_ref now.
        """
        src = code("src/src_rp/timelag_handle.f90")
        self.assertIn("msl = E2Col(j)%moist_ref", src)
        self.assertNotIn("E2Col(h2o)", src)
        self.assertNotIn("E2Col(h2o)", src,
                         "the water column must be resolved, not indexed")


class AFixtureTurnsThemOn(unittest.TestCase):
    def test_a_fixture_enables_the_statistics_files(self):
        """The defect survived because every fixture left out_st_* at 0, so
        the files were never written and never compared."""
        path = ROOT / FIXTURE
        self.assertTrue(path.exists(),
                        "base_n_gas_st.eddyflow must exist so the statistics "
                        "files are exercised by the harness")
        text = path.read_text(encoding="utf-8", errors="replace")
        self.assertRegex(text, r"(?m)^out_st_1=1$",
                         "the fixture must switch a statistics file on")


if __name__ == "__main__":
    unittest.main()
