"""One question, one answer: which water corrects this gas?

`E2Col(gas)%moist_ref` is that answer. `ResolveGasRef` fills it with the
project's explicit choice, failing that the H2O on the gas's own analyser,
failing that the first H2O anywhere - and the interface mirrors the same three
rules, so what the user picked is what the engine uses.

Three places consume it, and they used to disagree:

  - `MoistTerms` took sigma and rho_w from it unconditionally;
  - the water-flux covariance in `TimeLagHandle` declined the pairing whenever
    gas and hygrometer sat on different analysers;
  - `PointByPointToMixingRatio` declined it too, for the same stated reason.

So a gas whose own analyser carries no hygrometer was corrected with the other
one's water for the mean terms and with nothing at all for the flux term and the
dilution. Not a conservative choice: the correction ran, half-built, from a water
the other two halves held it did not share. It is the fault `TimeLagHandle`'s own
comment records for a second hygrometer - "the two halves of one term disagreed
about which water they meant" - arriving by a different route.

Observed on a CO2/H2O/COS project with COS on a MIRO and the hygrometer on an
LI-7200: `FH2O_CELL_COS` read `-9999` in every period while `COS_MOIST_RHOW`
carried the LI-7200's water density, and `cos_flux` sat 1.4% off the version that
applied the term.

Borrowing across analysers stays a compromise, so it is announced -
ExceptionHandler(106), once per gas - rather than hidden either way.
"""

from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]

TIMELAG = "src/src_rp/timelag_handle.f90"
DILUTION = "src/src_common/point_by_point_to_mixing_ratio.f90"
FLUXES23 = "src/src_rp/fluxes23_rp.f90"
E2SET = "src/src_common/define_e2_set.f90"
HANDLER = "src/src_common/exception_handler.f90"


def read(path):
    return (ROOT / path).read_text(encoding="utf-8")


def code(path):
    """Source with comment-only lines dropped.

    Every one of these files documents the retired same-analyser test in prose
    right where it used to stand, so a naive search matches the explanation
    rather than a live gate.
    """
    return "\n".join(ln for ln in read(path).splitlines()
                     if not ln.lstrip().startswith("!"))


class EveryConsumerAsksTheRecord(unittest.TestCase):
    def test_all_three_read_moist_ref(self):
        for path in (TIMELAG, DILUTION, FLUXES23):
            self.assertIn(
                "E2Col(gas)%moist_ref" if path != TIMELAG
                else "E2Col(j)%moist_ref",
                code(path),
                "%s must take the water from the gas's own record" % path)

    def test_neither_consumer_second_guesses_it(self):
        """The same-model `cycle` is what made the three disagree."""
        for path in (TIMELAG, DILUTION):
            body = code(path)
            self.assertNotIn(
                "%instr%model /= E2Col(msl)%instr%model) cycle",
                body,
                "%s is declining the pairing again; moist_ref already "
                "prefers a hygrometer on the gas's own analyser, and "
                "refusing what it resolved leaves the correction half-built"
                % path)


class TheCovarianceUsesTheWatersOwnLag(unittest.TestCase):
    """A gas and a hygrometer sharing an analyser share a tube, so the gas's
    lag is the water's. Down a different tube it is not - which is the true
    content of the objection the old gate raised, and it is answered by
    choosing the lag rather than by dropping the term."""

    def setUp(self):
        src = code(TIMELAG)
        start = src.index("Stats%h2ocov_tl = error")
        self.body = src[start: src.index("end do", start)]

    def test_the_lag_is_selected_by_analyser(self):
        self.assertIn("E2Col(j)%instr%model == E2Col(msl)%instr%model", self.body)
        self.assertIn("lagRow = RowLags(j)", self.body)
        self.assertIn("lagRow = RowLags(msl)", self.body)

    def test_the_covariance_is_taken_at_that_lag(self):
        self.assertIn("lagRow, Stats%h2ocov_tl(j)", self.body)
        self.assertNotIn(
            "RowLags(j), Stats%h2ocov_tl(j)", self.body,
            "the gas's lag must no longer be hard-wired: on a borrowed "
            "hygrometer it is the wrong series' lag")

    def test_a_gas_is_still_not_its_own_water(self):
        """That covariance is the water flux itself, not a cross term."""
        self.assertIn("if (j == msl) cycle", self.body)


class TheDilutionUsesTheGasOwnCell(unittest.TestCase):
    """Removing the same-analyser gate makes this path reachable for a
    cross-analyser pair, so where it reads cell conditions starts to matter.

    It asked `tc` and `pi` - `firstCell` and `firstCell + 3`, the *first* cell
    block - and gated them on the *water's* analyser. Both were the same thing
    while one global cell served every instrument. With per-instrument blocks
    the first is whichever analyser happens to hold cell record one, and the
    gate asks an unrelated instrument, so a gas could be converted with another
    analyser's cell temperature and pressure.
    """

    def setUp(self):
        self.body = code(DILUTION)

    def test_the_cell_block_comes_from_the_gas(self):
        self.assertIn("cellBase = E2Col(gas)%cell_ref", self.body)
        self.assertIn("Va(:) = Ru * Set(:, cellBase) / Set(:, cellBase + 3)",
                      self.body)

    def test_an_unresolved_reference_declines_rather_than_borrows(self):
        """A gas whose analyser owns no cell record has no cell molar volume.

        This coerced an unresolved reference to `firstCell`, which is another
        instrument's cell on a site where only some analysers have one. The
        fallback lives in DefineE2Set now, and only where the cell records name
        no analyser at all - the pre-record case, where the single cell is the
        site's. Every reader of cell_ref declines together.
        """
        self.assertNotIn(
            "if (cellBase < firstCell .or. cellBase > lastCell) cellBase = firstCell",
            self.body,
            "borrowing the first block hands one analyser's cell to gases "
            "measured in another")
        self.assertIn("cellBase >= firstCell .and. cellBase <= lastCell", self.body)

    def test_the_bounds_test_is_nested_not_chained(self):
        """Fortran does not promise to stop evaluating a compound condition
        once it is decided, so `cellBase >= firstCell .and. E2Col(cellBase)%...`
        reaches E2Col(0) on an unresolved reference. It killed the run in the
        WPL step the first time this was written."""
        self.assertNotIn(".and. E2Col(cellBase)%present", self.body,
                         "chaining the bounds test onto the presence test "
                         "lets E2Col(0) be evaluated")
        self.assertIn("if (cellBase >= firstCell .and. cellBase <= lastCell) then",
                      self.body)

    def test_the_water_analyser_no_longer_gates_the_cell(self):
        self.assertNotIn("E2Col(tc)%instr%model == E2Col(msl)%instr%model",
                         self.body,
                         "the cell belongs to the gas, not to its hygrometer")


class TheCellTermsComeFromTheGasOwnBlock(unittest.TestCase):
    """Cell temperature and cell pressure both belong to the analyser that
    measured the gas, and both are read from that gas's own block.

    The pressure covariance already asked `cellPressureSlot(gas)`. The
    temperature one asked `tc` - `firstCell`, the *first* block - and then
    admitted a gas only if its analyser matched that block's. One global cell
    made those the same thing; with per-instrument blocks the first is
    whichever analyser happens to hold cell record one, so every gas on any
    other analyser failed the test and lost the cell-temperature term of its
    WPL correction outright, reported as `H_CELL = -9999`.

    Seen on a two-analyser site: the MIRO owned cell record one, so the
    LI-7200's CO2 and H2O were the ones going without. v7.2.5 does the same -
    its single cell record is the MIRO's - so this is not a regression, it is
    the defect becoming expressible now that cells are per-instrument.
    """

    def setUp(self):
        src = code(TIMELAG)
        start = src.index("Stats%tc_cov_tl = error")
        self.body = src[start: src.index("end do", start)]

    def test_the_block_comes_from_the_gas(self):
        self.assertIn("cellBase = E2Col(j)%cell_ref", self.body)
        self.assertIn("ColTC(1:nrow) = Set(1:nrow, cellBase)", self.body)

    def test_an_unresolved_reference_declines_rather_than_borrows(self):
        """No cell record for this gas's analyser means no cell-temperature
        covariance - not the first analyser's cell temperature."""
        self.assertIn("if (cellBase < firstCell .or. cellBase > lastCell) cycle",
                      self.body)
        self.assertNotIn("cellBase = firstCell", self.body)

    def test_no_instrument_gate_decides_who_gets_a_cell_term(self):
        self.assertNotIn("E2Col(tc)%instr%model", self.body,
                         "matching against the first block's analyser denies "
                         "the term to every gas on any other")
        self.assertNotIn("Set(1:nrow, tc)", self.body,
                         "tc is a fixed slot, not this gas's cell")

    def test_the_column_is_read_inside_the_loop(self):
        """It varies per gas now. Filling ColTC once before the loop would
        give every gas whichever cell the last iteration happened to load."""
        self.assertLess(self.body.index("do j = firstGas, lastGas"),
                        self.body.index("ColTC(1:nrow)"))


class BorrowingIsAnnounced(unittest.TestCase):
    """Silent is the one thing it must not be. Both halves used to decline
    without a word, and the run looked ordinary."""

    def test_only_a_gas_that_borrows_is_reported(self):
        body = code(E2SET)
        self.assertIn("call ExceptionHandler(106)", body)
        self.assertIn("E2Col(j)%instr%model == E2Col(msl)%instr%model) cycle",
                      body)
        self.assertIn("GasSlotIsWater(j)) cycle", body,
                      "a hygrometer does not borrow from itself")

    def test_it_is_raised_once_per_run_not_once_per_period(self):
        """DefineE2Set runs for every averaging period.

        Unguarded this printed ninety-seven times on a single day of CH-LAE -
        eleven lines each - which buries the log the warning exists to be read
        in. The flag has to be `save`d: a local would reset on every call and
        change nothing.
        """
        body = code(E2SET)
        self.assertIn("logical, save :: crossWaterWarned", body)
        self.assertIn("if (crossWaterWarned(j)) cycle", body)
        self.assertIn("crossWaterWarned(j) = .true.", body)

    def test_it_names_the_gas_and_both_analysers(self):
        body = code(E2SET)
        self.assertIn("trim(GasOutputLabel(j))", body)
        self.assertIn("trim(E2Col(msl)%instr%model)", body)

    def test_the_handler_explains_the_compromise(self):
        text = read(HANDLER)
        self.assertIn("case(106)", text)
        block = text[text.index("case(106)"):]
        block = block[: block.index("end select")]
        self.assertIn("time lag", block)
        self.assertIn("own analyser", block,
                      "the message has to say how to remove the compromise")


if __name__ == "__main__":
    unittest.main()
