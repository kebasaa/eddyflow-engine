"""A gas slot may be named for its position, never for a species.

`co2 = 5, h2o = 6, ch4 = 7, n2o = 8, gas4 = 8` are gone, and so is the
`histCO2`/`histH2O`/`histCH4` spelling that replaced them. They read as species
and meant positions, and every defect this migration fixed began with that
reading: `E2Col(co2)` looks like "the carbon dioxide column" and is "the fifth
variable", which is carbon dioxide only when a project happens to declare it
first. The same four numbers survive as `histGas1..histGas4`, offsets from
`firstGas`, so a use says which position it means and nothing about species.

Where they legitimately remain:

  compatibility mode    fix_out_format promises the fixed EddyPro 7.x column
                        set, which is four gas blocks whatever the project
                        holds. Widening it would break the compatibility the
                        flag exists to provide.

  on-disk aliases       LegacySpectralVarTag, HistoricGasSlot,
                        TimelagOptGasLabel and GasSlotFromDynMDTag accept the
                        four historical spellings so files written before the
                        records keep being readable. These only ever *widen*
                        what a reader accepts.

  the interface's own   read_ini_fcc reads three month-grouping tables the
                        tables                interface exposes, for CO2, CH4
                        and the fourth gas; every other gas inherits CO2's.

  fallbacks             PrimaryWaterOutSlot and PrimaryCarbonOutSlot fall back
                        to the historical slot when a project describes no
                        water or no CO2 - which is what the gates they replace
                        evaluated in that case.

  dormant code          The LI-7550 analog-filter and analytic transfer
                        function paths are gated on flags hard-set false.
                        Pinned separately by test_li7550_dormant_static.

The flat ini layer used to be on this list. It is gone: read_ini_rp and
read_ini_fcc read every per-gas setting from its record, and a project without
records is refused rather than half-processed. That alone took read_ini_rp
from 76 occurrences to none.

What must not happen is a *new* `Set(:, histGas1)` or `E2Col(histGas2)` in a
processing path. So this pins the count per file: a number going up fails, a
number going down is the work continuing and the expectation should be lowered
to match. Deliberately a count and not a pattern - the legitimate uses are too
varied to describe by regex, and a whitelist of exact strings would have to be
rewritten every time a comment moved.
"""

from pathlib import Path
import re
import subprocess
import unittest


ROOT = Path(__file__).resolve().parents[1]

SLOT_NAMES = re.compile(r"\b(histGas1|histGas2|histGas3|histGas4)\b")

#: The spellings that must never come back. A bare `co2` as an integer
#: parameter is the thing this whole effort removed; as a *species string* it
#: is correct and common, so the check below strips string literals.
RETIRED_SPELLINGS = ("co2", "h2o", "ch4", "n2o", "gas4",
                     "w_co2", "w_h2o", "w_ch4", "w_gas4",
                     "histCO2", "histH2O", "histCH4", "histN2O")

#: Occurrences outside comments, per file. Generated from the tree; see the
#: module docstring for why each file still has any.
ALLOWED = {
    "src/src_common/define_all_var_set.f90": 6,
    "src/src_common/define_used_variables.f90": 1,
    "src/src_common/gas_slot_resolution.f90": 9,
    "src/src_common/m_common_global_var.f90": 4,
    "src/src_common/m_typedef.f90": 4,
    "src/src_fcc/cospectra_sorting_and_averaging.f90": 1,
    "src/src_fcc/init_out_files.f90": 1,
    "src/src_fcc/read_ini_fcc.f90": 5,
    "src/src_fcc/spectral_assessment_diagnostics.f90": 3,
    "src/src_rp/init_outfiles_rp.f90": 1,
    "src/src_rp/pwb_timelag_handle.f90": 4,
}


def tracked_sources():
    out = subprocess.run(["git", "ls-files", "src"], cwd=ROOT,
                         capture_output=True, text=True, check=True).stdout
    return [p for p in out.split()
            if p.endswith(".f90") or p.endswith(".inc")]


def slot_uses(rel):
    """Occurrences outside full-line comments."""
    text = (ROOT / rel).read_text(encoding="utf-8", errors="replace")
    return sum(len(SLOT_NAMES.findall(ln)) for ln in text.splitlines()
               if not ln.lstrip().startswith("!"))


def code_without_strings(rel):
    """Source with comments and string literals removed.

    A species name in quotes - `case ('co2')`, or a header literal listing
    `co2, ch4, 4th gas` - is data, not a slot. Only a bare identifier is.

    Quote state is carried across lines, because Fortran continues a string
    with a trailing `&` and out_raw_data opens one on one line and closes it
    two lines later. Stripping line by line left that text looking like code.
    """
    out, in_str = [], False
    for ln in (ROOT / rel).read_text(encoding="utf-8", errors="replace").splitlines():
        if not in_str and ln.lstrip().startswith("!"):
            continue
        kept = []
        for ch in (ln if in_str else ln.split("!")[0]):
            if ch == "'":
                in_str = not in_str
                continue
            if not in_str:
                kept.append(ch)
        out.append("".join(kept))
    return "\n".join(out)


class LegacySlotUseIsPinned(unittest.TestCase):
    def test_no_file_gains_a_historical_slot_reference(self):
        grew = []
        for rel in tracked_sources():
            n = slot_uses(rel)
            allowed = ALLOWED.get(rel, 0)
            if n > allowed:
                grew.append("%s: %d, expected at most %d" % (rel, n, allowed))
        self.assertFalse(
            grew,
            "these files gained a fixed gas-slot reference:\n  "
            + "\n  ".join(grew)
            + "\n\nA gas is addressed by its record. If the new use is one of "
              "the legitimate kinds - a compatibility mode, an on-disk alias, "
              "an interface table or a no-water fallback - say which in a "
              "comment and raise the count here.")

    def test_the_expectations_are_not_stale(self):
        """A count that has dropped is the work continuing; lower it, so the
        file cannot silently regain what it gave up."""
        stale = []
        for rel, allowed in sorted(ALLOWED.items()):
            if not (ROOT / rel).is_file():
                stale.append("%s: listed but does not exist" % rel)
                continue
            n = slot_uses(rel)
            if n < allowed:
                stale.append("%s: %d now, expectation still %d" % (rel, n, allowed))
        self.assertFalse(
            stale,
            "lower these expectations to match the tree:\n  " + "\n  ".join(stale))


class TheSpeciesSpellingsAreGone(unittest.TestCase):
    """`co2` as an integer is what this migration removed.

    Keeping the numbers is fine - the fixed output format and the on-disk
    readers need them. Keeping the *names* is not: they made a position read
    as a measurement, and that is what put water's molecular weight on a trace
    gas, gave methane carbon dioxide's flux column and left five subsystems
    silently correcting the wrong species.
    """

    def test_no_source_declares_or_uses_the_old_names(self):
        offenders = []
        for rel in tracked_sources():
            body = code_without_strings(rel)
            for name in RETIRED_SPELLINGS:
                if re.search(r"\b%s\b" % re.escape(name), body):
                    offenders.append("%s: %s" % (rel, name))
        self.assertFalse(
            offenders,
            "these are slot constants named for a species:\n  "
            + "\n  ".join(offenders)
            + "\n\nUse histGas1..histGas4 if one of the four historical "
              "positions is genuinely meant, or firstGas..lastGas and the gas "
              "record if a species is.")

    def test_the_historical_slots_are_positional(self):
        """Expressed as offsets from firstGas, not as numbers.

        Writing 5, 6, 7, 8 would be the same values and would drift the moment
        the anemometric block changed width; deriving them says these are the
        first four gas slots, which is all they ever were.
        """
        src = (ROOT / "src/src_common/m_typedef.f90").read_text(
            encoding="utf-8", errors="replace")
        self.assertIn("integer, parameter :: nHistoricGasSlots = 4", src)
        self.assertIn("integer, parameter :: histGas1 = firstGas", src)
        for n, name in enumerate(("histGas2", "histGas3", "histGas4"), start=1):
            self.assertIn(
                "integer, parameter :: %s = firstGas + %d" % (name, n), src,
                "%s must be an offset from firstGas, not a literal" % name)

    def test_the_dynamic_bounds_are_declared_from_the_capacity(self):
        """firstGas/lastGas are what everything else should iterate."""
        src = (ROOT / "src/src_common/m_typedef.f90").read_text(
            encoding="utf-8", errors="replace")
        self.assertIn("firstGas = NumAnemVar + 1", src)
        self.assertIn("lastGas  = GHGNumVar", src)


if __name__ == "__main__":
    unittest.main()
