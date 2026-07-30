"""Random uncertainty and the integral turbulence scale run over every gas.

Both files were bounded at the fourth gas slot and neither referenced
firstGas/lastGas at all. Because `Essentials` is never wholesale reset, the
slots past the fourth kept whatever was in them - so an 8-gas run reported
RANDUNC_HF = 0.00000 for gases 5+, a fabricated number where the sentinel
belongs, while the first four correctly read -9999.

Every loop here already guards on E2Col(var)%present with an `else` that
writes `error`, which is why widening the bound is safe: an unconfigured slot
reports "not computed" rather than being consulted at a zero default. That
guard is the thing to protect, so it is asserted alongside the bound.
"""

from pathlib import Path
import re
import unittest


ROOT = Path(__file__).resolve().parents[1]

FILES = (
    "src/src_rp/random_error_handle.f90",
    "src/src_rp/integral_turbulence_scale.f90",
)


def read(path):
    return (ROOT / path).read_text(encoding="utf-8")


class RandomUncertaintyCoversEveryGas(unittest.TestCase):
    def test_no_loop_is_bounded_at_the_fourth_gas(self):
        pattern = re.compile(r"\bu\s*[,:]\s*gas4\b")
        for path in FILES:
            hit = pattern.search(read(path))
            self.assertIsNone(
                hit,
                f"{path} still bounds a variable loop at gas4, so gases 5+ "
                f"keep whatever Essentials held from a previous period",
            )

    def test_loops_run_to_lastgas(self):
        for path in FILES:
            source = read(path)
            self.assertIn(
                "u, lastGas", source,
                f"{path} no longer iterates the full gas range")

    def test_absent_slots_report_not_computed(self):
        """The guard that makes widening safe.

        Without it, widening promotes every unconfigured slot from "never
        consulted" to "consulted at its zero default" - the same failure that
        wiped a whole gas's series when the absolute-limits loop was widened.
        """
        ru = read("src/src_rp/random_error_handle.f90")
        self.assertIn("E2Col(var)%present", ru)
        self.assertIn("Essentials%rand_uncer(var) = error", ru)
        self.assertIn("E2Col(gas_var)%present", ru)
        self.assertIn("Essentials%rand_uncer(gas_var) = error", ru)

        its = read("src/src_rp/integral_turbulence_scale.f90")
        self.assertIn("E2Col(var)%present", its)
        self.assertIn("ITS(var) = error", its)

    def test_the_arrays_are_wide_enough_for_the_bound(self):
        """A bound past an array's extent is a silent overrun, not a crash.

        -fbounds-check is on, so this would abort at runtime - but only on a
        project with more than four gases, which no default fixture has.
        """
        typedef = read("src/src_common/m_typedef.f90")
        self.assertIn("rand_uncer(E2NumVar)", typedef)
        self.assertIn("mahrt98_NR(GHGNumVar)", typedef)
        self.assertIn("ITS(E2NumVar)", read("src/src_common/m_common_global_var.f90"))


class RandomUncertaintySettingsReachTheEngine(unittest.TestCase):
    """The ru_* keys are Project tags and the interface must write them there.

    ParseIniFile is called with the section prefix 'Project', so EPPrjNTags are
    only matched inside sections whose name starts with Project. They have to
    be Project tags: RP and FCC both need ru_meth, and FCC sweeps only
    FluxCorrection*. The interface used to write all three into
    [RawProcess_RandomUncertainty_Settings], where nothing looked for them, so
    RUsetup%meth fell to its `case default` of 'none' and random uncertainty
    never ran for any project it had saved.

    This pins both halves of the agreement, so that moving either side is
    noticed.
    """

    def test_the_ru_keys_are_project_tags(self):
        source = read("src/src_common/m_common_global_var.f90")
        # Positional: write_processing_project_variables.f90 reads
        # EPPrjNTags(23), (24) and (25) by index, not by name.
        for index, tag in ((23, "ru_its_meth"), (24, "ru_meth"), (25, "ru_tlag_max")):
            self.assertIn(
                f"EPPrjNTags({index})%Label / '{tag}' /",
                source,
                f"{tag} moved out of slot {index} of the project tag table, "
                f"which its reader indexes positionally",
            )

    def test_project_tags_are_read_only_from_project_sections(self):
        for path in ("src/src_rp/read_ini_rp.f90", "src/src_fcc/read_ini_fcc.f90"):
            self.assertIn(
                "call ParseIniFile(PrjPath, 'Project', EPPrjNTags, EPPrjCTags",
                read(path),
                f"{path} changed which sections the project tags are swept from",
            )

    def test_the_rp_table_no_longer_duplicates_them(self):
        """Duplicates nothing read, and they are what made the keys look
        like RawProcess settings in the first place. Blanked, not deleted -
        the tables are positional."""
        source = read("src/src_rp/m_rp_global_var.f90")
        for tag in ("'ru_meth'", "'ru_its_meth'", "'ru_tlag_max'", "'ru_its_sec_factor'"):
            self.assertNotIn(f"%Label / {tag} /", source,
                             f"{tag} is back in the RP tag table; it is a "
                             f"Project tag and a duplicate there invites the "
                             f"interface to write it into a RawProcess group")

    def test_the_interface_writes_them_into_project(self):
        gui = ROOT.parent / "eddyflow-gui" / "src" / "ecinidefs.h"
        if not gui.is_file():
            self.skipTest("GUI checkout not present")
        source = gui.read_text(encoding="utf-8", errors="replace")
        self.assertIn("const auto INIGROUP_RAND_ERROR = INIGROUP_PROJECT;", source)
        self.assertIn("INIGROUP_RAND_ERROR_LEGACY", source,
                      "the legacy group name must survive so old files can "
                      "still be read and their stale copies cleared")


if __name__ == "__main__":
    unittest.main()
