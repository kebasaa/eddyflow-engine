"""The random errors are reported spectrally corrected, like their fluxes.

Langford et al. (2015) recommend that "any corrections applied to the
calculated fluxes (e.g. low or high frequency loss corrections) should also be
applied to the random error". With ``ru_spec_corr`` on - off by default, and
for a project without the key, because EddyPro does not do it and the defaults
reproduce EddyPro - each random error, and with it the Langford
limit of detection, is multiplied by its own flux's spectral correction factor
in ``ReportedRandomErrors`` (src_common/reported_random_errors.f90).

Three ways it could go wrong silently:

**Corrected twice.** RP's FLUXNET file is FCC's input. RP corrects only when
it writes the final outputs itself (FCC does not follow). FCC corrects unless
the project is one RP finishes on its own - then the file it reads is corrected
already - which it decides with RP's own test, before it changes any setting
that test reads.

**The wrong factor.** Tau takes w_u's, H w_ts's, LE and ET the designated
hygrometer's, each gas its own - the same BPCF%of the flux was multiplied by.

**A writer left reading the raw value.** Every write of a random error, and the
CEC significance test, goes through the reported copy.

Part of the EddyFlow engine's static checks.
"""

from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"


def code(*parts):
    text = SRC.joinpath(*parts).read_text(encoding="utf-8", errors="replace")
    return "\n".join(ln for ln in text.splitlines() if not ln.lstrip().startswith("!"))


ROUTINE = code("src_common", "reported_random_errors.f90")
DECODE = code("src_common", "write_processing_project_variables.f90")
GEN = (ROOT / "prj" / "gen_project_tags.py").read_text(encoding="utf-8")

RP_WRITERS = (("src_rp", "write_out_full.f90"), ("src_rp", "write_out_fluxnet.f90"))
FCC_WRITERS = (("src_fcc", "write_out_full_fcc.f90"), ("src_fcc", "write_out_fluxnet_fcc.f90"))


class TheSetting(unittest.TestCase):

    def test_the_key_sits_in_project(self):
        self.assertIn('14: "ru_spec_corr"', GEN)

    def test_it_is_off_unless_the_project_says_one(self):
        self.assertIn("RUsetup%spec_corr = .false.", DECODE)
        self.assertIn("if (EPPrjNTagFound(14)) RUsetup%spec_corr = nint(EPPrjNTags(14)%value) == 1",
                      DECODE)
        self.assertIn("logical :: spec_corr = .false.", code("src_common", "m_typedef.f90"))


class TheFactors(unittest.TestCase):

    def test_each_error_takes_its_own_flux_factor(self):
        self.assertIn("re_out(u) = Scaled(re_in(u), BPCF%of(w_u))", ROUTINE)
        self.assertIn("re_out(ts) = Scaled(re_in(ts), BPCF%of(w_ts))", ROUTINE)
        self.assertIn("re_out(var) = Scaled(re_in(var), BPCF%of(var))", ROUTINE)
        self.assertIn("le_out = Scaled(le_in, BPCF%of(wsl))", ROUTINE)
        self.assertIn("et_out = Scaled(et_in, BPCF%of(wsl))", ROUTINE)

    def test_it_can_be_switched_off_and_needs_the_caller_to_ask(self):
        self.assertIn("if (.not. (apply .and. RUsetup%spec_corr)) return", ROUTINE)

    def test_missing_values_are_not_invented(self):
        self.assertIn("if (re == error .or. factor == error .or. factor <= 0d0) then", ROUTINE)


class NothingIsCorrectedTwice(unittest.TestCase):

    def test_rp_corrects_only_when_it_writes_the_final_outputs(self):
        for parts in RP_WRITERS:
            src = code(*parts)
            self.assertIn("call ReportedRandomErrors(.not. EddyFlowProj%fcc_follows, Essentials%rand_uncer,",
                          src, parts)

    def test_fcc_corrects_unless_rp_already_did(self):
        for parts in FCC_WRITERS:
            self.assertIn("call ReportedRandomErrors(.not. RUsetup%rp_corrected, lEx%rand_uncer,",
                          code(*parts), parts)
        main = code("src_fcc", "eddyflow-fcc_main.f90")
        read = main.index("call ReadIniFCC('FluxCorrection')")
        decided = main.index("RUsetup%rp_corrected = .not. FccFollowsProject()")
        fallback = main.index("EddyFlowProj%hf_meth = 'moncrieff_97'")
        self.assertLess(read, decided)
        self.assertLess(decided, fallback)

    def test_rp_and_fcc_share_one_test(self):
        self.assertIn("if (FccFollowsProject()) then", code("src_rp", "eddyflow-rp_main.f90"))

    def test_no_writer_reads_the_raw_value_afterwards(self):
        for parts in RP_WRITERS + FCC_WRITERS:
            src = code(*parts)
            after = src[src.index("call ReportedRandomErrors("):]
            after = after[after.index("rep_et)"):]
            self.assertNotIn("Essentials%rand_uncer", after, parts)
            self.assertNotIn("lEx%rand_uncer", after, parts)


if __name__ == "__main__":
    unittest.main()
