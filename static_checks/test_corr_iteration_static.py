"""The iterative correction, and the three things that make it safe.

The spectral correction is evaluated at z/L; z/L comes from the corrected heat
flux; that flux is what the spectral correction produces. One pass leaves them
disagreeing. EddyUH closes the circle by repeating (``EddyUH.m:722-903``); this
makes the repetition optional.

Three properties carry the design:

1. **Off is one pass, and one pass is exactly what the block did before.** The
   loop bound is 1 unless the project asks otherwise, so byte-identity is
   structural rather than incidental.
2. **Nothing accumulates.** ``Fluxes1_rp`` rebuilds ``Flux1`` from ``Flux0``
   and ``BPCF``; ``Fluxes23_rp`` rebuilds ``Flux2``/``Flux3`` from ``Flux1``.
   Each pass corrects the same raw covariances afresh. If either ever started
   reading its own previous output, iterating would compound rather than
   converge - so that is asserted here, not assumed.
3. **RP and FCC stay twinned.** Both run the loop, and their full outputs must
   carry the same column in the same place.

Worth recording, because it contradicts what this repository said before:
**EddyUH has no convergence test.** Its loop is ``while indexITER <= 3`` with
the counter incremented at the top, so it runs four times, always. Its only
``break`` tests the urban footprint's roughness length, requires
``indexITER > 3`` - true on the final pass only - and lives in a branch that
``EddyUH_footprint.m:151`` makes unreachable. Its ``covsvar`` output is a
reported diagnostic that nothing reads.
"""

import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RP_MAIN = ROOT / "src" / "src_rp" / "eddyflow-rp_main.f90"
FCC_MAIN = ROOT / "src" / "src_fcc" / "eddyflow-fcc_main.f90"
FLUXES1 = ROOT / "src" / "src_rp" / "fluxes1_rp.f90"
FLUXES23 = ROOT / "src" / "src_rp" / "fluxes23_rp.f90"
HELPER = ROOT / "src" / "src_common" / "corr_iteration.f90"
DECODER = ROOT / "src" / "src_common" / "write_processing_project_variables.f90"
TAGS = ROOT / "src" / "src_common" / "m_common_global_var.f90"
RP_HDR = ROOT / "src" / "src_rp" / "init_outfiles_rp.f90"
RP_ROW = ROOT / "src" / "src_rp" / "write_out_full.f90"
FCC_HDR = ROOT / "src" / "src_fcc" / "init_out_files.f90"
FCC_ROW = ROOT / "src" / "src_fcc" / "write_out_full_fcc.f90"
EDDYUH = (ROOT.parent / "EddyUH_testing" / "EddyUH" / "EddyUH_1.7b_COS"
          / "EC_Software_FluxCalc" / "EddyUH.m")


def read(path):
    return path.read_text(encoding="utf-8", errors="replace")


def code(path):
    return re.sub(r"^\s*!.*$", "", read(path), flags=re.M)


def loop(path):
    src = read(path)
    i = src.index("corr_passes = 1")
    return src[i:src.index("end do", i)]


class OffIsOnePass(unittest.TestCase):

    def test_the_bound_is_one_unless_asked(self):
        for path in (RP_MAIN, FCC_MAIN):
            body = loop(path)
            self.assertIn("corr_passes = 1", body)
            self.assertIn("if (EddyFlowProj%corr_iter_meth)", body)
            self.assertIn("corr_passes = EddyFlowProj%corr_iter_max", body)

    def test_the_default_is_off(self):
        src = read(DECODER)
        i = src.index("EddyFlowProj%corr_iter_meth = .false.")
        block = src[i:i + 800]
        self.assertIn("EddyFlowProj%corr_iter_meth = .false.", block)
        #> EddyUH's own numbers, so switching it on reproduces EddyUH rather
        #> than something chosen here.
        self.assertIn("EddyFlowProj%corr_iter_max = 4", block)
        self.assertIn("EddyFlowProj%corr_iter_tol = 0d0", block)

    def test_a_zero_tolerance_never_exits_early(self):
        #> Which is EddyUH's behaviour: it runs every pass and tests nothing.
        #> A `>= tol` test would have made zero mean "exit immediately".
        for path in (RP_MAIN, FCC_MAIN):
            self.assertIn("EddyFlowProj%corr_iter_tol > 0d0", loop(path))

    def test_the_three_keys_are_guarded_reads(self):
        src = read(DECODER)
        for slot in (8, 9, 10):
            self.assertIn("if (EPPrjNTagFound(%d))" % slot, src)

    def test_a_pass_count_below_one_falls_back(self):
        #> Zero or negative would skip the correction entirely, which is not
        #> a setting anyone means - "no iteration" is what off already says.
        src = read(DECODER)
        self.assertIn("if (EddyFlowProj%corr_iter_max < 1) "
                      "EddyFlowProj%corr_iter_max = 4", src)


class NothingAccumulates(unittest.TestCase):

    def test_every_flux_level_is_cleared_at_entry(self):
        """The property that actually makes the loop safe.

        Not "no self-reference": Fluxes23_rp legitimately writes
        ``Flux3%gas(wsl) = Flux3%gas(wsl) * BPCF%of(wsl)``, applying the
        water correction to a value it computed a few lines earlier in the
        SAME call. That is fine. What would not be fine is carrying a value
        from the previous CALL, and the reset at the head of each routine is
        what rules that out - every level is rebuilt from scratch, so a
        second pass corrects the same raw covariances rather than the first
        pass's output.
        """
        for path, levels in ((FLUXES1, ["Flux1"]),
                             (FLUXES23, ["Flux2", "Flux3"])):
            body = code(path)
            body = body[body.index("subroutine "):]
            for level in levels:
                reset = "%s = errFlux" % level
                self.assertIn(reset, body,
                              "%s never clears %s" % (path.name, level))
                #> Before anything else touches that level.
                first_use = re.search(r"%s%%\w+[^=\n]*=" % level, body)
                self.assertLess(body.index(reset), first_use.start(),
                                "%s writes %s before clearing it"
                                % (path.name, level))

    def test_level_one_is_rebuilt_from_level_zero(self):
        self.assertIn("Flux1%gas(msl) = Flux0%gas(msl) * BPCF%of(msl)",
                      code(FLUXES1))

    def test_the_flux_stages_never_write_level_zero(self):
        #> The raw covariances each pass starts from. If a pass rewrote them,
        #> the second pass would be correcting the first pass's output.
        for path in (FLUXES1, FLUXES23):
            self.assertNotRegex(code(path), r"^\s*Flux0%\w+\s*=[^=]",
                                "%s assigns Flux0" % path.name)

    def test_neither_stage_keeps_state_between_calls(self):
        for path in (FLUXES1, FLUXES23):
            self.assertNotRegex(code(path), r"\bsave\b")

    def test_the_feedback_is_the_stability(self):
        #> Fluxes23 recomputes it and the next pass's correction reads it.
        #> That path existed already; only the repetition was missing.
        self.assertIn("Ambient%zL =", code(FLUXES23))
        self.assertIn("Ambient%zL", loop(RP_MAIN))

    def test_fcc_threads_the_stability_through_a_local(self):
        #> FCC's correction is handed lEx%Flux0%zL while Fluxes23 writes
        #> lEx%zL - two different fields, so the feedback does not close by
        #> itself the way RP's single Ambient%zL does.
        body = loop(FCC_MAIN)
        #> Seeded from the ex record's Level-0 value, then refreshed from
        #> what Fluxes23 produced.
        self.assertIn("iter_zL = lEx%Flux0%zL", body)
        self.assertIn("iter_zL = lEx%zL", body)
        #> And the correction reads the local, not the record - otherwise
        #> every pass would be handed the same starting stability and the
        #> loop would compute the same answer four times.
        call = body[body.index("call BandPassSpectralCorrections"):]
        call = call[:call.index("call Fluxes1")]
        self.assertIn("iter_zL", call)
        self.assertNotIn("Flux0%zL", call)


class TheTwoApplicationsStayTwinned(unittest.TestCase):

    def test_both_run_the_loop(self):
        for path in (RP_MAIN, FCC_MAIN):
            body = loop(path)
            self.assertIn("call BandPassSpectralCorrections", body)
            self.assertIn("WorstRelativeChange(prev_gas_flux", body)

    def test_the_column_is_written_only_when_the_loop_runs(self):
        #> A column of -9999 in every file that does not iterate would say
        #> nothing. Same rule the random-error columns follow.
        for path in (RP_HDR, RP_ROW, FCC_HDR, FCC_ROW):
            self.assertIn("if (EddyFlowProj%corr_iter_meth) then", read(path),
                          "%s writes the column unconditionally" % path.name)

    def test_the_column_sits_in_the_same_place_in_both(self):
        for hdr in (RP_HDR, FCC_HDR):
            src = read(hdr)
            self.assertLess(
                src.index("'corr_iter_dev'"),
                src.index("'corrected_fluxes_and_quality_flags,'"),
                "%s puts the column after the flux block" % hdr.name)

    def test_fcc_computes_its_own_convergence(self):
        #> Not read from the ex record. The two applications run their own
        #> loops over their own corrections, so RP's number would describe a
        #> different calculation.
        self.assertIn("lEx%corr_iter_dev = iter_dev", read(FCC_MAIN))
        self.assertNotIn("lEx%corr_iter_dev", read(ROOT / "src" / "src_common"
                                                   / "read_ex_record.f90"))


class TheConvergenceMeasure(unittest.TestCase):

    def test_it_is_the_worst_gas_not_each_one(self):
        src = read(HELPER)
        self.assertIn("if (dev > WorstRelativeChange) then", src)
        self.assertIn("do gas = firstGas, lastGas", src)

    def test_a_zero_or_missing_flux_is_skipped_not_counted(self):
        #> An error code counted as an infinite change, or a division by a
        #> flux that was exactly zero, would make a period with nothing to
        #> converge look like the worst in the run.
        src = read(HELPER)
        self.assertIn("if (before(gas) == error .or. after(gas) == error) cycle",
                      src)
        self.assertIn("if (before(gas) == 0d0) cycle", src)

    def test_it_is_a_percentage(self):
        self.assertIn("* 100d0", read(HELPER))


class WhatEddyUHActuallyDoes(unittest.TestCase):
    """Pinned because this repository asserted otherwise for a while."""

    @unittest.skipUnless(EDDYUH.is_file(), "EddyUH tree not beside this one")
    def test_its_loop_is_a_fixed_count(self):
        src = read(EDDYUH)
        self.assertIn("max_iterations = 3", src)
        self.assertIn("while indexITER <= max_iterations", src)

    @unittest.skipUnless(EDDYUH.is_file(), "EddyUH tree not beside this one")
    def test_its_only_break_is_the_urban_roughness_length(self):
        #> NOT a flux-convergence test. If this ever fails, EddyUH gained a
        #> real convergence criterion and the defaults here should follow it.
        src = read(EDDYUH)
        breaks = re.findall(r"^[^%\n]*\bbreak\b[^\n]*$", src, re.M)
        self.assertEqual(len(breaks), 1, "EddyUH's loop has more than one exit")
        self.assertIn("z_d_ave", breaks[0])
        self.assertIn("indexITER > 3", breaks[0])

    @unittest.skipUnless(EDDYUH.is_file(), "EddyUH tree not beside this one")
    def test_covsvar_is_reported_and_not_tested(self):
        src = read(EDDYUH)
        uses = re.findall(r"^[^%\n]*covsvar[^\n]*$", src, re.M)
        self.assertTrue(uses)
        for line in uses:
            self.assertNotRegex(line, r"\bif\b.*covsvar",
                                "EddyUH now tests covsvar: %s" % line.strip())


if __name__ == "__main__":
    unittest.main()
