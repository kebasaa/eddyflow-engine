"""Each hygrometer is WPL-corrected after its own path type, each gas after the
hygrometer it names, and nothing at all when WPL is off.

Three defects, all in `Fluxes23_rp` and its FCC twin `Fluxes23`:

  - Only the designated hygrometer got the evapotranspiration formula. A
    second one went through the trace-gas formula, which on a closed path has
    no (1 + mu*sigma) dilution of the flux itself and on an open path used the
    designated hygrometer's evapotranspiration for its own. Its numbered LE and
    ET therefore carried a different correction from the bare ones.

  - Every open-path gas took the designated hygrometer's evapotranspiration for
    its WPL water-vapour term, and the LI-7700 its E without WPL, whatever water
    the gas's record named. With a MIRO and an LI-7550 side by side, the MIRO's
    gases were corrected with the MIRO's humidity and the LI-7550's flux.

  - `Level2GasFlux` never read the WPL switch. Only the designated hygrometer
    and the point-by-point mixing ratios obeyed it, so with WPL off every gas
    was still corrected, open-path ones against an evapotranspiration the
    switch had left uncorrected. EddyPro 7 resets each gas to its level 1 flux
    when WPL is off; the 6.2.2 fork this engine came from did not.

And one found on the way: the numbered H of a second hygrometer was built from
Flux1%H, which carries no spectral factor, beside a bare H that does.
"""

from pathlib import Path
import re
import unittest


ROOT = Path(__file__).resolve().parents[1]

FLUXES_RP = "src/src_rp/fluxes23_rp.f90"
FLUXES_FCC = "src/src_fcc/fluxes23.f90"
READER = "src/src_common/read_ex_record.f90"
FLX_ROW_RP = "src/src_rp/write_out_fluxnet.f90"
FLX_ROW_FCC = "src/src_fcc/write_out_fluxnet_fcc.f90"
TYPEDEF = "src/src_common/m_typedef.f90"
FCC_MAIN = "src/src_fcc/eddyflow-fcc_main.f90"
E2SET = "src/src_common/define_e2_set.f90"

TWINS = (FLUXES_RP, FLUXES_FCC)


def read(path):
    return (ROOT / path).read_text(encoding="utf-8")


def code(path):
    """Source with comment-only lines dropped."""
    return "\n".join(ln for ln in read(path).splitlines()
                     if not ln.lstrip().startswith("!"))


def routine(path, name):
    body = code(path)
    start = re.search(r"^\s*subroutine %s\(" % name, body, re.M)
    assert start, "%s: no subroutine %s" % (path, name)
    end = body.index("end subroutine %s" % name, start.start())
    return body[start.start():end]


class TheWplSwitchReachesEveryGas(unittest.TestCase):

    def test_level2gasflux_reads_it(self):
        for path in TWINS:
            body = routine(path, "Level2GasFlux")
            self.assertIn("if (.not. EddyFlowProj%wpl) then", body, path)

    def test_off_means_the_level_one_flux(self):
        """No Va ratio, no LI-7700 multiplier: what EddyPro 7 does."""
        for path in TWINS:
            body = routine(path, "Level2GasFlux")
            off = body[body.index("if (.not. EddyFlowProj%wpl) then"):]
            off = off[:off.index("end if\n", off.index("return"))]
            self.assertIn("Flux2%gas(gas) = Flux1%gas(gas)", off, path)
            self.assertNotIn("Mul7700", off, path)
            self.assertNotIn("Va", off, path)

    def test_it_is_decided_before_any_correction(self):
        for path in TWINS:
            body = routine(path, "Level2GasFlux")
            self.assertLess(body.index("if (.not. EddyFlowProj%wpl) then"),
                            body.index("path_type == 'closed'"), path)

    def test_every_hygrometer_obeys_it_too(self):
        for path in TWINS:
            body = routine(path, "HygrometerWater")
            self.assertIn("if (EddyFlowProj%wpl) then", body, path)
            self.assertIn("E2 = E1", body, path)


class EachGasUsesTheWaterFluxOfTheHygrometerItNames(unittest.TestCase):

    def test_the_flux_is_passed_in_not_read_from_the_host(self):
        for path in TWINS:
            body = routine(path, "Level2GasFlux")
            self.assertIn("subroutine Level2GasFlux(gas, sigma_gas, rhow_gas, "
                          "e3_ref, enw_ref)", body, path)
            self.assertNotRegex(body, r"\bFlux3%E\b(?!_)", path)
            self.assertNotRegex(body, r"\bE_nowpl\b", path)

    def test_the_caller_chooses_it_per_gas(self):
        for path in TWINS:
            body = code(path)
            self.assertIn("call WaterFluxFor(msl, E3_ref, Enowpl_ref)", body, path)
            self.assertIn("call Level2GasFlux(msl, sigma_g, rhow_g, E3_ref, "
                          "Enowpl_ref)", body, path)

    def test_the_choice_is_the_gas_moisture_reference(self):
        rp = routine(FLUXES_RP, "WaterFluxFor")
        self.assertIn("ref = E2Col(gas)%moist_ref", rp)
        fcc = routine(FLUXES_FCC, "WaterFluxFor")
        self.assertIn("ref = lEx%moist_ref_at(gas)", fcc)
        for path, body in ((FLUXES_RP, rp), (FLUXES_FCC, fcc)):
            self.assertIn("e3 = E3w(ref)", body, path)
            self.assertIn("enw = Enowplw(ref)", body, path)
            #> The designated hygrometer is the fallback, and the default.
            self.assertLess(body.index("e3 = Flux3%E"),
                            body.index("e3 = E3w(ref)"), path)

    def test_a_biomet_reference_is_said_once_at_setup(self):
        """The fallback is a property of the project, so it is reported where
        the references are resolved - not from the flux routines, which run
        several times per period and keep no state."""
        body = code(E2SET)
        self.assertIn("logical, save :: borrowFluxNoted(GHGNumVar)", body)
        self.assertIn("E2Col(j)%moist_ref /= biometMoistRef) cycle", body)
        for path in TWINS:
            self.assertNotRegex(code(path), r"\bsave\b", path)


class EveryHygrometerGetsItsOwnEvapotranspiration(unittest.TestCase):

    def test_one_routine_for_all_of_them(self):
        for path in TWINS:
            body = code(path)
            self.assertIn("call HygrometerWater(wsl, Flux1%E,", body, path)
            other = routine(path, "OtherHygrometerWater")
            self.assertIn("call HygrometerWater(ws, E1,", other, path)

    def test_the_others_are_walked_from_the_shared_list(self):
        for path in TWINS:
            body = code(path)
            self.assertIn("call WaterOutSlots(w_slots, w_tags, n_w_slots)",
                          body, path)
            self.assertIn("call OtherHygrometerWater(w_slots(iw))", body, path)

    def test_no_hygrometer_takes_the_trace_gas_path(self):
        for path in TWINS:
            body = code(path)
            self.assertEqual(
                body.count("if (msl == wsl .or. isHygro(msl)) cycle"), 2,
                "%s: both the level 2 and the level 3 gas loops must skip "
                "every hygrometer" % path)

    def test_each_uses_its_own_level_three_heat(self):
        for path in TWINS:
            other = routine(path, "OtherHygrometerWater")
            self.assertIn("Flux3%H_at(ws)", other, path)

    def test_the_heat_is_computed_before_the_water(self):
        for path, heat in ((FLUXES_RP, "call PerHygrometerHeat_rp()"),
                           (FLUXES_FCC, "call PerHygrometerHeat()")):
            body = code(path)
            self.assertLess(body.index(heat),
                            body.index("call OtherHygrometerWater("), path)

    def test_their_own_air_and_cell(self):
        rp = routine(FLUXES_RP, "OtherHygrometerWater")
        for term in ("Ambient%sigma_at(ws)", "RHO%w_at(ws)",
                     "Ambient%RhoCp_at(ws)", "Ambient%Tcell_at(ws)",
                     "Ambient%Pcell_at(ws)", "cellPressureSlot(ws)"):
            self.assertIn(term, rp)
        fcc = routine(FLUXES_FCC, "OtherHygrometerWater")
        for term in ("lEx%sigma_at(ws)", "lEx%rhow_at(ws)",
                     "lEx%rhocp_at(ws)"):
            self.assertIn(term, fcc)

    def test_the_closed_path_spectral_factor_is_guarded(self):
        for path in TWINS:
            other = routine(path, "OtherHygrometerWater")
            self.assertIn("if (BPCF%of(ws) /= error) then", other, path)


class MissingLi7700MultipliersMeanAMissingFlux(unittest.TestCase):
    """They are error when the gas's water is the biomet RH; multiplied in,
    they reported methane at -8e8."""

    def test_both_twins_guard_them(self):
        for path, mul in ((FLUXES_RP, "Mul7700(gas)%A == error"),
                          (FLUXES_FCC, "lEx%Mul7700(gas)%A == error")):
            body = routine(path, "Level2GasFlux")
            self.assertIn(mul, body, path)
            self.assertLess(body.index(mul), body.index("%A * wpl"), path)


class TheNumberedHeatIsSpectrallyCorrected(unittest.TestCase):

    def test_both_twins_apply_the_sonic_factor(self):
        for path, name in ((FLUXES_RP, "PerHygrometerHeat_rp"),
                           (FLUXES_FCC, "PerHygrometerHeat")):
            body = routine(path, name)
            self.assertIn("Flux3%H_at(wslot) = Flux3%H_at(wslot) * BPCF%of(w_ts)",
                          body, path)


class FccKnowsWhichWaterEachGasNamed(unittest.TestCase):

    def test_rp_writes_the_reference(self):
        body = code(FLX_ROW_RP)
        self.assertIn("indx = E2Col(var)%moist_ref", body)
        self.assertIn("call AddIntDatumToDataline(indx, csv_row", body)

    def test_the_reader_takes_it_and_fcc_echoes_it(self):
        self.assertIn("integer :: moist_ref_at(GHGNumVar)", code(TYPEDEF))
        self.assertIn("lEx%moist_ref_at(gas) = moist_ref", code(READER))
        self.assertIn("AddIntDatumToDataline(lEx%moist_ref_at(gas)",
                      code(FLX_ROW_FCC))

    def test_an_older_ex_file_is_refused(self):
        body = code(FCC_MAIN)
        self.assertIn("index(header, '_MOIST_REF') <= 0", body)


if __name__ == "__main__":
    unittest.main()
