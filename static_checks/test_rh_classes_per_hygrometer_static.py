"""Every hygrometer's lag is classed by its own humidity, from its own period.

Two defects, found together while making a split PWB pre-pass reproduce a
serial one.

**The humidity was one period late.** The aggregate time-lag summary bins each
water lag by relative humidity. Under PWB its row was added straight after the
time-lag handling - before the period's FluxParams, which is what sets the
humidity. Stats%RH then still held the previous period's value, so every
period landed in the class of the one before it, and the first in the class of
whatever the pre-pass left behind. A serial pre-pass leaves the run's last
period there, a split one the end of its first piece, and on base_gappy_pwb
that put one period in a neighbouring class. The row now gets its humidity
after FluxParams (SetPwbTimelagSummaryRH).

**The pre-pass humidity depended on earlier periods.** The cache-generation
pre-pass shifted each gas by the streaming classifier's lag - carried, borrowed
or S2-anchored, all functions of the periods before - and took the humidity
from the shifted water. It now shifts by each gas's own-evidence lag
(pwb_raw_OwnRowLags), so what an assessment-only run writes is the same
however the walk was cut.

**Only the designated hygrometer was classed.** One table, one humidity: a
second hygrometer took a single window, its determinations were dropped
whenever the FIRST hygrometer had no humidity, and SetTimelags read only the
first one's classes. Now every closed-path hygrometer (WaterSlotClassed) has
its own table, binned by the biomet RH where the site has one and otherwise by
its own (PeriodWaterRH), written after the designated one's under
`..._for_<label>`, read back per slot, and looked up per slot. The designated
hygrometer's table keeps its title, columns and position, so older readers and
the interface's file check still find it.

Part of the EddyFlow engine's static checks.
"""

from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]


def read(path):
    return (ROOT / path).read_text(encoding="utf-8")


def code(path):
    return (chr(10)).join(ln for ln in read(path).splitlines()
                          if not ln.lstrip().startswith("!"))


def body_of(source, opener, closer):
    return source[source.index(opener):source.index(closer)]


MAIN = code("src/src_rp/eddyflow-rp_main.f90")
PWB = code("src/src_rp/pwb_timelag_handle.f90")
TLH = code("src/src_rp/timelag_handle.f90")
FLUX = code("src/src_rp/flux_params.f90")
OPT = code("src/src_rp/optimize_timelags.f90")
FIX = code("src/src_rp/fix_timelag_opt_dataset.f90")
ADD = code("src/src_rp/add_to_timelag_opt_dataset.f90")
WRITE = code("src/src_rp/writeout_timelag_optimization.f90")
READ = code("src/src_rp/read_timelag_opt_file.f90")
SET = code("src/src_rp/set_timelags.f90")
TYPES = code("src/src_common/m_typedef.f90")
GLOB = code("src/src_rp/m_rp_global_var.f90")


class TheHumidityIsThisPeriods(unittest.TestCase):

    def test_the_summary_row_takes_no_humidity_when_it_is_added(self):
        add = body_of(PWB, "subroutine AddPwbTimelagSummaryDataset",
                      "end subroutine AddPwbTimelagSummaryDataset")
        self.assertNotIn("Stats%RH", add)
        self.assertNotIn("RecordWaterRH", add)

    def test_it_is_filled_after_flux_params(self):
        added = MAIN.index("call AddPwbTimelagSummaryDataset(")
        flux = MAIN.index("call FluxParams(.true.)")
        filled = MAIN.index("call SetPwbTimelagSummaryRH(PwbTimelagOpt", flux)
        self.assertLess(added, flux)
        self.assertLess(flux, filled)

    def test_a_period_without_flux_params_is_gated_without_humidity(self):
        self.assertIn("PwbTimelagN, .false.)", MAIN)

    def test_the_pre_pass_records_after_flux_params_too(self):
        self.assertLess(MAIN.index("call FluxParams(.false.)"),
                        MAIN.index("call RecordPwbTimelagOptPeriod("))


class ThePrePassHumidityHasNoHistory(unittest.TestCase):

    def test_own_evidence_lags_are_captured_in_pass_one(self):
        self.assertIn("pwb_raw_OwnRowLags(j) = lPwbResult%row_lag", TLH)
        self.assertIn("pwb_raw_OwnRowLags(j) = mc_row", TLH)

    def test_the_pre_pass_apply_uses_them_and_production_does_not(self):
        self.assertIn("if (PwbCacheGenerate .and. InTimelagOpt) then", TLH)
        i = TLH.index("if (PwbCacheGenerate .and. InTimelagOpt) then")
        self.assertIn("RowLags(j) = pwb_raw_OwnRowLags(j)", TLH[i:i + 400])

    def test_the_global_exists(self):
        self.assertIn("pwb_raw_OwnRowLags(E2NumVar)", GLOB)


class EveryHygrometerHasItsOwnHumidity(unittest.TestCase):

    def test_one_humidity_per_slot(self):
        opt = body_of(TYPES, "type :: TimeLagOptType", "end type TimeLagOptType")
        self.assertIn("RH(E2NumVar)", opt)

    def test_biomet_first_then_the_hygrometers_own(self):
        rh = body_of(FLUX, "function PeriodWaterRH", "end function PeriodWaterRH")
        self.assertLess(rh.index("biomet%val(bRH)"), rh.index("Ambient%RH_at(slot)"))

    def test_the_designated_hygrometers_rh_at_is_stats_rh(self):
        #> What keeps a one-hygrometer project binned by the same number.
        self.assertIn("Ambient%RH_at(wsl)    = Stats%RH", FLUX)

    def test_each_slot_is_kept_on_its_own_humidity(self):
        self.assertIn("TimelagOpt(i)%RH(gas) == error", FIX)
        self.assertIn("toSet(actn(gas))%RH(gas) = TimelagOpt(i)%RH(gas)", FIX)

    def test_tlag_opt_gates_each_hygrometer_on_its_own_latent_heat(self):
        self.assertIn("le = Flux0%LE", ADD)
        self.assertIn("le = Flux0%gas(gas) * Ambient%lambda * MW_H2O * 1d-3", ADD)
        self.assertIn("rh = PeriodWaterRH(gas)", ADD)


class EveryClosedPathHygrometerIsClassed(unittest.TestCase):

    def test_the_rule(self):
        rule = body_of(OPT, "function WaterSlotClassed", "end function WaterSlotClassed")
        for test in ("TOSetup%h2o_nclass <= 1", "GasSlotIsWater(slot)",
                     "E2Col(slot)%present", "path_type == 'open'"):
            self.assertIn(test, rule)

    def test_the_open_path_override_is_gone(self):
        self.assertNotIn("TOSetup%h2o_nclass = 1", MAIN)

    def test_classes_are_per_slot(self):
        self.assertIn("toH2O(toMaxH2OClass, E2NumVar)", GLOB)
        self.assertIn("if (WaterSlotClassed(gas)) then", OPT)
        self.assertIn("toSet(i)%RH(gas)", OPT)
        self.assertIn("allocate(toH2On(TOSetup%h2o_nclass, E2NumVar))", MAIN)


class TheFileCarriesATablePerHygrometer(unittest.TestCase):

    def test_the_designated_table_keeps_its_title(self):
        self.assertIn("call WriteRhTable(wsl, 'H2O_timelag_determinations_as_a_function_of_relative_humidity')",
                      WRITE)

    def test_the_others_follow_named(self):
        self.assertIn("'H2O_timelag_determinations_as_a_function_of_relative_humidity_for_'", WRITE)
        self.assertLess(WRITE.index("call WriteRhTable(wsl,"),
                        WRITE.index("_relative_humidity_for_'"))

    def test_a_blank_line_ends_each_table(self):
        i = WRITE.index("_relative_humidity_for_'")
        self.assertIn("write(uto, '(a)')", WRITE[i - 200:i])

    def test_the_reader_reads_every_table(self):
        block = READ[READ.index("if (index(strg, 'H2O_timelag_determinations_as_a_function') /= 0) then"):]
        block = block[:block.index("end if" + chr(10) + "        end do")]
        self.assertNotIn(chr(10) + "                exit" + chr(10), block)
        self.assertIn("k = index(strg, '_for_')", block)
        self.assertIn("toH2O(nrows, slot)%def", block)
        self.assertIn("nrows > toMaxH2OClass", block)


class EachHygrometerLooksUpItsOwnClass(unittest.TestCase):

    def test_lookup_per_slot(self):
        self.assertIn("if (WaterSlotClassed(gas)) then", SET)
        self.assertIn("call LocalRhEstimate(lRH, gas)", SET)
        self.assertIn("toH2O(cls, gas)%max > toH2O(cls, gas)%min", SET)

    def test_the_local_estimate_takes_the_slot(self):
        self.assertIn("subroutine LocalRhEstimate(lRH, wsl)", SET)
        self.assertIn("integer, intent(in) :: wsl", SET)

    def test_an_empty_class_falls_back_to_the_plain_window_except_for_the_designated(self):
        self.assertIn("if (gas == wsl) classed = .true.", SET)
        self.assertIn("if (.not. classed) then", SET)


if __name__ == "__main__":
    unittest.main()
