"""A site with two hygrometers has two humidities, and said so about one.

`LE`, `ET` and `H` were single columns computed from the designated water
record. Everything downstream of humidity followed: one evapotranspiration, one
sensible heat, one Monin-Obukhov length. A second hygrometer produced a water
flux under its own name and then stopped - no latent heat, no
evapotranspiration, no sensible heat, no stability - so the disagreement that
is the whole reason for fielding two instruments could not be seen in any
quantity a reader uses.

Observed on CH-LAE, an LI-7200 and a MIRO MGA4 on one tower: the two report
humidities some twenty per cent apart, and only the LI-7200's reached LE.

Every hygrometer now gets a family. What has to hold:

  - The designated one keeps the bare FLUXNET spellings, and its entry is
    *assigned from* the scalars rather than recomputed beside them. Both
    directions would agree if the per-hygrometer loop were a perfect replay of
    the scalar path, and it nearly is - but the scalars have been through the
    Burba terms, the closed-path spectral correction and the storage chain, and
    only assignment makes agreement certain rather than likely.

  - Header and row walk the same `WaterOutSlots`, which returns slots and
    suffixes together. Getting the list right and the naming wrong shifts a
    file exactly as getting the list wrong does.

  - The FLUXNET block sits after the analyser block and before biomet. Not at
    the end: `ReadExRecord` finds the CEC descriptor by taking the last
    nCecFields of whatever remains, so a block appended after it is read as the
    descriptor. Not earlier either, because the fixed part is parsed by field
    position.

  - `tau` is per hygrometer and `u*` is not. u* comes from the wind covariances
    alone; tau is rho_a * u*^2 and carries the humidity through rho_a. A
    per-hygrometer u* column would be a copy of the one beside it, which is
    worse than absent - a reader would take the agreement for corroboration.

  - The humidity correction of H takes the *level 0* water flux, as the scalar
    does. Reaching for the corrected flux would make the numbered columns
    differ from the bare ones by more than the instruments do.

  - Air temperature is deliberately not per hygrometer. Ta comes from the sonic
    corrected for humidity, so two of them would give the site two air
    temperatures and, through molar volume and every WPL term, two values for
    every gas flux. That is the line that keeps this bounded.
"""

from pathlib import Path
import re
import unittest


ROOT = Path(__file__).resolve().parents[1]

RESOLUTION = "src/src_common/gas_slot_resolution.f90"
TYPEDEF = "src/src_common/m_typedef.f90"
FLUXES_FCC = "src/src_fcc/fluxes23.f90"
FLUXES_RP = "src/src_rp/fluxes23_rp.f90"
PARAMS = "src/src_rp/flux_params.f90"
READER = "src/src_common/read_ex_record.f90"
FLX_HDR = "src/src_rp/init_fluxnet_file_rp.f90"
FLX_ROW_RP = "src/src_rp/write_out_fluxnet.f90"
FLX_ROW_FCC = "src/src_fcc/write_out_fluxnet_fcc.f90"
FO_HDR = "src/src_fcc/init_out_files.f90"
FO_ROW = "src/src_fcc/write_out_full_fcc.f90"
BPCF = "src/src_common/bpcf_aux_subs.f90"
HANDLER = "src/src_common/exception_handler.f90"
FCC_MAIN = "src/src_fcc/eddyflow-fcc_main.f90"


def read(path):
    return (ROOT / path).read_text(encoding="utf-8")


def code(path):
    """Source with comment-only lines dropped.

    Every one of these files explains the retired single-hygrometer behaviour
    in prose right where it used to stand, so a naive search matches the
    explanation rather than live code.
    """
    return "\n".join(ln for ln in read(path).splitlines()
                     if not ln.lstrip().startswith("!"))


class OneListForSlotsAndNames(unittest.TestCase):
    def test_the_resolver_returns_both_together(self):
        body = code(RESOLUTION)
        self.assertIn("subroutine WaterOutSlots(slots, tags, nslots)", body)

    def test_every_writer_asks_it(self):
        for path in (FLX_HDR, FLX_ROW_RP, FLX_ROW_FCC, FO_HDR, FO_ROW):
            self.assertIn(
                "call WaterOutSlots(", code(path),
                "%s must take the hygrometer list from the shared resolver; "
                "deriving it again lets header and row disagree" % path)

    def test_the_suffix_names_the_hygrometer_not_the_count(self):
        """With the designated hygrometer second of three the families are
        H_1, H, H_3 - a gap where the bare name sits. A running count would
        make H_2 mean the first record, which changes meaning when a project
        is edited."""
        body = code(RESOLUTION)
        block = body[body.index("subroutine WaterOutSlots"):]
        block = block[: block.index("end subroutine WaterOutSlots")]
        self.assertIn("write(tags(i), '(a,i0)') '_', i", block)

    def test_only_records_naming_a_column_are_eligible(self):
        body = code(RESOLUTION)
        block = body[body.index("subroutine WaterOutSlots"):]
        block = block[: block.index("end subroutine WaterOutSlots")]
        self.assertIn("%col <= 0) cycle", block)
        self.assertIn("GasSlotIsWater(gas)) cycle", block)


class TheDesignatedEntryIsTheScalar(unittest.TestCase):
    """Assigned, not recomputed - in both executables."""

    def test_both_overwrite_the_entry_from_the_scalar(self):
        for path, lsym, zsym in ((FLUXES_FCC, "lEx%L", "lEx%zL"),
                                 (FLUXES_RP, "Ambient%L", "Ambient%zL")):
            body = code(path)
            self.assertIn("Flux3%H_at(wsl)   = Flux3%H", body, path)
            self.assertIn("Flux3%LE_at(wsl)  = Flux3%LE", body, path)
            self.assertIn("Flux3%ET_at(wsl)  = Flux3%ET", body, path)
            self.assertIn("Flux3%L_at(wsl)   = " + lsym, body, path)
            self.assertIn("Flux3%zL_at(wsl)  = " + zsym, body, path)

    def test_the_scalar_is_never_taken_from_the_entry(self):
        """The reverse assignment would move output already verified against
        v7.2.5, because the loop does not replay Burba, the spectral
        correction or the storage chain."""
        for path in (FLUXES_FCC, FLUXES_RP):
            body = code(path)
            for bad in ("Flux3%H = Flux3%H_at",
                        "Flux3%LE = Flux3%LE_at",
                        "Flux3%ET = Flux3%ET_at"):
                self.assertNotIn(bad, body,
                                 "%s: the bare column must not be recomputed "
                                 "from the per-hygrometer loop" % path)


class MomentumNotFrictionVelocity(unittest.TestCase):
    def test_tau_is_per_hygrometer(self):
        self.assertIn("real(kind = dbl) :: tau_at(GHGNumVar)", code(TYPEDEF))

    def test_ustar_is_not(self):
        """u* is Ambient%us, from the wind covariances. A numbered copy would
        be identical to the bare column and read as corroboration."""
        self.assertNotIn("ustar_at(GHGNumVar)", code(TYPEDEF))
        for path in (FLX_HDR, FO_HDR):
            self.assertNotIn("USTAR' // trim(w_tags", code(path))
            self.assertNotIn("u*' // trim(w_tags", code(path))


class TheSensibleHeatUsesTheLevelZeroFlux(unittest.TestCase):
    """As the scalar does: Flux2%H is built from Flux0%E, not from Flux3."""

    def test_both_twins_take_flux0(self):
        self.assertIn("e0_w = lEx%Flux0%gas(wslot) * MW_H2O * 1d-3",
                      code(FLUXES_FCC))
        self.assertIn("e0_w = Flux0%gas(wslot) * MW_H2O * 1d-3",
                      code(FLUXES_RP))

    def test_the_correction_consumes_it(self):
        for path in (FLUXES_FCC, FLUXES_RP):
            body = code(path)
            self.assertIn("* e0_w / rhoa_w", body, path)
            self.assertNotIn("* Flux3%E_at(wslot) / rhoa_w", body, path)

    def test_the_loop_variable_does_not_shadow_the_wind_index(self):
        """`w` is the module's wind-component parameter. A loop variable of
        that name silently retargets every Stats%Cov(w, ...) in the routine."""
        for path in (FLUXES_FCC, FLUXES_RP):
            body = code(path)
            start = body.index("PerHygrometerFluxes")
            block = body[start: start + 4000]
            self.assertNotIn("integer :: nw, iw, w\n", block, path)
            self.assertIn("wslot", block, path)


class TheFluxnetBlockSitsWhereItCanBeParsed(unittest.TestCase):
    """Between the analyser block and biomet: the CEC descriptor at the end is
    located by taking the last nCecFields, and the fixed part is positional."""

    def test_the_header_puts_it_before_biomet(self):
        body = code(FLX_HDR)
        self.assertLess(body.index("'NUM_WATER_FLUX'"),
                        body.index("'NUM_BIOMET_VARS'"))
        self.assertLess(body.index("'NUM_GAS_INSTR'"),
                        body.index("'NUM_WATER_FLUX'"))

    def test_both_rows_agree_with_the_header(self):
        rp = code(FLX_ROW_RP)
        self.assertLess(rp.index("n_w_flux, csv_row"),
                        rp.index("nbVars, csv_row"))
        fcc = code(FLX_ROW_FCC)
        self.assertLess(fcc.index("n_w_flux, csv_row"),
                        fcc.index("fluxnetChunks%s(6)"))

    def test_the_reader_steps_over_it_before_the_cec_descriptor(self):
        body = code(READER)
        self.assertIn("nWaterFluxFields = 7", body)
        self.assertLess(body.index("n_water_flux"),
                        body.index("lEx%cec%r_ET = error"))

    def test_the_block_is_self_describing(self):
        """A count first, so the reader needs no knowledge of the project."""
        self.assertIn("'NUM_WATER_FLUX'", code(FLX_HDR))
        for path in (FLX_ROW_RP, FLX_ROW_FCC):
            self.assertIn("AddIntDatumToDataline(n_w_flux", code(path))


class TheMoistureBlockCarriesTheWholeRegime(unittest.TestCase):
    """Two terms were enough for the WPL dilution. A per-hygrometer sensible
    heat needs the air that hygrometer implies, and the spectral corrections
    need its humidity."""

    def test_the_ex_record_carries_seven_fields_per_gas(self):
        self.assertIn("nGasMoistFields = 7", code(READER))

    def test_writer_and_header_carry_the_same_seven(self):
        hdr = code(FLX_HDR)
        for term in ("_MOIST_SLOT", "_MOIST_RHOW", "_MOIST_SIGMA",
                     "_MOIST_Q", "_MOIST_RHOA", "_MOIST_RHOCP", "_MOIST_RH"):
            self.assertIn(term, hdr)

    def test_rp_computes_them_per_hygrometer(self):
        body = code(PARAMS)
        for field in ("Ambient%Q_at", "Ambient%rho_a_at",
                      "Ambient%RhoCp_at", "Ambient%RH_at"):
            self.assertIn(field, body)

    def test_air_temperature_is_not_among_them(self):
        """Ta per hygrometer would give every gas flux two values."""
        self.assertNotIn("Ta_at(GHGNumVar)", code(TYPEDEF))
        self.assertNotIn("Ambient%Ta_at", code(PARAMS))

    def test_the_designated_entries_are_the_scalars(self):
        """The loop reproduces them to four places, not exactly - it reaches
        RHO%w_at where the scalar path reaches RHO%w. Four places is nothing
        until it decides which side of an RH-class boundary a period falls on,
        and then it moves a cutoff frequency and every flux behind it."""
        body = code(PARAMS)
        for field, scalar in (("Ambient%e_at(wsl)", "Ambient%e"),
                              ("Ambient%RH_at(wsl)", "Stats%RH"),
                              ("Ambient%rho_d_at(wsl)", "RHO%d"),
                              ("Ambient%rho_a_at(wsl)", "RHO%a"),
                              ("Ambient%Q_at(wsl)", "Ambient%Q"),
                              ("Ambient%RhoCp_at(wsl)", "Ambient%RhoCp")):
            self.assertTrue(
                re.search(re.escape(field) + r" *= *" + re.escape(scalar)
                          + r" *$", body, re.M),
                "%s must be assigned from %s, so a one-hygrometer project "
                "keeps the numbers it had" % (field, scalar))

    def test_sigma_at_is_left_alone(self):
        """It predates this work, the WPL dilution already reads it, and its
        designated entry has carried the loop's value all along. Assigning it
        here would move output verified against v7.2.5."""
        self.assertNotIn("Ambient%sigma_at(wsl) =", code(PARAMS))


class TheCutoffFollowsEachHygrometersOwnHumidity(unittest.TestCase):
    """The RH-to-cutoff relation describes how a tube and filter attenuate
    water at a given humidity - a property of the instrument and its own air.
    Every hygrometer was evaluated at the primary's RH, and the RH classes are
    ten points wide, so CH-LAE's two belong in different ones."""

    def setUp(self):
        body = code(BPCF)
        #> Each arm checked on its own text. Asked of the whole routine, a
        #> reverted arm still matches the other one's line and the check
        #> passes while half the fix is gone - which is what happened.
        self.iir = body[body.index("case('iir')"): body.index("case('sigma')")]
        self.sigma = body[body.index("case('sigma')"):]

    def test_the_exponential_fit_reads_the_gas_own_rh(self):
        self.assertIn("lRH = lEx%rh_at(gas)", self.iir)
        self.assertIn("if (lRH == error) lRH = lEx%RH", self.iir)

    def test_the_primary_rh_no_longer_sets_every_cutoff(self):
        #> Not `lRH = lEx%RH` on its own: that is the fallback for a
        #> hygrometer RP resolved no humidity for, and it is wanted. The
        #> defect was scaling the primary's RH directly into the fit.
        self.assertNotIn("lRH = lEx%RH * 1d-2", self.iir,
                         "evaluating one hygrometer's fit at another's "
                         "humidity answers the question for the wrong sample")

    def test_the_class_reads_the_gas_own_rh(self):
        self.assertIn("lRH = lEx%rh_at(gas)", self.sigma)
        self.assertIn("if (lRH == error) lRH = lEx%RH", self.sigma)

    def test_the_class_is_chosen_per_hygrometer(self):
        """Hygrometer outside, class inside. The other way round picks one
        class from the primary's RH and hands it to every instrument."""
        self.assertLess(self.sigma.index("do gas = firstGas, lastGas"),
                        self.sigma.index("do RH = RH10, RH90"))
        self.assertNotIn("lEx%RH > dfloat(RH)", self.sigma)


class AnOlderExFileIsRefusedRatherThanMisread(unittest.TestCase):
    """Widening the moisture records from three fields to seven made the ex
    file unreadable by position for anything written earlier - and unreadable
    in the worst way. A list-directed read of seven values from a three-field
    record does not fail; it continues into the next gas's fields and returns
    plausible numbers for the wrong slot. Both places that open the file used
    to skip the header unread, so nothing could notice.

    `ex_file=` is a project key: re-running FCC against a kept ex file is a
    supported workflow, not a corner case.
    """

    def test_neither_open_site_skips_the_header_blind(self):
        body = code(FCC_MAIN)
        opens = body.count("open(uex, file = AuxFile%ex")
        checks = body.count("call CheckExFileVintage()")
        self.assertEqual(opens, checks,
                         "every site that opens the ex file must check its "
                         "vintage; %d opens against %d checks"
                         % (opens, checks))

    def test_the_check_looks_for_a_column_only_this_version_writes(self):
        body = code(FCC_MAIN)
        self.assertIn("index(header, 'NUM_WATER_FLUX')", body)
        self.assertIn("subroutine CheckExFileVintage()", body)

    def test_it_stops_rather_than_warning(self):
        """Continuing would mix one gas's water terms into another's, which is
        not a degraded result but a wrong one."""
        self.assertIn("call ExceptionHandler(107)", code(FCC_MAIN))
        text = read(HANDLER)
        block = text[text.index("case(107)"):]
        block = block[: block.index("end select")]
        self.assertIn("stop 1", block)

    def test_the_message_says_how_to_recover(self):
        text = read(HANDLER)
        block = text[text.index("case(107)"):]
        block = block[: block.index("end select")]
        self.assertIn("Re-run EddyFlow-RP", block)
        self.assertIn("silently", block,
                      "the message has to say why this is fatal rather than "
                      "a warning: the misread does not announce itself")


class ASingleHygrometerProjectIsUnchanged(unittest.TestCase):
    """The whole widening is gated on a suffix being non-empty, and only the
    designated hygrometer has an empty one."""

    def test_every_emitter_skips_the_designated_entry(self):
        for path in (FLX_HDR, FLX_ROW_RP, FLX_ROW_FCC, FO_HDR, FO_ROW):
            self.assertTrue(
                re.search(r"len_trim\(w_tags\(\w+\)\) == 0\) cycle",
                          code(path)),
                "%s must skip the hygrometer that carries the bare names, or "
                "a one-hygrometer project gains a duplicate family" % path)

    def test_the_full_output_pairs_the_stability_name(self):
        """The bare column is spelled `(z-d)/L`, not zL. An unpaired name is
        the one thing a reader cannot resolve."""
        self.assertIn("',(z-d)/L' // trim(w_tags(k))", code(FO_HDR))


if __name__ == "__main__":
    unittest.main()
