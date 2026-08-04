import re
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def read(path):
    return (ROOT / path).read_text(encoding="utf-8")


class FluxnetHeaderGuard(unittest.TestCase):
    """Guards on the FLUXNET header and the positional reader that consumes it.

    Was a set of pytest-style module functions, which `unittest discover` - the
    runner this project actually uses - collects as zero. It had gone stale
    against the metadata split without anything reporting it.
    """

    def test_rp_fluxnet_custom_headers_are_sanitized_before_writing(self):
        source = read("src/src_rp/init_fluxnet_file_rp.f90")
        self.assertIn("usg(j) = SafeFluxnetCustomLabel(j)", source)
        self.assertIn("function SafeFluxnetCustomLabel", source)
        self.assertIn("function SanitizeFluxnetToken", source)
        self.assertIn("tmp = replace2(tmp, 'custom_', '')", source)
        self.assertIn("tmp = replace2(tmp, '_mean', '')", source)
        self.assertIn("case ('a':'z', '0':'9', '_', '-')", source)
        self.assertIn("CUSTOM_' // usg(i)(1:len_trim(usg(i)))", source)
        self.assertNotIn(
            "UserCol(j)%label(1:len_trim(UserCol(j)%label)) // '_'", source
        )

    def test_rp_fluxnet_flowrate_headers_are_model_numbered_and_sanitized(self):
        source = read("src/src_rp/init_fluxnet_file_rp.f90")
        self.assertIn("var_token = SanitizeFluxnetToken(UserCol(ordinal)%var)", source)
        self.assertIn(
            "model_token = CustomModelToken(UserCol(ordinal)%instr%model, var_token)",
            source,
        )
        self.assertIn("var_token = CustomVarToken(var_token)", source)
        self.assertIn(
            "case ('flowrate', 'co2', 'h2o', 'ch4', 'n2o', 'int_t', 'int_p')", source
        )
        self.assertIn("function CustomVarToken", source)
        self.assertIn("function CustomModelToken", source)
        self.assertIn("case ('cell_t', 'int_t_1', 'int_t_2')", source)
        self.assertIn("case ('co2', 'h2o', 'ch4', 'n2o')", source)
        self.assertIn(
            "clean_label = trim(var_token) // '_' // trim(model_token)", source
        )
        self.assertNotIn(
            "previous_tmp = SanitizeFluxnetToken(UserCol(i)%instr%model)", source
        )
        self.assertNotIn("flow_ordinal", source)

    def test_synthesized_custom_header_uses_explicit_custom_fallback(self):
        source = read("src/src_rp/init_fluxnet_file_rp.f90")
        self.assertIn(
            "clean_label = 'custom_' // trim(adjustl(ordinal_label))", source
        )
        self.assertNotIn("clean_label = trim(adjustl(ordinal_label))", source)

    def test_read_ex_record_rejects_bad_initial_read_before_fixed_skip(self):
        source = read("src/src_common/read_ex_record.f90")
        read_status_guard = (
            "if (read_status /= 0) then\n        call InvalidateRecord()"
        )
        fixed_skip = "ix = strCharIndex(dataline, ',', nMainFields)"
        self.assertLess(source.index(read_status_guard), source.index(fixed_skip))
        self.assertLess(
            source.index("if (lEx%fname == 'not_enough_data') then"),
            source.index(fixed_skip),
        )

    def test_read_ex_record_guards_every_comma_slice(self):
        """Every comma skip must reject a short line before advancing.

        This used to name six specific skips, which meant it went stale the
        moment the metadata block was split - the count it looked for no longer
        existed and the check errored rather than reporting anything useful.
        Enumerating them from the source instead cannot drift, and it covers
        the skips a widening adds without anyone remembering to list them.
        """
        source = read("src/src_common/read_ex_record.f90")
        self.assertIn("subroutine InvalidateRecord()", source)
        skips = [
            (m.start(), m.group(1))
            for m in re.finditer(
                r"ix = strCharIndex\(dataline, ',', (\w+)\)", source
            )
        ]
        self.assertGreaterEqual(
            len(skips), 15, "far fewer comma skips than expected: %s" % skips
        )
        for start, count in skips:
            with self.subTest(count=count):
                window = source[start : start + 140]
                self.assertIn("if (ix <= 0) then", window)
                self.assertIn("call InvalidateRecord()", window)

    def test_observed_header_continuation_shape_is_invalid_data_row(self):
        bad_line = (
            "2.41,-4.98,1.33,293.49,2,40,NAMEAN,CUSTOM_0,NA,NA,NA,NA,"
            "27.2397,69.99,0,42MEAN,CUSTOM_69.99,0,427.762,10.969,"
            "NUM_BIOMET_VARS,r_ET_cec,r_Fc_cec,CEC_N_VALID,CEC_N_O1"
        )
        fields = bad_line.split(",")
        self.assertLess(len(fields), 259)
        self.assertIn("NUM_BIOMET_VARS", fields)


class Vm97IsPerConfiguredGas(unittest.TestCase):
    """The VM97 flag chunk carries one field per configured gas, not four.

    Four sites state that width and they must agree, because FCC parses this
    chunk rather than echoing it: the header, RP's row writer, the reader, and
    FCC's re-emit. Two more carry the same assumption away from the row - the
    transposed string the flags are stored in, and the cospectra filter that
    reads it. A mismatch is not a crash; it shifts every later field and is
    consumed without complaint.
    """

    def test_header_names_each_gas_rather_than_listing_four(self):
        source = read("src/src_rp/init_fluxnet_file_rp.f90")
        self.assertIn(
            "trim(FluxnetLayoutTags(j)) // '_VM97_TEST'",
            source,
            "the VM97 header must be generated per configured gas",
        )
        for literal in (
            "CO2_VM97_TEST",
            "H2O_VM97_TEST",
            "CH4_VM97_TEST",
            "GS4_VM97_TEST",
        ):
            self.assertNotIn(
                literal,
                source,
                "%s is back as a header literal, so the header declares four "
                "gas columns whatever the project configures" % literal,
            )

    def test_the_transposed_flag_string_has_room_for_every_slot(self):
        """character(9) held a filler, u/v/w/ts and exactly four gases."""
        source = read("src/src_common/m_typedef.f90")
        self.assertIn(
            "character(FlagStrLen) :: vm_flags(8)",
            source,
            "vm_flags must be one character per variable; at character(9) a "
            "fifth gas is written past the end of the string",
        )

    def test_the_cospectra_filter_is_not_bounded_at_the_fourth_gas(self):
        """Otherwise gases 5+ are never filtered, whatever the tests said."""
        source = read("src/src_fcc/cospectra_qaqc.f90")
        block = source[source.index("filter_cosp_by_vm_flags") :]
        block = block[: block.index("Foken quality tests")]
        self.assertNotIn(
            "(2:9)", block, "the VM97 slices still stop at the fourth gas"
        )
        self.assertNotIn(
            "do i = co2, gas4",
            block,
            "the per-gas cospectra rejection still stops at the fourth gas",
        )
        self.assertIn("lastGas", block)

    def test_the_writers_bound_the_chunk_by_the_configured_gas_count(self):
        rp = read("src/src_rp/write_out_fluxnet.f90")
        fcc = read("src/src_fcc/write_out_fluxnet_fcc.f90")
        self.assertIn(
            "do var = u, ts + nFluxnetLayoutSlots",
            rp,
            "RP must emit one VM97 field per configured gas",
        )
        self.assertIn(
            "do var = u, ts + n_layout_gas",
            fcc,
            "FCC re-derives this chunk, so its loop has to match RP's width; "
            "echoing is not an option here",
        )

    def test_the_full_output_writers_agree_on_the_flag_width(self):
        """RP and FCC write the same column; it cannot be two widths.

        Under fcc_follows the FCC copy is the one that survives, so a narrower
        FCC write silently caps the full output at four gases - and a wider
        one carries filler digits for slots the project never configured.

        Both now cut to StatisticalFlagVars, which is also what builds the
        units row's legend, so the cell describes exactly the digits it holds.
        """
        rp = read("src/src_rp/write_out_full.f90")
        fcc = read("src/src_fcc/write_out_full_fcc.f90")
        for path, src in (("write_out_full.f90", rp),
                          ("write_out_full_fcc.f90", fcc)):
            self.assertIn("call StatisticalFlagVars(n_flag_vars", src,
                          "%s must cut the flag cells to the variables the "
                          "units row names" % path)
        self.assertIn(
            "lEx%vm_flags(i)(1 : 1 + n_flag_vars)", fcc,
            "vm_flags carries the leading filler at position one, so the "
            "slice is the same shape as RP's")
        self.assertNotIn(
            "CharHF%sr(2:FlagStrLen)", rp,
            "FlagStrLen is the array width, not the variable count: slicing "
            "to it emitted sixty-nine characters against a legend of eight")
        self.assertNotIn(
            "write(field_val, *) lEx%vm_flags(i)",
            fcc,
            "field_val is DatumLen and cannot hold a per-variable flag string; "
            "this aborts at runtime with 'End of record'",
        )

    def test_the_flag_legend_is_generated_not_spelled_out(self):
        """It was ten copies of a literal naming co2/h2o/ch4/gas4, in each of
        the four header branches - forty chances for one to drift away from
        what the row carries."""
        for path in ("src/src_rp/init_outfiles_rp.f90",
                     "src/src_fcc/init_out_files.f90"):
            src = read(path)
            self.assertNotIn("8u/v/w/ts/co2/h2o/ch4/", src,
                             "%s still spells the legend out" % path)
            self.assertIn("call StatisticalFlagVars(", src)
            self.assertIn("call TimelagFlagLegend(", src)


class Vm97ProducersCoverEveryGas(unittest.TestCase):
    """All eight VM97 tests must run for every gas, not just the first four.

    PackFlagString fills the positions beyond the count it is given with '9',
    which reads as "test not performed". Three producers passed gas4, so once
    the row carried a column for a fifth gas that column would have been a
    well-named, correctly-placed, permanently blank verdict - the failure mode
    this whole effort keeps running into.
    """

    PRODUCERS = (
        "src/src_rp/test_spike_detection_vickers_97.f90",
        "src/src_rp/test_spike_detection_mauder_13.f90",
        "src/src_rp/test_absolute_limits.f90",
        "src/src_rp/test_timelag.f90",
    )

    def test_no_producer_packs_only_the_historical_slots(self):
        for path in self.PRODUCERS:
            with self.subTest(producer=path):
                source = read(path)
                self.assertNotRegex(
                    source,
                    r"call PackFlagString\(.*,\s*gas4\s*,",
                    "packs only up to the fourth gas, so every later slot "
                    "reports 'not performed' whatever the test found",
                )
                self.assertRegex(
                    source, r"call PackFlagString\(.*,\s*GHGNumVar\s*,"
                )

    def test_the_time_lag_flags_are_not_an_integer(self):
        """The base-10 packing is what bounded TestTimeLag at four gases.

        `90000 + sum(hflags(j) * 10**(4 - j))` overflows a 32-bit integer past
        about nine variables, so the flag arrays could not grow - and because
        the loop was sized to them, neither could the test. A fifth gas's lag
        was never compared against its default while the flag column read as
        "not performed".

        The two row writers sliced the string from its end, which only found
        the right four digits because int2char right-aligned them there. They
        must slice from the first gas now, and against the same legend helper
        the header uses, or the cell and its units row disagree.
        """
        #> Comment lines dropped, or the note explaining the removed packing
        #> reads as the packing itself.
        source = "\n".join(
            ln for ln in read("src/src_rp/test_timelag.f90").splitlines()
            if not ln.lstrip().startswith("!"))
        self.assertNotIn("IntHF%tl", source)
        self.assertNotRegex(source, r"10\*\*\(4 - j\)")

        typedefs = read("src/src_common/m_typedef.f90")
        self.assertNotRegex(
            typedefs, r"type :: RSIntFlagType.*?integer :: tl.*?end type",
            "IntHF%tl is back; the per-variable flags must stay strings")
        for field in ("vm_tlag_hf", "vm_tlag_sf"):
            self.assertIn("character(FlagStrLen) :: %s" % field, typedefs,
                          "%s must grow with the gas count or the ex file "
                          "truncates it past eight gases" % field)

        for path in ("src/src_rp/write_out_full.f90",
                     "src/src_rp/write_out_fluxnet.f90"):
            writer = read(path)
            self.assertNotIn("CharHF%tl(FlagStrLen-3:FlagStrLen)", writer)
            self.assertIn(
                "CharHF%tl(firstGas + 1:firstGas + n_tl_vars)", writer,
                "%s must cut the cell to the gases the legend names" % path)
            self.assertIn("TimelagFlagLegend(n_tl_vars", writer)

    def test_the_spike_tests_evaluate_every_slot(self):
        for path in self.PRODUCERS[:2]:
            with self.subTest(producer=path):
                source = read(path)
                self.assertIn("do j = u, lastGas", source)
                self.assertNotIn("do j = u, gas4", source)

    def test_absolute_limits_is_one_loop_rather_than_four_copies(self):
        """Four copies is how gases past the fourth were skipped entirely."""
        source = read("src/src_rp/test_absolute_limits.f90")
        self.assertIn("do i = firstGas, lastGas", source)
        self.assertIn("integer :: hflags(GHGNumVar)", source)
        for slot in ("co2", "ch4", "gas4"):
            self.assertNotIn(
                "E2Col(%s)%%present" % slot,
                source,
                "a per-slot copy of the absolute-limits arm is back; the loop "
                "must cover every configured gas",
            )

    def test_a_gas_without_limits_is_not_tested_against_zero(self):
        """Absent limits mean "not performed", not "everything is invalid".

        Only the four historical slots take their limits from fixed project
        keys; past those they come from the per-gas records, and a project that
        names a gas without them leaves the pair at 0/0. Testing against that
        rejects every value - and because this routine filters on the same
        pass, it replaced the gas's entire series with the error code. The gas
        then reached the flux code empty and EliminateCorruptedVariables
        dropped it, so a fifth gas silently produced no concentrations or
        fluxes at all while still reporting screening flags.

        Widening the loop to every slot is what exposed this; the four-gas
        regression cannot see it, because four gases always have limits.
        """
        source = read("src/src_rp/test_absolute_limits.f90")
        body = source[source.index("do i = firstGas, lastGas") :]
        guard = "if (al%gas_max(i) <= al%gas_min(i)) then"
        self.assertIn(
            guard,
            body,
            "a gas whose limits were never configured is tested against 0/0, "
            "which wipes its whole time series",
        )
        # The guard has to precede the counting, or the damage is already done.
        self.assertLess(
            body.index(guard),
            body.index("Essentials%al_s(i) = count("),
            "the unconfigured-limits guard must run before the test itself",
        )

    def test_water_is_singled_out_by_species_not_by_slot_number(self):
        """The two values that differ between gases both belong to water.

        Water's mole fraction is reported in mmol mol-1 where the other gases
        use umol mol-1, so it needs a different molar-density scale and a
        different rough-outlier ceiling. That is a fact about the species.

        This assertion used to require `i == h2o`, which is the historical
        slot rather than the species - so it pinned the very thing its own
        docstring said was wrong. A second hygrometer sits well past slot 6
        and was given the trace-gas scale and the trace-gas ceiling, its
        readings compared against limits three orders of magnitude out.

        FilterDatasetForPhysicalThresholds consults the same
        al%gas_min/gas_max pair on the same pass, so both must settle water
        the same way or they disagree about one gas.
        """
        for path, opener in (
            ("src/src_rp/test_absolute_limits.f90", "do i = firstGas, lastGas"),
            ("src/src_rp/filter_dataset_for_physical_thresholds.f90",
             "do gas = firstGas, lastGas"),
        ):
            body = read(path)
            body = body[body.index(opener):]
            self.assertNotRegex(
                body, r"==\s*h2o\b",
                "%s must ask GasSlotIsWater, not compare against the h2o "
                "slot" % path)
            self.assertIn("GasSlotIsWater(", body,
                          "%s must single water out by species" % path)
        body = read("src/src_rp/test_absolute_limits.f90")
        body = body[body.index("do i = firstGas, lastGas"):]
        self.assertIn("dens_scale = StdVair", body)
        self.assertIn("dens_scale = StdVair * 1d3", body)


class FluxnetGasScalesAreBySpecies(unittest.TestCase):
    """No column scale may be chosen by which slot a gas occupies.

    The FLUXNET row carried a hard-coded x1000 applied only when the slot was
    ch4 or gas4. That names a position, not a species, so the same gas was
    reported in nmol mol-1 from slot 7 and in umol mol-1 from slot 9 - and a
    project with CO2 in slot four had it multiplied by a thousand.
    """

    WRITERS = (
        "src/src_rp/write_out_fluxnet.f90",
        "src/src_fcc/write_out_fluxnet_fcc.f90",
    )

    def test_no_gain_is_chosen_by_slot(self):
        for path in self.WRITERS:
            with self.subTest(writer=path):
                source = read(path)
                for line_no, line in enumerate(source.splitlines(), 1):
                    if "gain=" not in line:
                        continue
                    self.assertNotRegex(
                        line,
                        r"gas\s*==\s*(ch4|gas4)|gas\s*>=\s*ch4",
                        "line %d picks a scale from the slot" % line_no,
                    )
                # The instrument-geometry gains (cm, mm, l/min) are unit
                # conversions of metadata, not of gas quantities, and stay.
                gas_gains = [
                    l for l in source.splitlines()
                    if "gain=1d3" in l or "gain=1d6" in l
                ]
                self.assertTrue(
                    all("Instr%" in l or "gas_instr" in l for l in gas_gains),
                    "a literal per-gas gain is back: %s" % gas_gains,
                )

    def test_the_scale_comes_from_the_species_function(self):
        for path in self.WRITERS:
            with self.subTest(writer=path):
                source = read(path)
                self.assertIn("gain=FluxnetGasScale(", source)
                self.assertIn("gain=FluxnetGasAdvScale(", source)

    def test_the_reader_inverts_with_the_same_function(self):
        """A writer and an inverse derived separately drift silently.

        The inverse was a literal 1d-3 on two slots. If the writer's factor
        ever stops being 1000 for those slots - which is exactly what making
        it per-species does - a literal inverse un-scales by a number that was
        never applied, and nothing reports it.
        """
        source = read("src/src_common/read_ex_record.f90")
        block = source[source.index("!> Units adjustments") :]
        block = block[: block.index("!> Variances were actually read as")]
        # Every quantity the writer scales must be un-scaled here. Requiring
        # merely one use would pass with three of the four left as literals.
        # ...over every configured gas, not the historical four. Stopping at
        # gas4 left slot 5+ un-inverted, so FCC scaled it a second time and
        # N2O came out at 339273 nmol/mol instead of 339.
        self.assertIn("do gas = co2, ts + min(EddyFlowProj%gas_num", block)
        self.assertEqual(
            block.count("/ FluxnetGasScale(gas)"),
            4,
            "the inverse must cover all four scaled quantities - Flux0, "
            "rand_uncer, r and chi - each through the same function",
        )
        for slot in ("ch4", "gas4"):
            # Built by concatenation: '%r' in the Fortran text is a Python
            # conversion specifier if this goes through % formatting.
            self.assertNotIn(
                "lEx%r(" + slot + ") * 1d-3",
                block,
                "the inverse is a literal again; it must come from the same "
                "function the writer uses",
            )

    def test_the_species_rule_matches_the_fluxnet_bases(self):
        """CO2 umol, H2O mmol, everything else nmol - and nothing by slot."""
        source = read("src/src_common/gas4_output_units.f90")
        body = source[source.index("function FluxnetGasScale") :]
        body = body[: body.index("end function FluxnetGasScale")]
        self.assertIn("EddyFlowProj%gas(rec)%var", body)
        self.assertIn("case ('CO2', 'H2O')", body)
        for slot in ("ch4", "gas4", "co2 ", "h2o "):
            self.assertNotIn(
                "gas_slot == %s" % slot.strip(),
                body,
                "the scale is being chosen by slot again",
            )


class TrailingBlocksArePerConfiguredGas(unittest.TestCase):
    """No trailing block may still spell out a fixed four-gas quadruple.

    Each of these families is one column per gas. While they were written out
    as CO2/H2O/CH4/GS4 literals a fifth gas had no column at all, and the
    header and the row could disagree about the width without any error.
    """

    HEADER = "src/src_rp/init_fluxnet_file_rp.f90"

    #: family suffix -> shape. 'var' families cover u,v,w,ts then the gases;
    #: 'flux' families cover the flux names then the gases, and take their
    #: per-gas prefix from FluxnetFluxTag because CO2's flux column is FC.
    FAMILIES = {
        "_VM97_TEST": "var",
        "_LGD": "var",
        "_KID": "var",
        "_ZCD": "var",
        "_CORRDIFF": "flux",
        "_NSR": "flux",
        "_SS": "flux",
        "_SS_TEST": "flux",
        "_SSITC_TEST": "flux",
    }

    def test_no_family_is_written_as_a_four_gas_literal(self):
        source = read(self.HEADER)
        for suffix in sorted(self.FAMILIES):
            for slot in ("CO2", "CH4", "GS4"):
                with self.subTest(family=suffix, slot=slot):
                    self.assertNotIn(
                        "%s%s" % (slot, suffix),
                        source,
                        "%s%s is a header literal again, so this family is "
                        "back to a fixed four gases" % (slot, suffix),
                    )

    def test_every_family_is_generated_from_the_layout_list(self):
        source = read(self.HEADER)
        for suffix, shape in sorted(self.FAMILIES.items()):
            with self.subTest(family=suffix):
                if shape == "var":
                    # The three plain variable families share AddVariableFamily.
                    generated = (
                        "call AddVariableFamily('%s')" % suffix in source
                        or "trim(FluxnetLayoutTags(j)) // '%s'" % suffix in source
                    )
                else:
                    generated = (
                        "trim(FluxnetFluxTag(j)) // '%s'" % suffix in source
                    )
                self.assertTrue(
                    generated,
                    "%s is not generated per configured gas" % suffix,
                )

    def test_the_flux_prefix_keeps_carbon_dioxides_historical_name(self):
        """FC, not FCO2 - these columns are named for the flux.

        Deriving the prefix from the species tag alone would silently rename
        FC_SS to FCO2_SS and every other CO2 flux column with it.
        """
        source = read(self.HEADER)
        body = source[source.index("function FluxnetFluxTag") :]
        body = body[: body.index("end function FluxnetFluxTag")]
        self.assertIn("== co2", body)
        self.assertIn("tag = 'FC'", body)
        self.assertIn("tag = 'F' // trim(FluxnetLayoutTags(layout_index))", body)

    def test_the_per_gas_quality_types_are_slot_indexed(self):
        """QCType and ExType carried four named members apiece.

        Four scalars cannot hold a fifth gas, so the steady-state statistic and
        the Foken flags stopped at the fourth slot no matter how wide the row
        became.
        """
        source = read("src/src_common/m_typedef.f90")
        qc = source[source.index("type :: QCType") :]
        qc = qc[: qc.index("end type QCType")]
        ex = source[source.index("type ExType") :]
        ex = ex[: ex.index("end type ExType")]

        self.assertIn("integer :: w_gas(GHGNumVar)", qc)
        self.assertIn("integer :: gas(GHGNumVar)", qc)
        self.assertIn("real(kind = dbl) :: F_SS(GHGNumVar)", ex)

        # Scoped to the type bodies: w_co2 and friends also name the covariance
        # index constants, which are aliases like co2 itself and stay put. It
        # is only as *members* that they freeze a per-gas quantity at four.
        for member in ("w_co2", "w_h2o", "w_ch4", "w_gas4"):
            self.assertNotIn(
                ":: %s" % member,
                qc,
                "%s is back as a QCType member; four scalars cannot hold a "
                "fifth gas" % member,
            )
        for member in ("FC_SS", "FH2O_SS", "FCH4_SS", "FGS4_SS"):
            self.assertNotIn(
                ":: %s" % member,
                ex,
                "%s is back as an ExType member; the ex file could then carry "
                "the steady-state statistic for four gases only" % member,
            )


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
