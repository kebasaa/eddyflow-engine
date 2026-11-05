"""Carbonyl sulfide may sort the octants, experimentally, and its stomatal
component is flagged when it comes out as emission.

Zahn et al. (2022) sort moist ejections by the sign of c': air that lost CO2
passed through stomata. Stomata take carbonyl sulfide up the same way, so a
stated pairing may name a COS record as its octant gas instead. It is not the
default and it is not a drop-in: many soils take COS up as well, and soil
uptake carried in a moist ejection then sorts as stomatal. With CO2 octants a
COS target needs no such assumption, because its own sign never enters the
sort. Both are kept; the second is the one the interface calls experimental.

Before this, the resolver looked the octant gas up as CO2 and nothing else, and
a pairing that named any other record was silently given the designated CO2.

The flag: leaves are not known to emit COS while their stomata work, so a
stomatal component above zero is noise, soil uptake misattributed by COS
octants, or a senescent canopy. It is reported either way; cec_sign_<stem>
says which periods to look at.
"""

from pathlib import Path
import re
import unittest

from test_cec_reference_static import (NORMAL, apply_partition, carbon_status,
                                       target_sums)


ROOT = Path(__file__).resolve().parents[1]

RESOLVER = "src/src_common/gas_slot_resolution.f90"
HANDLER = "src/src_common/exception_handler.f90"
PARSER = "src/src_common/write_processing_project_variables.f90"
CEC = "src/src_common/m_cec.f90"
READER = "src/src_common/read_ex_record.f90"
ROW_WRITERS = ("src/src_rp/write_out_fluxnet.f90",
               "src/src_rp/write_out_fluxnet_only_biomet.f90",
               "src/src_fcc/write_out_fluxnet_fcc.f90")


def read(relative):
    return (ROOT / relative).read_text(encoding="utf-8", errors="replace")


def code(relative):
    """Source with whole-line comments removed."""
    return "\n".join(ln for ln in read(relative).splitlines()
                     if not ln.lstrip().startswith("!"))


def between(source, start, end):
    block = source[source.index(start):]
    return block[:block.index(end)]


def stomatal_sign_flag(stomatal, error=-9999):
    """What cec_sign_<stem> reports: 1 for emission, 0 for uptake or none."""
    if stomatal is None:
        return error
    return 1 if stomatal > 0 else 0


class TheOctantGasMayBeCarbonylSulfide(unittest.TestCase):
    def setUp(self):
        self.pairs = between(code(RESOLVER), "subroutine CecPairs(",
                             "end subroutine CecPairs")

    def test_a_stated_pairing_takes_co2_first_then_cos(self):
        co2 = "slot = CecSlotOfRecord(rec, 'CO2')"
        cos = "if (slot == 0) slot = CecSlotOfRecord(rec, 'COS')"
        self.assertIn(co2, self.pairs)
        self.assertIn(cos, self.pairs)
        self.assertLess(self.pairs.index(co2), self.pairs.index(cos))

    def test_any_other_species_falls_back_and_says_so(self):
        #> Silently partitioning against a gas the user did not choose is what
        #> this replaced. An unstated record is still resolved quietly - that
        #> is what not stating means - so only a record with a column warns.
        self.assertIn("if (CecRecordHasColumn(rec) .and. .not. octantWarned(i)) then",
                      self.pairs)
        self.assertIn("call ExceptionHandler(123)", self.pairs)
        self.assertLess(self.pairs.index("call ExceptionHandler(123)"),
                        self.pairs.index("slot = PrimaryCarbonSlot()"))
        self.assertIn("case(123)", read(HANDLER))
        self.assertIn("Warning(123)", read(HANDLER))

    def test_the_water_is_still_looked_up_as_water(self):
        self.assertIn("CecSlotOfRecord(EddyFlowProj%cec_pair(i)%water, 'H2O')",
                      self.pairs)

    def test_the_automatic_default_stays_co2(self):
        auto = between(code(RESOLVER), "subroutine AutoCecPairs(",
                       "end subroutine AutoCecPairs")
        self.assertIn("if (trim(adjustl(species)) /= 'CO2') cycle", auto)
        self.assertNotIn("'COS'", auto)

    def test_a_cos_pairing_carries_its_species_in_its_suffix(self):
        #> So E_cec_miro_mga4_6 and E_cec_miro_mga4_6_cos sit side by side
        #> rather than as a first and second occurrence of one partition. A CO2
        #> pairing's suffix is unchanged, which keeps every existing column.
        self.assertIn("if (len_trim(species) > 0 .and. trim(species) /= 'CO2') then",
                      self.pairs)
        self.assertIn("CecPairBaseTag = trim(CecPairBaseTag) // '_' // trim(species)",
                      self.pairs)


class TheColumnsFollowTheOctantGas(unittest.TestCase):
    def setUp(self):
        resolver = code(RESOLVER)
        self.columns = between(resolver, "subroutine CecOutputColumns(",
                               "end subroutine CecOutputColumns")
        self.rows = between(resolver, "subroutine CecRowValues(",
                            "end subroutine CecRowValues")
        self.flag = between(resolver, "logical function CecFlagsStomatalSign(",
                            "end function CecFlagsStomatalSign")

    def test_respiration_and_photosynthesis_are_co2_words(self):
        # The species is looked up once per target, outside the compound
        # condition, which gfortran need not evaluate in full.
        self.assertIn("k_species = CecSlotSpecies(slots(k))", self.columns)
        branch = ("else if (k == cecTargetCarbon .and. "
                  "trim(k_species) == 'CO2') then")
        self.assertIn(branch, self.columns)
        self.assertLess(self.columns.index(branch),
                        self.columns.index("call Emit('Reco_cec' // trim(tag)"))

    def test_any_other_octant_gas_is_named_like_an_extra(self):
        #> And still honours the pairing's choice of what to report.
        self.assertIn("if (k == cecTargetCarbon .and. pair%meth /= 1 "
                      ".and. pair%meth /= 3) cycle", self.columns)

    def test_only_carbonyl_sulfide_is_flagged(self):
        self.assertIn("CecFlagsStomatalSign = trim(CecSlotSpecies(slot)) == 'COS'",
                      self.flag)

    def test_the_header_and_the_row_flag_under_the_same_condition(self):
        """The header names the flag in the one branch both the non-CO2
        octant gas and the extras take; the row has those as two branches with
        identical bodies, so it asks twice. Water never gets one."""
        self.assertEqual(self.columns.count("if (CecFlagsStomatalSign(slots(k)))"), 1)
        self.assertIn("call Emit('cec_sign_' // trim(stem) // trim(tag), '[#]')",
                      self.columns)
        self.assertEqual(
            self.rows.count("if (CecFlagsStomatalSign(slots(k))) "
                            "call EmitSign(flux%comp(k)%stomatal)"), 2)
        water = between(self.rows, "if (k == cecTargetWater) then",
                        "else if (k == cecTargetCarbon) then")
        self.assertNotIn("EmitSign", water)

    def test_the_flag_is_a_count_with_the_error_code_for_nothing(self):
        self.assertIn("call EmitCount(merge(1, 0, stomatal > 0d0))", self.rows)
        self.assertIn("call EmitCount(nint(error))", self.rows)

    def test_the_flag_reads_the_component_not_the_total(self):
        self.assertEqual(stomatal_sign_flag(-3.0), 0)
        self.assertEqual(stomatal_sign_flag(0.0), 0)
        self.assertEqual(stomatal_sign_flag(0.4), 1)
        self.assertEqual(stomatal_sign_flag(None), -9999)


class ExtrasOnly(unittest.TestCase):
    """meth 4: a pairing kept only to sort the parcels its extras are split
    over, reporting neither of its own two gases."""

    def test_it_is_read_rather_than_reset(self):
        parser = read(PARSER)
        self.assertIn(".or. EddyFlowProj%cec_pair(i)%meth > 4) &", parser)
        self.assertNotIn(".or. EddyFlowProj%cec_pair(i)%meth > 3) &", parser)

    def test_it_wants_neither_of_the_pairs_own_gases(self):
        wanted = between(code(CEC), "logical function CecTargetIsWanted(",
                         "end function CecTargetIsWanted")
        self.assertIn("CecTargetIsWanted = meth == 1 .or. meth == 2", wanted)
        self.assertIn("CecTargetIsWanted = meth == 1 .or. meth == 3", wanted)
        self.assertIn("CecTargetIsWanted = .true.", wanted)

    def test_the_format_that_carries_it_is_one_this_engine_reads(self):
        version = re.search(r"MaxSupportedIniVer = '([\d.]+)'",
                            read("src/src_common/m_common_global_var.f90"))
        self.assertIsNotNone(version)
        self.assertGreaterEqual(tuple(int(p) for p in version.group(1).split(".")),
                                (5, 2, 0))


class TheEssentialsRowHoldsAFullPairing(unittest.TestCase):
    def test_the_buffer_is_as_wide_as_the_widest_pairing(self):
        """A pairing with four extras is six targets: 9 + 7 x 6 = 51 fields.

        The buffers were sized six per target, 45, and the emitters stop
        silently at the end of the buffer - so that pairing's row came out
        shorter than its header and every field after it was misread. Carbonyl
        sulfide with three companions is exactly the pairing that reaches it.
        """
        reader = read(READER)
        fixed = int(re.search(r"nCecPairFixedFields\s*=\s*(\d+)", reader).group(1))
        per_target = int(re.search(r"nCecTargetFields\s*=\s*(\d+)", reader).group(1))
        for path in ROW_WRITERS:
            source = read(path)
            for name in ("cec_values", "cec_is_int"):
                m = re.search(name + r"\(MaxNumCecTargets \* (\d+) \+ (\d+)\)", source)
                self.assertIsNotNone(m, "%s: %s is not sized" % (path, name))
                self.assertEqual((int(m.group(1)), int(m.group(2))),
                                 (per_target, fixed),
                                 "%s: %s cannot hold a full pairing" % (path, name))


class WhatTheOctantGasMeansForCos(unittest.TestCase):
    """The reason CO2 stays the default, and when COS octants agree with it.

    Four moist ejections: two off the canopy (CO2 drawn down), two off the
    ground (CO2 from respiration). The canopy takes COS up in both cases.
    """
    W = [1.0, 1.0, 1.0, 1.0]
    Q = [0.5, 0.5, 0.5, 0.5]
    C = [-1.0, -1.0, 1.0, 1.0]

    def split(self, sort_by, cos, total):
        f_o1, f_o2 = target_sums(self.W, sort_by, self.Q, cos)
        if f_o1 == 0:
            return f_o1, f_o2, (0.0, total)
        ratio = f_o1 / f_o2
        status = carbon_status(ratio)
        self.assertEqual(status, NORMAL)
        return f_o1, f_o2, apply_partition(status, ratio, total)

    def test_co2_octants_give_soil_uptake_to_the_ground(self):
        cos = [-2.0, -2.0, -1.0, -1.0]
        _, _, (nonstomatal, stomatal) = self.split(self.C, cos, -6.0)
        self.assertLess(nonstomatal, 0)
        self.assertLess(stomatal, 0)
        self.assertAlmostEqual(nonstomatal, -2.0)
        self.assertAlmostEqual(stomatal, -4.0)

    def test_cos_octants_give_soil_uptake_to_the_canopy(self):
        #> Every parcel lost COS, so every parcel sorts as stomatal: the soil's
        #> share is counted as plant uptake. This is the experimental caveat.
        cos = [-2.0, -2.0, -1.0, -1.0]
        f_o1, _, (nonstomatal, stomatal) = self.split(cos, cos, -6.0)
        self.assertEqual(f_o1, 0)
        self.assertEqual(nonstomatal, 0.0)
        self.assertAlmostEqual(stomatal, -6.0)

    def test_where_the_soil_emits_the_two_agree(self):
        cos = [-2.0, -2.0, 1.0, 1.0]
        _, _, by_co2 = self.split(self.C, cos, -2.0)
        _, _, by_cos = self.split(cos, cos, -2.0)
        self.assertAlmostEqual(by_co2[0], by_cos[0])
        self.assertAlmostEqual(by_co2[1], by_cos[1])
        self.assertGreater(by_co2[0], 0)
        self.assertLess(by_co2[1], 0)


if __name__ == "__main__":
    unittest.main()
