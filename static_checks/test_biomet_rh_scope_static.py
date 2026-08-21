"""The biomet humidity is a quantity of its own, not an overwrite.

When a biomet RH value is available, `FluxParams` takes the site scalars from
it. It also used to write that value straight into `Stats%chi/r/d` of the
*primary* hygrometer, so on a two-hygrometer site one instrument reported
biomet and the other reported itself - and which was which followed the primary
designation, a naming choice. On CH-LAE the tell was a mixing ratio of 19.9081
that followed the primary slot between two runs and matched neither instrument:
the LI-7200 read 17.1089 and the MIRO 16.354.

One variable was doing two jobs. It was the reported column of a hygrometer,
and it was the humidity fed to the drift correction, the LI-7700 multipliers
and, through `RHO%w_at`, the WPL ratio. Overwriting served the second and
wrecked the first.

They are separate now. The biomet value has its own three quantities, reported
as `h2o_biomet_*` so the number v7.2.5 put in `h2o_mixing_ratio` is still on the
file and the two can be compared. Every hygrometer reports what it measured.
Which humidity *corrects* a gas is `moist_ref`, and the user sets it per gas.

The split this preserves, which is main's: every mean moisture quantity may
come from biomet; every covariance comes from the instrument.
"""

from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]

PARAMS = "src/src_rp/flux_params.f90"
TYPEDEF = "src/src_common/m_typedef.f90"
FO_HDR = "src/src_fcc/init_out_files.f90"
FO_HDR_RP = "src/src_rp/init_outfiles_rp.f90"
FO_ROW = "src/src_fcc/write_out_full_fcc.f90"
FO_ROW_RP = "src/src_rp/write_out_full.f90"


def read(path):
    return (ROOT / path).read_text(encoding="utf-8")


def code(path):
    """Source with comment-only lines dropped.

    The routine explains the retired overwrite in prose right where it stood,
    so a naive search matches the explanation rather than live code.
    """
    return "\n".join(ln for ln in read(path).splitlines()
                     if not ln.lstrip().startswith("!"))


def biomet_branch():
    src = code(PARAMS)
    start = src.index("if (biomet%val(bRH) > 0d0")
    return src[start: src.index("elseif (wsl >= firstGas) then", start)]


class NoHygrometerIsOverwritten(unittest.TestCase):
    def setUp(self):
        self.block = biomet_branch()

    def test_the_branch_writes_no_hygrometer_concentration(self):
        """Not for the primary, not for any of them. A hygrometer's reported
        columns are its own measurement."""
        for field in ("Stats%chi(", "Stats%r(", "Stats%d("):
            self.assertNotIn(field, self.block,
                             "%s in the biomet branch replaces a measurement "
                             "with the site value" % field)

    def test_the_site_scalars_still_come_from_biomet(self):
        """That half is main's behaviour and stays: it is what the moist-air
        correction has always used."""
        for field in ("Stats%RH = biomet%val(bRH)", "Ambient%e =",
                      "Ambient%VPD =", "RHO%w ="):
            self.assertIn(field, self.block)


class TheBiometValueIsItsOwnQuantity(unittest.TestCase):
    def test_the_fields_exist(self):
        body = code(TYPEDEF)
        for field in ("chi_biomet", "r_biomet", "d_biomet"):
            self.assertIn("real(kind = dbl) :: %s" % field, body)

    def test_they_are_computed_by_the_formulas_the_overwrite_used(self):
        block = biomet_branch()
        self.assertIn("Ambient%chi_biomet = RHO%w * Ambient%Va / MW_H2O * 1d3",
                      block)
        self.assertIn("Ambient%r_biomet   = Ambient%chi_biomet", block)

    def test_the_molar_density_is_ambient_not_a_cell(self):
        """The overwrite divided by the primary analyser's cell volume. A site
        humidity is not measured in anybody's cell."""
        block = biomet_branch()
        self.assertIn("Ambient%d_biomet   = Ambient%chi_biomet / Ambient%Va",
                      block)
        self.assertNotIn("E2Col(wsl)%Va", block)

    def test_they_are_reset_outside_the_branch(self):
        """Ambient is a module global with no per-period reset. Left to the
        branch, a period without biomet humidity reports the previous one's -
        the trap RHO%w_at carries its own comment about."""
        src = code(PARAMS)
        reset = src.index("Ambient%chi_biomet = error")
        branch = src.index("if (biomet%val(bRH) > 0d0")
        self.assertLess(reset, branch,
                        "the reset must precede the branch, or a period with "
                        "no biomet humidity carries the last one's forward")


class TheyReachTheOutput(unittest.TestCase):
    """Both halves of the full output, and both executables."""

    def test_the_header_names_them(self):
        for path in (FO_HDR, FO_HDR_RP):
            body = code(path)
            for name in ("h2o_biomet_mole_fraction",
                         "h2o_biomet_mixing_ratio",
                         "h2o_biomet_molar_density"):
                self.assertIn(name, body, path)

    def test_the_row_writes_them(self):
        self.assertIn("lEx%chi_biomet", code(FO_ROW))
        self.assertIn("Ambient%chi_biomet", code(FO_ROW_RP))

    def test_the_ex_record_carries_them(self):
        """FCC computes nothing about humidity itself; without these in the
        record its full output would have three empty columns."""
        body = code("src/src_common/read_ex_record.f90")
        self.assertIn("lEx%chi_biomet, lEx%r_biomet, lEx%d_biomet", body)

    def test_an_older_ex_file_is_refused(self):
        """They sit in the fixed part of the row, so a file written before
        them parses three fields short from there on - silently."""
        self.assertIn("H2O_BIOMET_MOLE_FRACTION",
                      code("src/src_fcc/eddyflow-fcc_main.f90"))


if __name__ == "__main__":
    unittest.main()
