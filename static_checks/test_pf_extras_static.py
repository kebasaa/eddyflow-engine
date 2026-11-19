"""Planar-fit extras after innFLUX: moving windows and a tilt bootstrap.

innFLUX (Striednig et al. 2020, Sect. 2.3.1, Fig. 1; innFLUX_step0.m) fits the
planar fit of Wilczak et al. (2001) to the periods within +/-15 degrees of each
wind direction in 1-degree steps - 360 overlapping windows - and bootstraps the
tilt angles' uncertainty. Both are opt-in here, both off by default:

  pf_mw_num, pf_mw_width   RP.SNTags(440), (441)  windows instead of sectors
  pf_boot_n, pf_boot_seed  RP.SNTags(442), (443)  bootstrap of alpha and beta

When used, they are written into the planar-fit file, and a reused file
(pf_mode 0) applies its own windows whatever the project says. With both off
the file is written exactly as before.

Part of the EddyFlow engine's static checks.
"""

from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]


def read(*parts):
    return ROOT.joinpath(*parts).read_text(encoding="utf-8", errors="replace")


def code(*parts):
    return "\n".join(ln for ln in read(*parts).splitlines()
                     if not ln.lstrip().startswith("!"))


TAGS = read("src", "src_rp", "m_rp_global_var.f90")
TYPEDEF = code("src", "src_common", "m_typedef.f90")
COMMON = code("src", "src_common", "m_common_global_var.f90")
READER = code("src", "src_rp", "read_ini_rp.f90")
WRITER = code("src", "src_rp", "write_out_planar_fit.f90")
PFREAD = code("src", "src_rp", "read_planar_fit_file.f90")
SECTOR = code("src", "src_rp", "wind_sector.f90")
SORT = code("src", "src_rp", "sort_wind_by_sector.f90")
MAIN = code("src", "src_rp", "eddyflow-rp_main.f90")
BOOT = read("src", "src_common", "m_pfboot.f90")
PAR = code("src", "src_rp", "production_parallel.f90")


class TheKeys(unittest.TestCase):
    def test_slots(self):
        for idx, key in ((440, "pf_mw_num"), (441, "pf_mw_width"),
                         (442, "pf_boot_n"), (443, "pf_boot_seed")):
            self.assertRegex(TAGS, rf"SNTags\({idx}\)%Label\s*/\s*'{key}'")

    def test_off_unless_stated(self):
        self.assertIn("PFSetup%mw_num = 0", READER)
        self.assertIn("PFSetup%boot_n = 0", READER)

    def test_windows_replace_the_sectors(self):
        self.assertIn("PFSetup%num_sec = PFSetup%mw_num", READER)


class Capacity(unittest.TestCase):
    def test_rotations_sized_for_360_windows(self):
        self.assertIn("integer, parameter :: MaxNumPfSect = 360", TYPEDEF)
        self.assertIn("PFMat(3, 3, MaxNumPfSect)", COMMON)
        self.assertIn("PFb(3, MaxNumPfSect)", COMMON)
        self.assertIn("CtxPFMat(3, 3, MaxNumPfSect)", PAR)

    def test_user_sectors_keep_their_limit(self):
        self.assertIn("integer, parameter :: MaxNumWSect = 36", TYPEDEF)
        self.assertIn("if (PFSetup%mw_num == 0) then\n                    if (PFSetup%wsect_exclude(sec)) then",
                      MAIN)


class Windows(unittest.TestCase):
    def test_overlapping_membership(self):
        self.assertIn("dd = abs(modulo(twd - centre + 180d0, 360d0) - 180d0)", SORT)
        self.assertIn("if (dd < 0.5d0 * PFSetup%mw_width) then", SORT)

    def test_nearest_window_at_application(self):
        self.assertIn("n = modulo(nint(twdir * dble(PFSetup%mw_num) / 360d0), PFSetup%mw_num) + 1",
                      SECTOR)


class TheFile(unittest.TestCase):
    def test_moving_window_file_is_self_describing(self):
        self.assertIn("'Moving_window_sectors: '", WRITER)
        self.assertIn("index(line, 'Moving_window_sectors:')", PFREAD)
        self.assertIn("PFSetup%num_sec = nmw", PFREAD)

    def test_reused_file_wins_over_the_project(self):
        self.assertIn("the file''s are applied", PFREAD)

    def test_bootstrap_after_what_an_older_reader_reads(self):
        rot = WRITER.index("'Rotation matrices'")
        boot = WRITER.index("if (PFSetup%boot_n > 0) then")
        self.assertLess(rot, boot)

    def test_legacy_close_unchanged(self):
        """The legacy writer has always closed udf, not upf; left as it was."""
        legacy = WRITER[:WRITER.index("contains")]
        self.assertIn("close(udf)", legacy)


class Bootstrap(unittest.TestCase):
    def test_pure_and_seeded_per_sector(self):
        self.assertIn("pure subroutine BootstrapTilt", BOOT)
        self.assertIn("PFSetup%boot_seed + sec", MAIN)

    def test_innflux_angles(self):
        self.assertIn("alpha = asin(p31) * rad2deg", BOOT)
        self.assertIn("beta = asin(-p32 / sqrt(p32 ** 2 + p33 ** 2)) * rad2deg", BOOT)

    def test_credited(self):
        self.assertIn("Striednig et al. 2020", BOOT)
        self.assertIn("Wilczak et al. (2001)", BOOT)


if __name__ == "__main__":
    unittest.main()
