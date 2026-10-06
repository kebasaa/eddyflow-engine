"""The worker pool is shared by RP and FCC, so it may use only common modules.

It used to live in src_rp/prepass_parallel.f90, which uses m_rp_global_var
and the PWB modules - neither linked into FCC. Its generic part moved to
src_common/batch_pool.f90; RP's dump formats stayed where they were.
"""

import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
POOL = (ROOT / "src/src_common/batch_pool.f90").read_text(encoding="utf-8").replace("\r\n", "\n")
PREPASS = (ROOT / "src/src_rp/prepass_parallel.f90").read_text(encoding="utf-8").replace("\r\n", "\n")
PROD = (ROOT / "src/src_rp/production_parallel.f90").read_text(encoding="utf-8").replace("\r\n", "\n")
COMMON = {p.stem for p in (ROOT / "src/src_common").glob("*.f90")}


def module_files():
    """Module name -> file stem, for every module in the source tree."""
    out = {}
    for p in (ROOT / "src").rglob("*.f90"):
        for m in re.finditer(r"^\s*module\s+(\w+)\s*$", p.read_text(encoding="utf-8", errors="replace"),
                             re.I | re.M):
            out[m.group(1).lower()] = (p.parent.name, p.stem)
    return out


class ThePoolIsCommonCode(unittest.TestCase):

    def test_it_uses_only_modules_from_src_common_or_the_os_layer(self):
        where = module_files()
        for m in re.finditer(r"^\s*use\s+(\w+)", POOL, re.I | re.M):
            folder, _ = where[m.group(1).lower()]
            self.assertIn(folder, ("src_common", "win", "posix"), m.group(1))

    def test_rp_s_dump_formats_stay_with_rp(self):
        for name in ("WriteTlagBatchDump", "MergePwbBatchDumps", "WritePfBatchDump", "BatchMagic"):
            self.assertNotIn(name, POOL)
            self.assertIn(name, PREPASS)

    def test_rp_takes_the_pool_from_it(self):
        self.assertIn("    use m_batch_pool\n", PREPASS)
        self.assertIn("    use m_batch_pool, only:", PROD)
        for name in ("subroutine LaunchChunks(", "subroutine AppendBytes(",
                     "function WorkerRoot(", "subroutine RemoveStaleWorkerRoots("):
            self.assertIn(name, POOL)
            self.assertNotIn(name, PREPASS + PROD)

    def test_stale_worker_folders_of_every_kind_are_swept(self):
        self.assertIn("index(' pd pr fx sb ', ' ' // entry(under + 1:under + 2) // ' ') == 0", POOL)


if __name__ == "__main__":
    unittest.main()
