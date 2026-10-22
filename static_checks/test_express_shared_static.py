"""Express mode is applied where both RP and FCC read the project file.

`ConfigureForExpress` used to be RP's alone. It set WPL on, the spectral
correction, the QC and footprint methods and the output files, and FCC, which
reads the same project file and computes the final fluxes, took every one of
them as written. An express project carrying wpl_meth=0 therefore had RP
convert to mixing ratios point by point with WPL on, and FCC finish the fluxes
with WPL off: on base_ep_licor, LE came out ~22x off EddyPro's, which ignores
the 0 in express mode.

The settings FCC reads now live in `ConfigureProjectForExpress`, called at the
end of `WriteProcessingProjectVariables`, the reader both programs share, and
before FCC derives its spectral-assessment switches from them. RP's own routine
keeps only what lives in RP's types, and calls the shared one too.
"""

from pathlib import Path
import re
import unittest


ROOT = Path(__file__).resolve().parents[1]

SHARED = "src/src_common/configure_project_for_express.f90"
RP_EXPRESS = "src/src_rp/configure_for_express.f90"
READER = "src/src_common/write_processing_project_variables.f90"
FCC_INI = "src/src_fcc/read_ini_fcc.f90"

#> Shared-module objects RP's routine could set and FCC could then read.
COMMON_OBJECTS = ("EddyFlowProj", "Meth", "RUsetup", "Dir", "AuxFile")

#> Set by RP's routine alone on purpose, though FCC reads it. The metadata
#> file holds no value FCC computes, and FCC's writer lays it out differently
#> from RP's (no DOY, SI units, rounded coordinates): forcing it on in FCC
#> would only have FCC overwrite the file the express RP run already wrote.
RP_ONLY_ON_PURPOSE = {("EddyFlowProj", "out_md")}


def read(path):
    return (ROOT / path).read_text(encoding="utf-8")


def code(path):
    """Source with comments dropped."""
    out = []
    for ln in read(path).splitlines():
        if ln.lstrip().startswith("!"):
            continue
        out.append(ln.split("!", 1)[0] if "'" not in ln else ln)
    return "\n".join(out)


def routine(path, name):
    body = code(path)
    start = re.search(r"^\s*subroutine %s\b" % name, body, re.M | re.I)
    assert start, "%s: no subroutine %s" % (path, name)
    end = body.index("end subroutine %s" % name, start.start())
    return body[start.start():end]


def assigned(body):
    """{(object, field)} assigned in body, e.g. ('EddyFlowProj', 'wpl')."""
    return set(re.findall(r"^\s*(\w+)%(\w+)(?:\([^)]*\))?\s*=(?!=)",
                          body, re.M))


def fcc_code():
    return "\n".join(code(p.relative_to(ROOT).as_posix())
                     for p in sorted((ROOT / "src/src_fcc").glob("*.f90")))


class TheSharedRoutineCarriesWhatFccReads(unittest.TestCase):

    def test_it_fixes_wpl_and_the_corrections(self):
        got = assigned(routine(SHARED, "ConfigureProjectForExpress"))
        for need in (("EddyFlowProj", "wpl"), ("EddyFlowProj", "hf_meth"),
                     ("EddyFlowProj", "cosp_model"), ("Meth", "qcflag"),
                     ("Meth", "foot"), ("RUsetup", "meth"),
                     ("EddyFlowProj", "out_full"),
                     ("EddyFlowProj", "out_avrg_cosp"),
                     ("EddyFlowProj", "out_avrg_spec")):
            self.assertIn(need, got, "%s%%%s" % need)
        body = routine(SHARED, "ConfigureProjectForExpress")
        self.assertRegex(body, r"EddyFlowProj%wpl\s*=\s*\.true\.")

    def test_rp_leaves_nothing_fcc_reads_to_itself(self):
        """Anything RP's routine still sets on its own, FCC must not read.

        Otherwise the two programs run one project under two configurations
        again - which is the defect, whatever the field.
        """
        fcc = fcc_code()
        own = assigned(routine(RP_EXPRESS, "ConfigureForExpress"))
        for obj, field in sorted(own):
            if obj not in COMMON_OBJECTS or (obj, field) in RP_ONLY_ON_PURPOSE:
                continue
            self.assertNotRegex(
                fcc, r"\b%s%%%s\b" % (obj, field),
                "%s%%%s is set by RP's ConfigureForExpress alone but FCC "
                "reads it; move it to ConfigureProjectForExpress"
                % (obj, field))

    def test_rp_still_applies_it(self):
        body = routine(RP_EXPRESS, "ConfigureForExpress")
        self.assertRegex(body, r"call\s+ConfigureProjectForExpress\(\)")


class BothProgramsGetIt(unittest.TestCase):

    def test_the_shared_reader_calls_it_for_express(self):
        body = code(READER)
        main = body[:body.index("\ncontains\n")]
        self.assertRegex(
            main,
            r"if\s*\(EddyFlowProj%run_mode\s*==\s*'express'\)\s*"
            r"call\s+ConfigureProjectForExpress\(\)")

    def test_after_the_tags_it_overrides_are_read(self):
        """A tag read after the call would put the file's value back."""
        body = code(READER)
        main = body[:body.index("\ncontains\n")]
        call = main.index("call ConfigureProjectForExpress()")
        for obj, field in assigned(routine(SHARED,
                                           "ConfigureProjectForExpress")):
            later = re.search(r"^\s*%s%%%s\s*=(?!=)" % (obj, field),
                              main[call:], re.M)
            self.assertIsNone(later, "%s%%%s is read from the file after "
                              "the express override" % (obj, field))

    def test_fcc_reads_the_project_before_deriving_from_it(self):
        """WriteVariablesFCC turns hf_meth and out_avrg_* into FCC's
        spectral-assessment switches; it must see the express values."""
        body = routine(FCC_INI, "ReadIniFCC")
        self.assertLess(body.index("call WriteProcessingProjectVariables()"),
                        body.index("call WriteVariablesFCC()"))


if __name__ == "__main__":
    unittest.main()
