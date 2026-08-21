"""A path the project states, which is not there, must stop the run.

Five settings used to answer a missing path by quietly computing something
else: a missing planar fit file became double rotation, a missing time-lag file
became covariance maximisation, and a missing assessment file or co-spectra
directory became Moncrieff 1997. The run finished, the numbers were labelled as
the method that had been asked for, and nothing downstream could tell. The
first two change every flux in the run, not just the spectral correction.

It is not hypothetical. Every one of the 39 regression fixtures named at least
one path that does not exist, so the suite had been running covariance
maximisation in place of PWB and Moncrieff in place of Fratini - the two
methods most of the fixtures were written to exercise - and passing, because
the sweep checks that each row matches its header and a degraded run still
does.

The engine already stopped for three other missing inputs (the project file,
the alternative metadata file, the raw data directory). These checks hold the
remaining five to the same rule, and hold the message to being useful: an abort
that says "not found" and stops there just moves the dead end.
"""

import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
GUI = ROOT.parent / "eddyflow-gui"

ROUTINE = "src/src_common/abort_on_missing_path.f90"

#: file -> (the substitution that must not come back, the setting it guards)
SITES = {
    "src/src_rp/read_planar_fit_file.f90": ("Meth%rot = 'double_rotation'", "pf_file"),
    "src/src_rp/read_timelag_opt_file.f90": ("Meth%tlag = 'maxcov'", "to_file"),
    "src/src_fcc/read_ini_fcc.f90": (None, "sa_bin_spectra"),
    "src/src_fcc/read_spectral_assessment_file.f90": (None, "sa_file"),
}

#: The interface labels the remedies quote. A reworded option on the GUI side
#: would otherwise leave the engine giving directions that fit nothing on
#: screen, which is worse than saying nothing.
QUOTED_LABELS = (
    "Planar fit file not available",
    "Time lag file not available",
    "Spectral assessment file not available",
    "Full w/Ts cospectra files not available",
)


def read(rel, root=ROOT):
    return (root / rel).read_text(encoding="utf-8", errors="replace")


def code(rel):
    return "\n".join(ln for ln in read(rel).splitlines()
                     if not ln.strip().startswith("!"))


class TheRunStopsTests(unittest.TestCase):
    def test_the_routine_aborts_rather_than_returning(self):
        body = code(ROUTINE)
        self.assertIn("stop 1", body,
                      "a routine named Abort that returns is a warning")

    def test_every_site_calls_it(self):
        for rel, (_, setting) in SITES.items():
            self.assertIn("AbortOnMissingPath", code(rel),
                          f"{rel} no longer stops on a missing {setting}")

    def test_no_site_substitutes_a_method_instead(self):
        """This is the defect itself, and the line most likely to creep back.

        Scoped to the branch the open failed on. The same substitution appears
        legitimately where the file opened and then would not parse: that is a
        malformed file, a different diagnosis, and it keeps its old handling.
        """
        for rel, (substitution, setting) in SITES.items():
            if substitution is None:
                continue
            body = code(rel)
            #> The last else of the open block - "the file is not there".
            arm = body[body.rindex("\n    else\n"):]
            arm = arm[:arm.index("\n    end if")]
            self.assertIn("AbortOnMissingPath", arm,
                          f"{rel} does not stop when {setting} is missing")
            self.assertNotIn(
                substitution, arm,
                f"{rel} silently switches method again on a missing {setting}; "
                f"the whole point is that it stops")

    def test_the_fcc_does_not_demote_the_method_on_a_missing_directory(self):
        body = code("src/src_fcc/read_ini_fcc.f90")
        for guard in ("inquire(file = Dir%full", "inquire(file = Dir%binned"):
            arm = body[body.index(guard):]
            arm = arm[:arm.index("end if")]
            self.assertNotIn(
                "hf_meth = 'moncrieff_97'", arm,
                "a missing co-spectra directory demotes the method again "
                "instead of stopping")


class TheMessageHelpsTests(unittest.TestCase):
    def test_it_names_the_setting_and_the_path(self):
        body = code(ROUTINE)
        for piece in ("setting", "path"):
            self.assertIn(piece, body)
        self.assertIn("names a path that does not exist", read(ROUTINE))

    def test_every_call_passes_a_remedy(self):
        """Three arguments, always: stopping without saying what to change is
        half a fix."""
        for rel in SITES:
            for m in re.finditer(r"call AbortOnMissingPath\((.*?)\)\n",
                                 code(rel), re.S):
                self.assertGreaterEqual(
                    m.group(1).count(","), 2,
                    f"{rel} calls AbortOnMissingPath without a remedy")

    def test_the_remedy_survives_being_printed(self):
        """It was first passed as an array of lines, which compiled cleanly and
        printed nothing - an assumed-shape dummy needs an explicit interface,
        and these callers have none, so size() gave zero and the advice
        vanished. One string, wrapped in the routine."""
        body = code(ROUTINE)
        self.assertNotIn("remedy(:)", body,
                         "an assumed-shape remedy reaches these external "
                         "callers without an interface and prints nothing")
        self.assertIn("character(*), intent(in) :: remedy", body)

    @unittest.skipUnless(GUI.exists(), "eddyflow-gui not beside the engine")
    def test_each_quoted_label_still_exists_in_the_interface(self):
        """The remedies name options as the interface labels them. If a label
        is reworded the advice points at nothing, which is worse than generic
        advice - so the two are pinned together."""
        engine = "\n".join(read(rel) for rel in SITES)
        gui = "\n".join(
            (GUI / p).read_text(encoding="utf-8", errors="replace")
            for p in ("src/planarfitsettingsdialog.cpp",
                      "src/timelagsettingsdialog.cpp",
                      "src/advspectraloptions.cpp"))
        for label in QUOTED_LABELS:
            self.assertIn(label, engine,
                          f"no remedy quotes {label!r} any more")
            self.assertIn(label, gui,
                          f"the engine tells the user to choose {label!r}, "
                          f"which no longer appears in the interface")



class SharedFitArraysTests(unittest.TestCase):
    """xFit, yFit, zFit and ddum are one global set, shared by two routines
    that size them differently under the same "only if not already allocated"
    test. FitRHtoCutoff claimed them at a fixed 10 and then returned early on a
    project with no primary water without freeing them, so the next FitTFModels
    found them allocated, skipped its own sizing, and indexed past the end -
    a crash, not a wrong number.

    It only ever bit projects with no water, and those never reached it while a
    missing binned-spectra directory was quietly demoting every run to
    Moncrieff. Making that missing path fatal is what exposed it.
    """

    RH = "src/src_fcc/fit_rh_to_cutoff.f90"
    TF = "src/src_fcc/fit_tf_models.f90"

    def test_the_fixed_size_claim_comes_after_the_early_return(self):
        body = code(self.RH)
        guard = body.index("if (primary < firstGas) return")
        alloc = body.index("allocate(xFit(10))")
        self.assertLess(
            guard, alloc,
            "FitRHtoCutoff claims the shared fit arrays before the guard that "
            "returns without freeing them")

    def test_the_larger_consumer_resizes_rather_than_trusting(self):
        body = code(self.TF)
        self.assertIn(
            "if (size(xFit) < maxval(nlong)) deallocate(xFit)", body,
            "FitTFModels trusts an existing allocation again; whoever needs "
            "more has to ask for more, or a ten-element array left behind by "
            "the other routine caps it silently")


if __name__ == "__main__":
    unittest.main()
