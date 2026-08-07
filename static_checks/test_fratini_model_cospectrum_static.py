"""Fratini 2012 corrects every gas against the *measured w/T* cospectrum.

That substitution is the method. `BPCF_Fratini12` reads one column from the
full-cospectra file - `cov(w_ts)` - divides it by the sonic transfer function
to recover a "more theoretical" model cospectrum, and then integrates that one
curve against each gas's band-pass transfer function. The gas columns are
deliberately *not* imported: `wanted(firstGas:lastGas) = .false.`, with the
loc_var_present line left commented above it.

Converting the four hard-coded calls into a firstGas..lastGas loop moved the
first argument from `fullCospectra%of(w_ts)` to `fullCospectra%of(gas)`, which
looks like the same mechanical substitution as the nine sibling call sites and
is not. Everywhere else the cospectrum index and the corrected variable are the
same slot, so `w_co2 -> gas` is faithful. Here they are different by design.

The failure is silent in all three ways that matter:

  - `wanted` excludes the gas slots, so `ImportFullCospectra` leaves them at
    `error` rather than absent - the array is the right shape and full of
    sentinels.
  - `SpectralCorrectionFactors` sees `err_cnt == nfreq`, returns early and
    sets `BPCF%of(var) = error`. No exception, no log line.
  - the plausibility band below then reads `-9999 <= 0.8` as "implausible
    direct factor" and calls `CorrectionFactorsIbrom07`, which overwrites it
    with a perfectly reasonable number.

So every gas in every period fell back to Ibrom 2007 and the output looked
well-formed. It was caught by diffing a one-day CH-LAE project against v7.2.5:
mid-day `h2o_scf` had dropped from ~1.05-1.10 to ~1.02-1.03, uniformly, while
every uncorrected flux, covariance, time lag and binned cospectrum was
bit-identical - which is the signature of a correction factor that stopped
being computed rather than one that changed.
"""

from pathlib import Path
import re
import unittest


ROOT = Path(__file__).resolve().parents[1]

FRATINI = "src/src_common/bpcf_fratini_12.f90"
AUX = "src/src_common/bpcf_aux_subs.f90"


def read(path):
    return (ROOT / path).read_text(encoding="utf-8")


def code(path):
    """Source with comment-only lines dropped.

    Both the file and this defect are documented in prose right next to the
    code, and the retired call sites survive as comments, so a naive search
    matches the explanation rather than what compiles.
    """
    return "\n".join(ln for ln in read(path).splitlines()
                     if not ln.lstrip().startswith("!"))


def fratini_body():
    src = code(FRATINI)
    return src.split("subroutine BPCF_Fratini12")[1].split("end subroutine")[0]


class TheModelCospectrumIsTheSonicOne(unittest.TestCase):
    def test_every_gas_is_corrected_against_w_ts(self):
        body = fratini_body()
        self.assertRegex(
            body,
            r"call\s+SpectralCorrectionFactors\(\s*fullCospectra%of\(w_ts\)\s*,"
            r"\s*gas\s*,",
            "BPCF_Fratini12 must integrate the measured w/T cospectrum for "
            "every gas; only the corrected slot varies")

    def test_no_gas_is_corrected_against_its_own_cospectrum(self):
        body = fratini_body()
        self.assertNotRegex(
            body, r"SpectralCorrectionFactors\(\s*fullCospectra%of\(\s*gas\s*\)",
            "indexing the model cospectrum by gas reads a slot `wanted` "
            "excludes from the import - it is all error, and every gas then "
            "silently falls back to Ibrom 2007")

    def test_the_loop_still_reaches_every_gas(self):
        """The w_ts argument is the fix; the N-gas loop is what it was
        introduced for and must survive it."""
        body = fratini_body()
        self.assertRegex(
            body, r"do\s+gas\s*=\s*firstGas\s*,\s*lastGas",
            "the correction must be applied over every configured gas, not "
            "over the four historical slots")


class TheGasCospectraAreNotImported(unittest.TestCase):
    """What makes the wrong index fatal rather than merely redundant. If this
    ever changes, the assertion above stops being the whole story."""

    def test_only_the_sonic_cospectrum_is_requested(self):
        body = fratini_body()
        self.assertRegex(
            body, r"wanted\(w_ts\)\s*=\s*\.true\.",
            "the sonic cospectrum is the one column this method reads")
        self.assertRegex(
            body, r"wanted\(firstGas:lastGas\)\s*=\s*\.false\.",
            "the gas cospectra are deliberately not imported here")


class AnAllErrorCospectrumYieldsNoFactor(unittest.TestCase):
    """The short-circuit that turned a wrong index into a silent one. Pinned
    so the mechanism stays visible next to the check that depends on it."""

    def test_spectral_correction_factors_returns_error(self):
        src = code(AUX)
        body = src.split("subroutine SpectralCorrectionFactors")[1] \
                  .split("end subroutine")[0]
        self.assertRegex(
            body, r"if\s*\(\s*err_cnt\s*==\s*nfreq\s*\)",
            "a cospectrum of nothing but error codes must be recognised")
        self.assertRegex(
            body, r"BPCF%of\(var\)\s*=\s*error",
            "and must yield error rather than a computed factor")

    def test_the_plausibility_band_rejects_that_error(self):
        """`-9999 <= 0.8` is true, so the error propagates as a fallback
        rather than as a gap. This is why nothing in the output looked wrong."""
        body = fratini_body()
        self.assertRegex(
            body, r"min_bpcf_f12\(:\)\s*=\s*0\.8d0",
            "the lower bound is a plausibility floor, and it sits far above "
            "the error sentinel")
        self.assertRegex(
            body,
            r"BPCF%of\(gas\)\s*<=\s*min_bpcf_f12\(gas\)",
            "the band is tested per gas, and an error factor fails it")


if __name__ == "__main__":
    unittest.main()
