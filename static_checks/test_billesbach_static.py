"""Static checks for the Billesbach (2011) random-shuffle noise floor.

This routine sat in the tree for years marked ``\\todo Under development`` and
was never called. Reviving it meant repairing five separate faults, and each
of them is the kind that produces numbers rather than an error, so each is
pinned here.

**It was never reached.** ``RIN_Billesbach_11`` had no arm in the dispatch, so
selecting it was impossible. If the arm goes, the method silently becomes
unavailable again rather than failing.

**It wrote nowhere useful.** It put the shuffled covariance into
``Stats%Cov`` - the live covariance matrix - and nothing read the result. Two
faults in one line: no output, and the covariance matrix that ``Fluxes0_rp``
reads fifteen times a few hundred lines later was being overwritten with
noise. Every working method writes ``Essentials%rand_uncer(var)`` instead.

**Its indexing was transposed.** ``Set`` is declared ``Set(N, M)`` - records
down, variables across - and every other routine in the file indexes it
``Set(:, var)``. This one used ``Set(w, 1:N)``, which reads row ``w`` across
``N`` columns of an array that has about seventy, running off the end.

**Its shuffle was biased.** ``RandomBetween(min, max)`` computed
``int((max - min) * x + min)``, which never returns ``max``. Its only caller
is a Fisher-Yates loop where ``j`` must be able to equal ``i``; unable to, the
element at ``i`` always moves, and the loop draws from the permutations that
leave nothing in place rather than from all of them.

**Its generator was unseeded.** ``InitRandomSeed`` is also dead code, so the
sequence came from whatever gfortran seeds itself with - since GCC 7, the
operating system. The same project over the same data would publish a
different noise floor every run.
"""

import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "src"

RU = (SRC / "src_rp" / "random_error_handle.f90").read_text(
    encoding="utf-8", errors="replace")
PRJVARS = (SRC / "src_common" / "write_processing_project_variables.f90").read_text(
    encoding="utf-8", errors="replace")

#: EddyUH's realisation count, and the paper's. Fixed rather than settable.
NTIMES = 20


def body_of(subroutine):
    """The text between `subroutine <name>` and its `end subroutine`."""
    start = RU.index("subroutine %s(" % subroutine)
    end = RU.index("end subroutine %s" % subroutine, start)
    return RU[start:end]


class TheMethodIsReachable(unittest.TestCase):

    def test_the_dispatch_has_an_arm_for_it(self):
        self.assertRegex(
            RU,
            r"case\('billesbach_11'\)",
            "the dispatch no longer offers Billesbach, so selecting it falls "
            "through to `case default` and the method is unavailable",
        )

    def test_the_arm_calls_the_routine(self):
        self.assertIn(
            "call RU_Billesbach_11(", RU,
            "the arm exists but calls nothing - which is the state this "
            "routine spent years in",
        )

    def test_the_project_key_selects_it_as_four(self):
        #> Four, not three: three is already mahrt_98 and the interface maps
        #> its menu onto these numbers, so renumbering would change the
        #> method of every project that states one.
        self.assertRegex(
            PRJVARS,
            r"case\(4\)\s*\n(?:\s*!>[^\n]*\n)*\s*RUsetup%meth = 'billesbach_11'",
            "ru_meth = 4 no longer selects Billesbach",
        )

    def test_mahrt_did_not_move_off_three(self):
        self.assertRegex(
            PRJVARS,
            r"case\(3\)\s*\n\s*RUsetup%meth = 'mahrt_98'",
            "ru_meth = 3 is no longer Mahrt, so an existing project that "
            "states 3 has silently changed method",
        )


class TheOutputGoesWhereTheOthersGo(unittest.TestCase):

    def setUp(self):
        self.body = body_of("RU_Billesbach_11")

    def test_it_fills_rand_uncer(self):
        self.assertIn(
            "Essentials%rand_uncer(var)", self.body,
            "the routine computes a number and stores it nowhere, which is "
            "the fault it shipped with",
        )

    def test_it_never_writes_the_live_covariance_matrix(self):
        self.assertNotIn(
            "Stats%Cov", self.body,
            "writing Stats%Cov here overwrites the covariances Fluxes0_rp "
            "reads later in the same period - every flux would be computed "
            "from shuffled data",
        )

    def test_it_does_not_mutate_the_callers_array(self):
        #> RandomUncertaintyHandle declares Set intent(in) and passes it
        #> straight through. An intent(inout) here would not compile against
        #> that, but it did once, so state the requirement.
        self.assertRegex(
            self.body,
            r"real\(kind = dbl\), intent\(in\) :: Set\(N, M\)",
            "Set is no longer intent(in); the shuffle must work on a copy",
        )

    def test_the_indexing_is_records_down_variables_across(self):
        self.assertNotRegex(
            self.body,
            r"Set\(w, 1:N\)",
            "the transposed indexing is back: Set is Set(N, M), so the w "
            "column is Set(1:N, w)",
        )
        self.assertIn("Set(1:N, w)", self.body)
        self.assertIn("Set(1:N, var)", self.body)


class TheEstimatorIsBillesbachs(unittest.TestCase):

    def setUp(self):
        self.body = body_of("RU_Billesbach_11")

    def test_it_averages_twenty_realisations(self):
        self.assertRegex(
            self.body,
            r"ntimes = %d" % NTIMES,
            "the realisation count is no longer EddyUH's twenty",
        )

    def test_the_statistic_is_the_mean_of_absolute_covariances(self):
        #> Not a standard deviation. Billesbach takes |cov| per realisation
        #> and averages, which for a zero-mean Gaussian lands at sqrt(2/pi)
        #> of the scatter - a different number, deliberately.
        self.assertIn("dabs(cov)", self.body)
        self.assertRegex(self.body, r"acc / dble\(ntimes\)")

    def test_it_shuffles_the_scalar_not_the_wind(self):
        #> Shuffling w once and reusing it for every gas is cheaper and is
        #> what the dead version did, but it makes every gas's estimate a
        #> draw from one realisation.
        self.assertIn("call RandomShuffle(Set(1:N, var)", self.body)

    def test_it_skips_the_two_wind_components_that_have_no_flux_here(self):
        self.assertRegex(self.body, r"if \(var == v \.or\. var == w\) cycle")


class TheShuffleIsUnbiasedAndReproducible(unittest.TestCase):

    def test_random_between_can_return_its_upper_bound(self):
        body = body_of("RandomBetween") if "subroutine RandomBetween" in RU else None
        #> It is a function, not a subroutine, so slice it by hand.
        start = RU.index("integer function RandomBetween")
        end = RU.index("end function RandomBetween", start)
        body = RU[start:end]
        self.assertNotRegex(
            body,
            r"int\(\(max - min\) \* x \+ min\)",
            "the off-by-one is back: this form never returns max, which "
            "biases the Fisher-Yates loop that is its only caller",
        )
        self.assertIn("max - min + 1", body)

    def test_the_generator_is_seeded_once_from_a_fixed_value(self):
        self.assertIn("call SeedShuffleOnce()", RU,
                      "nothing seeds the generator, so the noise floor "
                      "changes from run to run")
        body = body_of("SeedShuffleOnce")
        self.assertIn("random_seed(put = seed)", body)
        self.assertRegex(
            body, r"logical, save :: seeded",
            "the seed is not latched, so it would be reset every period and "
            "successive periods would draw identical permutations",
        )

    def test_nothing_else_in_the_engine_draws_from_that_generator(self):
        #> Seeding here is safe only because the shuffle owns random_number.
        users = []
        for path in sorted((SRC).rglob("*.f90")):
            text = path.read_text(encoding="utf-8", errors="replace")
            if re.search(r"call random_number\(", text):
                users.append(path.name)
        self.assertEqual(
            users, ["random_error_handle.f90"],
            "something else now draws from the same generator, so seeding it "
            "for the shuffle would perturb that too: %s" % users,
        )


if __name__ == "__main__":
    unittest.main()
