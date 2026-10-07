"""A message raised while a progress line is open keeps its own line.

The engine opens progress lines - "   Absolute limits test.." - and ends them
with " Done." once the step is over. A warning raised in between was written
onto the open line, so its first line read "   Absolute limits test..
Warning(109)> One or more gases ...". The interface matches the progress text
first and stops there, so that first line never reached its warning panel.

Where a step can warn, the warning now either waits until the line is ended,
or the line is opened through m_log, whose message routines end an open line
before they write.
"""

import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def read(rel):
    return (ROOT / rel).read_text(encoding="utf-8", errors="replace").replace("\r\n", "\n")


def body(text, kind, name):
    i = text.index("%s %s(" % (kind, name))
    return text[i:text.index("end %s %s" % (kind, name), i)]


LOG = read("src/src_common/m_log.f90")


class MessagesEndAnOpenLine(unittest.TestCase):

    def test_both_message_routines_end_an_open_line_first(self):
        for name in ("LogSay", "LogSayList"):
            b = body(LOG, "subroutine", name)
            self.assertIn("call EndOpenLine()", b, name)
            self.assertLess(b.index("call EndOpenLine()"), b.index("write(*,"), name)

    def test_an_empty_logsay_only_ends_the_open_line(self):
        b = body(LOG, "subroutine", "LogSay")
        self.assertIn("if (LineOpen .and. len_trim(text) == 0) then", b)

    def test_open_and_end_keep_the_flag(self):
        self.assertIn("LineOpen = .true.", body(LOG, "subroutine", "LogOpenLine"))
        self.assertIn("LineOpen = .false.", body(LOG, "subroutine", "LogEndLine"))
        self.assertIn("LineOpen = .false.", body(LOG, "subroutine", "EndOpenLine"))


class StepsThatCanWarnUseThem(unittest.TestCase):

    def test_absolute_limits_warns_after_its_line_is_done(self):
        b = read("src/src_rp/test_absolute_limits.f90")
        self.assertLess(b.index("write(*,'(a)') ' Done.'"),
                        b.index("call ExceptionHandler(109)"))

    def test_borrowed_time_lags_are_announced_on_their_own_line(self):
        b = read("src/src_rp/timelag_handle.f90")
        self.assertIn("call LogOpenLine('  Compensating time-lags..')", b)
        self.assertIn("call LogEndLine(' Done.')", b)

    def test_main_program_lines_that_can_be_interrupted(self):
        main = read("src/src_rp/eddyflow-rp_main.f90")
        self.assertIn("call LogOpenLine(' Reading alternative metadata file: \"'", main)
        self.assertIn("call LogOpenLine(trim(SectorLine))", main)
        self.assertIn("call LogOpenLine('  Calculating CEC partitioning..')", main)

    def test_messages_that_bypassed_m_log_go_through_it(self):
        self.assertNotIn("write(*,*) ' Error(34)>", read("src/src_common/exception_handler.f90"))
        self.assertIn("call LogSay(trim(pairLine))", read("src/src_common/gas_slot_resolution.f90"))


class LinesAreEnded(unittest.TestCase):

    def test_storage_ends_its_line_when_the_periods_are_not_consecutive(self):
        b = read("src/src_rp/storage.f90")
        i = b.index("Stor%of(firstGas:lastGas)  = error")
        self.assertLess(b.index("call LogSay(' Done.')", i), b.index("return", i))

    def test_each_hygrometer_s_rh_fit_ends_its_own_line(self):
        b = body(read("src/src_fcc/fit_rh_to_cutoff.f90"), "subroutine", "FitRh2Fco")
        i = b.index("call LogSay('Done.')")
        self.assertLess(i, b.index("    end do\n\n    !> The primary's coefficients"))
        self.assertEqual(b.count("call LogSay('Done.')"), 1)


if __name__ == "__main__":
    unittest.main()
