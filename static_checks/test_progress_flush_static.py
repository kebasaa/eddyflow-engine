"""Progress reaches the console and the run log when it is written.

gfortran buffers a unit in blocks when it is not a terminal - the interface
reading the engine through a pipe, a pre-pass worker writing to its .out file.
On the Yatir run every worker's .out file stayed at 0 bytes for an hour and the
run log in out_path was not written after start-up, so the only visible
progress was the parent's own slice, arriving in bursts.

The run log was already flushed per LogSay line once it had a name. What was
not: the console, the no-advance writes the period loop uses for its per-period
"#", and the daily progress line DisplayProgress writes - which is what a long
pre-pass actually shows. All of them now go through LogFlush, which flushes
both. Content is unchanged; only when it arrives.

Part of the EddyFlow engine's static checks.
"""

from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]


def code(rel):
    text = (ROOT / rel).read_text(encoding="utf-8", errors="replace")
    return "\n".join(ln for ln in text.splitlines() if not ln.lstrip().startswith("!"))


LOG = code("src/src_common/m_log.f90")
PROGRESS = code("src/src_common/show_daily_advancement.f90")


def body(text, name):
    return text[text.index("subroutine %s(" % name):text.index("end subroutine %s" % name)]


class EveryProgressWriteIsFlushed(unittest.TestCase):

    def test_log_flush_flushes_the_console_and_the_named_log(self):
        b = body(LOG, "LogFlush")
        self.assertIn("flush(output_unit)", b)
        self.assertIn("if (Named) flush(ulog)", b)
        self.assertIn("use iso_fortran_env, only: iostat_end, output_unit", LOG)

    def test_each_log_routine_flushes_even_before_the_log_is_connected(self):
        for name in ("LogSay", "LogSayList", "LogSayNoAdv"):
            b = body(LOG, name)
            self.assertIn("call LogFlush()", b, name)
            self.assertNotIn("if (.not. Connected) return", b,
                             "%s must not return before flushing the console" % name)

    def test_the_daily_progress_line_is_flushed(self):
        self.assertIn("use m_log, only: LogFlush", PROGRESS)
        b = body(PROGRESS, "DisplayProgress")
        self.assertLess(b.index("end select"), b.index("call LogFlush()"))

    def test_log_flush_is_public(self):
        self.assertIn("LogIsOpen, LogFlush", LOG)


if __name__ == "__main__":
    unittest.main()
