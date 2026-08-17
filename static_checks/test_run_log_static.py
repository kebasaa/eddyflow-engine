"""The run log has to stay complete, and the console has to stay untouched.

A run leaves a dozen output files and, until now, no record of what the engine
said while producing them. The log is that record: the same filename pattern as
every other output, with a .log extension.

Fortran cannot duplicate the default output unit, so this is a tee, not a
redirect - a redirect would have silenced the console, and the interface reads
the engine's standard output to drive its progress display. Which means every
message site has to say things twice, and the failure mode is quiet: a message
added later reaches the console, the log misses it, and nobody notices until
they read a log to diagnose something and the line they needed is absent.

So the checks here are structural rather than exemplary. They assert that NO
console write is left without a log twin, anywhere in the engine, rather than
listing the ones that have them.
"""

import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC_DIRS = ("src/src_rp", "src/src_common", "src/src_fcc")

#: m_log is where the tee itself lives; its console writes are the tee.
EXEMPT = {"m_log.f90"}

CONSOLE = re.compile(r"\bwrite\s*\(\s*\*\s*,")
TWIN = re.compile(r"\bwrite\s*\(\s*ulog\s*,")
HELPER = re.compile(r"\bcall\s+LogSay(List|NoAdv)?\s*\(")


def sources():
    for d in SRC_DIRS:
        for p in sorted((ROOT / d).glob("*.f90")):
            if p.name not in EXEMPT:
                yield p


def statements(lines):
    """(index, block) per write statement, continuation lines included."""
    i = 0
    while i < len(lines):
        if CONSOLE.search(lines[i]):
            block = [i]
            while lines[block[-1]].rstrip().endswith("&") and block[-1] + 1 < len(lines):
                block.append(block[-1] + 1)
            yield block
            i = block[-1] + 1
        else:
            i += 1


class RunLogStaticTests(unittest.TestCase):
    def test_every_console_message_reaches_the_log(self):
        """Each write(*) is followed by its write(ulog) twin.

        Both forms count: the statement on its own line, and the one carrying a
        condition - `if (printout) write(*,...)`. The second kind was missed the
        first time round, and 53 messages went to the console alone; the ones
        with advance='no' beside them then ran two progress lines together in
        the log, which is how it was noticed.
        """
        orphans = []
        for path in sources():
            lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
            for block in statements(lines):
                after = block[-1] + 1
                if after >= len(lines) or not TWIN.search(lines[after]):
                    orphans.append(f"{path.name}:{block[0] + 1}  "
                                   f"{lines[block[0]].strip()[:70]}")
        assert not orphans, (
            "console messages with no log twin:\n  " + "\n  ".join(orphans[:20])
        )

    def test_the_twin_carries_the_same_condition(self):
        """`if (cond) write(*,...)` must be twinned by `if (cond) write(ulog,...)`,
        or the log gets a line the console did not, or misses one it did."""
        mismatched = []
        for path in sources():
            lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
            for block in statements(lines):
                head = lines[block[0]]
                after = block[-1] + 1
                if after >= len(lines):
                    continue
                cond = head[:head.index("write")].strip()
                twin = lines[after]
                if not TWIN.search(twin):
                    continue
                twin_cond = twin[:twin.index("write")].strip()
                if cond != twin_cond:
                    mismatched.append(f"{path.name}:{block[0] + 1}  "
                                      f"{cond!r} vs {twin_cond!r}")
        assert not mismatched, (
            "twins guarded differently from their originals:\n  "
            + "\n  ".join(mismatched[:20])
        )

    def test_the_helpers_write_with_the_format_they_replaced(self):
        """LogSay and LogSayList do NOT render alike: list-directed output opens
        a record with a blank, so `write(*,*) 'x'` prints " x" where
        `write(*,'(a)') 'x'` prints "x". Most of the engine's messages - every
        line of the exception handler among them - are the first kind, and
        collapsing them onto the wrong helper would shift every one of them by a
        column."""
        body = (ROOT / "src/src_common/m_log.f90").read_text(encoding="utf-8")
        say = body[body.index("subroutine LogSay(text)"):]
        say = say[:say.index("end subroutine LogSay")]
        assert "write(*, '(a)') text" in say and "write(ulog, '(a)') text" in say

        lst = body[body.index("subroutine LogSayList(text)"):]
        lst = lst[:lst.index("end subroutine LogSayList")]
        assert "write(*, *) text" in lst and "write(ulog, *) text" in lst

    def test_the_log_is_connected_before_the_first_message(self):
        """ulog has to be connected from the program's first statement.

        Connected late, every line up to that point writes to an unconnected
        unit and gfortran invents a fort.163 in whatever directory the run
        started in - which is exactly what happened, and it swallowed the whole
        preamble: the banner, the project file, and any exception raised while
        reading it.
        """
        for main in ("src/src_rp/eddyflow-rp_main.f90",
                     "src/src_fcc/eddyflow-fcc_main.f90"):
            body = (ROOT / main).read_text(encoding="utf-8", errors="replace")
            assert "call LogStart()" in body, main
            first_write = min(
                (body.index(m.group(0)) for m in CONSOLE.finditer(body)),
                default=len(body))
            assert body.index("call LogStart()") < first_write, (
                f"{main} says something before connecting the log"
            )

    def test_the_scratch_file_is_carried_over_not_dropped(self):
        """The preamble lives in the scratch file until the output path is
        known. LogInit copies it across; reading it back with trim() would lose
        a trailing blank, and " Executing EddyFlow " has one."""
        body = (ROOT / "src/src_common/m_log.f90").read_text(encoding="utf-8")
        init = body[body.index("subroutine LogInit"):]
        init = init[:init.index("end subroutine LogInit")]
        assert "status = 'scratch'" in body
        assert "advance = 'no', size = nchars" in init, (
            "the copy must take each record's own length, not trim() it"
        )
        assert "iostat_end" in init, (
            "end-of-record and end-of-file are both negative; told apart by "
            "sign alone, the copy stops at the first blank line"
        )

    def test_the_log_is_named_like_every_other_output(self):
        body = (ROOT / "src/src_common/init_run_log.f90").read_text(encoding="utf-8")
        for piece in ("Dir%main_out", "EddyFlowProj%id", "Log_FilePadding",
                      "Timestamp_FilePadding", "LogExt"):
            assert piece in body, piece


if __name__ == "__main__":
    unittest.main()
