"""Static checks for the ex-file (FLUXNET CSV) positional layout.

The ex file is navigated positionally: after each list-directed read,
ReadExRecord discards the consumed fields by skipping exactly that many
commas. Those counts used to be bare literals (263, 23, 12, 38, ...). A wrong
count does not raise an error - it shifts every following field by one and the
record is silently misread, so these are worth pinning down.

The counts are now named parameters expressed as sums of the groups they
cover. These checks assert (a) no bare literal crept back in, and (b) the
expressions still evaluate to the historical values for the 4-gas layout, so
the refactor stays provably equivalent for existing files.
"""

import re
import unittest
from pathlib import Path

SRC = Path(__file__).resolve().parent.parent / "src"
READER = SRC / "src_common" / "read_ex_record.f90"

#: Field counts of the shipped 4-gas layout. These are the literals that were
#: in the source before the parameters were introduced; they must not drift.
#:
#: A count stays here after it is widened to N gases. That is the whole point:
#: the layout must still describe the shipped file when the project configures
#: four gases, which is what makes the byte-identical regression a proof that a
#: conversion was faithful rather than merely plausible.
HISTORICAL = {
    "nMainFields": 263,
    "nNrexFields": 23,
    "nVmFields": 12,
    "nLgdFields": 38,
    "nSsItcFields": 9,
    "nSsTestFields": 9,
    "nLicorFields": 37,
    "nAgcFields": 3,
    "nWboostFields": 3,
    "nRotFields": 5,
    "nTlagMethFields": 4,
    "nMetaFixedFields": 22,
    "nMetaGasFields": 11,
    "nCecFields": 11,
}

#: Counts that have been converted from a compile-time `parameter` to a runtime
#: quantity, because the group they measure now carries one set per configured
#: gas. Each must be assigned from GAS_TERM, and must still reduce to its
#: HISTORICAL value when that term is 4.
#:
#: Add a name here in the same edit that demotes it in the reader. Leaving it
#: out is caught (it is no longer a parameter, so the parameter check fails);
#: adding it without demoting it is caught too (no runtime assignment exists).
RUNTIME_COUNTS = {
    "nNrexFields",
    "nVmFields",
    "nLgdFields",
    "nSsTestFields",
    "nSsItcFields",
    "nLicorFields",
}

#: The expression that stands for "how many gases this project configures".
#: Spelled exactly once here so a conversion that invents its own clamp - or
#: forgets to clamp at all - shows up as a failure.
GAS_TERM = "min(EddyFlowProj%gas_num, MaxNumGases)"

#: Counts whose value must track the gas count rather than being re-frozen at
#: its four-gas value. A `parameter` satisfies this symbolically, through
#: nExGas/nExVar/nExScal; a runtime count satisfies it through GAS_TERM.
GAS_DEPENDENT = (
    "nMainFields",
    "nNrexFields",
    "nVmFields",
    "nLgdFields",
    "nSsItcFields",
    "nSsTestFields",
    "nLicorFields",
)

#: Slot numbers from m_typedef.f90. Read rather than hard-coded so that a
#: change to the gas enumeration shows up here as a failure, not a surprise.
SLOT_NAMES = ("u", "ts", "co2", "gas4")


def _read_slots():
    text = (SRC / "src_common" / "m_typedef.f90").read_text(
        encoding="utf-8", errors="replace"
    )
    slots = {}
    for name in SLOT_NAMES:
        m = re.search(
            r"^\s*integer,\s*parameter\s*::\s*%s\s*=\s*(\d+)" % name,
            text,
            re.MULTILINE | re.IGNORECASE,
        )
        assert m, "could not find slot constant '%s' in m_typedef.f90" % name
        slots[name] = int(m.group(1))
    return slots


def _reader_text():
    return READER.read_text(encoding="utf-8", errors="replace")


def _join_continuations(text):
    """Reduce the parameter declarations to one line each.

    Only the declaration block is considered, so stripping '!' comments cannot
    disturb string literals elsewhere in the file. In Fortran a continuation
    '&' may be followed by a trailing comment, so comments go first.
    """
    start = text.find("integer, parameter :: nExGas")
    end = text.find("include 'interfaces_1.inc'", start if start >= 0 else 0)
    block = text[start:end] if start >= 0 and end > start else text
    block = re.sub(r"!.*", "", block)
    return re.sub(r"&\s*\n\s*&?", " ", block)


def _parse_parameters(text):
    """Evaluate the integer parameters declared in ReadExRecord.

    Handles Fortran line continuations and trailing '!<' comments, then
    evaluates each right-hand side in order so later parameters can refer to
    earlier ones.
    """
    joined = _join_continuations(text)
    env = dict(_read_slots())
    for name, expr in re.findall(
        r"^\s*integer,\s*parameter\s*::\s*(\w+)\s*=\s*([^\n!]+?)\s*$",
        joined,
        re.MULTILINE,
    ):
        expr = expr.replace("/ 2", "// 2")  # Fortran integer division
        try:
            env[name] = eval(expr, {"__builtins__": {}}, dict(env))  # noqa: S307
        except Exception:  # pragma: no cover - surfaced by the assertions below
            continue
    return env


#: GAS_TERM, tolerant of the whitespace a Fortran author might use.
GAS_TERM_RE = re.compile(
    r"min\s*\(\s*EddyFlowProj%gas_num\s*,\s*MaxNumGases\s*\)", re.IGNORECASE
)


def _runtime_assignment(text, name):
    """Right-hand side of `name = ...` in the body of ReadExRecord.

    A demoted count is declared bare (`integer :: nNrexFields`) and assigned
    where it is used, so it cannot be recovered from the declaration block the
    parameter parser reads.
    """
    m = re.search(
        r"^[ \t]*%s\s*=\s*((?:[^\n!&]|&[ \t]*\n[ \t]*&?)+)" % re.escape(name),
        text,
        re.MULTILINE,
    )
    if not m:
        return None
    return re.sub(r"&[ \t]*\n[ \t]*&?", " ", m.group(1)).strip()


def _gas_count_names(text):
    """Locals that hold the configured gas count.

    The clamp is worth computing once when several groups are sized from it,
    so a count may be written in terms of such a local rather than repeating
    the expression. Resolving one level of indirection lets that read well
    without weakening the check: the local must itself come from GAS_TERM.
    """
    return {
        m.group(1)
        for m in re.finditer(
            r"^[ \t]*(\w+)\s*=\s*min\s*\(\s*EddyFlowProj%gas_num\s*,\s*"
            r"MaxNumGases\s*\)",
            text,
            re.MULTILINE | re.IGNORECASE,
        )
    }


def _substitute_gas_count(rhs, text, value="4"):
    """Replace the gas count - however it is spelled - with `value`.

    Returns (expression, substituted). `substituted` is False when the
    expression does not depend on the gas count at all, which is the case a
    re-frozen literal produces.
    """
    out, n = GAS_TERM_RE.subn(value, rhs)
    for name in sorted(_gas_count_names(text), key=len, reverse=True):
        out, k = re.subn(r"\b%s\b" % re.escape(name), value, out)
        n += k
    return out, n > 0


class ExRecordLayout(unittest.TestCase):
    """The ex file is navigated by counting commas, so the counts are the format.

    These were pytest-style module functions, which the project's
    `unittest discover` runner collects as zero - so they had been silently
    absent for two landings of the very refactor they exist to guard, and three
    of them had gone stale in the meantime.
    """

    def test_reader_exists(self):
        self.assertTrue(
            READER.is_file(), "read_ex_record.f90 not found at %s" % READER
        )

    def test_no_bare_field_counts(self):
        """Every comma-skip must use a named count, not a literal.

        The sole exception is the per-custom-variable loop, which advances one
        field at a time; that 1 is a step, not a layout count.
        """
        text = _reader_text()
        offenders = [
            (n, lit)
            for n, line in enumerate(text.splitlines(), 1)
            for lit in re.findall(
                r"strCharIndex\(dataline,\s*','\s*,\s*(\d+)\s*\)", line
            )
            if lit != "1"
        ]
        self.assertFalse(
            offenders,
            "bare field counts in read_ex_record.f90 (use a named count): %s"
            % ", ".join("line %d: %s" % o for o in offenders),
        )

    def test_parameters_match_historical_layout(self):
        """Counts still declared `parameter` must describe the shipped file."""
        env = _parse_parameters(_reader_text())
        expected = {k: v for k, v in HISTORICAL.items() if k not in RUNTIME_COUNTS}
        missing = sorted(k for k in expected if k not in env)
        self.assertFalse(
            missing,
            "field-count parameters not found or unparsable: %s. If one was "
            "deliberately converted to a runtime count, add it to "
            "RUNTIME_COUNTS." % missing,
        )
        wrong = {k: (env[k], want) for k, want in expected.items() if env[k] != want}
        self.assertFalse(
            wrong,
            "ex-file field counts changed for the 4-gas layout; existing files "
            "would be misread. Got (value, expected): %s" % wrong,
        )

    def test_runtime_counts_reduce_to_the_historical_layout(self):
        """A widened count must still emit today's columns for four gases.

        This is the replacement for the guarantee `parameter` gave for free.
        Demoting a count to a runtime expression removes it from the check
        above; without this one, a group could be widened to something that
        never matched the shipped file and nothing would say so.
        """
        text = _reader_text()
        for name in sorted(RUNTIME_COUNTS):
            with self.subTest(count=name):
                self.assertIsNotNone(
                    re.search(
                        r"^[ \t]*integer[ \t]*::[ \t]*%s\b" % name, text, re.MULTILINE
                    ),
                    "%s is in RUNTIME_COUNTS but is not declared as a bare "
                    "`integer ::` in ReadExRecord" % name,
                )
                rhs = _runtime_assignment(text, name)
                self.assertIsNotNone(
                    rhs, "%s is never assigned; the comma skip would use "
                    "whatever was left on the stack" % name
                )
                four_gas, substituted = _substitute_gas_count(rhs, text)
                self.assertTrue(
                    substituted,
                    "%s must be sized from '%s' - directly, or through a local "
                    "assigned from it. Got '%s', which does not depend on the "
                    "configured gas count at all."
                    % (name, GAS_TERM, rhs),
                )
                four_gas = four_gas.replace("/ 2", "// 2")
                value = eval(  # noqa: S307 - fixed expression from our own source
                    four_gas, {"__builtins__": {}}, dict(_read_slots())
                )
                self.assertEqual(
                    value,
                    HISTORICAL[name],
                    "%s no longer reduces to the shipped 4-gas width (%s gives "
                    "%d, expected %d), so the byte-identical regression can no "
                    "longer prove the conversion faithful"
                    % (name, rhs, value, HISTORICAL[name]),
                )

    def test_counts_are_gas_dependent(self):
        """A gas-dependent count must actually reference the gas count.

        One silently re-hard-coded to its 4-gas value would pass both checks
        above while breaking as soon as the capacity grows.
        """
        text = _reader_text()
        joined = _join_continuations(text)
        for name in GAS_DEPENDENT:
            with self.subTest(count=name):
                if name in RUNTIME_COUNTS:
                    rhs = _runtime_assignment(text, name)
                    self.assertIsNotNone(rhs, "runtime count %s not assigned" % name)
                    _, substituted = _substitute_gas_count(rhs, text)
                    self.assertTrue(
                        substituted,
                        "%s no longer tracks the configured gas count" % name,
                    )
                    continue
                m = re.search(
                    r"integer,\s*parameter\s*::\s*%s\s*=\s*([^\n!]+)" % name, joined
                )
                self.assertIsNotNone(m, "parameter %s not found" % name)
                self.assertRegex(
                    m.group(1),
                    r"nExGas|nExVar|nExScal",
                    "%s no longer depends on the gas slot count; it will not "
                    "track a change in gas capacity" % name,
                )


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
