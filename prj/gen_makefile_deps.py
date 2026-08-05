#!/usr/bin/env python3
"""Regenerate the Makefile's per-object dependency block.

Companion to gen_project_tags.py and gen_metadata_tags.py, same contract: a
list the build depends on, kept by a script rather than by hand.

The list used to be hand-written, and it drifted the way hand-written lists
do. An object with no rule at all is the worst case, because make treats a
target with no prerequisites as up to date the moment the file exists: the
source is then edited, the object is not rebuilt, and the previous compile is
silently linked instead. The Makefile still carries a comment recording that
happening once to gas_slot_resolution.o. It happened again to
parse_month_grouping.o, which surfaced only as an undefined reference at link
time - and would not have surfaced at all had the edit not added a symbol.

Two kinds of dependency are emitted:

  the source itself   so an edit rebuilds the object. This is what a missing
                      entry costs.

  every `use`d module  named as the .o that produces it, because gfortran
                      writes the .mod and the .o together. Transitive
                      dependencies are left to make: m_rp_global_var.o already
                      depends on m_common_global_var.o, so a file that uses
                      only the former still rebuilds when the latter changes.

  every `include`d file  which the hand-written list had none of. Three files
                      - interfaces.inc, interfaces_1.inc and
                      version_and_date.inc - are included by more than fifty
                      sources between them, and editing one rebuilt nothing.

Paths are basenames throughout: the Makefile's VPATH covers all four source
directories and the object directory, which is how the existing rules already
resolve `m_typedef.f90` and `m_common_global_var.o` without a path.

Usage:  python gen_makefile_deps.py [--check]

--check exits non-zero if the Makefile is not what the generator would write.
"""

import argparse
import glob
import os
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MAKEFILE = ROOT / "prj" / "Makefile"

#: Where the compiler rules look for sources. Order is irrelevant - basenames
#: are unique across all four, which the generator asserts.
SOURCE_GLOBS = (
    "src/src_rp/*.f90",
    "src/src_rp/fft4/*.F",
    "src/src_fcc/*.f90",
    "src/src_common/*.f90",
)

BEGIN = "#> BEGIN GENERATED dependencies - edit gen_makefile_deps.py, not this block"
END = "#> END GENERATED dependencies"

#: `module foo`, but not `module procedure foo` and not `end module foo`.
MODULE_DECL = re.compile(r"(?im)^[ \t]*module[ \t]+([a-z][\w]*)[ \t]*$")
#: `use foo`, `use foo, only: bar`, `use, intrinsic :: iso_c_binding`.
MODULE_USE = re.compile(r"(?im)^[ \t]*use[ \t]+([a-z][\w]*)")
INCLUDE = re.compile(r"""(?i)^[ \t]*include[ \t]+['"]([^'"]+)['"]""", re.M)


def sources():
    out = []
    for pattern in SOURCE_GLOBS:
        out += sorted(glob.glob(str(ROOT / pattern)))
    seen = {}
    for path in out:
        stem = os.path.splitext(os.path.basename(path))[0]
        if stem in seen:
            raise SystemExit(
                "two sources share the basename %r (%s and %s). They compile "
                "to one object, so the build is already ambiguous."
                % (stem, seen[stem], path))
        seen[stem] = path
    return out


def read(path):
    return Path(path).read_text(encoding="utf-8", errors="replace")


def module_owners(paths):
    """module name (lowercased) -> the object that produces its .mod."""
    owners = {}
    for path in paths:
        stem = os.path.splitext(os.path.basename(path))[0]
        for m in MODULE_DECL.finditer(read(path)):
            owners[m.group(1).lower()] = stem
    return owners


def included(path, seen=None):
    """Every file `include`d by `path`, following nesting.

    interfaces.inc includes interfaces_1.inc, so a source naming only the
    first depends on both. Stopping at the first level would rebuild
    twenty-four objects on an edit to interfaces_1.inc and miss the
    twenty-five that reach it through interfaces.inc.
    """
    seen = set() if seen is None else seen
    for name in INCLUDE.findall(read(path)):
        base = os.path.basename(name)
        if base in seen:
            continue
        seen.add(base)
        #> Includes are written relative to the including file. Resolve
        #> against it first, then fall back to the common directory, which is
        #> where all three live.
        nested = Path(path).parent / name
        if not nested.is_file():
            nested = ROOT / "src" / "src_common" / base
        if nested.is_file():
            included(nested, seen)
    return seen


def dependencies(path, owners):
    stem = os.path.splitext(os.path.basename(path))[0]
    text = read(path)

    mods = set()
    for m in MODULE_USE.finditer(text):
        owner = owners.get(m.group(1).lower())
        #> An unknown module is intrinsic (iso_fortran_env and friends) and
        #> has no object to wait for.
        if owner and owner != stem:
            mods.add(owner + ".o")

    return ([os.path.basename(path)]
            + sorted(included(path)) + sorted(mods))


def block(paths, owners):
    lines = []
    for path in paths:
        stem = os.path.splitext(os.path.basename(path))[0]
        deps = dependencies(path, owners)
        lines.append("%s.o: \\" % stem)
        for i, dep in enumerate(deps):
            tail = "" if i == len(deps) - 1 else " \\"
            lines.append("\t%s%s" % (dep, tail))
    return lines


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()

    paths = sorted(sources(), key=lambda p: os.path.basename(p).lower())
    owners = module_owners(paths)

    raw = MAKEFILE.read_text(encoding="utf-8", errors="surrogateescape")
    nl = "\r\n" if "\r\n" in raw else "\n"
    if BEGIN not in raw or END not in raw:
        raise SystemExit(
            "marker pair missing from the Makefile. Add\n  %s\n  %s\n"
            "around the dependency list at the bottom." % (BEGIN, END))
    head, rest = raw.split(BEGIN, 1)
    _, tail = rest.split(END, 1)

    body = nl + nl.join(block(paths, owners)) + nl
    new = head + BEGIN + body + END + tail

    print("%d objects, %d modules, dependency block %d lines"
          % (len(paths), len(owners), body.count(nl)))

    if new == raw:
        print("Makefile dependencies up to date")
        return 0
    if args.check:
        print("\nMakefile dependency block is stale - "
              "re-run gen_makefile_deps.py")
        return 1
    MAKEFILE.write_text(new, encoding="utf-8", errors="surrogateescape")
    print("rewrote the dependency block")
    return 0


if __name__ == "__main__":
    sys.exit(main())
