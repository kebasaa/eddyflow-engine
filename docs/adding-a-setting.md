# Adding a setting

House rules for adding a source file or a project setting to the engine and its
interface. Each one is here because breaking it produces a *silent* wrong
result rather than an error — a stale binary, a key the engine never looks for,
a setting that re-points every later one. None of them is guessable from the
surrounding code.

Collected while porting the EddyUH capabilities; the comparison that drove
that work, and what became of each item, is in
[`eddyuh-comparison.md`](eddyuh-comparison.md).

## Every new source file

```bash
python prj/gen_makefile_deps.py
```

Then run `static_checks/test_makefile_deps_static.py` **before trusting any
build**. Without the regenerated rule, make treats the object as up to date the
moment the file exists and silently links the previous compile. This has bitten
`gas_slot_resolution.o` and `parse_month_grouping.o`.

The related trap, which no generator can catch: after editing a source file,
check that the binary is actually newer than it before drawing a conclusion
from a run. Make will decline to rebuild an object whose timestamp already
looks current, and the run then tests the previous compile.

## Every new ini key

1. `prj/gen_project_tags.py` → `FIXED_TAGS`, at a free index **below the record
   origin**. Above it silently re-indexes ~1600 per-gas slots;
   `test_ini_tag_collisions_static.py` asserts `rpGasOriginN == 425`.
2. `python prj/gen_project_tags.py` — the tag tables in
   `m_common_global_var.f90` and `m_rp_global_var.f90` are **generated**, not
   hand-edited, and the size parameters grow with them.
3. `m_typedef.f90` — the setup-type field.
4. Reader — a **literal default**, then an `if (…TagFound(n))` guarded
   override. Absent ≠ empty: an unguarded numeric read returns whatever the
   previous parse left in the `save`d array.
5. Consumers.
6. `m_eddypro_import.f90` — `defaultSect` / `defaultKey` / `defaultValue`, and
   bump `nProjectDefaults`. The value **must equal** the reader's default;
   `test_eddypro_import_static.py` checks it against the reader that owns it.
   Leave a key out only when *absent* is itself a decision no value reproduces
   — the file lists those and says why for each.
7. GUI: `ecinidefs.h`, `ecprojectstate.h`, `ecproject.h`, and **four** places in
   `ecproject.cpp` — `saveEcProject`, `loadEcProject`, `fuzzyCompare`,
   **`newEcProject`**. Copy `cec_*` as the model, **not** PWB, which is missing
   from `newEcProject` entirely.
8. Static check on both sides, and a CHANGELOG entry on both sides.

### Which INI group

Correctness, not cosmetics. FCC never sweeps `RawProcess_*`
(`ecinidefs.h`; `test_random_uncertainty_static.py` pins the `ru_*` case).

- A key **FCC reads** must be `[Project]`. Only a handful of numeric slots
  remain there — the blanks left by the keys the 5.0.0 record format retired.
- A key **only RP reads** goes in a `RawProcess_*` group, where there is room.

A setting applied at correction time — `hf_meth`, `cosp_model` — is read by
both applications, because RP runs the analytic corrections itself.

### Key naming

`test_ini_tag_collisions_static.py` forbids one label being a substring of
another in the same scope. `detlim_window` and `detlim_window_s` cannot
coexist; neither can `tlag_borrow` and `tlag_borrow_snr`. Rename in
`m_rp_global_var.f90` first — the generator refuses to overwrite a live tag.

### Per-gas keys are free

`gas_<i>_*` keys are append-only in the generator's `RP_GAS_NUMERIC` /
`GAS_NUMERIC` lists and consume none of the flat budget. **Append only, never
insert**, or every later setting silently re-points. The same holds for the
metadata generator's `ANTags` / `ACTags`, which are keyed by index: inserting
by hand produces duplicate indices and one of the pair is dropped without a
word.

### No `ini_version` bump

Additive keys need none. An older engine does not find the key and keeps its
default. Bump only for a change that alters the meaning of something an older
engine already reads.

## Fortran

**External subroutines have no interface block.** An argument-list mismatch
compiles clean and corrupts memory at run time — a logical written into a real
array, and a denormal in the output as the only symptom. When a signature
changes, a static check comparing the call sites' argument lists against it is
the only thing that catches it.

**`intent(out)` assigned only inside a branch returns whatever was on the
stack.** `CovMax` returned an undefined row lag this way until the
baseline-subtraction work went through it. Initialise before the branch.

**An error-coded value can win a comparison on magnitude.** `abs(-9999)` beats
every real covariance. Skip error codes explicitly rather than assuming they
lose.

## GUI

**`clicked`, not `toggled`,** for any gating checkbox. `refresh()` blocks the
*project's* signals, not the *widgets'*, so a `toggled` connection fires during
project load and rewrites a file the user only opened.

**A combo's row should be the ini value** where the engine's decode is a
contiguous run from zero. `setItemData` mappings and `index + 1` arithmetic
duplicate the engine's table, and the two drift.

**Sanitise on load.** A value from a later version lands on no row and
`setCurrentIndex` silently shows the first one — the interface then displays
something the file does not say. Clamp to the default, as the engine does.

**A control greyed by a setting on another page needs the other page named in
its tooltip.** Otherwise the user sees a disabled control with no way to find
out why.

## Verifying

Build from `prj/` in PowerShell, with `C:\Users\jonmuell\mingw64\bin` on
`PATH`:

```bash
mingw32-make rp; mingw32-make fcc
```

Then, in the order of what each actually catches:

- **`test_makefile_deps_static.py`** — first, after adding any file. A missing
  rule means the binary just tested is the previous one.
- **`test_ini_tag_collisions_static.py`** — runs `gen_project_tags.py --check`
  and asserts the record origin. Fails if a key landed above it.
- **`test_ex_record_layout_static.py`** — fails until `HISTORICAL` moves with
  the record.
- **Byte-identical regression** — `tests/regression/`, with every new feature
  defaulted off. The single most important check, and the only one that catches
  a column-layout fault: the static checks cannot, because a header and a row
  can be wrong in the same way and still agree with each other. To prove a
  default unchanged, build a reference binary with the new call sites reverted
  and diff `run.sh ref` against `run.sh chk`.
- **Round-trip** — open the project in the interface, save, confirm no key
  changed value. This is how `newEcProject` omissions surface.
- **Per feature, enabled** — confirm the new columns appear, are finite, and
  are ordered as the header says.

`tests/regression/README.md` describes the fixtures and what each is for.
Three traps it does not cover, each of which costs an hour: never edit
`sweep.sh` while it is running (bash reads a script incrementally, so inserting
a line shifts the byte offsets under it); never run two sweeps at once, since
they share `out_chk`; and a killed run leaves `eddyflow_rp.exe` alive holding
files open, after which every fixture fails with `Device or resource busy` or
`RP FAILED` — which looks exactly like a code regression.
