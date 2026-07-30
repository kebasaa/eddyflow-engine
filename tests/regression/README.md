# FLUXNET layout regression harness

Guards the per-gas column layout while it is converted from a fixed four
gases to one set per configured gas. The output format is positional and
FCC parses it by comma counts, so a group converted in the header but not
in the reader shifts every later field and is consumed without error.

## Running

    BASE=base_rec.eddyflow bash run.sh ref     # before the change
    BASE=base_rec.eddyflow bash run.sh chk     # after it
    diff -r out_ref out_chk

`run.sh` now normalises recursively, so `diff -r` on two run directories is
the whole comparison. It used to normalise only the top level, which left 21
of the 25 output files carrying the run timestamp in their names - every one
of them read as added *and* removed, burying whatever had really changed.

Two runs from an unchanged tree must diff clean. If they do not, the
normalisation is incomplete and every verdict from it is worthless; check
that before trusting a difference.

## Fixtures

| file | what it is for |
|---|---|
| `base.eddyflow` | pre-5.0.0 project, no gas records. Must abort with `Fatal error(99)`. |
| `base_rec.eddyflow` | the same selection as records: 4 gases, CH4 configured *without* a column, cell_t/int_p and the 7200 diagnostic as records. **Output must stay byte-identical** to the legacy reference - that is the proof a conversion is faithful. |
| `base_5gas.eddyflow` | adds N2O on column 9. Header and rows must agree; no duplicate column names. |
| `base_dup.eddyflow` | the same species in slots 4 and 5, to check the `_2` disambiguation. |
| `base_8gas.eddyflow` | eight gases across **two** analysers: the MIRO's four plus the LI-7200's CO2 and H2O, and a duplicate N2O. Exercises the capacity target (64 gases, 8 per instrument) and, because the two analysers measure the same species independently, catches slot cross-wiring: `CO2` and `CO2_2` must hold *different* real values, not the same one twice. |
| `base_neg.eddyflow` | `base_5gas` with `al_gas4_min` raised to 400, so COS fails the absolute-limits test. The negative fixture: exactly one gas's columns must move. Diff it against the `base_5gas` run, not against a reference. |
| `base_8gas_ru.eddyflow` | `base_8gas` with random uncertainty on. The only fixture that exercises `random_error_handle.f90` and `integral_turbulence_scale.f90` at all - every other one leaves `RUsetup%meth` at `none`, so those files run their `case('none')` arm and nothing else. Expect a real `RANDUNC_HF` for every gas that has a column, and `-9999` for one that does not. |

> **The `ru_*` keys are in the wrong section, and this fixture works around it.**
> `ru_meth`, `ru_its_meth` and `ru_tlag_max` are declared in `EPPrjNTags`, and
> `ParseIniFile` is called with the section prefix `'Project'`, so those tags
> are only ever matched inside `[Project*]`. The interface writes all three
> into `[RawProcess_RandomUncertainty_Settings]`, where nothing looks for them
> - so `RUsetup%meth` falls to its `case default` of `'none'` for every project
> the interface has saved, and random uncertainty has never actually run.
> `base_8gas_ru` repeats the keys under `[Project]` to get past it. Fixing the
> mismatch properly is a separate change: it turns the feature on for every
> project that asked for it, which moves output that has been `-9999` until now.

## Re-baselinings, and what each one accounted for

`out_ref` is regenerated only when a change is meant to move the output, and
only after every moved cell has been named. So far:

| change | delta against the previous reference |
|---|---|
| main record converted (`3511493`) | FLUXNET: exactly one column removed, `NUM_GAS_EXTRA`; every surviving cell byte-identical. `full_output`: 48 VM97 flag cells widened from 9 to 69 characters, each a pure extension of its old value with `'9'` (test not performed) padding |
| full output per gas | `full_output` **units row only**: 4 whitespace-only cells. Three unit labels and one flux label are `character(32)` and used to be concatenated unpadded, so those fields carried trailing blanks. Every gas now trims. No data cell moved, no column added or removed at four gases |

## The arithmetic cross-check

`nMainFields` in `read_ex_record.f90` should equal the header's main-record
width - the columns from `DOY_START` up to the first `*_NREX`. Measuring that
from a run and comparing against the formula is what localised a six-field
under-count that only appeared past four gases:

| gases | main-record width |
|---|---|
| 4 | 263 |
| 5 | 298 |
| 8 | 409 |

All three are confirmed against real runs. 64 gases gives 4133 by the same
formula, untested for want of a dataset.

## Traps that produced false results before

- **RP only.** FCC recomputes under `fcc_follows`; an RP-only run compares
  files nothing wrote. `run.sh` runs both.
- **Timestamp normalisation must recurse** into the per-period
  subdirectories, or every one reads as an added/removed file.
- **A fixture that is under-specified silently tests the old path.** An
  earlier `base_rec` carried gas records only, so cell temperature, cell
  pressure and the diagnostic were still arriving through the legacy tags
  and every "record path IDENTICAL" verdict rested on that.
- **"Nothing moved" is not a pass.** Perturbing `gas_N_mw` moves nothing
  for a gas reported as a ppb mole fraction - the molecular weight is only
  consulted for the `g_m3`/`mg_m3`/`ug_m3` arms. Use `ug_m3` to exercise it.
- The engine links `libgfortran-5.dll` dynamically; `PATH` needs the MinGW
  bin directory.

## A blocker found with `base_8gas_cell`: per-instrument cell T/P never arrives

`base_8gas_cell` is `base_8gas` with cell records on **both** analysers - the
MIRO's cell_t/int_p on columns 11/12 and the LI-7200's on 23/24. The two sets
carry genuinely different data (`T_CELL` 27.23 against 14.68, `PA_CELL` 0.0700
against 93.14 when each is configured alone), so a working per-instrument
resolution has to show up in the output.

It does not. Against `base_8gas`, adding the LI-7200's cell records moves
**one** column - `NUM_CUSTOM_VARS`, because columns 23/24 stop being custom
variables - and removes the two `CUSTOM_*` columns that carried them. No flux,
no `MV_AIR_CELL`, no cell value changes at all.

Localised with three temporary probes, all reverted:

| point | slot 10 (the LI-7200's CO2) |
|---|---|
| `ApplyCellDiagRecords` | record 3 → src 14 → **slot 73**, record 4 → src 15 → **slot 76**. Applied correctly |
| end of `DefineE2Set` | `cell_ref = 73`, i.e. block 2. Resolved correctly |
| `AirAndCellParameters` | `Stats%Mean(73)` and `Stats%Mean(76)` are **-9999** |

So the records reach the right slots and the gas points at the right block, but
the block carries no statistics by the time the physics reads it - and the
fallback in `AirAndCellParameters` ("no reading in this block") quietly hands
back instrument 1's scalars. Every gas then gets the MIRO's cell conditions.
Block 1 (slots 69-72) works, and those are exactly the slots that also have
legacy names (`tc`/`ti1`/`ti2`/`pi`).

**The window still to bisect** is `DefineE2Set` (eddyflow-rp_main.f90:651) to
`BasicStats` (:728, then :823 before `AirAndCellParameters` at :828), with
`FilterDatasetForDiagnostics`, `AdjustSonicCoordinates` and
`FilterDatasetForWindDirection` in between. Probe `E2Set(1, 73)` immediately
after `DefineE2Set` returns and again after `BasicStats`.

This has to be fixed before carrying per-instrument cell T/P across into FCC:
FCC reads scalar `lEx%Tcell`/`lEx%Pcell`, so widening the ex file would only
carry a value that is already wrong on the RP side.
