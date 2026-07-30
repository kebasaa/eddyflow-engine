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

## Per-instrument cell T/P: found, fixed, and what it exposed

`base_8gas_cell` is `base_8gas` with cell records on **both** analysers - the
MIRO's cell_t/int_p on columns 11/12 and the LI-7200's on 23/24.

**The bug.** `FilterDatasetForDiagnostics` selected the columns an analyser's
diagnostic may invalidate with `case (co2:gas4, pi:pe)`. Both halves stopped
meaning what they said once the slots widened. `pi:pe` was instrument 1's cell
pressure through to air pressure - three slots - but with one cell block per
instrument it spans instruments 2..8 entirely, so their cell *temperatures*
began being filtered on a rule instrument 1's `tc` has never been subject to.
The second analyser's cell temperature was wiped outright, `Stats%Mean` for its
block read `-9999`, and `AirAndCellParameters` quietly fell back to instrument
1's scalars for every gas. `co2:gas4` was the mirror image: a gas past the
fourth kept every record its analyser's diagnostic rejected.

Localised with probes at `ApplyCellDiagRecords` (records land in slots 73/76
correctly), the end of `DefineE2Set` (`cell_ref` resolves to 73 correctly), and
either side of the diagnostic filter (17474 values in, 0 out).

**The fix** names quantities rather than a slot span: every gas slot, one cell
pressure per instrument, and ambient T/P. With it, `base_8gas_cell` gives

    MV_AIR_CELL   MIRO gases 35.6776    LI-7200 gases 0.0253263

where both read 35.6776 before - the LI-7200's gases now use the LI-7200's own
cell conditions.

> ### The metadata calls an AGC percentage a diagnostic word
>
> Column 28 of this dataset is `LI72_AGC`, a gain percentage reading 93.33, and
> the `.metadata` declares it `col_28_variable=diag_72`. Read as a LI-7200
> diagnostic word that gives `ibits(93,5,4) = 2`, well under the 15 the test
> demands, so **every record is rejected**.
>
> It was inert only because the LI-7200's gases sat in slots 9+, which the old
> `co2:gas4` bound never reached. Closing that gap made it bite: the analyser's
> CO2 and H2O were rejected wholesale. That is the engine behaving correctly on
> bad metadata - a LI-7200 gas in one of the first four slots has always been
> rejected the same way - but it means **any project whose metadata mislabels a
> diagnostic column will now lose its gases past the fourth, where before they
> were silently kept.** Worth a release note.
>
> `base_8gas`, `base_8gas_ru` and `base_8gas_cell` therefore declare no
> diagnostic record; with one, their gas checks go inert. `base_rec` keeps its
> record untouched - all four of its gases are on the MIRO, so nothing matches
> and it stays the byte-identity anchor.

**Still open: FCC.** `src_fcc/fluxes23.f90` reads scalar `lEx%Tcell`/`lEx%Pcell`
and `lEx%cov_w(pi)` for every gas, so under `fcc_follows` the corrected fluxes
are still computed with instrument 1's cell conditions. `lEx%Vcell` already
crosses per gas; what is missing is Tcell/Pcell and the cell-pressure
covariance, which is a four-file lockstep change on the ex file.
