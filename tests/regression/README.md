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
| `base_5gas.eddyflow` | the minimal case past the historical four - one gas in slot 9 and nothing else changed, so a failure here is about crossing the boundary rather than about scale. Adds N2O on column 9. Header and rows must agree; no duplicate column names. |
| `base_dup.eddyflow` | the same species in slots 4 and 5, to check the `_2` disambiguation. |
| `base_n_gas.eddyflow` | currently eight gases across **two** analysers - the count is a property of the fixture, not of the test: the MIRO's four plus the LI-7200's CO2 and H2O, and a duplicate N2O. Exercises the capacity target (64 gases, with no per-instrument cap) and, because the two analysers measure the same species independently, catches slot cross-wiring: `CO2` and `CO2_2` must hold *different* real values, not the same one twice. |
| `base_neg.eddyflow` | `base_5gas` with `al_gas4_min` raised to 400, so COS fails the absolute-limits test. The negative fixture: exactly one gas's columns must move. Diff it against the `base_5gas` run, not against a reference. |
| `base_n_gas_ru.eddyflow` | `base_n_gas` with random uncertainty on. The only fixture that exercises `random_error_handle.f90` and `integral_turbulence_scale.f90` at all - every other one leaves `RUsetup%meth` at `none`, so those files run their `case('none')` arm and nothing else. Expect a real `RANDUNC_HF` for every gas that has a column, and `-9999` for one that does not. |
| `base_n_gas_sa.eddyflow` | `base_n_gas` with `hf_meth=ibrom_07` and `sa_mode=0`, reading the hand-authored assessment file `sa_n_gas_fitted.txt`. The only fixture where gases 5+ take a **fitted** transfer function rather than an analytic one. The file gives the first three gases `fc=1.00` and the rest `fc=0.05`, so the three resulting correction factors are far apart and cannot be confused: analytic 1.048, fitted-at-1.00 1.551, fitted-at-0.05 2.474. |
| `base_n_gas_sa_short.eddyflow` | the same, against `sa_n_gas_short.txt` - an assessment file carrying only the first two blocks, standing in for one written before the range widened. Must fall back to the analytic factors for **every** gas and say so, not correct the missing gases with a cut-off of zero. |
| `base_n_gas_bin.eddyflow` | `base_n_gas` reading the binned (co)spectra it wrote itself - `run.sh` rewrites the `SELF` token between RP and FCC. Every other fixture points `sa_bin_spectra` at a shared directory that predates the N-gas binned format and carries four gases, which makes those runs the **backward**-compatibility case and is why they are right to leave gases 5+ unassessed. This is the **forward** case: N2O, CO2_2, H2O_2 and N2O_2 go from `accepted periods=0` to `1` and from an ensemble count of 0 to 5. |

> **The `ru_*` keys reach the engine now; they never used to.**
> `ru_meth`, `ru_its_meth` and `ru_tlag_max` are declared in `EPPrjNTags`, and
> `ParseIniFile` is called with the section prefix `'Project'`, so those tags
> are only ever matched inside `[Project*]`. They have to be Project tags - RP
> and FCC both need `ru_meth`, and FCC sweeps only `FluxCorrection*`. The
> interface wrote all three into `[RawProcess_RandomUncertainty_Settings]`,
> where nothing looked for them, so `RUsetup%meth` fell to its `case default`
> of `'none'` and **random uncertainty had never actually run for any project
> the interface had saved**.
>
> The interface now writes them under `[Project]`, reads the legacy group as a
> fallback so an older file opens with its settings intact, and removes the
> stale copies on save. The engine is unchanged. The four duplicate `ru_*`
> slots in the RP tag table - which nothing read, and which are what made the
> keys look like RawProcess settings - are blanked.

> **The assessment could not be fitted from the fixtures themselves.**
> Fitting needs enough accepted half-hours to fill a class, and the regression
> dataset is three hours - too few for *any* gas, CO2 included. So the fitted
> path is exercised from the other end: `gen_sa.py` writes an assessment file
> declaring the parameters, and `sa_mode=0` feeds it back through the same
> reader the on-the-fly assessment writes for. That covers the half that was
> four-gas bounded and the half that decides whether a fitted transfer
> function can reach a fifth gas at all.
>
> Regenerate both files and both projects with:
>
> ```
> C:/Users/jonmuell/.platformio/python3/python.exe gen_sa.py
> ```

## What widening the assessment chain uncovered

The chain that fits transfer functions - bin, sort, ensemble, fit, write, read
back, apply - was bounded at the fourth gas end to end. Gases past it were
never assessed and fell through to an analytic transfer function, while the
output reported a correction factor either way and nothing said which gases
had actually been fitted.

Widening it turned up three defects that were inert *only* because the loops
stopped early. Each is the same failure class this effort keeps paying for:
widening a loop over gas slots promotes every unconfigured per-gas parameter
from never-consulted to consulted at its zero default.

1. **No month/class table past the fourth gas.** The interface exposes three
   month-grouping tables - CO2, CH4, the fourth gas - so a fifth gas kept
   class 0. That is not merely unfitted: 0 is not a valid `RegPar` index, and
   the assessment could never have fitted those gases either, because every
   month would have been written as `error`. Gases past the fourth now inherit
   CO2's grouping; the grouping bins the calendar, not the species.
2. **A phantom gas tested as present.** `var_present` was derived over every
   slot from `Flux0`, which is not reset between records, so a slot the
   project never declared held 0 rather than the error sentinel and read as
   present. It reached every `var_present`-gated loop. It surfaced as the
   *whole* spectral correction falling back to Moncrieff - the phantom has no
   class, so no cut-off to look up.
3. **Horst & Lenschow indexed the by-role instrument array by gas slot.**
   `lEx%instr` has one entry per role - CO2, H2O, CH4, the fourth gas, the
   sonic - and the code reached it as `gas - 3`, which runs off its end at the
   fifth gas and addresses an unrelated analyser before that. Under
   `-fbounds-check` this aborts the run; it is the per-gas `lEx%gas_instr`
   that carries the geometry. This is the last site of the "FCC passes
   instruments by role, not by slot" defect.

Two further things the widening made honest:

- **A short assessment file is detected rather than misread.** The reader used
  to skip the block header; a blind `read` of the next section's title line
  succeeds, so a file with too few blocks would have had its exponential-fit
  section consumed as transfer-function parameters. The header is now checked
  for `TFP`, and on a mismatch the peeked lines are put back so the sections
  below still parse.
- **An absent gas reads as unfitted, not as zero.** `RegPar` is zeroed before
  the read, and a cut-off of zero is not a missing value - it is an infinitely
  aggressive correction. It gave gases the short file did not carry a
  correction factor of 2.6. Gas blocks now start at `error`, which is what the
  readiness check keys on.

## Re-baselinings, and what each one accounted for

`out_ref` is regenerated only when a change is meant to move the output, and
only after every moved cell has been named. So far:

| change | delta against the previous reference |
|---|---|
| main record converted (`3511493`) | FLUXNET: exactly one column removed, `NUM_GAS_EXTRA`; every surviving cell byte-identical. `full_output`: 48 VM97 flag cells widened from 9 to 69 characters, each a pure extension of its old value with `'9'` (test not performed) padding |
| full output per gas | `full_output` **units row only**: 4 whitespace-only cells. Three unit labels and one flux label are `character(32)` and used to be concatenated unpadded, so those fields carried trailing blanks. Every gas now trims. No data cell moved, no column added or removed at four gases |
| per-gas cell conditions into FCC | FLUXNET row: **12 new columns**, `T_CELL_<tag>`, `PA_CELL_<tag>` and `W_PA_CELL_<tag>_COV` at four gases. Purely additive - every pre-existing column kept its value. `nMainFields` 263 -> 275 |
| FCC spectral assessment over every gas | `spectral_correction_diagnostics_adv.txt`, **one line**: the fourth gas's report line changed from `Gas 4:` to `COS:`. `GasName` returned the literal `Gas 4` for every slot past CH4, which on an eight-gas project named five different species identically; it now returns the species. No other file moved, and no flux value moved at four gases |
| spectral tags record-derived | `full_cospectra/*.csv` **header only**, two lines per file: `cov(w_gas4)` -> `cov(w_cos)` and `f_nat*cospec(w_gas4)` -> `...(w_cos)`. SpectralVarTags pinned slots 5-8 to co2/h2o/ch4/gas4, which names a position rather than a species. No numeric cell moved and no other file moved. Readers accept both spellings through LegacySpectralVarTag, so an existing full-cospectra directory still imports |
| N-gas binned (co)spectra | `binned_cospectra/*.csv` and `binned_ogives/*.csv`, **header line only**, 12 files: the fourth gas goes from `spec(COS)` to `spec(cos)` and the leading space of the list-directed write is gone. The writer named it from `SpecCol%label`, a third spelling - uppercase where every other file is lower. No data row moved and no other file moved. The reader matches case-insensitively, so the shared four-gas directory still imports |
| metadata header per gas | `eddyflow_*_metadata_*.csv`, **header line only**: the fourth gas's five run-together names become five columns. `cos_irga_tube_flowratecos_irga_kwcos_irga_ko...` was one field where the row writer emits six, so the header was 57 fields against 62 rows. No data row moved and no other file moved |

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

`base_n_gas_cell` is `base_n_gas` with cell records on **both** analysers - the
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
pressure per instrument, and ambient T/P. With it, `base_n_gas_cell` gives

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
> `base_n_gas`, `base_n_gas_ru` and `base_n_gas_cell` therefore declare no
> diagnostic record; with one, their gas checks go inert. `base_rec` keeps its
> record untouched - all four of its gases are on the MIRO, so nothing matches
> and it stays the byte-identity anchor.

**FCC: done, and it fixed a unit bug on the way.** Three per-gas groups were
added to the main record - `T_CELL_<tag>`, `PA_CELL_<tag>`, `W_PA_CELL_<tag>_COV`
- carried in **SI and read back unchanged**, the rule the `NUM_GAS_INSTR` block
already follows. On `base_n_gas_cell` they read

    MIRO gases     300.382 K   69.9983 Pa   cov 4.577e-4
    LI-7200 gases  287.827 K   94486.4 Pa   cov 0.576

> The scalars `lEx%Tcell`/`lEx%Pcell` are not only instrument 1's - they pass
> through the writer's `gain=1d0, offset=-273.15` and `gain=1d-3`, and the
> reader never inverts them. So FCC's two closed-path cell WPL terms, and the
> same three terms in the evapotranspiration block, were dividing by a
> temperature of about 27 instead of 300 and a pressure of 0.07 instead of 70.
> Latent on any project reporting mixing ratios, which is why it never showed.
> The per-gas columns replace them at every one of those sites.

## Spectral corrections reach every gas

Every gas now gets its own correction factor, computed from its own analyser:

    MIRO      CO2 / H2O / COS  1.04798     N2O  1.04877
    LI-7200   CO2_2 / H2O_2    1.04320

Four independent gates had to come out, and **each one alone was enough** to
leave `BPCF%of` at the error sentinel - so the symptom, `SCF = -9999`, looked
identical however many were fixed. That is why converting the `bpcf_*` maths
first changed nothing visible.

| gate | effect |
|---|---|
| `SetTransferFunctionsToValue` ran `u, gas4` | `BPTF` is `intent(out)`, so a skipped slot was left **undefined**; `SpectralCorrectionFactors` then found no usable band-pass value. This one made it fail under *every* method, analytic ones included - nothing to do with the cospectra file |
| `AnalyticLowPassTransferFunction`'s `case (co2, h2o, ch4, gas4)` | a gas past the fourth matched no arm and got no low-pass function |
| `CospectraMoncrieff97` copied the model to three named slots | slots 9+ kept the error value, so their cospectrum was "all error" |
| `RetrieveSensorParams` ran `co2, gas4` | path lengths and response time stayed at the sentinel. A response time of -9999 collapses the dynamic-response term: the LI-7200's gases came out at **5497** instead of 1.043 - plausible-looking only if nobody looks |

Plus the file itself: the reader matched columns against a compile-time table
whose slots 9+ were blank, and the writer's own name table was `character(4)`
and filled for eight slots. Both now come from `SpectralVarTags`, one helper,
so a name the two spell differently cannot arise. **The historical eight keep
their literal names, `gas4` included** - those strings are the shipped file
format, and renaming slot 8 to its species would change every existing
full-cospectra file.

> A gas only gets an in-situ correction if the project asks for its cospectrum.
> `base_n_gas` therefore sets `gas_N_out_full_cosp_w=1` for all eight; without
> it the column is not written and there is nothing to import.

**Still four-gas bounded, and not needed for the above:** the FCC spectral
*assessment* chain - `normalize_mean_spectra_cospectra.f90`,
`fit_tf_models.f90`, `fit_cospectral_models.f90`,
`read_spectral_assessment_file.f90` and the assessment file's own writer/reader
pair. That chain fits `RegPar` and `SA%class`, which the fully in-situ methods
consult; with them absent, `CorrectionFactorsIbrom07` and the Fratini fallback
skip a gas rather than invent a cutoff frequency, and the analytic path
supplies the factor instead. Widening it would let gases 5+ use a *fitted*
transfer function rather than an analytic one.

## Time lags: detection already worked, the PWB path did not

Under `maxcov`, which is what the fixtures use, detection was already correct
for gases past the fourth:

| gas | detected | used | window |
|---|---|---|---|
| N2O (slot 9) | 23.6 | 23.6 | [8, 25] |
| H2O_2 (slot 11) | 14.4 | 14.4 | [0, 20] |
| CO2_2 (slot 10) | 10.0 | **1.0** | [0, 10] |

CO2_2 is not a defect. Its detected lag lands exactly on `max_timelag`, is
rejected, and falls back to the nominal 1.0 - because the `.metadata` declares
that column's window as `[0, 10]` with a nominal of 1.0, while its own H2O on
the same tube detects 14.4. **The declared window is too narrow for the real
lag**; widening `col_21_max_timelag` is a site-configuration fix, not a code
one.

The PWB path *was* four-gas bounded, and is now widened: the loops in
`timelag_handle.f90`, `pwb_timelag_handle.f90`,
`adjust_timelag_opt_settings.f90`, `set_timelags.f90`, `cross_corr_test.f90`
and `optimize_timelags.f90`, the `gas < co2 .or. gas > gas4` range guards, and
the `GasLabel`/`GasIndexFromLabel` pair that names gases in the PWB time-lag
cache file. The historical four keep their literal labels there for the same
reason as in the cospectra file - the cache is read back by name, so renaming
one would orphan every cached lag.

> `donor_gas` and `GasLabel` were `character(8)`, which silently truncates a
> record-derived tag. A truncated label does not round-trip through
> `GasIndexFromLabel`, so the donor would read back as slot 0 and the
> instrument-sharing rule would drop it. Both are 32 now, and the compiler's
> truncation warning is what caught it.

> ### The widening broke the very thing it was meant to extend, and the 8-gas run caught it
>
> `SetTimelags`' `tlag_opt` branch reads `toPasGas`, the table filled from the
> time-lag optimisation file - which names the historical four. Widening the
> loop made gases 5+ take an **empty** entry, replacing the metadata's declared
> window with `[0, 0]`; every lag was then detected as zero and every flux past
> the fourth gas moved. That is precisely the "consulted at its zero default"
> failure recorded all over this document, reproduced by my own change.
>
> It is guarded now - the optimiser window is used only where it exists - and
> the guard is pinned by a check. Worth stating plainly: `base_rec` stayed
> **byte-identical through the whole episode**, because the four-gas path was
> never touched. Only the 8-gas fixture showed it. A widening gated solely on
> the byte-identity test would have shipped.

## The fabricated zero: three producers that stopped at the fourth gas

Three families reported **exactly `0.00000`** for every gas past the fourth,
in a file that declared a column for each of them. A zero is not a missing
value: a normalised cospectrum of zero, a kurtosis index of zero and a storage
term of zero are all claims about the data.

All three had the same shape - an N-gas *consumer* over a four-gas *producer*,
with the intervening array `intent(out)`, so the skipped slots were never
assigned and the writers' `/= error` guard passed whatever the stack held.

| family | producer that stopped at gas4 | symptom on `base_n_gas` |
|---|---|---|
| storage | `define_relative_separations.f90:43` | `n2o_strg`, `co2_2_strg`, `h2o_2_strg`, `n2o_2_strg` = `0.00000`; `LE_strg` and `ET` zero when the primary water sits past slot 8 |
| binned and full (co)spectra | `spectral_analysis.f90` - `AllCospectra`, `NormalizeCoSpectra`, `ExpAvrgCospectra`, the ogive binning | `spec(n2o)`, `spec(co2_2)`, `spec(h2o_2)`, `spec(n2o_2)` and their cospectra `0.00000` in all 50 bins |
| KID / ZCD, correlation difference | `kid.f90:45`, `fisher.f90:43,53,54` | `N2O_KID` … `N2O_2_KID` = `0.00000`, ZCDs `0` |

**Storage was the one with a causal chain worth following.** `Stor%of(gas)`
scales by `E2Col(gas)%Instr%height`, and `DefineRelativeSeparations` is what
turns the metadata's absolute separations into height above the anemometer.
It ran `co2, gas4`, so past slot 8 the height stayed zero and so did the
storage. Swapping records 2 and 5 (`base_h2o_late`) moved the zero with the
**slot**, not with the species - which is what identified it:

```
                h2o_strg     n2o_strg        LE_strg
base_n_gas      -0.651796     0.00000       -28.8760      (water at slot 6)
base_h2o_late    0.00000      0.0445455      -0.00000     (water at slot 9)
```

After the fix both fixtures agree by header name, and every configured gas
carries a real storage term.

> **The spectra were undefined, not merely wrong, and one accident proved it.**
> Changing `define_relative_separations.f90` - which touches no spectral code
> whatever - moved eight columns in the binned files, and *only* the eight
> belonging to gases 5+. Two runs of an unchanged tree had disagreed on them.
> They are deterministic now, and `AllCospectra`/`AllOgives` set every slot to
> `error` before filling the feasible ones, so a slot no loop reaches says
> "not performed" instead of carrying the previous period's memory.

**A second defect was sitting behind Fisher's bound.** It was passed
`E2Primes(:, 1:GHGNumVar)` - 68 columns - with `ncol = size(E2Primes, 2)`,
which is `E2NumVar`. The explicit-shape dummy therefore described half again
as many columns as were handed over, and `CorrelationMatrixNoError` read past
them. `-fbounds-check` does not catch this, and the four-gas loops never
reached far enough to trip it; widening them without fixing the count would
have walked straight in. The `KID` call site immediately above had it right.

**Gates.** `base_rec` byte-identical throughout. On `base_n_gas`: no
`spec()`/`cospec()` column may read `0.00000` for a configured gas, every
`<gas>_strg` and `*_KID`/`*_ZCD` is real or `-9999`, `base_h2o_late` matches
`base_n_gas` by header name, and two consecutive runs must diff clean. A gas
configured *without* a column - `base_rec`'s CH4 - must stay `-9999`
throughout, which is the check that the widening did not promote "absent" to
"zero".

Pinned by `static_checks/test_spectral_compute_bounds_static.py` (6 checks)
and `test_quality_test_bounds_static.py` (4), both negative-tested.

## What the water-displacement gate was actually measuring

`base_h2o_late` exists to prove that moving water off slot 6 changes nothing.
It was failing on five full-output columns - `h2o_var`, `n2o_var`,
`w/h2o_cov`, `w/n2o_cov` and `h2o_spikes` - and the five had three different
causes, only one of which was a defect.

**Four of them were a sixth fabricated zero.** `read_ex_record` copied
`stats%Cov` into `lEx%var` and `lEx%cov_w` over `u, gas4`, so every gas past
the fourth reported a variance and a w-covariance of exactly `0.00000`. The
data was in the file all along - `stats%Cov` is read over `u..lastCfg` - and
`write_out_full_fcc` already loops `firstGas..lastGas` over both. Only the
copy between them stopped early. A variance of zero says the series was
constant, which is a claim, not a gap.

It surfaced here because the swap moves H2O and N2O across the slot-8
boundary, so the two traded a real value for a zero:

```
                base_n_gas      base_h2o_late
h2o_var         0.274558E-01    0.00000        (h2o at slot 6 -> slot 9)
n2o_var         0.00000         0.206797E-06   (n2o at slot 9 -> slot 6)
```

Fixing it also lit up `co2_2_var`, `h2o_2_var`, `n2o_2_var` and their
covariances, which had been zero in *every* run of the 8-gas fixture.

The un-squaring immediately above it had the same bound, and the two had to
move together: variances are read from the file as standard deviations, so
widening only the copy would have replaced an obvious zero with a plausible
wrong number.

**The fifth was the fixture, not the code.** `h2o_spikes` differed in one
period out of six, 3 against 5. The legacy `al_h2o_min`/`al_h2o_max` keys are
keyed by *slot*: water at slot 6 is filtered against them before despiking,
water at slot 9 is not, so the two runs despike different inputs. Supplying
`gas_5_al_min=0.0` / `gas_5_al_max=40.0` - the per-gas record that says the
same thing about the species rather than about the position - makes the two
agree exactly.

That is the dual path behaving as designed, and it is worth stating plainly:
**a project saved by the interface always writes per-gas records, so it does
not depend on slot position; a hand-authored fixture that omits them does.**
`base_h2o_late` omits them deliberately, which is why the gate is quoted as
"195 of 196 columns" rather than "all".

The one column that then remains is `absolute_limits_hf`, and it is *correct*
that it differs: the flag string is one digit per slot, so moving a species
between slots moves its digit. Comparing it across the two fixtures compares
two different layouts.

> The general lesson repeats: this gate was written to test water resolution
> and spent most of its life failing for reasons that had nothing to do with
> water. Read what a red gate is actually telling you before fixing what you
> expected it to be telling you.
