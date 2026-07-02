# Engine Processing Variables Contract

This file documents the EddyFlow GUI project-file schema for processing multiple gas measurements, including duplicate gas names from the same or different IRGAs.

## INI Schema

The GUI writes all gas-processing selections in a repeated `[ProcessingVariables]` group. This group is mandatory. The legacy single-slot `[Project]` gas keys are not supported by the engine.

```ini
[ProcessingVariables]
count=N
row_1_id=co2_1_1
row_1_enabled=1
row_1_gas_col=5
row_1_gas=co2
row_1_irga=li7500_1
row_1_irga_index=1
row_1_gas_index=1
row_1_mw=44.0095
row_1_diff=0.13810
row_1_h2o_ref=h2o_1_1
row_1_cell_t=12
row_1_int_t_1=-1
row_1_int_t_2=-1
row_1_int_p=13
row_1_air_t=-1
row_1_air_p=-1
row_1_diag=20
```

Rows are 1-based. For row `i`, the full key list is:

- `row_i_id`: stable processing identifier, formatted as `gas_irgaIndex_gasIndex`, for example `co2_1_1`. If this key is absent, the engine synthesizes it from `row_i_gas`, `row_i_irga_index`, and `row_i_gas_index`. If it is present, it must match those fields exactly.
- `row_i_enabled`: `1` to process, `0` to ignore.
- `row_i_gas_col`: raw-file column number for this gas measurement; `-1` means unset.
- `row_i_gas`: normalized gas name such as `co2`, `h2o`, `ch4`, or a custom gas label.
- `row_i_irga`: stable GUI IRGA identifier.
- `row_i_irga_index`: numeric IRGA index used in the processing id.
- `row_i_gas_index`: duplicate counter for this `(gas, irga_index)` pair.
- `row_i_mw`: molecular weight in g/mol.
- `row_i_diff`: molecular diffusivity in air in cm2/s.
- `row_i_h2o_ref`: processing id of the H2O row used as reference; empty/self for H2O rows.
- `row_i_cell_t`: average cell temperature column, or `-1`.
- `row_i_int_t_1`: cell temperature in column, or `-1`.
- `row_i_int_t_2`: cell temperature out column, or `-1`.
- `row_i_int_p`: cell pressure column, or `-1`.
- `row_i_air_t`: ambient air temperature column, or `-1`.
- `row_i_air_p`: ambient pressure column, or `-1`.
- `row_i_diag`: IRGA diagnostics column for this gas row, or `-1`.

## Units

The GUI writes molecular constants in user-facing units:

- `row_i_mw`: g/mol.
- `row_i_diff`: cm2/s.

The engine converts these to existing internal units during project parsing:

- molecular weight: kg/mol.
- molecular diffusivity: m2/s.

## Removed Legacy Keys

The engine does not read gas selections from these old `[Project]` keys: `col_co2`, `col_h2o`, `col_ch4`, `col_gas4`, `col_int_t_1`, `col_int_t_2`, `col_int_p`, `col_air_t`, `col_air_p`, `col_cell_t`, `col_diag_75`, `col_diag_72`, `col_diag_77`, `gas_mw`, and `gas_diff`.

Projects without `[ProcessingVariables]` are invalid and must be updated before processing.

## Engine Identity Rule

`row_i_id` is the canonical identity for a gas measurement. It must always contain the gas name, numeric IRGA instrument id, and numeric gas-column id. Do not use `row_i_gas` as a unique key. Duplicate gas names and duplicate IRGAs are valid and must be carried independently through raw processing, corrections, diagnostics, and outputs.

Instrument model names may remain metadata in `row_i_irga`, but they are not part of core processed gas headers. Output names derive from `row_i_id`, for example `co2_1_1_flux`, `co2_1_2_mole_fraction`, `n2o_2_1_v-adv`, and `w/co2_1_1_cov`.

## Required Error Handling

Stop with a configuration error when:

- `count > MaxProcessingVariables`.
- `[ProcessingVariables]` is absent or `count <= 0`.
- An enabled row has an empty `row_i_id`.
- An enabled row has an empty `row_i_gas`.
- An enabled row has `row_i_gas_col <= 0`.
- An enabled row has `row_i_irga_index <= 0` or `row_i_gas_index <= 0`.
- A configured `row_i_id` does not match `row_i_gas`, `row_i_irga_index`, and `row_i_gas_index`.
- Two enabled rows have the same `row_i_id`.
- An H2O row has a non-empty `row_i_h2o_ref` that is not itself.
- A non-H2O `row_i_h2o_ref` points to a missing, disabled, or non-H2O row.
- A non-H2O row has an empty `row_i_h2o_ref` and there is not exactly one enabled H2O row.

Support column keys with `-1` are optional/unset and should use existing fallback behavior where scientifically valid.

## RP and FCC Handover

Both RP and FCC parse `[ProcessingVariables]` after the common `[Project]` tags. RP-to-FCC essentials and intermediate files persist the full processing-variable row list plus per-row results, using `processing_id` for gas-specific headers and products.

## Output Naming

Gas-specific output headers and intermediate products use `processing_id`. Examples:

- `co2_1_1_flux`
- `co2_1_2_flux`
- `co2_2_1_flux`
- `h2o_1_1_LE`
- `w/co2_1_1_cov`
- `co2_1_1_h2o_1_1_r_Fc_cec`
- `h2o_1_1_r_ET_cec`
