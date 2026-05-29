# Package Capability Inventory

## Audited date

2026-05-29

## Source-file inventory

### `src/TwoTimescaleResilience.jl`

- **Role:** Main module file and package entry point.
- **Exports:** Public API exported from the module; grouped by domain in the [Public API inventory](#public-api-inventory).
- **Dependencies:** Module root; includes package source files and manages exports. Dependency usage is mostly in included source files.
- **Tests covering it:** None directly; acts as module root.
- **Examples using it:** All examples that call `using TwoTimescaleResilience`.
- **Status:** core_implemented
- **Notes:** Central entry point for package capabilities. Avoid treating this file as owning model logic directly; it mainly wires together source files and exports.

### `src/amp_library.jl`

- **Role:** AmP species adapter and parsing module.
- **Exports:** `load_amp_species_library`, `amp_species_key`, `validate_amp_record`, `amp_record_to_deb_params`, `amp_species_deb_params`, `amp_species_profile`.
- **Dependencies:** `JSON`.
- **Tests covering it:** `test_amp_library.jl`.
- **Examples using it:** `ecotox_amp_multiaxis_response_calibrated_demo.jl`, `ecotox_amp_multispecies_multicompound_demo.jl`.
- **Status:** core_implemented_tested
- **Notes:** Bridges user-facing species names/keys with normalized AmP records and package-level species capacity parameters.

### `src/ECOTOXParser.jl`

- **Role:** ECOTOX offline parser / toxicity-library builder.
- **Exports:** `parse_ecotox_data`, `write_ecotox_library_json`, `build_ecotox_toxicity_library`, `build_ecotox_toxicity_library_multi`.
- **Dependencies:** `CSV`, `DataFrames`, `Statistics`, `JSON`.
- **Tests covering it:** `test_ecotox.jl` if parser functions are tested there; otherwise parser test coverage should be considered incomplete.
- **Examples using it:** None identified.
- **Status:** core_implemented, or core_implemented_tested if `test_ecotox.jl` directly exercises parser functions.
- **Notes:** Responsible for parsing raw/offline ECOTOX ASCII files and producing toxicity-library JSON. This is distinct from `src/ecotox_library.jl`, which is the runtime adapter for already-built ECOTOX libraries.

### `src/ecotox_library.jl`

- **Role:** ECOTOX runtime adapter and compound-memory runtime utilities.
- **Exports:** `load_ecotox_library`, `validate_ecotox_record`, `ecotox_active_stress`, `ecotox_effect_to_deb_axis`, `deb_axis_index`, `ecotox_record_to_deb_burden`, `ecotox_records_to_deb_burden`, `ecotox_records_to_deb_burden_stateful!`, `ecotox_burden_to_response`, `ecotox_filter_records`, `ecotox_records_for_taxon`, `load_compound_memory_library`, `validate_compound_memory_record`, `compound_retention`, `compound_bioaccumulation_factor`, `ecotox_default_retention`, `EcotoxExposureState`, `get_internal_burden`, `set_internal_burden!`, `reset_internal_burdens!`, `update_internal_burden!`.
- **Dependencies:** `JSON`, `CSV`.
- **Tests covering it:** `test_ecotox_library.jl`.
- **Examples using it:** `ecotox_amp_multispecies_multicompound_monthly_memory_demo.jl`, `ecotox_amp_multiaxis_response_calibrated_demo.jl`.
- **Status:** core_implemented_tested
- **Notes:** Runtime component for loading ECOTOX toxicity libraries, validating records, computing active stress, routing effect codes to DEB-informed axes, computing stateless/stateful burdens, and supporting compound memory behaviour.

### `src/compound_memory_warmup.jl`

- **Role:** Analytical compound-memory spin-up helper module.
- **Exports:** `analytical_initial_burden`, `background_for_target_burden`, `analytical_periodic_initial_burden`.
- **Dependencies:** None.
- **Tests covering it:** `test_analytical_warm_start.jl`.
- **Examples using it:** None identified.
- **Status:** core_implemented_tested
- **Notes:** Provides exact mathematical helpers for establishing background burden states without computational spin-up. These utilities initialise chemical memory `B_t`; they are not physiological condition memory `Z_t` and not DEBtox scaled damage `D_t`.

### `src/deb_axes.jl`

- **Role:** Core DEB-axis structures and mathematics.
- **Public API exported from module:** See `src/TwoTimescaleResilience.jl` and the [Public API inventory](#public-api-inventory).
- **Locally defined/internal helpers:** `default_pathogen_organic_deb_mapping`, `pairwise_axis_interaction`, if present.
- **Dependencies:** None.
- **Tests covering it:** `test_deb_axes.jl`, `test_deb_axes_grid.jl`.
- **Examples using it:** Core examples.
- **Status:** core_implemented_tested
- **Notes:** Contains foundational utilities for mapping burdens to physiological DEB-informed axes. Any internal pairwise-axis helper should not be interpreted as fitted synergism/antagonism or chemical interaction modelling. Current mixture behaviour is represented by explicit mixture-effect assumptions in `src/mixture_aggregation.jl`.

### `src/reduced_deb_response.jl`

- **Role:** Core reduced response calculations.
- **Exports:** `compute_adaptive_margin_response`, `compute_adaptive_margin_response_from_impairment` via the main module, if exported there.
- **Dependencies:** None.
- **Tests covering it:** `test_reduced_deb_response.jl`, and response-mode unit tests in `test_response_modes.jl` where applicable.
- **Examples using it:** Core workflows.
- **Status:** core_implemented_tested
- **Notes:** Computes adaptive margin `A_t`, restoring force `lambda_t`, and amplification factor `F_t`. Physiological condition memory `Z_t` is not an active validated model layer. Any placeholder or optional parameter should not be treated as a validated `Z_t` implementation.

### `src/mixture_aggregation.jl`

- **Role:** Mixture-effect mathematical aggregation models.
- **Exports:** `aggregate_deb_axis_burdens`, `mixture_contribution_diagnostics`, `aggregate_axis_mixture_effects`.
- **Dependencies:** None.
- **Tests covering it:** `test_mixture_aggregation.jl`.
- **Examples using it:** `mixture_effect_model_overlap_demo.jl`.
- **Status:** core_implemented_tested
- **Notes:** Implements mixture-effect assumptions such as toxic-unit summation, independent action, and grouped concentration-addition-then-independent-action. It does not implement synergism, antagonism, fitted interaction coefficients, or arbitrary interaction matrices.

### `src/vulnerability_feature_vectors.jl`

- **Role:** Threshold-free spatial vulnerability feature-vector construction and feature standardisation.
- **Exports:** `build_threshold_free_vulnerability_features`, `standardize_threshold_free_vulnerability_features`.
- **Dependencies:** `Statistics`.
- **Tests covering it:** `test_vulnerability_feature_vectors.jl`.
- **Examples using it:** None explicitly identified.
- **Status:** core_implemented_tested
- **Notes:** Builds continuous, threshold-free feature matrices from gridded response arrays and standardises/drops features for downstream clustering. Feature names/logic reject arbitrary threshold concepts such as `_gt_`, `_lt_`, `threshold`, `exceedance`, `above`, and `below` if those guards are present in source.

### `src/vulnerability_regime_clustering.jl`

- **Role:** Threshold-free vulnerability-regime clustering.
- **Exports:** `cluster_threshold_free_vulnerability_regimes`, `summarize_threshold_free_vulnerability_clusters`, `label_threshold_free_vulnerability_regimes`.
- **Dependencies:** `Statistics`.
- **Tests covering it:** `test_vulnerability_regime_clustering.jl`.
- **Examples using it:** None explicitly identified.
- **Status:** core_implemented_tested
- **Notes:** Clusters standardised threshold-free feature vectors into discrete vulnerability regimes. Cluster labels are relative regime descriptions, not safe/unsafe classes and not regulatory exceedance categories.

### `src/vulnerability_regime_outputs.jl`

- **Role:** Vulnerability-regime spatial output and export utilities.
- **Exports:** `vulnerability_regime_cluster_map`, `vulnerability_feature_maps`, `write_vulnerability_regime_netcdf`, `write_vulnerability_regime_summary_csvs`, `vulnerability_regime_output_bundle`.
- **Dependencies:** `NCDatasets`, `CSV`; include `DataFrames` only if the source actually imports/uses it.
- **Tests covering it:** `test_vulnerability_regime_outputs.jl`.
- **Examples using it:** None explicitly identified.
- **Status:** core_implemented_tested
- **Notes:** Builds cluster maps, feature maps, NetCDF outputs, and CSV summaries for threshold-free vulnerability regimes. Includes or should include checks against threshold/exceedance naming violations. If DataFrames is used in this source file, future cleanup may consider Tables-compatible rows to keep the core output layer lighter.

### `src/netcdf.jl`

- **Role:** NetCDF layer loading and normalisation utilities.
- **Exports:** `load_nc_layer`, `normalise_layer`.
- **Dependencies:** `NCDatasets`, `Statistics`.
- **Tests covering it:** `test_netcdf.jl`.
- **Examples using it:** ISIMIP and NetCDF raster examples, such as `examples/nc_real_raster_deb_axes_demo.jl`, if present.
- **Status:** core_implemented_tested
- **Notes:** Provides basic NetCDF layer loading and normalisation utilities. This is not by itself a complete stable real-raster ingestion workflow. Tests may be operationally heavier because NetCDF/raster dependencies can increase precompilation or IO cost.

## Public API inventory

The public API inventory should include only names exported from `src/TwoTimescaleResilience.jl`. If a function is defined but not exported, document it separately as an internal helper or example-support function.

- **DEB axis math / response:** `aggregate_deb_axis_burdens`, `axis_weights_for_species`, `compute_adaptive_margin_response`, `compute_adaptive_margin_response_from_impairment`.
- **AmP species adapter:** `load_amp_species_library`, `amp_species_key`, `validate_amp_record`, `amp_record_to_deb_params`, `amp_species_deb_params`, `amp_species_profile`.
- **ECOTOX parser/runtime:** `parse_ecotox_data`, `write_ecotox_library_json`, `build_ecotox_toxicity_library`, `build_ecotox_toxicity_library_multi`, `load_ecotox_library`, `validate_ecotox_record`, `ecotox_active_stress`, `ecotox_effect_to_deb_axis`, `deb_axis_index`, `ecotox_record_to_deb_burden`, `ecotox_records_to_deb_burden`, `ecotox_records_to_deb_burden_stateful!`, `ecotox_burden_to_response`, `ecotox_filter_records`, `ecotox_records_for_taxon`, if these are exported.
- **Compound memory:** `load_compound_memory_library`, `validate_compound_memory_record`, `compound_retention`, `compound_bioaccumulation_factor`, `ecotox_default_retention`, `EcotoxExposureState`, `get_internal_burden`, `set_internal_burden!`, `reset_internal_burdens!`, `update_internal_burden!`.
- **Analytical warm-up:** `analytical_initial_burden`, `background_for_target_burden`, `analytical_periodic_initial_burden`.
- **Mixture-effect aggregation:** `mixture_contribution_diagnostics`, `aggregate_axis_mixture_effects`, `aggregate_deb_axis_burdens`.
- **Threshold-free feature vectors:** `build_threshold_free_vulnerability_features`.
- **Standardisation:** `standardize_threshold_free_vulnerability_features`.
- **Clustering:** `cluster_threshold_free_vulnerability_regimes`, `summarize_threshold_free_vulnerability_clusters`, `label_threshold_free_vulnerability_regimes`.
- **Regime raster/NetCDF outputs:** `vulnerability_regime_cluster_map`, `vulnerability_feature_maps`, `write_vulnerability_regime_netcdf`, `write_vulnerability_regime_summary_csvs`, `vulnerability_regime_output_bundle`.
- **NetCDF layer utilities:** `load_nc_layer`, `normalise_layer`.
- **Exported plotting/example-support helpers, if exported:** `plot_scenario_comparison`, `plot_grid`, `plot_background_layers`, `plot_amplification_grid`, `synthetic_background_layers`, `run_synthetic_raster_demo`.
- **Internal/example-support helpers:** Any functions not exported from `src/TwoTimescaleResilience.jl` should not be listed as public API.
- **ISIMIP/spatial utilities:** Document exported names only. If utilities are internal or example-only, classify them as such.

## Data inventory

- `data/AmP_Species_Library.json`: Source library of AmP-derived physiological/DEB capacity records.
- `data/ECOTOX_Toxicity_Library.json`: Derived/runtime ECOTOX toxicity library used for stressor pressure records.
- `data/Compound_Memory_Library.csv`: Compound memory data defining retention and bioaccumulation/internal magnification values.
- `data/AmP_Species_Archetypes.csv` and `data/AmP_Species_Archetypes.json`: Derived species-archetype databases if present. Treat as derived artifacts generated from AmP response-capacity diagnostics, not raw AmP source data.

## Example inventory

### Fast / quick examples

- `synthetic_raster_demo.jl`
- `ecotox_amp_multiaxis_response_calibrated_demo.jl`
- `species_comparison_3x3_demo.jl`
- `mixture_effect_model_overlap_demo.jl`

### Medium examples

- `ecotox_amp_multispecies_multicompound_demo.jl`
- `ecotox_amp_monthly_memory_demo.jl`
- `ecotox_amp_multispecies_multicompound_monthly_memory_demo.jl`
- `ecotox_amp_multispecies_multicompound_3x3_grid_demo.jl`
- `ecotox_amp_synthetic_grid_mixture_demo.jl`

### Heavy / extended examples

- `isimip_moa_deb_3x3_demo.jl`
- `nc_real_raster_deb_axes_demo.jl`
- `nc_monthly_longterm_isimip_moa_deb_inspection.jl`
- `nc_monthly_longterm_isimip_moa_deb_inspection_real.jl`
- `nc_monthly_longterm_isimip_moa_deb_inspection_europe_geomakie.jl`
- `build_amp_species_archetype_database.jl`

Example heaviness should be rechecked when dependencies, grid sizes, or output generation change.

## Test inventory

### Fast / core tests

- `test_deb_axes.jl`
- `test_deb_axis_response.jl`
- `test_mixture_aggregation.jl`
- `test_deb_axes_grid.jl`
- `test_vulnerability_feature_vectors.jl`
- `test_vulnerability_regime_clustering.jl`
- `test_vulnerability_regime_outputs.jl`
- `test_ecotox_library.jl`
- `test_amp_library.jl`
- `test_analytical_warm_start.jl`

### Mixed tests

- `test_response_modes.jl`
  - Fast/core portion: response-mode unit tests.
  - Extended/generated-output portion: checks for generated multiaxis example CSVs. These output-dependent sections should be split or gated so missing generated files do not fail the fast suite.

### Extended / heavy / gated tests

- `test_synthetic_raster.jl`
- `test_netcdf.jl`
- `test_simulation.jl`
- `test_multistressor.jl`
- `test_deb_pipeline.jl`
- `test_examples.jl`
- `test_plotting.jl`

Some tests in this group may validate real core functionality but are operationally heavier because they trigger plotting, NetCDF, raster, example-output, or larger integration paths.

## Capability status table

| Capability | Status | Evidence | Notes |
| ---------- | ------ | -------- | ----- |
| AmP runtime adapter | core_implemented_tested | `src/amp_library.jl`, `test_amp_library.jl` | Loads AmP records and converts them to package species/capacity objects. |
| ECOTOX offline parser/library builder | core_implemented or core_implemented_tested | `src/ECOTOXParser.jl`; `test_ecotox.jl` if parser-tested | Rebuilds toxicity JSON from raw/offline ECOTOX files. |
| ECOTOX runtime adapter | core_implemented_tested | `src/ecotox_library.jl`, `test_ecotox_library.jl` | Runtime loading, validation, active stress, routing, burden, response. |
| Compound memory | core_implemented_tested | `src/ecotox_library.jl`, `test_ecotox_library.jl` | Stateful burden logic with `rho` and `K`. |
| Analytical compound-memory warm-up | core_implemented_tested | `src/compound_memory_warmup.jl`, `test_analytical_warm_start.jl` | Exact constant-background, inverse-target, and periodic-cycle helpers. |
| Mixture-effect assumptions | core_implemented_tested | `src/mixture_aggregation.jl`, `test_mixture_aggregation.jl` | Toxic-unit summation and independent-action style assumptions. |
| Grouped CA-then-IA | core_implemented_tested | `src/mixture_aggregation.jl`, `test_mixture_aggregation.jl` | Preferred grouped mixture-effect assumption when effect-code grouping exists. |
| Response modes | core_implemented_tested | `src/reduced_deb_response.jl`, `test_response_modes.jl` fast/unit sections | Generated-output sections of `test_response_modes.jl` should be gated or split. |
| Archetype database | core_implemented_tested, core_implemented, partial, or planned | `data/AmP_Species_Archetypes.csv`, `data/AmP_Species_Archetypes.json`, `examples/build_amp_species_archetype_database.jl`, `test_amp_species_archetypes.jl` if present | Derived species-archetype artifact; status depends on actual file/test presence. |
| Threshold-free features | core_implemented_tested | `src/vulnerability_feature_vectors.jl`, `test_vulnerability_feature_vectors.jl` | Continuous features only; no arbitrary exceedance features. |
| Standardisation | core_implemented_tested | `src/vulnerability_feature_vectors.jl`, `test_vulnerability_feature_vectors.jl` | Feature dropping and standardisation for clustering. |
| Clustering | core_implemented_tested | `src/vulnerability_regime_clustering.jl`, `test_vulnerability_regime_clustering.jl` | Deterministic threshold-free regime clustering. |
| Regime output bundle / NetCDF | core_implemented_tested | `src/vulnerability_regime_outputs.jl`, `test_vulnerability_regime_outputs.jl` | Writes vulnerability-regime cluster maps and selected feature rasters to NetCDF/CSV. |
| NetCDF layer utilities | core_implemented_tested | `src/netcdf.jl`, `test_netcdf.jl` | Basic layer loading/normalisation utilities; tests may be operationally heavier. |
| Synthetic raster examples | example_only | `examples/synthetic_raster_demo.jl` | Demonstration workflow, not source of core model equations. |
| Real external raster demonstration examples | example_only | `examples/nc_real_raster_deb_axes_demo.jl` | Demonstrates use of raster/NetCDF inputs if present. |
| Stable reusable real external raster ingestion pipeline | partial | `src/netcdf.jl`, raster examples | Basic utilities exist, but generalized robust real-raster ingestion is not yet the main stable workflow. |
| Physiological condition memory `Z_t` | not_implemented | No active validated `Z_t` model layer found. | Keep distinct from chemical memory `B_t`. |
| DEBtox scaled damage `D_t` | not_implemented | No source code evidence found. | Future TKTD extension only; do not conflate with `B_t` or `E_axis`. |
| Synergism/antagonism | not_implemented | `src/mixture_aggregation.jl` documents explicit mixture-effect assumptions. | TU, IA, and grouped CA-then-IA are not fitted interaction models. |
| Fitted interactions | not_implemented | No source evidence found for fitted interaction coefficients or arbitrary interaction matrices. | Current mixture layer is assumption-based. |

## Known architectural invariants

- **Compound memory vs physiological condition:** Maintain absolute separation of chemical memory `B_t` from physiological condition memory `Z_t` and DEBtox scaled damage `D_t`. `Z_t` and `D_t` are not implemented as active validated model layers.
- **Threshold-free definitions:** Threshold-free feature construction must not introduce arbitrary exceedance features. Naming or logic using patterns such as `_gt_`, `_lt_`, `threshold`, `exceedance`, `above`, or `below` should fail validation if the guards are present in source.
- **Mixture assumptions:** The package implements mathematical mixture-effect assumptions such as TU, IA, and grouped CA-then-IA, not arbitrary curve-tuned chemical interaction models.
- **Prohibited tuning parameters:** Mathematical tuning parameters named or behaving like `kappa`, `κ`, `gain`, `response_scale`, or `burden_to_margin_multiplier` should remain absent unless a separate, explicitly scoped model extension is introduced and reviewed.
