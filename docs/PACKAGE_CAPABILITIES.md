# Package Capability Inventory

## Audited date

February 22, 2025

## Source-file inventory

- `src/TwoTimescaleResilience.jl`:
  - **Role:** Main module file
  - **Exports:** Core functions, structs, and utilities for the package.
  - **Dependencies:** Imports sub-modules representing logical areas of the codebase.

- `src/amp_library.jl`:
  - **Role:** AmP species adapter and parsing module
  - **Exports:** `load_amp_species_library`, `amp_species_key`, `validate_amp_record`, `amp_record_to_deb_params`, `amp_species_deb_params`, `amp_species_profile`
  - **Dependencies:** `JSON`
  - **Notes:** Responsible for bridging user-facing species with normalized AmP records.

- `src/ecotox_library.jl`:
  - **Role:** ECOTOX parsing and empirical runtime module
  - **Exports:** `load_ecotox_library`, `validate_ecotox_record`, `ecotox_active_stress`, `ecotox_effect_to_deb_axis`, `deb_axis_index`, `ecotox_record_to_deb_burden`, `ecotox_records_to_deb_burden`, `ecotox_records_to_deb_burden_stateful!`, `ecotox_burden_to_response`, `ecotox_filter_records`, `ecotox_records_for_taxon`, `load_compound_memory_library`, `validate_compound_memory_record`, `compound_retention`, `compound_bioaccumulation_factor`, `ecotox_default_retention`, `EcotoxExposureState`, `get_internal_burden`, `set_internal_burden!`, `reset_internal_burdens!`, `update_internal_burden!`
  - **Dependencies:** `JSON`, `CSV`
  - **Notes:** Critical component for generating stress responses directly from empirical library data. Implements memory behavior.

- `src/compound_memory_warmup.jl`:
  - **Role:** Analytical compound memory spin-up helper module
  - **Exports:** `analytical_initial_burden`, `background_for_target_burden`, `analytical_periodic_initial_burden`
  - **Dependencies:** None
  - **Notes:** Utilizes exact mathematical solutions for establishing background burdens rather than requiring computational spin-up times.

- `src/deb_axes.jl`:
  - **Role:** Core DEB structures
  - **Dependencies:** None
  - **Notes:** Contains foundational math definitions mapping stress to physiological DEB axes.

- `src/reduced_deb_response.jl`:
  - **Role:** Core reduced mathematical calculations
  - **Dependencies:** None
  - **Notes:** Computes adaptive margin, restoring force, and amplification factor using derived stress.

- `src/mixture_aggregation.jl`:
  - **Role:** Mixture effect mathematical aggregation models
  - **Dependencies:** None
  - **Notes:** Provides options to implement IA, TU, and hybrid models when mapping aggregated effects.

- `src/vulnerability_feature_vectors.jl`:
  - **Role:** Threshold-free geographic vector analysis
  - **Exports:** `build_threshold_free_vulnerability_features`, `standardize_threshold_free_vulnerability_features`
  - **Dependencies:** `Statistics`
  - **Notes:** Contains strict invariants requiring names/logic to reject threshold concepts (like exceedance, above/below, limits, etc.).

- `src/vulnerability_regime_clustering.jl`:
  - **Role:** Implements vulnerability regime feature clustering
  - **Exports:** `cluster_threshold_free_vulnerability_regimes`, `summarize_threshold_free_vulnerability_clusters`, `label_threshold_free_vulnerability_regimes`
  - **Dependencies:** None
  - **Notes:** Analyzes vectors and builds regimes based on relative characteristics without strict failure logic.

- `src/vulnerability_regime_outputs.jl`:
  - **Role:** Maps spatial models and structures them into valid exports
  - **Exports:** `vulnerability_regime_cluster_map`, `vulnerability_feature_maps`, `write_vulnerability_regime_netcdf`, `write_vulnerability_regime_summary_csvs`, `vulnerability_regime_output_bundle`
  - **Dependencies:** `CSV`, `DataFrames`, `NCDatasets`
  - **Notes:** Acts as an output bundle builder integrating mapping features into export files.

## Public API inventory

- **DEB axis math:** `aggregate_deb_axis_burdens`, `axis_weights_for_species`, `compute_adaptive_margin_response`, `compute_adaptive_margin_response_from_impairment`
- **Species/AmP:** `load_amp_species_library`, `amp_species_key`, `validate_amp_record`, `amp_record_to_deb_params`, `amp_species_deb_params`, `amp_species_profile`
- **ECOTOX:** `load_ecotox_library`, `validate_ecotox_record`, `ecotox_active_stress`, `ecotox_effect_to_deb_axis`, `deb_axis_index`, `ecotox_record_to_deb_burden`, `ecotox_records_to_deb_burden`, `ecotox_records_to_deb_burden_stateful!`, `ecotox_burden_to_response`, `ecotox_filter_records`, `ecotox_records_for_taxon`
- **Compound memory:** `load_compound_memory_library`, `validate_compound_memory_record`, `compound_retention`, `compound_bioaccumulation_factor`, `ecotox_default_retention`, `EcotoxExposureState`, `get_internal_burden`, `set_internal_burden!`, `reset_internal_burdens!`, `update_internal_burden!`
- **Analytical warm-up:** `analytical_initial_burden`, `background_for_target_burden`, `analytical_periodic_initial_burden`
- **Mixture aggregation:** `mixture_contribution_diagnostics`, `aggregate_axis_mixture_effects`
- **Response modes:** (Included via standard evaluation workflows; see core equations)
- **Threshold-free features:** `build_threshold_free_vulnerability_features`
- **Standardisation:** `standardize_threshold_free_vulnerability_features`
- **Clustering:** `cluster_threshold_free_vulnerability_regimes`, `summarize_threshold_free_vulnerability_clusters`, `label_threshold_free_vulnerability_regimes`
- **Raster/NetCDF outputs:** `vulnerability_regime_cluster_map`, `vulnerability_feature_maps`, `write_vulnerability_regime_netcdf`, `write_vulnerability_regime_summary_csvs`, `vulnerability_regime_output_bundle`, `load_nc_layer`, `normalise_layer`

## Data inventory

- `AmP_Species_Library.json`: Source library of physiological DEB capabilities.
- `ECOTOX_Toxicity_Library.json`: Empirical data routing compound impact.
- `Compound_Memory_Library.csv`: Defines compound bioaccumulation and retention values.
- `AmP_Species_Archetypes.csv` and `AmP_Species_Archetypes.json`: Derived structures created by sorting parameter responses.

## Example inventory

- **Lightweight Demos:** `ecotox_amp_multiaxis_response_calibrated_demo.jl`, `synthetic_raster_demo.jl`, `ecotox_amp_multispecies_multicompound_demo.jl`
- **Memory Testing:** `ecotox_amp_monthly_memory_demo.jl`, `ecotox_amp_multispecies_multicompound_monthly_memory_demo.jl`
- **Spatial Applications:** `ecotox_amp_multispecies_multicompound_3x3_grid_demo.jl`, `ecotox_amp_synthetic_grid_mixture_demo.jl`
- **Heavy Data Extraction:** `isimip_moa_deb_3x3_demo.jl`, `nc_real_raster_deb_axes_demo.jl`

## Test inventory

- **Fast / Core Logic:** `test_deb_axes.jl`, `test_deb_axis_response.jl`, `test_mixture_aggregation.jl`, `test_deb_axes_grid.jl`
- **Current Development Focus:** `test_vulnerability_feature_vectors.jl`, `test_vulnerability_regime_clustering.jl`, `test_vulnerability_regime_outputs.jl`
- **Adapter Validation:** `test_ecotox_library.jl`, `test_amp_library.jl`
- **Extended / Heavy / Disabled tests:** `test_synthetic_raster.jl`, `test_netcdf.jl`, `test_simulation.jl`, `test_multistressor.jl`, `test_analytical_warm_start.jl`

## Capability status table

| Capability | Status |
| ---------- | ------ |
| AmP adapter | Implemented |
| ECOTOX parser | Implemented |
| ECOTOX runtime | Implemented |
| Compound memory | Implemented |
| Analytical warm-up | Implemented |
| Mixture effects | Implemented |
| Grouped CA-then-IA | Implemented |
| Response modes | Implemented |
| Archetype database | Implemented but example-only |
| Threshold-free features | Implemented |
| Standardisation | Implemented |
| Clustering | Implemented |
| Output bundle/NetCDF | Implemented |
| Real raster ingestion | Implemented but example-only (ISIMIP scripts) |
| Physiological $Z_t$ | Not implemented |
| DEBtox $D_t$ | Not implemented |
| Synergism/antagonism | Not implemented |
| Fitted interactions | Not implemented |

## Known architectural invariants

- Maintain absolute separation of "compound memory" ($B_t$) from "physiological condition memory" ($Z_t$) and DEBtox scaled damage ($D_t$).
- Maintain exact adherence to threshold-free features (`_gt_`, `threshold`, etc. strictly forbidden).
- Avoid arbitrary interaction coefficients. Implement strict "mixture-effect assumptions" instead of empirical curve tuning parameters (`kappa`, `κ`, `gain`, `response_scale`).