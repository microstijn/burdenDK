# Package Capability Inventory

## Audited date

2026-05-29

## Source-file inventory

- `src/TwoTimescaleResilience.jl`:
  - **Role:** Main module file
  - **Exports:** All public API, grouped into specific sub-domains.
  - **Dependencies:** None directly, imports all sub-modules.
  - **Tests covering it:** None directly; acts as module root.
  - **Examples using it:** All examples.
  - **Status:** core_implemented
  - **Notes:** Central entry point for all capabilities.

- `src/amp_library.jl`:
  - **Role:** AmP species adapter and parsing module
  - **Exports:** `load_amp_species_library`, `amp_species_key`, `validate_amp_record`, `amp_record_to_deb_params`, `amp_species_deb_params`, `amp_species_profile`
  - **Dependencies:** `JSON`
  - **Tests covering it:** `test_amp_library.jl`
  - **Examples using it:** `ecotox_amp_multiaxis_response_calibrated_demo.jl`, `ecotox_amp_multispecies_multicompound_demo.jl`
  - **Status:** core_implemented_tested
  - **Notes:** Responsible for bridging user-facing species with normalized AmP records.

- `src/ECOTOXParser.jl`:
  - **Role:** ECOTOX offline parser / toxicity-library builder
  - **Exports:** `parse_ecotox_data`, `write_ecotox_library_json`, `build_ecotox_toxicity_library`, `build_ecotox_toxicity_library_multi`
  - **Dependencies:** `CSV`, `DataFrames`, `Statistics`, `JSON`
  - **Tests covering it:** None explicitly (used to build data).
  - **Examples using it:** None.
  - **Status:** core_implemented
  - **Notes:** Responsible for parsing raw offline ECOTOX ASCII files and producing the library JSON.

- `src/ecotox_library.jl`:
  - **Role:** ECOTOX runtime adapter and compound-memory runtime utilities
  - **Exports:** `load_ecotox_library`, `validate_ecotox_record`, `ecotox_active_stress`, `ecotox_effect_to_deb_axis`, `deb_axis_index`, `ecotox_record_to_deb_burden`, `ecotox_records_to_deb_burden`, `ecotox_records_to_deb_burden_stateful!`, `ecotox_burden_to_response`, `ecotox_filter_records`, `ecotox_records_for_taxon`, `load_compound_memory_library`, `validate_compound_memory_record`, `compound_retention`, `compound_bioaccumulation_factor`, `ecotox_default_retention`, `EcotoxExposureState`, `get_internal_burden`, `set_internal_burden!`, `reset_internal_burdens!`, `update_internal_burden!`
  - **Dependencies:** `JSON`, `CSV`
  - **Tests covering it:** `test_ecotox_library.jl`
  - **Examples using it:** `ecotox_amp_multispecies_multicompound_monthly_memory_demo.jl`, `ecotox_amp_multiaxis_response_calibrated_demo.jl`
  - **Status:** core_implemented_tested
  - **Notes:** Critical component for generating stress responses directly from empirical library data. Implements compound memory behavior.

- `src/compound_memory_warmup.jl`:
  - **Role:** Analytical compound memory spin-up helper module
  - **Exports:** `analytical_initial_burden`, `background_for_target_burden`, `analytical_periodic_initial_burden`
  - **Dependencies:** None
  - **Tests covering it:** `test_analytical_warm_start.jl`
  - **Examples using it:** None
  - **Status:** core_implemented_tested
  - **Notes:** Utilizes exact mathematical solutions for establishing background burdens rather than requiring computational spin-up times.

- `src/deb_axes.jl`:
  - **Role:** Core DEB structures and mathematics
  - **Exports:** None (exports managed by TwoTimescaleResilience.jl)
  - **Internal Helpers:** `default_pathogen_organic_deb_mapping`, `pairwise_axis_interaction`
  - **Dependencies:** None
  - **Tests covering it:** `test_deb_axes.jl`, `test_deb_axes_grid.jl`
  - **Examples using it:** All core examples
  - **Status:** core_implemented_tested
  - **Notes:** Contains foundational math definitions mapping stress to physiological DEB axes. Any internal pairwise-axis helper (`pairwise_axis_interaction`) should not be interpreted as fitted synergism/antagonism or chemical interaction modelling. Current mixture behaviour is represented by explicit mixture-effect assumptions in `src/mixture_aggregation.jl`.

- `src/reduced_deb_response.jl`:
  - **Role:** Core reduced mathematical calculations
  - **Exports:** `compute_adaptive_margin_response`, `compute_adaptive_margin_response_from_impairment` (via main module)
  - **Dependencies:** None
  - **Tests covering it:** `test_reduced_deb_response.jl`
  - **Examples using it:** Core workflows
  - **Status:** core_implemented_tested
  - **Notes:** Computes adaptive margin ($A_t$), restoring force ($\lambda_t$), and amplification factor ($F_t$). Physiological condition memory $Z_t$ is not an active implemented model layer. Any placeholder or optional parameter should not be treated as a validated $Z_t$ implementation.

- `src/mixture_aggregation.jl`:
  - **Role:** Mixture effect mathematical aggregation models
  - **Exports:** `aggregate_deb_axis_burdens`, `mixture_contribution_diagnostics`, `aggregate_axis_mixture_effects`
  - **Dependencies:** None
  - **Tests covering it:** `test_mixture_aggregation.jl`
  - **Examples using it:** `mixture_effect_model_overlap_demo.jl`
  - **Status:** core_implemented_tested
  - **Notes:** Implements mixture-effect assumptions (`TU`, `IA`, and `grouped_CA_then_IA`). It does NOT implement synergism or antagonism.

- `src/vulnerability_feature_vectors.jl`:
  - **Role:** Threshold-free geographic vector analysis
  - **Exports:** `build_threshold_free_vulnerability_features`, `standardize_threshold_free_vulnerability_features`
  - **Dependencies:** `Statistics`
  - **Tests covering it:** `test_vulnerability_feature_vectors.jl`
  - **Examples using it:** None explicitly running
  - **Status:** core_implemented_tested
  - **Notes:** Contains strict invariants requiring names/logic to reject threshold concepts (e.g., `_gt_`, `_lt_`, `threshold`, `exceedance`).

- `src/vulnerability_regime_clustering.jl`:
  - **Role:** Implements vulnerability regime feature clustering
  - **Exports:** `cluster_threshold_free_vulnerability_regimes`, `summarize_threshold_free_vulnerability_clusters`, `label_threshold_free_vulnerability_regimes`
  - **Dependencies:** `Statistics`
  - **Tests covering it:** `test_vulnerability_regime_clustering.jl`
  - **Examples using it:** None explicitly running
  - **Status:** core_implemented_tested
  - **Notes:** Analyzes threshold-free vectors and builds discrete vulnerability regimes.

- `src/vulnerability_regime_outputs.jl`:
  - **Role:** Maps spatial models and structures them into valid exports
  - **Exports:** `vulnerability_regime_cluster_map`, `vulnerability_feature_maps`, `write_vulnerability_regime_netcdf`, `write_vulnerability_regime_summary_csvs`, `vulnerability_regime_output_bundle`
  - **Dependencies:** `NCDatasets`, `CSV`
  - **Tests covering it:** `test_vulnerability_regime_outputs.jl`
  - **Examples using it:** None explicitly running
  - **Status:** core_implemented_tested
  - **Notes:** Acts as an output bundle builder integrating mapping features into export files. Includes strict checks against threshold-free naming violations.

- `src/netcdf.jl`:
  - **Role:** Real and synthetic raster manipulation utilities
  - **Exports:** `load_nc_layer`, `normalise_layer`
  - **Dependencies:** `NCDatasets`, `Statistics`
  - **Tests covering it:** `test_netcdf.jl`
  - **Examples using it:** ISIMIP and nc real raster examples (e.g., `examples/nc_real_raster_deb_axes_demo.jl`)
  - **Status:** core_implemented_tested
  - **Notes:** Provides basic layer loading and normalization for real external raster workflows.

## Public API inventory

- **DEB axis math / response:** `aggregate_deb_axis_burdens`, `axis_weights_for_species`, `compute_adaptive_margin_response`, `compute_adaptive_margin_response_from_impairment`
- **AmP species adapter:** `load_amp_species_library`, `amp_species_key`, `validate_amp_record`, `amp_record_to_deb_params`, `amp_species_deb_params`, `amp_species_profile`
- **ECOTOX parser/runtime:** `load_ecotox_library`, `validate_ecotox_record`, `ecotox_active_stress`, `ecotox_effect_to_deb_axis`, `deb_axis_index`, `ecotox_record_to_deb_burden`, `ecotox_records_to_deb_burden`, `ecotox_records_to_deb_burden_stateful!`, `ecotox_burden_to_response`, `ecotox_filter_records`, `ecotox_records_for_taxon`
- **Compound memory:** `load_compound_memory_library`, `validate_compound_memory_record`, `compound_retention`, `compound_bioaccumulation_factor`, `ecotox_default_retention`, `EcotoxExposureState`, `get_internal_burden`, `set_internal_burden!`, `reset_internal_burdens!`, `update_internal_burden!`
- **Analytical warm-up:** `analytical_initial_burden`, `background_for_target_burden`, `analytical_periodic_initial_burden`
- **Mixture-effect aggregation:** `mixture_contribution_diagnostics`, `aggregate_axis_mixture_effects`
- **Threshold-free feature vectors:** `build_threshold_free_vulnerability_features`
- **Standardisation:** `standardize_threshold_free_vulnerability_features`
- **Clustering:** `cluster_threshold_free_vulnerability_regimes`, `summarize_threshold_free_vulnerability_clusters`, `label_threshold_free_vulnerability_regimes`
- **Raster/NetCDF outputs:** `vulnerability_regime_cluster_map`, `vulnerability_feature_maps`, `write_vulnerability_regime_netcdf`, `write_vulnerability_regime_summary_csvs`, `vulnerability_regime_output_bundle`
- **NetCDF layer utilities:** `load_nc_layer`, `normalise_layer`
- **Important internal helpers / example support (exported):** `plot_scenario_comparison`, `plot_grid`, `plot_background_layers`, `plot_amplification_grid`, `synthetic_background_layers`, `run_synthetic_raster_demo`
- **ISIMIP/spatial utilities:** Various internal unexported helpers, but relies on `netcdf.jl` tools.

## Data inventory

- `data/AmP_Species_Library.json`: Source library of physiological DEB capabilities.
- `data/ECOTOX_Toxicity_Library.json`: Empirical data routing compound impact.
- `data/Compound_Memory_Library.csv`: Defines compound bioaccumulation and retention values.
- `data/AmP_Species_Archetypes.csv` and `data/AmP_Species_Archetypes.json`: Derived databases mapping species onto generalized archetypes.

## Example inventory

- **Fast / Quick:**
  - `synthetic_raster_demo.jl`
  - `ecotox_amp_multiaxis_response_calibrated_demo.jl`
  - `species_comparison_3x3_demo.jl`
  - `mixture_effect_model_overlap_demo.jl`
- **Medium:**
  - `ecotox_amp_multispecies_multicompound_demo.jl`
  - `ecotox_amp_monthly_memory_demo.jl`
  - `ecotox_amp_multispecies_multicompound_monthly_memory_demo.jl`
  - `ecotox_amp_multispecies_multicompound_3x3_grid_demo.jl`
  - `ecotox_amp_synthetic_grid_mixture_demo.jl`
- **Heavy / Extended:**
  - `isimip_moa_deb_3x3_demo.jl`
  - `nc_real_raster_deb_axes_demo.jl`
  - `nc_monthly_longterm_isimip_moa_deb_inspection.jl`
  - `nc_monthly_longterm_isimip_moa_deb_inspection_real.jl`
  - `nc_monthly_longterm_isimip_moa_deb_inspection_europe_geomakie.jl`
  - `build_amp_species_archetype_database.jl`

## Test inventory

- **Fast / Core tests:**
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
- **Extended / Heavy / Gated tests (Temporarily disabled by default):**
  - `test_response_modes.jl` (May contain generated-output checks that should be gated)
  - `test_synthetic_raster.jl`
  - `test_netcdf.jl`
  - `test_simulation.jl`
  - `test_multistressor.jl`
  - `test_deb_pipeline.jl`
  - `test_examples.jl`
  - `test_plotting.jl`

## Capability status table

| Capability | Status | Evidence | Notes |
| ---------- | ------ | -------- | ----- |
| AmP runtime adapter | core_implemented_tested | `src/amp_library.jl`, `test_amp_library.jl` | |
| ECOTOX offline parser/library builder | core_implemented | `src/ECOTOXParser.jl` | Rebuilds JSON from ASCII. |
| ECOTOX runtime adapter | core_implemented_tested | `src/ecotox_library.jl`, `test_ecotox_library.jl` | |
| Compound memory | core_implemented_tested | `src/ecotox_library.jl`, `test_ecotox_library.jl` | Stateful burden logic. |
| Analytical compound-memory warm-up | core_implemented_tested | `src/compound_memory_warmup.jl`, `test_analytical_warm_start.jl` | |
| Mixture-effect assumptions | core_implemented_tested | `src/mixture_aggregation.jl`, `test_mixture_aggregation.jl` | TU, IA |
| Grouped CA-then-IA | core_implemented_tested | `src/mixture_aggregation.jl`, `test_mixture_aggregation.jl` | Preferred aggregation method. |
| Response modes | core_implemented_tested | `src/reduced_deb_response.jl`, implicit in tests | |
| Archetype database | core_implemented_tested | `data/AmP_Species_Archetypes.csv`, `data/AmP_Species_Archetypes.json`, `examples/build_amp_species_archetype_database.jl`, `test_amp_species_archetypes.jl` | |
| Threshold-free features | core_implemented_tested | `src/vulnerability_feature_vectors.jl`, `test_vulnerability_feature_vectors.jl` | |
| Standardisation | core_implemented_tested | `src/vulnerability_feature_vectors.jl`, `test_vulnerability_feature_vectors.jl` | |
| Clustering | core_implemented_tested | `src/vulnerability_regime_clustering.jl`, `test_vulnerability_regime_clustering.jl` | |
| Regime output bundle / NetCDF | core_implemented_tested | `src/vulnerability_regime_outputs.jl`, `test_vulnerability_regime_outputs.jl` | Spatial regime structures. |
| NetCDF layer utilities | core_implemented_tested | `src/netcdf.jl`, `test_netcdf.jl` | Basic loading/normalization. |
| Synthetic raster examples | example_only | `examples/synthetic_raster_demo.jl` | |
| Real external raster demonstration examples | example_only | `examples/nc_real_raster_deb_axes_demo.jl` | |
| Stable reusable real external raster ingestion pipeline | partial | `src/netcdf.jl` | Mostly demonstrated in scripts, missing generalized robust API. |
| Physiological $Z_t$ | not_implemented | Checks reject it; noted as explicitly not implemented in source. | |
| DEBtox scaled damage $D_t$ | not_implemented | No source code evidence found. | |
| Synergism/antagonism | not_implemented | `src/mixture_aggregation.jl` | Explicitly excluded. |
| Fitted interactions | not_implemented | | No arbitrary tuning parameters implemented. |

## Known architectural invariants

- **Compound memory vs Physiological condition:** Maintain absolute separation of "chemical memory" ($B_t$) from "physiological condition memory" ($Z_t$) and DEBtox scaled damage ($D_t$). $Z_t$ and $D_t$ are not implemented.
- **Threshold-free definitions:** Maintain exact adherence to threshold-free features. Naming or logic using patterns like `_gt_`, `_lt_`, `threshold`, `exceedance`, `above`, or `below` are strictly forbidden and assert failure.
- **Mixture assumptions:** The package implements mathematical mixture-effect assumptions (`TU`, `IA`, `grouped_CA_then_IA`), not arbitrary curve tuning interaction models (e.g., synergism/antagonism).
- **Prohibited tuning parameters:** Mathematical tuning parameters named or behaving like `kappa`, `κ`, `gain`, `response_scale`, or `burden_to_margin_multiplier` are strictly avoided.
