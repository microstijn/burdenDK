# Archetype Compound Memory Multi-Tranche Grid Demo

## Overview

The `archetype_compound_memory_multitranche_grid_demo.jl` is a deterministic simulation integrating real AmP species archetypes and actual compound-memory records evaluated via ECOTOX across continuous timeframes ("tranches"). This example focuses on **Threshold-Free Vulnerability Regimes** to assess progressive ecological exposure without inventing threshold-breach logics.

## Architecture

This multi-tranche approach runs directly on the foundations of the 10-year grid setup (`archetype_compound_memory_10yr_grid_demo.jl`), diverging strictly to implement:
- **Tranches vs. Continuous Simulation**: Continuous memory carries physiological response ($B_t$) over multiple chunks of continuous simulated time seamlessly, skipping re-warmup phases beyond month 1.
- **Tranche Transition Modelling**: Evaluates vulnerability changes by locking onto the *Tranche 1* standardisation baseline (mean/std vector metrics) and baseline clustered centroid positions, preventing random reclustering between discrete tranches.
- **Non-Uniform Temporal Shift**: Compounds inherit deterministic spatial properties per `behavior_profile` but modify sequentially over multiple decadal tranches using structured `tranche_trajectory_profile` multipliers across a [0.85, 1.25] bandwidth logic.
- **Vulnerability Transition Analytics**: Transitions natively summarize into deterministic cluster transition statistics.

## Standardizations & Restraints
- **No real raster ingestion**: This purely utilizes algorithmic grids rather than importing `.nc` or GeoTIFF datasets.
- **No parameter inflation**: Explicit variables like $\kappa$, gain parameters, burden to margins, physiological Z_t or DEBtox scaled damage (D_t) mechanisms are completely omitted.
- **No Synergisms**: Interactions strictly adhere to package `src/mixture_aggregation.jl` structures (predominantly `grouped_ca_then_ia_axis_effects`). No matrix configurations or antagonism modifiers were embedded.

## Fixed-Reference Clustering Feature Compatibility

The multi-tranche workflow relies on fixed-reference clustering to evaluate categorical regime shifts.
This involves standardizing later tranches using standardisation parameters derived strictly from the baseline (Tranche 1).

**Feature Requirements:**
- Fixed-reference clustering requires features to hold comparable meanings across tranches.
- Absolute simulation-month features such as `month_of_max_*` (which log absolute time indices like 120, 240, etc.) are entirely excluded from fixed-reference clustering because their standardisation shifts drastically and destroys comparability.
- If timing features are desired for clustering, they should be mathematically transformed into within-tranche month, month-of-year, or a cyclic seasonal phase.
- Baseline standardisation is applied only to the subset of retained comparable features.
- Later tranches are mathematically assigned to the baseline centroids using standardized nearest-neighbor Euclidean distance.
- Centroid-assignment diagnostics (`tranche_centroid_assignment_diagnostics.csv`) are written systematically to detect degeneracy and model failure states.

## Usage

```bash
# Execute with default overrides
julia --project=. examples/archetype_compound_memory_multitranche_grid_demo.jl

# Specify overrides specifically designed for 10-year (e.g. 40 years) scaling
TTR_GRID_NX=40 TTR_GRID_NY=30 TTR_N_TRANCHES=4 TTR_TRANCHE_LENGTH_YEARS=5 julia --project=. examples/archetype_compound_memory_multitranche_grid_demo.jl
```

## Plotting multi-tranche outputs

After running the simulation, you can generate visualizations using the dedicated plotting script:

```bash
julia --project=. examples/plot_archetype_compound_memory_multitranche_grid_demo.jl
```

To override the default output directory, use the `TTR_MULTITRANCHE_DEMO_OUTPUT_DIR` environment variable:

```bash
TTR_MULTITRANCHE_DEMO_OUTPUT_DIR=output/archetype_compound_memory_multitranche_grid_demo julia --project=. examples/plot_archetype_compound_memory_multitranche_grid_demo.jl
```

### Expected Figures

The script generates heatmaps and summary figures reflecting the continuous multi-tranche progression without referencing discrete safe/unsafe thresholds. Expected figures include:

- **`cluster_maps_by_tranche.png`**: Shows vulnerability regime geography across each tranche.
- **`cluster_transition_heatmap_T1_to_final.png`**: Shows categorical movement between regimes from Tranche 1 to the final tranche. Annotates transition probabilities natively via text percentiles.
- **`feature_change_heatmap_from_baseline.png`**: Explains continuous drivers underlying the regime transitions.
- **`cluster_area_fraction_heatmap.png`**: Shows relative regime expansion/contraction over time. Also integrates percentage annotations and automatically emits dominance warnings if one cluster occupies > 95% of a given tranche.
- **`cluster_area_delta_heatmap.png`**: (Optional) Displays the change in fraction of cells for each cluster relative to baseline.
- **`tranche_distance_heatmap.png`**: (Optional) Summarizes overall cluster-distribution differences between tranches.
- **`regime_intensity_delta_map_final.png`**: (Optional) Shows where final-tranche changes are spatially concentrated via interpretive centroid-derived scores.
- **`p95_F_delta_map_final.png`** / **`p95_Q_delta_map_final.png`**: (Optional) Maps showing continuous parameter deltas between the baseline and final tranche.

### Interpretation

- Cluster maps detail *relative vulnerability regimes* representing geographical states, not absolute safe/unsafe risk classes.
- Continuous metric maps evaluate parameter deltas spatially without introducing threshold or exceedance terminology.
