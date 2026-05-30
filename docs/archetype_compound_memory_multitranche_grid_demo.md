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

## Usage

```bash
# Execute with default overrides
julia --project=. examples/archetype_compound_memory_multitranche_grid_demo.jl

# Specify overrides specifically designed for 10-year (e.g. 40 years) scaling
TTR_GRID_NX=40 TTR_GRID_NY=30 TTR_N_TRANCHES=4 TTR_TRANCHE_LENGTH_YEARS=5 julia --project=. examples/archetype_compound_memory_multitranche_grid_demo.jl
```
