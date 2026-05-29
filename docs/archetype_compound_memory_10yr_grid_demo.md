# Archetype Compound Memory 10-year Grid Demo

This example script demonstrates an end-to-end integration of the `TwoTimescaleResilience` package's capabilities over a synthetic 10-year spatial domain. It utilizes realistic environmental pressure data matched to archetype physiological profiles while intentionally excluding computationally heavy or unsupported external elements.

## Features Demonstrated

*   **AmP Species Archetypes:** Loads authentic physiological profiles from the `AmP_Species_Library.json` and optionally the archetype database if generated. These are translated into limits on `DEBAxisParams` capacity.
*   **ECOTOX Toxicity Library & Compound Memory:** Filters real chemical stressor inputs from `Compound_Memory_Library.csv` with their corresponding memory characteristics ($K$ and $\rho$) matched deterministically against the `ECOTOX_Toxicity_Library.json`.
*   **Deterministic Spatial Stress Scenarios:** Computes deterministic concentrations ($C$) modeling distinct behavior profiles (e.g., persistent hotspots, point source plumes, pulse events) using 10-year monthly variations without introducing random noise.
*   **Analytical Warm-up:** Employs the `analytical_periodic_initial_burden` method to initialize baseline accumulated burden $B_0$.
*   **Mixture-Effect Assumptions:** Aggregates multi-stressor inputs strictly using the accepted rules: independent action (IA), toxic unit sum (TU), and the preferred grouped CA-then-IA model.
*   **Threshold-free Spatial Vulnerability Regimes:** Uses `build_threshold_free_vulnerability_features`, standardizes results, and clusters regimes. The clustering maps vulnerability via distinct response characteristics (e.g., growth dominated, upper envelope amplification) across spatial boundaries rather than defining arbitrary safety exceedances.

## Explicitly Excluded Concepts

To maintain mathematical soundness according to the current package maturity, this example strictly excludes:
*   Real external raster ingestion
*   Physiological condition memory carryover ($Z_t$)
*   DEBtox scaled damage models ($D_t$)
*   Curve-fitted interaction parameters (synergism / antagonism / kappa / gain)
*   Arbitrary threshold or exceedance features

## Running the Example

The demo operates as a callable script and will default its outputs to `output/archetype_compound_memory_10yr_grid_demo/`.

```bash
julia --project=. examples/archetype_compound_memory_10yr_grid_demo.jl
```

To run with environment overrides for dimensions, iteration scale, and testing properties:

```bash
TTR_GRID_NX=40 TTR_GRID_NY=30 TTR_N_YEARS=3 julia --project=. examples/archetype_compound_memory_10yr_grid_demo.jl
```
