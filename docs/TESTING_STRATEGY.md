# Testing Strategy for TwoTimescaleResilience

*Audited Date: 2026-05-29*

## Current Issue: Heavy Precompilation and Long Run Times

The framework integrates deeply with GIS mapping (`GeoMakie`), plotting (`CairoMakie`), NetCDF IO (`NCDatasets`), and large JSON datasets. This results in incredibly heavy dependency loads and precompilation times. Running `Pkg.test()` on the full suite during active development takes excessively long, disrupting iterative development cycles.

Furthermore, some generated-output example tests may fail fast unit tests merely because historical output reference files are missing on fresh environments. This shouldn't happen.

## Proposed Test Split

To solve this, the tests are conceptually split into multiple tiers. This ensures rapid feedback on core mathematical logic while preserving comprehensive integration capability.

### 1. Fast Default Tests
These evaluate the pure, dependency-light mathematical engine. They run in milliseconds and should be executed frequently.
- **Includes:** `test_deb_axes.jl`, `test_deb_axis_response.jl`, `test_mixture_aggregation.jl`, `test_ecotox_library.jl`, `test_vulnerability_feature_vectors.jl`, `test_analytical_warm_start.jl`

### 2. Extended Integration
Tests covering deeper adapters, full runtime execution pipelines, and time-series memory behaviors.
- **Includes:** `test_deb_pipeline.jl`, `test_simulation.jl`, `test_multistressor.jl`

### 3. Examples and Plotting
Tests that actually produce maps and complex terminal outputs.
- **Includes:** `test_examples.jl`, `test_plotting.jl`, `test_response_modes.jl` (mixed tests)

**Note on `test_response_modes.jl`:** This file currently mixes fast response-mode unit tests with extended generated-output regression checks for CSVs. Future cleanup should either split the fast component from the extended component, or gate the generated-output sections behind an environment variable.

### 4. Data Regeneration / Outputs
Tests that invoke NetCDF file construction or check historical outputs. Generated-output tests should not fail the fast suite merely because expected output files are absent.
- **Includes:** `test_netcdf.jl`, `test_synthetic_raster.jl`

## Environment Variable Gating

Currently, heavy tests are manually commented out in `test/runtests.jl`. Moving forward, the recommended strategy is to implement gating via environment variables so the CI server can opt-in to heavy validation, while local developers run the fast subset by default.

**Proposed environment variables (Not currently implemented):**
- `TTR_RUN_EXTENDED_TESTS=true`
- `TTR_RUN_EXAMPLE_TESTS=true`
- `TTR_RUN_OUTPUT_REGRESSION_TESTS=true`
- `TTR_MAKE_EXAMPLE_PLOTS=true`

## Recommended Commands for Fast Development

During active logic changes, do not run the full suite. Instead run individual files:

```bash
# Evaluate core math invariants
julia --project=. test/test_deb_axes.jl

# Evaluate vulnerability features
julia --project=. test/test_vulnerability_feature_vectors.jl
```

To run the full fast suite as currently structured:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Recommended Commands for Full Validation

Before pushing significant feature branches, simulating the CI environment and running the full integration suite is required.

```bash
TTR_RUN_EXTENDED_TESTS=true TTR_RUN_EXAMPLE_TESTS=true TTR_RUN_OUTPUT_REGRESSION_TESTS=true julia --project=. -e 'using Pkg; Pkg.test()'
```
*(Note: These variables will need to be checked in `test/runtests.jl` before this behaves automatically.)*