# Testing Strategy for TwoTimescaleResilience

## Current Issue

The TwoTimescaleResilience package spans from fast, isolated core logic to heavy integration tests invoking massive datasets, plotting libraries (CairoMakie, GeoMakie), and NetCDF generation. Currently, running `Pkg.test()` triggers full package precompilation for heavy dependencies, which often results in excessively long test runs or CI timeouts.

Because of this, `test/runtests.jl` currently has a significant portion of its test suite temporarily commented out.

## Proposed Split Strategy

To ensure rapid developer feedback without sacrificing robust end-to-end integration guarantees, the test suite should eventually be formally split into distinct categories.

### 1. Fast Default
- **Scope:** Core DEB axis math, basic response curves, mixture aggregation arithmetic, vulnerability vector calculations (standardization, simple clustering), and fast empirical parser logic.
- **Goal:** Execute in under 15 seconds. This should be the default execution path for standard `Pkg.test()`.

### 2. Extended
- **Scope:** Multi-stressor pipeline pipelines, complex condition buffering math, extended array validation, and heavy ECOTOX/AmP library queries.
- **Gated by:** `TTR_RUN_EXTENDED_TESTS=true`

### 3. Examples
- **Scope:** Ensures the demonstration scripts execute without crashing.
- **Gated by:** `TTR_RUN_EXAMPLE_TESTS=true`
- **Goal:** Do not fail fast unit tests if generated outputs are absent.

### 4. Plotting / GIS
- **Scope:** CairoMakie and GeoMakie exports, NetCDF writes (`NCDatasets`), spatial processing.
- **Gated by:** `TTR_MAKE_EXAMPLE_PLOTS=true`
- **Goal:** Keeps plotting dependency out of the core test precompilation path unless specifically requested.

### 5. Data Regeneration
- **Scope:** Archetype generation and data serialization testing.

## Environment Variable Gating

When formalizing the split, use Julia's `ENV` object within `runtests.jl` to conditionally execute files.

```julia
if get(ENV, "TTR_RUN_EXTENDED_TESTS", "false") == "true"
    @testset "Extended Integration Tests" begin
        include("test_isimip_deb_pipeline.jl")
        include("test_netcdf.jl")
    end
end
```

### Important Rule: Missing Outputs

Generated-output regression tests (like checking if an example produced a CSV) should be skipped or generate warnings, **not** fail the suite, if the execution of the example itself was gated out. Never let absent example outputs crash the fast core tests.