# AI Context for TwoTimescaleResilience

*Audited Date: 2026-05-29*

## Do first

- Inspect the actual code before assuming functionality is present. Look at files to determine API and state. Do not infer a stable real external raster ingestion pipeline from example scripts alone.
- Run targeted tests, not the full heavy suite, during development (see `docs/TESTING_STRATEGY.md`).
- Preserve the architecture: look closely at the data flow before editing components.
- Do not rewrite core math unless explicitly requested by the user.
- Consult `docs/ARCHITECTURE_GRAPH.md` for component relationships.

## Project identity

TwoTimescaleResilience is a framework for modelling the erosion of response capacity caused by accumulated environmental pressures, allowing researchers to evaluate vulnerability over spatial domains without relying on binary threshold-exceedance models.

## Current architecture

The framework uses a **capacity–pressure–memory** architecture. Capacity limits are defined by AmP species parameters; environmental pressures are informed by ECOTOX toxicity data; and temporal accumulation is modelled via a chemical memory component ($B_t$). The system evaluates multiple overlapping pressures using specific mixture-effect assumptions.

## Core files and what they do

- `src/TwoTimescaleResilience.jl`: Exports public APIs and controls logic domain.
- `src/deb_axes.jl` and `src/reduced_deb_response.jl`: Define the foundational impairment math, calculating $Q_t, A_t, \lambda_t,$ and $F_t$.
- `src/ecotox_library.jl`: Parses ECOTOX databases, manages empirical memory state $B_t$, and scales ambient exposures to active stress ($x_t$).
- `src/mixture_aggregation.jl`: Aggregates active stress using IA, TU, or grouped logic.
- `src/amp_library.jl`: Interacts with AmP library parameter stores.
- `src/vulnerability_feature_vectors.jl` and `src/vulnerability_regime_clustering.jl`: Build and cluster spatially explicit geographic vulnerability models without failure thresholds.

## What not to change

- Do not implement synergism or antagonism interactions. Do not call mixture-effect assumptions "interactions". Use standard mixture-effect assumptions.
- Do not implement arbitrary empirical curve tuning parameters such as `kappa`, `κ`, `gain`, `response_scale`, or `burden_to_margin_multiplier`.
- Do not implement threshold failures. Do not add thresholds to threshold-free features. Avoid naming or using variables like `_gt_`, `_lt_`, `threshold`, `exceedance`, `above`, or `below` in spatial feature/regime models.
- Maintain a strict boundary between chemical memory ($B_t$), physiological condition memory ($Z_t$), and DEBtox scaled damage ($D_t$). Keep $B_t$, $Z_t$, and $D_t$ distinct.
- Do not use DataFrames for lightweight internal data structures.
- Do not alter testing execution behavior arbitrarily.

## Common next tasks

- Implementing physiological condition memory ($Z_t$) (Future/Planned)
- Creating robust stable real-raster ingestion pipelines (Future/Planned)
- Expanding diagnostic outputs and visualizations.

## Test strategy

Full test suites may trigger massive precompilation tasks for dependencies like CairoMakie, GeoMakie, and NCDatasets. Avoid doing so unless specifically evaluating those modules. Run fast core mathematical logic manually via specific test files (e.g. `julia --project=. test/test_deb_axes.jl`). Refer to `docs/TESTING_STRATEGY.md` for details.

## Terms to use

- capacity
- pressure
- memory
- adaptive margin
- restoring force
- amplification
- mixture-effect assumptions
- threshold-free features
- vulnerability regimes

## Terms to avoid or use carefully

- attack/defence (reserved for manuscript framing, not code logic)
- interaction (when referring to TU/IA/grouped assumptions, use "mixture-effect assumption")
- damage / $D$ (do not use unless DEBtox scaled damage is actually implemented)
- kappa
- gain
- response_scale
