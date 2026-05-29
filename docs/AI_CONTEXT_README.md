# AI Context for TwoTimescaleResilience

## Do first

- Inspect actual code before assuming functionality (read files directly rather than relying on prompts or standard DEB assumptions).
- Run targeted tests during development, *not* the full heavy suite, to prevent timeouts from precompilation.
- Preserve the current architecture and invariants.
- Do not rewrite core math unless specifically requested.

## Project identity

TwoTimescaleResilience is a framework for mapping background chronic stress (like water quality) onto an organism's adaptive margin and restoring force to calculate its amplification of vulnerability to acute stressors.

## Current architecture

**Capacity--pressure--memory.**
The project evaluates multi-scale vulnerabilities: capacity provided by DEB/AmP capabilities, pressure exerted by ECOTOX empirical data arrays, and memory generated through analytical retention/bioaccumulation factors ($B_t$).

## Core files and what they do

- `src/deb_axes.jl` and `src/reduced_deb_response.jl`: Defines DEB axis structures and capacity margin logic.
- `src/amp_library.jl`: Parses and standardises the AmP library to produce DEB axis parameters.
- `src/ecotox_library.jl`: Filters ECOTOX data and converts it into DEB burden arrays, incorporating compound memory logic.
- `src/compound_memory_warmup.jl`: Fast mathematical solutions for initializing internal burdens without simulating spin-up.
- `src/mixture_aggregation.jl`: Combines specific pressures via mechanistic assumptions (IA, TU).
- `src/vulnerability_feature_vectors.jl`: Defines threshold-free continuous outputs for raster cells.
- `src/vulnerability_regime_clustering.jl`: Classifies those feature vectors.

## What not to change

- Do not add or combine interaction components: Synergism and antagonism are explicitly not implemented. Mixture effects must remain mathematically conservative (e.g. IA or TU).
- Do not add arbitrary physiological modifiers or scalar tuning variables: `kappa`, `κ`, `gain`, `response_scale`, `burden_to_margin_multiplier` are entirely forbidden.
- Do not conflate memory variables. $B_t$ is internal compound concentration and explicitly separate from unbuilt DEBtox scaled damage ($D_t$) or general physiological recovery condition ($Z_t$).
- Do not implement explicit threshold constraints (`_gt_`, `threshold`) within standard vulnerability output vectors.

## Common next tasks

- Integrating formal DEBtox $D_t$ scaled damage parameters. *(future)*
- Integrating formalized $Z_t$ physiological memory capabilities to model slow baseline health recoveries independent of chemical clearance. *(future)*
- More formalized tests integrating external GIS processes into automated suites. *(future)*

## Test strategy

The full suite is currently heavy as it invokes plotting packages, empirical test datasets, GIS systems, and NetCDF processing. Running a simple `Pkg.test()` will often time out or take a significant amount of time due to heavy precompilation of CairoMakie or NCDatasets.

**For tests:**
Use targeted test commands. E.g., `julia --project=. test/runtests.jl` with specific components disabled, or test individual files directly.

## Terms to use

- capacity
- pressure
- memory
- adaptive margin
- restoring force
- amplification
- mixture-effect assumptions
- threshold-free features

## Terms to avoid or use carefully

- "attack/defence" (only used previously in manuscript context)
- "interaction" when referring to TU/IA/grouped assumptions (these are not mathematically fitted interactions, they are mixture assumptions)
- "damage D" unless explicit DEBtox scaled damage is formally implemented.
- "kappa"
- "gain"
- "response_scale"