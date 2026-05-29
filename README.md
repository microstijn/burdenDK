# TwoTimescaleResilience / burdenDK

## One-sentence summary

TwoTimescaleResilience is a Julia framework for modelling background-conditioned vulnerability: chronic environmental pressure and retained compound burden narrow species-specific adaptive margin, reduce restoring force, and amplify the burden of later acute perturbations.

## What problem does this solve?

Regulatory threshold maps are useful but insufficient: chronic, cumulative, persistent, mixture-mediated pressure can erode response capacity before discrete endpoint failure. This package models adaptive-margin depletion and amplification rather than only threshold exceedance, using capacity–pressure–memory language to describe vulnerability.

## Conceptual architecture

- **Capacity** = AmP/DEB-informed species response capacity
- **Pressure** = ECOTOX-derived stressor evidence and concentration fields
- **Memory** = retained internal burden $B_t$ through $\rho$ and $K$
- **Response** = $Q_t$, $A_t$, $\lambda(A_t)$, $F_t$
- **Optional spatial/regime workflow** = threshold-free features $\rightarrow$ standardisation $\rightarrow$ clustering $\rightarrow$ raster/NetCDF outputs

**Core chain:**
$C_{j,t} \rightarrow B_{j,t} \rightarrow x_{j,t} \rightarrow E_{\text{axis}} \rightarrow Q_t \rightarrow A_t \rightarrow \lambda_t \rightarrow F_t$

## What this package is

The framework provides tools for multi-scale resilience and toxicokinetic modelling:

- **DEB-axis response math** (`src/deb_axes.jl`, `src/reduced_deb_response.jl`)
- **AmP species library loading and DEBAxisParams conversion** (`src/amp_library.jl`)
- **ECOTOX parsing and toxicity-library runtime loading** (`src/ecotox_library.jl`)
- **ECOTOX active stress and effect-code routing** (`src/ecotox_library.jl`)
- **Compound memory and analytical warm-up** (`src/compound_memory_warmup.jl`)
- **Mixture-effect aggregation** (`src/mixture_aggregation.jl`)
- **Response modes** (`src/mode_of_action.jl`)
- **Threshold-free vulnerability feature vectors** (`src/vulnerability_feature_vectors.jl`)
- **Feature standardisation** (`src/vulnerability_feature_vectors.jl`)
- **Vulnerability-regime clustering** (`src/vulnerability_regime_clustering.jl`)
- **Vulnerability-regime raster/NetCDF output helpers** (`src/vulnerability_regime_outputs.jl`, `src/netcdf.jl`)
- **Examples and diagnostics** (`examples/`)

## What this package is not

- **Not a full DEB implementation**
- **Not DEBkiss**
- **Not full DEBtox**
- Does not currently implement DEB reserve, structure, maturity, kappa allocation, growth/reproduction ODEs, starvation, GUTS survival, or DEBtox scaled damage.
- Does not implement synergism/antagonism/fitted interaction coefficients.
- Does not yet ingest real external raster products natively (unless audited example code explicitly does so, e.g. ISIMIP NetCDF processing).

## Core equations

**Memory:**
$$ B_t = \rho B_{t-1} + (1-\rho) K C_t $$

**Stress:**
$$ x = \max\left(0, \frac{B \text{ or } C - \text{NOEC}}{\text{EC50} - \text{NOEC}}\right) $$

**$E_{\text{axis}}$ mixture models:**
- TU (Toxic Unit sum)
- IA (Independent Action)
- Grouped CA then IA

**Response capacity mapping:**
$$ Q_t = \sum_a w_a E_a $$

$$ A_t = A_0 \max(10^{-6}, 1 - Q_t) \quad \text{for EC50/precomputed impairment path} $$

**Restoring force and amplification:**
$$ \lambda(A) $$
$$ F_t = \frac{\lambda(A_0)}{\lambda(A_t)} $$

**Analytical warm-up:**
- Constant background
- Periodic cycle

## Main source files

- `src/TwoTimescaleResilience.jl`: Main module definition and exports.
- `src/amp_library.jl`: Parses and loads `AmP_Species_Library.json` and creates DEBAxisParams.
- `src/compound_memory_warmup.jl`: Calculates initial burdens and analytical warm-up models.
- `src/deb_axes.jl`: Core structs and definitions for mapping stressors to DEB axes.
- `src/ecotox_library.jl`: Loads, validates, and routes ECOTOX data into empirical physiological burden vectors.
- `src/mixture_aggregation.jl`: Defines mixture effect arithmetic (`TU`, `IA`, `grouped_CA_then_IA`).
- `src/reduced_deb_response.jl`: Calculates adaptive margin, restoring force, and amplification factor.
- `src/vulnerability_feature_vectors.jl`: Threshold-free feature construction and standardisation for spatial model output arrays.
- `src/vulnerability_regime_clustering.jl`: Clusters threshold-free vulnerability features into discrete regimes.
- `src/vulnerability_regime_outputs.jl`: Provides functions to bundle spatial vulnerability features, clusters, and outputs.

## Data files

- `data/AmP_Species_Library.json`: Source library of AmP species parameters consumed by `amp_library.jl`.
- `data/ECOTOX_Toxicity_Library.json`: Source library of parsed ECOTOX empirical data consumed by `ecotox_library.jl`.
- `data/Compound_Memory_Library.csv`: Source library defining retention and bioaccumulation properties for compound memory consumed by `ecotox_library.jl`.
- `data/AmP_Species_Archetypes.csv` and `data/AmP_Species_Archetypes.json`: Derived databases mapping species onto generalized archetypes.

## Examples

- **Quick:**
  - `examples/synthetic_raster_demo.jl`: Fast demonstration of synthetic background generation and spatial metrics plotting.
  - `examples/ecotox_amp_multiaxis_response_calibrated_demo.jl`: Fast demonstration of ECOTOX and AmP integration via a multi-axis response.
- **Medium:**
  - `examples/ecotox_amp_multispecies_multicompound_demo.jl`: Runs ECOTOX metrics across multiple AmP species simultaneously.
  - `examples/ecotox_amp_multispecies_multicompound_monthly_memory_demo.jl`: Tests memory with time-series data for multi-species profiles.
- **Heavy / Extended:**
  - `examples/isimip_moa_deb_3x3_demo.jl`: Demonstrates ISIMIP spatial data manipulation and model outputs.
  - `examples/nc_real_raster_deb_axes_demo.jl`: Parses, calculates, and exports real-world vulnerability metrics to NetCDF.

*(Note: Do not run heavy examples in every test iteration; use environment variables or local execution for these.)*

## Tests

The testing framework evaluates units across fast logic and heavier empirical/output operations. Due to the high computational load and dependency precompilation, the test suite requires separating components.

**Current vs Proposed Strategy:**
Currently, the `test/runtests.jl` manually enables or disables subsets of tests. The proposed future split should segregate fast core logic, extended integration testing, example outputs, and plotting/GIS tools using test suites or specific environment gating.

**Environment variables used:**
- `TTR_RUN_EXTENDED_TESTS`: Used to optionally gate heavier integration tests.
- `TTR_RUN_EXAMPLE_TESTS`: Enables example-output tests.
- `TTR_MAKE_EXAMPLE_PLOTS`: Ensures plotting outputs are generated when executed.

## Quick start

```julia
using TwoTimescaleResilience

# 1. Load data
amp_lib = load_amp_species_library()
ecotox_lib = load_ecotox_library()

# 2. Extract specific AmP species parameters
params = amp_species_deb_params(amp_lib, "Daphnia magna")

# 3. Filter ECOTOX records to relevant taxa
records = ecotox_filter_records(ecotox_lib; taxon_class="Branchiopoda")

# 4. Compute stateless burden and response
concentrations = Dict("7647-14-5" => 2.5) # NaCl diagnostic concentration
burden = ecotox_records_to_deb_burden(concentrations, records)
response = ecotox_burden_to_response(burden, params)

println("Adaptive Margin: ", response.A)
println("Amplification Factor: ", response.amplification)

# 5. Stateful memory calculation
state = EcotoxExposureState()
stateful_burden = ecotox_records_to_deb_burden_stateful!(state, concentrations, records)
stateful_response = ecotox_burden_to_response(stateful_burden, params)
```

## Development notes

- Prefer adding adapters, tests, and documentation rather than rewriting existing architectural core math.
- Keep $B_t$ (chemical memory) explicitly distinct from future $Z_t$ (physiological condition) and DEBtox $D_t$ (scaled damage).
- Use `grouped_ca_then_ia_axis_effects` as the preferred default when grouping effects.
- Maintain strict "threshold-free" naming and behavior constraints within the clustering and standardisation models.
- **Do not introduce** mathematical tuning parameters named or behaving like `kappa`, `κ`, `gain`, `response_scale`, or `burden_to_margin_multiplier`.
