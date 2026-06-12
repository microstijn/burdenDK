# Getting started

[← Data & parameters](Data-and-Parameters.md) · next: [Limitations →](Limitations-and-Open-Questions.md)

## Julia version (read first)

**Use Julia 1.12.6.** The committed `Manifest.toml` is resolved under
`julia_version = "1.12.6"`. A machine default of the LTS (1.10.x) **cannot load
this project** — it fails while precompiling `PrecompileTools` with
`UndefVarError: StaticData not defined`, which is a version mismatch, not a code
bug. Run everything through the `release` channel:

```powershell
julia +release --project=. <script.jl>
```

`juliaup status` lists installed channels; `release` is 1.12.6 here. If you hit
the `StaticData` error, you are on the wrong Julia — switch to `+release`. See
also [`CLAUDE.md`](../../CLAUDE.md).

## Install / instantiate

```powershell
julia +release --project=. -e "using Pkg; Pkg.instantiate()"
```

## Quickstart

```julia
using TwoTimescaleResilience

# 1. Capacity: load AmP-derived species parameters
amp_lib = load_amp_species_library()
params  = amp_species_deb_params(amp_lib, "Daphnia magna")

# 2. Pressure: load ECOTOX records and filter to relevant taxa
ecotox_lib = load_ecotox_library()
records    = ecotox_filter_records(ecotox_lib; taxon_class = "Branchiopoda")

# 3. Stateless burden -> response
concentrations = Dict("7647-14-5" => 2.5)              # CAS => concentration
burden   = ecotox_records_to_deb_burden(concentrations, records)
response = ecotox_burden_to_response(burden, params)
println("Adaptive margin A_t   : ", response.A)
println("Amplification F_t     : ", response.amplification)

# 4. Stateful memory (carry burden across time steps)
state = EcotoxExposureState()
b2 = ecotox_records_to_deb_burden_stateful!(state, concentrations, records)
r2 = ecotox_burden_to_response(b2, params)
```

The point-level response API is `compute_adaptive_margin_response`; it now
defaults to the nondimensional `ec50_anchored_fractional_impairment` mode (see
[Pipeline §5](Pipeline.md#stage-5--response-e_axis--q--a---f)).

## Running the demos

Examples live in `examples/`. Lightweight ones first:

```powershell
julia +release --project=. examples/ecotox_amp_multiaxis_response_calibrated_demo.jl
julia +release --project=. examples/amp_kappa_collapse_diagnostic.jl     # read-only analysis
```

Heavier raster/plotting demos pull `CairoMakie` / `GeoMakie` / `NCDatasets` and
precompile slowly; run them deliberately.

## Regenerating the wiki figures

```powershell
julia +release --project=. examples/wiki_figures.jl   # writes docs/wiki/figures/*.png
```

## Testing

```powershell
# Fast core suite (default)
julia +release --project=. test/runtests.jl
```

Heavy / plotting / example / DynQual tests are gated behind environment flags so
the fast suite stays quick:

| Flag | Enables |
| --- | --- |
| `TTR_RUN_EXTENDED_TESTS=true` | extended grid/demo tests |
| `TTR_RUN_EXAMPLE_TESTS=true` | example-output tests |
| `TTR_RUN_PLOTTING_TESTS=true` | plotting tests |
| `TTR_RUN_DYNQUAL_TESTS=true` | DynQual helper tests |

See [Testing strategy](../TESTING_STRATEGY.md).

## Project conventions

These invariants are enforced by review and partly by code; don't violate them
without explicit sign-off (full list in [`CLAUDE.md`](../../CLAUDE.md)):

- **No arbitrary tuning knobs** (`gain`, `response_scale`, `burden_to_margin_*`,
  or κ used as a free knob). *(The `KA = 0.3·A0` knob has been removed — the recovery
  curve is now linear; see [Limitations §2](Limitations-and-Open-Questions.md).)*
- **No thresholds** in threshold-free spatial features (`_gt_`, `exceedance`, …).
- **Keep memory layers distinct:** `B_t` (chemical), `Z_t` (condition, opt-in),
  `D_t` (DEBtox, unimplemented).
- **Mixtures are assumptions, not fitted interactions** (TU / IA / grouped
  CA-then-IA only).
- Don't rewrite core math unless explicitly asked.
