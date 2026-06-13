# TwoTimescaleResilience / burdenDK

A Julia framework for **background-conditioned vulnerability**: chronic
environmental pressure and retained chemical burden slowly narrow a species'
*adaptive margin*, which weakens its *restoring force*, which *amplifies* the
burden of a later acute perturbation. It models capacity erosion and
amplification rather than threshold exceedance, using a **capacity – pressure –
memory** architecture (AmP/DEB capacity · ECOTOX pressure · compound-retention
memory) with continuous, **threshold-free** outputs.

![Adaptive margin erodes before any threshold is crossed](docs/wiki/figures/framework_margin_map.png)

> *The framework's output is the **adaptive margin** `A_t` (centre): chronic background
> pressure (left) erodes it across a broad area that mostly never crosses a regulatory limit
> (right). **Exceedance ≠ vulnerability** — the margin is depleted before any threshold is breached.
> Regenerate: `examples/framework_margin_map.jl`.*

## 📖 Documentation — start with the wiki

The full documentation lives in **[`docs/wiki/`](docs/wiki/Home.md)** (in-repo, so
it is versioned and reviewed with the code):

| | |
| --- | --- |
| **[Overview](docs/wiki/Overview.md)** | the idea — two timescales, capacity–pressure–memory, what it is / is not |
| **[How it works](docs/wiki/Pipeline.md)** | end-to-end pipeline, stage by stage |
| **[Model equations](docs/wiki/Equations.md)** | all the math in one place |
| **[Data & parameters](docs/wiki/Data-and-Parameters.md)** | AmP / ECOTOX / memory data + the offline mapping |
| **[Getting started](docs/wiki/Getting-Started.md)** | install, quickstart, demos, testing |
| **[Life stages & movement](docs/wiki/Life-Stages-and-Movement.md)** | stage-resolved capacity, occupancy-weighted movement, surface:volume TK — with the salmon worked example |
| **[External validation](docs/wiki/External-Validation.md)** | every external test in one place — scorecard, forest plot, the margin-first result |
| **[Limitations & open questions](docs/wiki/Limitations-and-Open-Questions.md)** | honest status — read this before trusting outputs |

## Quick start

> **Julia 1.12.6 required.** The default LTS (1.10.x) cannot load this project.
> Use the `release` channel — see [Getting started](docs/wiki/Getting-Started.md).

```powershell
julia +release --project=. -e "using Pkg; Pkg.instantiate()"
```

```julia
using TwoTimescaleResilience

amp_lib = load_amp_species_library()
params  = amp_species_deb_params(amp_lib, "Daphnia magna")

ecotox_lib = load_ecotox_library()
records    = ecotox_filter_records(ecotox_lib; taxon_class = "Branchiopoda")

burden   = ecotox_records_to_deb_burden(Dict("7647-14-5" => 2.5), records)
response = ecotox_burden_to_response(burden, params)

println("Adaptive margin A_t : ", response.A)
println("Amplification  F_t  : ", response.amplification)
```

## The core chain

$$ C_{j,t} \;\rightarrow\; B_{j,t} \;\rightarrow\; x_{j,t} \;\rightarrow\; E_{\text{axis}} \;\rightarrow\; Q_t \;\rightarrow\; \boxed{A_t} \;\rightarrow\; \lambda(A_t) \;\rightarrow\; F_t $$

concentration → memory → active stress → per-axis impairment → scalar load →
**adaptive margin (the product)** → restoring force → amplification. The
**adaptive-margin state** — relative depletion, the capacity-aware absolute margin,
and the axis composition — is the vulnerability signal; the amplification factor
`F` is a convenient *derived* scalar that falls out of the two-timescale algebra.
Full detail in [How it works](docs/wiki/Pipeline.md).

## Worked example — life stages & movement

![Atlantic salmon migration](docs/wiki/figures/salmon_migration.png)

The modelling unit can be a **life stage** and the target can **move**. An anadromous salmon carries
contaminant burden *between* regions (memory), recovers far faster as a small parr than as a
returning adult (stage-resolved capacity), and tracks the ambient water more tightly when small
(surface:volume toxicokinetics) — and at spawning it is *compliant yet eroded*. See
[Life stages & movement](docs/wiki/Life-Stages-and-Movement.md) for this and the GeoMakie migration
map.

## Status

| Capability | Status |
| --- | --- |
| DEB-axis response math, AmP adapter, ECOTOX runtime, compound memory | implemented, tested |
| Mixture models (TU / IA / grouped CA-then-IA) | implemented, tested |
| Threshold-free spatial features → clustering → NetCDF | implemented, tested |
| Stage-resolved capacity · mobile-target exposure · surface:volume TK | implemented, tested — **not yet validated** |
| Physiological condition memory `Z_t` | implemented, **opt-in, off by default** |
| Real-raster ingestion | partial (mostly examples) |
| External validation | recovery/margin layer corroborated; amplification `F` **null**; capacity weighting **untested** |
| DEBtox `D_t`, synergism/antagonism | not implemented (by design) |

> ⚠️ **Important caveat:** the amplification factor `F` is a **one-dimensional
> index** — by construction it carries one physiological number per species (now the
> energy investment ratio `g`). External validation found it carries **no external
> signal**, while the **recovery/margin layer is corroborated** — so the framework is
> read **margin-first**. See [External validation](docs/wiki/External-Validation.md)
> and [Limitations](docs/wiki/Limitations-and-Open-Questions.md).

## Not a full DEB model

This is a DEB-*informed*, physiologically-structured **vulnerability index**, not a
DEB / DEBkiss / DEBtox implementation: it borrows the process axes and parameters
but drops the dynamical reserve/structure/maturity/κ state equations. Mixtures are
combined by explicit assumptions, never fitted interaction coefficients.

## Contributing / maintaining

- Update the relevant [`docs/wiki/`](docs/wiki/Home.md) page in the same PR as any
  behaviour change; keep [Limitations](docs/wiki/Limitations-and-Open-Questions.md)
  current.
- Figures are regenerated by their scripts (`examples/framework_margin_map.jl` (landing),
  `wiki_figures.jl`, `validation_forest_plot.jl`, `salmon_migration_figure.jl`,
  `salmon_migration_map.jl`) — edit the script, not the images.
- Project invariants (no arbitrary knobs, no thresholds in spatial features, keep
  `B_t`/`Z_t`/`D_t` distinct, mixtures-as-assumptions) are listed in
  [`CLAUDE.md`](CLAUDE.md) and [Getting started](docs/wiki/Getting-Started.md).
