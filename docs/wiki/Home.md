# TwoTimescaleResilience (burdenDK) — Wiki

A Julia framework for modelling **background-conditioned vulnerability**: chronic
environmental pressure and retained chemical burden slowly narrow a species'
*adaptive margin*, which weakens its *restoring force*, which *amplifies* the
burden of a later acute perturbation.

It is built around a **capacity – pressure – memory** architecture, with capacity
derived from Add-my-Pet (AmP / DEB) parameters, pressure from EPA ECOTOX toxicity
data, and memory from a compound-retention recurrence. Outputs are continuous and
**threshold-free** — no safe/unsafe lines, no exceedance counts.

> **New here?** Read [Overview](Overview.md) → [How it works (pipeline)](Pipeline.md) → [Getting started](Getting-Started.md).

---

## Map of the wiki

| Page | What it covers |
| --- | --- |
| [Overview](Overview.md) | The idea: two timescales, capacity–pressure–memory, the Canguilhem framing, what the framework *is* and *is not*. |
| [How it works — the pipeline](Pipeline.md) | End-to-end data flow from concentrations to vulnerability regimes, stage by stage. |
| [Model equations](Equations.md) | All the math in one place: memory, stress, mixture models, margin, restoring force, amplification, the two-timescale ODE. |
| [Data & parameters](Data-and-Parameters.md) | AmP / ECOTOX / compound-memory data, the offline AmP→capacity mapping, and the proxy assumptions (ρ, K, effect codes). |
| [Getting started](Getting-Started.md) | Julia version, install/instantiate, quickstart, running demos, testing. |
| [Life stages & movement](Life-Stages-and-Movement.md) | Stage-resolved capacity (`λ_max=v_eff/L`, acceleration, maturation↔reproduction), occupancy-weighted movement, and surface:volume aquatic toxicokinetics — with the salmon worked example. |
| [Water-quality coupling & pMoA routing](Water-Quality-Coupling.md) | Connecting to monthly water-quality models (DynQual-type aggregate stressors): the knob-free pMoA stressor-routing table and the honest scope of the coupling. |
| [**External validation**](External-Validation.md) | Every external test in one place: the recovery/margin layer is corroborated (`k_M`, `R_i`, SFG, stress-on-stress, transplant), the amplification scalar is null, the capacity weighting is untested. Scorecard + forest plot. |
| [Limitations & open questions](Limitations-and-Open-Questions.md) | Honest status — read this before trusting outputs. |
| [Reproducibility](Reproducibility.md) | Exact commands + data provenance for the validation pipeline. |

**Validation bottom line:** external support lands on the **recovery/margin layer** — its *rate
endpoints* (`k_M`, `R_i` vs COMADRE), *state* (Scope for Growth), *function* (stress-on-stress), and a
*first sign of its dynamics* (transplant); the amplification scalar `g`/`F` predicts nothing —
consistent with the margin-first reframe. Full account:
[External validation](External-Validation.md) and the manuscript-ready
[synthesis note](../notes/external_validation_synthesis.md).

### Component reference (deeper dives)

- [Compound memory](../compound_memory.md) · [Mixture-effect models](../mixture_effect_models.md) · [Species archetypes](../species_archetypes.md)
- [Vulnerability feature vectors](../vulnerability_feature_vectors.md) · [Regime outputs](../vulnerability_regime_outputs.md) · [Tranche comparison](../vulnerability_tranche_comparison.md)
- [Architecture graph](../ARCHITECTURE_GRAPH.md) · [Package capabilities](../PACKAGE_CAPABILITIES.md) · [Testing strategy](../TESTING_STRATEGY.md)

---

## Status at a glance

| Capability | Status |
| --- | --- |
| DEB-axis response math (margin → restoring force → amplification) | implemented, tested |
| AmP species adapter + offline AmP→capacity mapping | implemented, tested |
| ECOTOX runtime + compound memory `B_t` | implemented, tested |
| Mixture-effect models (TU, IA, grouped CA-then-IA) | implemented, tested |
| Threshold-free spatial features → clustering → NetCDF outputs | implemented, tested |
| **Stage-resolved capacity, mobile-target exposure, surface:volume TK** | **implemented, tested — not yet externally validated** ([Life stages & movement](Life-Stages-and-Movement.md)) |
| Physiological condition memory `Z_t` | implemented, **opt-in, off by default** (not yet validated) |
| Stable real-raster ingestion | partial (mostly example scripts) |
| DEBtox scaled damage `D_t`, synergism/antagonism | not implemented (by design) |
| **External validation** | **recovery/margin layer corroborated**; amplification scalar `g`/`F` **null everywhere**; capacity weighting **untested** ([External validation](External-Validation.md)) |

See [Limitations & open questions](Limitations-and-Open-Questions.md) for the important caveats — in
particular, the amplification factor `F` is a **one-dimensional index** (it tracks the energy
investment ratio `g`), and external validation found it carries **no external signal**, while the
recovery/margin layer is corroborated — direct support for the margin-first reframe.

---

## Maintaining this wiki

This wiki lives **in the repository** (`docs/wiki/`), so it is versioned and reviewed in the same
pull requests as the code. When you change behaviour, update the relevant page in the same PR.
Conventions:

- One concept per page; link rather than duplicate.
- Keep [Limitations & open questions](Limitations-and-Open-Questions.md) honest and current — it is
  the page that protects the model's credibility.
- Figures live in [`docs/wiki/figures/`](figures/) and are regenerated by their scripts — edit the
  script, not the image: `examples/wiki_figures.jl` (mechanism), `examples/validation_forest_plot.jl`
  (forest plot), `examples/salmon_migration_figure.jl` and `examples/salmon_migration_map.jl`
  (life-stage + movement; the map is also the README landing figure).
- If you want a literal GitHub **Wiki tab**, mirror this folder into the `*.wiki.git` repo; the
  in-repo copy remains the source of truth.
