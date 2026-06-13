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
| [How it works — the pipeline](Pipeline.md) | End-to-end data flow from concentrations to vulnerability regimes, stage by stage, with links to component docs. |
| [Model equations](Equations.md) | All the math in one place: memory, stress, mixture models, margin, restoring force, amplification, the two-timescale ODE. |
| [Data & parameters](Data-and-Parameters.md) | AmP / ECOTOX / compound-memory data, the offline AmP→capacity mapping, and the proxy assumptions (ρ, K, effect codes). |
| [Getting started](Getting-Started.md) | Julia version, install/instantiate, quickstart, running demos, testing. |
| [Limitations & open questions](Limitations-and-Open-Questions.md) | Honest status — the κ-collapse, what is implemented vs deferred, where outputs lean on thin evidence. |
| [Water-quality coupling & pMoA routing](Water-Quality-Coupling.md) | Connecting to monthly water-quality models (DynQual-type aggregate stressors): the knob-free pMoA stressor-routing table, the `heavy_metals` field-vs-fitted dual route, and the honest scope of the coupling. |

### External validation (vs COMADRE)

The first **external** validation of the framework (2026-06-12) — does the model
predict independent demographic data from the COMADRE matrix database?

| Page | What it covers |
| --- | --- |
| [COMADRE external validation](COMADRE-External-Validation.md) | The scalar result: the DEB maintenance rate `k_M` predicts demographic recovery beyond pace-of-life + coarse phylogeny; GBIF species-name harmonisation. |
| [Phylogenetic PGLS](Phylogenetic-PGLS.md) (Idea A) | Real-phylogeny PGLS. A **dated TimeTree** (182 spp, 2026-06-13) **confirms** the `k_M`↔recovery rank signal survives pace + real phylogeny (β\*=0.22, p=0.011); the log-linear form is weak and Pagel's λ≈0 (phylogeny was never the confound). |
| [Per-axis resilience](Per-Axis-Resilience.md) (Idea B) | The multi-dimensional test: the DEB reproduction rate `R_i` specifically predicts demographic *compensation* (ρ=0.77) — the strongest external result. |
| [Margin validation](Margin-Validation.md) (SFG + SoS + GlobTherm) | The **margin state** (Scope for Growth, ρ=+0.41 estuary-scale, *no scale bridge*) **and its (static) function** (Stress-on-Stress survival, ρ=+0.39→+0.45 under confound control — the closest cross-sectional shadow of the two-timescale **amplification** claim; the *dynamics* themselves stay untested). GlobTherm bounds the capacity axis (recovery-specific, not general resilience). |
| [Reproducibility](Reproducibility.md) | Exact commands + data provenance for the validation pipeline (extractors, resolvers, analyses). |

**Bottom line:** external support lands on the **recovery/margin layer** — its *rate
endpoints* (`k_M`, `R_i` vs COMADRE), its *state* (Scope for Growth), its *function* (Stress-
on-Stress), and a *first sign of its dynamics* (transplant); the amplification scalar `g`/`F`
predicts nothing — consistent with the margin-first reframe.

> **One consolidated, manuscript-ready account of every external validation (with the honest
> scope of each):** [External validation synthesis](../notes/external_validation_synthesis.md).

### Component reference (deeper dives)

These pre-existing topic docs are the authoritative detail for each component:

- [Compound memory](../compound_memory.md) · [Mixture-effect models](../mixture_effect_models.md)
- [Species archetypes](../species_archetypes.md)
- [Vulnerability feature vectors](../vulnerability_feature_vectors.md) · [Regime outputs](../vulnerability_regime_outputs.md) · [Tranche comparison](../vulnerability_tranche_comparison.md)
- [Architecture graph](../ARCHITECTURE_GRAPH.md) · [Package capabilities](../PACKAGE_CAPABILITIES.md) · [Testing strategy](../TESTING_STRATEGY.md)

### Project / working notes

- [Source audit (2026-06-11)](../claude/TwoTimescaleResilience_source_audit_2026-06-11.md) — current state vs documentation, and the κ-collapse analysis.
- [Review & research agenda (2026-06-11)](../claude/TwoTimescaleResilience_review_and_agenda_2026-06-11.md).

---

## Status at a glance

| Capability | Status |
| --- | --- |
| DEB-axis response math (margin → restoring force → amplification) | implemented, tested |
| AmP species adapter + offline AmP→capacity mapping | implemented, tested |
| ECOTOX runtime + compound memory `B_t` | implemented, tested |
| Mixture-effect models (TU, IA, grouped CA-then-IA) | implemented, tested |
| Threshold-free spatial features → clustering → NetCDF outputs | implemented, tested |
| Physiological condition memory `Z_t` | **implemented, opt-in, off by default** (not yet validated) |
| Stable real-raster ingestion | partial (mostly example scripts) |
| DEBtox scaled damage `D_t`, synergism/antagonism | not implemented (by design) |
| External validation vs COMADRE | **recovery/margin layer corroborated** (`k_M`, `R_i`); amplification scalar null |
| External validation vs Scope for Growth | **margin state corroborated** where burden indexes exposure (ρ=+0.41; scale-attenuated, confound-bounded) — see [Margin validation](Margin-Validation.md) |
| External validation vs Stress-on-Stress | **static margin↔acute-resilience map corroborated** — modelled margin → acute-stress survival (ρ=+0.39→+0.45 confound-controlled, n=17); the *dynamics* (accumulate→erode→amplify) stay **untested** — see [Margin validation](Margin-Validation.md) |
| Single-trait `k_M`→toxicity (bounding) | **body-size-confounded** — maintenance→sensitivity replicates raw (n=310) but nulls under a size control; the distinctive content is the across-axis capacity weighting, not `k_M` alone — see [Margin validation §5](Margin-Validation.md) |
| Water-quality coupling (pMoA routing) | **engine validated; input routing is a declared pMoA assumption** — knob-free stressor→axis table; `heavy_metals` carries a field-vs-fitted dual route (metals→assim degrades all field anchors) — see [Water-quality coupling](Water-Quality-Coupling.md) |

See [Limitations & open questions](Limitations-and-Open-Questions.md) for the important caveats — in particular, the amplification factor is a **one-dimensional index** (it was the allocation fraction κ; after re-anchoring the recovery floor to the DEB maintenance rate constant it now tracks the energy investment ratio `g`). As of 2026-06-12 the framework has its **first external validation** against COMADRE (see [External validation](#external-validation-vs-comadre) above): the **recovery/margin layer is corroborated** (the DEB rates `k_M`↔recovery and `R_i`↔compensation), while the **amplification scalar `g`/`F` remains null** — direct support for the margin-first reframe.

---

## Maintaining this wiki

This wiki lives **in the repository** (`docs/wiki/`), so it is versioned and
reviewed in the same pull requests as the code it documents. When you change
behaviour, update the relevant page in the same PR. Conventions:

- One concept per page; link rather than duplicate.
- Keep [Limitations & open questions](Limitations-and-Open-Questions.md) honest and current — it is the page that protects the model's credibility.
- Figures live in [`docs/wiki/figures/`](figures/) and are regenerated by
  `examples/wiki_figures.jl` (`julia +release --project=. examples/wiki_figures.jl`) — edit the script, don't hand-edit images.
- If you want a literal GitHub **Wiki tab**, mirror this folder into the
  `*.wiki.git` repo; the in-repo copy remains the source of truth.
