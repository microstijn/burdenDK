# COMADRE External Validation (#1)

The first **external** anchor for the framework (all prior checks were internal).
Tests whether the model's recovery quantities predict independent demographic recovery
from the [COMADRE](https://compadre-db.org) animal matrix database, **beyond raw
pace-of-life** and beyond coarse phylogeny.

## Method
- **COMADRE recovery** = log₁₀ damping ratio `|λ₁|/|λ₂|` of each wild, unmanipulated
  projection matrix (rate of return to stable structure), aggregated per species.
- **Control** = Caswell generation time `log R₀ / log λ₁`, from the same matrices
  (independent of AmP).
- **Model quantities** (from AmP): `λ(A0)` pristine recovery rate, `λ_min = k_M`
  (maintenance rate constant), `g = λ_max/λ_min` (amplification axis).
- **Species matching** — harmonised name map (see below): **197 species** matched
  (193 with generation time).
- Three nested controls: raw Spearman; partial controlling generation time; and
  **Order-controlled** (within taxonomic Orders, group-mean-centred, *and* controlling
  generation time — a tree-free phylogenetic proxy; 174 species in 30 Orders).

## Result

| model quantity | raw ρ | \| gen. time | \| gen. time + Order |
| --- | --- | --- | --- |
| `λ(A0)` recovery rate | +0.362 ** | +0.173 * | +0.089 (n.s.) |
| **`λ_min = k_M`** | +0.406 ** | +0.264 ** | **+0.190 \*** |
| `g` (amplification) | −0.109 | −0.128 | −0.055 (n.s.) |

*\* p<0.05, ** p<0.01; n=197.* `k_M` predicts demographic recovery **beyond both
pace-of-life and coarse phylogeny**. `λ(A0)`'s broader signal was largely *between*-clade
(drops to n.s. under Order control). The `g` amplification axis is **not** corroborated
at any level — direct support for the margin-first reframe (the recovery/margin layer
validates; the `F`/`g` readout does not).

## Species name harmonisation
AmP↔COMADRE matching by exact string left 103 of 286 COMADRE species unmatched. A
[GBIF Backbone](https://www.gbif.org)-based resolver
(`scripts/resolve_comadre_amp_names.jl`) harmonises names in four stages:
1. exact match;
2. duplicated-genus-typo fix (`"Hydroprogne Hydroprogne caspia"` → `Hydroprogne caspia`);
3. trinomial → binomial (`"Pan troglodytes schweinfurthii"` → `Pan troglodytes`);
4. GBIF synonym/accepted-name resolution (`Rana catesbeiana` → `Lithobates catesbeianus`).

This recovered **15 species (183 → 197)**, all taxonomically verified (no false
positives). Multiple COMADRE names resolving to one AmP species are collapsed
(quantities averaged) to avoid pseudoreplication. The signal is **stable** under the
larger sample (`k_M` within-Order 0.190\* vs 0.200\* at n=183) — not an artifact of which
species happened to match exactly.

## Caveats
- **Scale bridge.** COMADRE recovery is *population* demographic; the model's `λ` is
  *individual-energetic*. Rank correlation is unit-invariant, but the
  individual→population bridge is argued (DEB-structured models), not 1:1.
- **Modest magnitudes.** Partial ρ ≈ 0.17–0.26 — corroboration, not strong prediction.
- **Specification-sensitive (important).** The gen-controlled `k_M` signal is
  **rank-based**: partial *Spearman* = +0.264, but the log-linear (Pearson) partial is
  only +0.04. Monotone-but-not-log-linear — see [Phylogenetic PGLS](Phylogenetic-PGLS.md).

## Verdict
The recovery-capacity framing has its first independent validation that survives both
pace-of-life and coarse-phylogeny controls, carried specifically by the DEB maintenance
rate constant `k_M`. The amplification scalar does not. Corroboration, not strong
prediction — but the cleanest external evidence to date, pointing at the margin/recovery
layer. Repo notes: `docs/notes/comadre_partial_validation.md`.
