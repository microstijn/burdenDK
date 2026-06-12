# External validation #1 — COMADRE demographic recovery

*Working note, 2026-06-12. The first **external** anchor for the framework (all
prior checks were internal). Tests whether the model's recovery quantities predict
independent demographic recovery from the COMADRE animal matrix database, **beyond
raw pace-of-life** — the generation-time-controlled partial correlation that
`external_anchor_scouting.md` identified as the only clean test. Reproduce:
`scripts/extract_comadre_recovery.jl` (extraction) →
`examples/comadre_partial_validation.jl` (analysis).*

## Method

- **COMADRE recovery** = log₁₀ damping ratio `|λ1|/|λ2|` of each wild, unmanipulated
  projection matrix (rate of return to stable structure), aggregated per species.
- **Control** = Caswell generation time `log R0 / log λ1`, computed from the same
  matrices (so the control is independent of AmP).
- **Model quantities** (from AmP): `λ(A0)` pristine recovery rate, `λ_min = k_M`
  (maintenance rate constant), and `g = λ_max/λ_min` (the amplification axis).
- Species matched by name (AmP key = COMADRE `SpeciesAccepted`). **183 species
  matched; 179 with generation time.** Partial Spearman controlling generation time.

## Result

| model quantity | raw ρ | partial ρ (\| gen. time) |
| --- | --- | --- |
| `λ(A0)` recovery rate | +0.358 ** | **+0.168 \*** |
| `λ_min = k_M` | +0.398 ** | **+0.256 \*\*** |
| `g` (amplification axis) | −0.031 | −0.115 (n.s.) |

(*p<0.05, **p<0.01, n=179/183.) References: ρ(COMADRE recovery, generation time)
= −0.402; ρ(`k_M`, generation time) = −0.557.

## What it means

1. **First external corroboration — and it's positive.** Species the model says
   recover faster (higher `λ(A0)`, higher `k_M`) genuinely have faster demographic
   return-to-equilibrium in COMADRE (raw ρ ≈ 0.36–0.40).
2. **The signal survives the circularity control.** After removing generation time,
   the partial correlations stay positive and significant (`k_M` ρ=0.256**, `λ(A0)`
   ρ=0.168*). So **the recovery/margin machinery predicts demographic recovery
   *beyond* what pace-of-life alone explains** — modest but real.
3. **The `g` amplification axis is *not* corroborated** (partial ρ=−0.115, n.s.) —
   orthogonal to demographic recovery, echoing the life-history check. The external
   data validates the **recovery/margin** half of the model, not the **`F`/`g`
   amplification scalar** — direct support for the margin-first reframe.

## Caveats (honest)

- **Scale bridge.** COMADRE recovery is *population* demographic (damping ratio per
  projection interval); the model's `λ` is *individual-energetic* (per day). Rank
  correlation is unit-invariant, but the individual→population bridge is assumed
  (defensible via DEB, not 1:1).
- **Modest magnitudes.** Partial ρ ≈ 0.17–0.26 — significant at n=179 but a small
  effect; this is corroboration, not strong prediction.
- **Matrix-quality filtering** (wild, unmanipulated, primitive matrices) and
  name-only species matching introduce noise; ~183/286 COMADRE species matched AmP.

## Verdict

The recovery-capacity framing has its **first independent, generation-time-robust
validation signal**; the amplification scalar does not. Next refinements if pursued:
(a) phylogenetic / taxonomic-rank controls, (b) a per-axis margin recovery metric
rather than `λ(A0)`, (c) sensitivity to the matrix-quality filters.

## Sources

- COMADRE database — [Salguero-Gómez et al. 2016, *J. Anim. Ecol.* (data CC-BY)](https://besjournals.onlinelibrary.wiley.com/doi/10.1111/1365-2656.12482) · [compadre-db.org](https://compadre-db.org)
- Demographic resilience framework — [Towards a Comparative Framework of Demographic Resilience (TREE)](https://www.sciencedirect.com/science/article/pii/S0169534720301312)
