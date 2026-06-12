# External check #5 — GlobTherm thermal tolerance (capacity coherence probe)

*Working note, 2026-06-13. The "glob first" step from `margin_validation_scouting.md`.
A **pre-registered coherence probe**, not a test of the margin's erosion mechanism
(the margin model does not encode thermal limits — that is the Scope-for-Growth
follow-up). Asks whether the AmP-derived **recovery-capacity axis** carries any
independent thermal-tolerance signal. Reproduce: `scripts/extract_amp_for_globtherm.jl`
→ `examples/globtherm_validation.jl`. n=664 AmP↔GlobTherm species matched.*

## Pre-registration (frozen before computing)
**H (general resilience):** species with greater recovery capacity (`k_M`, the
COMADRE-validated recovery floor; `λ_max`) tolerate **broader** thermal ranges
(`CTmax − CTmin`), beyond body-size and latitude (Janzen) confounds. **Null is
informative:** no signal ⇒ recovery capacity is stressor-*specific*; a clear positive
⇒ it generalises to thermal resilience. Circularity guard: AmP rates carry the
species' thermal regime, so control |latitude| + body size, and re-run on
**ectotherms only** (Class ≠ Mammalia/Aves) additionally controlling `T_typical`
(endotherm `T_typical` is the AmP default body temp ~310.65 K and cannot de-confound).

## Result
Spearman; partials control |latitude|, log body mass (+ `T_typical` for ectotherms).

| | CTmax | CTmin | breadth | breadth (partial) |
| --- | --- | --- | --- | --- |
| **All taxa (n=664)** | | | | |
| `k_M` | −0.136 ** | +0.368 ** | −0.453 ** | **−0.451 ** |
| `λ_max` | −0.010 | +0.327 ** | −0.247 ** | −0.307 ** |
| `g` | +0.139 ** | −0.236 ** | +0.366 ** | +0.381 ** |
| `E_m` (A0) | −0.100 * | +0.251 ** | −0.371 ** | −0.383 ** |
| **Ectotherms only (n=204; partial n=74)** | | | | |
| `k_M` | −0.249 ** | −0.078 | −0.196 | −0.276 * |
| `λ_max` | −0.158 * | −0.224 * | +0.112 | −0.285 * |
| `g` | +0.117 | −0.091 | +0.290 * | +0.129 |
| `E_m` (A0) | +0.038 | +0.100 | −0.133 | −0.021 |

## What it means — two honest, opposite-pointing findings
1. **The AmP capacity axis is externally *coherent*, not noise.** It correlates
   strongly (|ρ| up to 0.45) with an entirely independent physiological dataset
   (measured CTmax/CTmin). The AmP-derived quantities clearly carry real biological
   structure — they track the fast/slow and endotherm/ectotherm axes that also
   organise thermal tolerance. As a "is the AmP extraction sane against outside data"
   check, this passes.
2. **But the pre-registered general-resilience hypothesis is REFUTED.** Recovery
   capacity does **not** predict broader thermal tolerance; the sign is the opposite
   (higher `k_M`/`λ_max` → *narrower* breadth, −0.45 all-taxa, weak −0.28\* in the
   cleaner ectotherm subset; `g`/`E_m` null there). The all-taxa magnitude is mostly
   the endotherm/ectotherm contrast (fast homeotherms are thermal specialists); the
   ectotherm-only, regime-controlled signal is weak and marginal (n=74).

So **recovery capacity and thermal tolerance are separate physiological axes** — the
COMADRE-validated "recovery capacity" is *specific to demographic recovery*, not a
universal resilience currency, and if anything fast/high-capacity species are mild
thermal *specialists*.

## Verdict
GlobTherm does **not** externally validate the adaptive margin — as expected a priori
(the margin's currency is energetic capacity-under-pressure, not thermal limits). Its
value here is twofold and honest: (a) a *coherence* win — AmP capacity is sensibly,
strongly related to independent physiology; and (b) a *bounding* result — do **not**
sell the capacity axis as "general resilience"; it is recovery-specific. This cleanly
motivates the real test: **Scope for Growth**, the margin's own energetic currency,
measured along contamination gradients (next).

## Sources
- GlobTherm — [Bennett et al. 2018, *Scientific Data* 5:180022](https://www.nature.com/articles/sdata201822); data [Zenodo rec. 4976423](https://zenodo.org/records/4976423) (CC-BY), raw `GlobalTherm_upload_02_11_17.csv` (gitignored).
- Janzen (1967) thermal-breadth–latitude hypothesis; thermal specialist/generalist trade-off literature.
