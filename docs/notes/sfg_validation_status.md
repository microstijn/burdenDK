# Scope for Growth — the real margin test: RESULTS + design

## ✅ RESULTS (2026-06-13) — the adaptive margin tracks measured SFG
Data blocker resolved (the Widdows et al. 1995 SFG and contaminant tables were
provided manually and committed: `sfg_widdows1995_northsea.csv`,
`sfg_widdows1995_contaminants.csv`). Test run via `examples/sfg_margin_validation.jl`,
36 North Sea *Mytilus edulis* sites.

**Empirical (replicates Widdows):** SFG declines with hydrocarbon body burden —
`SFG ~ total_toxic_HC` Spearman ρ = **−0.467 \*\*** (PAH −0.467\*\*, THC −0.427\*\*); metals
mostly weak (Cu −0.33\*), organotins/organochlorines n.s. The gradient is real and
hydrocarbon-driven, exactly Widdows' finding.

**Model:** the **modeled adaptive margin tracks the measured energetic margin**,
`SFG ~ A_t` ρ = **+0.405 \*** (n=36). Comparison:

| predictor of SFG | \|ρ\| |
| --- | --- |
| best single contaminant (total hydrocarbons) | 0.467 |
| **model adaptive margin `A_t`** (MoA-routed) | **0.405** |
| naive equal-weight mean toxic-unit | 0.221 |

**Reading (honest):**
- **First external evidence for the adaptive margin *itself*** — not its rate endpoints
  (COMADRE) or a tangential capacity facet (GlobTherm). The modeled margin tracks an
  *independent, same-level* (individual-energetic) measurement along a *real pressure
  gradient*, with **no scale bridge**. The erosion mechanism (pressure → margin ↓ →
  energetic outcome ↓) is corroborated.
- **The model's structure earns its keep over naive aggregation** (0.41 vs 0.22): the
  mode-of-action routing sends hydrocarbons to the assimilation/feeding axis rather than
  diluting them with background metals — recovering the true signal.
- **It does not beat the single best contaminant (0.47)** — and shouldn't be expected
  to: with **one species** the AmP capacity weighting is constant across sites, so the
  model can at best *recover* the dominant-hydrocarbon signal, not exceed it. Testing the
  capacity dimension (the model's distinctive content) needs **across-species** SFG.

**Caveats (carried):** single-species → capacity weighting untested; tissue body burden
vs exposure-based potency (critical-body-residue mismatch) → a threshold-free
median-normalised relative-burden pressure was used, so per-contaminant potency is not
encoded; the contaminant→axis routing is a documented, approximate pMoA assignment.

**Verdict:** the adaptive margin has its first direct, same-level external corroboration.
Modest (ρ≈0.4) and single-species, but it validates the *thing we want to use* — the
margin state under pressure — not just adjacent rates. Next: an across-species SFG set
(or a second gradient) to test the capacity weighting the single-species design holds fixed.

---

# Design, feasibility (retained for reference)

*Working note, 2026-06-13. The "then SFG" step from `margin_validation_scouting.md`.
SFG is the **right** external anchor for the adaptive margin: SFG = energy absorbed −
energy respired = capacity beyond maintenance ≈ the margin in energetic terms, measured
at the **same organisational level** (individual energetics → no scale bridge) and
**independent** of AmP. This note pins the dataset, gives the test design, and records
the one blocker that stops a fully-autonomous run.*

## Why SFG is the test the margin actually needs
COMADRE validated the margin curve's *rate endpoints*; GlobTherm showed the capacity
axis is recovery-*specific* (not general resilience). Neither tested the margin's core
claim: **chronic pressure erodes the margin → an energetic outcome declines.** SFG along
a contamination gradient tests exactly that, with three advantages over COMADRE:
same level (no scale bridge), an independent physiological measurement, and a genuine
pressure gradient.

## Feasibility — the model side is READY
- **AmP coverage:** `Mytilus_edulis`, `M. galloprovincialis`, `M. californianus`,
  `M. trossulus` are all in AmP.
- **Pressure data:** the full EPA **ECOTOX** database is in-repo
  (`data/ecotox/…`, results/tests/species) — *Mytilus* and the relevant contaminants
  (PAHs, PCBs, metals, TBT) are extensively tested, so per-contaminant anchors are
  derivable.
- **Margin machinery:** the margin-first point API exists and is the default —
  `compute_adaptive_margin_response(...; response_mode="ec50_anchored_fractional_impairment")`
  → `A_t = A0·(1−Q)` and `λ(A_t)`.

So the model can produce a per-site margin prediction. The remaining engineering is
wiring per-site contaminant load → axis pressures → margin at the point level (a subset
of the de-scoped spatial-pipeline migration, but tractable point-wise).

## The one blocker — the empirical anchor is access-gated here
The per-site **SFG + contaminant** tables live in publisher PDFs/tables that this
environment cannot fetch (int-res.com returns 401 to automated access; ScienceDirect is
paywalled). This is the "literature-assembly step" flagged in the scouting note. It is a
~30-minute manual extraction with institutional access, not a research problem.

**Pinned candidate datasets (best first):**
1. **Widdows et al. 1995, *MEPS* 127:131–148** — North Sea / Langesundfjord transect,
   *Mytilus edulis* (in AmP), per-site SFG + PAH/PCB/metal tissue loads. The canonical
   gradient study. *(Open-access abstract; PDF access-gated here.)*
2. **Beiras et al. 2003 (Spanish coast, ~41 stations)** — large-scale SFG + contaminants
   in wild *M. galloprovincialis*. Best n if obtainable.
3. **Widdows & Page / Hamilton Harbour transplant** — *Mytilus*, TBT/PAH/PCB gradient.

A turnkey extraction template is committed at
`data/external/sfg_gradient_TEMPLATE.csv` (fill one row per site).

## Test design (pre-registered intent)
For species s at site i with measured contaminant vector `c_i`:
1. ECOTOX-anchor each contaminant → axis pressures (mode-of-action routing via
   `moa_deb_mapping`); mixture-aggregate (TU/IA, no fitted interactions).
2. `compute_adaptive_margin_response` → modeled margin `A_t,i` (and `λ(A_t,i)`).
3. **Primary test:** rank-correlate modeled margin across sites vs measured SFG.
4. **Model-specific test (the one that matters):** does the capacity-/axis-weighted
   margin predict SFG **better than raw summed contaminant load**? If yes, the model's
   structure (AmP capacity + per-axis routing) adds value beyond "more pollution = less
   SFG". If no, the model is only tracking total load — an honest negative.

Controls: body size, and (across-species version) generation time. Rank statistics
throughout (consistent with the rest of the programme).

## Honest weak points to state
- ECOTOX anchors for specific contaminant×*Mytilus* pairs may be sparse → surrogate
  species / QSAR fallback, documented per contaminant.
- Contaminant→axis (mode-of-action) routing is an approximate pMoA assignment.
- A single-contaminant or single-summed-load proxy would **not** test the model's
  specific content (it collapses to "load vs SFG"); the multi-contaminant,
  axis-routed version is the real test — hence the per-contaminant template columns.

## Decision / next step
The model side is buildable now; it waits only on the per-site table. Two ways forward:
1. **You extract** the Widdows 1995 (or Beiras) per-site SFG + contaminant table into
   `sfg_gradient_TEMPLATE.csv` (institutional access), then I build and run the
   ECOTOX-anchored margin-vs-SFG harness end-to-end.
2. **I build the harness now** against the template schema (parameterised, smoke-tested
   on synthetic rows, like the dated-tree pipeline) so it runs the instant the table is
   filled.

This is the genuine external test of the adaptive margin — worth doing right rather
than against a single-load proxy that would test nothing.

## Sources
- [Widdows et al. 1995, MEPS 127:131 (abstract)](https://www.int-res.com/abstracts/meps/v127/p131-148/) · ICES [TIMES 40 — SFG methods](https://repository.oceanbestpractices.org/handle/11329/667)
- [Beiras et al. — SFG large-scale Spanish coast survey](https://pubmed.ncbi.nlm.nih.gov/22885349/)
- EPA ECOTOX (in-repo `data/ecotox/`); AmP project.
