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

## ⚠️ RESULTS (2026-06-13) — second species (*M. galloprovincialis*): does NOT replicate (condition/age confound)
Second gradient added: **Albentosa et al. 2012**, *Sci. Total Environ.* 435–436:430–445,
doi:10.1016/j.scitotenv.2012.07.025 — the Iberian SMP survey of *Mytilus galloprovincialis*
(2007 & 2008, 39 site×survey records). Committed as `sfg_albentosa2012_iberia.csv` (SFG) +
`sfg_albentosa2012_contaminants.csv` (Table 6) + `sfg_albentosa2012_biometric_repeated.csv`
(CI/ST confounders for the 16 repeated-site records, Table 5); run via
`examples/sfg_margin_validation_albentosa2012.jl`. *(The 2026-06-13 handover and the first
draft mislabelled this "Beiras et al." — corrected here; the author is Albentosa.)*

**This is the cleanest possible negative: the authors' own paper diagnoses the failure.**
Albentosa et al.'s headline finding is that SFG here is dominated by **condition index**
(`SFG~CI` r=−0.617\*\*\*, R²=51.7%) and **age** (shell thickness, `SFG~ST` r=−0.465\*\*,
R²=26.4%); chemicals add only **16.95%** of variance, and *no* relationship exists with the
global pollution index (`SFG~CPI` r=0.092, n.s.). Their stepwise model keeps **Zn
*positive*** even with CI+ST included — tissue metal burden indexes something beneficial
(food/essential-element availability), not toxic stress. Only **DDTs and chlordanes**
behave as toxicants (inverse), weakly.

**Our harness reproduces every piece of that, independently:**
- *Empirical signs invert vs Widdows.* `SFG~Zn` ρ=**+0.59\*\*** (pooled), `SFG~Cu` +0.40\*,
  most metals positive; `SFG~CPI` ρ=+0.19 (n.s., matching their r=0.092). PAHs null
  (`SFG~PAH13` +0.12, vs Widdows' −0.47\*\*). Organochlorines the only "expected" sign
  (`DDTs` −0.27, `chlordanes` −0.34\*).
- *Axis diagnostic* (signed; toxic axis should be negative): metal route ρ=**+0.49\*\***
  (confound), organochlorine route ρ=−0.09 (toxic direction, weak — diluted by adding
  PCB7, which the authors also found n.s.), PAH route +0.12.
- *Margin* therefore anti-tracks SFG: `SFG~A_t` ρ=**−0.11** (pooled; per-survey −0.10/−0.07).
  The model is internally correct (pressure erodes margin) — the *pressure proxy* is wrong.
- *Confound control* (16 repeated-site records, partialling CI+ST): our `SFG~CI` ρ=−0.668\*\*
  and `SFG~ST` ρ=−0.491\* **reproduce the authors' coefficients** (a transcription
  cross-check). Controlling CI+ST does **not** rescue the margin (raw −0.45 → partial
  −0.39): the confound lives in the pressure proxy (burden∝food), not as an additive term.

**Reading — a BOUNDING result, not a tuning failure.** No re-routing was done to chase a
negative (that would be p-hacking against an invariant). The margin replicates **where
tissue burden indexes exposure** (Widdows North Sea, hydrocarbon-dominated, ρ=+0.41) and
fails **where burden is condition/food-confounded** (Iberian SMP metals, ρ=−0.11). The
boundary is exactly the long-standing *tissue-burden ≠ exposure* caveat — and the authors
say the same ("SFG biomarker needs corrective strategies to avoid effect of confounding
factors"). Concrete rule for the next gradients: **screen for condition/age/food
confounding (and prefer exposure-based pressure) before trusting tissue burden as the
pressure axis.** This does not overturn Widdows; it bounds it.

**Across-species capacity test:** still NOT testable. We now have 2 of the ~8–10 species'
gradients, but Widdows (+0.41) and Albentosa (−0.11) are not a clean capacity contrast —
the pressure proxy means different things in each (exposure vs confounded burden). The
capacity test needs gradients where burden→exposure holds across species. This survey is
**tissue-only** (no water/sediment concentrations), so it cannot be re-analysed on an
exposure basis to escape the confound.

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
2. **Albentosa et al. 2012, *STOTEN* 435–436:430–445 (Spanish coast, 41 stations)** —
   large-scale SFG + contaminants in wild *M. galloprovincialis*. **DONE** (see the ⚠️
   RESULTS section above): condition/age-confounded, does not replicate. *(Previously
   mislabelled "Beiras et al. 2003" here.)*
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
- [Albentosa et al. 2012, STOTEN 435–436:430–445 — SFG large-scale Spanish coast survey](https://doi.org/10.1016/j.scitotenv.2012.07.025)
- EPA ECOTOX (in-repo `data/ecotox/`); AmP project.
