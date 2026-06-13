# pMoA evidence digest — assembled papers (2026-06-13)

*Digest of the evidence set assembled in `~/Downloads/dat_1`, read against the pMoA stressor-routing
table (`data/pMoA_Stressor_Routing.csv`). The extracted fitted pMoAs are in
`data/external/debtox_fitted_pmoa_c7em00328e.csv`. **This digest records findings and one important
tension; it does NOT change the routing table or the codebase CAS-MoA routing — that needs a
decision (see §4).***

## 1. The files

| file | what it is | role |
| --- | --- | --- |
| `c7em00328e1.xlsx` | **fitted-pMoA compilation** (Kienzler/Jager-adjacent; ESPI 2018, doi 10.1039/c7em00328e). 49 DEBtox fits / 28 compounds | **Tier 1 gold** — the per-compound fitted pMoA |
| `jager2020.pdf` | Jager 2020, "Revisiting simplified DEBtox" (Ecol. Modelling 416:108904) | framework: the 4 pMoAs + the starvation cascade |
| `etc4531.pdf` | Kienzler et al. 2019, consensus MOA in the **EnviroTox** database (ETC 38:2294; open access) | Tier 2 — bulk CAS→MoA-class (narcosis vs specific) |
| `billoir2007.pdf` | Billoir et al. 2007, **DEBtox + matrix population models**, Daphnia/cadmium (Ecol. Modelling 203:204) | the individual→population bridge; pMoA-choice → population effect |
| `1-s2.0-S0273230014002803` | Bhatia et al. 2015, Cramer/TTC classification (Toxtree vs QSAR Toolbox) | tangential — classifier-disagreement caution, **not** pMoA |
| `rubach2011.pdf`, `…2010 Rubach.pdf` | Rubach traits (2011) + chlorpyrifos TK (2010) | Tier 4 trait-based / already used (synthesis §7b) |

## 2. The framework (jager2020) — and why it backs the ECOTOX diagnostic negative

Jager confirms toxicant stress is limited to **four metabolic processes: assimilation, maintenance
(somatic + maturity), growth, reproduction** (+ hazard). Crucially, he notes **toxicant-induced
starvation**: a toxicant hitting assimilation or maintenance leaves too little for the rest, so
growth and reproduction suffer downstream. That is exactly the **κ-rule cascade** that made our
ECOTOX effect-code diagnostic default to "reproduction" for almost every compound — independent
confirmation that *most-sensitive endpoint ≠ pMoA*, and that fitted pMoAs (which model the whole
budget) are required.

## 3. Fitted-pMoA findings (the data) — codes A=assim, M=maint, G=growth, R=repro, H=hazard

Dominant fitted pMoA per compound (`debtox_fitted_pmoa_c7em00328e.csv`):
- **Metals → ASSIMILATION-dominant:** cadmium **A×6** (A/M×1, G×1), copper **A×3** (G, M+R), mercury
  **A** (A/M), uranium **A×2** (A/M, M+G). Zinc mixed (M, A/M).
- **Organotins mixed:** tributyltin **A**, triphenyltin **M** (not reproduction).
- **PAHs mixed:** benzo(k)fluoranthene **A**, PAH-mixture **A**; fluoranthene **G+R**, pyrene **R**,
  nonylphenol **G+R**.
- **AChE/pesticides mixed:** aldicarb **M**, atrazine **M**, chlorpyrifos **R** (n=1; note this matches
  the reproduction signal our ECOTOX diagnostic threw — but both are weak).
- **Ocean acidification (pH) → M** ✓ (confirms the table's acidification→maintenance primary).

**Variability caveat (important).** The same compound shows *different* fitted pMoAs across
studies/species (cadmium A / A/M / G; copper A / G / M+R). The compilation itself flags that "in
several cases the identification of MoA is not very formal." So **pMoA is species/study/method-
dependent, not a fixed compound constant.** Use the *modal* pMoA with explicit confidence, and prefer
species-aware assignment where the model is species-resolved (AmP).

## 4. THE TENSION — fitted metals→assimilation vs the validated field routing → maintenance (decide before applying)

The biggest finding contradicts both the draft table **and the existing codebase CAS-MoA routing**:
- **DEBtox fits:** metals act via **assimilation** (gut/gill damage → reduced uptake).
- **Existing validated routing (`src/mode_of_action.jl`, used in SFG/SoS):** metals → **maintenance**,
  chosen *deliberately* to keep metals **off** the toxic (hydrocarbon→assimilation) axis, because in
  the mussel field data metals (As/Cd/Zn) are a **positive confound** (they track food/condition, not
  toxic exposure — "tissue burden ≠ exposure"). The routed margin beat naive load *partly because*
  metals were kept off assimilation.

So routing metals→assimilation per the fits would load a positive-confound burden onto the toxic
axis and **could break the SFG/SoS validation**. These are not reconcilable by fiat. Options:
1. **Keep them separate concerns (recommended framing):** the *fitted pMoA* (assimilation) is the
   true lab mechanism; the *field metals→maintenance* is a data-specific **confound-control**, not a
   pMoA claim. Use fitted pMoA for controlled/mechanistic application; retain confound-handling for
   field tissue-burden data. Document both explicitly.
2. **Species/context-dependence:** metals may genuinely act differently in mussel filter-feeders vs
   Daphnia; the compilation's variability supports this. Make the assignment species-aware.
3. **Settle empirically:** switch metals→assimilation and **re-run the SFG/SoS harnesses** to see
   whether the validation degrades. The data decides.

**Do not silently change `mode_of_action.jl` or the routing table.** This is a scientific call for the
project owner.

### 4-RESULT — Option 3 RAN (2026-06-13): metals→assimilation DEGRADES every field anchor

Re-routed metals from maintenance to assimilation in all four field harnesses (the one-line change:
fold `tu(metals…)` into `p_assim`, set `p_maint = fill(1.0, …)`) and re-ran. **Every dataset got
worse:**

| dataset | baseline (metals→maintenance) | variant (metals→assimilation) |
| --- | --- | --- |
| Widdows 1995 SFG | **+0.405\*** | +0.295 (loses significance) |
| Widdows 2002 SFG | +0.116 | **−0.137** (sign flips wrong) |
| Albentosa 2012 SFG | −0.114 | −0.206 (worse) |
| SoS DOME (raw) | +0.392 | +0.267 |
| SoS DOME (length + condition) | **+0.452** | +0.259 (control no longer strengthens) |

**Verdict: the fitted-pMoA metals→assimilation does NOT transfer to field tissue-burden data — it
degrades all five tests.** Clearest driver: in the Albentosa variant the metals are a demonstrable
**positive** confound (Cu **+0.40\***, Cd +0.29 *with* SFG) — tissue metal burden tracks food/condition,
not toxic exposure — so loading them onto the toxic axis imports that confound (and dilutes the PAH
signal that the assimilation axis is meant to carry).

**Resolution (option 1, now empirically backed): both assignments stay, for different data regimes.**
- **Field tissue-burden data → metals = maintenance** (a confound-control, *not* a pMoA claim; tissue
  burden ≠ exposure). **KEEP `src/mode_of_action.jl` and the field harnesses unchanged** — the
  validation depends on it.
- **Controlled / exposure-based data → metals = assimilation** (the fitted DEBtox pMoA, valid where
  burden = exposure; e.g. Viarengo, or an exposure-driven water-quality coupling).

So the routing table should carry metals with `pmoa_fitted = assimilation` *and* a
`field_confound_route = maintenance` note. **No change to the validated codebase.** (Reproduce: apply
the one-line metals→`p_assim` swap to any field harness and re-run; full numbers above.)

## 5. Bonus resources unlocked
- **EnviroTox consensus MOA (etc4531)** — open-access, consensus across 4 schemes, large CAS coverage
  (40% narcosis / 17% specific / 43% unclassified). Use to triage many compounds to a coarse MoA
  class; but it is *acute fish MoA*, coarser than DEBtox pMoA — a screen, not a final axis.
- **Billoir 2007** — a worked **DEBtox→Leslie-matrix→population-growth-rate** bridge (Daphnia/cadmium)
  that *quantifies how the pMoA choice changes the population outcome*. Directly supports the COMADRE
  individual→population argument, and shows why getting the routing right matters at population level.

## 6. Suggested next steps
1. **Decide the metals tension (§4)** — owner call; offer to run the empirical SFG/SoS re-test.
2. Add a `references` + `pmoa_fitted` column to `data/pMoA_Stressor_Routing.csv`; populate chemical
   rows from `debtox_fitted_pmoa_c7em00328e.csv` (raise `confidence` where a fit exists).
3. Pull the **EnviroTox** database (etc4531 supplement) for bulk CAS→MoA-class triage.
4. Keep the variability caveat visible: report modal pMoA + per-compound spread, not a single label.
