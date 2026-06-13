# Session handover (2026-06-13, session 3): negative controls, dated-tree PGLS, water-quality pMoA coupling, validation manuscript

*Cross-session handoff. Continues `validation_handover_margin_sos_dynamics_2026-06-13.md` (session 2,
which validated the adaptive margin across state/function/dynamics and ended on the open "powered
dynamic test" decision). This session **answered that decision (it's data-starved/size-confounded),
closed the real-phylogeny control, built the water-quality pMoA-coupling foundation, and drafted a
standalone validation manuscript.** Self-contained: nothing below should need re-deriving. Everything
is committed and pushed to `main` (through merge `af8d090`).*

---

## Part 0 — TL;DR (read first)

1. **The session 2 open decision is resolved — negatively, and that's the finding.** The proposed
   "powered dynamic test" (AmP `k_M` → GUTS `k_D`) was scoped and **run**: the single-trait
   maintenance→toxicity link is **body-size-confounded on both axes**.
   - *State axis* (`k_M`→acute LC50/EC50, ECOTOX, 4 AChE inhibitors, **n=310**): replicates Baas &
     Kooijman raw (ρ≈−0.27, all 4 chemicals) but **nulls under a body-size control** (partial ≈ −0.03).
     A robust null, not underpower.
   - *Rate axis* (`k_M`→Rubach 2010 chlorpyrifos elimination `k_out`, n=6): weak/wrong-sign, n.s.; and
     `k_out` is the *toxicokinetic* rate, not the thesis-relevant *toxicodynamic* recovery `k_r` (which
     is scarce, *G. pulex*-centric, chemical-specific → powered dynamic test stays data-starved).
   - **Consequence:** the model's distinctive leverage is the **across-axis capacity weighting**, NOT
     `k_M` as a scalar predictor. The capacity weighting remains **untested** (carried as an assumption).
2. **Real-phylogeny control — DONE.** A dated **TimeTree** (182 spp) PGLS: Pagel's λ≈0 again (phylogeny
   was never the confound). The `k_M`↔recovery signal **survives in rank form** under pace + dated
   phylogeny (β\*=0.221, **p=0.011**); the log-linear form nulls (p=0.96). So the headline anchor holds
   as a rank/monotone effect — exactly the scorecard's "rank-robust, magnitude-modest,
   specification-sensitive."
3. **Water-quality coupling foundation built.** Aggregate water-quality stressors (BOD, TDS, microplastics,
   …) are routed to DEB axes by **physiological mode of action (pMoA)**, not tuning weights — a declared
   mechanism (`data/pMoA_Stressor_Routing.csv` + `src/pmoa_stressor_routing.jl`). New wiki page
   `Water-Quality-Coupling.md`.
4. **THE metals finding (carry this).** Fitted DEBtox pMoAs put **metals on *assimilation*** (cadmium
   A×6, copper A×3, …), contradicting the field routing (metals→maintenance). **Tested empirically
   (option 3): metals→assimilation degrades ALL field anchors** (Widdows 1995 +0.41→+0.30; Widdows 2002
   +0.12→−0.14; DOME +0.39/+0.45→+0.27/+0.26). Resolution: **both routes coexist by data regime** —
   field measured burden → maintenance (confound-control, *keep `src/mode_of_action.jl` unchanged*);
   controlled/exposure-derived burden → assimilation (fitted pMoA). Documented in the table + digest.
5. **A standalone validation manuscript was drafted** (`docs/tex/external_validation_manuscript.tex`,
   compiles to 6 pp, with a forest-plot figure). It's the paper the framework paper flagged as priority
   #1. Needs: authors, journal, bib verification.
6. **Next actions (roadmap):** `docs/notes/where_next_roadmap_2026-06-13.md` — (A) finish & submit the
   manuscript [recommended first]; (B) test the capacity weighting [data-blocked]; (C) build the
   water-quality application [design-dependent].

---

## Part 1 — What this session did (chronological, all committed; 11 merges)

1. **GUTS `k_D` dynamic test — scoped + run** → size-confounded negative controls (state n=310, rate
   n=6). `docs/notes/guts_kd_dynamic_test_scoping.md`. Merge `1b85797`.
2. **External-validation synthesis write-up** — added §7b (the negative controls) and consolidated.
   Merge `f85e7ac`.
3. **Dated-tree PGLS** (TimeTree, 182 spp) — linear + rank; `k_M` rank-survives. Merges `f85e7ac`,
   `7b6c2c6`. Scripts `scripts/comadre_pgls_dated.jl`, `comadre_pgls_dated_rank.jl`.
4. **Water-quality coupling audit** — found the DynQual path uses heuristic proxy→axis weights (a
   prototype, not settled); the validated chain is tissue-burden→saturation→CAS-MoA.
5. **pMoA stressor-routing foundation** — table + loader + basis doc. Merge `65098c4`.
6. **ECOTOX effect-code pMoA diagnostic — NEGATIVE** (κ-cascade → endpoint ≠ pMoA). Kept on record.
   Merge `2f4a55e`. + `docs/notes/pmoa_evidence_to_gather.md`.
7. **Read assembled Tier-1 evidence** (`~/Downloads/dat_1`): fitted-pMoA compilation (ESPI 2018
   c7em00328e), Jager 2020, EnviroTox (Kienzler 2019), Billoir 2007. Extracted
   `data/external/debtox_fitted_pmoa_c7em00328e.csv`; `docs/notes/pmoa_evidence_digest.md`. Merge
   `b63b57b`.
8. **Metals routing tension → option 3 empirical test** → metals→assimilation degrades all field
   anchors; keep maintenance. Merge `53443eb` (digest §4-RESULT).
9. **Populated the routing table** (pmoa_fitted / field_confound_route / references columns) + **wiki
   refresh** (Phylogenetic-PGLS dated result, new Water-Quality-Coupling page, Margin-Validation §5,
   Home, Limitations, Data-and-Parameters). Merge `13f75a9`.
10. **Standalone validation manuscript** drafted (LaTeX, 6 pp). Merge `ba0efd6`. **Forest-plot figure**
    added. Merge `6c2cc8f` (`examples/validation_forest_plot.jl`).
11. **Strategic roadmap** `docs/notes/where_next_roadmap_2026-06-13.md`. Merge `af8d090`.

---

## Part 2 — Honest scorecard (current)

| anchor | level | result | status |
| --- | --- | --- | --- |
| COMADRE (`k_M`,`R_i`) | population | `k_M`↔recovery **+0.19–0.22\*** (rank; survives pace + **dated-tree** PGLS, β\*=0.221 p=0.011; linear nulls); `R_i`↔compensation **+0.77\*\*** | ✅ corroborated |
| GlobTherm (n=664) | physiology | capacity coherent (\|ρ\|≤0.45), recovery-specific (general-resilience refuted) | ✅ bounding |
| Scope for Growth ×3 | energetics | **+0.41 → +0.12 → −0.11** (scale-dependent) | ✅ / ◐ |
| Stress-on-Stress (DOME) | energetics | **+0.39 → +0.45** (confound-controlled) | ✅ static map |
| Viarengo 1995 | controlled dose | dose-response + additive mixture | ✅ controlled |
| Transplant (Veldhuizen) | dynamics | continued erosion static can't; Cd-alone **+0.90** | ◑ proof-of-concept |
| **`k_M`→toxicity (n=310)** | **cross-species** | **raw −0.27; nulls under body size (partial −0.03)** | **✅ bounding (negative control)** |
| amplification `g`/`F` | — | **null everywhere** | ✅ (margin-first) |

**The through-line:** rank-robust, magnitude-modest, **specification-sensitive**. A proper size/phylogeny
control nulls the single-trait stories; the `k_M` rate anchor survives only in rank form. The distinctive
content (capacity weighting) is untested.

---

## Part 3 — Model facts to carry (don't re-derive)

- **Julia 1.12.6 via `julia +release`** (default LTS cannot load the project).
- **Static API:** `compute_adaptive_margin_response`: `E=x/(1+x)` per axis → `Q=ΣwᵢEᵢ` → `A_t=A0·(1−Q)`.
  Pressure `x=conc/median` (threshold-free). Used by every SFG/SoS/ECOTOX correlation.
- **Dynamic API:** `simulate_deb_axis_response` integrates the slow erosion state, relaxing at
  `λ∈[k_M, k_M·g]`, time-constant `1/λ ≈ months` for *M. edulis* (unfitted).
- **MoA routing (FIELD):** hydrocarbons/PAH→assimilation; **metals→maintenance** (confound-control —
  see Part 4); organotin+organochlorine→reproduction; DBT→growth. `src/mode_of_action.jl`,
  `src/moa_deb_mapping.jl` — **unchanged this session and must stay so** (the validation depends on it).
- **pMoA routing (water-quality / aggregate stressors):** `data/pMoA_Stressor_Routing.csv` +
  `src/pmoa_stressor_routing.jl` — declared mechanism, not weights. Temperature is a **rate modifier**,
  not a pressure. `ρ` (tissue memory) only for accumulating stressors; ambient stressors erode via `λ`.
- **Dated-tree PGLS needs `Distributions`** in a throwaway env (script headers show the command). Tree
  at `data/external/comadre_amp_dated_tree.nwk` (TimeTree, 184 tips).
- **Repo workflow:** branch → commit → `git merge --no-ff` into `main` → push (each push needs explicit
  user authorization; the auto-classifier blocks autonomous pushes). The session used 11 such merges.

---

## Part 4 — The metals dual route (the crux — empirically settled)

- **Fitted DEBtox pMoA:** metals act via **assimilation** (gut/gill damage → reduced uptake; cadmium
  A×6, copper A×3, mercury A, uranium A in the ESPI 2018 compilation).
- **Field validated routing:** metals → **maintenance**, *deliberately*, because measured metal tissue
  burden is a **positive confound** (tracks food/condition, not exposure: tissue burden ≠ exposure).
- **Option 3 test (this session):** re-routing metals→assimilation **degraded every SFG/SoS anchor**
  (`docs/notes/pmoa_evidence_digest.md` §4-RESULT). So:
  - **field measured tissue burden → maintenance** (keep it);
  - **controlled / exposure-derived burden → assimilation** (the fitted pMoA — relevant if a
    water-quality coupling drives burden from *modelled* water concentration; that's an open judgment).
- The routing table carries BOTH (`primary_axis=maintenance`, `pmoa_fitted=assimilation`,
  `field_confound_route` note). **Do not silently switch `mode_of_action.jl`.**

---

## Part 5 — The validation manuscript (where it is)

- `docs/tex/external_validation_manuscript.tex` (+ committed PDF) — standalone, compiles
  (`pdflatex`+`bibtex`, 6 pp), self-contained embedded bib (`validation_refs`), matches the existing
  `article`/`natbib`/`plainnat` preamble. Figure: `validation_forest.pdf`
  (`examples/validation_forest_plot.jl` regenerates it + the wiki PNG).
- **Structure:** abstract → intro (margin-first logic) → methods (static+dynamic API, MoA routing,
  controls, rank stats) → results (7 anchors + scorecard table + forest plot) → discussion (the
  signature; tissue-burden≠exposure + metals lesson; capacity weighting as assumption; water-quality
  scope) → reproducibility.
- **Left to do:** author block (placeholder `(authors)`); **target journal** (Ecol. Modelling / ETC /
  STOTEN — drives length); **verify bib** page/volume fields (filled with plausible values); optional
  methods schematic + the `F`-vs-`g` figure.

---

## Part 6 — Open decisions & next actions (from the roadmap)

Full costed plan: `docs/notes/where_next_roadmap_2026-06-13.md`.
1. **(A) Finish & submit the manuscript** — recommended first; near-term, low-risk. Needs journal +
   authors + bib verification.
2. **(B) Test the capacity weighting** — the distinctive *untested* claim; **data-blocked** (needs
   across-species contaminant-gradient data with a shared outcome — largely doesn't exist). The central
   research question; scope a minimal multi-species SFG/SoS design.
3. **(C) Build the water-quality application** — foundation ready; **pin the target water-quality model
   + its outputs first** (DynQual: BOD/TDS/FC/temp), then build a clean coupling reading the pMoA table,
   decide the exposure-vs-field metals route, and **sensitivity-test** the routing (ranking stability =
   the defensibility statement).

**Owner decision points:** journal + authorship; which water-quality model; whether to invest in
purpose-built across-species data.

---

## Part 7 — Key new files (this session)

- **Docs/notes:** `guts_kd_dynamic_test_scoping.md`, `pmoa_evidence_to_gather.md`,
  `pmoa_evidence_digest.md`, `pmoa_stressor_routing_basis.md`, `where_next_roadmap_2026-06-13.md`.
- **Data:** `data/pMoA_Stressor_Routing.csv` (+ loader `src/pmoa_stressor_routing.jl`, wired into the
  module), `data/external/debtox_fitted_pmoa_c7em00328e.csv`, `state_axis_ecotox_amp_paired.csv`,
  `ecotox_acute_4chem.csv`, `rubach2010_kM_kout_paired.csv`, `comadre_amp_dated_tree.nwk`,
  `comadre_pgls_dated_results.txt`, `comadre_pgls_dated_rank_results.txt`,
  `ecotox_pmoa_inference_results.txt`, `comadre_species_for_timetree.txt`.
- **Scripts:** `comadre_pgls_dated.jl`, `comadre_pgls_dated_rank.jl`, `state_axis_ecotox_amp.jl`,
  `rubach2010_rate_axis.jl`, `extract_ecotox_acute.awk`, `ecotox_pmoa_inference.{awk,jl}`.
- **Manuscript / figure:** `docs/tex/external_validation_manuscript.{tex,pdf}`,
  `docs/tex/validation_forest.pdf`, `examples/validation_forest_plot.jl`,
  `docs/wiki/figures/validation_forest.png`.
- **Wiki updated:** `Phylogenetic-PGLS.md`, `Water-Quality-Coupling.md` (new), `Margin-Validation.md`,
  `Home.md`, `Limitations-and-Open-Questions.md`, `Data-and-Parameters.md`.
- **Source PDFs supplied by user** (NOT committed; in `~/Downloads/dat_1`): `c7em00328e1.xlsx` (fitted
  pMoA), `jager2020.pdf`, `etc4531.pdf` (EnviroTox), `billoir2007.pdf`, Rubach 2010/2011.

**Run:** harnesses `julia +release --project=. examples/<name>.jl`; PGLS needs a `Distributions`
throwaway env (script headers); manuscript `pdflatex→bibtex→pdflatex×2` in `docs/tex/`.

---

## Part 8 — Git state

All work **committed and pushed to `main`** (https://github.com/microstijn/burdenDK). 11 merges this
session, newest first: `af8d090` (roadmap), `6c2cc8f` (forest plot), `ba0efd6` (manuscript), `13f75a9`
(routing table + wiki), `53443eb` (metals option-3), `b63b57b` (fitted pMoA), `2f4a55e` (ECOTOX
diagnostic), `65098c4` (pMoA routing), `7b6c2c6` (rank PGLS), `f85e7ac` (synthesis + dated PGLS),
`1b85797` (GUTS scoping). Working tree clean except `.claude/` and one untracked PDF (Part 9).

---

## Part 9 — Loose ends

1. **`docs/notes/lambda_min_maintenance_rate.pdf`** — dropped by the user at the very start of session 1
   or 2; **still untracked and unread** across this session. Possibly relevant to the `λ_min=k_M` core
   (cf. `docs/notes/lambda_min_maintenance_rate.tex`). Decide whether to read/commit it.
2. **Tier-3 evidence gathering (user's task):** aggregate-stressor ecophysiology (microplastics,
   hypoxia, salinity, acidification, pathogens) for the routing rows with no DEBtox fit; EnviroTox
   consensus MOA for bulk CAS→class. List: `docs/notes/pmoa_evidence_to_gather.md`.
3. **GitHub Wiki-tab mirror** — `docs/wiki/` is the source of truth; the Wiki *tab* is a manual re-sync.
4. **(low) Regenerate `AmP_Species_Library.json`** to drop the vestigial `KA` field (cosmetic).
