# Session handover (2026-06-14, session 5): validation extensions, the validation-map figure, and the scale-free / no-calibration framing

*Cross-session handoff. Continues `session_handover_2026-06-13_session4_capabilities_figures_publication.md`
(session 4: built life-stage/movement/surface:volume capabilities, the salmon figures, consolidated the
wiki, made it margin-first, and set the publication-readiness plan). This session **extended the external
validation** (a bridge figure, a capacity-structure ablation, a response-curve robustness check, the
phenanthrene dynamics firm-up grown to two species, and a 5-fish cross-species capacity pilot), **rebuilt
the validation figure as an architecture↔forest "validation map,"** and **nailed the interpretive spine:
the model is validated as a zero-calibration, scale-free description of capacity erosion and recovery —
not of absolute capacity.** Everything below is committed + pushed to `main` (`0c51858`) except the
intentionally-local manuscript. Self-contained.*

---

## Part 0 — TL;DR (read first)

1. **The interpretive spine is now explicit (this was the missing story).** The framework validates three
   **scale-free** properties — a recovery **rate** (`k_M`), a **relative** margin state, and erosion
   **dynamics** — using **zero free parameters** (nothing fitted to outcomes). It does **not** claim the
   **absolute** margin level `A_0` (which cancels within-species and never enters a validated result), the
   **across-axis capacity weighting** (untested), or any single-number amplification scalar (null). The
   clean boundary: where absolute cross-species capacity *would* matter (clam > mussel baseline anoxia
   tolerance) the model honestly **fails**, exactly as the scale-free reading predicts.
2. **The validation figure is now a "validation map"** (`validation_architecture.png`): left = the model
   chain (capacity → margin state → recovery λ → function → dynamics → amplification) with TeX formulas,
   row-aligned + colour-linked to the forest on the right, so a reader sees *what part of the model each
   anchor validates*. A "scale-free claims" header + a grey "not established" band carry the framing.
   It **replaces** the plain forest in the manuscript (`fig:forest`).
3. **Dynamics firmed up to two species** (Dellali 2023, *M. galloprovincialis* + *R. decussatus*): each an
   independent 12-cell grid, ρ(erosion,LT50)=**−0.99 / −0.97**; control-normalised cross-species ρ=−0.91.
   Static map predicts flat; dynamic erosion reproduces the dose-ordered decline. Proof-of-concept, cleaner
   and now 2-species.
4. **New honest results:** a within-anchor **ablation** (routed margin beats naive load everywhere — the
   *structure* earns its keep, the *weighting* still doesn't); a **response-curve robustness** check
   (ranking stable, rank ρ≥0.99, under Hill/exp curves — no knob); a **DOME bridge figure** (relative
   ranking of a real monitoring network); a **5-fish capacity pilot** (n=5 — underpowered, confirms the
   open question).
5. **Manuscript:** 19 pp, compiles clean. Abstract, validation intro (4th commitment = zero free
   parameters + a "what is / is not claimed" paragraph), and the new figure caption all carry the
   scale-free / no-calibration stance. Still local/gitignored.

---

## Part 1 — What this session did (all committed + pushed to `main`; newest first)

| merge | what |
| --- | --- |
| `0c51858` | **validation-map figure** (architecture↔forest, TeX, scale-free header + not-established band) + no-calibration framing (abstract, validation intro, "what is/isn't claimed" paragraph, new caption) |
| `f98e4d5` | **2-species dynamics** (add *R. decussatus* clam) + **5-fish benzovindiflupyr capacity pilot**; forest gained the phenanthrene anchors |
| `16bea13` | phenanthrene dynamics **full 4×3 grid** (figure-digitised 15/21 d) → ρ=−0.99, n=12, non-trivial |
| `48e714b` | phenanthrene dynamics **firm-up** (Dellali 2023), endpoints PoC |
| `f19c12d` | **response-curve-form sensitivity** — ranking robust to the impairment curve |
| `d3a2938` | **DOME bridge figure** + capacity-weighting **ablation** + named the dynamics firm-up target |

Memories updated: [[framework-validation-single-paper]] (interface/response split, no-calibration decision,
the validation-map figure, the bridge/ablation/dynamics artifacts), [[always-store-decisions]].

---

## Part 2 — Model / validation facts to carry (don't re-derive)

- **Zero free parameters.** Every model quantity (`A_0`, axis weights `w_i`, `λ_min=k_M`, `λ_max=k_M g`)
  is fixed offline from AmP/DEB + ECOTOX; pressure is threshold-free (`E=x/(1+x)`); **nothing is fitted to
  the outcomes**. The reported ρ's are *uncalibrated*. A single fitted "fulcrum" would very likely raise
  them and is **deliberately omitted** (corroboration, not quantitative prediction). This is the headline
  honesty point now.
- **`A_0` is never validated in absolute terms, and never needs to be.** Within-species tests: `A_0` is a
  constant and **cancels** in the rank correlation (the dynamics harness notes the linear scale `q` is
  ratio-invariant). Cross-species recovery: a **rate** (`k_M`) is correlated, not a level. So validated =
  scale-free (rate / relative state / dynamics).
- **Two-species dynamics (Dellali 2023, *Animals* 13(1):151).** Constant waterborne phenanthrene
  (~10/45/89 µg/L, Table 1; flat over time) → assimilation axis. LT50 from text/Table 3 (7 d, 28 d) +
  Figure-4-digitised (15, 21 d). Each species independent 12-cell grid; ρ(erosion,LT50)=−0.99 (mussel),
  −0.97 (clam). The clam's **higher baseline anoxia tolerance** (~13 vs ~9 d) is anoxia physiology (shell
  closure, anaerobic metabolism — the paper's own explanation), **not** the contaminant margin; raw pooled
  ρ=−0.82 conflates it, control-normalised ρ=−0.91. *That gap is the untested across-species capacity
  question, made concrete.*
- **5-fish benzovindiflupyr pilot (Nickisch Born Gericke et al. 2022, ETC 41(7):1732, SI; acute data
  Ashauer 2013).** Per-species 96 h LC50 extracted (carp 3.5 → bluegill 28.5 µg/L). `k_M`→sensitivity
  −0.90 raw, −0.92 | size — but **n=5** and AmP only exposes **structural** length `L_m` (a weak
  cross-species size proxy, δ_M varies), so it **cannot** overturn the powered n=310 control. The
  distinctive axis weighting (`alpha_maint`) does **not** predict sensitivity (single chemical = one axis).
  **Confirms the open question stands.** SDHI → maintenance routing.
- **Ablation:** routed/structured margin beats a naive equal-weight load at every field anchor (W1995
  +0.41/+0.22; W2002 +0.12/+0.005; DOME +0.39/+0.32) — supports the *operative structure* (routing,
  impairment, aggregation), **not** the across-species *weighting*.
- **Response-curve robustness:** swapping `E=x/(1+x)` for Hill (h=0.5,2) or `1−e^{-x}` (all threshold-free,
  half-saturating at the reference — **no knob**) barely moves the corroboration (ρ within ±0.03) and
  leaves the ranking near-identical (rank ρ≥0.99). The relative ranking doesn't depend on the curve.
- **DOME bridge figure:** 17-station ICES DOME network ranked by modelled margin, coloured by survival-in-
  air; ρ=+0.39 raw (reproduces the SoS anchor). Validated machinery, **no water concentration** — the
  licensed relative use, the meeting point of water-quality monitoring and ecological response.
- **Julia 1.12.6 via `julia +release`.** All harnesses run; figures need CairoMakie. **TeX in figures
  works** via MathTeXEngine (already a Makie dep — no dependency added; use `L"..."`). Manifest carries
  `LaTeXStrings` + `MathTeXEngine`.

---

## Part 3 — Next validation steps (ranked by leverage)

### 1. The across-axis capacity weighting — THE central open question (highest value)
The distinctive content is the per-species **across-axis weighting** `w_i` (how burden on each DEB axis
erodes margin). It is **untested**: every individual-level anchor is single-species (weighting held
constant), the one cross-species single-trait link (`k_M`→toxicity) is size-confounded, and the 5-fish
pilot is single-MoA (one axis) and n=5. **A real test needs multiple modes of action × the same species
set**, so that the model's *axis* weighting (not just overall capacity) predicts *which species is more
vulnerable to which MoA*.
- **Concrete build:** assemble a **multi-MoA × species matrix from EPA ECOTOX** — chemicals with clear
  energetic pMoA spanning axes (assimilation: PAHs/narcotics; maintenance: metals/uncouplers/SDHIs;
  reproduction: organotins/EDCs; growth: chitin-synthesis inhibitors) × a shared species set, with acute
  LC50 (or, better, a shared sublethal endpoint). Test: does the AmP axis weighting predict the
  *species × MoA* sensitivity pattern **beyond body size**? The §7b ECOTOX pipeline
  (`scripts/extract_ecotox_acute.awk`, `scripts/state_axis_ecotox_amp.jl`) and the benzovindiflupyr
  pipeline (`scripts/extract_benzovindiflupyr_fish.jl`, `examples/benzovindiflupyr_capacity_probe.jl`) are
  the starting harnesses — generalise to multiple MoA.
- **Watch:** body size dominates acute LC50 (the n=310 lesson); a clean size control is essential, and AmP
  **structural** `L_m` is a poor cross-species size proxy (use physical length `L_m/δ_M` or ultimate
  weight `W_w` — needs δ_M, which is in the AmP raw records, not the current `auxiliary_metrics`).

### 2. A calibration / quantitative-prediction study — the "fulcrum" (the user's point)
The validation is deliberately uncalibrated. A **separate** study should introduce **one** global
calibration (e.g. a per-outcome slope, or a single global gain on the response) and quantify how much it
raises the fits and whether the model then makes *calibrated quantitative* predictions. This answers "is
ρ≈0.2–0.45 a structural ceiling or just missing scale?" and is the natural complement that the
no-calibration stance sets up. Keep it **out of Paper 1** (which is the conservative corroboration); it is
a Paper-2 / methods contribution. Do **not** smuggle a fulcrum into Paper 1.

### 3. Powered dynamics — beyond proof-of-concept
The phenanthrene 2-species result is a clean PoC but single-PAH, half figure-digitised, and dose-uniform
in fractional erosion. To make it powered: (a) **mine more reported SoS LT50 tables** from the
Mediterranean biomonitoring school (Dellali/Boufahja/Khessiba/Banni often report LT50 in tables across
species/time); (b) a **clean toxicodynamic recovery-rate (`k_r`) dataset** across species (GUTS-proper;
scarce, *G. pulex*-centric — Nyman 2012). Lower priority than (1) — the dynamics are already 2-species
corroborated.

### 4. (Paper 2) the water-quality coupling + the new capabilities
Validate the **designed interface** (exposure filter + compound memory `B_t`) against spatial fields
(DynQual / monitoring), and the **stage / movement / surface:volume** capabilities against stage- or
movement-resolved contaminant–outcome data. Big application step; explicitly deferred.

### What NOT to chase
- GUTS **direct-survival** datasets for the *dynamics* — they don't give the discriminating static-vs-
  dynamic contrast (survival declining under constant exposure is trivially expected). They are useful for
  the **capacity** question (step 1), not the dynamics.
- "Fixing" the clam/mussel absolute gap by adding an anaerobic-metabolism axis or a tolerance covariate —
  that risks a tuning knob (violates the invariants). Document it as the scope boundary instead.

---

## Part 4 — Where things live (new/changed this session)

- **Validation figure (manuscript):** `examples/validation_architecture_forest.jl` →
  `docs/wiki/figures/validation_architecture.png` (+ `docs/tex/validation_architecture.pdf`, build
  artifact). The plain forest (`examples/validation_forest_plot.jl` → `validation_forest.png`) is retained
  but **no longer the manuscript figure**.
- **Bridge figure:** `examples/dome_margin_ranking_figure.jl` → `docs/wiki/figures/dome_margin_ranking.png`.
- **Dynamics (2 species):** `examples/sos_dynamic_validation_dellali.jl`;
  data `data/external/sos_dellali2023_phenanthrene.csv` (species column; provenance in header).
- **Fish capacity pilot:** `scripts/extract_benzovindiflupyr_fish.jl` (parses the SI xlsx — an xlsx is a
  zip of XML; unzip + point `DUMP` at it) → `data/external/benzovindiflupyr_fish_survival.csv`;
  `examples/benzovindiflupyr_capacity_probe.jl`.
- **Response-curve robustness:** `examples/response_curve_sensitivity.jl`.
- **Synthesis (the working record):** `docs/notes/external_validation_synthesis.md` — §7-bis (2-species
  dynamics), §7c (fish pilot), §9b (ablation + bridge + response-curve), §9 (no-calibration bullet),
  §2 scorecard rows.
- **Wiki:** `docs/wiki/External-Validation.md` (scorecard + dynamics paragraph + bridge + response-curve).
- **Manuscript (local/gitignored):** `docs/tex/twotimescaleresilience_framework_paper.tex` — abstract
  (scale-free/no-calibration sentence), validation intro (4th commitment + "what is/isn't claimed"
  paragraph), `fig:forest` swapped to the architecture map + new caption, dynamics §, scorecard, future-
  work. **19 pp, compiles clean** (`pdflatex → bibtex → pdflatex ×2`; `rm ttr_framework_refs.bib` to pick
  up bib edits — `filecontents*` won't overwrite it).
- **Raw SI** for the fish data (`etc5348-sup-0004`) is **not** in the repo (raw downloads gitignored); the
  derived CSV is committed.

---

## Part 5 — Git state

All committed + pushed to `main` (https://github.com/microstijn/burdenDK), through `0c51858`. `main` in
sync with origin. Working tree carries only the intentionally-local items: the modified `.gitignore` (the
manuscript-local rule, uncommitted by design), the gitignored manuscript `.tex`/`.pdf`, and figure-PDF
build artifacts (`validation_architecture.pdf`, `dome_margin_ranking.pdf`, `salmon_*.pdf`). Workflow:
branch → commit → `git merge --no-ff` → push (**push needs explicit user authorization**; the auto-mode
classifier enforced this once this session when a push wasn't authorized).

---

## Part 6 — Open decisions / loose ends

1. **Owner decisions still blocking submission (from session 4, unchanged):** target journal (the
   interface framing favours *Ecological Modelling*, or the bridge venues *Environmental Modelling &
   Software* / STOTEN), and **authorship** (replace the `(authors)` placeholder).
2. **Bibliography verification** — folded validation refs use `and others` + plausible-but-unverified
   vol/pages (checklist in the `.bib` header). Verify before submission. The new `dellali2023phenanthrene`
   ref is verified.
3. **The figure swap is in the local manuscript only.** If a co-author setup is created (Overleaf/private
   repo per session 4), carry the architecture figure + its caption.
4. **`examples/wiki_figures.jl` is stale** (still the removed `KA` λ form) — carried from earlier sessions.
5. **`docs/notes/lambda_min_maintenance_rate.pdf`** — still untracked/unread (carried from sessions 1–4).
6. **Body-size proxy for cross-species capacity:** AmP `auxiliary_metrics` exposes structural `L_m` but not
   `δ_M`/`W_w`; the cross-species capacity tests (step 1, 5-fish pilot) need a physical-size proxy, which
   means pulling `δ_M`/`W_w` from the AmP raw records into the library (a small `AmP_Translator.jl` add).
