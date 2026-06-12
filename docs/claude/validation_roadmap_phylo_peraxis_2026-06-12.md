# Validation roadmap — phylogeny + per-axis recovery (session handoff)

*Working note, 2026-06-12. Cross-session handoff. The current session validated the
framework against external data (COMADRE) for the first time and reframed the model
to be margin-first. This document is the briefing for a **fresh session** to carry
the validation further along two lines: (A) a real phylogeny + PGLS, and (B) a
per-axis margin recovery metric mapped to the demographic-resilience decomposition.
It is deliberately extensive and self-contained so no context has to be re-derived.*

---

## Part 0 — Where we are (read this first)

### The one-paragraph state of the model
`burdenDK` / TwoTimescaleResilience maps chronic environmental pressure → erosion of
a physiological **adaptive margin** `A_t` → weakened **restoring force** `λ(A_t)` →
amplified burden of a later acute event, `F = λ(A0)/λ(A_t)`. Over this session we (1)
found `F` collapses to a single parameter, fixed the cause (`λ_min` was
mis-normalized: `[p_M]/[E_m]` → the textbook maintenance rate constant `k_M =
[p_M]/[E_G]`, which makes the timescale ratio the **energy investment ratio `g`**
instead of the allocation fraction `κ`); (2) showed even `g` reduces to ~`1/{p_Am}`
and is **not** a real biological vulnerability axis; (3) **reframed the model so the
adaptive-margin *state* is the product and `F` is a derived "happy analytical
accident"**; (4) added capacity-aware margin features + a margin-first clustering
standardiser; and (5) ran the **first external validation** against COMADRE.

### The headline external result (the thing to build on)
Matching AmP species to the COMADRE animal matrix database (now **197 species**
after GBIF name harmonisation — see Part 3, ✅) and predicting demographic recovery
(log damping ratio `|λ1|/|λ2|`):

| model quantity | raw ρ | \| gen. time | \| gen. time + Order |
| --- | --- | --- | --- |
| `λ(A0)` recovery rate | +0.362 ** | +0.173 * | +0.089 (n.s.) |
| **`λ_min = k_M`** | +0.406 ** | +0.264 ** | **+0.190 \*** |
| `g` (amplification) | −0.109 | −0.128 | −0.055 (n.s.) |

*(n=197, 193 with gen. time; was 183/179 — the signal is stable under harmonisation.)*

**`k_M` (the DEB maintenance rate constant) predicts demographic recovery beyond
both pace-of-life and coarse phylogeny.** The amplification scalar does not. This is
the cleanest external evidence to date and it lands on the **recovery/margin layer**,
not the `F`/`g` readout — consistent with the margin-first reframe.

### Files you will touch / reuse (all on `main`)
- `scripts/extract_comadre_recovery.jl` — **standalone** COMADRE extraction. Needs
  `RData.jl` + `DataFrames` in a *throwaway* env (NOT the project Manifest):
  `julia -e 'using Pkg; Pkg.activate(mktempdir()); Pkg.add(["RData","DataFrames"])'`.
  Downloads `COMADRE_v.4.26.4.0.RData` (CC-BY) and writes
  `data/external/comadre_recovery.csv` (committed; raw `.RData` is gitignored).
  Currently outputs: `species, n_matrices, mean_log_damping, mean_generation_time,
  class, order, family`.
- `examples/comadre_partial_validation.jl` — the analysis (project env). Matching,
  raw/partial/Order-controlled Spearman. **This is the file to extend for Idea A/B.**
- `docs/notes/comadre_partial_validation.md` — the result writeup + caveats.
- `docs/notes/external_anchor_scouting.md` — why COMADRE, what other anchors exist,
  the circularity trap.
- `docs/notes/g_lifehistory_check.md`, `feature_redundancy_check.md`,
  `mixture_assumption_sensitivity.md` — the internal diagnostics.
- AmP capacity parameters live in `data/AmP_Species_Library.json`; the mapping is in
  `src/AmP_Translator.jl` (offline; `auxiliary_metrics` now includes `k_M`, `E_G`,
  `g`, `L_m`, `p_Am`, `p_M`).

### Environment reminder (do not skip)
**Use Julia 1.12.6 via `julia +release`.** The default LTS (1.10.x) cannot load the
project (`UndefVarError: StaticData`). See `CLAUDE.md`.

### The honest current caveats (what these two ideas address)
1. The phylogenetic control is a **taxonomic-rank proxy** (Order group-mean-centering),
   not a real tree / PGLS. → **Idea A.**
2. The model side is a **single scalar recovery rate** (`λ(A0)` or `k_M`), which throws
   away the per-axis margin structure that the margin-first reframe says is the point.
   → **Idea B.**

---

## Part 1 — Idea A: real phylogeny + PGLS

### Why
Related species are not statistically independent: a correlation across species can be
driven entirely by a few clades. The current Order group-mean-centering removes
*between-Order* mean differences but is coarse (27 Orders, treats within-Order
structure as flat, ignores branch lengths). A **phylogenetic generalized least
squares (PGLS)** or **phylogenetic independent contrasts (PIC)** analysis uses the
actual tree (topology + branch lengths) as the covariance structure, which is the
field standard for comparative analyses. If `k_M`'s partial signal (ρ=0.20*) survives
proper PGLS, the result is publishable-grade; if it vanishes, we learn the Order
proxy was too lenient.

### Data — getting a tree for the 183 matched species
Options, best first:
1. **Open Tree of Life (OTL)** synthetic tree via API. The `rotl` R package
   (`tnrs_match_names` → `tol_induced_subtree`) is the canonical route; there is also
   a REST API (`https://api.opentreeoflife.org/v3/`) usable from Julia via `HTTP.jl`
   + `JSON`. OTL gives topology but **branch lengths are not dated** — you'd need to
   either (a) use `compute.brlen` (Grafen) for an ultrametric proxy, or (b) graft on
   dated branch lengths.
2. **TimeTree.org** — dated trees; can upload a species list and download a Newick
   with divergence times. Manual but gives real branch lengths. Best for a defensible
   dated tree.
3. **Clade-specific published timetrees** — VertLife (birds/mammals/squamates/
   amphibians/fish: `birdtree.org`, `vertlife.org`), the Open Tree dated "DaTeD"
   variants. Highest quality where coverage exists, but COMADRE animals span many
   phyla (corals, insects, molluscs, fish, birds, mammals…), so a single source won't
   cover all 183 — likely a **mixed-source supertree** or restrict to the
   best-covered subclade (e.g. chordates) for a clean first pass.

**Recommendation for the first pass:** restrict to **vertebrates** (or even just
mammals+birds, where dated trees are excellent and COMADRE coverage is densest), get
a dated tree from VertLife/TimeTree, and run a clean PGLS there. Then, separately,
do an OTL-based all-taxa PGLS with Grafen branch lengths as a coarser robustness
check. Report both.

### Tooling
- **Julia** phylo ecosystem is usable but thinner: `Phylo.jl` (trees, can read
  Newick), `PhyloNetworks.jl` (has `phylolm` = PGLS!), `GLM.jl`, `StatsModels.jl`.
  `PhyloNetworks.phylolm` does PGLS directly given a tree + data frame. This keeps
  everything in Julia.
- **R fallback** (no R currently installed in this env — would need `winget install`
  or similar): `ape` (read tree, PIC), `caper`/`phylolm`/`nlme::gls(correlation=
  corPagel)` for PGLS with estimated λ (Pagel's lambda). R is the field standard and
  has the most mature tooling; if Julia phylo proves painful, export the matched
  table to CSV and do PGLS in R.

### Concrete pipeline
1. Extend `examples/comadre_partial_validation.jl` (or a new
   `examples/comadre_pgls_validation.jl`) to write the matched table to CSV:
   `species, k_M, lambda_A0, g, comadre_log_damping, generation_time, class, order,
   family` (the model + COMADRE columns already computed there).
2. Fetch a tree for those species (OTL via `HTTP`/`JSON`, or TimeTree/VertLife
   download). Save Newick to `data/external/comadre_amp_tree.nwk` (gitignore raw
   downloads, commit the final pruned Newick if small).
3. Prune/match tree tips to the species list (handle name mismatches: OTL TNRS, or a
   GBIF/Catalogue-of-Life synonym pass — see Part 3).
4. PGLS: `comadre_log_damping ~ k_M + generation_time` with phylogenetic covariance
   (Pagel's λ estimated). Report the partial effect (coefficient, t, p) of `k_M`
   controlling generation time *and* phylogeny. Repeat for `lambda_A0` and `g`.
   Also report Pagel's λ (how much phylogenetic signal is in the residuals).
5. Robustness: PIC version; Grafen vs dated branch lengths; vertebrate-only vs
   all-taxa.

### What confirms / refutes
- **Confirm:** `k_M` retains a significant positive partial effect under PGLS. Then
  "the DEB maintenance rate constant predicts demographic recovery independently of
  pace-of-life and phylogeny" — a strong, defensible external-validation sentence.
- **Refute:** `k_M`'s effect drops to n.s. under PGLS. Then the Order-proxy signal was
  phylogenetic after all; report it honestly and lean on the *internal* coherence of
  the margin-first model instead.

### Pitfalls / cautions
- **Coverage**: not all 183 species will be in any one tree; report n actually used.
- **Branch lengths**: undated OTL topology with Grafen lengths is a weak covariance;
  prefer dated trees where possible and say which you used.
- **Polytomies**: OTL synthetic trees have many; PGLS tolerates them but PIC needs
  resolution (randomly resolve + average, or use `multi2di`).
- **Name matching is the silent killer** (see Part 3); a 183→tree match can silently
  drop a third of species.
- **`k_M` units**: it is a per-day rate from AmP; COMADRE damping ratio is per
  projection interval (often per year). Rank/PGLS on monotone transforms is fine, but
  for PGLS use a sensible transform (log `k_M`, log generation time) and check
  residual normality.

---

## Part 2 — Idea B: per-axis margin recovery metric

### Why
The whole point of the margin-first reframe is that the **per-axis** margin state
(which DEB process is impaired — assimilation / maintenance / growth / reproduction)
carries information that the scalar `F` and even scalar `λ` throw away. So far the
external validation used only *scalar* model quantities (`λ(A0)`, `k_M`, `g`). Idea B
asks: **does an axis-resolved recovery description predict demographic recovery — and
specifically its components — better than a single scalar?** This is the test that
would most directly vindicate (or not) the margin-first claim.

### The key opportunity: align 4 DEB axes ↔ 3 demographic-resilience components
COMADRE matrices decompose into **`matU`** (survival/growth transitions) and
**`matF`** (fertility). The **demographic resilience** framework
(Capdevila/Salguero-Gómez; the `Rage` R package has the functions) decomposes
resilience into three components:
- **Resistance** — ability to resist decline after a disturbance (≈ short-term, driven
  by survival/structure → `matU`).
- **Compensation** — ability to over-perform / boom after disturbance (≈ driven by
  reproduction/fertility → `matF`).
- **Recovery time** — time to return to stable structure (≈ damping ratio, what we
  already used).

The model's four DEB axes map naturally onto this:
- **maintenance + growth axes** ↔ survival/structure ↔ **resistance** (`matU`).
- **reproduction axis** ↔ fertility ↔ **compensation** (`matF`).
- **assimilation axis** ↔ overall energy throughput ↔ scales both.

So the rich hypothesis is: **the model's maintenance/growth capacity predicts
demographic *resistance*, and its reproduction capacity predicts demographic
*compensation*** — a structured, multi-dimensional correspondence, not a single
correlation. If that pattern holds, it is far stronger evidence than a scalar ρ.

### Defining a per-axis model recovery metric (options)
The model does not currently expose an axis-resolved recovery rate; you must define
one. Options, in increasing ambition:

- **B1 — axis weights as a capacity profile.** The κ-rule weights
  `w = [½, κ/4, κ/4, (1−κ)/2]` already are an axis-resolved capacity allocation.
  Test whether `w_reproduction` predicts compensation, `w_maintenance + w_growth`
  predicts resistance, etc. Cheapest; but the weights are mostly κ, so this largely
  re-tests κ per component.
- **B2 — axis-localised recovery sensitivity.** For each axis `a`, impair only that
  axis to a reference level, compute the resulting `λ` drop:
  `s_a = λ(A0) − λ(A0(1 − w_a·E_ref))`. This 4-vector `(s_A, s_M, s_G, s_R)` is the
  per-axis recovery sensitivity. Correlate each `s_a` with the matching demographic
  component. Moderate effort; uses the existing response functions.
- **B3 — DEB rate constants per process.** Pull the *process-specific* DEB rates from
  AmP (`allStat` has many: `k_M` maintenance, `k_J` maturity maintenance, `v`
  conductance, `kap_R` reproduction efficiency, `R_i` ultimate reproduction rate,
  ages `a_b`/`a_p`). Map each to a demographic component and test. Most faithful to
  DEB; requires extending `AmP_Translator.jl` (or a side-extraction) to surface these
  per species. **This is the most promising** because `k_M` already validated, and
  reproduction-side rates (`R_i`, `kap_R`) are the natural predictors of compensation.

**Recommendation:** do **B3** for the reproduction-side (predict compensation from
reproduction rates) and keep `k_M` for resistance/recovery, i.e. a small structured
test: resistance↔`k_M`, compensation↔reproduction-rate, recovery↔damping (already
done). Report the 3×3 correlation matrix (model component × demographic component) and
look for the diagonal being strongest.

### Data — the demographic components
- Compute **resistance / compensation / recovery** per COMADRE matrix from `matA`,
  `matU`, `matF`. The `Rage` R package (`Rage::...`) has these; in Julia you'd
  implement them from the matrices (the damping ratio is already done; resistance and
  compensation are transient-dynamics metrics — first-timestep amplification /
  attenuation, reactivity, maximal amplification `Kreiss`/`rho_max`). References:
  Stott et al. transient bounds; Capdevila et al. demographic resilience.
- Extend `scripts/extract_comadre_recovery.jl` to also output, per species:
  reactivity / first-timestep attenuation (resistance proxy) and maximal amplification
  (compensation proxy), alongside the damping ratio (recovery).

### What confirms / refutes
- **Confirm (strong):** the model component × demographic component correlation matrix
  has a clear diagonal — maintenance/growth capacity ↔ resistance, reproduction
  capacity ↔ compensation — surviving pace-of-life + phylogeny controls. That would be
  a genuinely novel, multi-dimensional external validation of the *margin-state*
  (not just a scalar).
- **Refute:** the matrix is diffuse / off-diagonal; the axes don't map to demographic
  components. Then the per-axis structure, while conceptually appealing, isn't
  externally distinguishable — important to know before the manuscript leans on it.

### Pitfalls
- Transient metrics (reactivity, amplification) are sensitive to the chosen initial
  stage distribution and to matrix imprimitivity — filter and document carefully.
- The axis↔component mapping is a *hypothesis*; pre-register it (write it down before
  computing) to avoid fishing.
- B3 requires faithfully extracting reproduction-side DEB rates; check units and that
  they vary independently of `k_M` (else it's circular again).

---

## Part 3 — Residual refinements (smaller, do alongside)

- **Species name harmonisation. ✅ DONE (2026-06-12).** `scripts/resolve_comadre_amp_names.jl`
  (standalone, throwaway HTTP+JSON env) harmonises COMADRE→AmP names: exact →
  duplicated-genus-typo fix → trinomial→binomial → GBIF Backbone synonym/accepted
  resolution. Writes the committed map `data/external/comadre_amp_namemap.csv`
  (`comadre_species,amp_key,method`), which `comadre_partial_validation.jl` now reads
  (collapsing pseudoreplicated AmP keys by averaging). Recovered **15 species
  (183→197 matched, 193 with gen. time)**; the `k_M` signal is **stable** under the
  larger sample (within-Order partial 0.190* vs 0.200*). The 89 still-unresolved are
  genuinely absent from AmP (≈39 congeners, ≈50 whole genera/clades — corals,
  sponges, molluscs). **The name map is the prerequisite for tree matching (Part 1).**
- **Matrix-quality filter sensitivity.** Re-run with stricter/looser filters (composite
  vs individual matrices, min study duration, primitivity tolerance) and confirm the
  `k_M` signal is robust.
- **Scale-bridge formalisation.** State the individual-DEB → population-matrix link
  explicitly (DEB-IPM / DEB-structured population models; Smallegange, de Roos). The
  damping ratio is a population quantity; `k_M` an individual rate. The bridge is
  defensible but should be argued, not assumed, in the manuscript.
- **Effect sizes & multiple testing.** Report confidence intervals (bootstrap over
  species) and correct for the several model quantities tested.

---

## Part 4 — The big modelling questions still open (not validation, but adjacent)

These are unresolved and will recur; a fresh session should know them:
1. **`KA = 0.3·A0`** is still an undocumented constant (the last no-knob violation).
   Secondary to the λ-bounds but unjustified. See `docs/claude/...source_audit...md`
   Findings 2–3c.
2. **Should `F` exist at all?** The margin-first reframe says the margin state is the
   product; `F` is a lossy scalar. The validation supports this (recovery layer
   validates, `F`/`g` doesn't). A future decision: demote/remove `F` from the headline
   entirely.
3. **Follow-up A (engineering):** the grid / ECOTOX / ISIMIP pipelines still use the
   inert *raw* subtractive margin, so the spatial maps don't yet reflect the
   nondimensional margin-first model. Flagged repeatedly; deliberately de-scoped.
4. **The Kooijman questions** (from the source audit): are the four DEB axes a
   legitimate vulnerability coordinate system without the dynamics; is the
   AmP-parameter extraction legitimate once removed from the dynamical model.

---

## Appendix — reproducibility quick reference

```powershell
# 1. (once) throwaway env for COMADRE extraction
julia +release -e 'using Pkg; Pkg.activate(mktempdir()); Pkg.add(["RData","DataFrames"]); using RData, DataFrames'
# then run the extractor with that env active (it downloads COMADRE + writes the CSV)
julia +release --project=<that-env> scripts/extract_comadre_recovery.jl

# 2. the validation analysis (project env)
julia +release --project=. examples/comadre_partial_validation.jl

# 3. internal diagnostics (project env)
julia +release --project=. examples/amp_kappa_collapse_diagnostic.jl
julia +release --project=. examples/amp_lambda_structure_comparison.jl
julia +release --project=. examples/amp_g_lifehistory_check.jl
julia +release --project=. examples/feature_redundancy_check.jl
julia +release --project=. examples/mixture_assumption_sensitivity.jl
```

**Data:** `data/external/comadre_recovery.csv` (committed, 286 species, with
taxonomy). Raw `COMADRE_v.4.26.4.0.RData` (~907 KB, CC-BY) is gitignored — re-fetch
via the extractor.

**Recommended first move for the next session:** Part 3 name-harmonisation (cheap,
boosts every downstream n), then **Idea A restricted to vertebrates with a dated tree**
(cleanest, most defensible), then **Idea B3** (the structured axis↔component test —
the highest-upside result). Keep every step generation-time-controlled; the
circularity trap is the recurring danger.
