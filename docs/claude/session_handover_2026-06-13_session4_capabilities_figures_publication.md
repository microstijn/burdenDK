# Session handover (2026-06-13, session 4): life stages, movement, surface:volume, publication figures, wiki consolidation — and the publication-readiness plan

*Cross-session handoff. Continues `session_handover_2026-06-13_session3_pmoa_and_validation_writeup.md`
(session 3, which closed the validation, built the water-quality pMoA foundation, and drafted the
standalone validation manuscript — then resolved that the framework should be published as ONE
framework+validation paper). This session **built three new model capabilities (life stages, movement,
surface:volume TK), turned them into publication figures, consolidated the wiki, made the README
margin-first, and now sets the publication-readiness plan.** Everything below is committed and pushed
to `main` (`b2680e9`), except the intentionally-local manuscript. Self-contained.*

---

## Part 0 — TL;DR (read first)

1. **The publication is one paper: framework + external validation.** It is drafted, compiles (~15 pp),
   margin-first, conceptual-modelling venue: `docs/tex/twotimescaleresilience_framework_paper.tex`
   (**LOCAL / gitignored on purpose** — see [[framework-validation-single-paper]]).
2. **Three new model capabilities shipped this session — additive, tested, but NOT validated:**
   stage-resolved capacity, mobile-target (movement) exposure, surface:volume aquatic toxicokinetics.
   They are flagged in the manuscript (Conceptual Scope + Limitations) and the wiki as
   *implemented-but-not-validated*. They are **Paper-2 material**, not Paper-1 content.
3. **Publication figures built** (margin-first): salmon trajectory + GeoMakie migration map (Hovmöller).
   The salmon **map is now the README landing hero** (replacing the amplification-centric concept
   figure — `F` is null and off-message).
4. **Wiki consolidated** 13→11 pages: 5 validation pages → one `External-Validation.md`; new
   `Life-Stages-and-Movement.md`; Home streamlined; README margin-first and current.
5. **The plan (Part 3):** get Paper 1 submission-ready (journal, authors, bib, SI, submission
   artifacts, Zenodo). Owner decisions: **target journal + authorship** (blocking).

---

## Part 1 — What this session did (all committed + pushed; newest first)

| commit | what |
| --- | --- |
| `b2680e9` (merge) | docs/wiki-consolidation → main |
| `9e27218` | README hero → salmon migration map (drop the grid figure) |
| `c67c7aa` | (superseded) grid margin-map hero + salmon worked example |
| `3da28a7` | **wiki consolidation** (5 validation pages → 1) + currency pass |
| `d79229f` (merge) | feat/migration-visuals → main |
| `9b053e7` | salmon migration **map** figure (encounter-collapsed map + Hovmöller) |
| `ccc97f1` | salmon migration **trajectory** figure |
| `ff1ba3d` (merge) | feat/surface-volume-tk → main |
| `32e635a` | **surface:volume** aquatic TK (gated, length-dependent retention `ρ(L)=ρ_ref^(L_ref/L)`) |
| `0dd345a` (merge) | feat/stage-resolved-capacity → main |
| `6917b69` | **mobile-target** exposure: occupancy-weighted migration (salmon POC) |
| `6fedeef` | **stage-resolved DEB capacity** (life-stage integration) |

Memories updated: [[life-stage-capacity-implemented]] (incl. the **A₀-stays-density decision** and the
surface:volume/migration POC), [[framework-validation-single-paper]], [[always-store-decisions]].

---

## Part 2 — Model facts to carry (don't re-derive)

- **Stage-resolved capacity (additive).** `λ_max = v_eff(L)/L` (younger/smaller recovers faster);
  `λ_min = k_M` and `A₀ = E_m` (a reserve **density**) stay stage-invariant; metabolic acceleration
  `v_eff(L)=v·s(L)` (abj/ssj); the `(1−κ)` axis is relabelled maturation→reproduction at puberty
  (weight unchanged). Opt-in API: `deb_params_at_length`, `deb_params_for_stage(:juvenile/:adult)`;
  `DEBStageProfile`, `v_eff_at_length`, `deb_stage_profile`. The AmP library carries an additive
  `ontogeny` block; regenerated with **zero core-value change** (only the vestigial `KA` dropped).
  Whole-organism path unchanged. **Gotcha:** *M. edulis* is **abj**; adult-at-`L_i` == whole-organism
  relies on `L_i = s_M·L_m`.
- **A₀ decision (user-confirmed):** `A₀` is a density → stage-invariant; the stage effect lives in `λ`,
  NOT `A₀`. Do **not** make `A₀` absolute (∝L³) without sign-off (breaks the threshold-free design).
  The DEB-correct lever for small-stage vulnerability is the **exposure side** (surface:volume) — done.
- **Movement.** `occupancy_weighted_exposure(region_concs, occupancy)` → `C=Σ_g π_g C_g`, fed into the
  existing compound-memory recurrence (a single region = the resident case).
- **Surface:volume aquatic TK.** `surface_volume_retention(ρ_ref,L,L_ref)=ρ_ref^(L_ref/L)` and the
  **gated** `waterborne_stage_retention(...; waterborne)` (inert for terrestrial/air-breathing/dietary
  targets). Rate/lag reading (both uptake & elimination ∝1/L): small stages equilibrate fast, large
  stages lag. No new knob.
- **Validation (unchanged):** recovery/margin layer corroborated (`k_M`, `R_i`, SFG, stress-on-stress,
  transplant); amplification `g`/`F` **null everywhere**; the **across-axis capacity weighting is
  untested** (the central open question). Full account: `docs/wiki/External-Validation.md` +
  `docs/notes/external_validation_synthesis.md`.
- **Julia 1.12.6 via `julia +release`.** Fast suite green (**18870** tests). Figures need CairoMakie;
  the map needs GeoMakie (`GeoMakie.coastlines()`, no external files).

---

## Part 3 — Publication-readiness plan (THE focus)

### Scope decision (recommended)
**Paper 1 = the framework + external-validation manuscript** (already drafted). Keep it tight and
validated. The new **stage/movement/surface:volume** work is unvalidated → **Paper 2** (extensions +
the water-quality spatial application), later. In Paper 1 they stay as flagged future-work (already
written into Conceptual Scope + Limitations). *Open for the owner to override (fold-in vs two papers).*

### Paper 1 — to submission-ready
1. **Owner decisions (blocking):**
   - **Target journal** — recommend *Ecological Modelling* (conceptual modelling; publishes the
     DEBtox/GUTS/AOP–DEB lineage the paper cites). Drives length/format. Alternatives: MEE (length
     limits), ETC/STOTEN (compress the philosophy).
   - **Authorship** — author list, affiliations, CRediT roles, corresponding author, ORCIDs. Replace
     the `(authors)` placeholder.
2. **Bibliography verification** — the folded validation refs use `and others` + plausible-but-unverified
   page/volume fields; verify against DOIs before submission.
3. **Title + abstract** — finalize; ensure the margin-first / honest-validation framing is consistent
   end-to-end (the validation is *corroboration, not strong prediction*).
4. **Figures** — forest plot ✓ in-paper. Optional: a methods schematic. **Graphical abstract:** the
   margin-first salmon **map** (`docs/wiki/figures/salmon_migration_map.png`) is the strongest candidate.
   Decide whether to include *one* brief, explicitly-unvalidated "extensions" illustration (a salmon
   figure) or leave stage/movement as text-only future work (recommended: text-only in Paper 1).
5. **Reproducibility / data-availability statement** — point to the public GitHub repo
   (https://github.com/microstijn/burdenDK) + the synthesis note; cut a **tagged release + Zenodo DOI**
   at submission so the code is citable.
6. **Supplementary information** — fold the framework-comparison supplement (from
   `docs/tex/twotimescales_intro_memory_updated.tex`) and the full derivation
   (`docs/tex/framework_derivation_complete.tex`) into SI; the validation reproducibility detail lives
   in `docs/notes/external_validation_synthesis.md`.
7. **Submission artifacts** — cover letter; highlights / graphical abstract; the journal's author
   checklist.
8. **Collaboration logistics** — the manuscript is **local/gitignored**, which does not scale to
   co-authors. Before circulating, move it to **Overleaf or a private repo** (keep the in-repo figures
   as the source of truth; regenerate via the figure scripts).

### Paper 2 (future, not now)
Stage-resolved capacity + mobile-target exposure + surface:volume TK + the water-quality pMoA spatial
coupling (DynQual). Needs **validation data** (stage- and movement-resolved contaminant–outcome series)
or framing as a methods/application paper. Foundation already in the repo (the salmon POC + figures +
`data/pMoA_Stressor_Routing.csv`).

---

## Part 4 — Open decisions (owner)

1. **Target journal** (blocks formatting/length). Recommend *Ecological Modelling*.
2. **Authorship** (blocks submission).
3. **Scope** — Paper 1 only (recommended) vs fold the new capabilities in vs two papers.
4. **Whether a salmon figure goes in Paper 1** (recommend: text-only future-work; map as graphical
   abstract).
5. **Where the manuscript lives for co-authoring** (Overleaf vs private repo).

---

## Part 5 — Where things live

- **Manuscript (local, gitignored):** `docs/tex/twotimescaleresilience_framework_paper.tex` (+ PDF).
  §3 math is current (threshold-free `E=x/(1+x)`, `A_t=A0(1−Q)`, linear `λ`, `λ_min=k_M`,`λ_max=k_M·g`).
- **Validation:** `docs/wiki/External-Validation.md` (consolidated) + `docs/notes/external_validation_synthesis.md` (full).
- **New-capability docs:** `docs/wiki/Life-Stages-and-Movement.md`.
- **Figures (tracked PNGs; PDFs are build artifacts):** `docs/wiki/figures/salmon_migration.png`,
  `salmon_migration_map.png`, `validation_forest.png`. Scripts: `examples/salmon_migration_scenario.jl`
  (shared), `salmon_migration_demo.jl` (table), `salmon_migration_figure.jl` (trajectory),
  `salmon_migration_map.jl` (map). `examples/validation_forest_plot.jl`, `examples/wiki_figures.jl`.
- **Code (new):** `src/deb_axes.jl` (stage API), `src/amp_library.jl` (ontogeny + stage wrappers),
  `src/AmP_Translator.jl` (ontogeny block), `src/movement_exposure.jl`, `src/exposure_filters.jl`
  (surface:volume). Tests: `test/test_amp_lifestage.jl`, `test/test_movement_exposure.jl`.
- **Memories:** `~/.claude/.../memory/` — `framework-validation-single-paper`,
  `life-stage-capacity-implemented`, `always-store-decisions`.

---

## Part 6 — Git state

All committed and pushed to `main` (https://github.com/microstijn/burdenDK), through `b2680e9`. `main`
is in sync with origin. Working tree carries only the intentionally-local items: the modified
`.gitignore` (the rule that keeps the manuscript local — **uncommitted by design**) and the gitignored
manuscript `.tex`/`.pdf` + figure-PDF build artifacts. Repo workflow: branch → commit →
`git merge --no-ff` into `main` → push (push needs explicit user authorization).

---

## Part 7 — Loose ends

1. **`.gitignore` change is uncommitted** — it is what keeps the manuscript local. Decide whether to
   commit the rule (documents intent) or move it to `.git/info/exclude` (zero repo footprint). Until
   then the manuscript is local simply because it was never committed.
2. **`docs/notes/lambda_min_maintenance_rate.pdf`** — still untracked/unread (carried from sessions 1–3).
3. **`examples/wiki_figures.jl` is stale** — it still uses the removed `KA` half-saturation `λ` form for
   its illustrative restoring-force figure. Not the landing figure anymore, but regenerate/fix it before
   relying on those wiki concept images.
4. **GitHub Wiki-tab mirror** — `docs/wiki/` is the source of truth; the Wiki *tab* is a manual re-sync.
