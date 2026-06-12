# Validation programme — session handover (2026-06-13)

*Cross-session handoff. The previous handover
(`validation_roadmap_phylo_peraxis_2026-06-12.md`) is now largely **executed**; this
document records what this (very long) session did, the **conceptual pivot** it produced,
the current honest state of the model + its external validation, and the **one open
decision** for the next session. Self-contained: no context should need re-deriving.*

---

## Part 0 — TL;DR (read this first)

1. **Model cleanup (Part 4) — DONE.** The last no-knob violation `KA = 0.3·A0` is
   **removed**: the recovery curve `λ(A)` is now **linear** in the relative margin,
   `λ(A) = λ_min + (λ_max−λ_min)·clamp(A/A0,0,1)`, so pristine margin → `λ_max` exactly and
   `Fmax = λ_max/λ_min = g`. The amplification scalar `F` is **demoted** (docs/framing only;
   code kept as a labelled diagnostic). Fast test suite green (18,808).
2. **Remaining roadmap (Part 3 + Idea A/B follow-ups) — DONE.** Bootstrap CIs +
   Benjamini-Hochberg, matrix-quality filter sensitivity, scale-bridge note, positive-`a_p`
   resolution, and an Idea-A **dated-tree PGLS pipeline** (ready; needs one manual tree
   download — the single remaining manual step there).
3. **THE PIVOT.** A user question exposed that COMADRE validated the margin curve's
   **rate endpoints** (`k_M`, `R_i`) — **not the adaptive margin itself** (its state, or its
   erosion under chronic pressure). The user's actual deliverable is **the adaptive margin**.
   So the programme turned to validating *the margin*.
4. **Margin validation — two external checks run this session:**
   - **GlobTherm (thermal tolerance), n=664:** the AmP capacity axis is *coherent* with
     independent physiology (|ρ|≤0.45, not noise) **but** recovery capacity is **NOT** general
     resilience — it does not predict thermal-tolerance breadth (refuted; bounding result).
   - **Scope for Growth (Widdows et al. 1995, 36 *Mytilus edulis* sites): the real test.**
     The **modeled adaptive margin tracks measured SFG, ρ=+0.405\*** — the **first direct,
     same-level (no scale bridge) external corroboration of the margin itself**. MoA-routed
     margin beats naive aggregation (0.41 vs 0.22) but not the single best contaminant
     (0.47) — expected, since one species fixes the capacity weighting.
5. **OPEN DECISION (where we stopped):** to test the **capacity weighting** (which the
   single-species SFG test holds constant), we need an **across-species SFG set**. A clean
   *multi-species × chemical-contaminant* SFG dataset wasn't found ready-made; the choice is
   (a) **chemical route** — assemble more single-species contaminant-SFG gradients (next:
   Beiras/Spanish *M. galloprovincialis* survey), or (b) **temperature route** — a
   multi-species SFG×temperature study (adjacent facet, partly overlaps GlobTherm).
   **Recommended (a).** Awaiting the user's steer + the Beiras table. See Part 5.

---

## Part 1 — Model state after this session (what changed in `src/`)

### Linear recovery curve (KA removed) — commit `6074895`
- `src/deb_axes.jl`: `DEBAxisParams.KA` field **removed**; `restoring_force_from_margin`
  and `restoring_force_from_margin_and_axes` now linear:
  `λ_min + (λ_max−λ_min)·clamp(A/A0, 0, 1)`.
- `src/background.jl`: same change; `BackgroundParams.KA` removed.
- `src/amp_library.jl`: `KA` no longer required/read (legacy JSON `KA` ignored).
- `src/AmP_Translator.jl`: stops computing/writing `KA`. **NB: the committed
  `data/AmP_Species_Library.json` still carries a now-ignored vestigial `KA` field** — NOT
  regenerated, because Julia `JSON.print` key ordering is non-deterministic and a full
  regen would produce a spurious whole-file reorder diff. A future clean regen will drop it.
- Tests updated (`test_amp_translator_identities` → `Fmax==g`; `test_amp_library`,
  `test_amp_species_profile` [also fixed a stale pre-re-anchoring `λ_min`],
  `test_amp_species_archetypes`, `test_deb_axes`, `test_ecotox_library`); archetype builder
  dropped `KA`/`A0_over_KA` columns.
- Docs: `CLAUDE.md` invariant note updated; wiki Equations/Pipeline/Limitations/
  Data-and-Parameters updated; `F` reframed as a derived diagnostic that is null externally.

### `F` demotion
Docs/framing only — `F`/`amplification_from_margin`/`amplification_factor` and the `F_t`
feature dimension are **retained** (labelled diagnostic). The margin **state** is the product.

### Prior context that still holds (from the 2026-06-12 handover)
- `λ_min` was earlier re-anchored to `k_M = [p_M]/[E_G]` (fixing the κ-collapse); the
  timescale ratio is the **energy investment ratio `g`**, not `1/κ`.
- **Use Julia 1.12.6 via `julia +release`** (default LTS cannot load the project).

---

## Part 2 — External validation: the honest scorecard

| test | what it validates | result | status |
| --- | --- | --- | --- |
| COMADRE scalar | `k_M` (recovery floor) ↔ demographic recovery | +0.19\* beyond pace + Order | ✅ |
| COMADRE per-axis (Idea B) | `R_i` (reproduction) ↔ compensation | **+0.77\*\*** (beyond pace+size) | ✅ |
| COMADRE PGLS (Idea A) | `k_M` ↔ recovery under real phylogeny | OTL+Grafen tree too weak (Pagel λ≈0); rank-vs-linear fragility | ⚠️ dated tree TODO |
| amplification `g`/`F` | (the scalar readout) | **null in every test** | ✅ (margin-first) |
| GlobTherm | capacity axis = general resilience? | **refuted** (capacity is recovery-specific) | ✅ bounding |
| **SFG (Widdows)** | **the adaptive margin itself** | **modeled margin ↔ measured SFG +0.40\*** | ✅ **first direct margin evidence** |

**The big picture:** external support lands on the **margin/recovery layer**; the
amplification scalar never validates. The SFG result is the first to corroborate the
**margin state under pressure** (the thing the user wants to use) at the margin's own
organisational level (individual energetics → **no scale bridge**, unlike COMADRE).

### What is STILL not validated (be honest in the manuscript)
- The **capacity weighting** of the margin (A0 / κ-rule axis weights): the single-species
  SFG test holds it constant. → needs across-species SFG (Part 5).
- A **real phylogeny** control for the scalar `k_M` result (dated tree, Part 4 below).
- Everything is **rank-based / modest** (ρ≈0.2–0.4) and **specification-sensitive**
  (`k_M`↔recovery is robust in ranks, absent in log-linear regression).

---

## Part 3 — Part 3 robustness + Idea B follow-up (commit `d00bc15`)

- **Effect sizes + multiple testing** (`examples/comadre_bootstrap_effectsizes.jl`):
  bootstrap (resample-over-species, seed 20260612) 95% CIs + Benjamini-Hochberg over 7
  tests. Every positive finding survives BH with a CI excluding 0; `g` is the lone null
  (CI spans 0). (Self-contained t/incomplete-beta p-values, validated against the known t.)
- **Matrix-quality filter sensitivity** (`scripts/comadre_filter_sensitivity.jl`):
  `ρ(k_M, recovery | gen)` stable 0.18–0.33 across 6 COMADRE filter variants — not a
  filter artifact.
- **Scale-bridge note** (`docs/notes/comadre_scale_bridge.md`): argues the individual-DEB →
  population-matrix link via DEB-IPM / PSPM (licenses the rank-correlational tests).
- **Positive-`a_p` RESOLVED** (`scripts/comadre_ap_diagnostic.jl`): the pre-registered
  *negative* `a_p`→compensation came out positive because `a_p` is pace-loaded
  (ρ(`a_p`,gen)=+0.50 — the naive intuition lives in pace). The residual-after-pace signal
  is genuinely positive, is **not** a fecundity proxy (ρ(`a_p`,`R_i`)=−0.13) and **not** a
  matrix-dimension artifact (ρ(dim,comp)≈0; survives a dimension control). Within a pace
  class, delayed maturity independently predicts compensation — a reproduction-*timing* axis.
- Writeups: `docs/notes/comadre_robustness_effectsizes.md`, updated
  `docs/notes/comadre_peraxis_validation.md`.

---

## Part 4 — Idea A dated-tree PGLS (commit `27c9e0a`) — pipeline READY, one manual step

`scripts/comadre_pgls_dated.jl` parses a **dated** Newick (real branch lengths → VCV;
Pagel's λ by ML) and re-runs the recovery models. Smoke-tested end-to-end. **190/198
matched species are vertebrates**, so a VertLife/TimeTree vertebrate timetree covers
~everything. **The one blocker:** no reliable public *dated-tree* API (datelife
unreachable; VertLife/TimeTree behind UIs). A human drops one dated Newick into
`data/external/comadre_amp_dated_tree.nwk` (the script prints the species list + path),
then the PGLS runs unattended. This is the genuine phylogenetic test the OTL+Grafen pass
(Pagel λ≈0) could not provide. See `docs/notes/comadre_pgls_validation.md`.

---

## Part 5 — THE OPEN FRONTIER: validate the margin's *capacity* via across-species SFG

### Where we are
- Single-species SFG (Widdows, *M. edulis*) validated the **erosion mechanism + MoA
  routing** but holds the AmP **capacity** constant. To test capacity we need SFG across
  **multiple AmP species**.
- **AmP coverage is not the constraint** — 15 bivalves are in AmP: `Mytilus_galloprovincialis`,
  `Perna_perna`, `Perna_viridis`, `Mya_arenaria`, `Cerastoderma_edule`, `Macoma_balthica`,
  `Magallana_gigas` (=*Crassostrea gigas*), `Ostrea_edulis`, `Mercenaria_mercenaria`,
  `Ruditapes_philippinarum`, `Ruditapes_decussatus`, `Ensis_directus`, `Spisula_subtruncata`,
  `Modiolus_modiolus`, `Pinna_nobilis`.
- A ready-made **multi-species × chemical-contaminant** SFG table (the ideal) was **not
  found**; multi-species SFG studies mostly use **temperature**.

### The decision (pending the user)
- **(a) Chemical route [recommended].** Extend Widdows cross-species by assembling more
  single-species **contaminant-SFG gradient** studies for AmP bivalves; test whether AmP
  capacity predicts the cross-species **sensitivity** (SFG-decline slope per unit
  contaminant). Tests the margin in its native pressure. **Next species: Beiras et al.
  Spanish "Mussel Watch" SFG + contaminant survey of *M. galloprovincialis* (~40 stations)**
  — the user was asked to pull that table (same shape as the Widdows CSVs). Caveat: a
  well-powered capacity correlation needs ~8–10 species' gradients → ongoing lit-assembly;
  each added study has institutional-access/extraction cost (same wall as Widdows).
- **(b) Temperature route.** A multi-species SFG×temperature study (e.g. van Erkom Schurink
  & Griffiths 1992 — 4 South African mussels, only 2 in AmP; or Sarà's Mediterranean DEB
  work; NE-Atlantic DEBIB set — but that is DEB-*modelled*, circular-risk). Readily
  available but tests the **thermal/energetic** facet, not the chemical margin, and partly
  overlaps GlobTherm (though sublethal SFG-under-warming ≠ lethal CTmax, so still distinct).

### Honest design notes for whoever builds it
- The capacity test = does AmP capacity (A0, κ-weights) predict cross-species margin
  resilience? With single-species gradients, n is the number of *species* (low) — needs
  several species for power.
- **Circularity guard:** validate against an *independent* outcome (measured SFG), never
  the EC50/NOEC that feeds the model's pressure. Avoid DEB-*modelled* energy budgets
  (e.g. DEBIB) as the "outcome" — that is AmP-like and circular.
- Body burden vs exposure-based potency is a real mismatch (critical body residue); the
  Widdows harness used a threshold-free **median-normalised relative burden** as the
  pressure proxy (so per-contaminant potency is not encoded). Document whatever choice.

---

## Part 6 — Files, commands, reproducibility

### New this session
| file | role |
| --- | --- |
| `examples/comadre_bootstrap_effectsizes.jl` | bootstrap CIs + BH (Part 3) |
| `scripts/comadre_filter_sensitivity.jl` | matrix-quality filter sensitivity (RData env) |
| `scripts/comadre_ap_diagnostic.jl` | positive-`a_p` resolution (RData env) |
| `scripts/comadre_pgls_dated.jl` | dated-tree PGLS (Distributions env) — needs the tree |
| `scripts/extract_amp_for_globtherm.jl` | AmP↔GlobTherm matched table (MAT env) |
| `examples/globtherm_validation.jl` | thermal coherence probe |
| `examples/sfg_margin_validation.jl` | **the SFG margin test** (project env) |
| `data/external/globtherm_amp_matched.csv` | committed (raw `GlobalTherm.csv` gitignored) |
| `data/external/sfg_widdows1995_northsea.csv` | SFG outcome (36 sites) |
| `data/external/sfg_widdows1995_contaminants.csv` | contaminant pressure (same sites) |
| `data/external/sfg_gradient_TEMPLATE.csv` | turnkey template for new SFG gradients |

### Notes (all in `docs/notes/`)
`comadre_robustness_effectsizes.md`, `comadre_scale_bridge.md`,
`margin_validation_scouting.md` (the anchor-scouting index for the margin),
`globtherm_validation.md`, `sfg_validation_status.md` (SFG **results** + design).

### Commands
```powershell
# margin SFG test (project env)
julia +release --project=. examples/sfg_margin_validation.jl
# thermal coherence (project env)
julia +release --project=. examples/globtherm_validation.jl
# COMADRE effect sizes (project env)
julia +release --project=. examples/comadre_bootstrap_effectsizes.jl
# throwaway-env extractors (see each script header for the Pkg.add line):
#   RData+DataFrames: comadre_filter_sensitivity.jl, comadre_ap_diagnostic.jl
#   MAT:              extract_amp_for_globtherm.jl
#   Distributions:    comadre_pgls_dated.jl
```
All raw downloads gitignored (`*.RData`, `GlobalTherm.csv`); derived CSVs committed.

### Wiki
`docs/wiki/` is the **source of truth**; the GitHub Wiki tab is a mirror (re-mirror via
clone → copy → rewrite `[[links]]`/`../` → push to `master`). **NOT yet mirrored this
round:** the GlobTherm + SFG margin-validation results live only in `docs/notes/` — a
follow-up should add a wiki "Margin validation" page (GlobTherm + SFG) and re-sync.

---

## Part 7 — Prioritised next actions
1. **Across-species SFG capacity test (Part 5)** — get the user's route decision (a/b);
   if (a), obtain the Beiras *M. galloprovincialis* table, run a 2-species test
   (`examples/sfg_margin_validation.jl` generalises), keep adding species.
2. **Dated-tree PGLS (Part 4)** — one manual VertLife/TimeTree download → runs unattended.
3. **Mirror GlobTherm + SFG to the wiki** (new "Margin validation" page).
4. (Lower) regenerate `AmP_Species_Library.json` to drop the vestigial `KA` field.
