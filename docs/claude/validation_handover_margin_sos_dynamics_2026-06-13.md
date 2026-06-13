# Margin validation — session handover (2026-06-13, session 2: SoS + dynamics)

*Cross-session handoff. Continues `validation_handover_2026-06-13.md` (session 1, which pivoted
the programme to validating the **adaptive margin** and ran the first SFG test). This session
**extended the margin validation across four organisational levels** (state → function →
controlled mechanics → dynamics), wrote a consolidated manuscript-ready synthesis, and ended on
**one open design decision** (a circularity-safe *powered* dynamic test). Self-contained: nothing
below should need re-deriving. Everything described is committed and pushed to `main`
(through merge `97fcc34`).*

---

## Part 0 — TL;DR (read first)

1. **The margin now has external support at four levels, each honestly scoped** (full account in
   `docs/notes/external_validation_synthesis.md`, the new authoritative write-up):
   - **state** — Scope for Growth (3 studies): **+0.41** (Widdows 1995, estuary) → **+0.12**
     (Widdows 2002, basin) → **−0.11** (Albentosa 2012, condition-confounded). Scale-dependent.
   - **function** — Stress-on-Stress (DOME, 17 UK stations): **+0.39 → +0.45** under body-size +
     condition control (controlling the confound *strengthens* it — the opposite of Albentosa).
   - **controlled mechanics** — Viarengo 1995: monotone dose-response (potency Cu>PAH>PCB) + a
     mixture that is **additive, no antagonism** (the model's TU/IA rules bracket it).
   - **dynamics** — Veldhuizen 1991 transplant: the model's erosion *dynamics* reproduce
     continued erosion the static map can't; firmed up single-contaminant (Cd-alone, **+0.90**).
2. **Crucial honesty correction made this session:** the SFG/SoS/GlobTherm tests are **static
   point-API maps** (`compute_adaptive_margin_response`), **not** simulations of the dynamics
   (`B_t` accumulation, the erosion ODE). The wiki + notes were rewritten to say so; "amplification"
   was demoted to "cross-sectional shadow." Only the transplant test (and its firm-up) touches the
   actual dynamics.
3. **The amplification scalar `g`/`F` is null in every test** — the margin-first prediction,
   confirmed again.
4. **THE OPEN DECISION (where we stopped):** a *powered* dynamic test is hard because off-the-shelf
   mussel-SoS data lacks dense post-plateau sampling, and the obvious powered route (DEBtox/DEB-TKTD
   datasets) is **circular** (DEB validating DEB). The proposed circularity-safe + powered design is
   the **toxicodynamic twin of the COMADRE test: does AmP `k_M` predict the GUTS dominant rate
   constant `k_D`** (fitted to *raw* survival, non-DEB), **controlling body size + phylogeny**? The
   user was asked whether to **scope the `k_D`-database ↔ AmP-species overlap** — *no answer yet*.
   See Part 5.

---

## Part 1 — What this session added (chronological, all committed)

Building on session 1 (which delivered Widdows 1995 SFG = +0.41):

1. **SFG #2 — Albentosa 2012** (*STOTEN* 435–436:430), *M. galloprovincialis*, Iberian SMP, 39
   site×survey. **margin~SFG = −0.11** (bounding negative). The authors themselves show SFG is
   dominated by **condition index** (r=−0.62***) and **age**; metals (As/Zn/Cd) are a *positive*
   confound; controlling CI/age does **not** rescue it. *Tissue-burden ≠ exposure* biting hard.
   (Also corrected a misattribution: session 1 called this "Beiras" — it is **Albentosa**.)
   Files: `data/external/sfg_albentosa2012_{iberia,contaminants,biometric_repeated}.csv`,
   `examples/sfg_margin_validation_albentosa2012.jl`.

2. **ICES DOME probe — is there an SFG *database*?** **Dead end.** SFG is a recognised ICES method
   (TIMES 40) but has been **dropped from active OSPAR/DOME monitoring** (checked the full 2024
   OSPAR CEMP bundle, figshare 27211422 — no `SFG` determinand). Capacity-route scouting: non-mussel
   SFG is **temperature-** not contaminant-driven → the chemical capacity test stays data-blocked.

3. **SFG #3 — Widdows 2002** (*Mar. Environ. Res.* 53:327), *M. edulis*, **Irish Sea, 23 sites**.
   **margin~SFG = +0.12** — weak/partial, *as the paper predicts* (single-marker correlations
   collapse at basin scale). PAH/assimilation axis −0.27, metals +0.25 (confound again). **SFG was
   figure-digitised from Fig 2A** (no numeric table) → rank-only. Files:
   `data/external/sfg_widdows2002_{irishsea,contaminants}.csv`,
   `examples/sfg_margin_validation_widdows2002.jl`.

4. **`docs/wiki/Margin-Validation.md`** — new wiki page consolidating SFG + SoS + GlobTherm; linked
   from `Home.md` (validation table, status row, bottom line).

5. **SoS (Stress-on-Stress) — DOME static** *(the headline new result)*. Survival-in-air (days)
   under emersion is the closest *outcome* to the margin's purpose (resilience to an acute hit).
   ICES DOME 2024 OSPAR CEMP, *M. edulis*, **17 UK stations, 2012–2022**, co-located contaminants +
   body size — the multi-station, exposure-paired, QA'd, open dataset SFG lacked.
   **survival~margin = +0.39 → +0.45** under length + condition control. **The decisive contrast
   with Albentosa:** here controlling the confound *strengthens* the signal (genuine erosion, not a
   health-proxy artifact); metals near-dead (`p_maint` +0.09); PAH −0.43, PCB −0.48*. Files:
   `scripts/extract_dome_sos.jl`, `data/external/sos_dome_ukcemp.csv`,
   `examples/sos_margin_validation_dome.jl`.

6. **SoS temporal (within-station).** Station-year panel (nearest-year ±2 matching). The clean
   within-station fixed-effects test is **+0.15 (n.s.)** — erosion *over time* underpowered (thin
   panel, stable within-station burden). Pooled +0.28 (PAH −0.33*, PCB −0.34*); station-level
   QC-cleaned **+0.62**. Files: `scripts/extract_dome_sos_yearly.jl`,
   `data/external/sos_dome_ukcemp_yearly.csv`, `examples/sos_temporal_validation_dome.jl`.

7. **Honest-framing correction (important).** Realised + documented that all of §§1–6 use the
   **static point API** — they validate the response curve + MoA routing + capacity weighting, NOT
   the dynamics. Added a "what these tests exercise / don't" scope box to the wiki; softened the
   "amplification" language to "static margin↔acute-resilience map." (This was a user-prompted
   correction — keep this honesty.)

8. **Dynamic transplant proof-of-concept — Veldhuizen 1991** (*ACET* 21:497–504). The **first** test
   that runs the **dynamics** (`simulate_deb_axis_response`). Clean *M. edulis* transplanted to a
   Western Scheldt gradient, stress indices at **2.5 & 5 months**. *Discriminating feature:* **Cd
   plateaus by 2.5 mo, yet SoS keeps dropping to 5 mo** → a static map predicts no further erosion;
   the dynamic erosion state rises **+33%** (still rising because `1/λ = 1/k_M ≈ months`, **unfitted**)
   and matches; ρ(erosion,survival)=−1 (n=4). Files: `data/external/sos_veldhuizen1991_transplant.csv`,
   `examples/sos_dynamic_validation_veldhuizen.jl`.

9. **`docs/notes/external_validation_synthesis.md`** — consolidated, manuscript-ready account of
   **everything externally validated** (scorecard, per-anchor §§3–7, the `g`/`F` null, cross-cutting
   caveats, validated-vs-open, per-anchor reproducibility, references). Linked from `Home.md`.
   **This is the document to update as the single source of truth going forward.**

10. **Viarengo 1995 — controlled validation** (*Mar. Environ. Res.* 39:245). First **controlled-
    exposure** test of the **impairment curve + mixture model** (field tests can't isolate single
    contaminants/mixtures). *M. galloprovincialis*, 3-day exposure, SoS LT50.
    **(A)** monotone dose-response on each MoA axis, potency **Cu>PAH>PCB**.
    **(B)** Cu+DMBA mixture worse than either component (**no antagonism**); the model's *own* rules
    (`aggregate_axis_mixture_effects`: TU/CA & IA) **bracket** observed LT50 (CA 5.25/3.89, IA
    5.14/3.57 vs obs 5.0/3.0) → corroborates "mixtures are additive assumptions, not fitted
    interactions." Files: `data/external/sos_viarengo1995_doseresponse.csv`,
    `examples/sos_mixture_validation_viarengo.jl`.

11. **Dynamic firm-up — Veldhuizen 1991** (*ACET* **20**:259–265, single-contaminant). Removes the
    Cd-vs-PCB confound of the transplant: Cd **alone** erodes SoS progressively (lab LT50
    10.7→9.5→7.6; semi-field 9.3→8.6), modelled margin tracks it **ρ=+0.90 (n=5)**; PCB **alone**
    erodes SoS with a **delayed onset**. So the transplant's continued erosion is **not** a PCB
    artifact. Does **not** add a clean *constant-burden* test (burden rises throughout; near-plateau
    10-mo LT50 figure-only). Files: `data/external/sos_veldhuizen1991_singlecontaminant.csv`,
    `examples/sos_dynamic_firmup_veldhuizen_singlecontaminant.jl`.

---

## Part 2 — Honest scorecard (current)

| anchor | level | result | status |
| --- | --- | --- | --- |
| COMADRE (`k_M`, `R_i`) | population demography | `k_M`↔recovery **+0.19\***, `R_i`↔compensation **+0.77\*\*** | ✅ rate endpoints |
| GlobTherm (n=664) | individual physiology | capacity coherent (\|ρ\|≤0.45) but **recovery-specific** | ✅ bounding |
| SFG ×3 | individual energetics (no scale bridge) | **+0.41 → +0.12 → −0.11** (scale-dependent) | ✅ / ◐ |
| SoS DOME (static) | individual energetics | **+0.39 → +0.45** (confound-controlled), +0.62 QC | ✅ static map |
| Viarengo 1995 (controlled) | individual, controlled dose | dose-response + additive mixture, no antagonism | ✅ controlled |
| Dynamics (transplant + single-contaminant) | individual, *over time* | continued erosion static can't; Cd-alone **+0.90** (de-confounded) | ◑ proof-of-concept |
| amplification scalar `g`/`F` | — | **null everywhere** | ✅ (margin-first) |

**Two recurring threads to remember:** (a) **metals are a positive confound** (As/Cd/Zn) across
SFG and SoS — the MoA routing exists to keep them off the toxic (hydrocarbon→assimilation) axis,
which is why the routed margin beats naive load; (b) results are **rank-robust, magnitude-modest,
specification-sensitive** (the programme's signature, cf. COMADRE `k_M`).

---

## Part 3 — Model facts to carry (don't re-derive)

- **Julia 1.12.6 via `julia +release`** (default LTS cannot load the project — `StaticData`
  precompile error = wrong Julia).
- **Static API:** `compute_adaptive_margin_response(axis_pressures, params)` →
  `E=x/(1+x)` per axis → `Q=ΣwᵢEᵢ` → `A_t=A0·(1−Q)`. **This is what every SFG/SoS correlation uses.**
- **Dynamic API:** `simulate_deb_axis_response(t, pulse_axes, A0, params, q; y0, dt)` integrates the
  slow erosion state `y` (rises toward `q·cost/λ` with time-constant `1/λ`). Used only by the
  transplant harness. `λ = restoring_force_from_margin(A0−cost)`.
- **MoA routing (used throughout):** hydrocarbons/PAH→assimilation; metals→maintenance;
  organotin+organochlorine(PCB,DDT)→reproduction; DBT→growth. Pressure = **median-normalised
  relative tissue burden** (`tu = conc/median`), threshold-free.
- **M. edulis params** (from `data/AmP_Species_Library.json`): `A0=181.43`, `λ_min=k_M=0.001127/day`,
  `λ_max=0.014637/day`, `g=12.99`, `alpha_axes=[0.0055, 0.240, 0.997, 0.0035]`. The erosion timescale
  `1/λ ≈ 68–887 days ≈ months` — this is why the dynamic transplant test works, and it is **unfitted**.
- **Mixture rules** (`src/mixture_aggregation.jl`, `aggregate_axis_mixture_effects`):
  `axis_toxic_unit_sum` (CA: sum x then `X/(1+X)`), `independent_action_axis_effects` (IA:
  `1−∏(1−Eᵢ)`), `grouped_ca_then_ia_axis_effects`. Returns `E_assimilation/maintenance/growth/reproduction`.
- **DEBtox scaled damage `D_t` is still UNIMPLEMENTED** (CLAUDE.md invariant). The dynamic erosion is
  the `simulate_deb_axis_response` `y`-state, *not* a DEBtox damage variable.
- Raw DOME download `data/external/ospar_cemp2024_*.csv` is **gitignored** (re-fetch via the extract
  scripts' header URLs); derived CSVs are committed (the repo's standard pattern).

---

## Part 4 — The circularity principle (carry this — it frames the open decision)

A DEB-derived model is **not** inherently tautological. Circularity arises **only** when a
DEB/DEBtox quantity sits on the **outcome** side, or the same parameters appear on both sides.
- ✅ **Non-circular (what we did):** raw measured outcome (SoS survival, mortality, condition,
  energy reserves) + measured tissue burden + `k_M` from **clean-condition AmP life-history** (never
  fit to the toxicity data). COMADRE is the template: DEB `k_M` predicted an independent demographic
  recovery rate, controlling pace.
- ❌ **Circular (avoid):** a **DEBtox/DEB-TKTD–fitted** curve used as "data"; or organism parameters
  drawn from AmP on *both* sides. **Rule: never let a DEB-derived quantity touch the outcome side.**

This is why the DEBtox/DEB-TKTD route for the powered dynamic test was **rejected**.

---

## Part 5 — THE OPEN DECISION: a circularity-safe *powered* dynamic test

**Problem.** A clean powered dynamic test needs dense *post-plateau* sampling, which off-the-shelf
mussel-SoS data lacks; the only powered datasets with that design are DEBtox/DEB-TKTD — which are
circular (Part 4).

**Proposed design (the toxicodynamic twin of COMADRE).** The dynamic claim reduces to one statement:
*the margin's erosion/recovery timescale ≈ the maintenance timescale `1/k_M`* (the model's erosion
state relaxes at `λ ∈ [k_M, k_M·g]`). There is an independent, **non-DEB** measurement of exactly
that: **GUTS's dominant rate constant `k_D`** — the rate at which the toxicodynamic damage state
follows exposure, **fitted to raw survival time-courses** by a phenomenological (non-DEB) model,
available for **hundreds of species×chemical calibrations** (openGUTS, EFSA TKTD case studies,
published `k_D` compilations).

> **Test:** across species, does **AmP `k_M`** (clean life-history) predict **GUTS `k_D`** (raw
> survival), **controlling body size + phylogeny**?

- **Non-circular:** `k_D` from raw survival via a non-DEB model; `k_M` from independent life-history.
- **Powered:** far more `k_D` values than we ever had SoS sites.
- **On-thesis:** tests the exact dynamic prediction (erosion timescale = maintenance timescale),
  out-of-sample.
- **The one real caveat:** both `k_M` and `k_D` scale allometrically with size, so — exactly as
  COMADRE controlled pace-of-life — the test must ask whether `k_M` predicts `k_D` **beyond body
  size** (the residual-after-size signal is the non-trivial result). Phylogeny control too.

**Fallback if `k_D`↔AmP overlap is thin:** extract the **observed** erosion timescale (time-to-plateau
of a raw effect time-course) for a handful of AmP species, compare to their independent `1/k_M`.
Same non-circular logic, less power.

**Immediate next step (was offered to the user, not yet answered):** *scope feasibility* — find the
best public `k_D` compilation (openGUTS / EFSA 2018 TKTD opinion datasets / a `k_D` meta-analysis)
and check how many of those species are in AmP (so `k_M` is available). That decides whether the
powered version is on the table. **Start here next session.**

Useful leads already found: GUTS / DEB-TKTD reviews — Jager et al. (GUTS); EFSA 2018 TKTD opinion
(*EFSA Journal* 16:5377, case-study datasets); openGUTS project (software + ring-test datasets);
EST 2021 "Predicting Mixture Effects over Time with GUTS" (PMC7893709).

---

## Part 6 — Other open items (lower priority than Part 5)

1. **Across-species capacity weighting** (the model's distinctive content) — needs across-species
   contaminant-gradient data that **largely does not exist** (SFG corpus mussel-dominated; DOME
   dropped SFG; non-mussel SFG is temperature-driven). Likely needs purpose-designed data.
2. **Dated-tree PGLS** (real-phylogeny control for the COMADRE `k_M` result) — pipeline built &
   smoke-tested (`scripts/comadre_pgls_dated.jl`); needs **one** VertLife/TimeTree dated-Newick
   download into `data/external/comadre_amp_dated_tree.nwk`, then runs unattended.
3. **GitHub Wiki-tab mirror** — `docs/wiki/` (now incl. `Margin-Validation.md`) is the source of
   truth; the Wiki *tab* is a separate manual re-sync (clone `*.wiki.git`, copy, rewrite
   `[[links]]`/`../`, push to `master`).
4. **(low) Regenerate `AmP_Species_Library.json`** to drop the vestigial `KA` field (session-1 item;
   non-deterministic JSON key order makes a clean regen worthwhile but cosmetic).

---

## Part 7 — Files, commands, reproducibility

**New data (this session)** — all in `data/external/`:
`sfg_albentosa2012_*.csv`, `sfg_widdows2002_*.csv`, `sos_dome_ukcemp.csv`,
`sos_dome_ukcemp_yearly.csv`, `sos_viarengo1995_doseresponse.csv`,
`sos_veldhuizen1991_transplant.csv`, `sos_veldhuizen1991_singlecontaminant.csv`.
(Raw DOME `ospar_cemp2024_*.csv` gitignored — re-fetch via `scripts/extract_dome_sos*.jl` headers.)

**New harnesses** — `examples/`:
`sfg_margin_validation_albentosa2012.jl`, `sfg_margin_validation_widdows2002.jl`,
`sos_margin_validation_dome.jl`, `sos_temporal_validation_dome.jl`,
`sos_mixture_validation_viarengo.jl`, `sos_dynamic_validation_veldhuizen.jl`,
`sos_dynamic_firmup_veldhuizen_singlecontaminant.jl`.
**Extractors** — `scripts/extract_dome_sos.jl`, `scripts/extract_dome_sos_yearly.jl` (base Julia,
read the gitignored raw DOME csv).

**Docs** — `docs/notes/external_validation_synthesis.md` (the synthesis — keep current),
`docs/notes/sos_validation_status.md` (SoS detail), `docs/notes/sfg_validation_status.md`,
`docs/notes/globtherm_validation.md`, `docs/wiki/Margin-Validation.md`, `docs/wiki/Home.md`.

**Run any harness:** `julia +release --project=. examples/<name>.jl`
**Run an extractor:** `julia +release scripts/extract_dome_sos.jl` (no `--project`; base Julia).

**Source PDFs supplied by the user** (in `~/Downloads`, NOT committed — paywalled): `albentosa2012.pdf`,
`widdows2002.pdf`, `viarengo1995.pdf`, `BF01183870.pdf` (Veldhuizen transplant, ACET 21),
`veldhuizen-tsoerkan1991 (1).pdf` (Veldhuizen single-contaminant, ACET 20).

---

## Part 8 — Git state

All work is **committed and pushed to `main`** (https://github.com/microstijn/burdenDK). Key recent
commits (newest first):
- `97fcc34` Merge — Viarengo controlled mixture validation + single-contaminant dynamic firm-up
- `5b1c235` Merge — SoS temporal + dynamics proof-of-concept + honest framing + validation synthesis
- `c3483c1` Merge — DOME dead end + Widdows 2002 (partial replication)
- `ba0749a` Merge — Margin-Validation wiki page + Stress-on-Stress amplification test
- `c6be02e` Merge — SFG margin test #2 (Albentosa, bounding negative)

Working tree clean; no unpushed branches. (Note: pushing to `main` requires explicit user
authorization each time — the auto-classifier blocks autonomous pushes.)

---

## Part 9 — Prioritised next actions

1. **Scope the powered dynamic test (Part 5)** — find a public GUTS `k_D` compilation and measure its
   species overlap with AmP. If viable: build `k_M` (AmP) vs `k_D` (GUTS), controlling size +
   phylogeny — the non-circular, powered, on-thesis dynamic test. *This is the live thread.*
2. **Dated-tree PGLS (Part 6.2)** — one VertLife/TimeTree download → runs unattended.
3. **Mirror the wiki tab (Part 6.3)** — now that `Margin-Validation.md` exists.
4. Keep `external_validation_synthesis.md` as the single source of truth; update it as results land.
