# pMoA stressor-routing table — basis and rationale

*Companion to `data/pMoA_Stressor_Routing.csv`. This routes **aggregate water-quality / model
stressors** (microplastics, BOD, salinity, pathogens, …) onto the four DEB margin axes
(assimilation, maintenance, growth, reproduction) by **physiological mode of action (pMoA)**, not by
fitted weights — so the mapping is a **declared mechanistic assumption** (like the mixture rules),
not a tuning knob. Drafted 2026-06-13; per-stressor citations are slots to fill/verify, not yet
audited.*

## Why pMoA, not weights

The earlier DynQual prototype mapped each proxy to axes with free `w_A/w_M/w_G/w_R` weight vectors —
i.e. tuning knobs, which the project's no-arbitrary-knobs invariant forbids. The pMoA framework
removes the knob: each stressor is assigned to **the metabolic process it physiologically taxes**,
stated as a falsifiable mechanism. The residual judgment (the assignment itself) is then
**sensitivity-tested**, not fitted (see §4).

## 1. The three pillars the assignments rest on

1. **DEBtox / DEB-TKTD pMoA taxonomy** (Kooijman; Jager et al.). Every stressor is assumed to act on
   *one primary* process: assimilation, (somatic) maintenance, growth cost, reproduction cost, or
   hazard. These are the model's four axes (+ hazard = the acute survival / Stress-on-Stress endpoint,
   handled outside the margin). This is the theoretical skeleton — established, not invented.
2. **Internal consistency with the existing CAS-MoA routing** (`src/mode_of_action.jl`,
   `src/moa_deb_mapping.jl`): hydrocarbons/PAH→assimilation; metals→maintenance;
   organotin+organochlorine→reproduction; DBT→growth. The chemical rows in the table **inherit these
   axes unchanged**; aggregate stressors slot in by the same logic. The table extends the existing
   routing, never contradicts it.
3. **Stressor-class ecophysiology** — the documented dominant sublethal mechanism per stressor
   (carried per-row with a `confidence` grade and a citation slot to fill).

## 2. Schema (`data/pMoA_Stressor_Routing.csv`)

| column | meaning |
| --- | --- |
| `primary_axis` | the committed pMoA assignment (the conservative default uses primary only) |
| `secondary_axis`, `axis_split` | optional coarse mechanistic prior (e.g. `0.7/0.3`) — **not fitted**; sensitivity-test it, or drop to primary-only |
| `tissue_retention_rho` | monthly tissue-memory ρ — **only for accumulating stressors**; ambient conditions get 0 (see §3) |
| `ambient_condition` | `yes` = an environmental state (salinity, pH, hypoxia), not a tissue burden |
| `pmoa_basis` | which pillar/category the routing rests on |
| `confidence` | high / medium / low (evidence strength for the assignment) |
| `mechanism_notes` | the one-line falsifiable mechanism + consistency note |

## 3. Two structural rules encoded in the table

- **Temperature is NOT a pressure axis.** It is the Arrhenius clock on *all* DEB rates — including the
  erosion/recovery rate `λ` and `k_M`, hence the whole two-timescale logic. Route it through the DEB
  temperature correction, never as a burden (`primary_axis = NONE`, `pmoa_basis = rate-modifier`).
  Treating temperature as a maintenance pressure double-counts and corrupts the timescale.
- **`ρ` (tissue memory) ≠ margin erosion.** Tissue memory applies only to stressors with a tissue
  reservoir (microplastics, metals, lipophilic organics → high ρ). *Ambient* stressors (salinity, pH,
  hypoxia, ammonia) get ρ=0 — the current condition is the burden — yet a **sustained** ambient
  stressor still erodes the margin over months via the slow `λ`-dynamics, not via accumulation. The
  two memory layers stay distinct, per the invariant.

## 4. How to defend it without external validation of the mapping

You do not (yet) have data validating "BOD→maintenance" directly. The honest defensibility claim is:
- **engine** (saturation `E=x/(1+x)`, margin, recovery dynamics) — externally validated (SFG/SoS/COMADRE);
- **routing** — declared mechanism (pMoA), literature-anchored, in the same epistemic class as the
  mixture rules;
- **robustness** — **sensitivity-test the assignment**: perturb each stressor's axis (and the coarse
  splits / ρ) across its plausible alternatives and show the **site/month vulnerability rankings are
  stable**. Given the model is rank-robust everywhere else, they almost certainly are — and that
  stability *is* the defensibility statement for the water-quality use, standing in for the external
  validation you cannot get.

## 5. Per-row notes worth flagging

- **microplastics** — primary assimilation (feeding impairment) is the best-supported sublethal MoA;
  secondary maintenance (inflammation/oxidative) is a coarse prior. High ρ (particles retained).
- **nutrients/eutrophication** — flagged `indirect`, `low`: not a direct toxicant; it acts through
  secondary hypoxia and food-web shifts. Consider modelling it *via* its BOD/hypoxia consequence
  rather than as a primary stressor.
- **pathogens** — `low` confidence; the immune-energetics → axis mapping is the least-settled row.
- **AChE pesticides** — chronic route maintenance; **acute route is hazard** (the GUTS survival
  endpoint), which is exactly the ECOTOX-tested class in synthesis §7b (where single-trait `k_M`→LC50
  was size-confounded — a caution that this row's *cross-species* strength is bounded).
- **pharmaceuticals_generic** — a placeholder default only; must be reassigned per known API MoA.

## 6. Open / to fill
- Replace the citation slots with verified references per row (I have not audited these).
- Decide the conservative default: primary-axis-only (fewest assumptions) vs. primary+secondary split.
- Wire the table into the pipeline as the single source of stressor routing (replacing inline weights),
  mirroring how `Compound_Memory_Library.csv` is loaded.
