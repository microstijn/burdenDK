# External Validation of the Adaptive Margin — Synthesis (2026-06-13)

*A consolidated, manuscript-ready account of every external validation run against the
TwoTimescaleResilience / burdenDK framework, with the honest scope of each. Self-contained;
per-anchor detail lives in the cited `docs/notes/*` and `docs/wiki/*` files.*

---

## 1. What is being validated, and the strategy

The framework models **background-conditioned vulnerability**: chronic environmental pressure
and retained chemical burden slowly narrow a species' **adaptive margin** `A_t`, which weakens
its **restoring force** `λ(A_t)`, which **amplifies** the burden of a later acute perturbation.
Capacity is derived offline from Add-my-Pet (AmP/DEB) parameters; pressure from EPA ECOTOX;
memory from a compound-retention recurrence. The recovery curve is linear in relative margin,
`λ(A) = λ_min + (λ_max−λ_min)·clamp(A/A0,0,1)`, with `λ_min` anchored to the DEB maintenance
rate constant `k_M = [p_M]/[E_G]` and `λ_max/λ_min = g` (the energy-investment ratio).

**The margin-first reframe (central to reading the results).** The validated *product* is the
**margin / recovery layer** — its rate endpoints, state, function, and (now, tentatively)
dynamics. The one-dimensional **amplification scalar** `g`/`F` is a derived diagnostic and is
**null in every external test**; that null is itself a prediction the data confirm.

**Strategy.** External anchors at *multiple organisational levels* (population demography →
individual energetics); **rank statistics throughout** (unit-invariant, robust); **pre-registered
signs**; and honest reporting of nulls and bounds. Two anchors require an individual→population
**scale bridge** (COMADRE) and are argued, not 1:1; the SFG/SoS anchors are at the margin's *own*
organisational level (**no scale bridge**).

---

## 2. The scorecard

| # | anchor | organisational level | what it tests | headline result | status |
| --- | --- | --- | --- | --- | --- |
| 1 | **COMADRE** (matrix demography) | population | recovery **rate endpoints** `k_M`, `R_i` | `k_M`↔recovery **+0.19\*** (beyond pace+phylo); `R_i`↔compensation **+0.77\*\*** | ✅ corroborated |
| 2 | **GlobTherm** (thermal tolerance) | individual physiology | is the capacity axis *general* resilience? | coherent (\|ρ\|≤0.45) but **recovery-specific** (general-resilience refuted) | ✅ bounding |
| 3 | **Scope for Growth** (3 studies) | individual energetics | the margin **state** under pressure | **+0.41** (estuary) → +0.12 (basin) → −0.11 (confounded) | ✅ / ◐ scale-dependent |
| 4 | **Stress-on-Stress** (ICES DOME) | individual energetics | the margin **function** (acute resilience), static | **+0.39 → +0.45** (confound-controlled), +0.62 QC-cleaned | ✅ static map |
| 4b | **Viarengo 1995** (controlled SoS) | individual, controlled dose | the **impairment curve + mixture model** | monotone dose-response (potency Cu>PAH>PCB); mixture additive, no antagonism (TU/IA bracket) | ✅ controlled |
| 5 | **Transplant + single-contaminant time-courses** (Veldhuizen 1991) | individual, *over time* | the margin **dynamics** (sustained-burden erosion) | dynamics reproduce continued erosion the static map can't (n=4); Cd-alone erodes SoS progressively, ρ(margin,LT50)=+0.90 (PCB confound removed) | ◑ proof-of-concept, de-confounded |
| — | **amplification scalar `g`/`F`** | — | the 1-D readout | **null everywhere** (−0.05…−0.13) | ✅ (margin-first prediction) |

*\* p<0.05, ** p<0.01, partial/within-Order where noted. Magnitudes are modest and rank-based
throughout (see §9).*

---

## 3. Recovery rate endpoints — COMADRE (population demography)

The first external anchor: do the model's recovery quantities predict independent demographic
recovery from the [COMADRE](https://compadre-db.org) matrix database, **beyond pace-of-life and
beyond coarse phylogeny**? Recovery = log₁₀ damping ratio `|λ₁|/|λ₂|` per species (n=197 matched
via a GBIF-backbone name resolver). Three nested controls: raw, partial on Caswell generation
time, and within-Order group-mean-centred (tree-free phylogeny proxy).

| model quantity | raw ρ | \|gen | \|gen+Order |
| --- | --- | --- | --- |
| `λ(A0)` recovery rate | +0.362\*\* | +0.173\* | +0.089 n.s. |
| **`λ_min = k_M`** | +0.406\*\* | +0.264\*\* | **+0.190\*** |
| `g` (amplification) | −0.109 | −0.128 | −0.055 n.s. |

- **`k_M` predicts demographic recovery beyond pace *and* coarse phylogeny.** `λ(A0)`'s broader
  signal is largely between-clade (n.s. under Order control). `g` is null at every level.
- **Per-axis (Idea B) — the strongest single result:** the DEB **reproduction rate `R_i`**
  specifically predicts demographic **compensation** ρ = **+0.77\*\*** (beyond pace + size).
- **Maturation timing (`a_p`) resolved:** the pre-registered *negative* `a_p`→compensation came
  out *positive* because `a_p` is pace-loaded (ρ(`a_p`,gen)=+0.50); the residual-after-pace signal
  is genuinely positive, is not a fecundity proxy (ρ(`a_p`,`R_i`)=−0.13) and survives a
  matrix-dimension control — a reproduction-*timing* axis.
- **Robustness:** bootstrap (resample-over-species) `k_M`→recovery ρ=+0.264, **95% CI [+0.14,
  +0.38]**, survives Benjamini-Hochberg over 7 tests; `g` is the **only** test whose CI spans 0.
  Partial ρ holds 0.18–0.33 across 6 COMADRE matrix-quality filters (not a filter artifact).
- **Specification-sensitive (important):** the gen-controlled `k_M` signal is **rank-based**
  (partial Spearman +0.264); the log-linear partial is only +0.04 — monotone but not log-linear.
- **Open:** a *real-phylogeny* control. The OTL+Grafen tree gives Pagel λ≈0 (too weak to test);
  a **dated-tree PGLS** pipeline is built and smoke-tested, awaiting one VertLife/TimeTree download.

*Detail:* `docs/wiki/COMADRE-External-Validation.md`, `Per-Axis-Resilience.md`,
`Phylogenetic-PGLS.md`; notes `comadre_*`.

---

## 4. Capacity coherence and its bound — GlobTherm

A pre-registered probe (n=664 AmP↔GlobTherm species): does the AmP recovery-capacity axis carry
a *general* thermal-tolerance signal (broader `CTmax−CTmin`), beyond body size and latitude? Two
honest, opposite-pointing findings:
- ✅ **Coherence:** the AmP capacity axis correlates strongly (|ρ| up to **0.45**) with an
  entirely independent physiological dataset — the offline AmP→capacity extraction carries real
  biological structure, not noise.
- ❌ **General-resilience refuted:** higher `k_M`/`λ_max` → *narrower* thermal breadth (−0.45
  all-taxa; weak −0.28\* in the cleaner ectotherm subset). Recovery capacity and thermal tolerance
  are **separate axes** — the COMADRE-validated recovery capacity is *specific to demographic
  recovery*, **not** a universal resilience currency.

This bounds the claim (do not sell capacity as "general resilience") and motivates testing the
margin's *own* currency — energetics under contaminant pressure. *Detail:* `globtherm_validation.md`.

---

## 5. The margin **state** — Scope for Growth (individual energetics, no scale bridge)

SFG = energy absorbed − energy respired ≈ capacity beyond maintenance ≈ the adaptive margin in
energetic terms, measured at the **same organisational level** (no scale bridge) and independent
of AmP. Method: per-site tissue burden → mode-of-action axis routing → median-normalised relative
burden → margin-first point API → rank-correlate modelled margin `A_t` vs measured SFG.

| study | gradient | n | ρ(margin, SFG) |
| --- | --- | --- | --- |
| **Widdows et al. 1995** (North Sea, *M. edulis*) | estuary/regional, hydrocarbon | 36 | **+0.41\*** (beats naive 0.22; ≈ best single 0.47) |
| **Widdows et al. 2002** (Irish Sea, *M. edulis*) | basin-scale, hydrocarbon | 23 | **+0.12** (PAH axis −0.27; beats naive 0.005) |
| **Albentosa et al. 2012** (Iberia, *M. galloprovincialis*) | condition/food-confounded | 39 | **−0.11** (the authors: SFG~CI −0.62\*\*\*) |

**A coherent, scale-dependent picture.** The margin tracks SFG **strongly where tissue burden
indexes exposure** (Widdows 1995), **weakly at basin scale** where burden decouples from exposure
(Widdows 2002 — the paper's own point), and is **confound-flipped where condition/food dominates**
(Albentosa — controlling CI/age does not rescue it; a bounding result, not a tuning failure).
**Across all three, metals behave as a positive confound** (As/Cd/Zn correlate *positively* with
SFG); the mode-of-action routing's job is to keep them off the toxic (hydrocarbon→assimilation)
axis, which is why the routed margin consistently beats naive equal-weight load.
*Detail:* `sfg_validation_status.md`. *(Widdows 2002 SFG was figure-digitised → rank-only.)*

---

## 6. The margin **function** — Stress-on-Stress, static map (ICES DOME)

Where SFG is the margin *state*, **stress-on-stress** (survival-in-air, days, under emersion/
anoxia) is the closest *outcome* to the margin's purpose: resilience to an acute hit. Data: ICES
DOME 2024 OSPAR CEMP (open, CC BY 4.0), *M. edulis*, **17 UK stations**, 2012–2022, co-located
contaminants + body size — the multi-station, exposure-paired, QA'd dataset SFG lacked.
**Pre-registered positive prediction.**

| test | ρ (n=17) |
| --- | --- |
| survival ~ modelled margin `A_t` | **+0.39** |
| survival ~ margin \| body length | +0.40 |
| **survival ~ margin \| length + condition** | **+0.45** |

Axis diagnostic: PAH/assimilation **−0.43**, metals/maintenance +0.09 (confound near-dead here),
PCB/reproduction **−0.48\***; routed margin (0.39) beats naive load (0.32). **The decisive contrast
with Albentosa:** there, controlling condition could not rescue the margin (the confound *was* the
signal); here, partialling size + condition **strengthens** the margin signal (0.39→0.45) — genuine
margin erosion, not a health-proxy artifact. QC-cleaned/nearest-year aggregation gives +0.62.

> **Scope (applies to §§3–6).** These are **static, cross-sectional maps**: the *already-
> accumulated* burden is passed through the instantaneous point API (`compute_adaptive_margin_
> response`: burden → bounded impairment `E=x/(1+x)` → `Q` → `A_t=A0·(1−Q)`) and `A_t` is
> correlated with the outcome. They validate the **response-curve shape, MoA routing, and AmP
> capacity weighting** — not the time dynamics (`B_t` accumulation, the erosion ODE), which are §7.

**Honest power:** n=17 → two-sided n.s. (p≈0.1); under the pre-registered one-sided prediction the
controlled result is marginally significant (p≈0.04). A within-station temporal pass (station-year
panel) is **positive but underpowered** (+0.15 n.s.; pooled +0.28 with PAH −0.33\*, PCB −0.34\*).
*Detail:* `sos_validation_status.md`.

### 6b. Controlled cross-check — impairment curve + mixture model (Viarengo 1995)
The only **controlled-exposure** validation: *M. galloprovincialis*, 3-day single-contaminant and
binary-mixture exposures, survival-in-air LT50 (Viarengo et al. 1995). **(A)** LT50 falls
monotonically with dose for Cu (maintenance), DMBA/PAH (assimilation), Aroclor/PCB (reproduction),
ρ=−1 each, potency **Cu>PAH>PCB** — consistent with the saturating per-axis impairment `E=x/(1+x)`.
**(B)** the Cu+DMBA mixture is worse than either component (real effect, **no antagonism**) and the
model's own mixture rules (`aggregate_axis_mixture_effects`, CA/TU & IA) **bracket** the observed
LT50 (CA 5.25/3.89, IA 5.14/3.57 vs obs 5.0/3.0) — corroborating the **"mixtures are additive
assumptions, not fitted interactions"** invariant; mild supra-additive excess unresolved (n=2,
LT50 rounding). Static (single timepoint). *Detail:* `sos_validation_status.md`.

---

## 7. The margin **dynamics** — transplant proof-of-concept (the only test that runs the dynamics)

The **first** test that exercises the erosion *dynamics* (`simulate_deb_axis_response`), not the
static map. Data: **Veldhuizen-Tsoerkan et al. 1991** (ACET 21:497–504) — clean *M. edulis*
transplanted to a Western Scheldt contamination gradient (4 sites), stress indices at **2.5 and
5 months**.

**The discriminating feature:** Cd accumulates fast and **plateaus by 2.5 mo**, yet SoS survival
keeps **dropping** 2.5→5 mo at contaminated sites. A static burden→margin map predicts ~no further
erosion after plateau; a dynamic model integrating erosion under *sustained* burden predicts
continued erosion **iff** its timescale `1/λ` ≈ months. For *M. edulis*, `λ_min=k_M=0.00113/day ⇒
1/λ ≈ 68–887 d ≈ months`, **unfitted** — it falls out of the maintenance rate. Result: the dynamic
erosion state rises **~33% from 2.5→5 mo** (matching the continued SoS decline); the static map
gives ~0 extra erosion and **cannot** reproduce it. Cross-sectionally ρ(erosion, survival)=−1 (n=4).

**Honest scope:** n=4 sites × 2 times, figure-digitised; ρ=−1 trivial at n=4; the model predicts a
*uniform fractional* continued-erosion (cost≪A0 ⇒ `λ≈λ_max`) where the observed absolute drops are
similar; and **PCB roughly doubles 2.5→5 mo**, so reality cannot fully exclude a PCB (vs
time-integration) cause — only the *model* shows time-integration suffices. A **proof-of-concept**,
not a powered validation.

**Firm-up (Veldhuizen-Tsoerkan et al. 1991, ACET 20:259–265) — PCB confound removed.** The
companion study exposed *M. edulis* to Cd **or** PCB **separately**: **Cd alone erodes SoS
progressively** (lab LT50 10.7→9.5→7.6 over 0/2/4 wk; semi-field 9.3→8.6 over 3/6 mo), and
**modelled margin tracks it ρ=+0.90** (n=5); **PCB alone** also erodes SoS, with a **delayed
onset** (no effect at 3 mo, effect at 6 mo). So the transplant's continued erosion is **not** a
PCB artifact — a single toxicant suffices, time-/dose-dependently. It does **not** add a clean
*constant-burden* continued-erosion test (burden rises through every measured point). The dynamic
claim is thus **de-confounded and reinforced**, still short of a *powered* test. *Detail:*
`sos_validation_status.md` (dynamic).

---

## 8. The consistent null — the amplification scalar `g`/`F`

Across **every** anchor the one-dimensional amplification scalar predicts nothing: COMADRE
(−0.05…−0.13, the lone bootstrap CI spanning 0), GlobTherm (`g`/`E_m` null in the clean ectotherm
subset), and every SFG/SoS test (the *margin state* `A_t` is the predictor, never `g`/`F`). This
is the **margin-first prediction confirmed**: external support is for the margin/recovery layer,
not for a scalar amplification readout.

---

## 9. Cross-cutting honest assessment

- **Rank-robust, magnitude-modest, specification-sensitive.** Effects are ρ≈0.2–0.45 (COMADRE,
  SFG, SoS) — corroboration, not strong prediction — and several are robust in *ranks* but weak in
  log-linear form (COMADRE `k_M`). Report as monotone tendencies.
- **Scale bridge.** COMADRE needs an individual→population bridge (argued via DEB-structured
  models). SFG/SoS do **not** — they are at the margin's own level; that is their methodological
  strength.
- **The capacity weighting is still untested.** Every SFG/SoS study is single-species, so the AmP
  capacity (A0, κ-rule axis weights — the model's distinctive content) is held constant; the tests
  validate the *erosion mechanism + MoA aggregation*, not the weighting. Testing it needs
  across-species contaminant-gradient data, which is largely absent (the SFG corpus is
  mussel-dominated; ICES DOME no longer holds SFG; non-mussel SFG is temperature- not
  contaminant-driven).
- **Tissue burden ≠ exposure.** A threshold-free median-normalised relative burden is the pressure
  proxy; where burden tracks food/condition rather than exposure, it fails (Albentosa).
- **Metals as a positive confound** recur across SFG and SoS (As/Cd/Zn) — the MoA routing exists
  precisely to keep them off the toxic axis.
- **The dynamics are only a proof-of-concept** (§7); `B_t`/the erosion ODE have one small,
  qualitatively-positive, confound-limited test.

---

## 10. Validated vs open — one-paragraph summary

External evidence supports the **margin/recovery layer** at four levels: its **rate endpoints**
(COMADRE `k_M`, `R_i`), its **state** under pressure (SFG, where burden indexes exposure), its
**function** as acute-stress resilience (SoS, with the condition confound *strengthening* the
signal), and a **first positive sign of its dynamics** (the transplant continued-erosion result on
the model's own unfitted timescale). The **amplification scalar `g`/`F` is null throughout**,
exactly as the margin-first reframe predicts. **Open:** a real-phylogeny control (dated-tree PGLS,
one download away); the **capacity weighting** (needs across-species data that largely does not
exist); and a **powered dynamic test** (a denser single-contaminant time series, in progress).

---

## 11. Reproducibility

| anchor | data (provenance) | script / harness |
| --- | --- | --- |
| COMADRE | COMADRE `.RData` (gitignored); derived CSV committed | `scripts/extract_comadre_recovery.jl`, `examples/comadre_*`, `scripts/comadre_*` |
| COMADRE robustness | — | `examples/comadre_bootstrap_effectsizes.jl`, `scripts/comadre_filter_sensitivity.jl`, `comadre_ap_diagnostic.jl` |
| PGLS (pending tree) | dated Newick (manual) | `scripts/comadre_pgls_dated.jl` |
| GlobTherm | `GlobalTherm.csv` (gitignored); `globtherm_amp_matched.csv` committed | `scripts/extract_amp_for_globtherm.jl`, `examples/globtherm_validation.jl` |
| SFG (Widdows 1995) | `sfg_widdows1995_*.csv` | `examples/sfg_margin_validation.jl` |
| SFG (Widdows 2002) | `sfg_widdows2002_*.csv` (SFG figure-digitised) | `examples/sfg_margin_validation_widdows2002.jl` |
| SFG (Albentosa 2012) | `sfg_albentosa2012_*.csv` | `examples/sfg_margin_validation_albentosa2012.jl` |
| SoS static | `sos_dome_ukcemp.csv` (raw DOME gitignored) | `scripts/extract_dome_sos.jl`, `examples/sos_margin_validation_dome.jl` |
| SoS temporal | `sos_dome_ukcemp_yearly.csv` | `scripts/extract_dome_sos_yearly.jl`, `examples/sos_temporal_validation_dome.jl` |
| SoS controlled (Viarengo) | `sos_viarengo1995_doseresponse.csv` | `examples/sos_mixture_validation_viarengo.jl` |
| Dynamics (transplant) | `sos_veldhuizen1991_transplant.csv` (digitised) | `examples/sos_dynamic_validation_veldhuizen.jl` |
| Dynamics firm-up (single-contaminant) | `sos_veldhuizen1991_singlecontaminant.csv` | `examples/sos_dynamic_firmup_veldhuizen_singlecontaminant.jl` |

All rank statistics; pre-registered signs; raw downloads gitignored, derived CSVs committed.
Run via `julia +release --project=. <harness>` (Julia 1.12.6).

---

## 12. Key references

- COMADRE: Salguero-Gómez et al. 2016, *J. Anim. Ecol.* (database). AmP: Add-my-Pet / DEB.
- GlobTherm: Bennett et al. 2018, *Sci. Data* 5:180022.
- SFG: Widdows et al. 1995, *MEPS* 127:131; Widdows et al. 2002, *Mar. Environ. Res.* 53:327;
  Albentosa et al. 2012, *STOTEN* 435–436:430.
- SoS: ICES DOME 2024 OSPAR CEMP (figshare 27211422); Veldhuizen-Tsoerkan et al. 1991,
  *Arch. Environ. Contam. Toxicol.* 21:497–504 (transplant) and 20:259–265 (controlled, firm-up).
