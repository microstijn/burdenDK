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
| 1 | **COMADRE** (matrix demography) | population | recovery **rate endpoints** `k_M`, `R_i` | `k_M`↔recovery **+0.19–0.22\*** (rank; survives pace + **dated-tree** PGLS); `R_i`↔compensation **+0.77\*\*** | ✅ corroborated |
| 2 | **GlobTherm** (thermal tolerance) | individual physiology | is the capacity axis *general* resilience? | coherent (\|ρ\|≤0.45) but **recovery-specific** (general-resilience refuted) | ✅ bounding |
| 3 | **Scope for Growth** (3 studies) | individual energetics | the margin **state** under pressure | **+0.41** (estuary) → +0.12 (basin) → −0.11 (confounded) | ✅ / ◐ scale-dependent |
| 4 | **Stress-on-Stress** (ICES DOME) | individual energetics | the margin **function** (acute resilience), static | **+0.39 → +0.45** (confound-controlled), +0.62 QC-cleaned | ✅ static map |
| 4b | **Viarengo 1995** (controlled SoS) | individual, controlled dose | the **impairment curve + mixture model** | monotone dose-response (potency Cu>PAH>PCB); mixture additive, no antagonism (TU/IA bracket) | ✅ controlled |
| 5 | **Transplant + single-contaminant time-courses** (Veldhuizen 1991) | individual, *over time* | the margin **dynamics** (sustained-burden erosion) | dynamics reproduce continued erosion the static map can't (n=4); Cd-alone erodes SoS progressively, ρ(margin,LT50)=+0.90 (PCB confound removed) | ◑ proof-of-concept, de-confounded |
| 5-bis | **Phenanthrene time×dose** (Dellali 2023, §7-bis) | individual, *over time* | the margin **dynamics**, controlled | static flat vs dynamic erosion across 4 wk × 3 doses; ρ(erosion,LT50)=**−0.99** (n=12, non-trivial inversions); single PAH, flat control | ◑ proof-of-concept, cleaner/complete |
| 5b | **Single-trait `k_M`→toxicity** (ECOTOX LC50 n=310; Rubach `k_out` n=6) | individual, cross-species | does maintenance predict toxic response *beyond body size*? | raw maintenance↔sensitivity ρ≈−0.27 (all 4 chemicals) **nulls under a size control** (partial −0.03); rate axis weak/n.s. | ✅ bounding (size-confounded) |
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
- **Real-phylogeny control — DONE (2026-06-13); the rank signal survives it.** A dated **TimeTree**
  (182 spp, real branch lengths) over the COMADRE-matched species, with Pagel's λ estimated by ML.
  Two specifications, and the contrast *is* the result:
  - *Linear (log-linear) PGLS:* `k_M`→recovery is significant alone (β\*=0.30, p=0.014) but **nulls
    under a generation-time covariate** (β\*=0.009, p=0.96). The log-linear form is weak (as already
    known, +0.04) and Pagel's λ≈0 (the damping-ratio trait carries little phylogenetic signal).
  - *Rank PGLS (phylogenetic Spearman):* rank-transform `y`, predictor and generation time, same dated
    VCV. **`k_M`→recovery survives: β\*=0.221, p=0.011** under the generation-time control — barely
    below the non-phylogenetic Spearman partial (+0.264\*\*), i.e. the dated-tree correction only mildly
    attenuates it. `λ(A0)` drops to n.s. (β\*=0.10); `g` is null (β\*=−0.05).
  **Net: the `k_M`↔recovery signal is a genuine rank/monotone effect that survives *both* pace-of-life
  *and* a real dated-phylogeny control — only its log-linear form is weak.** This is the strongest
  control combination the programme has applied, and the headline anchor holds in ranks.
  *Output:* `comadre_pgls_dated_results.txt` (linear), `comadre_pgls_dated_rank_results.txt` (rank);
  tree `comadre_amp_dated_tree.nwk` (TimeTree, 184 tips; list `comadre_species_for_timetree.txt`).

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

### 7-bis. Phenanthrene firm-up — a cleaner, two-species controlled dynamic replicate (Dellali et al. 2023)
`examples/sos_dynamic_validation_dellali.jl`; data `data/external/sos_dellali2023_phenanthrene.csv`.
**Dellali et al. 2023** (*Animals* 13(1):151) exposed **both** *M. galloprovincialis* and *R. decussatus*
to **constant** waterborne phenanthrene (Table 1: ~10/45/89 µg/L, flat over time) and measured survival-
in-air LT50 at 7/15/21/28 d. A **cleaner** dynamic test than the transplant: **single PAH** (→ assimilation
axis, same routing as Widdows/DOME), a **flat control**, a **3-level dose gradient** — no metal/PCB confound,
exposure genuinely constant (so the static map's prediction is unambiguous).

- **The discriminator.** Constant exposure ⇒ a static burden→margin map predicts LT50 **flat** in time.
  Observed LT50 instead **declines progressively and dose-dependently** (mussel: control 8.6→8.15 d (5%, flat),
  WC1 7.42→4.31 d (42%), WC2 6.4→2.77 d (57%); clam: control 13→11.7 d (10%), WC1 9.49→4.73 d (50%),
  WC2 7.67→3.85 d (50%), over 7→28 d). The **dynamic** erosion keeps rising over the four weeks
  (1/λ ≫ 28 d) with **dose-ordered** magnitude; **static change = 0** for every treatment, both species.
- **Two independent 12-cell replicates:** ρ(dynamic erosion, LT50) = **−0.99** (mussel) and **−0.97**
  (clam). Unlike the transplant's small-n ρ, these are **not monotone by construction** — they carry genuine
  dose×time inversions (e.g. WC1@28 d < WC2@15 d), so the erosion state is a **non-trivial** predictor.
- **Cross-species (n=24):** raw pooled ρ = **−0.82**, but the clam is far more anoxia-tolerant at baseline
  (control ~13 vs ~8.6 d) — the paper ascribes this to **anoxia physiology** (tighter shell closure,
  anaerobic metabolism), **not** the contaminant margin. **Control-normalising** (the contaminant effect
  proper) restores ρ = **−0.91**. The absolute cross-species tolerance gap is the **untested capacity-
  weighting** question (§9), not a dynamics claim.
- **Honest scope.** 7 d & 28 d text/Table-3 reported (max→7 d, min→28 d); **15 d & 21 d figure-digitised**
  from the Figure 4 species panels (±~0.5 d → rank-reliable). Fractional erosion is dose-uniform
  (cost≪A0 ⇒ λ≈λ_max), so the model captures the **trajectory shape + dose-ordering**, not the dose×time
  magnitude interaction.
- **Net:** a cleaner, complete, **two-species** dynamic proof-of-concept corroborating the transplant — the
  static map cannot produce the temporal decline; the unfitted-`k_M` dynamics can, in both species.
  Powered validation still needs reported, multi-species series spanning a wider capacity range.

## 7b. Direct test of the maintenance-timescale claim — single-trait `k_M`→toxicity is body-size-confounded

§7 rests on the prediction that the margin's erosion/recovery timescale is set by the maintenance
rate (`1/λ ≈ 1/k_M`). That prediction was tested **directly and out-of-sample**: does AmP `k_M`
predict an independent, **non-DEB** toxicity rate or threshold, *controlling body size*? Two axes,
both reducing to a **body-size confound** — a bounding result in the spirit of §4 (GlobTherm).

- **State axis (`k_M` → acute sensitivity), n=310.** Replicating Baas & Kooijman (2015): species-level
  acute LC50/EC50 for their four AChE inhibitors (chlorpyrifos, malathion, carbaryl, carbofuran),
  pulled from **raw EPA ECOTOX** and matched to AmP. The maintenance↔sensitivity link **replicates raw**
  (higher `k_M`/`[p_M]` → more sensitive; ρ≈−0.27 pooled, **same sign in all four chemicals**) but
  **vanishes under a body-size control** (partial ρ ≈ −0.03; −0.01…+0.10 per chemical). At n=310 this
  is a **robust null, not underpower**: small species are both faster-`k_M` and more sensitive, so
  `k_M` carries no signal *beyond* size.
- **Rate axis (`k_M` → elimination rate), n=6.** Against measured chlorpyrifos elimination constants
  `k_out` (Rubach et al. 2010, 15 freshwater arthropods, radiotracer — non-DEB), AmP overlap is 6
  species and `k_M`↔`k_out` is weak/wrong-signed (ρ≈−0.5, n.s.). `k_out` is anyway the **toxicokinetic**
  (chemical-clearance) rate, not the thesis-relevant **toxicodynamic recovery** rate `k_r`; the clean
  `k_r` (GUTS-proper; e.g. Nyman et al. 2012, *G. pulex*/propiconazole `k_r`=1.0–2.3 d⁻¹) is **scarce,
  *G. pulex*-centric and chemical-specific**, leaving a *powered* dynamic rate test **data-starved**.

**Reading.** Like GlobTherm (§4), this **bounds rather than refutes**: the model's distinctive leverage
is **not** "`k_M` predicts toxic response" — that is a body-size story — but the **across-axis capacity
weighting** (§9), which every single-trait *and* single-species test here holds constant. Both axes
stand as honest negative controls. *Detail:* `guts_kd_dynamic_test_scoping.md`.

---

## 7c. Cross-species capacity pilot — 5 fish, benzovindiflupyr (single chemical)

`examples/benzovindiflupyr_capacity_probe.jl`; data `data/external/benzovindiflupyr_fish_survival.csv`
(extracted by `scripts/extract_benzovindiflupyr_fish.jl` from the **Nickisch Born Gericke et al. 2022**
*ETC* 41(7):1732–1741, doi:10.1002/etc.5348 SI; underlying acute survival = Ashauer et al. 2013). The
first **reported, multi-species** single-chemical cross-species set wired in. Benzovindiflupyr is an SDHI
(complex-II inhibitor) → energetic cost → the **maintenance** axis. Per-species acute (96 h) LC50: carp
3.5 < fathead 4.8 < trout 8.9 < sheepshead 26 < bluegill 28.5 µg/L (~8×).

- **n=5 → pilot, direction only** (only |ρ|>~0.9 reaches p<0.05). A single chemical exercises **one** axis,
  so this is a **§7b-style** single-trait cross-species probe (new chemical class/taxon), **not** the
  across-axis weighting.
- `k_M`→sensitivity = **−0.90** raw (higher `k_M` → more sensitive), **−0.92** partialling body size — here
  `k_M` *survives* the size control, the **opposite** of the powered §7b (n=310). **But**: n=5, and the only
  size proxy AmP exposes is **structural** length `L_m`, which does **not** rank physical body size cleanly
  (δ_M varies; ρ(`L_m`,LC50)=−0.10) — so this **cannot** overturn the n=310 result; a noisy pilot with a
  weak size control.
- The framework's **distinctive axis weighting** (`alpha_maint`) does **not** predict sensitivity
  (ρ=+0.10, partial −0.10) — as expected for a single-MoA set. `A0` is a size proxy (ρ(`A0`,`L_m`)=+0.90).
- **Net:** the cross-species GUTS pipeline now runs reproducibly, but this single-chemical pilot **confirms
  the open question stands** (§9): single-trait `k_M` is suggestive-but-confounded/underpowered, and the
  real across-axis **weighting** test needs **multiple MoA × species** data — which this is not.

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
  log-linear form. The cleanest case is COMADRE `k_M`: the rank effect **survives pace + a dated-tree
  PGLS** (phylo-Spearman β\*=0.22, p=0.011) while its log-linear form nulls under the same controls.
  Report as monotone tendencies, not linear effects.
- **Scale bridge.** COMADRE needs an individual→population bridge (argued via DEB-structured
  models). SFG/SoS do **not** — they are at the margin's own level; that is their methodological
  strength.
- **The capacity weighting is still untested — and is carried as an assumption.** Every SFG/SoS study
  is single-species, so the AmP capacity (A0, κ-rule axis weights — the model's distinctive content)
  is held constant; the tests validate the *erosion mechanism + MoA aggregation*, not the weighting.
  Testing it needs across-species contaminant-gradient data, which is largely absent (the SFG corpus
  is mussel-dominated; ICES DOME no longer holds SFG; non-mussel SFG is temperature- not
  contaminant-driven). Pending such data it is stated as a **model assumption** — like the mixture
  rules — not a validated result.
- **Single-trait maintenance is a body-size story (§7b).** A well-powered cross-species test (n=310)
  shows the raw `k_M`→sensitivity correlation is *fully accounted for by body size* (partial ≈ 0). The
  model's leverage is therefore the across-axis *weighting*, not `k_M` as a scalar predictor — the
  strongest form yet of the "specification-sensitive" signature (a proper size control nulls it).
- **Tissue burden ≠ exposure.** A threshold-free median-normalised relative burden is the pressure
  proxy; where burden tracks food/condition rather than exposure, it fails (Albentosa).
- **Metals as a positive confound** recur across SFG and SoS (As/Cd/Zn) — the MoA routing exists
  precisely to keep them off the toxic axis.
- **The dynamics are only a proof-of-concept** (§7); `B_t`/the erosion ODE have one small,
  qualitatively-positive, confound-limited test.

---

## 9b. The structure earns its keep — within-anchor ablation (routed vs naive)

Every field anchor admits a direct ablation: replace the MoA-routed, bounded, AmP-weighted margin
with a **naive equal-weight load index** (mean median-normalised burden over all contaminants, no
routing, no capacity structure) and ask which better tracks the outcome. Consolidated across the
three field anchors with a contaminant gradient:

| anchor (outcome) | naive equal-weight load | routed margin |
| --- | --- | --- |
| Widdows 1995, estuary (SFG) | +0.22 | **+0.41\*** |
| Widdows 2002, basin (SFG) | +0.005 | **+0.12** |
| Stress-on-Stress, DOME (survival-in-air) | +0.32 | **+0.39** (+0.45 \| size+cond.) |

**The routed margin beats the naive index at every anchor** — the response *structure* is not
decorative. **Precise scope:** this supports the **operative structure** (MoA routing, the saturating
`E=x/(1+x)`, axis aggregation) over an unstructured index; it does **not** test the across-**species**
capacity weighting `w_i`, which is held constant within each single-species anchor and remains the
open question (§9, §10). A single hand-picked axis (e.g. PAH/assimilation alone, 0.43–0.47) can edge
the routed margin, but that axis is not known a priori; the routed margin is the **choice-free**
aggregate that consistently improves on naive. *This is the strongest internal support for the
structure available without across-species gradient data, and it is now consolidated as
Table (`tab:ablation`) in the manuscript (`\subsection{Does the structured margin earn its keep?}`).*

### Response-curve-form robustness — the ranking does not depend on `E=x/(1+x)`
`examples/response_curve_sensitivity.jl`. The licensed use is a *relative ranking*, so its
defensibility rests on ranking stability — tested here against the **impairment curve form** (in
addition to the MoA routing already tested in the Discussion). The engine applies `E(x)=x/(1+x)`;
feeding a transformed pressure `x' = E_alt(x)/(1−E_alt(x))` reproduces any alternative curve *through
the real engine* (Hill-`h`: `x'=x^h`; saturating-exp: `x'=exp(x)−1`). All alternatives are
threshold-free, monotone, bounded, and **half-saturate at the reference `x=1`** — shape probes, **no
half-saturation knob reintroduced** (cf. the removed `K_A`).

| anchor | margin↔outcome ρ (baseline → Hill h=2 / h=0.5 / 1−e^−x) | ranking stability (rank ρ vs baseline) |
| --- | --- | --- |
| Widdows 1995 SFG (n=36) | +0.405 → +0.402 / +0.407 / +0.406 | **0.988–0.997** |
| DOME stress-on-stress (n=17) | +0.392 → +0.426 / +0.392 / +0.414 | **0.990–1.000** |

**Neither the corroboration nor the licensed ranking depends on the specific impairment form** (ρ moves
within ±0.03; rank ρ ≥ 0.99). Manuscript: Discussion `\paragraph{The relative ranking is robust to the
response-curve form.}`.

### Bridge figure — the licensed relative use, on validated ground
`examples/dome_margin_ranking_figure.jl` → `docs/wiki/figures/dome_margin_ranking.png` (+
`docs/tex/dome_margin_ranking.pdf`). The 17-station ICES DOME network **ranked by modelled adaptive
margin** `A_t/A0`, coloured by the independent survival-in-air outcome (warmer = more resilient
clusters toward larger retained margin; in-script `ρ(margin, survival)=+0.392` raw, reproducing the
SoS anchor). Built **entirely from validated SoS machinery**; **no modelled water concentration enters**
— it demonstrates the relative, mechanistically-structured ranking the framework licenses (Discussion
"licensed use"), not the spatial coupling. This is the figure that speaks to the water-quality reader
on already-defended ground (manuscript `fig:ranking`).

## 10. Validated vs open — one-paragraph summary

External evidence supports the **margin/recovery layer** at four levels: its **rate endpoints**
(COMADRE `k_M`, `R_i`), its **state** under pressure (SFG, where burden indexes exposure), its
**function** as acute-stress resilience (SoS, with the condition confound *strengthening* the
signal), and a **first positive sign of its dynamics** (the transplant continued-erosion result on
the model's own unfitted timescale). The **amplification scalar `g`/`F` is null throughout**,
exactly as the margin-first reframe predicts. A **direct cross-species test** of the maintenance
claim (§7b, n=310) finds the single-trait `k_M`→toxicity signal is **body-size-confounded** —
bounding, like GlobTherm. The **real-phylogeny control is now done** (dated TimeTree PGLS, §3): the
`k_M`↔recovery signal **survives in rank form** under pace + a real dated-phylogeny correction
(phylogenetic Spearman β\*=0.221, p=0.011); only its log-linear form is weak. **Open:** the **capacity weighting**, which the single-species corpus cannot test and is therefore carried as
a **model assumption** (like the mixture rules), pending across-species gradient data that largely
does not exist; and a **powered dynamic test**, which this session found is **data-starved** for the
*toxicodynamic recovery rate* `k_r` (scarce, chemical-specific; the abundant rate/threshold endpoints
reduce to body size, §7b) — **but a well-matched *erosion-dynamics* firm-up dataset exists and is now
partially run (§7-bis):** Dellali et al. 2023 (*Animals* 13(1):151, doi:10.3390/ani13010151) report
**weekly survival-in-air LT50 of *M. galloprovincialis* under graded phenanthrene (~10/45/89 µg/L) at 7,
15, 21, 28 d, with a flat control** — a controlled, *time- and concentration-resolved* single-PAH
(→ assimilation axis) survival-in-air series. The **full 4-week × 3-dose grid now replicates the
transplant** (static flat; dynamic erosion rises, dose-ordered; ρ=−0.99, n=12; §7-bis), with the 15/21 d
points figure-digitised; powered validation still needs reported (not figure-read), multi-species series. (Companion: the openGUTS / Nyman *G. pulex*–propiconazole `k_r`
data remain the toxicodynamic-rate route, but are *G. pulex*-centric.)

---

## 11. Reproducibility

| anchor | data (provenance) | script / harness |
| --- | --- | --- |
| COMADRE | COMADRE `.RData` (gitignored); derived CSV committed | `scripts/extract_comadre_recovery.jl`, `examples/comadre_*`, `scripts/comadre_*` |
| COMADRE robustness | — | `examples/comadre_bootstrap_effectsizes.jl`, `scripts/comadre_filter_sensitivity.jl`, `comadre_ap_diagnostic.jl` |
| PGLS (DONE) | `comadre_amp_dated_tree.nwk` (TimeTree, 184 tips; list `comadre_species_for_timetree.txt`); results `comadre_pgls_dated_results.txt` (linear), `comadre_pgls_dated_rank_results.txt` (rank) | `scripts/comadre_pgls_dated.jl` (linear), `scripts/comadre_pgls_dated_rank.jl` (rank/phylo-Spearman) — both need `Distributions` env |
| GlobTherm | `GlobalTherm.csv` (gitignored); `globtherm_amp_matched.csv` committed | `scripts/extract_amp_for_globtherm.jl`, `examples/globtherm_validation.jl` |
| SFG (Widdows 1995) | `sfg_widdows1995_*.csv` | `examples/sfg_margin_validation.jl` |
| SFG (Widdows 2002) | `sfg_widdows2002_*.csv` (SFG figure-digitised) | `examples/sfg_margin_validation_widdows2002.jl` |
| SFG (Albentosa 2012) | `sfg_albentosa2012_*.csv` | `examples/sfg_margin_validation_albentosa2012.jl` |
| SoS static | `sos_dome_ukcemp.csv` (raw DOME gitignored) | `scripts/extract_dome_sos.jl`, `examples/sos_margin_validation_dome.jl` |
| SoS temporal | `sos_dome_ukcemp_yearly.csv` | `scripts/extract_dome_sos_yearly.jl`, `examples/sos_temporal_validation_dome.jl` |
| SoS controlled (Viarengo) | `sos_viarengo1995_doseresponse.csv` | `examples/sos_mixture_validation_viarengo.jl` |
| Dynamics (transplant) | `sos_veldhuizen1991_transplant.csv` (digitised) | `examples/sos_dynamic_validation_veldhuizen.jl` |
| Dynamics firm-up (single-contaminant) | `sos_veldhuizen1991_singlecontaminant.csv` | `examples/sos_dynamic_firmup_veldhuizen_singlecontaminant.jl` |
| Maintenance→toxicity, state axis (§7b) | EPA ECOTOX ASCII (gitignored); `ecotox_acute_4chem.csv`, `state_axis_ecotox_amp_paired.csv` committed | `scripts/extract_ecotox_acute.awk` (extract), `scripts/state_axis_ecotox_amp.jl` (analysis) |
| Maintenance→toxicity, rate axis (§7b) | Rubach 2010 Table 2 (transcribed); `rubach2010_kM_kout_paired.csv` committed | `scripts/rubach2010_rate_axis.jl` |

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
- Maintenance→toxicity (§7b): Baas & Kooijman 2015, *Ecotoxicology* 24:657 (metabolic rate↔sensitivity);
  Rubach et al. 2010, *Environ. Toxicol. Chem.* 29:2225 (15-species chlorpyrifos toxicokinetics);
  Nyman et al. 2012, *Ecotoxicology* 21:1828 (GUTS-proper `k_r`); US-EPA ECOTOX (ASCII release 03/2026).
