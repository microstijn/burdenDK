# TwoTimescaleResilience — Source Audit (Addendum)

*Working note, 2026-06-11. Companion to `TwoTimescaleResilience_review_and_agenda_2026-06-11.md`. This is the result of agenda task #1 ("Tidy and audit"): verify that the documented self-description matches the actual `src/`, and check the review's specific technical claims against the real files. Append-only: this addendum refines the review where the code disagrees with it; it does not overwrite the review.*

---

## How to read this

The review note is conceptually sound and its numbers are correct, but two of its claims are mis-located against the current tree, and one of its central recommendations turns out to be **already implemented**. The sharpest scientific finding (the Fmax–κ collapse and the hidden `0.3` gain) survives the audit fully intact. Findings are ordered by how much they change the agenda.

---

## Finding 1 — The mapping math is in `AmP_Translator.jl`, not `amp_library.jl`

The review analyses "`amp_library.jl`" for the DEB compound-parameter derivations. That attribution is wrong for the current tree:

- **`src/amp_library.jl`** is *only* a JSON adapter: `load_amp_species_library`, `validate_amp_record`, `amp_record_to_deb_params`, etc. It reads precomputed `A0`, `alpha_axes`, and `lambda_bounds` straight out of `data/AmP_Species_Library.json`. It contains none of the `{p_Am, p_M, κ, v}` math.
- **`src/AmP_Translator.jl:45–65`** is where the derivation actually lives. Every numeric claim in the review checks out against it verbatim:
  - `E_m = p_Am / v`; `A_0 = E_m`
  - `alpha_A = 1/E_m`; `alpha_M = p_M/(κ·p_Am)` (= 1/L_m); `alpha_G = κ`; `alpha_R = 1−κ`
  - `L_m = κ·p_Am/p_M`; `lambda_max = v/L_m`; `lambda_min = p_M/A_0`
  - `K_A = A_0 * 0.3` — the hidden gain, confirmed at `AmP_Translator.jl:65`.

**Two facts the review missed about this file:**
1. It is **not `include`d** in the module (`TwoTimescaleResilience.jl` never includes it). It is a standalone offline data-generation script.
2. It depends on **`MAT`** (reads `data/allStat.mat`).

**Consequence for the agenda:** "fix the mapping" means *regenerate the JSON via `AmP_Translator.jl`*, not edit live module code. The α-axes and `KA` reaching the model are frozen artifacts in the JSON.

## Finding 2 — The review's proposed fix already exists as a second response mode ⚠️

This is the finding that most changes the work program. The review's "fix direction" —

> `A_t = A0·(1 − Σ_i w_i · g_i(s_i))`, each `g_i` a bounded dimensionless impairment (the `x/(1+x)` shape), weights `w_i` tied to the α split

— is **already implemented** as the `ec50_anchored_fractional_impairment` response mode in `src/deb_axes.jl` (note: the orientation docs misattribute these functions to `reduced_deb_response.jl`, which is actually a small file holding only `event_cost_from_axes` / `simulate_reduced_deb_response` — another stale-docs artifact):

- `compute_adaptive_margin_response(...; response_mode="ec50_anchored_fractional_impairment")` (`src/deb_axes.jl`) computes `A_t = A0 * max(floor, 1 − Q_t)`.
- `Q_t = Σ wᵢ·Eᵢ`, with `Eᵢ = xᵢ/(1+xᵢ)` from `ec50_anchored_fractional_impairment` (`src/deb_axes.jl`) — the exact bounded shape the review proposes.
- `wᵢ` = (as-audited) normalized α-axes via `axis_weights_for_species` (`src/deb_axes.jl`).
- `compute_adaptive_margin_response_from_impairment` (`src/deb_axes.jl`) is the precomputed-impairment entry to the same form.

The **dimensional-heterogeneity and numerical-inertness critique therefore applies only to the *other* mode**, `raw_margin_subtraction` (`deb_adaptive_margin`, `deb_axes.jl:91`: `A_t = A0 − Σαᵢ·sᵢ`). Both modes coexist; the calibrated demo runs and compares both side by side (`examples/ecotox_amp_multiaxis_response_calibrated_demo.jl:219`).

**Consequence for the agenda:** task #3 ("implement the nondimensional A_t") should be **rewritten** as:
- (a) **Canonicalize.** Decide whether `raw_margin_subtraction` is retired, kept as a diagnostic, or kept as the default. Right now both are first-class, which is itself a source of the "documentation ahead of code" ambiguity.
- (b) Audit the weights: `wᵢ` = *normalized* α is α-derived but not explicitly the "κ-rule split" the review imagines; because `alpha_A = 1/E_m` is tiny and `alpha_G = κ`, `alpha_R = 1−κ` are O(1), the normalized weights are in practice dominated by the κ / (1−κ) / (1/L_m) terms. Worth a written justification either way.

## Finding 3 — The Fmax–κ collapse and the `0.3` knob are real and mode-independent ✅ (survives)

`Fmax` is a property of `restoring_force_from_margin` (`deb_axes.jl:97`) and the JSON `KA`/λ-bounds, **not** of the A_t mapping — so it is identical in *both* response modes. With `KA = 0.3·A0`:

```
λ(A0) = λ_min + (λ_max − λ_min) · A0/(KA+A0) = λ_min + (λ_max − λ_min)/1.3
Fmax  = λ(A0)/λ_min = 1 + (λ_max/λ_min − 1)/1.3 = 1 + (1/κ − 1)/1.3
```

using the exact identity `λ_max/λ_min = 1/κ` (`λ_max = v/L_m`, `λ_min = p_M/A_0` ⇒ ratio `= 1/κ`). The `0.3` is an unjustified constant that only sets the `1/1.3` factor and violates the project's own no-knob invariant. **This is the finding that genuinely still needs a fix**, and it is untouched by Finding 2.

**Verified against real data:** the existing `test/test_amp_library.jl` pins `Abatus_cordatus` at `alpha_axes[3]` (κ) `= 0.77712`, `KA = 461.996 = 0.3 × A0 (1539.987)`, and `λ_max/λ_min = 0.011569/0.008990 = 1.2868 = 1/0.77712`. The identities hold exactly on shipped data.

## Finding 3b — Empirical confirmation: amplification is *functionally* κ-only

Run `examples/amp_kappa_collapse_diagnostic.jl` (read-only, whole library, 7,335 usable species). It does not just confirm the collapse — it shows it is stronger than "rank-locked":

| quantity | Spearman(·, κ) | max\|actual − predict-from-κ-alone\| |
| --- | --- | --- |
| Fmax (A_t→0) | −1.000 | 1.8×10⁻¹⁵ |
| F at realized erosion Q=0.2 | −1.000 | 4.4×10⁻¹⁶ |
| F at realized erosion Q=0.5 | −1.000 | 4.4×10⁻¹⁶ |
| F at realized erosion Q=0.8 | −1.000 | 4.4×10⁻¹⁶ |

`F(κ, Q) = [1 + (1/κ−1)/(1+c)] / [1 + (1/κ−1)(1−Q)/(c+1−Q)]` with `c = 0.3` reproduces the realized F to **machine epsilon** at every erosion level. So:

- **A0 contributes nothing.** It spans **2.7 → 2.86×10⁶ J/cm³ (six orders of magnitude)** across the library and cancels out of F entirely (because `KA ∝ A0` *and* `A_t ∝ A0`). The "capacity" signal the framework advertises is numerically inert in the output.
- **Three of the four α-axes contribute nothing** to F beyond their κ content.
- **The mechanism is near-off for the median species.** κ has median **0.928** (most species cluster near κ=1), giving median Fmax **1.059** — half the library can amplify at most ~6%. F only fires for the low-κ tail (min κ=0.069 → Fmax up to 11.4). Note `max κ = 1.000`, exactly the boundary the missing κ≥1 guard (Finding 7) protects.

**Also confirms Finding 2's weight caveat (D2):** the normalized-α assimilation weight has median `w_A = 0.00004`, and **100% of species have `w_A < 0.01`**. The assimilation axis is dead library-wide — an assimilation-targeting toxicant is ignored for every species.

### What this proves about the fix path

The diagnostic settles which proposed fixes touch the headline collapse and which do not:

- **Nondimensionalizing A_t does *not* break it.** The table above already uses the nondimensional operating point (`A_t = (1−Q)·A0`); F is still pure-κ. The collapse is entirely in the λ-curve, not the margin mapping.
- **Reweighting the axes does *not* break it.** It changes *which* `Q` a species reaches under a given stress (and fixes the dead assimilation axis — a real defect worth fixing), but F-as-a-function-of-Q stays pure-κ.
- **Reparametrizing `KA` onto normalized margin does *not* break it.** Working in `a = A_t/A0` with a fixed half-saturation is algebraically identical to `KA ∝ A0` — that proportionality *is* the cancellation. It only makes the knob honest.

**The collapse breaks only by severing `KA ∝ A0` or the `λ_max/λ_min = 1/κ` identity.** That is a scientific fork, now concrete:

1. **Keep `KA ∝ A0`** (recovery set by *relative* margin depletion) ⇒ **vulnerability ≡ κ**. Defensible if relative depletion is the right physics; then own it as a *result*, state it in the manuscript, and run the ecology backwards-check (is "low κ = more vulnerable" correct or inverted?).
2. **Make `KA` absolute / differently-scaled** (recovery set by *absolute* margin) ⇒ A0's six-order range re-enters and **vulnerability ≡ reserve capacity A0** (bigger buffer = more resilient). A different, possibly more intuitive, model.
3. **Re-anchor the λ-bounds** to two independently-varying DEB rates so their ratio is no longer `1/κ` but reflects the genuine slow/fast timescale separation.

This is exactly the first Kooijman question in agenda task #6, and the diagnostic shows it cannot be dodged: as built, the model asserts "chronic-then-acute vulnerability = allocation fraction κ, full stop." Whether that is a finding or a bug is a modelling decision, not a refactor — it must be made before the headline result is trusted.

## Finding 3c — Structural comparison: the κ-lock is in the λ-bounds, not `KA`

`examples/amp_lambda_structure_comparison.jl` (read-only, 7,335 species) re-runs the library under three single-lever variants of the restoring-force structure, to make the Finding 3b fork concrete. Each variant changes exactly one lever from the current baseline.

**F correlations at Q = 0.5 erosion:**

| structure | ρ(F, κ) | ρ(F, A0) | ρ(F, L_m) | reshuffle vs S0 |
| --- | --- | --- | --- | --- |
| S0 baseline (`KA = 0.3·A0`) | **−1.000** | −0.15 | 0.21 | 1.000 |
| S1 `KA` absolute (`0.3·median A0`) | −0.901 | −0.40 | 0.27 | 0.901 |
| S2 `λ_min` absolute floor (`median λ_min`) | **+0.176** | −0.06 | **−0.811** | **−0.176** |

(The full-erosion ceiling `Fmax` tells the same story; under S2 the max `Fmax` explodes to ~574 and 44% of species fall below the floor and cannot amplify.)

**The non-obvious result: the κ-lock is anchored in the λ-bounds, not in `KA`.**

1. **Making `KA` absolute (S1) is a *weak* lever.** F stays ~90% κ-ranked, A0 enters only weakly, the ranking barely moves. Structural reason: at the ceiling, `Fmax = λ(A0)/λ_min ≤ λ_max/λ_min = 1/κ` for *any* `KA` — `KA` only slides F *within* `[1, 1/κ]`; it cannot escape the κ-bounded ceiling. So the Finding 3b option (2) framing ("`KA` absolute → vulnerability ≡ A0") **overstated `KA`'s power**; it does not meaningfully decouple F from κ.
2. **Re-anchoring `λ_min` (S2) is the *real* lever.** F decouples from κ entirely (ρ → +0.18) and becomes size/rate-driven (ρ(F, L_m) = −0.81), full ranking reshuffle — but a crude flat floor is violent (44% non-amplifying, `Fmax` to ~574), so option (3) needs a *principled* re-anchoring, not a constant.

**Consequence — the fork is really a decision about the λ-bounds.** Specifically: should the slow-recovery floor `λ_min` stay reserve-normalized (`p_M/A0`, which is exactly what forces `λ_max/λ_min ≡ 1/κ`) or be anchored to a size/rate quantity? `KA` is a secondary shape knob. This sharpens the Kooijman question from "relative-vs-absolute margin" to **"what physiological quantity sets the slow recovery rate `λ_min`?"** Caveat: S1's specific magnitude depends on the chosen absolute `KA` (here `0.3·median A0`); the κ-bounded-ceiling argument, however, holds for any `KA`.

## Finding 4 — Physiological condition memory `Z_t` IS implemented (contradicts the docs)

`docs/PACKAGE_CAPABILITIES.md` and `docs/AI_CONTEXT_README.md` state `Z_t` is `not_implemented` / "no math." That is now **false**. `src/condition_buffer.jl` is a full `Z_t` layer:

- `ConditionBufferParams` (ρ_A, ρ_M, ρ_E, δ_Z, ω_Z, Z bounds)
- `condition_buffer_derivative`, `update_condition_buffer`, `simulate_condition_buffer`
- `adaptive_margin_with_buffer`: `A_t = A0 + ω_Z·Z − Σαᵢ·sᵢ`

It is wired into `deb_amplification_pipeline` (`deb_axes.jl:307`, via `buffer_grid`/`buffer_params`) and into `restoring_force_from_margin_and_axes` (`use_buffer_recovery_factor`, `beta_Z`). It is **off by default** (`beta_Z=0`, `use_buffer_recovery_factor=false`), so "not the default path" is fair — but "not implemented" must be corrected, and the `B_t`/`Z_t`/`D_t` separation invariant should now read "`Z_t` exists as an opt-in layer; keep it distinct from `B_t`; `D_t` remains unimplemented."

## Finding 5 — The documented inventory is stale

All three reference docs are dated **2026-05-29**. The actual `src/` contains **17 files absent from the inventory**:

`AmP_Translator.jl`, `ascii.jl`, `background.jl`, `condition_buffer.jl`, `default_isimip_moa.jl`, `exposure_filters.jl`, `grids.jl`, `isimip_deb_pipeline.jl`, `isimip_event_response.jl`, `metrics.jl`, `moa_deb_mapping.jl`, `mode_of_action.jl`, `plotting.jl`, `pulses.jl`, `simulation.jl`, `species_defaults.jl`, `vulnerability_tranche_comparison.jl`.

The module export list likewise carries names not in the documented Public API (tranche-comparison utilities, `simulate_two_timescale`/`PulseStressor`/pulse helpers, `compute_metrics`/`trapezoid_auc`, `background_index`/`BackgroundParams`, ASCII grid IO, etc.). The capability inventory should be regenerated from the actual module rather than hand-patched.

## Finding 6 — Archetype DB: present, not missing ✅

The review flagged "confirm before relying on files." Confirmed present: `data/AmP_Species_Archetypes.csv`, `data/AmP_Species_Archetypes.json`, and `test/test_amp_species_archetypes.jl`. This item needs no action beyond keeping its status truthful.

## Finding 7 — Minor: no κ ≥ 1 guard at the source

`AmP_Translator.jl:38–40` rejects `kap <= 0` but not `kap >= 1`. Since `λ_max/λ_min = 1/κ` assumes κ ∈ (0,1), an out-of-range κ would silently violate `λ_max > λ_min`. The guard belongs *here*, at generation, not in `validate_amp_record` (which never sees κ — it only sees the derived `A0`/α/λ). Risk is low (AmP κ are physical), but it is cheap insurance and was explicitly called for in the review.

---

## Revised next actions (supersedes agenda task #3 framing)

1. **Reconcile the docs** — *done.* Refreshed the `2026-05-29` dates, listed the 17 missing source files, corrected the `Z_t` status to `implemented_opt_in`.
2. **Pin current behavior with regression tests** — *done.* `test/test_amp_translator_identities.jl` pins `KA == 0.3·A0`, `λ_max/λ_min == 1/κ` (κ = `alpha_axes[3]`), and `Fmax == 1 + (1/κ − 1)/1.3` across the shipped library (29,347 assertions pass). These will (intentionally) fail when the `0.3`/λ-bound structure changes — that failure is the signal that the refactor changed as-built behavior.
3. **Quantify the collapse** — *done.* `examples/amp_kappa_collapse_diagnostic.jl` (Finding 3b) is the empirical baseline: F is functionally κ-only to machine epsilon; A0 and three α-axes are inert; the assimilation weight is dead library-wide.
4. **Two unambiguous fixes that do NOT depend on the scientific fork** — *implemented on branch `fix/margin-response-d1-d2` (point-level API only; grid/ECOTOX/ISIMIP pipelines still on the raw margin, deferred to a follow-up). They do **not** change the κ-collapse — that is Finding 3b's structural fork.*
   - **D1** — `compute_adaptive_margin_response` now defaults to the nondimensional `ec50_anchored_fractional_impairment` mode; `raw_margin_subtraction` retained as a documented diagnostic. Fixes margin inertness in the point API.
   - **D2** — `axis_weights_for_species` default is now `kappa_rule_assimilation_led`: dimensionless weights `w = [1/2, κ/4, κ/4, (1−κ)/2]` with `κ = α_G/(α_G+α_R)`. Assimilation now carries 0.5 instead of ~0.00004; the legacy `normalized_alpha_axes` remains as a diagnostic method.
5. **Resolve the structural fork (Finding 3b) — this gates the headline result.** Decide among: (1) keep `KA ∝ A0` and own "vulnerability ≡ κ" as a result + ecology backwards-check; (2) make `KA` absolute so "vulnerability ≡ reserve capacity A0"; (3) re-anchor the λ-bounds. This is review task #6's first Kooijman question and **cannot be decided from the code**. Do not implement the λ/`KA` change until it is answered. The `0.3` removal and κ≥1 guard ride along with whichever resolution is chosen.

The single ordering question from the review is now answerable in two halves. The *dimensional/inertness* half is largely already in the tree (Finding 2) and finished off by D1/D2. The *κ-collapse* half (Finding 3 + 3b) is real, structural, and reduces to one modelling decision — what physiological quantity sets recovery capacity. Everything downstream of that decision (manuscript framing, "mechanistic" vs "physiologically structured," whether low-κ-as-vulnerable is defensible) waits on it.
