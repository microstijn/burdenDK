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

— is **already implemented** as the `ec50_anchored_fractional_impairment` response mode in `src/reduced_deb_response.jl`:

- `compute_adaptive_margin_response(...; response_mode="ec50_anchored_fractional_impairment")` (`reduced_deb_response.jl:599`) computes `A_t = A0 * max(floor, 1 − Q_t)`.
- `Q_t = Σ wᵢ·Eᵢ`, with `Eᵢ = xᵢ/(1+xᵢ)` from `ec50_anchored_fractional_impairment` (`reduced_deb_response.jl:465`) — the exact bounded shape the review proposes.
- `wᵢ` = normalized α-axes via `axis_weights_for_species` (`reduced_deb_response.jl:482`).
- `compute_adaptive_margin_response_from_impairment` (`reduced_deb_response.jl:649`) is the precomputed-impairment entry to the same form.

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

1. **Reconcile the docs** (done alongside this audit): refresh the `2026-05-29` dates, add the 17 missing source files, and correct the `Z_t` status from `not_implemented` to `implemented_opt_in`.
2. **Pin current behavior with regression tests** (agenda #2, done alongside this audit): characterization tests for `λ_max/λ_min == 1/κ` (κ = `alpha_axes[3]`), `KA == 0.3·A0`, and `Fmax == 1 + (1/κ − 1)/1.3`, run across the shipped library. These will (intentionally) fail when the `0.3` knob is removed — that failure is the signal that the refactor changed as-built behavior.
3. **Canonicalize the response mode** before touching the margin equation (replaces "implement the nondimensional A_t"). The nondimensional form exists; decide its status relative to `raw_margin_subtraction`.
4. **Break the Fmax–κ collapse / remove the `0.3`** in `AmP_Translator.jl` + `restoring_force_from_margin`, then regenerate the JSON and re-run the response-capacity diagnostics. This is the one load-bearing math change the audit confirms is still outstanding.
5. Add the κ ≥ 1 guard at `AmP_Translator.jl`.

The single ordering question from the review is unchanged: *does the framework still produce meaningful, defensible amplification once the Fmax–κ collapse is broken?* The audit narrows the surface needed to answer it — the dimensional fix is largely already in the tree; the κ-collapse is not.
