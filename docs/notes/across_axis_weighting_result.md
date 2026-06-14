# Across-axis capacity-weighting test — result (2026-06-14)

**The central open validation question** (session-5 handover §3.1): does the model's per-species
across-axis weighting predict species × mode-of-action sensitivity *beyond body size*? Answer from a
powered multi-MoA × shared-species test built from the local EPA ECOTOX dump:
**not corroborated — and the conclusion is robust to endpoint choice on the powered (acute) side.**
A clean axis-matched *sublethal* test is what the framework really needs, but ECOTOX cannot supply it;
that requires a **curated apical-EC50 database** (scoped below).

## What the model predicts

Operative weights `w = (1/2, κ/4, κ/4, (1−κ)/2)` over (assim, maint, growth, repro), `κ = alpha_axes[3]`
(DEB allocation fraction). Assimilation is fixed → no cross-species signal. Discriminating contrast =
**maintenance** (`w_M = κ/4`) vs **reproduction** (`w_R = (1−κ)/2`); `w_M/w_R = κ/(2(1−κ))` rises in κ.
So **higher κ → relatively more maintenance-sensitive**. Test statistic: within-species, chemical-centered
`Δ = resid_R − resid_M` (`Δ>0` ⇒ relatively maintenance-sensitive); differencing removes the species-level
overall sensitivity **and body size**. Prediction: **ρ(κ, Δ) > 0**. Harness
`examples/across_axis_weighting_capacity_test.jl`; panel/tiers `data/ecotox_multimoa_panel.csv`
(core = defensible pMoA; stratum = AChE/DDT/PAH, contested).

## Result — PRIMARY: acute apical mortality (LC50 | MOR)

`examples/across_axis_weighting_endpoint_sensitivity.jl` (the clean, endpoint-homogeneous endpoint):

| panel | n | ρ(κ, Δ) | perm p |
| --- | --- | --- | --- |
| CORE | 19 | **−0.22** | 0.36 (null) |
| CORE+STRATUM | 88 | **−0.22** | 0.043 (significant, **wrong-signed**) |

The predicted sign (ρ>0) does not appear. **Endpoint-robust:** the broader pooled set
(LC50/EC50/IC50, any effect) gives the same answer — CORE ρ=−0.28 (p=0.15, n=27); CORE+STRATUM ρ=−0.20
(p=0.044, n=101); partial ρ(κ,Δ | log Ww) ≈ −0.29 / −0.17 (κ–size collinearity ≈ +0.03, so genuinely
size-independent). So the negative result is not an artifact of endpoint pooling.

**Reading.** The negative pattern is a **taxonomic target-site artifact**: low-κ arthropods (cladocerans,
*Chironomus*, mysids) are relatively maintenance-sensitive because of insecticide/uncoupler target-site
physiology — which the DEB *energetic* weighting does not (and should not) capture. The CORE (no AChE) is
a non-significant null; significance appears only when the contested AChE/DDT chemicals enter, and is also
wrong-signed. **Conclusion: the across-axis weighting is not corroborated by acute cross-species data.**

## Axis-matched SUBLETHAL test — attempted, ECOTOX too thin

`examples/across_axis_weighting_axismatched.jl`. To make each chemical engage its pMoA, match the ECOTOX
effect group to the axis (reproduction→REP; maintenance→GRO/PHY/DVP/MPH). Outcome:
**underpowered and endpoint-incoherent.** CORE n=**8** (ρ=−0.29, p=0.50); CORE+STRATUM n=24
(ρ=**+0.14**, p=0.51 — sign flips toward the prediction but far from significant). The apical sublethal
**EC50s** that would properly engage each pMoA are nearly absent — `GRO` EC50 = 6 sp, `REP` EC50 = 8 sp,
`PHY` EC50 = 2 sp — so the test could only run by pooling **NOEC/LOEC** (censored, concentration-spacing-
dependent thresholds that must not be treated as potencies). Verdict: **ECOTOX cannot support a clean
axis-matched weighting test.**

## "Not all endpoints are equal" — the policy this exposes

Three heterogeneity axes: **statistic** (LC50/EC50 potency ✓ vs NOEC/LOEC censored ✗); **effect group**
(apical MOR/GRO/REP/DVP ✓ vs sub-organismal biomarkers ENZ/BCM/GEN/CEL/HRM ✗); **duration/direction**.
DECISION (2026-06-14): **acute LC50|MOR is the primary reported endpoint** (powered, homogeneous, robust);
pooled acute = robustness check; the NOEC/LOEC axis-matched attempt is documented as the ECOTOX limit.

## Required next deliverable — a CURATED APICAL-EC50 database  (PARKED → Paper-2)

**Consolidated 2026-06-14:** a retrieval attempt established this can't be assembled from existing summary
data — local ECOTOX apical-EC50 → M∩R = 2 species; the cross-species chronic literature is NOEC-dominated;
the DEBtox-fitted route gives clean labels but n≈5. So the clean test = whole-budget DEBtox refits across a
designed panel, **parked as Paper-2** (`curated_apical_ecx_plan.md`). The text below is the design spec.

The ECOTOX → axis pMoA *mapping* is a core framework component, so a proper test of its cross-species
consequence (the weighting) must use endpoints that engage the pMoA. ECOTOX can't, so this needs a
hand-curated database:

- **Object:** per (species, chemical) a single coherent **ECx/EC50 on an APICAL endpoint, matched to the
  chemical's DEB axis** — maintenance→respiration/metabolic-rate EC50; growth→somatic-growth EC50;
  reproduction→reproduction-output EC50; assimilation→feeding-rate EC50. Exclude NOEC/LOEC, biomarkers,
  acute lethality. Species ∈ AmP (for κ); body-size covariate `Ww_i` already in hand.
- **Candidate sources (ranked):** (1) DEBtox-fitted dose-response compilations — clean pMoA *by
  construction*, report ECx (e.g. the ESPI2018 / `c7em00328e` set already referenced); (2) OECD standard
  chronic tests — Daphnia reproduction (OECD 211), fish early-life-stage growth (210/212), mysid life-cycle
  — which cover exactly our deep backbone (Daphnia, Pimephales, Oncorhynchus, Americamysis); (3) published
  chronic EC50 / SSD compilations for the clean EDCs (EE2, nonylphenol, BPA) and uncouplers (PCP);
  (4) the ECOTOX apical-EC50 subset as a seed (n≈6–8 per axis — too thin alone).
- **Binding constraint (again):** the shared-species × multi-axis intersection, now on apical EC50 —
  expected small, hence *curated* (hand-assembled), not mined. The existing harness (chemical-centering →
  Δ vs κ, size partial) is reused unchanged once the matrix exists.

## Paper-1 use

Report the acute result as the honest, powered, endpoint-robust **negative** that delimits the validated
claims (rate / relative state / dynamics), exactly as the validation map's grey "not established" band
states — extending the n=310 single-axis null and the n=5 pilot. The curated apical-EC50 axis-matched test
is the dedicated-data complement (Paper-2 / methods), flagged as the open route.

Artifacts: `scripts/extract_ecotox_multimoa.awk`, `scripts/ecotox_multimoa_coverage.jl`,
`scripts/ecotox_multimoa_chem_breakdown.jl`, `scripts/extract_amp_size_proxies.jl`,
`examples/across_axis_weighting_capacity_test.jl`, `…_axismatched.jl`, `…_endpoint_sensitivity.jl`;
data `data/external/ecotox_multimoa_extract.csv`, `…_core_candidates.csv`,
`across_axis_weighting_core_results.csv`, `amp_size_proxies.csv`.
