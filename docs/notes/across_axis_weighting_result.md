# Across-axis capacity-weighting test — result (2026-06-14)

**The central open validation question** (session-5 handover §3.1): does the model's per-species
across-axis weighting predict species × mode-of-action sensitivity *beyond body size*? Answer from a
powered multi-MoA × shared-species test built from the local EPA ECOTOX dump:
**not corroborated — the core test is a null (and wrong-signed); the apparent signal in the larger set
is most parsimoniously a taxonomic target-site artifact, not the energetic weighting.**

## What the model predicts

Operative weights `w = (1/2, κ/4, κ/4, (1−κ)/2)` over (assim, maint, growth, repro), `κ = alpha_axes[3]`
(DEB allocation fraction). Assimilation is fixed → no cross-species signal. The discriminating contrast
is **maintenance** (`w_M = κ/4`) vs **reproduction** (`w_R = (1−κ)/2`); `w_M/w_R = κ/(2(1−κ))` rises
monotonically in κ. So a **higher-κ** species should be **relatively more sensitive to maintenance-axis**
chemicals and less to reproduction-axis chemicals.

## Design (uncalibrated, size-controlled by construction)

`examples/across_axis_weighting_capacity_test.jl`. Per (species, chemical): median log10 LC50/EC50/IC50
(µg/L). **Chemical-center** to remove the chemical main effect; per species `Δ = resid_R − resid_M`
(`Δ>0` ⇒ relatively more maintenance-sensitive). Within-species differencing removes the species-level
overall sensitivity **and body size**. Prediction: **ρ(κ, Δ) > 0**. Tests: permutation Spearman, size
partial (log Ww_i), OLS slope. Panel `data/ecotox_multimoa_panel.csv`: **core** = defensible pMoA
(maintenance: uncouplers PCP/DNP + DEBtox-fitted-M pyridine/aldicarb/organotins; reproduction: estrogenic
EDCs EE2/estradiol/BPA/nonylphenol); **stratum** = AChE pesticides + DDT + PAHs (contested pMoA).

## Result

| panel | n (M∩R) | ρ(κ, Δ) | perm p | partial ρ \| log Ww | OLS slope (t) |
| --- | --- | --- | --- | --- | --- |
| CORE | 27 | **−0.28** | 0.15 (null) | −0.29 | −1.97 (t=−1.97, df=25) |
| CORE+STRATUM | 101 | **−0.20** | 0.043 | −0.17 | −0.99 (t=−2.00, df=99) |

κ–size collinearity is negligible (Spearman ρ(κ, log Ww) = +0.03 core), so the result is genuinely
size-independent — the partial barely moves the estimate.

## Reading (honest)

- **The predicted sign (ρ>0) does not appear.** Core is a non-significant null with a *wrong-signed* point
  estimate; the only "significant" result needs the contested AChE/DDT stratum and is also wrong-signed.
- **The negative pattern is a taxonomic/target-site artifact.** The low-κ species (cladocerans *Daphnia*/
  *Ceriodaphnia*, *Chironomus*, mysids, *Oryzias*) carry positive Δ = relatively maintenance-sensitive —
  but that is arthropod/insect sensitivity to the maintenance-class chemicals (insecticides, uncouplers),
  a **target-site** effect the DEB energetic weighting does not (and should not) capture. The core (no
  AChE) is null; the significance appears only when AChE insecticides enter (stratum) — i.e. the signal is
  the insecticide-sensitivity confound, not corroboration.
- **Bounding caveat — endpoint–axis mismatch.** Pooled *acute* LC50/EC50 does not exercise an EDC's
  *reproductive* pMoA (reproductive effects are chronic/sublethal); the reproduction arm may be acting via
  acute narcosis. So this is "not corroborated, partly expected from endpoint mismatch," **not a hard
  falsification** of the weighting.

## Bottom line + next

Confirms the **scale-free / no-calibration framing**: the across-axis weighting is **not established**, now
shown by a *powered cross-MoA test* (n=27 core / 101) — extending the n=310 single-axis null and the n=5
pilot. **Paper-1 use:** report as the honest, powered negative that delimits the validated claims (rate /
relative state / dynamics), exactly as the validation map's grey "not established" band already states.

**Decisive refinement (before any stronger claim):** repeat with **axis-matched sublethal endpoints** —
maintenance: respiration / growth-rate EC50; reproduction: reproduction EC50 — so each chemical actually
exercises its DEB axis. ECOTOX codes these effect groups (REP/GRO/PHY); the extractor already captures the
`effect` column, so this is an endpoint filter + per-axis endpoint map, not new data. A clean physical-size
covariate is already in hand (`Ww_i`, `scripts/extract_amp_size_proxies.jl`).

Artifacts: `scripts/extract_ecotox_multimoa.awk`, `scripts/ecotox_multimoa_coverage.jl`,
`scripts/ecotox_multimoa_chem_breakdown.jl`, `scripts/extract_amp_size_proxies.jl`,
`examples/across_axis_weighting_capacity_test.jl`; data `data/external/ecotox_multimoa_extract.csv`,
`…_core_candidates.csv`, `across_axis_weighting_core_results.csv`, `amp_size_proxies.csv`.
