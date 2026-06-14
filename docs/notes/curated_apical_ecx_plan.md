# Curated apical-EC50 database — plan + the power/label tension (2026-06-14)

> **STATUS (consolidated 2026-06-14): PARKED as Paper-2 future work.** A retrieval attempt concluded that
> **assembling a clean apical-EC50 cross-species matrix from existing summary data is not feasible**:
> local ECOTOX apical-EC50 → M∩R = 2 species (`scripts/ecotox_apical_ecx_seed.jl`); the cross-species
> chronic literature (Caldwell EE2 SSD over 26 species; PCP datasets) is **NOEC-dominated** — the
> excluded statistic; the DEBtox-fitted route gives clean labels but only 5 AmP species with ≥2 axes.
> **Conclusion:** the only clean, powered route is **whole-budget DEBtox refits across a designed
> species × MoA panel** (re-analysing raw dose–response–time data, not assembling summary ECx) — a real
> dedicated study. Paper-1 reports the powered ECOTOX **negative** (`across_axis_weighting_result.md`).
> Everything below is the design spec for that future study.

Goal: a clean, powered cross-species test of the **across-axis capacity weighting** (κ-driven
maintenance-vs-reproduction). The acute ECOTOX test was powered (n=27/101) but not clean (acute
endpoint + contested pMoA → not corroborated, robustly). A clean test needs (i) **defensible pMoA axis
labels** and (ii) **cross-species apical sublethal ECx** (NOT NOEC/LOEC — see
`across_axis_weighting_result.md`). The hard finding below: those two requirements pull against each other
in the available data.

## What the DEBtox-fitted source gives (axis labels — gold, but sparse)

`scripts/extract_debtox_fitted_pmoa_species.jl` parses the ESPI2018 compilation
(`c7em00328e1.xlsx`, doi 10.1039/c7em00328e) → `data/external/debtox_fitted_pmoa_species_c7em00328e.csv`:
species-resolved, **cascade-aware** pMoA (Jager 2020: the most-sensitive apical endpoint ≠ pMoA, because
the κ-rule starvation cascade pushes most damage downstream onto reproduction; only a whole-budget DEBtox
fit recovers the true axis).

- **45 AmP-resident (species, compound) fitted pairs across 13 AmP species.**
- **Only 5 species have ≥2 distinct fitted axes** (Daphnia magna [A/G/M/R], C. elegans [A/G/M],
  Folsomia candida [A/M/R], Danio rerio [A/M], Lumbricus rubellus [A/M]) — and fits are
  **assimilation-skewed** (the metals→A result). Maintenance/reproduction fits are 1–3 compounds each.
- The species are DEBtox-study organisms (nematode, springtail, earthworms, snail, Daphnia, mussels,
  urchin, zebrafish), only **Daphnia + Danio** overlap the ECOTOX powered backbone.

**Verdict:** the fitted route supplies clean *labels* but, used as the data itself, gives n≈5 — back to
the pilot. It cannot power a cross-species test alone. (Potencies aren't even in the compilation; they're
in the papers behind the DOI column.)

## The way through — labels (DEBtox) × power (cross-species apical chronic ECx)

Combine: take the **clean fitted axis label** for a chemical from the compilation (or unambiguous
mechanism where no fit exists), then source its **cross-species apical chronic ECx** from richer data.

- **Reproduction arm (richest apical data):** the **EDC literature** — 17α-ethinylestradiol, 17β-estradiol,
  nonylphenol, bisphenol-A — has abundant *cross-species chronic reproduction* ECx (fish, Daphnia). Plus
  fitted-R PAHs (pyrene, fluoranthene, 3,4-dichloroaniline R/H).
- **Maintenance arm:** **uncouplers** PCP, 2,4-DNP (clean mechanism) for chronic growth/respiration ECx;
  plus fitted-M pyridine, aldicarb, atrazine, triphenyltin, zinc(earthworm).
- **Endpoint:** apical sublethal **ECx** — reproduction-output EC50 (R), growth/respiration EC50 (M).
  Exclude NOEC/LOEC, biomarkers, acute lethality.
- **Sources to retrieve (user can pull):** OECD chronic test dossiers (Daphnia reproduction OECD 211, fish
  early-life-stage 210/212 — cover Daphnia/Pimephales/Oncorhynchus/Americamysis), EFSA/ECHA dossiers,
  published chronic SSDs for the EDCs/uncouplers, and the DEBtox papers' fitted ECx/NEC (DOIs in the
  shopping list printed by the parser).

## The deeper caveat (state it, don't paper over it)

The κ-rule cascade also undercuts the *observable*: a maintenance chemical's strongest apical effect is on
reproduction (downstream), so a reproduction-EC50 does not cleanly isolate the reproduction axis. The
weighting's true consequence is on the **integrative margin**, not a single apical trait. So the rigorous
test may ultimately need **whole-budget DEBtox refits across a species panel** (or a population-rate
integrative endpoint à la Billoir 2007), not summary ECx — i.e. a dedicated study (Paper-2+), not a mine.
The curated apical-ECx DB is the best *feasible* approximation and the honest next increment; its limits
are real and should be reported.

## Schema (curated_apical_ecx.csv, when populated)

`species, amp_key, kappa, chemical, cas, axis, pmoa_label_source, endpoint, ecx_type, ecx_value_ugL,
exposure_days, source_ref, doi, confidence, notes` — one apical ECx per (species, chemical); axis from the
fitted compilation/mechanism; reuse the Δ-vs-κ harness (`across_axis_weighting_capacity_test.jl`) unchanged.

## Status / next

Backbone (axis labels) built + AmP-crossed. Next is a **data hunt** (user-assisted): assemble cross-species
chronic reproduction ECx for the EDCs and chronic growth/respiration ECx for the uncouplers — the
best-bet powered, clean M-vs-R contrast. Recommend starting there (EDCs vs uncouplers) before the broader
fitted-paper ECx extraction (which mostly yields assimilation + n≈5).
