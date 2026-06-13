# Stress-on-Stress (survival-in-air) — the margin↔acute-resilience map: RESULTS

## ✅ RESULTS (2026-06-13) — accumulated burden → (static) modeled margin → acute-stress survival
The strongest external support for the adaptive margin so far, and the closest to the
framework's *acute-resilience* claim. Where Scope for Growth validated the margin **state**
(the energetic budget), **stress-on-stress (SoS)** probes what the margin is **for**: the
capacity to withstand an **acute perturbation**. SoS = survival-in-air (days) under
emersion/anoxia.

> **Scope of the test — read this.** This is a **static, cross-sectional map**, not a
> simulation of the two-timescale dynamics. We take the *already-accumulated* field tissue
> burden and pass it through the instantaneous point API
> (`compute_adaptive_margin_response`: burden → per-axis impairment → `A_t = A0·(1−Q)`), then
> correlate `A_t` with measured survival across sites. It does **not** run the dynamics: the
> `B_t` accumulation memory, the slow-timescale erosion ODE, and the amplification of an acute
> pulse are **not** invoked (those functions exist —`simulate_deb_axis_response`,
> `pulse_deb_axes_timeseries`— but the harness never calls them). So this validates the
> **pressure→margin→outcome mapping** (response-curve shape, MoA routing, AmP capacity
> weighting) — the *cross-sectional shadow* of the two-timescale story — **not** the
> accumulate→erode→amplify *sequence*. The dynamic claim is tested only weakly by the
> within-station temporal pass below, and remains essentially unvalidated (see end).

**Data — the multi-station, exposure-paired, QA'd dataset SFG never had.** ICES DOME 2024
OSPAR CEMP (figshare 27211422, CC BY 4.0), *Mytilus edulis*, **17 UK stations** (one SURVT
station has no contaminants), 2012–2022. Station-level medians via
`scripts/extract_dome_sos.jl` → `data/external/sos_dome_ukcemp.csv`; test in
`examples/sos_margin_validation_dome.jl`. Real pressure gradient (PAH ~100-fold, Pb 15×, PCB
14× across stations) and real outcome spread (median survival 6–13.5 d). Mussels are
standard-size (length 4.4–5.9 cm) so the size confound is naturally narrow. Same species as
the validated Widdows 1995 SFG result → coherent.

**Pre-registered (sign frozen before running):** modeled margin `A_t` (contaminant burden →
MoA routing → margin) correlates **positively** with SoS survival; MoA-routed margin beats
naive load; survives body-size control.

### Result
| test | ρ | n |
| --- | --- | --- |
| **survival ~ margin `A_t`** | **+0.39** | 17 |
| survival ~ margin \| length | +0.40 | 16 |
| **survival ~ margin \| length + dry-weight (condition)** | **+0.45** | 16 |

| axis diagnostic (negative = behaves as toxicant) | ρ |
| --- | --- |
| `p_assim` (2–3 ring PAH — toxic hydrocarbons) | **−0.43** |
| `p_maint` (metals Cd/Cu/Hg/Pb/Zn) | +0.09 (near-null) |
| `p_repro` (PCB) | **−0.48 \*** |

Baselines: best single axis (PAH) |ρ|=0.43; **naive mean toxic-unit |ρ|=0.32**; model margin
|ρ|=0.39 (and 0.45 confound-controlled). Empirical: `survival~PAH` −0.43, `~Cu` −0.41,
`~PCB(SCB7)` −0.83\*\* (n=6, striking but small), metals otherwise weak/mixed (`Hg` +0.31,
`Zn` +0.25 — the familiar weak metal positive-confound).

### Reading — why this is *stronger* than the SFG line
1. **Closest to the acute-resilience claim.** The outcome *is* resilience to an acute
   perturbation (not the energetic budget SFG measures) — so the static margin→survival map is
   the cross-sectional shadow of the amplification story. (It does not *simulate* the
   accumulate→erode→amplify sequence — see the scope note above.)
2. **Confound control STRENGTHENS it (0.39 → 0.45).** The decisive contrast with Albentosa:
   there, controlling condition could not rescue the margin (the confound *was* the signal);
   here, partialling size + condition makes the margin signal *stronger* — strong evidence
   this is genuine margin erosion, not a health-proxy artifact. The metal confound that
   dominated SFG (As/Zn/Cd positive) is near-dead here (`p_maint` +0.09).
3. **Mechanism is the validated one.** The PAH/assimilation axis carries the toxic signal
   (−0.43), exactly Widdows' hydrocarbon-narcosis mechanism; PCB/reproduction also negative.
4. **Structure earns its keep.** Routed margin (0.39) > naive equal-weight load (0.32).
5. **Clean data.** Multi-station, exposure-paired, QA'd, open — no figure-digitizing, no
   per-paper grind, no scale/confound apology.

### Honest caveats
- **Power.** n=17 stations → the *two-sided* test is n.s. (p≈0.1). Under the **pre-registered
  one-sided** prediction the size-controlled result is marginally significant (t≈1.9, p≈0.04
  one-sided). Suggestive-to-moderate, not a slam dunk — but positive, confound-robust, and
  mechanistically coherent across every axis.
- **Station-level cross-sectional.** SoS (2012–2022) and contaminants are pooled to station
  medians (years pooled) — station-typical pressure vs station-typical survival.
- **Single species** (*M. edulis*) — the AmP **capacity weighting** is still held constant
  (untestable without across-species data; see `sfg_validation_status.md`).
- PCB axis rests partly on a sparse `SCB7` (nd→neutral); the PAH/metal axes are fully populated.

### Temporal analysis (2026-06-13) — DONE: static claim robust, dynamic claim underpowered
Built a station-YEAR panel (`scripts/extract_dome_sos_yearly.jl` → `sos_dome_ukcemp_yearly.csv`,
nearest-year ±2 contaminant matching since burden changes slowly; one Cu≈1046 mg/kg outlier
station-year QC-dropped) → `examples/sos_temporal_validation_dome.jl`. 36 usable station-years.
Three readouts:

| readout | ρ | note |
| --- | --- | --- |
| **(A) within-station fixed effects** (de-mean by station; the clean dynamic test) | **+0.15** | n.s. (11 stations, 30 st-yr, df=19) — erosion *over time* underpowered |
| **(B) pooled station-years** (n=36; pseudoreplicated) | **+0.28** | PAH axis −0.33\*, PCB axis −0.34\* (toxic axes significant) |
| **(C) station-level, QC-cleaned** (n=17) | **+0.62 \*\*** | nearest-year + outlier-drop *strengthens* vs the +0.39 raw-median |

**Honest reading.** The **within-station** design — which removes *every* fixed between-station
confound (size regime, population, hydrography), the cleanest control possible — is **positive
but not significant** (+0.15). The panel is too thin for the dynamic test: 11 multi-year
stations, mostly 2-year spans, and tissue burden is fairly stable within a station over these
windows, so there is little within-station pressure variation to drive a year-to-year signal.
So the **dynamic erosion-over-time claim is not established with power** here — directionally
consistent, no more.

What the temporal pass *did* establish: (i) the **cross-sectional** margin↔acute-resilience
link is **robustly positive across aggregations** (+0.39 raw-median → +0.62 QC-cleaned/
nearest-year), magnitude specification-sensitive but sign-stable (the programme's recurring
"rank-robust, magnitude-sensitive" pattern, cf. COMADRE `k_M`); and (ii) in the larger panel
the **toxic axes (PAH, PCB) carry significant signal** (−0.33\*, −0.34\*) while metals remain
the weak positive confound. Net: the **static** margin→acute-resilience result is solid; the
**dynamic** (temporal erosion) claim needs longer station time series than DOME currently
offers (and ideally stations with a real within-station contaminant trend).

## Where SoS sits in the validation programme
| layer | anchor | result |
| --- | --- | --- |
| recovery rate endpoints | COMADRE (`k_M`, `R_i`) | corroborated |
| capacity coherence (bounding) | GlobTherm | recovery-specific, not general resilience |
| margin **state** | Scope for Growth | +0.41 (estuary) → +0.12 (basin) → −0.11 (confounded) |
| margin **function (acute resilience), static map** | **Stress-on-Stress** | **+0.39 / +0.45 controlled — burden→margin→acute-survival (not the dynamics)** |
| margin **dynamics** (accumulate→erode→amplify) | within-station temporal (proxy) | **+0.15 n.s. — essentially unvalidated** |

## Sources
- ICES DOME 2024 OSPAR CEMP biota (figshare [27211422](https://ices-library.figshare.com/articles/dataset/Data_and_results_for_the_2024_OSPAR_CEMP_assessment/27211422), CC BY 4.0); SURVT = stress-on-stress survival.
- SoS method/biomarker background: Viarengo et al.; Eertman et al. 1993; ICES biological-effects monitoring.
- Mechanism (hydrocarbon narcosis on feeding/energetics): Widdows et al. 1995, 2002.
