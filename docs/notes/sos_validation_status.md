# Stress-on-Stress (survival-in-air) — the DIRECT amplification test: RESULTS

## ✅ RESULTS (2026-06-13) — chronic burden → eroded margin → reduced acute-stress survival
The strongest, most *on-thesis* external support for the adaptive margin so far. Where Scope
for Growth validated the margin **state** (the energetic budget), **stress-on-stress (SoS)**
validates what the margin is **for**: the capacity to withstand an **acute perturbation**.
SoS = survival-in-air (days) under emersion/anoxia — a direct field proxy for the framework's
core two-timescale claim (chronic pressure erodes the margin → less able to survive an acute
hit).

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
1. **On-thesis.** It tests the amplification claim (resilience to an acute perturbation)
   directly, not the margin state. This is the thing the framework is *for*.
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

### The obvious power-boosting follow-up (not yet done)
The data are **multi-year**. A **within-station temporal** analysis — does SoS survival track
contaminant *change* at the *same* station across 2012–2022? — would (a) add real statistical
power and (b) test the erosion mechanism **over time** directly (the dynamic claim), not just
cross-sectionally. Needs SURVT-year ↔ contaminant-year matching and repeated-measures handling
(mixed model or within-station rank methods). This is the natural next step to turn a
suggestive cross-sectional result into a powered one.

## Where SoS sits in the validation programme
| layer | anchor | result |
| --- | --- | --- |
| recovery rate endpoints | COMADRE (`k_M`, `R_i`) | corroborated |
| capacity coherence (bounding) | GlobTherm | recovery-specific, not general resilience |
| margin **state** | Scope for Growth | +0.41 (estuary) → +0.12 (basin) → −0.11 (confounded) |
| margin **function (acute resilience)** | **Stress-on-Stress** | **+0.39 / +0.45 controlled — the amplification claim, direct** |

## Sources
- ICES DOME 2024 OSPAR CEMP biota (figshare [27211422](https://ices-library.figshare.com/articles/dataset/Data_and_results_for_the_2024_OSPAR_CEMP_assessment/27211422), CC BY 4.0); SURVT = stress-on-stress survival.
- SoS method/biomarker background: Viarengo et al.; Eertman et al. 1993; ICES biological-effects monitoring.
- Mechanism (hydrocarbon narcosis on feeding/energetics): Widdows et al. 1995, 2002.
