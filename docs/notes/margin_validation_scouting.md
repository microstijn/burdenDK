# External validation of the *adaptive margin* — scouting (#4)

*Working note, 2026-06-12. The COMADRE work (#3) validated the model's **rate
constants** (`k_M`↔recovery, `R_i`↔compensation) — the *endpoints* of the margin→
restoring-force curve — but **not the adaptive margin itself**: the state `A_t`, its
erosion under chronic pressure, and the curve shape were never tested (COMADRE has no
chronic-pressure axis). This note scouts data that could validate the margin's actual
claim. Companion to `external_anchor_scouting.md` (which targeted recovery capacity).*

## The claim that needs an anchor
The adaptive margin is a **state**, and its load-bearing, still-untested assertion is:

> chronic pressure erodes the margin (`A_t = A₀ − Σ αₐ sₐ`) → the restoring force
> `λ(A_t)` falls → a measurable energetic/fitness outcome declines.

So the right anchor is **a chronic-pressure gradient + an independent outcome**, not the
cross-species baseline demography COMADRE provides.

## Two boundaries that decide usability
1. **Circularity.** The model takes ECOTOX `EC50/NOEC` as the *pressure input*
   (`src/ecotox_library.jl`). Therefore we must **not** validate against `EC50/NOEC`
   itself — that is circular. The outcome must be an *independent* measurement (an
   energetic flux, growth, reproduction, or survival), not the toxicity anchor that
   sets `sₐ`.
2. **Level.** COMADRE forced an individual→population **scale bridge** (assumed, not
   measured). An anchor measured at the **individual-energetic level** removes that
   assumption entirely — a strict improvement.

## Candidate anchors (verified reachable), ranked

### 1. Scope for Growth (SFG) along contamination gradients — *best concept match* — ✅ DONE (2026-06-13)
**Run** (`sfg_validation_status.md`, Widdows et al. 1995, 36 North Sea *Mytilus edulis*
sites). The **modeled adaptive margin tracks measured SFG** (ρ=+0.405\*) — the first
direct, same-level (no scale bridge) external corroboration of the *margin itself*. The
MoA-routed margin beats naive aggregation (0.41 vs 0.22) but not the single best
contaminant (0.47) — expected, since one species holds the capacity weighting constant.
Erosion mechanism corroborated; capacity dimension needs an across-species SFG set next.

SFG = energy absorbed − energy respired = net energy available for growth/reproduction
(Widdows). **This is, almost literally, an energetic adaptive margin** — capacity
beyond maintenance — and it has been measured in mussels (*Mytilus edulis*,
*M. galloprovincialis*, all in AmP) along pollution gradients for decades (Mussel
Watch / BEQUALM). A published **SFG-vs-DEB comparison for *Mytilus edulis*** already
works out the DEB↔SFG mapping.
- **Why it's the right test:** same organisational level (no scale bridge); SFG is an
  *independent physiological measurement*, not AmP-derived and not the EC50; it is a
  *gradient* (pressure → erosion), so it tests the margin's actual mechanism.
- **Test design:** for a mussel multi-site study, predict the margin from AmP capacity
  + the site contaminant pressure, and rank-correlate against measured SFG across sites.
  A within-species, across-sites proof-of-concept (not a broad cross-species sweep).
- **Cost/caveat:** SFG data live in papers/tables, not one download — a small
  **literature-assembly** step (pick 1–2 good multi-site mussel studies). Narrow taxa
  (mussels), so a focused validation, not a sweep.

### 2. DEBtox / GUTS sublethal time-series — *most downloadable, partly circular*
Chronic-toxicant growth/reproduction (DEBtox/DEBkiss) and survival (openGUTS) time-
series for workhorse species (*Daphnia magna* is in AmP). Directly tests erosion →
life-history outcome. openGUTS ships example datasets in a documented format.
- **Caveat:** shares the DEB framework, so partial circularity — but the *toxicity
  response* data is empirical and independent of the AmP **capacity** fit, and our
  margin is a *reduced form*, not full DEBtox. A defensible, accessible second test.

### 3. GlobTherm thermal tolerance — *cheap, but tests a different facet* — ✅ DONE (2026-06-13)
**Run** (`globtherm_validation.md`, n=664): the AmP capacity axis is strongly,
coherently related to independent thermal data (|ρ| up to 0.45 — *not noise*), but the
pre-registered general-resilience hypothesis is **refuted** — recovery capacity does
not predict broader thermal tolerance (sign is opposite/weak). **Bounding result:** the
capacity axis is recovery-*specific*, not a universal resilience currency; thermal
tolerance is a separate axis. Confirms GlobTherm doesn't validate the margin — on to SFG.

2,133 species CTmax/CTmin, one CSV download (Dryad doi:10.5061/dryad.1cv08), broad taxa
overlapping AmP. AmP `allStat` carries Arrhenius thermal parameters, so DEB predicts a
thermal window.
- **Caveat:** validates the *thermal* AmP parameters / a capacity boundary — **not the
  margin's erosion dynamic** (which uses the `A₀/α/λ` capacity, not the thermal params).
  Useful as an independent *capacity* check (the `A₀` leg the `F`-cancellation hid), but
  tangential to the margin itself. Runnable immediately if a quick external signal is
  wanted.

## Recommendation
**Pursue Scope for Growth as the primary anchor.** It is the only option that tests the
margin's actual claim (energetic margin ↓ under chronic pressure) at the margin's own
level, with an independent measurement, no scale bridge, and an existing DEB↔SFG
bridge. The cost is a small literature-assembly step rather than a one-click download.
Use **DEBtox/GUTS** as the pragmatic, downloadable second test (flagging the partial
circularity), and **GlobTherm** only as a cheap orthogonal capacity check.

**First concrete step (cheap, no harness yet):** pin one well-cited mussel SFG-gradient
study whose species is in AmP (e.g. Widdows et al. Venice Lagoon for
*M. galloprovincialis*; or a transplant-gradient study), extract the per-site SFG +
contaminant load, and draft the margin-prediction design against it. Then decide.

## The decision for the PI / Kooijman (per project norm — don't build to a circular anchor)
> *Is Scope for Growth along a contamination gradient an acceptable independent anchor
> for the adaptive margin (margin ↓ ⇒ SFG ↓, controlling body size/temperature) — given
> SFG is an energetic flux and the margin a capacity state — or does validating the
> margin require direct chronic-then-acute recovery data (the `F` claim, which already
> failed against COMADRE)?*

## Sources
- Scope for Growth — [Physiological energetics of *Mytilus edulis*: SFG (MEPS)](https://www.int-res.com/articles/meps/46/m046p113.pdf) · [SFG of mussels from the Venice Lagoon (Mar. Poll. Bull.)](https://www.sciencedirect.com/science/article/abs/pii/0141113696000037) · [SFG vs DEB models for *Mytilus edulis* (J. Sea Res.)](https://www.sciencedirect.com/science/article/abs/pii/S1385110111000438)
- DEBtox / GUTS — [openGUTS](https://openguts.info/) · [Revisiting simplified DEBtox models (Ecol. Modelling)](https://www.sciencedirect.com/science/article/abs/pii/S0304380019304120) · [DEBtox (Wikipedia)](https://en.wikipedia.org/wiki/DEBtox)
- Thermal tolerance — [GlobTherm (Scientific Data)](https://www.nature.com/articles/sdata201822) · [Dryad dataset doi:10.5061/dryad.1cv08](https://datadryad.org/dataset/doi:10.5061/dryad.1cv08)
