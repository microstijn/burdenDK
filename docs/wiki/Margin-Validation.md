# Margin Validation (Scope for Growth + GlobTherm)

The COMADRE validation ([page](COMADRE-External-Validation.md)) corroborated the recovery
curve's **rate endpoints** — the DEB maintenance rate `k_M` (↔ demographic recovery) and the
reproduction rate `R_i` (↔ compensation). It did **not** test the *adaptive margin itself* —
its state, or its erosion under chronic pressure — which is the framework's actual product.
This page covers the two external checks that target the **margin**, at the margin's own
organisational level (**individual energetics → no scale bridge**, unlike COMADRE):

1. **GlobTherm** (thermal tolerance) — a *coherence/bounding* probe of the capacity axis.
2. **Scope for Growth (SFG)** — the real test: does the modelled margin track an *independent,
   same-level* energetic outcome along real contaminant gradients?

---

## 1. GlobTherm — capacity coherence (a bounding result)

A pre-registered probe (`examples/globtherm_validation.jl`, n=664 AmP↔GlobTherm species):
does the AmP recovery-capacity axis carry a general thermal-tolerance signal (broader
`CTmax − CTmin`), beyond body-size and latitude?

Two honest, opposite-pointing findings:
- ✅ **Coherence:** the AmP capacity axis correlates strongly (|ρ| up to **0.45**) with an
  entirely independent physiological dataset (measured CTmax/CTmin) — the AmP extraction
  carries real biological structure, not noise.
- ❌ **General resilience REFUTED:** recovery capacity does **not** predict broader thermal
  tolerance; if anything higher `k_M`/`λ_max` → *narrower* breadth (−0.45 all-taxa, weak
  −0.28\* in the cleaner ectotherm subset). Recovery capacity and thermal tolerance are
  **separate axes** — the COMADRE-validated "recovery capacity" is *specific to demographic
  recovery*, not a universal resilience currency.

**Role here:** GlobTherm does not validate the margin (the margin's currency is energetic
capacity-under-pressure, not thermal limits). It is a *coherence win* + a *bounding rule*:
do not sell the capacity axis as "general resilience." This motivates the SFG test — the
margin's own energetic currency. Detail: `docs/notes/globtherm_validation.md`.

---

## 2. Scope for Growth — the margin's own currency

SFG = energy absorbed − energy respired ≈ capacity beyond maintenance ≈ **the adaptive
margin in energetic terms**, measured at the **same organisational level** (individual
energetics) and **independent** of AmP. Along a contaminant gradient, SFG tests the margin's
core claim directly: *chronic pressure erodes the margin → an energetic outcome declines.*

**Method** (`examples/sfg_margin_validation*.jl`). For each site, per-contaminant tissue
burden → DEB-axis pressures via a documented mode-of-action (pMoA) routing
(hydrocarbons→assimilation, metals→maintenance, organotin/organochlorine→reproduction),
aggregated as a threshold-free **median-normalised relative burden** (`tu = conc/median`);
then the margin-first point API `compute_adaptive_margin_response` → modelled margin `A_t`;
then **rank**-correlate `A_t` against measured SFG across sites. Single species per study, so
the AmP **capacity weighting is held constant** — this tests the *erosion mechanism + MoA
aggregation*, not the capacity weighting (which would need across-species SFG; see Caveats).

### Results — three gradients, a scale-dependent picture

| study | gradient | n | ρ(margin, SFG) | reading |
| --- | --- | --- | --- | --- |
| **Widdows et al. 1995** (North Sea, *M. edulis*) | estuary/regional, hydrocarbon-dominated | 36 | **+0.41 \*** | ✅ tracks SFG (beats naive 0.22; ≈ best single contaminant 0.47) |
| **Widdows et al. 2002** (Irish Sea, *M. edulis*) | basin-scale, hydrocarbon/industrial | 23 | **+0.12** | ◐ right direction, attenuated by scale (beats naive 0.005) |
| **Albentosa et al. 2012** (Iberia, *M. galloprovincialis*) | condition/food-confounded | 39 | **−0.11** | ✗ confound-flipped |

*\* |t|>2.0. Widdows 2002 SFG was figure-digitized (rank only). See per-study notes.*

The three studies cohere into a **scale-dependent, confound-sensitive** story:
- **Strong where tissue burden indexes exposure** — Widdows 1995, a fine-scale
  hydrocarbon-dominated gradient: the modelled margin reproduces the measured energetic
  gradient (ρ=+0.41), the **first direct same-level corroboration of the margin itself**, and
  the MoA routing beats naive equal-weight aggregation (0.41 vs 0.22).
- **Weak at basin scale** — Widdows 2002: same group/method/species, but over the large Irish
  Sea contaminants stop co-varying with any single marker (the *paper's own point*: r≈−0.9 in
  small estuaries → weak at basin scale). The margin still tracks SFG in the right direction
  (+0.12) and still beats the naive load (0.005).
- **Flipped under condition confounding** — Albentosa: the authors show SFG is dominated by
  **condition index** (`SFG~CI` r=−0.62\*\*\*) and **age**, with chemicals adding only ~17% of
  variance. Tissue burden here indexes *food/condition*, not toxic exposure, so the
  burden-driven margin anti-tracks SFG. Our harness reproduces their coefficients
  (`SFG~CI` ρ=−0.67\*\*) and controlling CI+age does not rescue it — a **bounding result**, not
  a tuning failure (routing was held fixed).

**One consistent thread:** in *every* dataset, **metals behave as a positive confound**
(As/Cd/Zn correlate *positively* with SFG — both Widdows 2002 and Albentosa report this). The
MoA routing's job is precisely to keep these off the toxic axis and route the genuine signal
(hydrocarbons) to assimilation; that is why the routed margin consistently beats naive
equal-weight load.

Detail + reproducibility: `docs/notes/sfg_validation_status.md`.

---

## Verdict

The **margin/recovery layer** has external support at two levels now: its *rate endpoints*
(COMADRE: `k_M`, `R_i`) and the *margin state under pressure itself* (SFG). The SFG evidence
is the first **same-level, no-scale-bridge** corroboration of the adaptive margin — strong
where tissue burden indexes exposure (ρ=+0.41), attenuating predictably with spatial scale,
and honestly bounded where condition/food confounds dominate. Modest and qualified, but it
validates *the thing the framework is for* — the margin state under chronic pressure — not
just adjacent rates. The amplification scalar `g`/`F` remains null throughout, consistent with
the [margin-first reframe](Limitations-and-Open-Questions.md).

## Caveats (carried)
- **Capacity weighting untested.** Single-species designs hold the AmP capacity constant;
  testing the κ-rule axis weights needs *across-species* SFG. That data is largely absent —
  the SFG corpus is mussel-dominated, ICES DOME no longer holds SFG, and non-mussel SFG is
  temperature- not contaminant-driven (see `sfg_validation_status.md`).
- **Tissue burden ≠ exposure.** A threshold-free median-normalised relative burden is the
  pressure proxy, so per-contaminant potency is not encoded; where burden tracks
  food/condition rather than exposure, the proxy fails (Albentosa).
- **pMoA routing** is a documented, approximate assignment.
- **Rank statistics throughout** (consistent with the rest of the validation programme);
  Widdows 2002 SFG is figure-digitized (±~1 J/g/h), rank-only.

## Sources
- SFG: [Widdows et al. 1995, *MEPS* 127:131](https://www.int-res.com/abstracts/meps/v127/p131-148/) ·
  [Widdows et al. 2002, *Mar. Environ. Res.* 53:327](https://doi.org/10.1016/S0141-1136(01)00120-9) ·
  [Albentosa et al. 2012, *STOTEN* 435–436:430](https://doi.org/10.1016/j.scitotenv.2012.07.025)
- GlobTherm: [Bennett et al. 2018, *Scientific Data* 5:180022](https://www.nature.com/articles/sdata201822)
- Repo notes: `docs/notes/sfg_validation_status.md`, `docs/notes/globtherm_validation.md`.
