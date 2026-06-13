# Water-Quality Coupling & pMoA Stressor Routing

How the framework connects to **monthly water-quality models** (e.g. DynQual: BOD, TDS, fecal
coliform, …) — and the honest scope of that coupling. This is newer than the validation work; the
*engine* is validated, the *input routing* is a declared mechanistic assumption.

> **One-line status.** The margin engine, the saturating response, and the recovery dynamics are
> externally validated. The step that maps an aggregate water-quality stressor to a DEB axis is a
> **declared physiological-mode-of-action (pMoA) assignment** — not a tuning knob, not yet externally
> validated for the aggregate stressors. Treat water-quality outputs as **relative, mechanistically
> structured vulnerability indices**, defensible by ranking-stability, not as absolute risk.

---

## 1. The input-type problem (why this needs care)

The validated chain wants **compound-resolved tissue burden**: water concentration of a named CAS →
bioaccumulation recurrence → CAS mode-of-action routing → saturating impairment `E=x/(1+x)` → margin.
This is what the SFG/SoS/Viarengo validations exercise.

But monthly water-quality models mostly emit **aggregate stressor indices** — BOD, salinity/TDS,
fecal coliform, microplastics, nutrients — *not* per-CAS concentrations. Those aren't tissue burdens
of anything, so **no validated routing exists for them**; any "BOD → maintenance" mapping is a
judgment call. This is an **input-type mismatch**, not a code-quality issue: aggregate models force a
mapping step the compound-resolved validation never covered.

**Two regimes, two honesties:**
- **Compound-resolved input** (a fate model emitting CAS concentrations) → inherits the full validated
  chain.
- **Aggregate-index input** (DynQual-type) → the stressor→axis step is a **declared assumption**;
  restrict claims to relative rankings and sensitivity-test the mapping.

---

## 2. The fix: pMoA routing, not tuning weights

An early DynQual prototype mapped proxies to axes with free `w_A/w_M/w_G/w_R` weight vectors — exactly
the kind of tuning knob the project's *no-arbitrary-knobs* invariant forbids. The replacement assigns
each stressor to **the metabolic process it physiologically taxes** (DEBtox pMoA), as a falsifiable
mechanism. Same epistemic category as the mixture rules: *assumptions, not fits.*

**Table:** `data/pMoA_Stressor_Routing.csv` (loader `src/pmoa_stressor_routing.jl`; basis
`docs/notes/pmoa_stressor_routing_basis.md`). 15 stressor classes → primary (± secondary) DEB axis,
tissue-memory `ρ`, fitted-pMoA evidence, and references.

| stressor | primary axis | note |
| --- | --- | --- |
| microplastics, suspended sediment | **assimilation** | feeding/gut impairment |
| BOD/hypoxia, salinity, pH, ammonia, pathogens | **maintenance** | metabolic/osmotic/immune cost |
| heavy metals | **maintenance** (field) / **assimilation** (fitted) | dual route — see §3 |
| hydrocarbons/PAH | **assimilation** | narcosis/feeding |
| endocrine (PCB/DDT/TBT) | **reproduction** | (fitted organotins are mixed — caution) |
| temperature | **none — rate modifier** | Arrhenius on `k_M`/`λ`, *never* a pressure |

Two structural rules: **temperature is the clock, not a pressure** (it rescales every rate, including
the erosion timescale `1/λ`); and **`ρ` (tissue memory) applies only to accumulating stressors** —
ambient stressors (salinity, pH) get `ρ=0` but still erode the margin over months via the slow
`λ`-dynamics.

---

## 3. The `heavy_metals` dual route (the crux — empirically settled)

The fitted DEBtox pMoAs put **metals on *assimilation*** (gut/gill damage → reduced uptake; cadmium
A×6, copper A×3, mercury A, uranium A in the compilation). But the validated field routing puts
**metals on *maintenance*** — deliberately, because in mussel field data metal tissue burden is a
**positive confound** (it tracks food/condition, not exposure: *tissue burden ≠ exposure*).

These were tested against each other (re-routing metals→assimilation in all four field harnesses):

| dataset | metals→maintenance (baseline) | metals→assimilation (variant) |
| --- | --- | --- |
| Widdows 1995 SFG | **+0.41\*** | +0.30 (loses sig.) |
| Widdows 2002 SFG | +0.12 | **−0.14** (sign flips) |
| Albentosa 2012 SFG | −0.11 | −0.21 |
| SoS DOME (raw / +condition) | **+0.39 / +0.45** | +0.27 / +0.26 |

**metals→assimilation degrades every anchor.** Resolution: **both routes stand, chosen by how the
burden is obtained** —
- **field measured tissue burden → maintenance** (confound-control; validated, *keep it*);
- **controlled / exposure-derived burden → assimilation** (the fitted pMoA; valid where burden =
  exposure, e.g. a coupling driven by *modelled* water concentration).

`src/mode_of_action.jl` and the field harnesses are **unchanged** — the validation depends on the
maintenance route. Detail: `docs/notes/pmoa_evidence_digest.md` §4-RESULT.

---

## 4. What did NOT work: inferring pMoA from ECOTOX effect codes

A tempting in-house shortcut — read a stressor's pMoA off which ECOTOX effect endpoint (growth /
reproduction / feeding / physiology) responds at the lowest concentration — was prototyped and is a
**documented negative** (`scripts/ecotox_pmoa_inference.*`,
`data/external/ecotox_pmoa_inference_results.txt`). It does **not** work: the DEB **κ-rule cascade**
makes *reproduction* the most-sensitive endpoint for almost any pMoA (reproduction is funded by what's
left after maintenance + growth), so 4/5 test compounds "infer" reproduction regardless of mechanism;
the comparison is also taxon-confounded. **Most-sensitive endpoint ≠ pMoA.** Use *fitted* DEBtox
pMoAs instead.

---

## 5. Evidence base & status

- **Fitted pMoAs:** `data/external/debtox_fitted_pmoa_c7em00328e.csv` (28 compounds, from the ESPI 2018
  compilation) populate the chemical rows. Caveat: pMoA is **species/study-dependent** (the same
  compound fits as A, A/M, *and* G) — carry the modal pMoA + spread, not a hard label.
- **Framework:** Jager 2020 (the four pMoAs + the starvation cascade that explains §4).
- **To gather** (shopping list in `docs/notes/pmoa_evidence_to_gather.md`): EnviroTox consensus MOA
  (bulk CAS→class), and the aggregate-stressor ecophysiology (microplastics, hypoxia, salinity, …)
  for the rows with no DEBtox fit.

**Defensibility plan:** the routing is a declared mechanism; validate it by **sensitivity-testing** —
perturb each axis/split/`ρ` over its plausible alternatives and show the site/month vulnerability
**rankings are stable**. Given the model is rank-robust everywhere else, they almost certainly are,
and that stability *is* the defensibility statement for the water-quality use.

## Sources / repo
- `docs/notes/pmoa_evidence_digest.md`, `pmoa_stressor_routing_basis.md`, `pmoa_evidence_to_gather.md`.
- `data/pMoA_Stressor_Routing.csv`, `src/pmoa_stressor_routing.jl`.
