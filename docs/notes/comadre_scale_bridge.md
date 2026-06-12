# The individual→population scale bridge (roadmap Part 3)

*Working note, 2026-06-12. The COMADRE validation correlates an **individual-energetic**
model rate (`k_M`, `R_i`, `r_B`, per day, from a DEB energy budget) with a
**population-demographic** quantity (damping ratio, reactivity, attenuation, per
projection interval, from a matrix population model). These live at different
organisational levels and in different units. A correlation is therefore not
self-evidently meaningful — the bridge must be **argued, not assumed**. This note
states the argument and its caveats so the manuscript can cite it rather than gloss it.*

## The apparent mismatch
- **Model side.** `λ(A)` and its bounds are individual physiological rates: `k_M =
  [p_M]/[E_G]` is the somatic maintenance rate constant (1/d), `λ_max = v/L_m` a
  reserve-mobilisation rate (1/d), `R_i` the ultimate reproduction rate (#/d). They
  describe how fast one individual's state variables relax/turn over.
- **Demographic side.** The damping ratio `|λ₁|/|λ₂|` is the asymptotic rate at which a
  *population's stage structure* returns to equilibrium after perturbation; reactivity
  and attenuation are first-timestep *population* transients. These are properties of
  the projection matrix **A**, per census interval.

Per-day individual rate vs per-interval population rate is not a 1:1 identity.

## Why the bridge is principled
The link is exactly what **physiologically-structured population models (PSPMs)** and
**DEB-structured population models** formalise: population dynamics are *derived from*
the individual energy budget.

1. **DEB-IPM (Smallegange, Caswell, et al.).** The Dynamic Energy Budget Integral
   Projection Model builds the population kernel (survival, growth, reproduction of the
   structuring trait) directly from the standard DEB individual model. The matrix /
   kernel entries are functions of the same DEB primitives (`p_Am, κ, v, [p_M], [E_G]`)
   that define `k_M`, `λ_max`, `R_i`. So `matU` (survival/growth) inherits the
   individual maintenance/growth timescale and `matF` (fertility) inherits the
   individual reproduction rate — *by construction of the model class*, not by analogy.
2. **PSPM (de Roos & Persson).** Continuous physiologically-structured models show the
   population return rate (the demographic analogue of the damping ratio) is set by the
   slowest individual-state relaxation timescale. When somatic maintenance/turnover
   (`k_M`) is rate-limiting for individual state, it is rate-limiting for the
   population's return to stable structure → a **monotone** individual→population map.
3. **Reproduction → fertility transient.** Reactivity (max column sum of **Â**) is
   dominated by the `matF` entries; in a DEB-IPM those entries are the individual
   reproduction output. Hence the observed `R_i` ↔ reactivity association is the
   model-class prediction surfacing in independent field matrices (see the per-axis
   note for the partly-mechanical caveat on that specific cell).

## What this licenses, and what it does not
- **Licenses:** a *monotone* (hence rank-correlational) association between an
  individual DEB rate and the matching population-demographic quantity. Rank/Spearman
  is unit-invariant, so the per-day vs per-interval difference and any monotone
  temperature/size rescaling do not invalidate the test — they are exactly what rank
  statistics absorb. This is why every COMADRE test here is rank-based.
- **Does not license:** a quantitative 1:1 equality of rates, or a claim that the
  individual rate is the *only* driver of the population quantity. The bridge is "the
  individual rate is a monotone determinant of the population rate, other things equal",
  with pace-of-life (generation time) and phylogeny as the obvious "other things" —
  which is precisely why the headline statistics are *partial* (gen-controlled) and
  why we tested phylogeny ([Phylogenetic PGLS](comadre_pgls_validation.md)).

## Caveats to state in the manuscript
- **Temperature.** AmP rates are at a reference temperature; COMADRE matrices are at
  field temperatures. Arrhenius rescaling is monotone per species, so rank tests
  survive, but absolute-rate comparisons would not.
- **Rate-limiting assumption.** The PSPM argument requires the focal individual rate to
  be (close to) rate-limiting for the population return; where another process dominates,
  the monotone map weakens — consistent with the modest partial ρ (~0.2–0.3) on the
  maintenance/recovery side vs the strong reproduction→compensation map.
- **Aggregation.** COMADRE matrices vary in structure/dimension; the bridge is to the
  *asymptotic/transient* population descriptors, not to any single matrix entry.

## Sources
- Smallegange, I.M., Caswell, H., et al. — Dynamic Energy Budget Integral Projection
  Model (DEB-IPM). *Methods Ecol. Evol.*
- de Roos, A.M. & Persson, L. — *Population and Community Ecology of Ontogenetic
  Development* (physiologically-structured population models).
- Caswell, H. (2001) *Matrix Population Models* — damping ratio, transient dynamics.
- Stott, Townley, Hodgson (2011) — transient bounds (reactivity/attenuation).
