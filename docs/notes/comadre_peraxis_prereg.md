# Pre-registration — per-axis DEB rates ↔ demographic resilience components (Idea B)

*Working note, 2026-06-12. **Written BEFORE computing the correlation matrix**, to
avoid fishing (roadmap Part 2 explicitly requires this). Records the axis↔component
mapping, the metrics, the controls, and the confirm/refute criteria. The analysis
(`examples/comadre_peraxis_validation.jl`) must not alter this mapping after seeing
results; any deviation is logged in the results writeup.*

## The question
The scalar validation (`k_M` ↔ recovery) used one model rate against one demographic
quantity. The margin-first reframe claims the **per-axis** structure carries
information. Idea B tests whether **different DEB process rates predict different
demographic-resilience components** — a structured (diagonal) correspondence, not a
single correlation.

## Model side — DEB process rates (AmP, reference T), per species
Extracted (side-extraction, `scripts/extract_amp_reproduction_rates.jl`) from
`allStat.mat`; coverage >99.6% of 7335 species.

| axis | DEB rate | symbol | units |
| --- | --- | --- | --- |
| maintenance | somatic maintenance rate constant | `k_M` | 1/d |
| reproduction | ultimate reproduction rate | `R_i` | #/d |
| growth | von Bertalanffy growth rate | `r_B` | 1/d |

(Secondary/intensive: `kap_R` reproduction efficiency; `k_J` maturity maintenance.)
**Known confound:** `R_i` is absolute fecundity → strongly body-size dependent. It is
therefore tested as a *partial* (controlling generation time, and body mass `Ww_i`
as a secondary control), never raw. Rank (Spearman) is the primary statistic —
Idea A showed these signals are monotone-but-not-log-linear; log-linear is reported
alongside for transparency.

## Demographic side — resilience components (COMADRE `matA`)
Standardised matrix Â = A / λ₁ (removes asymptotic growth). Stott et al. (2011)
first-timestep bounds (over all initial structures — avoids the initial-vector
sensitivity the roadmap warns about):

| component | metric | definition | direction |
| --- | --- | --- | --- |
| recovery | damping ratio | \|λ₁\|/\|λ₂\| | higher = faster return (DONE) |
| compensation | reactivity (P̄₁) | max column sum of Â | higher = more amplification |
| resistance | attenuation (P̲₁) | min column sum of Â | higher (→1) = resists decline |

Aggregated per species as mean over wild/unmanipulated matrices (same filter as the
damping-ratio extraction). Stored as log₁₀.

## Pre-registered hypothesis (the diagonal)
Partial correlations **controlling generation time** (all components and rates are
pace-of-life-loaded; the test is residual-after-pace). Predicted sign in brackets.

| model rate ↓ \ component → | resistance | compensation | recovery |
| --- | --- | --- | --- |
| `k_M` (maintenance) | (−, fast↓resist) | | **(+) strongest** |
| `R_i` (reproduction) | | **(+) strongest** | |
| `r_B` (growth) | | | (+) |

**Core predictions (must hold to confirm):**
1. `R_i` is the strongest model-rate partial predictor of **compensation**, and `R_i`
   predicts compensation more strongly than it predicts recovery or resistance.
2. `k_M` remains the strongest model-rate partial predictor of **recovery**
   (replicates the established scalar result on the larger per-axis table).
3. **Diagonal dominance:** in the rate×component partial-correlation matrix, the
   largest |ρ| in each *column* sits on the hypothesised diagonal cell.

**Confirm (strong):** diagonal dominance survives the generation-time control (and,
secondarily, body mass + the taxonomic-Order proxy). This would be genuinely novel
multi-dimensional external support for the margin *state*, not a scalar.

**Refute:** the matrix is diffuse / off-diagonal — the axes do not map to distinct
demographic components, so the per-axis structure is not externally distinguishable.
Important to know before the manuscript leans on it.

## Caveats fixed in advance
- `R_i` size-confound → always partial on generation time (+ body mass secondary).
- Transient bounds need primitive, non-negative matrices → reuse the existing
  matrix-quality filter; document n.
- The mapping above is a *hypothesis*; this file freezes it pre-computation.
