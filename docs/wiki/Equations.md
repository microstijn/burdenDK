# Model equations

[← Pipeline](Pipeline.md) · next: [Data & parameters →](Data-and-Parameters.md)

All symbols in one place. Source of truth: [`src/deb_axes.jl`](../../src/deb_axes.jl)
(response math), [`src/AmP_Translator.jl`](../../src/AmP_Translator.jl)
(offline capacity mapping), [`src/ecotox_library.jl`](../../src/ecotox_library.jl)
(memory & stress). The manuscript derivations are in `docs/tex/`.

## 1. Capacity mapping (offline, AmP → parameters)

From AmP compound parameters `{p_Am, p_M, κ, v, k_M, E_G}`:

$$ E_m = \frac{p_{Am}}{v}, \qquad A_0 = E_m, \qquad L_m = \frac{\kappa\,p_{Am}}{p_M} $$

$$ \alpha_A = \frac{1}{E_m},\quad \alpha_M = \frac{p_M}{\kappa\,p_{Am}}=\frac{1}{L_m},\quad \alpha_G = \kappa,\quad \alpha_R = 1-\kappa $$

$$ \lambda_{\max} = \frac{v}{L_m}, \qquad \lambda_{\min} = \min\!\left(k_M,\ \lambda_{\max}\right),\quad k_M = \frac{[p_M]}{[E_G]}, \qquad K_A = 0.3\,A_0 $$

> **`λ_min` was re-anchored** to the DEB somatic maintenance rate constant `k_M`
> (it used to be `p_M/A_0 = [p_M]/[E_m]`). The timescale ratio is now the energy
> investment ratio $\lambda_{\max}/\lambda_{\min} = g = [E_G]/(\kappa[E_m])$ instead of
> the artifact $1/\kappa$, so amplification tracks `g`, not κ. The `0.3` in `K_A`
> remains an unjustified constant. Full story:
> [the tex note](../notes/lambda_min_maintenance_rate.tex) ·
> [Limitations §1](Limitations-and-Open-Questions.md).

## 2. Memory (chemical burden `B_t`)

$$ B_t = \rho\,B_{t-1} + (1-\rho)\,K\,C_t,\qquad 0 \le \rho < 1 $$

`ρ` = retention/carryover, `K` = bioaccumulation/internal magnification, `C_t` =
ambient concentration. Analytical steady states (constant, inverse-target,
periodic) are in [`compound_memory_warmup.jl`](../../src/compound_memory_warmup.jl).
Detail: [Compound memory](../compound_memory.md).

## 3. Active stress

$$ x = \max\!\left(0,\ \frac{(B\text{ or }C) - \text{NOEC}}{\text{EC50} - \text{NOEC}}\right) $$

NOEC acts as a no-effect floor; the EC50–NOEC span sets the slope.

## 4. Per-axis impairment and mixtures

Each compound's stress is routed to a DEB axis by its effect code, then bounded
to an impairment. The per-axis bounded shape is

$$ E_a = \frac{x_a}{1 + x_a} \in [0,1) $$

Multiple compounds on the same axis are combined by one **mixture-effect
assumption** (no fitted interactions):

- **TU** (concentration addition within a shared target): $x_{\text{tot}} = \sum_j x_j$
- **IA** (distinct targets): $E_{\text{tot}} = 1 - \prod_j (1 - E_j)$
- **grouped CA-then-IA** (default): TU within an effect-code group, IA across groups.

See [`mixture_aggregation.jl`](../../src/mixture_aggregation.jl) ·
[Mixture-effect models](../mixture_effect_models.md).

## 5. Scalar load `Q_t` and axis weights

$$ Q_t = \sum_{a \in \{A,M,G,R\}} w_a\,E_a,\qquad \sum_a w_a = 1 $$

The default weights are **κ-rule, assimilation-led** (dimensionless, grounded in
the DEB allocation split):

$$ w = \left[\tfrac12,\ \tfrac{\kappa}{4},\ \tfrac{\kappa}{4},\ \tfrac{1-\kappa}{2}\right], \qquad \kappa = \frac{\alpha_G}{\alpha_G + \alpha_R} $$

Assimilation gates the reserve buffer (the margin itself), so it carries the
largest share; the somatic fraction κ is split equally over maintenance and
growth; reproduction carries `1−κ`. The legacy `normalized_alpha_axes` weighting
is retained as a diagnostic only (it gave assimilation ≈ 0 weight; see
[Limitations](Limitations-and-Open-Questions.md)). Code:
[`axis_weights_for_species`](../../src/deb_axes.jl).

## 6. Adaptive margin `A_t`

**Canonical (nondimensional) mode** — the default:

$$ A_t = A_0\,\max(10^{-6},\ 1 - Q_t) $$

**Diagnostic (raw subtraction) mode** — retained for comparison only, numerically
inert because `A0 ≫ Σ α·s`:

$$ A_t = A_0 - \sum_a \alpha_a\,s_a $$

## 7. Restoring force and amplification

$$ \lambda(A) = \lambda_{\min} + (\lambda_{\max}-\lambda_{\min})\,\frac{A_+}{K_A + A_+},\qquad A_+ = \max(A,0) $$

$$ \boxed{\,F_t = \dfrac{\lambda(A_0)}{\lambda(A_t)}\,} $$

`F_t = 1` ⇒ no amplification; `F_t > 1` ⇒ chronic pressure has weakened recovery.

## 8. Why `F ∝ 1/λ` — the two-timescale argument

For an acute event profile `C_event(τ)` on the fast time `τ`, with the slow margin
`A` held fixed, the induced burden obeys

$$ \frac{dy}{d\tau} = -\lambda(A)\,y + \beta\,C_{\text{event}}(\tau),\qquad y(0)=0 $$

Integrating the displacement over the event (a Fubini swap of the time integrals)
gives total response burden $\propto 1/\lambda(A)$. Standardising against pristine
conditions yields the amplification factor in §7 as a pure ratio of restoring
forces — independent of the event's shape. This separation (slow `A`, fast `y`)
is what avoids per-cell ecosystem Jacobians over large rasters.

> **`F` is a derived readout, not the primary output.** The closed form above is a
> *happy analytical accident* of the two-timescale separation. The product is the
> **adaptive-margin state** (§5–§6): relative depletion `Q_t`, absolute margin `A_t`
> (which carries capacity `A_0`), and the axis composition `E_a`. `F` collapses all
> of that into one scalar that is capacity-blind and one-dimensional — useful, but
> secondary. See [Limitations §1](Limitations-and-Open-Questions.md).

## 9. Optional layers

- **Physiological condition memory `Z_t`** (opt-in, off by default):
  [`condition_buffer.jl`](../../src/condition_buffer.jl) integrates a condition
  buffer and can shift the margin via `A_t = A0 + ω_Z·Z − Σ α·s`. Kept strictly
  distinct from chemical memory `B_t`.
- **DEBtox scaled damage `D_t`**: not implemented (by design).

## Symbol table

| Symbol | Meaning |
| --- | --- |
| `C_t` | ambient concentration at time `t` |
| `B_t` | retained internal chemical burden (memory) |
| `ρ`, `K` | retention factor, bioaccumulation/magnification factor |
| `x`, `x_a` | active stress (per axis) |
| `E_a` | bounded per-axis impairment ∈ [0,1) |
| `w_a` | axis weight (κ-rule) |
| `Q_t` | scalar impairment load ∈ [0,1] |
| `A_0`, `A_t` | baseline / current adaptive margin |
| `α_a` | per-axis sensitivity (assimilation, maintenance, growth, reproduction) |
| `κ` | DEB somatic allocation fraction = `α_G/(α_G+α_R)` |
| `λ_min`, `λ_max`, `K_A` | restoring-force bounds and half-saturation |
| `λ(A)` | restoring force at margin `A` |
| `F_t` | amplification factor `λ(A_0)/λ(A_t)` |
