# Mixture-Effect Models in TwoTimescaleResilience

The `TwoTimescaleResilience` framework incorporates a mixture-effect layer designed to combine already computed compound-level burdens into bounded fractional axis impairments ($E_a$). 

**Important Conceptual Note:** This framework is DEB-informed, but it is **not a full DEBtox or DEBkiss implementation**.
- It does not implement scaled damage ($D$), damage repair dynamics, $k_d$, GUTS survival, starvation/food-limitation modules, DEBkiss growth/reproduction ODEs, or specific $\kappa$ allocation rules.
- The four axes—assimilation, maintenance, growth, and reproduction—are used simply as coarse, DEBtox-compatible physiological process axes.
- The mixture-effect layer operates **before** the adaptive margin calculation, outputting bounded fractional axis impairments.
- The response layer consumes $E_a$ and remains:
  $$ Q_t = \sum_a w_a E_a $$
  $$ A_t = A_0 \max(10^{-6}, 1 - Q_t) $$
  $$ \lambda_t = \lambda(A_t) $$
  $$ F_t = \frac{\lambda(A_0)}{\lambda(A_t)} $$
- The legacy "raw" response mode remains completely unchanged and independent of this layer.
- The EC50-anchoring formula used within this layer remains unchanged.

These are **axis-level mixture-effect models**, not fitted interaction models. They do **not** implement synergism or antagonism, and they do not rely on fitted interaction coefficients or arbitrary interaction matrices.

## Supported Models

The framework supports three mixture-effect models for aggregating compound-level burdens $x_{j,a}$ (compound $j$ on axis $a$):

### 1. `axis_toxic_unit_sum`
This model computes the exact same-axis sum of active stress units, functioning as a toxic-unit summation prior to nonlinear response conversion.
$$ X_a = \sum_j x_{j,a} $$
$$ E_a = \frac{X_a}{1 + X_a} $$

### 2. `independent_action_axis_effects`
This model assumes that each compound acts on the axis entirely independently. It first computes the fractional impairment for each compound, then combines them according to independent action probability logic.
$$ E_{j,a} = \frac{x_{j,a}}{1 + x_{j,a}} $$
$$ E_a = 1 - \prod_j (1 - E_{j,a}) $$

### 3. `grouped_ca_then_ia_axis_effects`
This is a hybrid grouped mixture-effect model that applies concentration addition within defined effect groups, and independent action across groups.
**Grouping Rule:** The current implementation groups by DEB axis + `effect_code`. If `effect_code` is missing, it is assigned to the string `"unknown_effect_code"`.
Note: `effect_code` is used as an empirical grouping proxy, not a confirmed physiological mode of action.

Within effect group $g$ on axis $a$:
$$ X_{a,g} = \sum_{j \in g} x_{j,a} $$
$$ E_{a,g} = \frac{X_{a,g}}{1 + X_{a,g}} $$

Across groups:
$$ E_a = 1 - \prod_g (1 - E_{a,g}) $$

---

**Distinction:**
- `mixture_method`: Refers to the legacy aggregation label (e.g. `additive_axis_burden`, `axis_toxic_unit_sum`).
- `mixture_effect_model`: Refers specifically to the new axis-level effect-combination models (`axis_toxic_unit_sum`, `independent_action_axis_effects`, `grouped_ca_then_ia_axis_effects`).