# TwoTimescaleResilience.jl Math Audit

## 1. Exposure and Stressor Handling

### 1.1 Exposure Filter Application
- **Concept Name:** Exposure Filter Application
- **Reconstructed Equation:**
  $$E_i = h_i \cdot R_i$$
  where $R_i$ is the raw raster/background value for stressor $i$, and $h_i$ is the filter multiplier.
- **File Path & Line Number:** `src/exposure_filters.jl:10`
- **Julia Implementation:**
  ```julia
  function apply_exposure_filter(values::AbstractVector{<:Real}, filter::ExposureFilter)
      return filter.multipliers .* values
  end
  ```
- **Notes/Observations:** Simple element-wise scaling. Serves as a contact, habit, or use filter $H$. Default aquatic filter multipliers are `1.0`. Human filters have placeholder values `[0.10, 0.25, 0.50, 0.60, 0.35, 0.40, 0.30]`.

### 1.2 Background Index Calculation
- **Concept Name:** Background Stress Index Calculation
- **Reconstructed Equation:**
  $$B = \max\left( \sum_{i} w_i v_i + \sum_{j < k} \Gamma_{jk} v_j v_k, 0 \right)$$
  where $w_i$ are weights, $v_i$ are the values of background stressors, and $\Gamma_{jk}$ is an optional interaction matrix.
- **File Path & Line Number:** `src/background.jl:59`
- **Julia Implementation:**
  ```julia
  function background_index(stressors::Vector{BackgroundStressor}; interaction=nothing)
      B_additive = sum(w .* vals)
      B_interaction = 0.0
      if interaction !== nothing
          B_interaction = pairwise_sum(interaction, vals)
      end
      B = B_additive + B_interaction
      return max(B, 0.0)
  end
  ```
- **Notes/Observations:** By design constraint, the interaction sum only evaluates the strict upper triangular elements ($j < k$). The code also correctly computes only non-negative background stress.


## 2. Physiological Mapping (DEB-like Axes)

### 2.1 Mode of Action Vector
- **Concept Name:** Mapping Burdens to Modes of Action
- **Reconstructed Equation:**
  $$m_r = \sum_{i} U_{ri} v_i + \sum_{j < k} \Gamma_{jk}^{(r)} v_j v_k$$
  where $U$ is a weight matrix mapping $n$ stressors to $m$ modes of action, and $\Gamma^{(r)}$ are optional interaction matrices for each mode $r$. The result can be clamped.
- **File Path & Line Number:** `src/mode_of_action.jl:12`
- **Julia Implementation:**
  ```julia
  modes = mapping.U * values
  if mapping.interactions !== nothing
      for r in 1:n_modes
          for j in 1:(n_vars-1)
              for k in (j+1):n_vars
                  modes[r] += int_mat[j, k] * values[j] * values[k]
              end
          end
      end
  end
  ```
- **Notes/Observations:** Linear mapping with optional strict upper-triangular quadratic pairwise interaction terms. Results can optionally be clamped to be non-negative and $\le 1.0$.

### 2.2 DEB Axes Mapping
- **Concept Name:** Mapping Variables to DEB Axes Costs
- **Reconstructed Equation:**
  $$s_a = \sum_{i} W_{ai} v_i + \sum_{j < k} \Gamma_{jk}^{(a)} v_j v_k$$
  where $W$ maps variables to DEB axes (Assimilation, Maintenance, Growth, Reproduction).
- **File Path & Line Number:** `src/deb_axes.jl:30`
- **Julia Implementation:**
  ```julia
  s = W * values
  if !isnothing(mapping.interactions)
      for a in 1:4
          s[a] += pairwise_axis_interaction(mapping.interactions[a], values)
      end
  end
  ```
- **Notes/Observations:** Mathematical structure is identical to Mode of Action mapping. Returns a NamedTuple with axes: `assimilation`, `maintenance`, `growth`, `reproduction`.


## 3. State Variables and Margin

### 3.1 Adaptive Margin Background Penalty
- **Concept Name:** Background Stress Penalty ($\phi$)
- **Reconstructed Equation:**
  If saturating (Hill function):
  $$\phi(B) = \frac{\alpha B^{h}}{K_B^{h} + B^{h}}$$
  Otherwise (linear):
  $$\phi(B) = \alpha B$$
- **File Path & Line Number:** `src/background.jl:15`
- **Julia Implementation:**
  ```julia
  function phi_background(B::Real, params::BackgroundParams)
      if params.use_saturating_phi
          num = params.alpha * (B ^ params.hill)
          den = (params.KB ^ params.hill) + (B ^ params.hill)
          return num / den
      else
          return params.alpha * B
      end
  end
  ```

### 3.2 Adaptive Margin
- **Concept Name:** Adaptive Margin Formulation (Scalar background)
- **Reconstructed Equation:**
  $$A(B) = A_0 - \phi(B)$$
- **File Path & Line Number:** `src/background.jl:25`

### 3.3 Adaptive Margin from DEB Axes
- **Concept Name:** Adaptive Margin Formulation (DEB axes formulation)
- **Reconstructed Equation:**
  $$A(s) = A_0 - \sum_a \alpha_a s_a$$
  where $s = [s_{A}, s_{M}, s_{G}, s_{R}]^T$ and $\alpha$ are the scaling costs.
- **File Path & Line Number:** `src/deb_axes.jl:72`
- **Julia Implementation:**
  ```julia
  function deb_adaptive_margin(axes, params::DEBAxisParams)
      s = _deb_axes_to_vector(axes)
      alpha = collect(params.alpha_axes)
      return params.A0 - sum(alpha .* s)
  end
  ```

### 3.4 Buffer Derivative
- **Concept Name:** Condition Buffer (Reserve Dynamics)
- **Reconstructed Equation:**
  $$\frac{dZ}{dt} = \rho_A(1 - s_A) - \rho_M(1 + s_M) - \rho_E C_{event} - \delta_Z Z$$
- **File Path & Line Number:** `src/condition_buffer.jl:14`
- **Julia Implementation:**
  ```julia
  dZdt = params.rho_A * (1.0 - sA) - params.rho_M * (1.0 + sM) - params.rho_E * C_event - params.delta_Z * Z
  ```
- **Notes/Observations:** Provides buffer dynamics, dependent on assimilation and maintenance DEB axes $s_A$ and $s_M$, as well as an event cost $C_{event}$.

### 3.5 Adaptive Margin with Buffer
- **Concept Name:** Adaptive Margin Formulation (incorporating Z buffer)
- **Reconstructed Equation:**
  $$A(s, Z) = A_0 + \omega_Z Z - \sum_a \alpha_a s_a$$
- **File Path & Line Number:** `src/condition_buffer.jl:59`
- **Julia Implementation:**
  ```julia
  function adaptive_margin_with_buffer(axes, Z, deb_params::DEBAxisParams, buffer_params::ConditionBufferParams)
      return deb_params.A0 + buffer_params.omega_Z * Z - sum(alpha .* s)
  end
  ```


## 4. Resilience and Restoring Force

### 4.1 Restoring Force $\lambda$
- **Concept Name:** Restoring Force Calculation
- **Reconstructed Equation:**
  $$\lambda(A) = \lambda_{min} + (\lambda_{max} - \lambda_{min}) \frac{[A]_+}{K_A + [A]_+}$$
  where $[A]_+ = \max(A, 0)$.
- **File Path & Line Number:** `src/background.jl:29` and `src/deb_axes.jl:78`
- **Julia Implementation:**
  ```julia
  function restoring_force(B::Real, params::BackgroundParams)
      A_val = adaptive_margin(B, params)
      A_pos = positive_part(A_val)
      diff = params.lambda_max - params.lambda_min
      return params.lambda_min + diff * A_pos / (params.KA + A_pos)
  end
  ```
- **Notes/Observations:** A saturating Michaelis-Menten style function dependent on the positive part of the adaptive margin.

### 4.2 Restoring Force with Axes Penalty
- **Concept Name:** Restoring Force $\lambda$ with DEB Axes & Buffer Penalty
- **Reconstructed Equation:**
  $$\lambda(A, s, Z) = \min \left( \lambda_{max}, \max \left( \lambda_{min}, \lambda(A) \cdot \exp \left( - \sum_a \beta_a s_a \right) \cdot \exp(\beta_Z Z) \right) \right)$$
- **File Path & Line Number:** `src/deb_axes.jl:83`
- **Julia Implementation:**
  ```julia
  penalty = 1.0
  if params.use_axis_recovery_penalty
      penalty *= exp(-sum(beta .* s))
  end
  if params.use_buffer_recovery_factor && Z !== nothing
      penalty *= exp(params.beta_Z * Float64(Z))
  end
  lambda = base_lambda * penalty
  return clamp(lambda, params.lambda_min, params.lambda_max)
  ```


## 5. Perturbation Events and Amplification

### 5.1 Amplification Factor
- **Concept Name:** Amplification Factor Calculation
- **Reconstructed Equation:**
  $$\mathcal{F} = \frac{\lambda(A_0)}{\lambda(A)}$$
  or $\frac{\lambda(B=0)}{\lambda(B)}$.
- **File Path & Line Number:** `src/background.jl:38` and `src/deb_axes.jl:108`
- **Julia Implementation:**
  ```julia
  function amplification_factor(B::Real, params::BackgroundParams)
      l_0 = restoring_force(0.0, params)
      l_B = restoring_force(B, params)
      return l_0 / l_B
  end
  ```

### 5.2 Pulse Exposure Matrix
- **Concept Name:** Simple Rectangular Pulse Implementation
- **Reconstructed Equation:**
  $$P_{ij} = \begin{cases} A_j & \text{if } t_i \in [t_{start}^{(j)}, t_{end}^{(j)}] \\ 0 & \text{otherwise} \end{cases}$$
- **File Path & Line Number:** `src/pulses.jl:10`

### 5.3 TKTD Burden Matrix Update
- **Concept Name:** Toxicokinetic Exact Burden Update
- **Reconstructed Equation:**
  Assuming $P_n$ is constant over $[t_n, t_{n+1})$, the exact update solves $dD/dt = -k D + P$:
  $$D_{n+1} = D_n e^{-k \Delta t} + P_n (1 - e^{-k \Delta t})$$
- **File Path & Line Number:** `src/pulses.jl:35`
- **Julia Implementation:**
  ```julia
  for j in 1:n_pulses
      k_j = pulses[j].k
      exp_term = exp(-k_j * dt)
      D[i+1, j] = D[i, j] * exp_term + P[i, j] * (1.0 - exp_term)
  end
  ```

### 5.4 Acute Mixture Burden
- **Concept Name:** Aggregation of TKTD Burdens
- **Reconstructed Equation:**
  $$M = \sum_j w_j D_j + \sum_{j < k} \Gamma_{jk} D_j D_k$$
- **File Path & Line Number:** `src/pulses.jl:50`

### 5.5 Two Timescale Simulation Response
- **Concept Name:** Response Dynamics to Pulses
- **Reconstructed Equation:**
  $$ \frac{dy}{dt} = -\lambda y + q M(t) $$
  Solved with exact update assuming piecewise constant $M$:
  $$y_{n+1} = y_n e^{-\lambda \Delta t} + \frac{q M_n}{\lambda} (1 - e^{-\lambda \Delta t})$$
- **File Path & Line Number:** `src/simulation.jl:46`
- **Julia Implementation:**
  ```julia
  y[i+1] = y[i] * exp_term + (q * M[i] / lam) * (1.0 - exp_term)
  ```

### 5.6 Metrics: Area Under Curve (AUC)
- **Concept Name:** Trapezoidal Numerical Integration for AUC
- **Reconstructed Equation:**
  $$\text{AUC} = \sum_{i=1}^{n-1} \frac{y_i + y_{i+1}}{2} (t_{i+1} - t_i)$$
- **File Path & Line Number:** `src/metrics.jl:1`
- **Julia Implementation:**
  ```julia
  function trapezoid_auc(t::Vector{Float64}, y::Vector{Float64})
      n = length(t)
      auc = 0.0
      for i in 1:(n-1)
          dt = t[i+1] - t[i]
          auc += (y[i] + y[i+1]) / 2.0 * dt
      end
      return auc
  end
  ```
- **Notes/Observations:** Metric used for computing AUC of burden $M$ and response $y$.
