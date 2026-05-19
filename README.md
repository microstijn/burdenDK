# TwoTimescaleResilience.jl

A Julia package implementing a two-timescale multi-stressor DEB-TKTD resilience framework.

## Concept

This package models how persistent background water-quality stressors reduce the "adaptive margin" and "restoring force" of an ecosystem or organism, thereby amplifying the response to additional short perturbation pulses.

The framework is divided into two layers:
1. **Slow Background Layer:** Represents long-term or persistent stressors (e.g., monthly averages) which dictate the baseline health and resilience of the system.
2. **Fast Perturbation Layer:** Represents acute perturbation pulses (e.g., short chemical spills, acute heatwaves) modeled via Toxicokinetic-Toxicodynamic (TKTD) equations.

**Note:** The framework does *not* claim that monthly global water-quality data resolve acute contaminant pulses. The pulse layer serves as a standardized perturbation probe to evaluate the background-conditioned restoring force.

## Key Equations

The central output of the framework is the amplification factor:
$$F(B) = \frac{\lambda(0)}{\lambda(B)}$$

When perturbations are simulated, the area under the curve (AUC) of the response coordinate $y$ is approximately:
$$\text{AUC}_y(B) \approx \frac{q}{\lambda(B)} \text{AUC}_{M_P}$$
which leads to the central analytical result that for identical pulse burdens under $B = 0$ and $B > 0$:
$$\frac{\text{AUC}_y(B)}{\text{AUC}_y(0)} \approx \frac{\lambda(0)}{\lambda(B)} = F(B)$$

## Installation

```julia
using Pkg
Pkg.add("TwoTimescaleResilience")
```

## Minimal Scalar Example

```julia
using TwoTimescaleResilience

params = BackgroundParams()
B = background_index([BackgroundStressor("TDS", 1.0, 0.5)])

A = adaptive_margin(B, params)
lam = restoring_force(B, params)
F = amplification_factor(B, params)

println("Background burden: $B")
println("Restoring force: $lam")
println("Amplification factor: $F")
```

## Multi-Stressor Simulation Example

```julia
using TwoTimescaleResilience

bg_stressors = [BackgroundStressor("s1", 1.0, 0.5)]
pulses = [PulseStressor("p1", 10.0, 5.0, 10.0, 1.0, 1.0)]

params = BackgroundParams()
sim = simulate_two_timescale("Demo", bg_stressors, pulses, params; tmax=50.0, dt=0.1)

metrics = compute_metrics(sim)
println("Response AUC: ", metrics.auc_y)
```

## Raster Amplification and CairoMakie Example

```julia
using TwoTimescaleResilience
using CairoMakie

# Generate synthetic grid and export outputs
params = BackgroundParams()
Bgrid, Agrid, lambdagrid, Fgrid = run_synthetic_raster_demo(params; output_dir="output")

# Or plot directly
fig = plot_amplification_grid(Fgrid)
```

## Testing

To run the package tests:
```julia
using Pkg
Pkg.test("TwoTimescaleResilience")
```
