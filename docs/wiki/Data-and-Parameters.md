# Data & parameters

[← Equations](Equations.md) · next: [Getting started →](Getting-Started.md)

The framework draws on three data sources (capacity, pressure, memory) plus a
derived archetype database. This page describes each, the **offline mapping** that
turns raw AmP parameters into capacity numbers, and the **proxy assumptions** that
carry weight on thin evidence.

## Data files (`data/`)

| File | Role |
| --- | --- |
| `AmP_Species_Library.json` | **Capacity.** Per-species `A0`, `alpha_axes`, `lambda_bounds` (`λ_min`, `λ_max`), plus auxiliary `L_m`, `p_Am`, `p_M`, `k_M`, `E_G`, `g`. Derived offline from AmP. (Legacy records may carry a now-ignored `KA` field.) |
| `ECOTOX_Toxicity_Library.json` | **Pressure.** Parsed ECOTOX records (NOEC/EC50, effect codes, taxa). Built from raw ASCII by [`ECOTOXParser.jl`](../../src/ECOTOXParser.jl). |
| `Compound_Memory_Library.csv` | **Memory.** Per-compound retention `ρ` and bioaccumulation `K`. |
| `AmP_Species_Archetypes.csv` / `.json` | **Derived.** Response-capacity archetype labels (see [Species archetypes](../species_archetypes.md)). Built from AmP diagnostics; does not modify AmP. |
| `allStat.mat` | Raw AmP `allStat` dump, input to the offline translator. |

## The offline AmP → capacity mapping

This is the most load-bearing step and a frequent source of confusion:

- **`src/AmP_Translator.jl`** reads `data/allStat.mat`, applies the mapping in
  [Model equations §1](Equations.md#1-capacity-mapping-offline-amp--parameters),
  and writes `data/AmP_Species_Library.json`. It is a **standalone script**,
  depends on `MAT`, and is **not `include`d** in the module.
- **`src/amp_library.jl`** only *loads and validates* that JSON at runtime
  (`load_amp_species_library`, `amp_species_deb_params`, …). It contains none of
  the `{p_Am, p_M, κ, v}` math.

**Therefore:** to change how capacity is derived (e.g. to address the κ-collapse),
edit `AmP_Translator.jl` and **regenerate the JSON** — editing `amp_library.jl`
does nothing to the mapping. The α-axes and λ-bounds reaching the model are frozen
artifacts in the JSON.

```powershell
# Regenerate the capacity library after editing AmP_Translator.jl
julia +release --project=. src/AmP_Translator.jl
```

## Proxy assumptions (keep these loud)

Several inputs are **modelling proxies**, not measured quantities. They are honest
defaults, but they carry real weight and must stay visible in any results:

- **`ρ` (retention) and `K` (bioaccumulation)** are class-level defaults in
  `Compound_Memory_Library.csv`, not per-compound measured kinetics. The route to
  defensible values is one-compartment toxicokinetics (BCF/BAF, elimination
  constants).
- **ECOTOX `effect_code` as mode-of-action.** The routing of a compound to a DEB
  axis is a physiological-mode-of-action (pMoA) assignment proxy. This is a known,
  formalised problem in the DEBtox literature; the routing is defensible but
  approximate.
- **Recovery-curve shape.** The restoring force is now **linear** in the relative
  margin `A/A0` (no half-saturation constant); the old `KA = 0.3·A0` knob was removed
  — see [Limitations §2](Limitations-and-Open-Questions.md).

## Parameter reference

| Parameter | Source | Meaning |
| --- | --- | --- |
| `A0` | AmP (`E_m = p_Am/v`) | baseline adaptive margin (reserve density) |
| `alpha_axes` | AmP (`1/E_m`, `1/L_m`, `κ`, `1−κ`) | per-axis sensitivities |
| `λ_max` | AmP (`v/L_m`) | fast recovery-rate bound |
| `λ_min` | AmP (`min(k_M, λ_max)`, `k_M = [p_M]/[E_G]`) | slow recovery floor = somatic maintenance rate constant ([why](Limitations-and-Open-Questions.md)) |
| `k_M`, `E_G`, `g` | AmP (`auxiliary_metrics`) | maintenance rate constant, cost of structure, energy investment ratio |
| *(KA removed)* | — | the recovery curve is now linear in `A/A0`; no half-saturation constant |
| `ρ`, `K` | `Compound_Memory_Library.csv` | memory kinetics |
| `NOEC`, `EC50` | ECOTOX | stress anchors |
| effect code | ECOTOX | axis routing (pMoA proxy) |

For the runtime API that consumes these, see [Getting started](Getting-Started.md)
and [`PACKAGE_CAPABILITIES.md`](../PACKAGE_CAPABILITIES.md).
