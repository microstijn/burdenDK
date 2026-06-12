# Reproducibility

All external-validation artifacts are committed; raw downloads are gitignored and
re-fetched by the scripts. **Use Julia 1.12.6 via `julia +release`** (the default LTS
cannot load the project — see `CLAUDE.md`).

## Environment for the standalone extractors
The extractors/resolvers depend on packages that are **not** in the project Manifest
(by design — they are offline pre-processing). Run each in a throwaway environment:

```powershell
# COMADRE extraction (RData + DataFrames)
julia +release -e 'using Pkg; Pkg.activate(mktempdir()); Pkg.add(["RData","DataFrames"])'
# GBIF name resolution + OTL tree fetch (HTTP + JSON)
julia +release -e 'using Pkg; Pkg.activate(mktempdir()); Pkg.add(["HTTP","JSON"])'
# AmP reproduction-rate side-extraction (MAT)
julia +release -e 'using Pkg; Pkg.activate(mktempdir()); Pkg.add("MAT")'
# PGLS (Distributions, for the t p-value)
julia +release -e 'using Pkg; Pkg.activate(mktempdir()); Pkg.add("Distributions")'
```

## Pipeline

```powershell
# 1. COMADRE extraction -> data/external/comadre_recovery.csv
#    (damping=recovery, generation time, + reactivity=compensation, attenuation=resistance)
julia +release --project=<rdata-env> scripts/extract_comadre_recovery.jl

# 2. GBIF name harmonisation -> data/external/comadre_amp_namemap.csv
julia +release --project=<http-env>  scripts/resolve_comadre_amp_names.jl

# 3. AmP reproduction rates -> data/external/amp_reproduction_rates.csv
#    (k_M, R_i, r_B, kap_R, k_J, Ww_i, a_p; does NOT regenerate the library JSON)
julia +release --project=<mat-env>   scripts/extract_amp_reproduction_rates.jl

# --- analyses (project env) ---
# scalar validation (raw / pace / Order-controlled)
julia +release --project=. examples/comadre_partial_validation.jl
# per-axis matrix + robustness (Idea B)
julia +release --project=. examples/comadre_peraxis_validation.jl

# --- phylogeny (Idea A) ---
julia +release --project=. scripts/export_comadre_matched_table.jl   # matched table
julia +release --project=<http-env> scripts/fetch_comadre_tree.jl    # OTL tree
julia +release --project=<dist-env> scripts/comadre_pgls.jl          # PGLS
```

## Data provenance
- **COMADRE** `COMADRE_v.4.26.4.0.RData` (~907 KB, CC-BY; Salguero-Gómez et al. 2016,
  *J. Anim. Ecol.*) — gitignored, re-fetched by the extractor.
- **AmP** `data/allStat.mat` (Add-my-Pet collection) — the source of all DEB parameters.
- **GBIF Backbone** and **Open Tree of Life** — queried live via their public APIs.

## Committed derived tables (`data/external/`)
| file | contents |
| --- | --- |
| `comadre_recovery.csv` | per-species recovery, gen. time, reactivity, attenuation, taxonomy |
| `comadre_amp_namemap.csv` | COMADRE→AmP harmonised names (+ method) |
| `comadre_amp_matched.csv` | matched model+COMADRE table (PGLS input) |
| `comadre_amp_tree.nwk` / `_tree_map.csv` | OTL induced subtree + tip→species map |
| `amp_reproduction_rates.csv` | per-species DEB process rates (k_M, R_i, r_B, a_p, …) |

See also the cross-session handoff
`docs/claude/validation_roadmap_phylo_peraxis_2026-06-12.md`.
