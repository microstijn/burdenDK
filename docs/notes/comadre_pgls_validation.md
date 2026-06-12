# External validation #1b — PGLS (real phylogeny) on COMADRE recovery

*Working note, 2026-06-12. Idea A of the validation roadmap. Replaces the
taxonomic-Order proxy (`comadre_partial_validation.jl`) with a **real phylogeny**:
an Open Tree of Life induced subtree + phylogenetic generalized least squares
(PGLS) with Pagel's λ estimated by ML. Asks: does `k_M` predict COMADRE
demographic recovery beyond pace-of-life **and** phylogenetic non-independence?
Reproduce: `scripts/export_comadre_matched_table.jl` →
`scripts/fetch_comadre_tree.jl` → `scripts/comadre_pgls.jl`.*

## Method

- **Tree.** OTL TNRS-matched all 197 harmonised AmP↔COMADRE species; the synthetic
  induced subtree placed 197/197 (topology only — OTL has no dated branch lengths).
  188/197 re-joined to the model table (9 were remapped by OTL to a different
  accepted taxon and dropped).
- **Branch lengths.** Grafen (1989) ultrametric proxy (node height = descendant-tip
  count − 1, power ρ=1), since the topology is undated.
- **PGLS.** Pure-Julia GLS (no phylo-package dependency): phylogenetic VCV from the
  Grafen tree, Pagel's λ estimated by ML (golden-section), `comadre_log_damping ~
  log(predictor) (+ log generation_time)`, standardized predictors. Subsetting the
  full-tree VCV to the generation-time subset is exact under BM (marginalization).

## Result

| | β\* alone | p | β\* \| gen. time | p | Pagel λ |
| --- | --- | --- | --- | --- | --- |
| `λ(A0)` | +0.246 | 0.044 \* | −0.087 | 0.58 | ≈0.01–0.05 |
| **`k_M`** | **+0.297** | **0.013 \*** | −0.004 | 0.98 | ≈0.00–0.04 |
| `g` | −0.087 | 0.49 | −0.156 | 0.21 | ≈0.02–0.03 |

(β\* = standardized partial slope per SD of log predictor; n=188, 185 with gen. time.)

## What it means — two distinct findings

**1. The undated OTL+Grafen tree carries ~no phylogenetic signal; this PGLS ≈ OLS.**
ML Pagel's λ ≈ 0 for every model; the logL profile (k_M+gen model) peaks at λ≈0.1
(−349.9) and falls monotonically to full Brownian motion at λ=1 (−467.5). So the
Grafen-on-undated-topology covariance is essentially uninformative — the analysis
**cannot meaningfully adjudicate phylogenetic non-independence**. A *dated* tree
(VertLife/TimeTree; needs a manual download — see roadmap Part 1) is the genuine
phylogenetic test; this all-taxa pass is only the coarse robustness check.

**2. `k_M` predicts recovery alone, but not under a log-linear generation-time
control — and this is a rank-vs-linear effect, not phylogeny.** `k_M` is a
significant positive predictor on its own (β\*=0.30, p=0.013), but vanishes once
log generation time is included (β\*≈0). Diagnostics on the same 193 species:

| partial `k_M ~ recovery \| gen` | value |
| --- | --- |
| rank (Spearman) partial | **+0.264** (the headline result) |
| log-linear (Pearson) partial | +0.04 |

(log k_M vs log gen r = −0.55 — moderate, *not* damaging collinearity; gen itself
gives −0.215 rank vs −0.19 log-linear, i.e. stable.) So the "k_M beyond
pace-of-life" signal is **monotone but not log-linear** — it lives in the ranks and
a linear specification does not see it. The Order-proxy partial Spearman (0.190\*)
and this log-linear PGLS (≈0) are *both right* about different functional forms.

## Verdict

**Neither a clean confirmation nor a clean refutation — and that is informative.**
(a) The cheap all-taxa OTL+Grafen tree is too weak to test phylogeny (Pagel λ≈0),
so the dated-tree follow-up is now clearly the necessary next step, not optional.
(b) The `k_M`-beyond-pace-of-life signal is **specification-sensitive**: robust under
rank-based partial correlation (0.26\*\*) but absent under log-linear regression.
The honest headline must state this. It does not overturn the rank result — both
hold under their own functional-form assumption — but it bounds how strongly the
manuscript can lean on a single number. Next: (i) a dated vertebrate tree for a real
phylogenetic test; (ii) characterise the rank-vs-linear gap (influential species?
nonlinearity?); (iii) the per-axis metric (Idea B), which may be less
specification-fragile than the scalar.

## Sources

- Open Tree of Life synthetic tree / induced subtree API — [opentreeoflife.org](https://opentreeoflife.org), [api.opentreeoflife.org/v3](https://github.com/OpenTreeOfLife/germinator/wiki/Synthetic-tree-API-v3)
- Grafen, A. (1989) The phylogenetic regression. *Phil. Trans. R. Soc. B* 326:119–157.
- Pagel, M. (1999) Inferring the historical patterns of biological evolution. *Nature* 401:877–884.
