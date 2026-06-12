# Phylogenetic PGLS (Idea A)

Related species are not statistically independent. The [COMADRE External Validation](COMADRE-External-Validation.md)
used a taxonomic-Order proxy for phylogeny; this replaces it with a **real phylogeny**
(Open Tree of Life) and **phylogenetic generalized least squares (PGLS)** with Pagel's
λ estimated by maximum likelihood.

## Method
- **Tree.** OTL TNRS matched all 197 harmonised species; the synthetic induced subtree
  placed 197/197 (topology only — OTL has no dated branch lengths). 188 re-joined to the
  model table (9 remapped by OTL to a different accepted taxon).
- **Branch lengths.** Grafen (1989) ultrametric proxy (undated topology).
- **PGLS.** Pure-Julia GLS (no phylo-package dependency): phylogenetic VCV from the
  Grafen tree, Pagel's λ by ML, `comadre_log_damping ~ log(predictor) (+ log gen. time)`,
  standardized predictors.

## Result

| | β\* alone | p | β\* \| gen. time | p | Pagel λ |
| --- | --- | --- | --- | --- | --- |
| `λ(A0)` | +0.246 | 0.044 \* | −0.087 | 0.58 | ≈0.01–0.05 |
| **`k_M`** | **+0.297** | **0.013 \*** | −0.004 | 0.98 | ≈0.00–0.04 |
| `g` | −0.087 | 0.49 | −0.156 | 0.21 | ≈0.02–0.03 |

## Two distinct findings
**1. The undated OTL+Grafen tree carries ~no phylogenetic signal — this PGLS ≈ OLS.**
ML Pagel's λ ≈ 0 for every model; the logL profile peaks at λ≈0.1 (−349.9) and falls
monotonically to full Brownian motion at λ=1 (−467.5). The Grafen-on-undated-topology
covariance is essentially uninformative, so the analysis **cannot meaningfully
adjudicate phylogenetic non-independence**. A *dated* tree (VertLife/TimeTree) is the
genuine phylogenetic test — the necessary follow-up.

**2. `k_M` predicts recovery alone, but not under a log-linear gen.-time control — and
this is a rank-vs-linear effect, not phylogeny.** `k_M` is significant on its own
(β\*=0.30, p=0.013) but vanishes once log generation time is added (β\*≈0). On the same
193 species:

| `k_M ~ recovery \| gen` | value |
| --- | --- |
| rank (Spearman) partial | **+0.264** (the headline) |
| log-linear (Pearson) partial | +0.04 |

(log `k_M` vs log gen r = −0.55 — moderate, *not* damaging collinearity.) The signal is
**monotone but not log-linear** — it lives in the ranks and a linear specification does
not see it.

## Verdict
**Neither a clean confirmation nor a clean refutation — and that is informative.**
(a) The cheap all-taxa OTL+Grafen tree is too weak to test phylogeny, so a dated tree is
now clearly the necessary next step. (b) The `k_M`-beyond-pace signal is
**specification-sensitive** (robust under rank correlation, absent under log-linear
regression); both hold under their own functional-form assumption, but this bounds how
strongly a single number can be leaned on. Repo notes:
`docs/notes/comadre_pgls_validation.md`.

## Sources
- Open Tree of Life synthetic tree / induced-subtree API.
- Grafen (1989) *Phil. Trans. R. Soc. B* 326:119–157; Pagel (1999) *Nature* 401:877–884.
