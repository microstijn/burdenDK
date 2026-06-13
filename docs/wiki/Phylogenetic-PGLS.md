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
covariance is essentially uninformative, so this pass **cannot meaningfully adjudicate
phylogenetic non-independence**. The *dated*-tree test that resolves it is now **done** —
see "Dated-tree PGLS — DONE" below (λ≈0 again; the rank signal survives).

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

## Dated-tree PGLS — DONE (2026-06-13): the genuine phylogenetic test

A dated **TimeTree** (184 tips over the 197 COMADRE-matched species, real branch lengths) was
obtained and the PGLS re-run on the **182** matched species, with Pagel's λ by ML
(`scripts/comadre_pgls_dated.jl` linear, `scripts/comadre_pgls_dated_rank.jl` rank). The
**linear vs rank contrast is the result**:

| `k_M` → recovery, + generation-time control, dated VCV | β\* | p |
| --- | --- | --- |
| **linear** (log-linear) PGLS | 0.009 | 0.96 |
| **rank** PGLS (phylogenetic Spearman) | **0.221** | **0.011** |

- **Pagel's λ ≈ 0 *again* on the dated tree** (0.00–0.13): the damping-ratio trait carries little
  phylogenetic signal, so PGLS ≈ OLS and the earlier within-Order proxy was adequate — **phylogeny
  was never the confound.**
- **In rank form, `k_M`→recovery SURVIVES** the generation-time control *and* the dated-tree
  correction (β\*=0.221, p=0.011) — barely below the non-phylogenetic Spearman partial (+0.264). The
  *linear* form nulls (p=0.96), confirming the effect is **monotone, not log-linear** (`λ(A0)` drops
  to n.s.; `g` null). Results: `data/external/comadre_pgls_dated_results.txt` (linear),
  `comadre_pgls_dated_rank_results.txt` (rank); tree `comadre_amp_dated_tree.nwk`.

## Verdict (updated)
**The `k_M`↔recovery anchor holds.** It is a genuine **rank/monotone** effect that survives the
strongest control combination applied anywhere in the programme — **pace-of-life *and* a real dated
phylogeny, simultaneously** — while its log-linear form is weak. Phylogeny neither rescues nor
refutes it (λ≈0). So the honest reading is exactly the scorecard's: **rank-robust, magnitude-modest,
specification-sensitive** — corroboration, reported as a monotone tendency, not a linear effect. Repo
notes: `docs/notes/comadre_pgls_validation.md`, `external_validation_synthesis.md` §3.

## Sources
- Open Tree of Life synthetic tree / induced-subtree API.
- Grafen (1989) *Phil. Trans. R. Soc. B* 326:119–157; Pagel (1999) *Nature* 401:877–884.
