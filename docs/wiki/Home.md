# burdenDK / TwoTimescaleResilience — Wiki

`burdenDK` maps chronic environmental pressure → erosion of a physiological
**adaptive margin** `A_t` → weakened **restoring force** `λ(A_t)` → amplified burden
of a later acute event. Over 2026 the model was reframed to be **margin-first** (the
adaptive-margin *state* is the product; the amplification scalar `F` is a derived
readout), and validated for the first time against **external** data.

This wiki documents the external-validation programme against the **COMADRE** animal
matrix database. For the running design notes see `docs/notes/` and the cross-session
handoff `docs/claude/validation_roadmap_phylo_peraxis_2026-06-12.md` in the repo.

## The headline
The **DEB maintenance rate constant `k_M`** predicts independent demographic recovery,
and — newly — the **DEB reproduction rate `R_i`** specifically predicts the demographic
**compensation** component. The amplification scalar `g`/`F` predicts nothing. All
external support lands on the **margin/recovery layer**, consistent with the
margin-first reframe.

## Pages
- **[COMADRE External Validation](COMADRE-External-Validation.md)** — the scalar result (`k_M` ↔ recovery), GBIF
  species name harmonisation, and the honest specification-sensitivity caveat.
- **[Phylogenetic PGLS](Phylogenetic-PGLS.md)** (Idea A) — a real Open-Tree-of-Life phylogeny + PGLS;
  what it could and could not adjudicate.
- **[Per-Axis Resilience](Per-Axis-Resilience.md)** (Idea B) — the per-axis test: different DEB process
  rates predict different demographic-resilience components.
- **[Reproducibility](Reproducibility.md)** — exact commands and data provenance.

## One-table summary

| test | model quantity | demographic quantity | result |
| --- | --- | --- | --- |
| scalar (#1) | `k_M` | recovery (damping ratio) | +0.19\* beyond pace + Order |
| phylogeny (Idea A) | `k_M` | recovery | rank-robust, log-linear-fragile; tree too weak to test phylogeny |
| per-axis (Idea B) | `R_i` | compensation (reactivity) | **+0.77\*\*** beyond pace + size |
| per-axis (Idea B) | `g`/`F` | any | null (as in every test) |

*\* p<0.05, ** p<0.01. Last reviewed 2026-06-12.*
