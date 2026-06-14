# Paper 1 — where we are still open to critique (red-team, 2026-06-14)

*A deliberate, adversarial inventory after consolidating the across-axis weighting result. Goal: name
the sharpest objections a hostile-but-fair reviewer (DEB theorist, quantitative ecotoxicologist, or
statistician) would raise, the current defence, and the residual exposure. Ranked by how much each
threatens the paper's headline. Companion to `external_validation_synthesis.md` (the evidence) — this
file is the *weaknesses*.*

The paper's headline: a **zero-calibration, scale-free** description of capacity erosion & recovery,
validated as **monotone corroboration** at multiple organisational levels, honestly disclaiming absolute
capacity, the across-axis weighting, and a 1-D amplification scalar.

---

## TIER 1 — threatens the headline

### C1. "Zero free parameters" overclaims — structural choices are de-facto parameters, and at least one (metals routing) was selected with reference to anchor performance.
**Bite.** The honesty headline is "nothing fitted to the outcomes." But the model embeds fixed *choices*
that function as parameters: the κ-rule axis weights (why assimilation = 1/2?), the MoA→axis **routing
table**, `A0 = E_m`, the **linear** `λ(A)`, `λ_min=k_M`, `λ_max=v/L_m`. Worse, the routing of **metals**
was set to *maintenance* after a re-run showed metals→assimilation **degrades all five field anchors**
(`pmoa_evidence_digest.md` §4-RESULT). That is a structural decision **informed by the validation
outcomes** — exactly what "zero calibration" is meant to exclude.
**Defence.** The metals→maintenance choice has an *independent* principled basis (field tissue-burden ≠
exposure; metals track food/condition — a positive confound), not merely "it fit better"; signs were
pre-registered; the curves are DEB-derived, not tuned.
**Residual exposure.** A reviewer who reads the digest will call the metals routing a hidden
researcher-degree-of-freedom and the "zero parameters" claim misleading. **Mitigation:** (i) soften to
"no per-outcome calibration" + an explicit **table of structural choices**; (ii) report the metals
routing as a transparent **two-way sensitivity** (both routings, both anchor results), not just the
chosen one; (iii) state the routing was fixed *a priori on mechanism* and the §4 re-run was a *check*,
if that is in fact the order — and if it is not, say so.

### C2. The model's nominal mechanism (amplification `g`/`F`) is null in **every** test, reframed as a "confirmed margin-first prediction" — an unfalsifiability risk; and the scale-free boundary may be post-hoc.
**Bite.** The framework is *named and motivated* by amplification (background narrows margin → amplifies
a later perturbation). Yet the amplification scalar predicts nothing anywhere (ρ ≈ −0.05…−0.13).
Re-describing "our central scalar fails" as "the margin-first prediction is confirmed" reads, to a sharp
reviewer, as moving the goalposts. Relatedly: the **scale-free boundary** (validated = rate/state/dynamics;
disclaimed = absolute capacity, weighting) — was it drawn *before* or *after* seeing which tests passed?
The clam/mussel "honest failure predicted by the scale-free reading" is only compelling if the prediction
preceded the result.
**Defence.** The reframe is mathematically grounded: the within-species rank tests are invariant to `A0`
and to the scalar, so the *margin* is the only identifiable object — the scalar's null is a real,
arguably-predictable consequence, not a save. Nulls are reported, not buried.
**Residual exposure.** "Keeps the wins, explains the losses." **Mitigation:** state explicitly which
predictions were **pre-registered** vs reframed; give an explicit **falsification criterion** (what
result would have refuted the validated claims?); present the amplification-null as a *finding with an
interpretation*, not a predicted success.

### C3. With the across-axis weighting now **tested and not corroborated**, what distinctive, validated content remains — and does the elaborate DEB-axis machinery earn its keep over a simple capacity proxy?
**Bite.** The model's *distinctive* content was the per-species κ-driven axis weighting. This session
tested it (powered, n=27/101 ECOTOX) → **not corroborated** (and the clean apical-EC50 version is
data-starved). The amplification scalar is null. So the validated residue is: `k_M`→recovery (modest,
rank-only) and margin↔SFG/SoS (mussel, where burden indexes exposure). A reviewer asks: *could a simple
body-size / single-DEB-covariate model reproduce these without the axis machinery?* If yes, the
distinctive apparatus is not validated content.
**Defence.** The **within-anchor ablation** shows the *routed, structured* margin beats a naive
equal-weight load at every field anchor — the *operative structure* (routing + saturating impairment +
aggregation) earns its keep, even though the cross-species *weighting* does not.
**Residual exposure.** The ablation defends *structure*, not the *distinctive cross-species weighting*,
and only on mussels. **Mitigation:** foreground the weighting-negative; delimit cleanly **validated
mechanism/structure** vs **not-corroborated weighting** vs **absolute capacity**; add a head-to-head vs a
body-size-only / single-best-DEB-covariate baseline at each anchor (not just naive load).

---

## TIER 2 — substantial, addressable

### C4. The validated effects are modest, rank-only, and specification-sensitive; the strongest control barely survives, on a trait with ~no phylogenetic signal.
ρ ≈ 0.2–0.45; the flagship `k_M`→recovery is **β\*≈0.22 in ranks** and **nulls log-linearly**; Pagel's
λ≈0 means the PGLS ≈ OLS (so "survives phylogeny" corrects for little). A quantitative reviewer: the model
explains little variance and makes no calibrated prediction. **Defence:** uncalibrated monotone
corroboration is the stated, conservative claim; rank stats + nested controls + bootstrap CIs.
**Mitigation:** report variance-explained honestly; frame as monotone tendency; note the deliberate
no-fulcrum stance (the Paper-2 calibration study).

### C5. The strongest positives may be the most tautological.
SFG ≈ "absorbed − respired energy" ≈ "capacity beyond maintenance" ≈ the margin's own definition; `R_i`
(DEB reproduction rate) → demographic **compensation** is close to "reproduction drives population
growth." **Defence:** the modelled margin is computed from *contaminant burden via routing + AmP
capacity*, using **no** SFG data — a genuine out-of-sample prediction, not a restatement; `R_i` is the
AmP-derived rate, independent of the matrices. **Mitigation:** make the data-independence explicit; test
whether a trivial "reproduction predicts compensation" baseline does as well.

### C6. Researcher degrees of freedom across many anchors and choices.
Several headline numbers improve along a chain of *defensible* choices (clam/mussel −0.82→−0.91 control-
normalised; SoS 0.39→0.45→0.62 with successive controls/QC). Signs are pre-registered, but the full
pipeline is not. **Defence:** each step is mechanistically motivated and the failures are reported.
**Mitigation:** a pre-registered analysis protocol; a "forking-paths" sensitivity (how far do results
move across defensible choices?) — the response-curve robustness does this for one axis; extend it.

### C7. DEB-derivation legitimacy (domain-expert scrutiny).
A Kooijman-school reviewer will contest the identifications: "adaptive margin = reserve density `E_m`",
the two-timescale `λ_min=k_M`, `λ_max=v/L_m`, the κ-rule weighting derivation, the linear `λ(A)`.
**Mitigation:** a formal DEB appendix deriving/justifying each mapping and pre-empting the obvious
objections; be explicit about where this *reinterprets* vs *applies* standard DEB.

---

## TIER 3 — manageable, name them

### C8. Taxonomic narrowness of the energetics anchors.
The margin's *own-level* validation (SFG/SoS) is almost entirely **Mytilus** (+ one clam). The cross-taxon
tests (fish, n=310 ECOTOX, the weighting) largely nulled. So "the margin works" is concentrated in one
genus's field physiology. COMADRE/GlobTherm are taxonomically broad but at other levels. **Mitigation:**
scope generality claims to each anchor's level; the DEBtox-refit panel (Paper-2) broadens taxa.

### C9. The dynamics evidence is thin and partly expected.
n=4 transplant + 2-species phenanthrene (half figure-digitised). The static-vs-dynamic contrast is the
novel part, but "more dose/time → more erosion → lower LT50" is close to expected, and the unfitted-
timescale argument (`1/λ ≈ months`) spans a wide band (68–887 d) that "months" fits easily.
**Mitigation:** powered, reported (not digitised) multi-species LT50 series; tighten the timescale
prediction to a sharper interval.

---

## The one-line version
The honest **negatives** (amplification-null, weighting-not-corroborated, single-trait `k_M` size-
confounded) are the paper's integrity — but they also hollow out the *distinctive* validated content,
leaving a **modest, rank-only, mussel-concentrated, uncalibrated** core whose "zero-parameter" framing has
one real crack (metals routing). The defensible thesis to lead with is narrow and structural: *a
mechanistically-structured margin beats an unstructured load index at predicting energetic resilience,
on its own organisational level, with no per-outcome calibration* — and everything beyond that is
honestly bounded. Future test that could add distinctive validated content: the **DEBtox-refit species ×
MoA panel** (`curated_apical_ecx_plan.md`).
