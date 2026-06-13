# pMoA stressor-routing — evidence to gather

*Shopping list to firm up `data/pMoA_Stressor_Routing.csv`. Prioritised by yield. Drafted 2026-06-13.
Context: the in-house ECOTOX effect-code shortcut was prototyped and found **confounded** (the DEB
κ-rule cascade makes reproduction the most-sensitive endpoint for almost any pMoA; taxon mix differs
across endpoints). So **fitted DEBtox pMoAs are the gold standard**, not raw endpoint sensitivities.*

## Tier 1 — fitted-pMoA compilations (the gold standard; replaces a-priori guesses)

These report which energetic process a toxicant *taxes*, from DEBtox model fits (not endpoint
sensitivity). Use these to set the chemical rows and upgrade `confidence`.

1. **"Physiological modes of action across species and toxicants"** — *Environ. Sci.: Processes &
   Impacts* 2018 (doi 10.1039/c7em00328e). The single most on-target source: a cross-species×toxicant
   compilation of fitted pMoAs. **Get this first.**
2. **Jager, "Revisiting simplified DEBtox models for analysing ecotoxicity data"** (2019/2020). *You
   already have the PDF in `~/Downloads`.* The canonical pMoA definitions — cite for the framework.
3. **DEBtox / DEB-TKTD per-compound fits** (each gives a fitted pMoA for one or more compounds):
   - Daphnia metals Cu/Ni/Zn, *combined* pMoAs — *ETC* 43(2):338, 2024 (doi 10.1002/etc.5793).
   - Pb in *Lymnaea stagnalis* — *Environ. Pollut.* 2023 (PMID 37369157).
   - C. elegans Cd / fluoranthene / atrazine, pMoA↔gene expression — PMC2857823 (also: pMoA class is
     method-dependent — a caution).
   - openGUTS / Jager & Ashauer GUTS e-book worked examples (also feed the dynamic `k_r` side).

## Tier 2 — chemical MoA-class databases (broad routing for many CAS at once)

For mapping many compounds to a coarse MoA class (then to an axis), rather than one-by-one:
4. **MOAtox** (Barron et al. 2015, *ETC* 34:2370) — curated aquatic mode-of-action database (~1000s
   of chemicals: narcosis, AChE inhibition, reactivity, etc.). Maps CAS → broad MoA; bridge MoA→axis.
5. **Verhaar scheme** (Verhaar et al. 1992) + its automated implementations (e.g. Toxtree) — narcosis
   vs reactive vs specific-acting classes; the standard a-priori MoA classifier.
6. **US-EPA ECOTOX** (*already on disk*). Use it NOT for pMoA (confounded — see header) but for
   axis-specific endpoint **sensitivity rankings within a species**, and to flag
   reproduction-disproportionate (endocrine) stressors via the growth-vs-reproduction ratio.

## Tier 3 — aggregate stressors with no DEBtox fit (the genuinely new rows)

These are not single chemicals, so no DEBtox pMoA exists — gather mechanistic ecophysiology instead.

7. **Microplastics → assimilation.** Strong support that MP cut clearance/feeding (→ assimilation),
   with growth/reproduction as the DEB downstream:
   - Manila clam (*Ruditapes*) — *Environ. Pollut.* 2021 (S0269749121020844).
   - *Mytilus* MP accumulation + DEB — *Ocean Science* 16:927, 2020.
   - Mediterranean cultured mussels — *Environ. Pollut.* 2024 (S0269749124017676).
   (Reported: clearance ↓8–25% → growth ↓6–16% → reproduction ↓7–19%.)
8. **Organic pollution / hypoxia (BOD) → maintenance.** Aerobic-scope / metabolic-cost-of-hypoxia
   literature (e.g. metabolic index / oxygen- and capacity-limited thermal tolerance; Pörtner;
   Deutsch et al. metabolic index). Gather one quantitative aerobic-scope-vs-DO source.
9. **Salinity / osmotic → maintenance.** Osmoregulation-energetics: cost of ion pumping as a
   maintenance term (euryhaline fish/invertebrate energetics reviews).
10. **Acidification / low pH → maintenance + growth.** Ocean-acidification energetics / calcification
    cost in calcifiers (DEB-OA studies in bivalves/echinoderms).
11. **Pathogens / immune → maintenance (+ reproduction trade-off).** Eco-immunology cost-of-immunity
    literature; least-settled row, keep `confidence = low`.
12. **Nutrients/eutrophication** — flagged indirect; prefer to model via its BOD/hypoxia consequence
    (row 8) rather than gather direct evidence.

## Tier 4 — trait-based parallel (alignment + local access)

13. **Rubach / Van den Brink (WUR) "species traits × stressor mode of action"** — *IEAM* 7:172 (2011)
    and the chlorpyrifos-traits paper (PMC3431471). Same conceptual programme at the species-trait
    level; Wageningen line → likely local access and people to corroborate the routing against.

## Databases to pull (one-off)
- MOAtox (download — supplementary of Barron 2015).
- ESPI c7em00328e (the pMoA compilation — supplementary table is the prize).
- add-my-pet / AmP (have it) — for the DEB capacity side.
- openGUTS / EFSA 2018 TKTD datasets — for the dynamic `k_r` side (separate thread).

## How to use what you gather
- **Chemical rows** → set `primary_axis` from the *fitted* pMoA (Tier 1) or MoA class (Tier 2); mark
  `confidence = high` where a fit exists.
- **Aggregate rows** → keep the mechanistic assignment (Tier 3), `confidence = medium/low`, cite the
  ecophysiology.
- **Everything** → record the source in a `references` column (add it to the CSV), then
  **sensitivity-test** the assignment and report ranking stability as the defensibility statement.
