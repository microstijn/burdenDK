# Where next — strategic roadmap (2026-06-13)

*Written after the validation programme reached a manuscript-ready state and the water-quality pMoA
coupling foundation was laid. Everything below is committed/pushed; the wiki is current. This is a
planning document: the genuine forward directions, honestly costed, with a recommended sequence.*

## Where we are (one paragraph)

The **external-validation programme is effectively complete and drafted** — a standalone 6-page
manuscript (`docs/tex/external_validation_manuscript.tex`, with the forest-plot figure) built from the
consolidated synthesis. Support lands on the **margin/recovery layer** at four organisational levels;
the amplification scalar `g/F` is null throughout; the one positive rate anchor (`k_M`) survives pace +
a real dated phylogeny in rank form; two well-powered **negative controls** bound the maintenance
claim. The **water-quality coupling foundation** exists (knob-free pMoA routing table + loader +
fitted-pMoA evidence + the metals dual-route settled). The model's **distinctive content — the
across-species capacity weighting — remains untested** and is carried as an assumption.

## The three live directions

### A. Finish & submit the validation manuscript  *(near-term · low-risk · high-value · recommended first)*
- **State:** complete first draft, compiles, forest plot in. Placeholder author block.
- **What's left:** author list/affiliations; **pick a target journal** (drives length/scope — candidates:
  *Ecological Modelling*, *Environmental Toxicology & Chemistry*, *Science of the Total Environment*);
  **verify the bibliography** page/volume fields (I filled plausible values); optionally expand Methods
  (the AmP→capacity extraction; the COMADRE individual→population scale-bridge argument) and add 1–2
  figures (the `F`-vs-`g` figure already exists in `docs/wiki/figures/`; a methods schematic of the
  capacity–pressure–memory chain would help); cover-letter framing.
- **Buys:** communicates the done science and locks the validated baseline before anything is built on
  top. **Blocker:** none — your call on journal + authors.

### B. Test the capacity weighting — the distinctive, untested claim  *(medium-term · data-blocked)*
- **Why it matters:** every individual-level test is single-species, so the AmP capacity weighting (A0,
  κ-rule axis weights — the model's novel content) is held constant; and the single-trait `k_M`→toxicity
  link is body-size-confounded (synthesis §7b). So *the thing that makes the model distinctive has never
  been tested.* Until it is, it is honestly an assumption, like the mixture rules.
- **What it'd take:** across-species contaminant-gradient data with a *shared* outcome (SFG or
  survival-in-air) — which "largely does not exist" (the SFG corpus is mussel-dominated; ICES DOME
  dropped SFG; non-mussel SFG is temperature-driven). Realistic routes: (1) a **purpose-built
  experiment** (several species, a common contaminant gradient, SFG/SoS) — the clean test; (2) **mine a
  multi-species mesocosm dataset** if one with co-located burden + a common energetic/survival endpoint
  exists; (3) the cross-species GUTS `k_r` route (see D) — currently data-starved.
- **Buys:** the only validation of the model's distinctive mechanism. **Blocker:** *data generation*,
  not analysis. This is the central scientific question, and it is a data problem.

### C. Build the water-quality application  *(medium-term · design-dependent · the north-star deliverable)*
- **State:** pMoA routing foundation ready (`data/pMoA_Stressor_Routing.csv`, `src/pmoa_stressor_routing.jl`);
  the pipeline that *uses* it is **not designed** (the DynQual prototype was an explicit tryout).
- **What it'd take:** (1) **pin the target water-quality model + its exact outputs** (e.g. DynQual:
  BOD/TDS/fecal-coliform/temperature — aggregate stressors, not CAS) — this decides everything; (2)
  build a **clean coupling that reads the pMoA table** (replacing inline weights), with temperature as a
  rate-modifier not a pressure; (3) decide the **exposure-vs-field metals route** for the
  modelled-burden case (the burden is exposure-derived, so the *fitted assimilation* route may apply —
  a judgment to make explicitly); (4) **sensitivity-test the routing** (perturb axis/split/ρ → show
  site/month ranking stability) — this *is* the defensibility statement for the application.
- **Buys:** the actual instrument for "relative vulnerability statements about monthly water quality."
  **Blocker:** the input-model decision; then it's an engineering build.

## Smaller / parallel threads
- **Powered dynamic test (`k_r`):** the on-thesis rate test is data-starved (GUTS-proper `k_r` is
  *G. pulex*-centric, chemical-specific). Revisit only if a multi-species `k_r` compilation surfaces.
- **Evidence gathering (your task):** Tier-3 aggregate-stressor ecophysiology (microplastics, hypoxia,
  salinity, acidification, pathogens) for the routing rows with no DEBtox fit; EnviroTox consensus MOA
  for bulk CAS→class. List: `docs/notes/pmoa_evidence_to_gather.md`.
- **Loose ends:** the `lambda_min_maintenance_rate.pdf` dropped at session start is still **unread**
  (possibly relevant to the `λ_min=k_M` core); the GitHub **Wiki-tab mirror** is a manual re-sync if you
  want the wiki public.

## Recommended sequencing
1. **Ship the manuscript (A).** It's done science; the marginal value is in communicating it. Decide
   journal + authors, verify the bib, add the methods schematic, submit.
2. **In parallel, fix the water-quality input model** so C can start against a stable target rather than
   a moving one.
3. **Frame B as the strategic research question.** It needs purpose-built across-species data; scope a
   *minimal* multi-species SFG/SoS design now so the data can be gathered while A and C proceed.

## Decision points that are yours to make
- **Journal + authorship** for the validation paper.
- **Which water-quality model** the coupling targets (and whether its metals burden is treated as
  exposure-derived → fitted assimilation, or field-like → maintenance confound-control).
- **Whether to invest in purpose-built across-species data** to test the capacity weighting — the only
  route to validating the model's distinctive content.
