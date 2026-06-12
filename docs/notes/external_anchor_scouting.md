# External anchor — can we even find validation data? (#3 scouting)

*Working note, 2026-06-12. Before building any validation harness, this records the
"think hard + look around" pass: what kind of data would actually validate the
framework, whether it exists, and the trap that makes most of it unusable. No code
was built. Companion to `g_lifehistory_check.md` and `feature_redundancy_check.md`.*

## The framing that decides everything

The framework computes **recovery capacity** — chronic pressure erodes the adaptive
margin, which weakens the restoring force, which amplifies a later acute event. It
is derived **entirely from AmP physiological (capacity) parameters** and uses *no*
species toxicity data on the capacity side. So the right external anchor is an
**independent, species-resolved measure of recovery / resilience**, not of
toxicodynamic *susceptibility*.

This split is not ours to invent — **SPEAR** (Species At Risk) already formalises it:
it flags taxa that are *both* "physiologically sensitive to toxicants" *and* "cannot
recover quickly following a pulse disturbance (not resilient)" — two separate trait
sets. **Our model is the recovery/resilience half only.** Consequence:

- The **abundant** data — ECOTOX EC50/NOEC, species sensitivity distributions (SSDs),
  fitted DEBtox tolerance — measures *susceptibility*, the **wrong half**. A
  correlation (or its absence) against it would not validate recovery capacity.
- The data we **need** — independent species-resolved *recovery* — is the scarce kind.

## What exists, ranked (findability × conceptual match)

1. **COMADRE / COMPADRE demographic resilience — strongest real lead.** A comparative
   framework over ~910 animal/plant populations decomposes resilience into
   *resistance*, *compensation*, and **recovery time**, from publicly downloadable
   matrix population models. Genuine, quantitative, cross-species *recovery* data,
   independent of toxicity, overlapping AmP taxa; DEB bridges individual → population.
   Caveat: *population* recovery vs our *individual-energetic* recovery (a defensible
   scale bridge, not 1:1).
2. **freshwaterecology.info / Tachet traits — taxon-resolved resilience traits.**
   20,000+ European freshwater organisms with "resistance/resilience potential" and
   life-cycle traits; downloadable (registration); overlaps AmP aquatic invertebrates.
   Coarser (categorical), but directly recovery-flavoured.
3. **SPEAR — field-validated, but a community index.** Its species-level recovery
   traits (generation time, recolonisation) are the right comparison and it is
   validated against pesticide effects worldwide; the underlying species trait table
   is the usable artifact.
4. **Multi-stressor amplification meta-analyses — validate the *premise*, not the
   ranking.** ~142 three-stressor interactions across 38 papers; synergy/amplification
   is real but quantified at population/community level, not species-resolved
   physiology. Good for "does chronic-then-acute amplify at all?".
5. **Direct chronic-then-acute experiments — the ideal, but scattered.**
   Cross-sensitization / carry-over is a *real, documented* phenomenon (HPA priming in
   vertebrates; amphipod predation-risk carry-over) but exists as single-species
   mechanistic studies, not a ready cross-species ranking. Assembling one is a
   literature-synthesis project.
6. **ECOTOX / SSD / DEBtox tolerance — the wrong half** (susceptibility). Useful only
   as a *contrast* ("does capacity-vulnerability correlate with susceptibility?"), not
   as validation.

## The circularity trap (the reason this is hard)

The vulnerability index is built **entirely from AmP life-history parameters**. Every
recovery-flavoured anchor above (COMADRE recovery time, Tachet resilience, SPEAR
generation-time) is *also* life-history-driven. So a raw rank correlation risks the
**same partial-circularity** already seen in `g_lifehistory_check.md` (where
`g ≈ 1/{p_Am}` by construction). The only genuinely clean tests are:

- **direct chronic-then-acute amplification data** (scarce, item 5), or
- a **partial correlation controlling for generation time / pace-of-life** against
  COMADRE recovery — i.e. *does the margin model predict recovery beyond what raw
  life history already does?*

## Verdict

**Anchorable data does exist** — COMADRE demographic recovery is the best public,
quantitative, independent option, with Tachet/freshwaterecology and SPEAR as
secondary trait-based checks. But two honest caveats stand:

1. **Scale mismatch** — population (COMADRE) vs individual-energetic (our model)
   recovery; defensible via DEB but not 1:1.
2. **Circularity** — the real test is a generation-time-controlled *partial*
   correlation, not a raw rank correlation. The *ideal* direct chronic-then-acute
   data is real-but-scattered and would need a literature synthesis.

## The precise question for Kooijman / the PI

> *Is COMADRE demographic recovery time a legitimate external anchor for an
> AmP-derived recovery-capacity index — controlling for generation time — or do we
> need direct chronic-then-acute amplification data, which would be a
> literature-synthesis effort?*

Until that is answered we should **not** build a validation harness or tune anything
to a possibly-circular anchor.

## Sources

- SPEAR — [UFZ effect monitoring](https://www.ufz.de/index.php?en=38122) · [SPEARpesticides (Wikipedia)](https://en.wikipedia.org/wiki/SPEARpesticides)
- Demographic resilience / COMADRE — [Towards a Comparative Framework of Demographic Resilience (TREE)](https://www.sciencedirect.com/science/article/pii/S0169534720301312) · [Life history mediates demographic resilience trade-offs (PMC)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9314072/)
- Freshwater traits — [freshwaterecology.info](https://www.freshwaterecology.info/about.php) · [EPA Freshwater Biological Traits DB](https://www.epa.gov/risk/freshwater-biological-traits-database-traits)
- Multi-stressor meta-analysis — [Three-stressor higher-order interactions (PubMed)](https://pubmed.ncbi.nlm.nih.gov/36572303/)
- Carry-over / cross-sensitization — [Adaptive response to chronic predation risk in amphipods (Oecologia)](https://link.springer.com/article/10.1007/s00442-020-04594-z)
- AmP × TKTD — [Standard DEB animal model for TKTD (Jager)](https://www.sciencedirect.com/science/article/pii/S030438002200285X) · [The AmP project (PLOS Comp Biol)](https://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1006100)
