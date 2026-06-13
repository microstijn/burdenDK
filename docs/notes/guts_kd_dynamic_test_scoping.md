# Scoping: the circularity-safe *powered* dynamic test (AmP `k_M` ↔ GUTS `k_D`)

*Scoping the open decision in `validation_handover_margin_sos_dynamics_2026-06-13.md` Part 5.
**Started as a feasibility memo; this session also RAN two of the tests** (rate axis §5d, state axis
§5f). §§1–4 are the original scoping; §§5a–5f are the empirical session log; the revised bottom line
is the box immediately below. Read that first.*

> ## REVISED BOTTOM LINE (after running it, 2026-06-13)
> **The single-trait maintenance→toxicity tests are body-size-confounded on *both* axes — that is the
> session's finding.**
> - **State axis** (`k_M`→sensitivity, ECOTOX LC50, **n=310**, §5f): the Baas & Kooijman signal
>   *replicates raw* (maintenance↑→more sensitive, r≈−0.27, all 4 chemicals) but **nulls under a body-
>   size control** (partial ≈ −0.03). Well-powered → a robust null, not underpower.
> - **Rate axis** (`k_M`→Rubach `k_out`, **n=6–7**, §5d): weak/wrong-sign, n.s.; and `k_out` is anyway
>   the *toxicokinetic* (elimination) rate, not the thesis-relevant *toxicodynamic* recovery rate.
> - **Clean TD rate `k_r`** is the right endpoint but is **scarce, *G. pulex*-centric, and chemical-
>   specific** (§5e) → the *powered* dynamic test stays data-starved (confirms the handover's worry).
>
> **Consequence:** the distinctive, defensible content is **not** "`k_M` predicts toxic response"
> (that's a size story) but the **across-axis capacity weighting + MoA routing** — which still needs
> purpose-built across-species gradient data (handover Part 6.1). Neither axis is worth writing up as a
> positive result; both are valuable **honest negative controls** that bound the maintenance claim.

## 1. The question and the *original* verdict (pre-data)

**The test (Part 5).** Across species, does **AmP `k_M`** (clean life-history maintenance rate)
predict the **GUTS dominant rate constant `k_D`** (fitted to *raw* survival time-courses by a
non-DEB model), **controlling body size + phylogeny**? It is the toxicodynamic twin of the COMADRE
`k_M`↔recovery test: non-circular (`k_D` never touches a DEB quantity), powered (more `k_D` values
than we ever had SoS sites), on-thesis (tests the dynamic prediction that the erosion timescale
≈ the maintenance timescale `1/k_M`).

**Original verdict (pre-data): GO — conditional.** The overlap is real and the dynamic range is
excellent, but `k_D` must be **hand-assembled from ~6–10 source papers** (no single large public
`k_D` database), and the result is **at the mercy of a strong body-size collinearity** (see §4).
*This session then ran it and the body-size collinearity won — see the revised bottom line above.*

## 2. The two sides of the overlap

**AmP side — abundant.** `data/AmP_Species_Library.json` holds **7,335 species**, each with `k_M`
and `L_m` (body size) directly. Coverage of the standard GUTS test-species universe is **excellent**:
of 30 canonical test species checked, **23 are in AmP**; the 7 misses are oligochaetes/rotifers or
genus-only entries (Lumbriculus, Brachionus, Enchytraeus, Sialis, Notonecta, Chironomus dilutus,
Neocaridina). Vertebrate/amphibian test species are also present (Poecilia, Clarias, Lithobates,
Xenopus, Salmo, Oncorhynchus kisutch, Carassius, Ictalurus, Anguilla, …). **AmP is not the
binding constraint.**

**GUTS side — the binding constraint, and it is scattered, not centralised.** There is no
"hundreds of species" `k_D` database. `k_D` lives per-paper:
- **Ashauer et al. 2016** (*Sci Rep* 6:29178, open access / PMC4933929) — **6 species**: *Gammarus
  pulex* (4 pesticides) + a malathion SSD over *Poecilia reticulata, Rana sylvatica, Clarias
  gariepinus, Rana catesbeiana, Pimephales promelas*. Also explicitly relates `k_D` to body size —
  directly relevant to §4.
- **Nickisch Born Gericke et al. 2022** (*ETC* 41:1732, PMC9328144) — **5 fish** (rainbow trout,
  fathead minnow, common carp, + 2) calibrated GUTS-RED-IT for benzovindiflupyr.
- **Focks et al. 2018** (*Ecotoxicology* 27:992, PMC6132984) — three neonicotinoids × several
  aquatic **macroinvertebrates**.
- **openGUTS ring-test datasets** (openguts.info) + **EFSA 2018 TKTD opinion** (*EFSA J* 16:5377)
  case studies — a handful of standard species (Daphnia magna, Gammarus, fish).
- **Jager & Ashauer "Modelling survival under chemical stress" (2018 e-book)** — worked `k_D`
  values for the canonical species.

Realistic **distinct-species pool with a published `k_D` AND an AmP `k_M`: ~20–35.** That is the
`n` for the regression — modest but in COMADRE's regime (COMADRE `k_M`↔recovery ran at comparable n).

## 3. The assembled candidate roster (already on disk)

`data/external/guts_kd_amp_candidates.csv` — **25 standard GUTS test species confirmed in AmP**,
each with AmP `k_M` and `L_m`, and **blank `guts_kd_perday` / `guts_source` / `chemical` columns to
fill** from the §2 papers. This is the skeleton the next session populates.

Dynamic range across the 25 (the reason the test can work at all):
- `k_M`: **0.00113 → 2.06 /day — ~1,820× spread** (Mytilus edulis → Caenorhabditis elegans).
- `L_m`: **0.005 → 16.0 cm — ~3,200× spread** (C. elegans → Clarias gariepinus).
- Phylogenetically broad: cladocerans, amphipods, isopods, mysids, insects (Diptera, Collembola),
  gastropods, bivalves, nematodes, annelids, decapods, teleosts, amphibians.

## 4. The make-or-break caveat (read before building)

**Body-size collinearity is strong.** Across the 25 candidates,
**cor(log `k_M`, log `L_m`) = −0.668** — smaller animals have faster maintenance, as expected from
pace-of-life. Because the test *must* control body size (both `k_M` and a TK-dominated `k_D` scale
allometrically), partialling out size will **eat a large share of the `k_M` signal**. This is
exactly the COMADRE pace-of-life situation: the **size-partialled residual** is the non-trivial
result, and it is the **lower-power** version of the test. A raw `k_M`↔`k_D` correlation that
vanishes under a size control would be the *expected* null, not a surprise — so the size control is
doing the real scientific work and must be reported as the headline, not an afterthought.

**Conceptual caveat — what `k_D` actually is.** The GUTS-RED *dominant* rate constant lumps
toxicokinetics (uptake/elimination) with damage recovery into one one-compartment rate; for many
chemicals it is **TK/elimination-dominated**, which tracks surface-to-volume (size), *not*
maintenance metabolism per se. So `k_D ≈ 1/k_M` is a genuine **hypothesis**, not a near-identity.
That cuts both ways and keeps the test honest:
- a **positive** size-partialled `k_M`→`k_D` link would be a strong, surprising, on-thesis result
  (maintenance pace sets the toxicodynamic damage timescale beyond mere size);
- a **null** is equally publishable — it says the damage timescale is set by toxicokinetics/size,
  not by maintenance pace, which *bounds* the dynamic claim rather than refuting the margin.

Either way it is non-circular and informative. This is the right next test.

## 5. Concrete next steps (for the session that builds it)

1. **Populate `guts_kd_amp_candidates.csv`** `k_D` columns from the §2 sources (start with Ashauer
   2016 + the openGUTS / EFSA worked values — the cleanest, most-cited). Record source + chemical
   per row; one species may carry several `k_D` (per chemical) → decide aggregation (geometric mean
   per species, or keep chemical as a random effect).
2. **Regression** (COMADRE template, `scripts/comadre_pgls*.jl` is the pattern): `log k_D ~ log k_M
   + log body_size`, then add phylogeny via PGLS using the dated tree machinery already built
   (`scripts/comadre_pgls_dated.jl`; reuse the VertLife/TimeTree download route of Part 6.2). Report
   the **size-partialled** `k_M` coefficient as the headline.
3. **Fallback if filled-in `n` < ~15** (Part 5 fallback): extract the **observed** erosion timescale
   (time-to-plateau of a raw effect time-course) for a handful of AmP species and compare to their
   independent `1/k_M`. Same non-circular logic, less power.

## 5b. Data-access status — the `k_D` numbers are SI-gated (blocker)

Attempted to harvest actual `k_D` values online (2026-06-13). **Every GUTS paper puts `k_D` in a
supplementary file, not the main text**, and those are unreachable with the available web tools:
- **PMC supplementary downloads** are now behind a **JavaScript proof-of-work challenge**
  (`cloudpmc-viewer-pow`) — `curl`/WebFetch get the "Preparing to download…" interstitial, not the file.
- **Publisher SI mirrors** (Springer `static-content`) return **403** (Akamai bot-wall).
- No **machine-readable compiled `k_D` dataset** (Zenodo/Dryad/GitHub CSV) was findable.
- Main-text **HTML tables** *are* reachable (that's how the species lists above were confirmed), but
  they carry only the model *structure*, not the fitted `k_D` numbers.

**`k_D` numbers were NOT filled in — deliberately. Do not approximate them; a wrong `k_D` poisons
the regression.** The scaffold (`guts_kd_amp_candidates.csv`) now carries, per species, the
**source paper(s) + chemical(s)** that report a `k_D`, with the `guts_kd_perday` column blank.

**Shopping list (drop into `~/Downloads`, mirroring the paywalled-PDF workflow in the handover
Part 7), ordered by yield:**
1. **Ashauer et al. 2016** *Sci Rep* 6:29178 — **SI PDF** `srep29178-s1.pdf` (Tables S14–S17: `k_D`
   for *Gammarus pulex* ×4 pesticides + the 5-species malathion SSD). Highest yield, broadest size span.
2. **Focks et al. 2018** *Ecotoxicology* 27:992 — **SI** ESM4 / Table S2 (`k_D` for *Asellus
   aquaticus, Caenis horaria, Chaoborus obscuripes, Cloeon dipterum, Plea minutissima* × neonics).
3. **Nickisch Born Gericke et al. 2022** *ETC* 41:1732 — **SI** Tables S1–S6 (`k_D` for 5 fish ×
   benzovindiflupyr; gives rainbow trout / fathead minnow / common carp).
4. **openGUTS** ring-test parameter outputs + **EFSA 2018 TKTD opinion** (*EFSA J* 16:5377)
   case-study `k_D` (Daphnia magna, Gammarus, fish) — for the canonical anchors.

Once any of these are local, filling the CSV is mechanical (read table → paste `k_D` + model variant).

## 5c. Upgrade from reading the Ashauer 2016 main text (PDF supplied 2026-06-13)

The main PDF (`~/Downloads/srep29178.pdf`) does **not** carry the `k_D` tables — those are in the SI
(Table **S14** = *G. pulex* ×4 pesticides; Tables **S16/S17** = the 5-species malathion SSD; **S18**
= the synthetic study). So the SI PDF is still required for the numbers. But the main text delivers
two things worth more than the four `G. pulex` values:

**(i) Our test is a recognised, explicitly-stated open problem — not our invention.** Discussion, p. 8:
> "Species traits such as metabolic rate, which scales with size, correlated with the dominant rate
> constant for a small set of chemicals … and the dominant rate constant could be related to the
> size of three different species. … predictions of species sensitivity based on phylogeny could be
> refined to predict GUTS parameters as proxy for sensitivity. Such new statistical models should
> then be tested on a wide range of species and chemicals."

That is precisely the AmP-`k_M` → GUTS-`k_D` (controlling size + phylogeny) test, named as the field's
own next step. Strong positioning; also means it **must be cited as prior art**, and we must say what
we add (AmP `k_M` specifically; many species; explicit phylogeny control; thesis framing as
maintenance-timescale).

**(ii) Better, richer data leads than the scattered GUTS SIs** (and the closest prior attempt):
- **Baas & Kooijman 2015** *Ecotoxicology* 24:657 — checked (abstract). **Resolves the novelty
  question in our favour.** They regressed AmP **specific somatic maintenance `[p_M]`** (same
  add-my-pet source as our `k_M = [p_M]/[E_G]`) against the **NEC** (no-effect *concentration*; a
  survival **threshold/potency**), NECs from US-EPA ECOTOX, 4 pesticides (chlorpyrifos, malathion,
  carbofuran, carbaryl). **That is the *threshold* side, not the *rate* side.** NEC is a
  concentration (potency); `k_D` is a rate (timescale). The two-timescale thesis is exactly this
  state-vs-dynamics split: **Baas & Kooijman validated the margin/state axis (maintenance→potency);
  our `k_M`→`k_D` test is the untested dynamic/rate twin.** Net effect: AmP maintenance is *already a
  proven predictor of toxic response* (de-risk), and the rate axis is open (novelty). Note: their
  outcome (NEC from ECOTOX) is non-DEB → non-circular, and the repo already holds
  `data/ECOTOX_Toxicity_Library.json`, so a maintenance→NEC *cross-check* is even runnable in-repo as
  a bonus replication of the state side.
- **Gergs, Kulkarni & Preuss 2015** *Environ. Pollut.* 206:449 — body-size-dependent TK/TD; `k_D`
  vs size for three species.
- **Rubach et al. 2010** *ETC* 29:2225 (toxicokinetics of **15 freshwater arthropod species**,
  chlorpyrifos) + **Rubach 2011** *Ecotoxicology* 21:2088 (traits → sensitivity). A **single-chemical,
  many-species TK-rate** dataset — likely the **highest-yield, most size/phylogeny-controlled** source
  of a `k_D`-analog (uptake/elimination rate constant) available. Check AmP overlap of those 15
  arthropods next.

**Revised source priority:** (1) Baas & Kooijman 2015 — prior-art check; (2) Rubach 2010 — the rich
15-species TK set; (3) Ashauer 2016 SI (S14/S16/S17); (4) Focks 2018, Gericke 2022 SIs.

## 5d. First real paired run — Rubach 2010 (the data is in; result + a key caveat)

Rubach 2010 SI wasn't needed — **Table 2 (main text) gives elimination rate constants `k_out`** for
all 15 species, and Table 1 gives fresh weight (allometry control) + EC50 (sensitivity). Extracted
and paired against AmP `k_M`. Overlap: **6 clean species-level matches** (*Asellus aquaticus, Cloeon
dipterum, Culex pipiens, Daphnia magna, Gammarus pulex, Notonecta maculata*) + *Procambarus*
(sp.→clarkii proxy). *Chaoborus* dropped — two AmP congeners differ 20× in `k_M`, proxy unreliable.
Paired data: `data/external/rubach2010_kM_kout_paired.csv`. First runnable `k_M`↔rate test in the repo.

**Result (n=6, and n=7 with the proxy):**
| | cor(log k_M, log k_out) | partial, controlling log weight |
| --- | --- | --- |
| n=6 species | **−0.52** | **−0.42** |
| n=7 +proxy | **−0.55** | **−0.58** |

Underpowered (n=6 → nothing significant; |r|<0.81 is n.s.), but the sign is **weakly NEGATIVE** —
the *opposite* of the naïve thesis direction (faster maintenance → faster dynamics), and it does
**not** collapse under the size control.

**The key caveat this exposes — and it reshapes the Part 5 design.** Rubach's `k_out` is the
**toxicokinetic elimination rate** (chemical clearance, from radiotracer), *not* the **toxicodynamic
damage-recovery rate**. The two-timescale thesis predicts maintenance pace sets the **TD
recovery/erosion** timescale (`λ ∈ [k_M, k_M·g]`), i.e. GUTS-*proper* `k_r` (damage recovery) — **not**
chemical elimination `k_e`. Rubach himself (p.9) finds TK is governed by lipophilicity/biotransformation,
explaining only ~28% of sensitivity via `k_out`; elimination did **not** track size or lipid simply.
So the null/negative here is the *expected* "off-target endpoint" result, not evidence against the thesis.

**Consequence for the powered test (important, revises Part 5):** the GUTS-RED **dominant rate
constant `k_D` is frequently TK(elimination)-dominated** — and TK is about chemical partitioning, not
organism metabolism. If `k_D ≈ k_e` (true for *slow*-eliminating chemicals like chlorpyrifos, `t95`
up to 143 d here), then `k_D` will *not* track `k_M` either, for the same reason `k_out` doesn't. **The
thesis-relevant quantity is specifically the TD recovery rate `k_r` (GUTS-proper), or `k_D` only where
it is TD-limited (fast-eliminating chemicals).** That is a *stricter* data requirement than Part 5
assumed: we need `k_r`-resolved fits (GUTS-proper), not just any GUTS-RED `k_D`.

**Revised recommendation:**
1. Treat Rubach `k_out` as the **TK-axis control** (it *should* be null vs `k_M` — and is). Good
   two-timescale hygiene: it shows the pipeline discriminates TK from TD.
2. For the real (TD) test, target **GUTS-proper `k_r`** or `k_D` from **fast-eliminating** compounds
   (where `k_D` is recovery-limited). Sources: Ashauer 2016 SI (check whether `k_D` there is TK- or
   TD-limited per chemical), GUTS-proper fits (Nyman 2012 propiconazole; openGUTS), EFSA 2018.
3. The Baas & Kooijman maintenance→NEC (threshold/state-axis) in-repo cross-check remains the cleanest
   *state*-side confirmation and needs no new data.

## 5e. Retrieving GUTS-proper `k_r` (the TD rate) — what came back

Chased the TD recovery rate `k_r` (the on-target endpoint identified in §5d).

**Reliable (open PMC, internally consistent):**
- **Nyman et al. 2012** (*Ecotoxicology* 21:1828, PMC3431474) — *G. pulex* × propiconazole, GUTS-proper:
  **`k_r` = 2.3 /d (SD model), 1.0 /d (IT model)**; killing rate `k_k`=0.0051 g nmol⁻¹d⁻¹; threshold
  z=311.6 nmol/g; TK `k_e`=6.9 /d. Here **`k_r` < `k_e`** → recovery is TD-limited, so this is a *clean*
  `k_r` (not elimination-contaminated). One solid `k_r` point; *G. pulex* AmP `k_M`=0.0726/d.

**Could NOT reliably retrieve (tool limitation — do not trust auto-extracted numbers):**
- **Ashauer 2013** (*ETC* 32:954, fish, PMC3615168) and **Gericke 2022** (*ETC* 41:1732, fish, PMC9328144):
  WebFetch + WebSearch **cross-contaminated the two** — both ETC fish-GUTS papers — and attributed
  *benzovindiflupyr* (registered ~2015, **anachronistic for a 2013 paper**) and the same 5-fish set to
  both. Numbers returned are untrustworthy. **These two need the actual PDFs read directly.** Note both
  are likely **GUTS-RED** (report a single lumped dominant rate `k_e`/`k_D`, not a separated `k_r`).

**Two structural facts that constrain the powered TD test (carry these):**
1. **`k_r` is chemical-specific, not a pure species trait.** Nyman's `k_r` is *the propiconazole value*
   for *G. pulex*. `k_M` is a pure species property. So a cross-species `k_M`↔`k_r` regression must
   **control chemical** (or use one chemical across many species — which for `k_r` does not exist).
2. **GUTS-proper `k_r` is overwhelmingly *G. pulex*-centric.** The TD-resolved literature is dominated
   by one species across chemicals; multi-*species* `k_r` (let alone single-chemical multi-species) is
   **not available off the shelf**. This confirms the handover's Part 5 worry: the *powered* dynamic
   test is genuinely data-starved on the TD axis. Realistic assembled `k_r` set ≈ 5–8 species,
   chemical-confounded — COMADRE-scale at best, and messier.

**Bottom line of the retrieval pass:** the *clean* on-target rate (`k_r`) exists but is scarce and
chemical-confounded; the *abundant* rate (Rubach `k_out`, GUTS-RED `k_D`) is TK-contaminated and
doesn't track `k_M` (§5d). The dynamic test is therefore **real but low-power**; the highest-value
move may be the **state-axis** (`k_M`→NEC, Baas & Kooijman replication, runnable in-repo) rather than
forcing the rate axis.

## 5f. State-axis cross-check RAN (in-repo, n=310) — maintenance→sensitivity is a body-size story

Replicated the Baas & Kooijman *state*-axis idea in-repo, no downloads. Pulled species-level acute
lethality (LC50/MOR + EC50/IMM, 24–96 h, water-conc units → µg/L) from the **raw ECOTOX ASCII dump**
(`data/ecotox/ecotox_ascii_03_12_2026/`) for B&K's **four chemicals** (chlorpyrifos, malathion,
carbaryl, carbofuran), aggregated to a per-species LC50 (median), and matched to AmP `k_M`/`[p_M]`/`L_m`.
**310 chemical×species pairs** matched to AmP (42–95 per chemical). Pipeline (committed, reproducible):
`data/ecotox/.../` → `scripts/extract_ecotox_acute.awk` → `data/external/ecotox_acute_4chem.csv` →
`scripts/state_axis_ecotox_amp.jl` → `data/external/state_axis_ecotox_amp_paired.csv`
(rate axis: `scripts/rubach2010_rate_axis.jl`).

**Result (log LC50 = potency; higher LC50 = LESS sensitive):**
| chemical (n) | cor(LC50, k_M) | cor(LC50, `[p_M]`) | cor(LC50, size L_m) | **partial(LC50, k_M ∣ size)** |
| --- | --- | --- | --- | --- |
| Chlorpyrifos (86) | −0.27 | −0.28 | +0.36 | **−0.01** |
| Carbaryl (87) | −0.35 | −0.32 | +0.47 | **−0.03** |
| Malathion (95) | −0.12 | −0.13 | +0.13 | **−0.03** |
| Carbofuran (42) | −0.19 | −0.16 | +0.42 | **+0.10** |
| **POOLED (310)** | **−0.27** | **−0.27** | **+0.35** | **−0.03** |

**Two clean findings:**
1. **The raw Baas & Kooijman signal replicates** — higher AmP maintenance (`k_M` *and* their `[p_M]`)
   → lower LC50 → **more sensitive** (pooled r ≈ −0.27, n=310, same sign in all 4 chemicals). At face
   value the state axis shows a real maintenance↔sensitivity link.
2. **…but it is entirely body-size-mediated.** Larger `L_m` → higher LC50 (+0.35), and the **partial
   correlation of maintenance with sensitivity controlling size collapses to ≈ 0** (−0.03 to +0.10
   across chemicals; pooled −0.03). With n=310 this is a *robust null*, not underpower: **maintenance
   adds nothing to sensitivity beyond body size.** Small animals are both higher-`k_M` and more
   sensitive; `k_M` rides the size correlation but carries no independent signal.

**Interpretation (the through-line of the whole session).** *Both* AmP-maintenance↔toxic-response
axes are dominated by the same body-size confound: the **rate** axis (§5d, Rubach `k_out`: weak/wrong-
sign, n.s.) **and** the **state** axis (here: real raw signal, but zero after size). The recurring
programme signature (rank-ish at face value, magnitude-modest, **specification-sensitive** — cf.
COMADRE `k_M`) holds again, in its strongest form yet: a well-powered size control **nulls** the
single-trait maintenance story. This does **not** refute the margin model (which uses across-axis
capacity weighting + MoA routing, not raw maintenance), but it **bounds** any "maintenance pace as
master predictor of toxic response" claim — the distinctive content has to live in the *across-axis
capacity weighting*, not in `k_M` alone.

**Caveats (honest):** LC50 is a cruder sensitivity proxy than B&K's curated NEC (multi-source,
multi-duration aggregation adds noise); `L_m` is a DEB output collinear with `k_M` by construction, so
the size control is conservative and may *over*-remove; and the LC50~size effect may partly proxy a
fish-vs-invertebrate taxonomic split (a phylogeny/taxon control is the natural next refinement). But
the direction is unambiguous and consistent across 4 chemicals.

## 6. Files
- `data/external/amp_kM_roster.csv` — all 7,335 AmP species with `k_M`, `L_m`, `g`, `A0` (the full
  AmP side; regenerate via the inline script if needed).
- `data/external/guts_kd_amp_candidates.csv` — the 25-species curated skeleton with blank `k_D`.
- `data/external/rubach2010_kM_kout_paired.csv` — **the first real paired set** (AmP `k_M` vs Rubach
  `k_out`, 6 species + 1 proxy); the TK-axis control (see §5d).
- `data/external/ecotox_acute_4chem.csv` — 4,473 species-level acute LC50/EC50 records (4 chemicals),
  extracted from the raw ECOTOX ASCII dump via `/tmp/extract_ecotox.awk`.
- `data/external/state_axis_ecotox_amp_paired.csv` — **310 chemical×species** LC50 pairs matched to AmP
  `k_M`/`p_M`/`L_m`; the state-axis result (§5f). The raw maintenance↔sensitivity signal that *nulls*
  under a body-size control.
- Sources to mine for `k_D`: Ashauer 2016 (PMC4933929), Nickisch Born Gericke 2022 (PMC9328144),
  Focks 2018 (PMC6132984), openguts.info ring-test, EFSA 2018 TKTD opinion (*EFSA J* 16:5377).
