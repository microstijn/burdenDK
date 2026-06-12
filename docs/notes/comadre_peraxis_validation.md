# External validation #1c — per-axis DEB rates ↔ demographic resilience (Idea B)

*Working note, 2026-06-12. Idea B of the roadmap. Tests the pre-registered mapping
(`comadre_peraxis_prereg.md`, frozen before computing): do distinct DEB process
rates predict distinct COMADRE demographic-resilience components? This is the test
of the **margin-first multi-dimensional** claim, vs. the scalar `k_M`↔recovery
result. Reproduce: `scripts/extract_comadre_recovery.jl` (now also reactivity +
attenuation) → `scripts/extract_amp_reproduction_rates.jl` →
`examples/comadre_peraxis_validation.jl`. n=197 species (193 with generation time).*

## Metrics
- **Model rates** (AmP, ref T, log): `k_M` maintenance, `R_i` ultimate reproduction,
  `r_B` von Bertalanffy growth.
- **Demographic** (COMADRE `matA`, Â=A/λ₁, Stott et al. 2011, log₁₀):
  recovery = damping ratio; compensation = reactivity (max col sum); resistance =
  attenuation (min col sum).

## Result — gen-controlled partial rank correlations

| rate \ component | resistance | compensation | recovery |
| --- | --- | --- | --- |
| `k_M` (maint.) | +0.313 ** | −0.300 ** | +0.260 ** |
| `R_i` (repro.) | −0.597 ** | **+0.773 \*\*** | −0.092 |
| `r_B` (growth) | +0.524 ** | −0.590 ** | +0.326 ** |

(Raw correlations show the same pattern; R_i→compensation is +0.739 raw, +0.773
gen-controlled, **+0.775 controlling generation time AND body mass `Ww_i`**.)

**Pre-registered prediction outcomes:**
- **P1 PASS** — `R_i`'s strongest component is **compensation** (+0.773).
- **P1b PASS** — compensation's strongest rate is **`R_i`** (vs −0.30 `k_M`, −0.59 `r_B`).
- **P2b PASS** — recovery's strongest rate is `r_B` (+0.326; within pred. k_M-or-r_B).
- **P2 fail** — `k_M`'s strongest component is resistance (+0.313) just above recovery
  (+0.260); both significant and positive, but not the predicted column-max.

## What it means
1. **The per-axis structure is real and strong — the reproduction axis maps to the
   compensation component.** `R_i`→compensation (+0.77, robust to pace and size) is
   the single strongest association in the entire COMADRE validation, far exceeding
   any scalar recovery correlation. The multi-dimensional margin-first claim has
   genuine external support: which DEB process is impaired predicts *which kind* of
   demographic resilience, not just an undifferentiated rate.
2. **`k_M`/`r_B` sit at the persistence end.** Maintenance and growth rates predict
   **resistance** and **recovery** positively and **compensation** negatively — the
   slow/persistent pole. `k_M`'s recovery signal (+0.26**) replicates the scalar
   result on this table, though resistance edges it (P2 "fail" is marginal).
3. **The dominant structure is a resistance↔compensation trade-off, and the rates
   sort along it.** Resistance and compensation are near-mirror images (`R_i`: −0.60
   / +0.77; `r_B`: +0.52 / −0.59). Recovery (damping) is the more independent third
   component. So the data resolve a coherent fast-reproduction (compensation) vs
   slow-maintenance/growth (resistance) axis, with recovery partly orthogonal —
   exactly a multi-axis picture, not one scalar.

## The honest caveat (must accompany the headline)
**The `R_i`↔compensation link is partly mechanical.** Reactivity (max column sum of
Â) is dominated by the fertility (`matF`) entries, and `R_i` is a fecundity rate —
both are "reproductive output," so a correlation is expected by construction. What is
*not* tautological: `R_i` is predicted from the DEB **energy budget** ({p_Am, κ, v,
…}, lab/reference conditions) while COMADRE reactivity comes from **independent field
projection matrices**. A 0.77 rank match across aphids→whales, holding pace and body
size, therefore validates that **DEB-predicted reproduction tracks field
reproductive demography** — a meaningful cross-dataset correspondence, but it is a
convergent measurement of fecundity, not proof of an emergent abstract "margin." The
manuscript must frame it as the former. (The `k_M`/`r_B`↔resistance/recovery cells
are less mechanically coupled and so carry more independent weight, but are weaker.)

## Verdict
**The pre-registered core passed: the reproduction axis specifically and robustly
predicts the compensation component** — the first multi-dimensional external result,
and stronger than anything scalar. Tempered by the mechanical-coupling caveat for
that specific cell. The maintenance/growth↔resistance/recovery side is directionally
as predicted but weaker and noisier. Next: (i) a less fecundity-mechanical
reproduction predictor (e.g. `kap_R`, or reproduction *timing* `a_p`) to test whether
the reproduction↔compensation map holds beyond raw fecundity; (ii) carry the dated
phylogeny (Idea A TODO) through to this matrix; (iii) decide which cells are
manuscript-grade given the coupling audit.

## Sources
- Stott, I., Townley, S., Hodgson, D.J. (2011) A framework for studying transient
  dynamics of population projection matrices. *Methods Ecol. Evol.* 2:946–955.
- Capdevila, P. et al. (2020) Towards a comparative framework of demographic
  resilience. *Trends Ecol. Evol.* 35:776–786.
- COMADRE — Salguero-Gómez et al. 2016, *J. Anim. Ecol.* (CC-BY).
