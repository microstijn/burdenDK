# Per-Axis Resilience (Idea B)

The scalar tests ([COMADRE External Validation](COMADRE-External-Validation.md)) used one model rate against one
demographic quantity. The margin-first reframe claims the **per-axis** structure (which
DEB process is impaired) carries information a scalar throws away. This is the
multi-dimensional test — **pre-registered** before computing (repo:
`docs/notes/comadre_peraxis_prereg.md`) to avoid fishing.

## Mapping (pre-registered)
**Model rates** (AmP, reference T, log): `k_M` maintenance, `R_i` ultimate reproduction,
`r_B` von Bertalanffy growth.
**Demographic components** (COMADRE `matA`, Â = A/λ₁, Stott et al. 2011 first-timestep
bounds, log₁₀):

| component | metric | definition |
| --- | --- | --- |
| recovery | damping ratio | \|λ₁\|/\|λ₂\| |
| compensation | reactivity (P̄₁) | max column sum of Â |
| resistance | attenuation (P̲₁) | min column sum of Â |

## Result — generation-time-controlled partial rank correlations (n=197)

| rate \ component | resistance | compensation | recovery |
| --- | --- | --- | --- |
| `k_M` (maint.) | +0.313 ** | −0.300 ** | +0.260 ** |
| `R_i` (repro.) | −0.597 ** | **+0.773 \*\*** | −0.092 |
| `r_B` (growth) | +0.524 ** | −0.590 ** | +0.326 ** |

**The reproduction axis maps to the compensation component.** `R_i` → compensation
(+0.773, robust to pace and to body mass `Ww_i`: +0.775) is the **single strongest
association in the entire COMADRE validation** — far exceeding any scalar correlation.
`k_M`/`r_B` sit at the persistence pole (predict resistance/recovery, anti-predict
compensation). The dominant structure is a **resistance↔compensation trade-off** with
recovery more independent — a genuine multi-axis picture.

Pre-registered predictions: R_i↔compensation **PASS** (both directions); recovery's
strongest rate is `r_B` **PASS**; `k_M`'s strongest component was resistance not
recovery (marginal "fail", both significant).

## Is it more than raw fecundity? (robustness)
The strongest cell is **partly mechanical**: reactivity is fertility-driven and `R_i`
is a fecundity count, so a correlation is partly built-in. Tested with reproduction
descriptors that are *not* fecundity magnitude (gen-controlled, n=193):

| predictor → compensation | partial rank ρ |
| --- | --- |
| `a_p` (timing) \| gen | +0.359 ** |
| `R_i/Ww_i` (mass-specific) \| gen | +0.628 ** |
| `R_i` \| gen, `a_p` | +0.793 ** |
| `a_p` \| gen, `R_i` | +0.442 ** |

(`kap_R` not usable: AmP default 0.95 for ~97% of species.) **This substantially allays
the mechanical concern** — size-free fecundity and (non-mechanical) reproduction timing
*each independently* predict compensation. The map is a real, multi-faceted
reproduction-strategy signal, not reducible to the `R_i`↔reactivity identity. *Logged
pre-reg deviation:* `a_p`→compensation came out positive (+0.36), not the predicted
negative — interpretation flagged for follow-up.

## The honest framing
`R_i` is predicted from the DEB **energy budget** (lab/reference) while COMADRE
reactivity is **independent field** demography. A 0.77 match across aphids→whales,
holding pace and size, validates that **DEB-predicted reproduction tracks field
reproductive demography** — meaningful cross-dataset corroboration, but a convergent
measurement of reproductive strategy, not proof of an abstract "margin." Repo notes:
`docs/notes/comadre_peraxis_validation.md`.

## Sources
- Stott, Townley, Hodgson (2011) *Methods Ecol. Evol.* 2:946–955.
- Capdevila et al. (2020) *Trends Ecol. Evol.* 35:776–786.
