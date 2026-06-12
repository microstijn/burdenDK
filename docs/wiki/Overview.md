# Overview — what are we doing?

[← Home](Home.md) · next: [How it works →](Pipeline.md)

## The problem

Regulatory risk assessment is dominated by **thresholds**: a concentration is
compared to a no-effect or effect level, and a cell/sample is flagged if it
exceeds it. That is useful but incomplete. Chronic, cumulative, persistent,
mixture-mediated pressure can **erode an organism's capacity to recover** long
before any single endpoint is exceeded. A population that looks "safe" on every
threshold map can still be one heat-wave or spill away from collapse because its
buffering capacity has been quietly spent.

TwoTimescaleResilience models that erosion directly. Instead of asking *"is the
concentration above a line?"* it asks *"how much has chronic pressure narrowed
this species' room to absorb the next shock, and by how much does that amplify
the shock's impact?"*

## The two timescales

The framework's name is its central idea. Two processes run at very different
speeds, and separating them is what makes the model tractable.

1. **Slow — margin erosion.** Over months to years, background pressure and
   retained chemical burden reduce a species' **adaptive margin** `A` (its
   physiological room to manoeuvre). This is the "memory" timescale.
2. **Fast — perturbation response.** When an acute event hits, the organism's
   recovery is governed by a **restoring force** `λ` that depends on the *current*
   margin. A depleted margin means a weaker restoring force and a larger,
   longer-lasting burden from the same event.

Because the slow variable is approximately constant during a fast event, the
burden of an acute perturbation works out to be proportional to `1/λ` (a Fubini
argument over the event; see [Model equations](Equations.md)). The headline
output is the **amplification factor**

$$ F = \frac{\lambda(A_0)}{\lambda(A_t)} $$

— how much larger the acute burden is for a chronically-stressed organism
(margin `A_t`) than for a pristine one (margin `A_0`). `F = 1` means no
amplification; `F > 1` means chronic pressure has made the same shock worse.

This two-timescale separation is also what lets the framework scale to large
spatial rasters without integrating full ecosystem Jacobian matrices per cell.

## Capacity – pressure – memory

The model is organised around three sources of information:

- **Capacity** — *what the organism can absorb.* Derived from **Add-my-Pet (AmP)**
  Dynamic Energy Budget parameters: the baseline margin `A0`, the per-axis
  sensitivities, and the recovery-rate bounds `λ_min`, `λ_max`.
- **Pressure** — *what the environment applies.* Derived from **EPA ECOTOX**
  toxicity data (NOEC / EC50) plus ambient concentration fields, routed to four
  DEB-informed physiological axes (assimilation, maintenance, growth,
  reproduction) by the chemical's effect code.
- **Memory** — *what persists.* A compound-retention recurrence carries internal
  burden `B_t` forward in time via a retention factor `ρ` and a
  bioaccumulation/magnification factor `K`.

These feed the response chain `C → B → x → E_axis → Q → A → λ → F`, detailed in
[How it works](Pipeline.md).

## The Canguilhem framing

In the manuscript, capacity–pressure–memory is read through Georges Canguilhem's
*The Normal and the Pathological*. The load-bearing distinction is between
**normativity** — the capacity to *establish new viable norms* — and
**homeostasis** — return to a fixed setpoint. The model's mathematics (a margin
that narrows, reducing the capacity to re-establish function) matches
"capacity to re-establish viable norms" better than "defence of a fixed normal".
An older **defence–attack–memory** vocabulary survives only as code shorthand;
the scientific framing is capacity–pressure–memory.

## What this framework *is*

- A **physiologically-structured vulnerability index**, DEB-*informed*.
- A way to turn AmP capacity + ECOTOX pressure + exposure memory into a
  continuous, threshold-free amplification field over species and space.
- A spatial workflow that clusters that field into **vulnerability regimes**
  (relative descriptions, not safe/unsafe classes).

## What this framework is *not*

- **Not a full DEB / DEBkiss / DEBtox implementation.** It borrows DEB's process
  *axes* and parameters, but drops the dynamical reserve/structure/maturity/κ
  state equations. The four axes are coarse, DEBtox-compatible coordinates for
  *where* a stressor acts — not a conserved-flux energy budget. (For this reason
  "mechanistic" overclaims; "physiologically structured" is the honest register.)
- **Not a synergism/antagonism model.** Mixtures are combined with explicit
  *assumptions* (TU, IA, grouped CA-then-IA), never fitted interaction
  coefficients. See [Mixture-effect models](../mixture_effect_models.md).
- **Not threshold-based.** The spatial layer forbids exceedance features by
  construction.

For the important caveats — including the fact that the amplification factor is a
one-dimensional index (the allocation fraction κ originally; the energy investment
ratio `g` after the recovery-floor fix) and is not yet externally validated — see
[Limitations & open questions](Limitations-and-Open-Questions.md).
