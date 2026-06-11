# TwoTimescaleResilience — Review & Research Agenda

*Working note, 2026-06-11. Captures the current review of the framework, the specific finding about the AmP→capacity mapping, and a forward-looking study/work program. Intended as a `docs/` reference and as briefing material for Claude Code sessions.*

---

## Part 1 — What we have so far

### The framework, in one breath
TwoTimescaleResilience (`burdenDK`) is a species-aware Julia framework that maps chronic background stress onto erosion of a physiological *adaptive margin*, models how that erosion weakens a *restoring force*, and outputs an *amplification factor* for the burden of a later acute event. The two-timescale separation (slow margin erosion vs. fast perturbation response) is what lets it avoid Jacobian ecosystem matrices over large rasters; the Fubini argument supplies burden ∝ 1/λ. Capacity comes from Add-my-Pet (AmP), pressure from ECOTOX, and exposure history from a compound-memory recurrence with retention ρ and bioaccumulation K. The manuscript framing is Canguilhemian capacity–pressure–memory; defence–attack–memory is retained only as code shorthand.

### Assessment of the project as documented
The conceptual core is genuinely good, and the project discipline is unusually strong: the standing invariants (no κ, no gain/`burden_to_margin` knobs, threshold-free outputs, no safe/unsafe language) protect the model's honesty in exactly the places it could otherwise inflate effects; the Canguilhem reframing is substantive rather than cosmetic, because "capacity to re-establish viable norms" matches what the math does (margin narrowing) better than "defence of a fixed normal"; and the append-don't-delete addendum structure keeps the context file truthful over time.

The main risks identified:

- **Documentation ahead of code.** Several sections hedge between "done," "implemented/documented," and "planned/prompted." The species-archetype database is explicitly "confirm before relying on files," yet later treated as part of a complete Stage 1–4 pipeline. The project's own recommended next step is a source audit — that instinct is correct and should be acted on first.
- **DEB-adjacency is the sharpest scientific exposure** (expanded below).
- **Proxies carrying weight on thin evidence:** ρ and K are class-default modeling proxies; `effect_code` is a stand-in for mode of action in the grouped CA-then-IA model. Both are flagged honestly in the YAML but must stay loud in the manuscript.
- **The mixture-model demo has a quiet null result:** in the large calibrated demo all three mixture models are numerically identical because no axis/month has more than one contributing compound. The capability is only actually exercised by the small overlap demo.
- **Minor:** top-level `last_updated` (2026-05-28) lags the newest addendum and filename (2026-06-10).

### The DEB question — can we use DEB this way?
Verdict: **yes, as DEB-*informed* phenomenology, not as DEB.** Three layers, very different footing:

- **Axes (strong ground).** Using assimilation / maintenance / growth / reproduction as a coordinate system for *where* a stressor acts is exactly what DEBtox does via physiological mode of action. Borrowing the target classification has direct precedent.
- **Dynamics (given up — fine if stated).** Dropping reserve / structure / maturity / κ-allocation means no conserved-flux fidelity. Acceptable for a vulnerability index, but the word "mechanistic" in the current `scientific_identity` overclaims; "physiologically structured" is the honest register.
- **Parameters (the load-bearing step).** AmP parameters do not automatically keep their DEB meaning once removed from the dynamical system. The whole "AmP-derived capacity" claim stands on the explicit map from {p_Am, p_M, κ, v} to A0 and the α-axes.

### The mapping — analysis of `amp_library.jl`
The pure DEB compound-parameter derivations are textbook-correct, but the step that builds the framework's quantities has a dimensional problem with real numerical consequences.

**Correct and worth keeping**
- `E_m = p_Am / v` — standard maximum reserve density. ✓
- `L_m = κ·p_Am / p_M` — standard maximum structural length. ✓
- `A0 = E_m` — conceptually the *best* available choice, since reserve genuinely is the buffer in DEB.
- The identity `lambda_max / lambda_min = 1/κ` holds exactly, so the comment that κ ∈ (0,1) guarantees `max > min` is right, and removing the safety catch is justified.

**Problems, in priority order**
1. **Dimensionally heterogeneous α-axes.** `alpha_A = 1/E_m` (cm³/J), `alpha_M = p_M/(κ·p_Am)` = 1/L_m (1/cm), `alpha_G = κ` (dimensionless), `alpha_R = 1−κ` (dimensionless). Dotting these with a dimensionless stress vector and subtracting from A0 (J/cm³) does not combine. This is where a quantitative referee stops.
2. **Numerical inertness.** A0 is O(thousands) while every α·s term is O(1) or smaller for realistic stress, so A_t ≈ A0, λ(A_t) ≈ λ(A0), and F ≈ 1 for almost every species. The diagnostics confirm it: median `min_x_collapse` ≈ 17,700 (max ≈ 3×10⁶), while stress x is normally O(1). The amplification mechanism barely fires.
3. **Fmax collapses to κ alone.** With `KA = 0.3·A0`, the half-saturation term is `A0/(KA+A0) = 1/1.3` for every species, so `Fmax = 1 + (1/κ − 1)/1.3`. Everything else cancels. Verified against the table (Thalia democratica → κ = 0.10 exactly). "High amplification potential" therefore means "low κ" — a strong claim needing ecological defence, possibly backwards.
4. **`KA = 0.3·A0` is a hidden gain parameter.** The 0.3 is unjustified and only sets the 1/1.3 constant; it violates the project's own no-knob invariants.
5. **Minor:** no `κ ≥ 1` guard in the validity check; `A0 = E_m` inherits body-size scaling through {p_Am}'s zoom factor z, though Fmax does not.

**Fix direction.** Nondimensionalize per axis: `A_t = A0·(1 − Σ_i w_i · g_i(s_i))`, with each `g_i` a bounded dimensionless impairment (the existing `E_axis = x/(1+x)` shape) and weights `w_i` defensibly tied to the κ-rule split. This restores dimensional coherence, keeps A_t responsive at realistic x, removes the rogue 0.3, and breaks the κ-only collapse in Fmax.

### Tooling
Claude Code via the native installer (no Node.js; runs natively on Windows or under WSL2) plus the VS Code extension. Before use: move the working copy off the OneDrive-synced path (or pause sync) to avoid conflict copies, and commit a clean baseline so every agent change is a reviewable, revertable diff. High-value move: distill the YAML's invariants into a repo-root `CLAUDE.md` so the guardrails load automatically each session.

---

## Part 2 — What to study, and what to do

This is a launching point, not an exhaustive plan. The aim of the next phase is to convert "I think the framework is sound" into "I have verified it, and I can defend every place it leans on DEB." Study and action are interleaved on purpose — the reading directly feeds specific fixes.

### What to study

**DEB foundations, focused on the parameter question.** Kooijman's DEB theory: the κ-rule, reserve dynamics, the compound parameters {p_Am}, [p_M], v, κ, what [E_m] means, and specifically its body-size scaling through the zoom factor z. The question to hold throughout: *is reserve density a legitimate buffer proxy, and what would A actually do if you let it move dynamically?*

**DEBtox / DEBkiss and physiological mode of action.** Your routing of ECOTOX `effect_code` to axes is the pMoA assignment problem, already formalized in this literature (Jager and collaborators). This tells you whether the routing is defensible and how others justify it — directly relevant to making the `effect_code`-as-MoA proxy honest in the manuscript.

**Mixture toxicity theory.** Concentration Addition vs. Independent Action, the conditions each assumes, and two-stage / grouped models (which you already implement as grouped CA-then-IA). Enough to defend the null-mixture-model framing and to state precisely why these are not synergism/antagonism.

**Toxicokinetics for the memory model.** One-compartment TK, BCF/BAF, elimination rate constants — the route from ρ and K as class defaults to measured, defensible parameters.

**Canguilhem, *The Normal and the Pathological*.** Enough to make the normativity framing precise: the distinction between normativity (capacity to set new norms) and homeostasis (return to a setpoint) is the load-bearing philosophical claim and should not be decorative.

**Nondimensionalization craft.** Buckingham Pi and the practice of building a dimensionally coherent erosion equation — the immediate technical skill for the A_t fix.

### What to do (sequenced)

1. **Tidy and audit.** Fix the metadata date, then run the source audit: verify every function in the documentation inventory exists in `src/` with the documented signature, and flag documented-but-missing items (the archetype DB first). Make the self-description true before building on it.
2. **Pin current behavior with tests.** Before changing anything, write regression tests for the DEB identities: `lambda_max/lambda_min == 1/κ` and `Fmax == 1 + (1/κ − 1)/1.3`. Add the missing `κ ≥ 1` guard. This documents the as-built model so the refactor is auditable.
3. **Implement the nondimensional A_t** on a branch, with the explicit constraints "F must respond at realistic x" and "Fmax must not be a pure function of κ" written into the task. Re-run the response-capacity diagnostics and confirm both.
4. **Check the body-size confound empirically.** Correlate A0 and Fmax against a size variable (L_m or ultimate weight) across the species set; report the result either way.
5. **Sensitivity analysis.** Quantify how much outputs move under ρ/K choices, the KA choice, and mixture-model choice. This is reviewer-proofing for the proxy assumptions.
6. **Kooijman outreach** with three precise questions, not a general blessing: (a) are the four process axes defensible as a vulnerability coordinate system *without* the dynamics; (b) is the A0/α-axes map a legitimate use of AmP parameters once extracted from the model; (c) which DEB notation must be avoided to not imply unimplemented state variables.
7. **Build the 3×3×12 monthly-memory diagnostic with deliberate same-axis overlap**, so it actually exercises the mixture models — directly answering the null-result gap, not repeating it.
8. **Hold the deferrals.** Physiological condition memory Z_t, DEBtox scaled damage D_t, and real raster ingestion stay deferred until the margin equation, mixture diagnostics, and clustering pipeline are stable — consistent with the existing invariants.

### The single question that orders all of this
Does the framework still produce meaningful, defensible amplification once the margin equation is dimensionally coherent and the Fmax–κ collapse is broken? Everything in the program above either protects that result or tests whether it survives. If it does, you have a publishable mechanistic-vulnerability framework. If it doesn't, you have found that out cheaply, before the manuscript hardened around it.
