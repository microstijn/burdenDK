# CLAUDE.md — burdenDK / TwoTimescaleResilience

Guidance for Claude Code (and humans) working in this repository.

## Julia version (read first)

**Use Julia 1.12.6 for this project.** The committed `Manifest.toml` is resolved
under `julia_version = "1.12.6"`. The machine default is the LTS (1.10.11), which
**cannot load this project** — it fails while precompiling `PrecompileTools` with
`UndefVarError: StaticData not defined` (a version-mismatch symptom, not a code bug).

Run everything through the `release` channel:

```powershell
julia +release --project=. <script.jl>          # run a script
julia +release --project=. test/runtests.jl     # run the test suite
julia +release --project=. -e 'using Pkg; Pkg.instantiate()'
```

`juliaup status` shows the available channels; `release` is 1.12.6 here. If you
see the `StaticData` precompile error, you are on the wrong Julia — switch to
`+release`, don't try to "fix" the package.

## Testing

- Prefer targeted test files over the full suite during development; heavy
  dependencies (CairoMakie, GeoMakie, NCDatasets) trigger large precompilation.
  Example: `julia +release --project=. test/test_amp_translator_identities.jl`
- The fast suite is `test/runtests.jl`; extended/plotting/example tests are gated
  behind `TTR_RUN_EXTENDED_TESTS`, `TTR_RUN_EXAMPLE_TESTS`, `TTR_RUN_PLOTTING_TESTS`,
  `TTR_RUN_DYNQUAL_TESTS` env vars (default off).
- See `docs/TESTING_STRATEGY.md`.

## Where the model actually lives

- Core math: `src/deb_axes.jl`, `src/reduced_deb_response.jl` (`A_t`, `λ_t`, `F_t`).
- AmP → capacity parameters are **precomputed offline** by `src/AmP_Translator.jl`
  (a standalone script, depends on `MAT`, reads `data/allStat.mat`, writes
  `data/AmP_Species_Library.json`). It is **not** `include`d in the module.
  `src/amp_library.jl` only *loads/validates* that JSON. To change the
  `{p_Am,p_M,κ,v} → A0/α-axes/KA` mapping you edit `AmP_Translator.jl` and
  regenerate the JSON — editing `amp_library.jl` does nothing to the mapping.
- Orientation docs: `docs/AI_CONTEXT_README.md`, `docs/PACKAGE_CAPABILITIES.md`,
  `docs/ARCHITECTURE_GRAPH.md` (re-audited 2026-06-11 — see
  `docs/claude/TwoTimescaleResilience_source_audit_2026-06-11.md`).

## Invariants — do not violate without explicit sign-off

- **No arbitrary tuning knobs.** Do not add parameters named/behaving like
  `kappa`/`κ` (as a free knob), `gain`, `response_scale`, or
  `burden_to_margin_multiplier`. (Note: the existing `KA = 0.3·A0` in
  `AmP_Translator.jl` is a known violation flagged for removal — see the audit.)
- **No thresholds in threshold-free features.** Avoid `_gt_`, `_lt_`, `threshold`,
  `exceedance`, `above`, `below` in spatial feature/regime code.
- **Keep memory layers distinct:** chemical memory `B_t`, physiological condition
  memory `Z_t` (`src/condition_buffer.jl`, opt-in/off by default), and DEBtox
  scaled damage `D_t` (unimplemented). Do not conflate them.
- **Mixtures are assumptions, not interactions.** TU / IA / grouped CA-then-IA only;
  no synergism, antagonism, or fitted interaction coefficients. Don't call them
  "interactions."
- Don't use DataFrames for lightweight internal data structures.
- Don't rewrite core math unless explicitly asked.
