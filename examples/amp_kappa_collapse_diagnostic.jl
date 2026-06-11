# ===========================================================================
# AmP kappa-collapse diagnostic
#
# Read-only diagnostic that quantifies, across the whole AmP species library,
# how much of the amplification factor F is determined by the somatic
# allocation fraction kappa (= alpha_axes[3]) ALONE under the current
# parameterization.
#
# Background. With the as-built mapping (see src/AmP_Translator.jl) the
# restoring-force curve uses lambda_max/lambda_min == 1/kappa and KA == 0.3*A0.
# Because KA is proportional to A0, the species scale A0 cancels out of
#
#     lambda(A) = lambda_min + (lambda_max - lambda_min) * A / (KA + A)
#
# and F(A_t) = lambda(A0)/lambda(A_t) becomes a function of (kappa, erosion)
# only. This script measures that collapse empirically and is the evidence
# behind the "relative-vs-absolute margin" fork in
# docs/claude/TwoTimescaleResilience_source_audit_2026-06-11.md.
#
# It is intentionally dependency-light (no plotting / DataFrames) so it runs
# fast. Run with the project's Julia (1.12.x):
#
#     julia +release --project=. examples/amp_kappa_collapse_diagnostic.jl
# ===========================================================================

using TwoTimescaleResilience
using Statistics
using Printf

const OUTPUT_DIR = joinpath(@__DIR__, "..", "output", "amp_kappa_collapse_diagnostic")

# Realized erosion levels Q (fraction of pristine margin lost): A_t = (1-Q)*A0.
const Q_LEVELS = (0.2, 0.5, 0.8)

# Spearman rank correlation via Pearson correlation of ordinal ranks.
function ordinalrank(v::AbstractVector)
    p = sortperm(v)
    r = similar(p)
    r[p] = 1:length(v)
    return Float64.(r)
end
spearman(x, y) = cor(ordinalrank(x), ordinalrank(y))

# Closed-form F predicted from kappa ALONE (KA = c*A0). If the realized F
# matches this to ~machine epsilon, then A0 and the other alpha-axes
# contribute nothing to F beyond kappa.
function F_predicted_from_kappa(kappa::Real, Q::Real; c::Real = 0.3)
    r = 1 / kappa
    num = 1 + (r - 1) / (1 + c)
    den = 1 + (r - 1) * (1 - Q) / (c + 1 - Q)
    return num / den
end

function run_diagnostic(; library = load_amp_species_library())
    kappa = Float64[]
    A0    = Float64[]
    wA    = Float64[]                                   # normalized assimilation weight
    Fmax  = Float64[]                                   # F at full erosion (A_t -> 0)
    F_Q   = Dict(Q => Float64[] for Q in Q_LEVELS)      # F at realized erosion Q

    for (_key, record) in library
        local p
        try
            p = amp_record_to_deb_params(record)
        catch
            continue
        end
        a = p.alpha_axes
        (all(isfinite, a) && p.A0 > 0 && p.lambda_min > 0) || continue
        push!(kappa, a[3])
        push!(A0, p.A0)
        push!(wA, a[1] / sum(a))
        push!(Fmax, amplification_from_margin(0.0, p))
        for Q in Q_LEVELS
            push!(F_Q[Q], amplification_from_margin((1 - Q) * p.A0, p))
        end
    end

    return (; kappa, A0, wA, Fmax, F_Q, n = length(kappa))
end

function report(d; io::IO = stdout)
    @printf(io, "Species with usable records: %d\n\n", d.n)

    println(io, "=== D3: is amplification a pure function of kappa? ===")
    @printf(io, "%-30s %14s %22s\n", "quantity", "Spearman(.,kappa)", "max|actual - pred(kappa)|")
    @printf(io, "%-30s %14.6f %22.3e\n", "Fmax (A_t -> 0)",
            spearman(d.Fmax, d.kappa),
            maximum(abs.(d.Fmax .- F_predicted_from_kappa.(d.kappa, 1.0))))
    for Q in Q_LEVELS
        @printf(io, "%-30s %14.6f %22.3e\n", "F at realized erosion Q=$Q",
                spearman(d.F_Q[Q], d.kappa),
                maximum(abs.(d.F_Q[Q] .- F_predicted_from_kappa.(d.kappa, Q))))
    end
    println(io, "\n  Spearman = -1 => F is rank-identical to kappa.")
    println(io, "  max|actual - pred(kappa)| ~ 0 => A0 and the other alpha-axes")
    println(io, "  contribute NOTHING to F beyond kappa (the KA proportional-to-A0 cancellation).")

    println(io, "\n=== D2: assimilation weight collapse (normalized alpha_A) ===")
    @printf(io, "median w_A = %.5f   mean w_A = %.5f\n", median(d.wA), mean(d.wA))
    @printf(io, "fraction of species with w_A < 0.01: %.1f%%\n", 100 * count(<(0.01), d.wA) / d.n)
    @printf(io, "fraction of species with w_A < 0.05: %.1f%%\n", 100 * count(<(0.05), d.wA) / d.n)

    println(io, "\n=== Context: kappa, A0, Fmax spread ===")
    @printf(io, "kappa : min %.3f  median %.3f  max %.3f\n", minimum(d.kappa), median(d.kappa), maximum(d.kappa))
    @printf(io, "A0    : min %.1f  median %.1f  max %.3g\n", minimum(d.A0), median(d.A0), maximum(d.A0))
    @printf(io, "Fmax  : min %.3f  median %.3f  max %.3f\n", minimum(d.Fmax), median(d.Fmax), maximum(d.Fmax))
end

if abspath(PROGRAM_FILE) == @__FILE__
    d = run_diagnostic()
    report(d)
    try
        mkpath(OUTPUT_DIR)
        open(joinpath(OUTPUT_DIR, "kappa_collapse_summary.txt"), "w") do io
            report(d; io = io)
        end
        @printf("\nSummary written to %s\n", joinpath(OUTPUT_DIR, "kappa_collapse_summary.txt"))
    catch err
        @warn "Could not write summary file" exception = (err, catch_backtrace())
    end
end
