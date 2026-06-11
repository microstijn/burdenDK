# ===========================================================================
# AmP lambda/KA structure comparison (fork-support diagnostic)
#
# READ-ONLY. Does not change any model default. It re-runs the whole AmP
# library under three single-lever variants of the restoring-force structure
# and reports what the amplification factor F correlates with under each, so
# the "relative-vs-absolute margin" fork in
# docs/claude/TwoTimescaleResilience_source_audit_2026-06-11.md (Finding 3b)
# becomes a concrete comparison rather than an abstract question.
#
# Restoring force (all variants):
#     lambda(A) = lambda_min + (lambda_max - lambda_min) * A / (KA + A)
#     F(Q)      = lambda(A0) / lambda((1-Q)*A0)      # Q = fraction of margin lost
#
# Variants (each changes exactly ONE lever from the current baseline, so the
# attribution is clean):
#   S0  baseline (current)         KA = 0.3*A0          [expect F driven by kappa]
#   S1  KA absolute                KA = 0.3*median(A0)  [re-introduces A0 / reserve capacity]
#   S2  lambda_min absolute floor  lambda_min = median(lambda_min)
#                                                       [ratio driven by the species'
#                                                        real rate lambda_max = v/L_m, not kappa]
#
# kappa is alpha_axes[3]; the species size variable is L_m (auxiliary_metrics).
#
# Run:  julia +release --project=. examples/amp_lambda_structure_comparison.jl
# ===========================================================================

using TwoTimescaleResilience
using Statistics
using Printf

const C_KA = 0.3                       # the (current) KA-to-A0 ratio
const Q_LEVELS = (0.5, 1.0)            # 0.5 = moderate erosion; 1.0 = full (Fmax)
const OUTPUT_DIR = joinpath(@__DIR__, "..", "output", "amp_lambda_structure_comparison")

lambda_curve(A, lambda_min, lambda_max, KA) =
    lambda_min + (lambda_max - lambda_min) * A / (KA + A)

function F_amp(Q, A0, lambda_min, lambda_max, KA)
    A_t = (1 - Q) * A0
    return lambda_curve(A0, lambda_min, lambda_max, KA) /
           lambda_curve(A_t, lambda_min, lambda_max, KA)
end

function ordinalrank(v::AbstractVector)
    p = sortperm(v); r = similar(p); r[p] = 1:length(v); return Float64.(r)
end
spearman(x, y) = cor(ordinalrank(x), ordinalrank(y))

function collect_species(library)
    kappa = Float64[]; A0 = Float64[]; Lm = Float64[]
    lmin = Float64[]; lmax = Float64[]
    for (_k, rec) in library
        haskey(rec, "auxiliary_metrics") || continue
        haskey(rec["auxiliary_metrics"], "L_m") || continue
        local p
        try
            p = amp_record_to_deb_params(rec)
        catch
            continue
        end
        k = p.alpha_axes[3]
        L = Float64(rec["auxiliary_metrics"]["L_m"])
        (all(isfinite, (k, p.A0, L, p.lambda_min, p.lambda_max)) &&
         0 < k < 1 && p.A0 > 0 && L > 0 && 0 < p.lambda_min < p.lambda_max) || continue
        push!(kappa, k); push!(A0, p.A0); push!(Lm, L)
        push!(lmin, p.lambda_min); push!(lmax, p.lambda_max)
    end
    return (; kappa, A0, Lm, lmin, lmax, n = length(kappa))
end

# Compute F vectors for one structure across all species, at a given Q.
function F_for_structure(d, Q; ka_of, lmin_of)
    F = Vector{Float64}(undef, d.n)
    for i in 1:d.n
        KA   = ka_of(d.A0[i])
        lmin = lmin_of(d.lmax[i], d.lmin[i])
        F[i] = F_amp(Q, d.A0[i], lmin, d.lmax[i], KA)
    end
    return F
end

function report(d; io::IO = stdout)
    median_A0   = median(d.A0)
    median_lmin = median(d.lmin)

    # lever definitions
    ka_baseline(A0)  = C_KA * A0
    ka_absolute(A0)  = C_KA * median_A0
    lmin_keep(lmax, lmin) = lmin
    # absolute floor, clamped below lambda_max so the curve stays well-defined
    lmin_floor(lmax, lmin) = min(median_lmin, 0.999 * lmax)

    structures = (
        ("S0 baseline  (KA = 0.3*A0)",        ka_baseline, lmin_keep),
        ("S1 KA absolute (0.3*median A0)",     ka_absolute, lmin_keep),
        ("S2 lambda_min absolute floor",       ka_baseline, lmin_floor),
    )

    @printf(io, "Species compared: %d\n", d.n)
    @printf(io, "median A0 = %.1f   median lambda_min = %.5f\n", median_A0, median_lmin)
    n_clamped = count(i -> median_lmin >= 0.999 * d.lmax[i], 1:d.n)
    @printf(io, "S2 floor exceeds lambda_max for %d species (%.1f%%) -> those cannot amplify under S2\n\n",
            n_clamped, 100 * n_clamped / d.n)

    F0_ref = F_for_structure(d, 0.5; ka_of = ka_baseline, lmin_of = lmin_keep)  # reshuffle baseline

    println(io, "What does amplification F correlate with, under each structure? (Q = 0.5 erosion)")
    @printf(io, "%-32s %10s %10s %10s %14s %12s\n",
            "structure", "rho(F,kappa)", "rho(F,A0)", "rho(F,L_m)", "reshuffle vs S0", "median F")
    for (name, ka_of, lmin_of) in structures
        F = F_for_structure(d, 0.5; ka_of = ka_of, lmin_of = lmin_of)
        @printf(io, "%-32s %10.3f %10.3f %10.3f %14.3f %12.4f\n",
                name,
                spearman(F, d.kappa), spearman(F, d.A0), spearman(F, d.Lm),
                spearman(F, F0_ref), median(F))
    end

    println(io, "\nFull-erosion ceiling Fmax (Q -> 1):")
    @printf(io, "%-32s %10s %10s %10s %12s %12s\n",
            "structure", "rho(Fmax,k)", "rho(Fmax,A0)", "rho(Fmax,Lm)", "median Fmax", "max Fmax")
    for (name, ka_of, lmin_of) in structures
        F = F_for_structure(d, 1.0; ka_of = ka_of, lmin_of = lmin_of)
        @printf(io, "%-32s %10.3f %10.3f %10.3f %12.3f %12.3f\n",
                name,
                spearman(F, d.kappa), spearman(F, d.A0), spearman(F, d.Lm),
                median(F), maximum(F))
    end

    println(io, "\nReading the table:")
    println(io, "  rho near -1 in the kappa column  => F is (still) a kappa ranking  [vulnerability = allocation]")
    println(io, "  rho moving toward +-1 in the A0   => reserve capacity now drives vulnerability")
    println(io, "  reshuffle vs S0 << 1              => the species ranking genuinely changes")
    println(io, "\nThis is decision-support for the Finding 3b fork; it changes no model default.")
end

if abspath(PROGRAM_FILE) == @__FILE__
    d = collect_species(load_amp_species_library())
    report(d)
    try
        mkpath(OUTPUT_DIR)
        open(joinpath(OUTPUT_DIR, "lambda_structure_comparison.txt"), "w") do io
            report(d; io = io)
        end
        @printf("\nSummary written to %s\n", joinpath(OUTPUT_DIR, "lambda_structure_comparison.txt"))
    catch err
        @warn "Could not write summary file" exception = (err, catch_backtrace())
    end
end
