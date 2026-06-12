# ===========================================================================
# Validation for the lambda_min = k_M fix (see docs/notes/lambda_min_maintenance_rate.tex).
#
# READ-ONLY. Confirms that re-anchoring the slow recovery floor to the DEB
# somatic maintenance rate constant k_M = [p_M]/[E_G] breaks the kappa-collapse:
# the amplification factor F should now track the energy investment ratio g
# = (v/L_m)/k_M, not the allocation fraction kappa.
#
#   julia +release --project=. examples/amp_lambda_min_validation.jl
# ===========================================================================

using TwoTimescaleResilience
using Statistics
using Printf

function ordinalrank(v::AbstractVector)
    p = sortperm(v); r = similar(p); r[p] = 1:length(v); return Float64.(r)
end
spearman(x, y) = cor(ordinalrank(x), ordinalrank(y))

function collect_metrics(lib)
    kappa = Float64[]; gval = Float64[]
    Fhalf = Float64[]; Fmax = Float64[]
    ratio = Float64[]                      # lambda_max / lambda_min
    clamped = 0; total = 0
    for (_k, rec) in lib
        haskey(rec, "auxiliary_metrics") && haskey(rec["auxiliary_metrics"], "g") || continue
        local p
        try
            p = amp_record_to_deb_params(rec)
        catch
            continue
        end
        a = p.alpha_axes
        (all(isfinite, a) && p.A0 > 0 && 0 < a[3] < 1 && 0 < p.lambda_min <= p.lambda_max) || continue
        total += 1
        push!(kappa, a[3])
        push!(gval, Float64(rec["auxiliary_metrics"]["g"]))
        push!(Fhalf, amplification_from_margin(0.5 * p.A0, p))
        push!(Fmax, amplification_from_margin(0.0, p))
        push!(ratio, p.lambda_max / p.lambda_min)
        p.lambda_min == p.lambda_max && (clamped += 1)
    end
    return (; kappa, gval, Fhalf, Fmax, ratio, clamped, total)
end

(; kappa, gval, Fhalf, Fmax, ratio, clamped, total) = collect_metrics(load_amp_species_library())

@printf("Species: %d   (clamped g<1, no separation: %d = %.1f%%)\n\n",
        total, clamped, 100 * clamped / total)

println("Correlations (Spearman):")
@printf("  rho(F@0.5, kappa) = %+.3f      <- was -1.000 before the fix\n", spearman(Fhalf, kappa))
@printf("  rho(F@0.5, g)     = %+.3f      <- F should now track g\n", spearman(Fhalf, gval))
@printf("  rho(g,     kappa) = %+.3f      <- g must NOT be just kappa for the fix to add info\n",
        spearman(gval, kappa))

println("\nTimescale-separation ratio lambda_max/lambda_min:")
unclamped = ratio[ratio .> 1.0 + 1e-9]
@printf("  among unclamped species: should equal g.  rho(ratio, g) = %+.3f\n",
        spearman(ratio, gval))
@printf("  max |ratio - g| (unclamped) = %.3e\n",
        isempty(unclamped) ? NaN : maximum(abs.(ratio[ratio .> 1.0 + 1e-9] .- gval[ratio .> 1.0 + 1e-9])))

println("\nAmplification magnitudes:")
@printf("  F@0.5 : median %.4f   max %.3f\n", median(Fhalf), maximum(Fhalf))
@printf("  Fmax  : median %.4f   max %.3f   (was median 1.059, max 11.4)\n",
        median(Fmax), maximum(Fmax))
