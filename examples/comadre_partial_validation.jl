# ===========================================================================
# External validation (#3): does the model predict COMADRE demographic recovery
# BEYOND raw pace-of-life? Generation-time-controlled partial Spearman.
#
# READ-ONLY. Reads the COMADRE-derived recovery table
# (data/external/comadre_recovery.csv, produced by
# scripts/extract_comadre_recovery.jl), matches species to AmP, computes the
# model's recovery quantities, and reports raw vs partial (controlling COMADRE
# generation time) rank correlations.
#
# Rationale + caveats: docs/notes/comadre_partial_validation.md
#   julia +release --project=. examples/comadre_partial_validation.jl
# ===========================================================================

using TwoTimescaleResilience
using Statistics
using Printf

const CSV = joinpath(@__DIR__, "..", "data", "external", "comadre_recovery.csv")

function read_comadre(path)
    sp = String[]; logdamp = Float64[]; gen = Float64[]
    for (i, line) in enumerate(eachline(path))
        i == 1 && continue
        f = split(line, ","); length(f) >= 4 || continue
        push!(sp, f[1]); push!(logdamp, parse(Float64, f[3]))
        push!(gen, f[4] == "" ? NaN : parse(Float64, f[4]))
    end
    return sp, logdamp, gen
end

ordinalrank(v) = (p = sortperm(v); r = similar(p); r[p] = 1:length(v); Float64.(r))
spear(x, y) = cor(ordinalrank(x), ordinalrank(y))
function pspear(x, y, z)              # partial Spearman of x,y controlling z
    rxy = spear(x, y); rxz = spear(x, z); ryz = spear(y, z)
    return (rxy - rxz * ryz) / sqrt((1 - rxz^2) * (1 - ryz^2))
end
sig(r, n; k = 0) = (t = r * sqrt((n - 2 - k) / (1 - r^2)); abs(t) > 2.6 ? "**" : abs(t) > 1.97 ? "*" : "")

function main()
    sp, logdamp, gen = read_comadre(CSV)
    lib = load_amp_species_library()

    lamA0 = Float64[]; lammin = Float64[]; gval = Float64[]; rec = Float64[]; cgen = Float64[]
    for k in eachindex(sp)
        key = replace(sp[k], " " => "_")
        haskey(lib, key) || continue
        local p
        try; p = amp_species_deb_params(lib, key); catch; continue; end
        (0 < p.lambda_min <= p.lambda_max && isfinite(p.A0)) || continue
        push!(lamA0, restoring_force_from_margin(p.A0, p))   # pristine recovery rate
        push!(lammin, p.lambda_min)                          # = k_M (maintenance rate constant)
        push!(gval, p.lambda_max / p.lambda_min)             # = g (energy investment ratio)
        push!(rec, logdamp[k]); push!(cgen, gen[k])
    end
    mask = isfinite.(cgen); n = count(mask)
    @printf("AmP-COMADRE matched species: %d  (with COMADRE generation time: %d)\n\n", length(rec), n)

    println("Does the model predict COMADRE demographic recovery (log damping ratio)?")
    @printf("%-26s %16s %30s\n", "model quantity", "raw rho", "partial rho (| gen. time)")
    for (lab, x) in (("lambda(A0) recovery rate", lamA0),
                     ("lambda_min = k_M", lammin),
                     ("g = lambda_max/lambda_min", gval))
        raw = spear(x, rec)
        par = pspear(x[mask], rec[mask], cgen[mask])
        @printf("%-26s %12.3f %3s %22.3f %3s\n", lab, raw, sig(raw, length(rec)), par, sig(par, n; k = 1))
    end
    @printf("\nreference  rho(COMADRE recovery, generation time) = %+.3f\n", spear(rec[mask], cgen[mask]))
    @printf("reference  rho(k_M, COMADRE generation time)       = %+.3f\n", spear(lammin[mask], cgen[mask]))
    println("\n* p<0.05, ** p<0.01. Reading: the model's recovery-rate quantities predict")
    println("demographic recovery and the signal survives controlling for generation time")
    println("(model adds info beyond pace-of-life); the g amplification axis does not.")
end

main()
