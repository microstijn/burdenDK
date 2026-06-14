# ===========================================================================
# Cross-species CAPACITY probe — 5 fish, benzovindiflupyr (single chemical).
# Data: Nickisch Born Gericke et al. 2022 (ETC 41(7):1732-1741) SI -> per-species acute LC50
# (data/external/benzovindiflupyr_fish_survival.csv via scripts/extract_benzovindiflupyr_fish.jl).
#
# QUESTION (the capacity-weighting open question, in miniature). Across 5 fish exposed to ONE
# chemical, does the model's per-species AmP capacity predict relative sensitivity (LC50) -- and,
# the decisive control from the n=310 negative control (manuscript Sec. validation), does it survive
# a BODY-SIZE control? Benzovindiflupyr is an SDHI (mitochondrial complex-II inhibitor) -> energetic
# cost -> the MAINTENANCE axis; so the framework-relevant capacity is the maintenance-axis weighting.
#
# HONEST SCOPE UP FRONT: n=5 species. This is a PILOT, not a powered test -- at n=5 only |rho|>~0.9
# reaches p<0.05. A single chemical exercises ONE axis, so this probes single-axis cross-species
# capacity (a Sec.-7b-style extension to a new chemical class/taxon), NOT the full across-axis
# weighting. Read it as direction + whether body size dominates, exactly as the n=310 control warned.
#
#   julia +release --project=. examples/benzovindiflupyr_capacity_probe.jl
# ===========================================================================

using TwoTimescaleResilience
using Statistics, Printf

const FILE = joinpath(@__DIR__, "..", "data", "external", "benzovindiflupyr_fish_survival.csv")
const FISH = ["Cyprinodon_variegatus", "Cyprinus_carpio", "Lepomis_macrochirus",
              "Oncorhynchus_mykiss", "Pimephales_promelas"]

ordinalrank(v) = (p = sortperm(v); r = similar(p); r[p] = 1:length(v); Float64.(r))
spear(x, y) = cor(ordinalrank(x), ordinalrank(y))
function partial_spear(x, y, z)            # Spearman of x,y controlling z (rank-linear residuals)
    rx, ry, rz = ordinalrank(x), ordinalrank(y), ordinalrank(z)
    Z = hcat(ones(length(rz)), rz); resid(v) = v .- Z * (Z \ v)
    return cor(resid(rx), resid(ry))
end

# --- per-species acute LC50 from the long survival CSV ---
function lc50_table()
    rows = NamedTuple[]
    for line in eachline(FILE)
        (startswith(line, "#") || startswith(line, "species,")) && continue
        f = split(line, ","); length(f) == 6 || continue
        push!(rows, (sp = String(f[1]), sample = String(f[2]), time = parse(Float64, f[4]),
                     conc = parse(Float64, f[5]), nsurv = parse(Float64, f[6])))
    end
    lc = Dict{String,Float64}()
    for sp in FISH
        ac = filter(r -> r.sp == sp && startswith(r.sample, "acute"), rows)
        tmax = maximum(r.time for r in ac); concs = sort(unique(r.conc for r in ac))
        frac = [ (n0 = sum(r.nsurv for r in ac if r.conc == c && r.time == 0.0; init = 0.0);
                  nf = sum(r.nsurv for r in ac if r.conc == c && r.time == tmax; init = 0.0);
                  n0 > 0 ? nf / n0 : NaN) for c in concs ]
        l = NaN
        for i in 1:length(concs)-1
            if isfinite(frac[i]) && isfinite(frac[i+1]) && frac[i] >= 0.5 >= frac[i+1] && concs[i] > 0
                t = (frac[i] - 0.5) / (frac[i] - frac[i+1])
                l = exp(log(concs[i]) + t * (log(concs[i+1]) - log(concs[i]))); break
            end
        end
        lc[sp] = l
    end
    return lc
end

function main()
    lc = lc50_table()
    lib = load_amp_species_library()
    sp = String[]; LC50 = Float64[]; A0 = Float64[]; aMaint = Float64[]; kM = Float64[]; Lm = Float64[]
    println("=== per-species: sensitivity (LC50) vs AmP capacity + body size ===")
    @printf("  %-22s %9s %12s %10s %9s %8s\n", "species", "LC50", "A0", "alpha_maint", "k_M", "L_m(cm)")
    for s in FISH
        p = amp_species_deb_params(lib, s); am = lib[s]["auxiliary_metrics"]
        push!(sp, s); push!(LC50, lc[s]); push!(A0, p.A0); push!(aMaint, p.alpha_axes[2])
        push!(kM, am["k_M"]); push!(Lm, am["L_m"])
        @printf("  %-22s %9.2f %12.3g %10.4f %9.4f %8.2f\n", s, lc[s], p.A0, p.alpha_axes[2], am["k_M"], am["L_m"])
    end

    println("\n=== rank correlations with LC50 (lower LC50 = MORE sensitive), n=5 ===")
    @printf("  rho( body size L_m , LC50 )          = %+.2f   [+ = bigger fish less sensitive]\n", spear(Lm, LC50))
    @printf("  rho( k_M , LC50 )                    = %+.2f   [Sec.7b raw: higher k_M -> more sensitive => -]\n", spear(kM, LC50))
    @printf("  rho( A0 (capacity) , LC50 )          = %+.2f   [+ = more capacity, less sensitive]\n", spear(A0, LC50))
    @printf("  rho( alpha_maint (weighting) , LC50 )= %+.2f   [framework axis weighting; - = more vuln]\n", spear(aMaint, LC50))
    println("  -- controlling body size (the decisive Sec.-7b control) --")
    @printf("  partial rho( k_M , LC50 | L_m )            = %+.2f\n", partial_spear(kM, LC50, Lm))
    @printf("  partial rho( A0 , LC50 | L_m )             = %+.2f\n", partial_spear(A0, LC50, Lm))
    @printf("  partial rho( alpha_maint , LC50 | L_m )    = %+.2f\n", partial_spear(aMaint, LC50, Lm))
    @printf("  collinearity check  rho( A0 , L_m ) = %+.2f ;  rho( k_M , L_m ) = %+.2f\n", spear(A0, Lm), spear(kM, Lm))

    println("\nReading: at n=5 only |rho|>~0.9 is significant -- treat as DIRECTION, not proof. The")
    println("question is whether any capacity term predicts sensitivity BEYOND body size (the n=310")
    println("negative control showed single-trait k_M does NOT). The dimensionless axis WEIGHTING")
    println("(alpha_maint), unlike A0/k_M, is not a size proxy -- its partial-on-size is the cleanest")
    println("read of the framework's distinctive content. Single chemical = single axis; a real")
    println("across-axis weighting test needs multiple MoA x species. See docs/notes/.")
end

main()
