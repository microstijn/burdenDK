# ===========================================================================
# Build the AmP <-> GlobTherm matched table for the thermal-tolerance coherence
# check (external validation #5, "glob first"). GlobTherm is an INDEPENDENT
# physiological dataset (measured CTmax/CTmin), so it is a clean external probe of
# whether the AmP-derived capacity axis carries any thermal-tolerance signal.
#
# NOTE ON SCOPE (honest): the adaptive-margin model does not predict thermal limits
# directly (it uses the energetic capacity A0/k_M/lambda, not the Arrhenius thermal
# params). So this is a COHERENCE / general-resilience probe, not a test of the
# margin's erosion mechanism (that is the Scope-for-Growth follow-up). The thermal
# regime T_typical is extracted so the analysis can de-confound it (AmP rates carry
# the species' thermal regime, which also drives CTmax).
#
# STANDALONE — needs MAT (throwaway env). Self-contained.
#   julia +release -e 'using Pkg; Pkg.activate(mktempdir()); Pkg.add("MAT")'
#   julia +release --project=<that-env> scripts/extract_amp_for_globtherm.jl
#
# Reads:  data/allStat.mat, data/external/GlobalTherm.csv (Zenodo rec 4976423, CC-BY;
#         Bennett et al. 2018 Sci. Data)
# Writes: data/external/globtherm_amp_matched.csv (committed):
#   species,k_M,E_m,lambda_max,g,T_A,T_typical,Ww_i,CTmax,CTmin,lat,class
# ===========================================================================

using MAT, Statistics

const MATF = joinpath(@__DIR__, "..", "data", "allStat.mat")
const GT   = joinpath(@__DIR__, "..", "data", "external", "GlobalTherm.csv")
const OUT  = joinpath(@__DIR__, "..", "data", "external", "globtherm_amp_matched.csv")

# --- GlobTherm: Genus_species => (CTmax, CTmin, lat, class) ----------------
pf(x) = (x = strip(x); x == "" ? NaN : something(tryparse(Float64, x), NaN))
gt = Dict{String, NTuple{4, Any}}()
for (i, line) in enumerate(eachline(GT))
    i == 1 && continue
    f = split(line, ","); length(f) >= 43 || continue
    key = string(strip(f[1]), "_", strip(f[2]))
    gt[key] = (pf(f[4]), pf(f[23]), pf(f[17]), strip(f[43]))   # Tmax, tmin, lat_max, Class
end
println("GlobTherm species: ", length(gt))

# --- AmP: pull capacity + thermal + size, match to GlobTherm --------------
file = matopen(MATF); allStat = read(file, "allStat"); close(file)
ok(d, k) = haskey(d, k) && d[k] isa Number && isfinite(d[k])

n = 0
open(OUT, "w") do io
    println(io, "species,k_M,E_m,lambda_max,g,T_A,T_typical,Ww_i,CTmax,CTmin,lat,class")
    for sp in sort(collect(keys(allStat)))
        haskey(gt, sp) || continue
        d = allStat[sp]
        (ok(d, "k_M") && ok(d, "v") && ok(d, "L_m") && ok(d, "E_m") && ok(d, "g")) || continue
        lambda_max = d["v"] / d["L_m"]
        T_A  = ok(d, "T_A") ? string(d["T_A"]) : ""
        Ttyp = ok(d, "T_typical") ? string(d["T_typical"]) : ""
        wwi  = ok(d, "Ww_i") ? string(d["Ww_i"]) : ""
        cm, cn, lat, cls = gt[sp]
        println(io, sp, ",", d["k_M"], ",", d["E_m"], ",", lambda_max, ",", d["g"], ",",
                T_A, ",", Ttyp, ",", wwi, ",",
                isnan(cm) ? "" : cm, ",", isnan(cn) ? "" : cn, ",", isnan(lat) ? "" : lat, ",", cls)
        global n += 1
    end
end
println("AmP-GlobTherm matched: ", n, " -> ", OUT)
