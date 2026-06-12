# ===========================================================================
# Side-extraction of per-species DEB *process rates* from AmP allStat.mat, for
# the per-axis validation (Idea B). Surfaces the reproduction-/growth-side rates
# that the main library JSON does not carry, WITHOUT regenerating/altering the
# committed data/AmP_Species_Library.json (the {p_Am,p_M,kap,v} -> A0/axes/KA
# mapping is untouched). Read-only on allStat.mat.
#
# STANDALONE — needs MAT (+ Statistics), NOT part of the project Manifest:
#   julia +release -e 'using Pkg; Pkg.activate(mktempdir()); Pkg.add("MAT")'
#   julia +release --project=<that-env> scripts/extract_amp_reproduction_rates.jl
#
# Writes data/external/amp_reproduction_rates.csv (committed):
#   species,k_M,R_i,r_B,kap_R,k_J,Ww_i
#     k_M  somatic maintenance rate constant   (1/d)  — maintenance axis
#     R_i  ultimate reproduction rate          (#/d)  — reproduction axis
#     r_B  von Bertalanffy growth rate         (1/d)  — growth axis
#     kap_R reproduction efficiency            (-)    — intensive reproduction
#     k_J  maturity maintenance rate constant  (1/d)
#     Ww_i ultimate wet weight                 (g)    — body-size control
# All at AmP reference temperature. See docs/notes/comadre_peraxis_prereg.md.
# ===========================================================================

using MAT, Statistics

const MATF = joinpath(@__DIR__, "..", "data", "allStat.mat")
const OUT  = joinpath(@__DIR__, "..", "data", "external", "amp_reproduction_rates.csv")

file = matopen(MATF); allStat = read(file, "allStat"); close(file)

okpos(d, k) = haskey(d, k) && d[k] isa Number && isfinite(d[k]) && d[k] > 0
okany(d, k) = haskey(d, k) && d[k] isa Number && isfinite(d[k])   # kap_R may be in (0,1]

n = 0
open(OUT, "w") do io
    println(io, "species,k_M,R_i,r_B,kap_R,k_J,Ww_i")
    for sp in sort(collect(keys(allStat)))
        d = allStat[sp]
        # require the core rates; allow Ww_i/kap_R to be blank if missing
        (okpos(d, "k_M") && okpos(d, "R_i") && okpos(d, "r_B") && okpos(d, "k_J")) || continue
        kapR = okany(d, "kap_R") ? string(d["kap_R"]) : ""
        wwi  = okpos(d, "Ww_i") ? string(d["Ww_i"]) : ""
        println(io, sp, ",", d["k_M"], ",", d["R_i"], ",", d["r_B"], ",", kapR, ",", d["k_J"], ",", wwi)
        global n += 1
    end
end
println("wrote $n species -> ", OUT)
