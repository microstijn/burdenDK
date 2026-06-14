# extract_amp_size_proxies.jl -- physical body-size proxies per AmP species, as a SIDECAR to the
# committed AmP_Species_Library.json (which predates these fields; a deliberate key-sorted regen
# would fold them in via the AmP_Translator.jl add). Reads data/allStat.mat directly.
#
# Provides the cross-species BODY-SIZE CONTROL for the across-axis capacity-weighting test. The
# n=310 negative control showed acute LC50 is body-size dominated; structural L_m is a poor
# cross-species size proxy (delta_M varies). Ultimate WET WEIGHT Ww_i [g] is mass-based and the
# standard acute-LC50 size covariate. Keyed by allStat species name (== AmP library key).
#
# Run (needs MAT; ~20 s to load the 23 MB .mat):
#   julia +release --project=. scripts/extract_amp_size_proxies.jl
using MAT

const OUT = "data/external/amp_size_proxies.csv"
_posf(x) = (x isa Float64 && isfinite(x) && x > 0) ? x : NaN

f = matopen("data/allStat.mat")
allStat = read(f, "allStat"); close(f)

n = Ref(0)
open(OUT, "w") do io
    println(io, "species,Ww_i_g,Wd_i_g,L_i_cm")
    for sp in sort(collect(keys(allStat)))
        rec = allStat[sp]
        rec isa AbstractDict || continue
        Ww = _posf(get(rec, "Ww_i", NaN))
        Wd = _posf(get(rec, "Wd_i", NaN))
        Li = _posf(get(rec, "L_i", NaN))
        isnan(Ww) && isnan(Wd) && isnan(Li) && continue
        println(io, join([sp,
            isnan(Ww) ? "" : round(Ww, sigdigits = 6),
            isnan(Wd) ? "" : round(Wd, sigdigits = 6),
            isnan(Li) ? "" : round(Li, sigdigits = 6)], ","))
        n[] += 1
    end
end
println("wrote $OUT  ($(n[]) species with >=1 size proxy)")
