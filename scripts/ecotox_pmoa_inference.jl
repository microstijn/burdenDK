# ecotox_pmoa_inference.jl — DIAGNOSTIC (negative result, kept on record).
#
# Aggregates the records from scripts/ecotox_pmoa_inference.awk into per-compound, per-axis
# median effect concentrations, and prints the "most-sensitive sublethal axis". The point of
# keeping this is the NEGATIVE finding it demonstrates:
#
#   "most-sensitive endpoint" does NOT equal "pMoA". The DEB kappa-rule cascade routes any
#   energetic stress into reduced reproduction first (reproduction is funded by the leftover
#   after maintenance + growth), so reproduction "wins" for nearly every compound regardless
#   of the process actually taxed. Effect-code sensitivity is also taxon-confounded. Only an
#   endocrine, reproduction-DISPROPORTIONATE stressor (e.g. bisphenol-A) gives a clean signal.
#
# => Do not route stressors from raw ECOTOX endpoint sensitivities. Use fitted DEBtox pMoAs
#    (docs/notes/pmoa_evidence_to_gather.md, Tier 1).
#
# Run (from repo root, after the awk step has written the records file):
#   julia +release --project=. scripts/ecotox_pmoa_inference.jl
using Statistics
const REC = "data/external/ecotox_pmoa_records.csv"
const OUT = "data/external/ecotox_pmoa_inference_results.txt"
isfile(REC) || error("$REC not found — run scripts/ecotox_pmoa_inference.awk first (see its header).")
name = Dict("7440439"=>"Cadmium","7440508"=>"Copper","2921882"=>"Chlorpyrifos","1912249"=>"Atrazine","80057"=>"Bisphenol A")
prior = Dict("7440439"=>"maintenance(table)","7440508"=>"maintenance(table)","2921882"=>"maintenance(table)","1912249"=>"not in table","80057"=>"reproduction(expected)")
axes = ["assimilation","maintenance","growth","reproduction","survival"]
d = Dict{Tuple{String,String,String},Vector{Float64}}()
for ln in eachline(REC)
    f = split(ln, ","); length(f) == 4 || continue
    push!(get!(d, (String(f[1]),String(f[2]),String(f[3])), Float64[]), parse(Float64, f[4]))
end
open(OUT, "w") do io
    for sink in (stdout, io)
        println(sink, "ECOTOX effect-code pMoA inference — DIAGNOSTIC (negative result).")
        println(sink, "most-sensitive sublethal endpoint != pMoA: the DEB kappa-cascade makes")
        println(sink, "reproduction win for almost any stress. See script headers.\n")
        for cas in ["7440439","7440508","2921882","1912249","80057"]
            println(sink, "===== ", name[cas], "   [a priori: ", prior[cas], "] =====")
            println(sink, rpad("axis",14), rpad("n_chron",8), rpad("med_chronic_ugL",16), rpad("n_acute",8), "med_acute_ugL")
            best = ("", Inf)
            for ax in axes
                ch = get(d,(cas,ax,"chronic"),Float64[]); ac = get(d,(cas,ax,"acute"),Float64[])
                mc = isempty(ch) ? NaN : median(ch); ma = isempty(ac) ? NaN : median(ac)
                println(sink, rpad(ax,14), rpad(length(ch),8), rpad(isnan(mc) ? "-" : string(round(mc,sigdigits=3)),16),
                        rpad(length(ac),8), isnan(ma) ? "-" : string(round(ma,sigdigits=3)))
                if ax != "survival" && !isnan(mc) && mc < best[2]; best = (ax, mc); end
            end
            println(sink, "  --> most-sensitive sublethal axis (chronic) = ", best[1]=="" ? "n/a" : uppercase(best[1]),
                    "   (confounded — see header)\n")
        end
    end
end
println("\nWrote ", OUT)
