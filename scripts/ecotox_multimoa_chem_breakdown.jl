# ecotox_multimoa_chem_breakdown.jl -- per-chemical AmP coverage, to see WHICH chemicals carry
# each axis (matters because the high-coverage chemicals -- AChE pesticides, DDT, atrazine, metals
# -- have the most CONTESTED pMoA, while the clean-pMoA chemicals are sparse).
# Run after the extractor + coverage scan.  julia +release --project=. scripts/ecotox_multimoa_chem_breakdown.jl
import JSON

const POTENCY = Set(["LC50", "EC50", "IC50"])
lib = JSON.parsefile("data/AmP_Species_Library.json")
amp = Set{String}()
for (sp, v) in lib
    (v isa AbstractDict && haskey(v, "alpha_axes") && length(v["alpha_axes"]) == 4) || continue
    k = float(v["alpha_axes"][3]); (isfinite(k) && 0 < k < 1) && push!(amp, sp)
end
genus_species(l) = (t = split(l, "_"); length(t) >= 2 ? string(t[1], "_", t[2]) : String(l))

name = Dict{String,String}(); axis = Dict{String,String}()
for ln in eachline("data/ecotox_multimoa_panel.csv")
    f = split(ln, ","); (length(f) < 4 || f[1] == "cas") && continue
    name[f[1]] = f[2]; axis[f[1]] = f[4]
end

sp_of = Dict{String,Set{String}}()   # cas -> {amp species}
for ln in eachline("data/external/ecotox_multimoa_extract.csv")
    startswith(ln, "cas,") && continue
    f = split(ln, ","); length(f) >= 11 || continue
    String(f[7]) in POTENCY || continue
    key = genus_species(String(f[5])); (key in amp) || continue
    push!(get!(sp_of, String(f[1]), Set{String}()), key)
end

for ax in ["maintenance", "reproduction", "growth", "assimilation"]
    println("\n=== $ax : distinct AmP species per chemical ===")
    rows = [(get(name, c, c), length(s)) for (c, s) in sp_of if get(axis, c, "") == ax]
    for (nm, n) in sort(rows; by = x -> -x[2])
        println("  ", rpad(nm, 28), n)
    end
end
