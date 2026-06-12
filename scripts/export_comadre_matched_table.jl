# ===========================================================================
# Export the AmP-COMADRE matched table for the phylogenetic (PGLS) validation
# (Idea A). Reuses the same matching as examples/comadre_partial_validation.jl
# (harmonised name map + pseudoreplication collapse) and writes one row per
# matched AmP species with the model quantities, the COMADRE recovery/pace
# quantities, and taxonomy — the input table for tree fetching + PGLS.
#
#   julia +release --project=. scripts/export_comadre_matched_table.jl
#
# Writes data/external/comadre_amp_matched.csv (committed):
#   amp_key,comadre_species,k_M,lambda_A0,g,comadre_log_damping,generation_time,
#   class,order,family
# ===========================================================================

using TwoTimescaleResilience
using Statistics

const REC = joinpath(@__DIR__, "..", "data", "external", "comadre_recovery.csv")
const MAP = joinpath(@__DIR__, "..", "data", "external", "comadre_amp_namemap.csv")
const OUT = joinpath(@__DIR__, "..", "data", "external", "comadre_amp_matched.csv")

# namemap: comadre_species => amp_key (resolved only)
namemap = Dict{String, String}()
for line in Iterators.drop(eachline(MAP), 1)
    f = split(line, ","); length(f) >= 2 && f[2] != "" && (namemap[String(f[1])] = String(f[2]))
end

# COMADRE rows grouped by resolved AmP key
struct Row; logdamp::Float64; gen::Float64; class::String; order::String; family::String; sp::String; end
grouped = Dict{String, Vector{Row}}()
for line in Iterators.drop(eachline(REC), 1)
    f = split(line, ","); length(f) >= 7 || continue
    sp = String(f[1]); key = get(namemap, sp, nothing); key === nothing && continue
    gen = f[4] == "" ? NaN : parse(Float64, f[4])
    push!(get!(grouped, key, Row[]),
          Row(parse(Float64, f[3]), gen, String(f[5]), String(f[6]), String(f[7]), sp))
end

lib = load_amp_species_library()
n = 0
open(OUT, "w") do io
    println(io, "amp_key,comadre_species,k_M,lambda_A0,g,comadre_log_damping,generation_time,class,order,family")
    for (key, rows) in sort(collect(grouped); by = first)
        haskey(lib, key) || continue
        local p
        try; p = amp_species_deb_params(lib, key); catch; continue; end
        (0 < p.lambda_min <= p.lambda_max && isfinite(p.A0)) || continue
        kM = p.lambda_min
        lamA0 = restoring_force_from_margin(p.A0, p)
        g = p.lambda_max / p.lambda_min
        gens = [r.gen for r in rows if isfinite(r.gen)]
        gt = isempty(gens) ? "" : string(round(mean(gens); digits = 4))
        ld = round(mean(r.logdamp for r in rows); digits = 5)
        r1 = rows[1]
        sps = join((r.sp for r in rows), ";")   # all COMADRE names that collapsed here
        println(io, key, ",", sps, ",", kM, ",", lamA0, ",", g, ",", ld, ",", gt,
                ",", r1.class, ",", r1.order, ",", r1.family)
        global n += 1
    end
end
println("wrote $n matched species -> ", OUT)
