# ecotox_apical_ecx_seed.jl -- the ECOTOX-derived SEED of the curated apical-EC50 database.
# Authoritative, local, reproducible: keeps ONLY proper apical EC50 (no NOEC/LOEC), endpoint matched
# to axis (reproduction -> REP; maintenance -> GRO/PHY), for the panel's M/R chemicals, intersected
# with AmP. Shows exactly what the local dump yields so the external (OECD/SSD) gap is explicit.
#   julia +release --project=. scripts/ecotox_apical_ecx_seed.jl
import JSON
using Statistics, Printf

const EXTRACT = "data/external/ecotox_multimoa_extract.csv"
const PANEL   = "data/ecotox_multimoa_panel.csv"
const OUT     = "data/external/curated_apical_ecx_ecotox_seed.csv"
const MATCH   = Dict("maintenance" => Set(["GRO", "PHY"]), "reproduction" => Set(["REP"]))

ax_of = Dict{String,String}(); tier_of = Dict{String,String}(); name_of = Dict{String,String}()
for ln in eachline(PANEL)
    f = split(ln, ","); (length(f) < 5 || f[1] == "cas") && continue
    name_of[f[1]] = f[2]; ax_of[f[1]] = f[4]; tier_of[f[1]] = f[5]
end
lib = JSON.parsefile("data/AmP_Species_Library.json")
kappa_of = Dict{String,Float64}()
for (sp, v) in lib
    (v isa AbstractDict && haskey(v, "alpha_axes") && length(v["alpha_axes"]) == 4) || continue
    k = float(v["alpha_axes"][3]); (isfinite(k) && 0 < k < 1) && (kappa_of[sp] = k)
end
gs(l) = (t = split(l, "_"); length(t) >= 2 ? string(t[1], "_", t[2]) : String(l))

acc = Dict{Tuple{String,String},Vector{Float64}}()   # (species, cas) -> [log10 EC50]
for ln in eachline(EXTRACT)
    startswith(ln, "cas,") && continue
    f = split(ln, ","); length(f) >= 11 || continue
    f[7] == "EC50" || continue                         # apical POTENCY only -- no NOEC/LOEC
    cas = String(f[1]); haskey(ax_of, cas) || continue
    ax = ax_of[cas]; haskey(MATCH, ax) || continue
    String(f[8]) in MATCH[ax] || continue              # endpoint matched to axis
    sp = gs(String(f[5])); haskey(kappa_of, sp) || continue
    v = tryparse(Float64, f[11]); (v === nothing || v <= 0) && continue
    push!(get!(acc, (sp, cas), Float64[]), log10(v))
end

rows = [(sp = k[1], cas = k[2], chem = name_of[k[2]], axis = ax_of[k[2]], tier = tier_of[k[2]],
         kappa = kappa_of[k[1]], ec50 = 10^median(v), n = length(v)) for (k, v) in acc]
open(OUT, "w") do io
    println(io, "species,kappa,chemical,cas,axis,tier,EC50_ugL,n_records,endpoint")
    for r in sort(rows; by = x -> (x.axis, x.chem, x.sp))
        println(io, join([r.sp, round(r.kappa, digits=3), r.chem, r.cas, r.axis, r.tier,
            round(r.ec50, sigdigits=4), r.n, r.axis == "reproduction" ? "REP" : "GRO/PHY"], ","))
    end
end

spM = Set(r.sp for r in rows if r.axis == "maintenance")
spR = Set(r.sp for r in rows if r.axis == "reproduction")
mr  = sort(collect(intersect(spM, spR)); by = s -> -kappa_of[s])
println("apical EC50 cells: ", length(rows), "   maintenance sp: ", length(spM),
        "   reproduction sp: ", length(spR), "   M∩R: ", length(mr))
println("\nper-chemical apical-EC50 species coverage:")
bychem = Dict{String,Set{String}}()
for r in rows; push!(get!(bychem, "$(r.axis):$(r.chem)", Set{String}()), r.sp); end
for (c, s) in sort(collect(bychem); by = x -> -length(x[2])); @printf("  %-32s %d sp\n", c, length(s)); end
println("\nM∩R species (the only ones usable for the apical-EC50 weighting test):")
for s in mr; @printf("  %-26s kappa=%.3f\n", s, kappa_of[s]); end
println("\nwrote $OUT")
