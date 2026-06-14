# ===========================================================================================
# ENDPOINT-SENSITIVITY diagnostic for the across-axis weighting test.
# "Not all endpoints are equal." The committed acute test pooled endpoint in {LC50,EC50,IC50}
# across ALL effect groups. This re-runs the within-species kappa-contrast Delta under several
# endpoint POLICIES to see whether the conclusion depends on which endpoints we admit:
#   (a) LC50 | MOR        -- clean apical acute mortality (the gold-standard acute endpoint)
#   (b) EC50 | IMM        -- clean acute immobilization
#   (c) apical acute      -- (a)+(b)
#   (d) pooled acute      -- LC50/EC50/IC50 any effect  (the COMMITTED acute test)
#
#   julia +release --project=. examples/across_axis_weighting_endpoint_sensitivity.jl
# ===========================================================================================
import JSON
using Statistics, Printf, Random

const EXTRACT = "data/external/ecotox_multimoa_extract.csv"
const PANEL   = "data/ecotox_multimoa_panel.csv"
const MIN_SPC = 3
Random.seed!(20260614)
ordinalrank(v) = (p = sortperm(v); r = similar(p); r[p] = 1:length(v); Float64.(r))
spearman(x, y) = cor(ordinalrank(x), ordinalrank(y))
function perm_p(x, y; nperm = 50000)
    obs = abs(spearman(x, y)); c = 0
    for _ in 1:nperm; abs(spearman(x, Random.shuffle(y))) >= obs - 1e-12 && (c += 1); end
    return (c + 1) / (nperm + 1)
end

ax_of = Dict{String,String}(); tier_of = Dict{String,String}()
for ln in eachline(PANEL)
    f = split(ln, ","); (length(f) < 5 || f[1] == "cas") && continue
    ax_of[f[1]] = f[4]; tier_of[f[1]] = f[5]
end
lib = JSON.parsefile("data/AmP_Species_Library.json")
kappa_of = Dict{String,Float64}()
for (sp, v) in lib
    (v isa AbstractDict && haskey(v, "alpha_axes") && length(v["alpha_axes"]) == 4) || continue
    k = float(v["alpha_axes"][3]); (isfinite(k) && 0 < k < 1) && (kappa_of[sp] = k)
end
genus_species(l) = (t = split(l, "_"); length(t) >= 2 ? string(t[1], "_", t[2]) : String(l))

rows = NamedTuple[]
for ln in eachline(EXTRACT)
    startswith(ln, "cas,") && continue
    f = split(ln, ","); length(f) >= 11 || continue
    cas = String(f[1]); haskey(ax_of, cas) || continue
    ax_of[cas] in ("maintenance", "reproduction") || continue
    sp = genus_species(String(f[5])); haskey(kappa_of, sp) || continue
    v = tryparse(Float64, f[11]); (v === nothing || v <= 0) && continue
    push!(rows, (sp = sp, cas = cas, axis = ax_of[cas], ep = String(f[7]), ef = String(f[8]), val = v))
end

function delta_rho(keep, allowed)
    acc = Dict{Tuple{String,String},Vector{Float64}}()
    for r in rows
        (tier_of[r.cas] in allowed && keep(r)) || continue
        push!(get!(acc, (r.sp, r.cas), Float64[]), log10(r.val))
    end
    cell = Dict(k => median(v) for (k, v) in acc)
    chem_sp = Dict{String,Vector{String}}()
    for ((sp, cas), _) in cell; push!(get!(chem_sp, cas, String[]), sp); end
    chems = [c for (c, s) in chem_sp if length(unique(s)) >= MIN_SPC]
    yc = Dict{Tuple{String,String},Float64}()
    for c in chems
        sps = unique(chem_sp[c]); m = mean(cell[(s, c)] for s in sps)
        for s in sps; yc[(s, c)] = cell[(s, c)] - m; end
    end
    spM = Dict{String,Vector{Float64}}(); spR = Dict{String,Vector{Float64}}()
    for ((s, c), yv) in yc; push!(get!(ax_of[c] == "maintenance" ? spM : spR, s, Float64[]), yv); end
    k = Float64[]; d = Float64[]
    for s in keys(spM)
        haskey(spR, s) || continue
        push!(k, kappa_of[s]); push!(d, mean(spR[s]) - mean(spM[s]))
    end
    length(k) < 6 && return (n = length(k), rho = NaN, p = NaN)
    return (n = length(k), rho = spearman(k, d), p = perm_p(k, d))
end

policies = [
    ("a) LC50|MOR  (apical mortality)", r -> r.ep == "LC50" && r.ef == "MOR"),
    ("b) EC50|IMM  (immobilization)",   r -> r.ep == "EC50" && r.ef == "IMM"),
    ("c) apical acute (a+b)",           r -> (r.ep == "LC50" && r.ef == "MOR") || (r.ep == "EC50" && r.ef == "IMM")),
    ("d) pooled acute LC50/EC50/IC50 [COMMITTED]", r -> r.ep in ("LC50", "EC50", "IC50")),
]
println("ENDPOINT-SENSITIVITY of rho(kappa, Delta)   [model predicts rho>0]")
for (lab, tiers) in (("CORE", Set(["core"])), ("CORE+STRATUM", Set(["core", "stratum"])))
    println("\n=== $lab ===")
    @printf("  %-46s %5s %8s %8s\n", "endpoint policy", "n", "rho", "perm p")
    for (name, keep) in policies
        res = delta_rho(keep, tiers)
        if isnan(res.rho)
            @printf("  %-46s %5d   (too thin)\n", name, res.n)
        else
            @printf("  %-46s %5d %+8.3f %8.4f\n", name, res.n, res.rho, res.p)
        end
    end
end
