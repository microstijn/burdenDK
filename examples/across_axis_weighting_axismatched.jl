# ===========================================================================================
# ACROSS-AXIS WEIGHTING -- AXIS-MATCHED SUBLETHAL endpoints (refinement of the acute-pooled test).
#
# The acute test (examples/across_axis_weighting_capacity_test.jl) pooled LC50/EC50 -- acute
# lethality does NOT exercise an EDC's reproductive pMoA (it is narcosis). This version matches the
# observed ENDPOINT to the chemical's DEB axis, so each chemical actually engages its pMoA:
#   maintenance chemical  ->  growth / physiology(respiration) endpoint  (somatic energy drain)
#   reproduction chemical ->  reproduction endpoint
# then re-runs the within-species kappa-contrast  Delta = resid_R - resid_M ;  prediction rho>0.
#
# Reuses data/external/ecotox_multimoa_extract.csv (already carries the `effect` column + NOEC/LOEC).
# FEASIBILITY FIRST: sublethal axis-matched data is sparse; the script prints effect-code coverage
# per axis before attempting the test, and reports honestly if coverage is too thin.
#
#   julia +release --project=. examples/across_axis_weighting_axismatched.jl
# ===========================================================================================
import JSON
using Statistics, Printf, Random, LinearAlgebra

const EXTRACT = "data/external/ecotox_multimoa_extract.csv"
const PANEL   = "data/ecotox_multimoa_panel.csv"
const SIZES   = "data/external/amp_size_proxies.csv"
const MIN_SPC = 3
Random.seed!(20260614)

# axis -> matched effect codes (ECOTOX `effect` group). Reproduction = REP; maintenance is a
# somatic energy drain -> growth (GRO) + physiology/respiration (PHY). MPH/DVP folded into growth.
const MATCH = Dict("maintenance" => Set(["GRO", "PHY", "DVP", "MPH"]),
                   "reproduction" => Set(["REP"]))
const SUBLETHAL = Set(["EC50", "NOEC", "LOEC", "IC50"])   # exclude LC50 (lethal) for the matched test

ordinalrank(v) = (p = sortperm(v); r = similar(p); r[p] = 1:length(v); Float64.(r))
spearman(x, y) = cor(ordinalrank(x), ordinalrank(y))
function partial_spearman(x, y, z)
    rx, ry, rz = ordinalrank(x), ordinalrank(y), ordinalrank(z)
    Z = hcat(ones(length(rz)), rz); resid(v) = v .- Z * (Z \ v)
    return cor(resid(rx), resid(ry))
end
function perm_p(x, y, stat; nperm = 50000)
    obs = abs(stat(x, y)); c = 0
    for _ in 1:nperm; abs(stat(x, Random.shuffle(y))) >= obs - 1e-12 && (c += 1); end
    return (c + 1) / (nperm + 1)
end

# panel + AmP
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
ww_of = Dict{String,Float64}()
for ln in eachline(SIZES)
    f = split(ln, ","); (length(f) < 2 || f[1] == "species") && continue
    w = tryparse(Float64, f[2]); (w !== nothing && w > 0) && (ww_of[f[1]] = w)
end
genus_species(l) = (t = split(l, "_"); length(t) >= 2 ? string(t[1], "_", t[2]) : String(l))

# read extract: keep (species, cas, axis, tier, endpoint, effect, value)
rows = NamedTuple[]
for ln in eachline(EXTRACT)
    startswith(ln, "cas,") && continue
    f = split(ln, ","); length(f) >= 11 || continue
    cas = String(f[1]); haskey(ax_of, cas) || continue
    ax = ax_of[cas]; (ax in ("maintenance", "reproduction")) || continue
    sp = genus_species(String(f[5])); haskey(kappa_of, sp) || continue
    v = tryparse(Float64, f[11]); (v === nothing || v <= 0) && continue
    push!(rows, (sp = sp, cas = cas, axis = ax, ep = String(f[7]), ef = String(f[8]), val = v))
end

# --- FEASIBILITY: effect-code coverage per axis (AmP-resident, sublethal endpoints) ---
println("=== sublethal effect-code coverage (distinct AmP species), per axis ===")
for ax in ("maintenance", "reproduction")
    println("  [$ax]")
    bytab = Dict{Tuple{String,String},Set{String}}()   # (effect, endpoint) -> species
    for r in rows
        (r.axis == ax && r.ep in SUBLETHAL) || continue
        push!(get!(bytab, (r.ef, r.ep), Set{String}()), r.sp)
    end
    for ((ef, ep), sps) in sort(collect(bytab); by = x -> -length(x[2]))
        matched = ef in MATCH[ax] ? "  <-- axis-matched" : ""
        @printf("    effect=%-5s endpoint=%-5s  %3d sp%s\n", ef, ep, length(sps), matched)
    end
end

# --- build axis-matched cells: (species, cas) -> median log10(value) over matched effect+sublethal ---
function matched_cell(allowed_tiers)
    cell = Dict{Tuple{String,String},Float64}(); acc = Dict{Tuple{String,String},Vector{Float64}}()
    for r in rows
        tier_of[r.cas] in allowed_tiers || continue
        (r.ep in SUBLETHAL && r.ef in MATCH[r.axis]) || continue
        push!(get!(acc, (r.sp, r.cas), Float64[]), log10(r.val))
    end
    for (k, v) in acc; cell[k] = median(v); end
    return cell
end

function run_matched(allowed, label)
    cell = matched_cell(allowed)
    chem_sp = Dict{String,Vector{String}}()
    for ((sp, cas), _) in cell; push!(get!(chem_sp, cas, String[]), sp); end
    chems = [c for (c, s) in chem_sp if length(unique(s)) >= MIN_SPC]
    yc = Dict{Tuple{String,String},Float64}()
    for c in chems
        sps = unique(chem_sp[c]); m = mean(cell[(s, c)] for s in sps)
        for s in sps; yc[(s, c)] = cell[(s, c)] - m; end
    end
    sp_M = Dict{String,Vector{Float64}}(); sp_R = Dict{String,Vector{Float64}}()
    for ((s, c), yv) in yc
        push!(get!(ax_of[c] == "maintenance" ? sp_M : sp_R, s, Float64[]), yv)
    end
    recs = NamedTuple[]
    for s in keys(sp_M)
        haskey(sp_R, s) || continue
        push!(recs, (sp = s, kappa = kappa_of[s], Delta = mean(sp_R[s]) - mean(sp_M[s]),
                     logWw = haskey(ww_of, s) ? log10(ww_of[s]) : NaN))
    end
    sort!(recs; by = r -> -r.kappa)
    println("\n", "="^78)
    println("[$label]  chemicals(>=3sp)=", length(chems),
            "   maint-cells=", count(c -> ax_of[c[2]] == "maintenance", keys(yc)),
            "   repro-cells=", count(c -> ax_of[c[2]] == "reproduction", keys(yc)),
            "   species (M&R) = ", length(recs))
    if length(recs) < 6
        println("  -> TOO THIN for a powered test (n=", length(recs), "). Honest verdict: ECOTOX",
                " axis-matched sublethal coverage cannot support this test.")
        for r in recs; @printf("     %-26s kappa=%.3f  Delta=%+.3f\n", r.sp, r.kappa, r.Delta); end
        return recs
    end
    k = [r.kappa for r in recs]; d = [r.Delta for r in recs]
    rho = spearman(k, d); p = perm_p(k, d, spearman)
    @printf("  PRIMARY  rho(kappa, Delta) = %+.3f   (perm p=%.4f, n=%d)   [predict rho>0]\n", rho, p, length(recs))
    has = [r for r in recs if isfinite(r.logWw)]
    if length(has) >= 5
        prho = partial_spearman([r.kappa for r in has], [r.Delta for r in has], [r.logWw for r in has])
        @printf("  SIZE     partial rho(kappa, Delta | log Ww) = %+.3f  (n=%d)\n", prho, length(has))
    end
    for r in recs; @printf("     %-26s kappa=%.3f  Delta=%+.3f  logWw=%6.2f\n", r.sp, r.kappa, r.Delta, r.logWw); end
    return recs
end

println("\n", "#"^78, "\nAXIS-MATCHED SUBLETHAL TEST")
run_matched(Set(["core"]), "CORE (defensible pMoA)")
run_matched(Set(["core", "stratum"]), "CORE + STRATUM")
