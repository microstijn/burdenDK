# ===========================================================================================
# ACROSS-AXIS CAPACITY-WEIGHTING TEST -- the central open validation question.
#
# Does the model's per-species kappa-driven axis weighting predict species x MoA sensitivity
# BEYOND body size? Operative weights w = (1/2, kappa/4, kappa/4, (1-kappa)/2) over
# (assim, maint, growth, repro), kappa = alpha_axes[3]. Assimilation is fixed (no cross-species
# signal); the discriminating contrast is MAINTENANCE (w_M = kappa/4) vs REPRODUCTION
# (w_R = (1-kappa)/2). Ratio w_M/w_R = kappa/(2(1-kappa)) increases monotonically in kappa.
#
# DESIGN (uncalibrated, size-controlled by construction):
#   1. Per (species, chemical): median log10 LC50/EC50/IC50 (ug/L), from the local ECOTOX dump.
#   2. CHEMICAL-CENTER: y~_{s,c} = y_{s,c} - mean_s y_{.,c}  (removes the chemical main effect;
#      y~ high = species s is relatively MORE RESISTANT to chemical c).
#   3. Per species: resid_M = mean y~ over maintenance chemicals; resid_R = mean over reproduction.
#      Delta(s) = resid_R - resid_M. Differencing removes the species main effect AND body size.
#      Delta > 0  <=>  relatively more SENSITIVE to maintenance than reproduction.
#   4. PREDICTION: rho(kappa, Delta) > 0. Test: permutation Spearman, size partial, OLS slope.
#
# Panel tiers (data/ecotox_multimoa_panel.csv): core = defensible pMoA; stratum = AChE/DDT/PAH
# (contested pMoA) reported as a robustness arm. Inputs produced by
# scripts/extract_ecotox_multimoa.awk + scripts/extract_amp_size_proxies.jl.
#
#   julia +release --project=. examples/across_axis_weighting_capacity_test.jl
# ===========================================================================================
import JSON
using Statistics, Printf, Random, LinearAlgebra

const EXTRACT  = "data/external/ecotox_multimoa_extract.csv"
const PANEL    = "data/ecotox_multimoa_panel.csv"
const SIZES    = "data/external/amp_size_proxies.csv"
const POTENCY  = Set(["LC50", "EC50", "IC50"])
const MIN_SPC  = 3          # min AmP species per chemical for stable centering
const AXES2    = ("maintenance", "reproduction")
Random.seed!(20260614)

# ---- rank helpers ----
ordinalrank(v) = (p = sortperm(v); r = similar(p); r[p] = 1:length(v); Float64.(r))
spearman(x, y) = cor(ordinalrank(x), ordinalrank(y))
function partial_spearman(x, y, z)             # Spearman of x,y controlling z (rank-linear)
    rx, ry, rz = ordinalrank(x), ordinalrank(y), ordinalrank(z)
    Z = hcat(ones(length(rz)), rz); resid(v) = v .- Z * (Z \ v)
    return cor(resid(rx), resid(ry))
end
function perm_p(x, y, stat; nperm = 50000)     # two-sided permutation p for a correlation stat
    obs = abs(stat(x, y)); c = 0
    for _ in 1:nperm
        abs(stat(x, Random.shuffle(y))) >= obs - 1e-12 && (c += 1)
    end
    return (c + 1) / (nperm + 1)
end
# OLS: returns (coefs, tstats, df) for y ~ [1 X...]
function ols(y, Xcols...)
    X = hcat(ones(length(y)), Xcols...)
    n, p = size(X)
    b = X \ y
    resid = y .- X * b
    s2 = sum(abs2, resid) / (n - p)
    XtXi = inv(X' * X)
    se = sqrt.(max.(diag(XtXi) .* s2, 0.0))
    return b, b ./ se, n - p
end

# ---- panel: cas -> (axis, tier, name) ----
ax_of = Dict{String,String}(); tier_of = Dict{String,String}(); name_of = Dict{String,String}()
for ln in eachline(PANEL)
    f = split(ln, ","); (length(f) < 5 || f[1] == "cas") && continue
    name_of[f[1]] = f[2]; ax_of[f[1]] = f[4]; tier_of[f[1]] = f[5]
end

# ---- AmP: kappa + size ----
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

# ---- cells: (species, cas) -> [log10 potency] (maintenance/reproduction, AmP-resident) ----
vals = Dict{Tuple{String,String},Vector{Float64}}()
for ln in eachline(EXTRACT)
    startswith(ln, "cas,") && continue
    f = split(ln, ","); length(f) >= 11 || continue
    String(f[7]) in POTENCY || continue
    cas = String(f[1]); haskey(ax_of, cas) || continue
    ax_of[cas] in AXES2 || continue
    sp = genus_species(String(f[5])); haskey(kappa_of, sp) || continue
    v = tryparse(Float64, f[11]); (v === nothing || v <= 0) && continue
    push!(get!(vals, (sp, cas), Float64[]), log10(v))
end
cell = Dict((k => median(v)) for (k, v) in vals)   # (species,cas) -> median log10 potency

# ---- run the test for a given tier set ----
function run_test(allowed::Set{String}, label::String)
    # chemicals in scope, with >= MIN_SPC species (for stable centering)
    chem_sp = Dict{String,Vector{String}}()
    for ((sp, cas), _) in cell
        tier_of[cas] in allowed || continue
        push!(get!(chem_sp, cas, String[]), sp)
    end
    chems = [c for (c, sps) in chem_sp if length(unique(sps)) >= MIN_SPC]
    # chemical-center
    yc = Dict{Tuple{String,String},Float64}()
    for c in chems
        sps = unique(chem_sp[c]); m = mean(cell[(s, c)] for s in sps)
        for s in sps; yc[(s, c)] = cell[(s, c)] - m; end
    end
    # per-species Delta = mean_R(yc) - mean_M(yc)
    sp_M = Dict{String,Vector{Float64}}(); sp_R = Dict{String,Vector{Float64}}()
    for ((s, c), yv) in yc
        tgt = ax_of[c] == "maintenance" ? sp_M : sp_R
        push!(get!(tgt, s, Float64[]), yv)
    end
    recs = NamedTuple[]
    for s in keys(sp_M)
        haskey(sp_R, s) || continue
        push!(recs, (sp = s, kappa = kappa_of[s], Delta = mean(sp_R[s]) - mean(sp_M[s]),
                     logWw = haskey(ww_of, s) ? log10(ww_of[s]) : NaN,
                     nM = length(sp_M[s]), nR = length(sp_R[s])))
    end
    sort!(recs; by = r -> -r.kappa)

    println("\n", "="^78)
    println("[$label]   tiers=", join(sort(collect(allowed)), "+"),
            "   chemicals used=", length(chems), "   species (M&R)=", length(recs))
    k = [r.kappa for r in recs]; d = [r.Delta for r in recs]
    rho = spearman(k, d); p = perm_p(k, d, spearman)
    @printf("  PRIMARY   rho(kappa, Delta) = %+.3f   (perm p = %.4f, n=%d)\n", rho, p, length(recs))
    @printf("            [prediction: rho > 0  =>  higher kappa = more maintenance-sensitive]\n")
    # size control
    has = [r for r in recs if isfinite(r.logWw)]
    if length(has) >= 5
        kk = [r.kappa for r in has]; dd = [r.Delta for r in has]; ww = [r.logWw for r in has]
        prho = partial_spearman(kk, dd, ww)
        @printf("  SIZE      partial rho(kappa, Delta | log Ww) = %+.3f   (n=%d with Ww)\n", prho, length(has))
        @printf("            collinearity rho(kappa, log Ww) = %+.2f\n", spearman(kk, ww))
        b, t, df = ols(dd, kk, ww)
        @printf("  OLS       Delta ~ kappa + logWw :  slope_kappa=%+.3f  t=%+.2f (df=%d)\n", b[2], t[2], df)
    end
    b, t, df = ols(d, k)
    @printf("  OLS       Delta ~ kappa          :  slope=%+.3f  t=%+.2f (df=%d)\n", b[2], t[2], df)
    return recs
end

println("ACROSS-AXIS CAPACITY-WEIGHTING TEST (maintenance vs reproduction)")
println("cells (species x chemical) in scope: ", length(cell))
core_recs = run_test(Set(["core"]), "CORE (defensible pMoA)")
run_test(Set(["core", "stratum"]), "CORE + STRATUM (incl. AChE/DDT/PAH)")

# ---- per-species core table + save ----
println("\n=== CORE per-species (sorted by kappa) ===")
@printf("  %-30s %6s %8s %9s %s\n", "species", "kappa", "Delta", "log10Ww", "nM/nR")
for r in core_recs
    @printf("  %-30s %6.3f %+8.3f %9.3g  %d/%d\n", r.sp, r.kappa, r.Delta, r.logWw, r.nM, r.nR)
end
open("data/external/across_axis_weighting_core_results.csv", "w") do io
    println(io, "species,kappa,Delta,log10_Ww,n_maint,n_repro")
    for r in core_recs
        println(io, join([r.sp, round(r.kappa, digits = 4), round(r.Delta, digits = 4),
            isfinite(r.logWw) ? round(r.logWw, digits = 4) : "", r.nM, r.nR], ","))
    end
end
println("\nwrote data/external/across_axis_weighting_core_results.csv")
println("\nReading: Delta>0 = relatively more maintenance-sensitive. A POSITIVE rho(kappa,Delta)")
println("that SURVIVES the size partial corroborates the across-axis weighting beyond body size.")
println("A null (rho~0) is the honest result the n=310/n=5 work anticipated -- report either way.")
