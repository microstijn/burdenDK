# ===========================================================================
# Resolve the positive-a_p surprise (roadmap Idea B follow-up). The per-axis
# pre-registration predicted age-at-puberty a_p -> compensation NEGATIVE (later
# puberty -> less compensation); it came out POSITIVE and survived controlling
# generation time and R_i. This script tests the alternative explanations:
#   (1) pace conflation     -> control generation time
#   (2) fecundity proxy     -> control R_i
#   (3) body-size proxy     -> control Ww_i
#   (4) MATRIX-DIMENSION artifact -> control mean matrix dimension
#       (late-puberty species have more juvenile stages -> bigger matrices;
#        if reactivity scales with dimension, a_p<->compensation could be an
#        artifact). Computed directly from the matrices here.
#
# STANDALONE — RData.jl + DataFrames throwaway env. Self-contained (joins R_i,
# a_p, Ww_i, k_M via the committed CSVs).
#   julia +release --project=<rdata-env> scripts/comadre_ap_diagnostic.jl
# ===========================================================================

using RData, DataFrames, LinearAlgebra, Statistics, Printf

const RAW = joinpath(@__DIR__, "..", "data", "external", "COMADRE_v.4.26.4.0.RData")
const MAP = joinpath(@__DIR__, "..", "data", "external", "comadre_amp_namemap.csv")
const RATES = joinpath(@__DIR__, "..", "data", "external", "amp_reproduction_rates.csv")

d = load(RAW); c = d["comadre"]; md = c["metadata"]; mats = c["mat"]
getcol(col, i) = (hasproperty(md, col) ? (v = md[i, col]; v === missing ? "" : string(v)) : "")
clean(A) = any(ismissing, A) ? nothing : Float64.(A)

# returns (log10 damping, reactivity=max col sum of A/lambda1, dimension) or NaNs
function metrics(Araw)
    A = clean(Araw); A === nothing && return (NaN, NaN, NaN)
    (size(A, 1) == size(A, 2) && size(A, 1) >= 2 && all(isfinite, A) && all(>=(0), A)) || return (NaN, NaN, NaN)
    mods = sort(abs.(eigvals(A)); rev = true)
    (mods[1] > 0 && mods[2] > 0 && mods[1] > mods[2]) || return (NaN, NaN, NaN)
    cs = vec(sum(A ./ mods[1]; dims = 1))
    return (log10(mods[1] / mods[2]), maximum(cs), Float64(size(A, 1)))
end
function gentime(Uraw, Fraw, Araw)
    U = clean(Uraw); F = clean(Fraw); A = clean(Araw)
    (U === nothing || F === nothing || A === nothing) && return NaN
    try
        k = size(U, 1); N = inv(I(k) - U)
        R0 = maximum(abs.(eigvals(F * N))); lam = maximum(abs.(eigvals(A)))
        (R0 > 0 && lam > 0 && !isapprox(lam, 1.0; atol = 1e-6)) || return NaN
        return log(R0) / log(lam)
    catch; return NaN; end
end

namemap = Dict{String, String}()
for line in Iterators.drop(eachline(MAP), 1)
    f = split(line, ","); length(f) >= 2 && f[2] != "" && (namemap[String(f[1])] = String(f[2]))
end
rates = Dict{String, NTuple{4, Float64}}()   # key => (R_i, r_B, Ww_i, a_p)
for line in Iterators.drop(eachline(RATES), 1)
    f = split(line, ","); length(f) >= 8 || continue
    pf(x) = x == "" ? NaN : parse(Float64, x)
    rates[String(f[1])] = (pf(f[3]), pf(f[4]), pf(f[7]), pf(f[8]))
end

# aggregate per species: recovery, compensation(reactivity), gen, dimension
agg = Dict{String, Vector{NTuple{4, Float64}}}()
for i in 1:nrow(md)
    sp = getcol(:SpeciesAccepted, i)
    (sp == "" || occursin(" sp", sp) || occursin("_sp", sp)) && continue
    cap = getcol(:MatrixCaptivity, i); (cap == "W" || cap == "Wild") || continue
    tr = getcol(:MatrixTreatment, i); (tr == "Unmanipulated" || tr == "") || continue
    dr, reac, dim = metrics(mats[i]["matA"]); isfinite(dr) || continue
    gt = gentime(mats[i]["matU"], mats[i]["matF"], mats[i]["matA"])
    push!(get!(agg, replace(sp, "," => " "), NTuple{4, Float64}[]),
          (dr, reac > 0 ? log10(reac) : NaN, gt, dim))
end

aP = Float64[]; Ri = Float64[]; Ww = Float64[]; comp = Float64[]; gen = Float64[]; dim = Float64[]
for (spc, vs) in agg
    key = get(namemap, spc, nothing); key === nothing && continue
    haskey(rates, key) || continue
    r = rates[key]; isfinite(r[4]) || continue   # need a_p
    cv = [v[2] for v in vs if isfinite(v[2])]; isempty(cv) && continue
    gv = [v[3] for v in vs if isfinite(v[3])]
    push!(aP, r[4]); push!(Ri, r[1]); push!(Ww, r[3])
    push!(comp, mean(cv)); push!(gen, isempty(gv) ? NaN : mean(gv))
    push!(dim, mean(v[4] for v in vs))
end

ordinalrank(v) = (p = sortperm(v); r = similar(p); r[p] = 1:length(v); Float64.(r))
spear(x, y) = cor(ordinalrank(x), ordinalrank(y))
function prank(x, y, ctrls::Vector)
    n = length(x); Z = hcat(ones(n), (ordinalrank(c) for c in ctrls)...)
    res(v) = (vr = ordinalrank(v); vr .- Z * (Z \ vr)); cor(res(x), res(y))
end
finite(xs...) = [all(isfinite, t) for t in zip(xs...)]
lg(v) = log.(v)

println("species with a_p + compensation: ", length(aP), "\n")
println("--- where does a_p sit? (raw Spearman) ---")
m = finite(aP, gen, Ri, Ww, dim, comp)
@show round(spear(aP[m], gen[m]); digits=3)   # pace loading
@show round(spear(aP[m], Ri[m]); digits=3)    # fecundity co-variation
@show round(spear(aP[m], dim[m]); digits=3)   # matrix-dimension co-variation
@show round(spear(dim[m], comp[m]); digits=3) # does dimension drive compensation?
@show round(spear(aP[m], comp[m]); digits=3)  # raw a_p -> compensation

println("\n--- a_p -> compensation, progressively controlled (partial rank) ---")
aPm, Rim, Wwm, dimm, compm, genm = aP[m], Ri[m], Ww[m], dim[m], comp[m], gen[m]
@printf("%-40s %+.3f\n", "| gen",                 prank(lg(aPm), compm, [genm]))
@printf("%-40s %+.3f\n", "| gen, R_i",            prank(lg(aPm), compm, [genm, lg(Rim)]))
@printf("%-40s %+.3f\n", "| gen, R_i, mass",      prank(lg(aPm), compm, [genm, lg(Rim), lg(Wwm)]))
@printf("%-40s %+.3f\n", "| gen, R_i, dimension", prank(lg(aPm), compm, [genm, lg(Rim), dimm]))
@printf("%-40s %+.3f\n", "| gen, R_i, mass, dim", prank(lg(aPm), compm, [genm, lg(Rim), lg(Wwm), dimm]))
println("\nIf a_p->compensation survives the matrix-dimension control, it is a real")
println("reproduction-timing/strategy signal, not a stage-count artifact.")
