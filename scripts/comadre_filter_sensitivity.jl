# ===========================================================================
# Matrix-quality filter sensitivity (roadmap Part 3): is the k_M -> demographic
# recovery signal robust to how COMADRE matrices are filtered? Re-derives the
# generation-time-controlled partial rank correlation rho(k_M, recovery | gen)
# under several matrix-quality filter variants and reports them side by side.
#
# STANDALONE — needs RData.jl + DataFrames (same throwaway env as
# scripts/extract_comadre_recovery.jl). Self-contained: joins AmP k_M via the
# committed CSVs (no project/library dependency).
#   julia +release --project=<rdata-env> scripts/comadre_filter_sensitivity.jl
#
# Reads the gitignored raw COMADRE_v.4.26.4.0.RData (re-fetched by the extractor),
# data/external/comadre_amp_namemap.csv and data/external/amp_reproduction_rates.csv.
# ===========================================================================

using RData, DataFrames, LinearAlgebra, Statistics

const RAW = joinpath(@__DIR__, "..", "data", "external", "COMADRE_v.4.26.4.0.RData")
const MAP = joinpath(@__DIR__, "..", "data", "external", "comadre_amp_namemap.csv")
const RATES = joinpath(@__DIR__, "..", "data", "external", "amp_reproduction_rates.csv")

d = load(RAW); c = d["comadre"]; md = c["metadata"]; mats = c["mat"]
getcol(col, i) = (hasproperty(md, col) ? (v = md[i, col]; v === missing ? "" : string(v)) : "")
clean(A) = any(ismissing, A) ? nothing : Float64.(A)

function damping_ratio(Araw)
    A = clean(Araw); A === nothing && return NaN
    (size(A, 1) == size(A, 2) && size(A, 1) >= 2 && all(isfinite, A) && all(>=(0), A)) || return NaN
    mods = sort(abs.(eigvals(A)); rev = true)
    (mods[1] > 0 && mods[2] > 0 && mods[1] > mods[2]) || return NaN
    return mods[1] / mods[2]
end
function generation_time(Uraw, Fraw, Araw)
    U = clean(Uraw); F = clean(Fraw); A = clean(Araw)
    (U === nothing || F === nothing || A === nothing) && return NaN
    try
        k = size(U, 1); N = inv(I(k) - U)
        R0 = maximum(abs.(eigvals(F * N))); lam = maximum(abs.(eigvals(A)))
        (R0 > 0 && lam > 0 && !isapprox(lam, 1.0; atol = 1e-6)) || return NaN
        return log(R0) / log(lam)
    catch; return NaN; end
end

# AmP k_M per resolved species (via namemap -> amp_key -> rates)
namemap = Dict{String, String}()
for line in Iterators.drop(eachline(MAP), 1)
    f = split(line, ","); length(f) >= 2 && f[2] != "" && (namemap[String(f[1])] = String(f[2]))
end
kM_of = Dict{String, Float64}()
for line in Iterators.drop(eachline(RATES), 1)
    f = split(line, ","); length(f) >= 2 || continue
    kM_of[String(f[1])] = parse(Float64, f[2])
end

ordinalrank(v) = (p = sortperm(v); r = similar(p); r[p] = 1:length(v); Float64.(r))
function pspear(x, y, z)
    rx, ry, rz = ordinalrank(x), ordinalrank(y), ordinalrank(z)
    rab = cor(rx, ry); rac = cor(rx, rz); rbc = cor(ry, rz)
    (rab - rac * rbc) / sqrt((1 - rac^2) * (1 - rbc^2))
end
spear(x, y) = cor(ordinalrank(x), ordinalrank(y))

# evaluate one filter variant: keep(i) is a predicate on metadata row i.
function run_variant(keep)
    rows = Dict{String, Vector{NTuple{2, Float64}}}()
    nmat = 0
    for i in 1:nrow(md)
        sp = getcol(:SpeciesAccepted, i)
        (sp == "" || occursin(" sp", sp) || occursin("_sp", sp)) && continue
        keep(i) || continue
        dr = damping_ratio(mats[i]["matA"]); isfinite(dr) || continue
        gt = generation_time(mats[i]["matU"], mats[i]["matF"], mats[i]["matA"])
        spc = replace(sp, "," => " ")
        push!(get!(rows, spc, NTuple{2, Float64}[]), (log10(dr), gt)); nmat += 1
    end
    kM = Float64[]; rec = Float64[]; gen = Float64[]
    for (spc, vs) in rows
        key = get(namemap, spc, nothing); key === nothing && continue
        haskey(kM_of, key) || continue
        gvals = [v[2] for v in vs if isfinite(v[2])]
        push!(kM, kM_of[key]); push!(rec, mean(v[1] for v in vs))
        push!(gen, isempty(gvals) ? NaN : mean(gvals))
    end
    m = isfinite.(gen)
    return (nmat = nmat, nsp = length(kM), nsp_gen = count(m),
            raw = spear(log.(kM), rec), partial = pspear(log.(kM[m]), rec[m], gen[m]))
end

base(i) = (cap = getcol(:MatrixCaptivity, i); (cap == "W" || cap == "Wild")) &&
          (tr = getcol(:MatrixTreatment, i); (tr == "Unmanipulated" || tr == ""))
matdim(i) = (s = getcol(:MatrixDimension, i); s == "" ? 0 : (v = tryparse(Int, s); v === nothing ? 0 : v))
composite(i) = getcol(:MatrixComposite, i)

variants = [
    ("baseline (wild, unmanipulated)", base),
    ("+ individual matrices only",     i -> base(i) && composite(i) == "Individual"),
    ("+ composite (Mean/Pooled) only", i -> base(i) && composite(i) in ("Mean", "Pooled", "Seasonal")),
    ("+ dimension >= 3",               i -> base(i) && matdim(i) >= 3),
    ("+ dimension == 2 only",          i -> base(i) && matdim(i) == 2),
    ("ALL captivity (looser)",         i -> (tr = getcol(:MatrixTreatment, i); tr == "Unmanipulated" || tr == "")),
]

println("Sensitivity of rho(k_M, recovery | generation time) to matrix-quality filters:\n")
println(rpad("variant", 36), rpad("n_mat", 8), rpad("n_sp", 7), rpad("raw rho", 10), "partial rho | gen")
for (lab, keep) in variants
    r = run_variant(keep)
    println(rpad(lab, 36), rpad(r.nmat, 8), rpad(r.nsp, 7),
            rpad(round(r.raw; digits = 3), 10), round(r.partial; digits = 3), "  (n=", r.nsp_gen, ")")
end
println("\nStable partial rho across variants => the k_M signal is not an artifact")
println("of a particular matrix-quality filter.")
