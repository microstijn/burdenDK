# ===========================================================================
# Extract a per-species demographic-recovery table from the COMADRE animal
# matrix database, for the external validation in
# examples/comadre_partial_validation.jl.
#
# STANDALONE — needs RData.jl + DataFrames, NOT part of the project Manifest.
# Run it in a throwaway environment so it does not touch the project deps:
#
#   julia -e 'using Pkg; Pkg.activate(mktempdir()); Pkg.add(["RData","DataFrames"])'  # once
#   # then run this file with that environment active.
#
# It (1) downloads COMADRE_v.4.26.4.0.RData (open access, CC-BY; cite
# Salguero-Gomez et al. 2016, J. Anim. Ecol.), (2) computes per wild,
# unmanipulated matrix the damping ratio (|lambda1|/|lambda2|, recovery), the
# Caswell generation time (log R0 / log lambda1), and the Stott et al. (2011)
# first-timestep transient bounds on Ahat=A/lambda1 — reactivity (compensation)
# and attenuation (resistance) — for the per-axis validation (Idea B), and (3)
# writes the species-aggregated table to data/external/comadre_recovery.csv
# (committed; the raw .RData is gitignored).
# ===========================================================================

using RData, DataFrames, LinearAlgebra, Statistics, Downloads

const URL = "https://compadre-db.org/Data/Download/COMADRE_v.4.26.4.0.RData"
const RAW = joinpath(@__DIR__, "..", "data", "external", "COMADRE_v.4.26.4.0.RData")
const OUT = joinpath(@__DIR__, "..", "data", "external", "comadre_recovery.csv")

isfile(RAW) || (mkpath(dirname(RAW)); Downloads.download(URL, RAW))

d = load(RAW); c = d["comadre"]; md = c["metadata"]; mats = c["mat"]
getcol(col, i) = (v = md[i, col]; v === missing ? "" : string(v))
clean(A) = any(ismissing, A) ? nothing : Float64.(A)

# Per-matrix transient/asymptotic metrics, all from matA:
#   damping ratio |lambda1|/|lambda2| (recovery), and the Stott et al. (2011)
#   standardised first-timestep bounds on Ahat = A/lambda1 — reactivity P̄1 (max
#   column sum, = compensation: max one-step amplification) and attenuation P̲1
#   (min column sum, = resistance: max one-step decline). Returns (dr, reac, att).
function matrix_metrics(Araw)
    A = clean(Araw); A === nothing && return (NaN, NaN, NaN)
    (size(A, 1) == size(A, 2) && size(A, 1) >= 2 && all(isfinite, A) && all(>=(0), A)) || return (NaN, NaN, NaN)
    mods = sort(abs.(eigvals(A)); rev = true)
    (mods[1] > 0 && mods[2] > 0 && mods[1] > mods[2]) || return (NaN, NaN, NaN)
    cs = vec(sum(A ./ mods[1]; dims = 1))
    return (mods[1] / mods[2], maximum(cs), minimum(cs))
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

logornan(x) = (isfinite(x) && x > 0) ? log10(x) : NaN
rows = Dict{String, Vector{NTuple{4, Float64}}}(); kept = 0
taxo = Dict{String, NTuple{3, String}}()   # species => (class, order, family)
for i in 1:nrow(md)
    sp = getcol(:SpeciesAccepted, i)
    (sp == "" || occursin(" sp", sp) || occursin("_sp", sp)) && continue
    cap = getcol(:MatrixCaptivity, i); (cap == "W" || cap == "Wild") || continue
    tr = getcol(:MatrixTreatment, i); (tr == "Unmanipulated" || tr == "") || continue
    dr, reac, att = matrix_metrics(mats[i]["matA"]); isfinite(dr) || continue
    gt = generation_time(mats[i]["matU"], mats[i]["matF"], mats[i]["matA"])
    push!(get!(rows, sp, NTuple{4, Float64}[]), (log10(dr), gt, logornan(reac), logornan(att)))
    global kept += 1
    haskey(taxo, sp) || (taxo[sp] = (getcol(:Class, i), getcol(:Order, i), getcol(:Family, i)))
end

clean_field(s) = replace(s, "," => " ")
meanstr(xs) = isempty(xs) ? "" : string(round(mean(xs); digits = 5))
open(OUT, "w") do io
    println(io, "species,n_matrices,mean_log_damping,mean_generation_time,class,order,family,",
                "mean_log_reactivity,mean_log_attenuation")
    for (sp, vs) in sort(collect(rows); by = first)
        dvals = [v[1] for v in vs]; gvals = [v[2] for v in vs if isfinite(v[2])]
        rvals = [v[3] for v in vs if isfinite(v[3])]; avals = [v[4] for v in vs if isfinite(v[4])]
        gt = isempty(gvals) ? "" : string(round(mean(gvals); digits = 4))
        cl, or, fa = get(taxo, sp, ("", "", ""))
        println(io, clean_field(sp), ",", length(vs), ",", round(mean(dvals); digits = 5), ",", gt,
                ",", clean_field(cl), ",", clean_field(or), ",", clean_field(fa),
                ",", meanstr(rvals), ",", meanstr(avals))
    end
end
println("kept $kept matrices across ", length(rows), " species -> ", OUT)
