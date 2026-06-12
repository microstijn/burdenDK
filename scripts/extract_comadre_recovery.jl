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
# unmanipulated matrix the damping ratio (|lambda1|/|lambda2|, the rate of
# return to stable structure) and the Caswell generation time (log R0 / log
# lambda1), and (3) writes the species-aggregated table to
# data/external/comadre_recovery.csv (committed; the raw .RData is gitignored).
# ===========================================================================

using RData, DataFrames, LinearAlgebra, Statistics, Downloads

const URL = "https://compadre-db.org/Data/Download/COMADRE_v.4.26.4.0.RData"
const RAW = joinpath(@__DIR__, "..", "data", "external", "COMADRE_v.4.26.4.0.RData")
const OUT = joinpath(@__DIR__, "..", "data", "external", "comadre_recovery.csv")

isfile(RAW) || (mkpath(dirname(RAW)); Downloads.download(URL, RAW))

d = load(RAW); c = d["comadre"]; md = c["metadata"]; mats = c["mat"]
getcol(col, i) = (v = md[i, col]; v === missing ? "" : string(v))
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

rows = Dict{String, Vector{NTuple{2, Float64}}}(); kept = 0
taxo = Dict{String, NTuple{3, String}}()   # species => (class, order, family)
for i in 1:nrow(md)
    sp = getcol(:SpeciesAccepted, i)
    (sp == "" || occursin(" sp", sp) || occursin("_sp", sp)) && continue
    cap = getcol(:MatrixCaptivity, i); (cap == "W" || cap == "Wild") || continue
    tr = getcol(:MatrixTreatment, i); (tr == "Unmanipulated" || tr == "") || continue
    dr = damping_ratio(mats[i]["matA"]); isfinite(dr) || continue
    gt = generation_time(mats[i]["matU"], mats[i]["matF"], mats[i]["matA"])
    push!(get!(rows, sp, NTuple{2, Float64}[]), (log10(dr), gt)); global kept += 1
    haskey(taxo, sp) || (taxo[sp] = (getcol(:Class, i), getcol(:Order, i), getcol(:Family, i)))
end

clean_field(s) = replace(s, "," => " ")
open(OUT, "w") do io
    println(io, "species,n_matrices,mean_log_damping,mean_generation_time,class,order,family")
    for (sp, vs) in sort(collect(rows); by = first)
        dvals = [v[1] for v in vs]; gvals = [v[2] for v in vs if isfinite(v[2])]
        gt = isempty(gvals) ? "" : string(round(mean(gvals); digits = 4))
        cl, or, fa = get(taxo, sp, ("", "", ""))
        println(io, clean_field(sp), ",", length(vs), ",", round(mean(dvals); digits = 5), ",", gt,
                ",", clean_field(cl), ",", clean_field(or), ",", clean_field(fa))
    end
end
println("kept $kept matrices across ", length(rows), " species -> ", OUT)
