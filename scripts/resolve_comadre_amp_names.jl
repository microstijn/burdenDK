# ===========================================================================
# Harmonise COMADRE species names to AmP species-library keys, so the external
# validation (examples/comadre_partial_validation.jl) is not limited to exact
# string matches. Recovers species lost to (a) COMADRE duplicated-genus typos
# ("Genus Genus species"), (b) trinomial subspecies ("Genus species subsp"),
# and (c) taxonomic synonymy/spelling, via the GBIF Backbone matching API.
#
# STANDALONE — needs HTTP.jl + JSON, NOT part of the project Manifest. Run it in
# a throwaway environment so it does not touch the project deps (same pattern as
# scripts/extract_comadre_recovery.jl):
#
#   julia +release -e 'using Pkg; Pkg.activate(mktempdir()); Pkg.add(["HTTP","JSON"])'  # once
#   # then run this file with that environment active:
#   julia +release --project=<that-env> scripts/resolve_comadre_amp_names.jl
#
# Reads:  data/AmP_Species_Library.json        (AmP keys, underscore-separated)
#         data/external/comadre_recovery.csv   (COMADRE species, space-separated)
# Writes: data/external/comadre_amp_namemap.csv (committed)
#           columns: comadre_species,amp_key,method
#         One row per COMADRE species. amp_key empty => unresolved. method is one
#         of: exact, dedup_genus, trinomial, gbif_accepted, gbif_synonym, "".
#
# GBIF Backbone taxonomy (CC-BY): GBIF.org. The /species/match endpoint returns
# the accepted name for synonyms; /species/{key}/synonyms enumerates the cluster.
# ===========================================================================

using HTTP, JSON

const AMP   = joinpath(@__DIR__, "..", "data", "AmP_Species_Library.json")
const COM   = joinpath(@__DIR__, "..", "data", "external", "comadre_recovery.csv")
const OUT   = joinpath(@__DIR__, "..", "data", "external", "comadre_amp_namemap.csv")

amp_keys = Set(keys(JSON.parsefile(AMP)))
und(name::AbstractString) = replace(strip(name), " " => "_")
ampget(name) = (k = und(name); k in amp_keys ? k : nothing)

# COMADRE species list (col 1, skip header).
comadre = String[]
for (i, line) in enumerate(eachline(COM))
    i == 1 && continue
    push!(comadre, String(split(line, ",")[1]))
end

# --- GBIF helpers ---------------------------------------------------------
function gbif_json(url)
    try
        r = HTTP.get(url; readtimeout = 20, retries = 2)
        return JSON.parse(String(r.body))
    catch err
        @warn "GBIF request failed" url err
        return nothing
    end
end

# canonical binomial candidates GBIF offers for a queried name
function gbif_candidates(name)
    cands = String[]
    q = HTTP.escapeuri(name)
    m = gbif_json("https://api.gbif.org/v1/species/match?kingdom=Animalia&name=$q")
    m === nothing && return cands, nothing
    # accepted species binomial (present when matched at/below species rank)
    for f in ("species", "canonicalName")
        haskey(m, f) && m[f] isa AbstractString && push!(cands, m[f])
    end
    # key of the accepted taxon, for the synonym cluster
    acckey = get(m, "acceptedUsageKey", get(m, "speciesKey", get(m, "usageKey", nothing)))
    return unique(cands), acckey
end

function gbif_synonym_names(acckey)
    acckey === nothing && return String[]
    s = gbif_json("https://api.gbif.org/v1/species/$acckey/synonyms?limit=200")
    (s === nothing || !haskey(s, "results")) && return String[]
    out = String[]
    for r in s["results"]
        haskey(r, "canonicalName") && r["canonicalName"] isa AbstractString && push!(out, r["canonicalName"])
    end
    return out
end

# --- resolution pipeline --------------------------------------------------
function resolve(sp)
    # (a) exact
    k = ampget(sp); k !== nothing && return k, "exact"
    toks = split(sp)
    # (b) duplicated leading genus: "Genus Genus species" -> "Genus species"
    if length(toks) >= 3 && toks[1] == toks[2]
        k = ampget(join(toks[2:end], " ")); k !== nothing && return k, "dedup_genus"
    end
    # (c) trinomial subspecies -> binomial (after de-dup, so use first two
    #     *distinct* leading tokens)
    base = (length(toks) >= 3 && toks[1] == toks[2]) ? toks[2:end] : toks
    if length(base) >= 3
        k = ampget(join(base[1:2], " ")); k !== nothing && return k, "trinomial"
    end
    # (d) GBIF synonymy — query the cleaned binomial if we have one, else raw
    query = length(base) >= 2 ? join(base[1:2], " ") : sp
    cands, acckey = gbif_candidates(query)
    for c in cands
        k = ampget(c); k !== nothing && return k, "gbif_accepted"
    end
    for c in gbif_synonym_names(acckey)
        k = ampget(c); k !== nothing && return k, "gbif_synonym"
    end
    return nothing, ""
end

counts = Dict{String, Int}()
open(OUT, "w") do io
    println(io, "comadre_species,amp_key,method")
    for sp in comadre
        k, method = resolve(sp)
        counts[method == "" ? "unresolved" : method] = get(counts, method == "" ? "unresolved" : method, 0) + 1
        println(io, sp, ",", k === nothing ? "" : k, ",", method)
    end
end

# --- report ---------------------------------------------------------------
resolved = length(comadre) - get(counts, "unresolved", 0)
println("COMADRE species: ", length(comadre))
println("resolved to AmP: ", resolved, "  (exact ", get(counts, "exact", 0),
        " + recovered ", resolved - get(counts, "exact", 0), ")")
for m in ("exact", "dedup_genus", "trinomial", "gbif_accepted", "gbif_synonym", "unresolved")
    haskey(counts, m) && println("  ", rpad(m, 14), counts[m])
end
# collision check: >1 COMADRE name mapping to the same AmP key (pseudoreplication)
seen = Dict{String, Vector{String}}()
for line in Iterators.drop(eachline(OUT), 1)
    f = split(line, ","); f[2] == "" && continue
    push!(get!(seen, f[2], String[]), f[1])
end
coll = [(k, v) for (k, v) in seen if length(v) > 1]
if !isempty(coll)
    println("\nWARNING: ", length(coll), " AmP key(s) claimed by >1 COMADRE species (pseudoreplication):")
    for (k, v) in coll; println("  ", k, " <- ", join(v, " ; ")); end
end
println("\nwrote -> ", OUT)
