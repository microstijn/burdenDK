# ===========================================================================
# Fetch an Open Tree of Life induced subtree for the AmP-COMADRE matched species
# (Idea A — phylogenetic/PGLS validation). STANDALONE — needs HTTP.jl + JSON,
# NOT part of the project Manifest. Run in a throwaway env:
#
#   julia +release -e 'using Pkg; Pkg.activate(mktempdir()); Pkg.add(["HTTP","JSON"])'
#   julia +release --project=<that-env> scripts/fetch_comadre_tree.jl
#
# Pipeline: (1) read the matched table's amp_key column -> species names;
# (2) OTL TNRS match_names -> OTT ids (Animals context, fuzzy on); (3) OTL
# tree_of_life/induced_subtree -> Newick topology, dropping ott_ids the
# synthetic tree reports as unknown/broken and retrying.
#
# OTL gives TOPOLOGY ONLY (no dated branch lengths) — the PGLS step assigns
# Grafen (1989) ultrametric branch lengths and estimates Pagel's lambda.
#
# Reads:  data/external/comadre_amp_matched.csv
# Writes: data/external/comadre_amp_tree.nwk        (Newick, committed)
#         data/external/comadre_amp_tree_map.csv    (amp_key,ott_label,matched)
# Open Tree of Life (opentreeoflife.org), synthetic tree; CC0/CC-BY per source.
# ===========================================================================

using HTTP, JSON

const TBL  = joinpath(@__DIR__, "..", "data", "external", "comadre_amp_matched.csv")
const NWK  = joinpath(@__DIR__, "..", "data", "external", "comadre_amp_tree.nwk")
const TMAP = joinpath(@__DIR__, "..", "data", "external", "comadre_amp_tree_map.csv")
const OTL  = "https://api.opentreeoflife.org/v3"

post(path, body) = HTTP.post("$OTL/$path", ["Content-Type" => "application/json"],
                             JSON.json(body); readtimeout = 60, retries = 2)

# amp_key -> "Genus species"
amp_keys = String[]
for line in Iterators.drop(eachline(TBL), 1)
    push!(amp_keys, String(split(line, ",")[1]))
end
names = [replace(k, "_" => " ") for k in amp_keys]
println("species to place: ", length(names))

# --- 1. TNRS match -> ott_id per name -------------------------------------
r = post("tnrs/match_names",
         Dict("names" => names, "do_approximate_matching" => true, "context_name" => "Animals"))
res = JSON.parse(String(r.body))
key_ott = Dict{String, Int}()      # amp_key => ott_id
for item in res["results"]
    ms = item["matches"]; isempty(ms) && continue
    nm = item["name"]; ott = ms[1]["taxon"]["ott_id"]
    idx = findfirst(==(nm), names)
    idx === nothing || (key_ott[amp_keys[idx]] = ott)
end
println("TNRS matched: ", length(key_ott), " / ", length(names))

# --- 2. induced subtree, dropping ott_ids not in the synthetic tree -------
ott_key = Dict(v => k for (k, v) in key_ott)   # ott_id => amp_key (1:1 enough)
ids = collect(values(key_ott))
local newick = ""
for attempt in 1:6
    try
        rr = post("tree_of_life/induced_subtree", Dict("ott_ids" => ids))
        global newick = JSON.parse(String(rr.body))["newick"]
        break
    catch e
        e isa HTTP.Exceptions.StatusError || rethrow()
        body = JSON.parse(String(e.response.body))
        bad = Int[]
        for fld in ("unknown", "broken", "node_ids_not_in_tree", "unknown_ids")
            if haskey(body, fld)
                v = body[fld]
                v isa AbstractDict ? append!(bad, parse.(Int, replace.(collect(keys(v)), "ott" => ""))) :
                v isa AbstractVector ? append!(bad, [x isa Integer ? x : parse(Int, replace(string(x), "ott" => "")) for x in v]) : nothing
            end
        end
        isempty(bad) && (println("induced_subtree failed, no droppable ids: ", body); rethrow())
        filter!(x -> !(x in bad), ids)
        println("attempt $attempt: dropped ", length(bad), " ott_ids not in tree; retrying with ", length(ids))
    end
end
isempty(newick) && error("could not obtain induced subtree")

open(NWK, "w") do io; println(io, newick); end

# --- 3. record which amp_keys made it onto the tree -----------------------
placed = Set(ids)
open(TMAP, "w") do io
    println(io, "amp_key,ott_id,placed")
    for k in amp_keys
        ott = get(key_ott, k, nothing)
        println(io, k, ",", ott === nothing ? "" : ott, ",",
                (ott !== nothing && ott in placed))
    end
end
nplaced = count(in(placed), values(key_ott))
println("placed on tree: ", nplaced, " / ", length(names), " -> ", NWK)
