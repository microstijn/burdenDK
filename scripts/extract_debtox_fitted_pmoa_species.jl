# extract_debtox_fitted_pmoa_species.jl
# Species-resolved fitted-pMoA backbone for the CURATED APICAL-EC50 database (across-axis weighting
# test). Parses the ESPI2018 DEBtox compilation (Kienzler/Jager-adjacent, doi 10.1039/c7em00328e),
# which carries the CASCADE-AWARE pMoA per (compound, SPECIES) from full DEBtox fits -- the correct
# axis assignment (Jager 2020: the most-sensitive apical endpoint != pMoA, because of the kappa-rule
# starvation cascade). The compilation gives the AXIS; fitted potencies (ECx/NEC) live in the
# original papers (the DOI column) -> the shopping list printed below.
#
# The source xlsx is gitignored (raw download). Unzip it first, then run:
#   unzip -o ~/Downloads/dat_1/c7em00328e1.xlsx -d /tmp/debtox_xlsx
#   julia +release --project=. scripts/extract_debtox_fitted_pmoa_species.jl
import JSON
using Printf

const DIR = get(ENV, "DEBTOX_XLSX_DIR", "/tmp/debtox_xlsx")
const OUT = "data/external/debtox_fitted_pmoa_species_c7em00328e.csv"

# --- shared strings (index -> text) ---
ss_raw = read(joinpath(DIR, "xl", "sharedStrings.xml"), String)
shared = String[]
for m in eachmatch(r"<si>(.*?)</si>"s, ss_raw)
    txt = join((t.captures[1] for t in eachmatch(r"<t[^>]*>([^<]*)</t>", m.captures[1])), "")
    push!(shared, txt)
end

# --- sheet1 cells -> (row, col_letter) => value ---
sheet = read(joinpath(DIR, "xl", "worksheets", "sheet1.xml"), String)
cells = Dict{Tuple{Int,String},String}()
for m in eachmatch(r"<c r=\"([A-Z]+)(\d+)\"([^>]*?)>(.*?)</c>"s, sheet)
    col, row, attrs, body = m.captures[1], parse(Int, m.captures[2]), m.captures[3], m.captures[4]
    vm = match(r"<v>([^<]*)</v>", body); vm === nothing && continue
    val = occursin("t=\"s\"", attrs) ? shared[parse(Int, vm.captures[1]) + 1] : String(vm.captures[1])
    cells[(row, col)] = strip(val)
end
getc(r, c) = get(cells, (r, c), "")

# --- reconstruct table (cols B..G from row 4); forward-fill merged Ref/Species/Scientific ---
norm(s) = replace(strip(replace(s, r"\s+" => " ")), " " => "_")
A2axis = Dict("A" => "assimilation", "M" => "maintenance", "G" => "growth", "R" => "reproduction", "H" => "hazard")
primary_axis(p) = (c = match(r"[AMGRH]", p); c === nothing ? "" : A2axis[c.match])

rows = NamedTuple[]
ref = spc = sci = ""
for r in 4:60
    global ref, spc, sci
    b, c, d = getc(r, "B"), getc(r, "C"), getc(r, "D")
    isempty(b) || (ref = b); isempty(c) || (spc = c); isempty(d) || (sci = d)
    compound, pmoa, doi = getc(r, "E"), getc(r, "F"), getc(r, "G")
    (isempty(compound) || isempty(pmoa)) && continue
    push!(rows, (ref = ref, species_common = spc, scientific = sci, amp_key = norm(sci),
                 compound = compound, pmoa = pmoa, axis = primary_axis(pmoa), doi = doi))
end
println("parsed fitted rows: ", length(rows))

# --- AmP kappa ---
lib = JSON.parsefile("data/AmP_Species_Library.json")
kappa_of = Dict{String,Float64}()
for (sp, v) in lib
    (v isa AbstractDict && haskey(v, "alpha_axes") && length(v["alpha_axes"]) == 4) || continue
    k = float(v["alpha_axes"][3]); (isfinite(k) && 0 < k < 1) && (kappa_of[sp] = k)
end
inamp(k) = haskey(kappa_of, k)

open(OUT, "w") do io
    println(io, "ref,species_common,scientific_name,amp_key,in_amp,kappa,compound,pmoa,primary_axis,doi")
    for r in rows
        println(io, join([r.ref, r.species_common, r.scientific, r.amp_key, inamp(r.amp_key),
            inamp(r.amp_key) ? round(kappa_of[r.amp_key], digits = 3) : "",
            r.compound, r.pmoa, r.axis, r.doi], ","))
    end
end
println("wrote $OUT")

# --- AmP-resident summary + shopping list ---
amp_rows = [r for r in rows if inamp(r.amp_key)]
sps = sort(unique(r.amp_key for r in amp_rows))
println("\n=== AmP-resident fitted (species, compound) pairs: ", length(amp_rows),
        "  across ", length(sps), " AmP species ===")
for s in sps
    rs = [r for r in amp_rows if r.amp_key == s]
    axset = sort(unique(r.axis for r in rs))
    @printf("  %-28s kappa=%.3f  axes=%-30s  %d compound(s)\n",
            s, kappa_of[s], join(axset, "/"), length(rs))
end

# species with >=2 distinct discriminating axes (the multi-MoA targets)
disc = ("maintenance", "reproduction", "growth", "assimilation")
multi = [s for s in sps if length(unique(r.axis for r in amp_rows if r.amp_key == s && r.axis in disc)) >= 2]
println("\n  AmP species with >=2 distinct fitted axes: ", length(multi), " -> ", join(multi, ", "))

println("\n=== SHOPPING LIST (fitted ECx/NEC to retrieve), grouped by paper DOI ===")
bydoi = Dict{String,Vector{NamedTuple}}()
for r in amp_rows; push!(get!(bydoi, r.doi, NamedTuple[]), r); end
for (doi, rs) in sort(collect(bydoi); by = x -> -length(x[2]))
    println("  ", rs[1].ref, "  (", doi, ")")
    for r in rs
        @printf("      %-26s  %-22s  pMoA=%-4s -> %s\n", r.amp_key, r.compound, r.pmoa, r.axis)
    end
end
