# Representative ECOTOX CAS list for ISIMIP/DynQual-style water quality stressors.
#
# DynQual/ISIMIP water quality variables are aggregate indicators:
# - WT: water temperature, no CAS
# - BOD / organic pollution: aggregate oxygen-demanding organic load, no single CAS
# - FC / pathogen pollution: biological indicator, no CAS
# - TDS / salinity: mixture of dissolved ions/salts, represented here by common salts/ions
# - Nutrient/Chemical/Plastic placeholders: represented by common nutrient, metal, organic, and polymer proxies

using Pkg
try
    Pkg.activate(joinpath(@__DIR__, ".."))
catch err
    @warn "Could not activate project" exception=(err, catch_backtrace())
end

using CSV
using DataFrames

include(joinpath(@__DIR__, "..", "src", "ECOTOXParser.jl"))
using .ECOTOXParser

println("\n--- Building ECOTOX toxicity library ---")

# Adjust this if your files are inside validation/
ecotox_dir = joinpath(@__DIR__, "..", "data", "ecotox", "ecotox_ascii_03_12_2026")

results_path = joinpath(ecotox_dir, "results.txt")
tests_path   = joinpath(ecotox_dir, "tests.txt")
species_path = joinpath(ecotox_dir, "validation", "species.txt")

# If your files are actually in validation/, uncomment this block:
# ecotox_dir = joinpath(@__DIR__, "..", "data", "ecotox", "ecotox_ascii_03_12_2026", "validation")
# results_path = joinpath(ecotox_dir, "results.txt")
# tests_path   = joinpath(ecotox_dir, "tests.txt")
# species_path = joinpath(ecotox_dir, "species.txt")

@assert isfile(results_path) "Missing results file: $results_path"
@assert isfile(tests_path)   "Missing tests file: $tests_path"
@assert isfile(species_path) "Missing species file: $species_path"

# Start with a small explicit CAS list.
# CAS may be hyphenated or digits-only.
cas_list = [
    # -------------------------------------------------------------------------
    # Salinity / TDS proxies
    # -------------------------------------------------------------------------
    "7647-14-5",    # Sodium chloride
    "10043-52-4",  # Calcium chloride
    "7786-30-3",   # Magnesium chloride
    "7757-82-6",   # Sodium sulfate
    "7778-18-9",   # Calcium sulfate

    # -------------------------------------------------------------------------
    # Nutrient / eutrophication / oxygen-demand related proxies
    # -------------------------------------------------------------------------
    "7664-41-7",   # Ammonia
    "14798-03-9",  # Ammonium ion
    "14797-55-8",  # Nitrate ion
    "14797-65-0",  # Nitrite ion
    "14265-44-2",  # Orthophosphate / phosphate ion
    "7723-14-0",   # Phosphorus

    # -------------------------------------------------------------------------
    # Organic pollution / BOD representative compounds
    # BOD itself has no CAS; these are proxy oxygen-demanding or wastewater organics.
    # -------------------------------------------------------------------------
    "50-00-0",     # Formaldehyde
    "67-56-1",     # Methanol
    "64-17-5",     # Ethanol
    "64-19-7",     # Acetic acid
    "108-95-2",    # Phenol

    # -------------------------------------------------------------------------
    # General chemical toxicity proxies / common aquatic toxicants
    # Useful for the Chemical placeholder layer and ECOTOX attack calibration.
    # -------------------------------------------------------------------------
    "7440-43-9",   # Cadmium
    "7440-50-8",   # Copper
    "7440-66-6",   # Zinc
    "7439-92-1",   # Lead
    "7439-97-6",   # Mercury
    "7440-38-2",   # Arsenic
    "7440-02-0",   # Nickel
    "7440-47-3",   # Chromium

    # -------------------------------------------------------------------------
    # Plastic / polymer proxies
    # Plastic is not a single chemical; these are representative polymers.
    # ECOTOX coverage may be sparse or inconsistent.
    # -------------------------------------------------------------------------
    "9002-88-4",   # Polyethylene
    "9003-07-0",   # Polypropylene
    "9003-53-6",   # Polystyrene
    "9002-86-2",   # Polyvinyl chloride
    "25038-59-9"   # Polyethylene terephthalate
]

println("CAS list:")
for cas in cas_list
    println("  ", cas, " -> ", normalize_cas(cas), " -> ", hyphenate_cas(cas))
end

output_path = joinpath(
    @__DIR__,
    "..",
    "data",
    "ECOTOX_Toxicity_Library.json"
)

summary = build_ecotox_toxicity_library_multi(
    results_path,
    tests_path,
    species_path,
    cas_list;
    output_path = output_path,
    skip_empty = true
)



println("\nRows in ECOTOX toxicity summary: ", nrow(summary))
println("Output written to: ", output_path)

if nrow(summary) > 0
    println("\nFirst rows:")
    show(first(summary, min(10, nrow(summary))); allcols=true, allrows=true)
    println()
else
    println("\nNo records found for the selected CAS list.")
end

println("\nDone.")

using JSON
records = JSON.parsefile(joinpath(@__DIR__, "..", "data", "ECOTOX_Toxicity_Library.json"))

length(records)
first(records)

usable = filter(r ->
    haskey(r, "NOEC_median") &&
    haskey(r, "EC50_median") &&
    r["NOEC_median"] !== nothing &&
    r["EC50_median"] !== nothing &&
    r["EC50_median"] > r["NOEC_median"],
    records
)

length(usable)
first(usable)