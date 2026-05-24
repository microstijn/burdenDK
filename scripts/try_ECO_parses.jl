include(joinpath(@__DIR__, "..", "src", "ECOTOXParser.jl"))

using .ECOTOXParser

base = joinpath(@__DIR__, "..", "data", "ecotox", "ecotox_ascii_03_12_2026")

results_path = joinpath(base, "results.txt")
tests_path   = joinpath(base, "tests.txt")
species_path = joinpath(base, "validation", "species.txt")

target_cas = "50000"  # example cadmium; replace as needed

summary = build_ecotox_toxicity_library(
    results_path,
    tests_path,
    species_path,
    target_cas;
    output_path = joinpath("data", "ECOTOX_Toxicity_Library_$(replace(target_cas, "-" => "_")).json")
)

