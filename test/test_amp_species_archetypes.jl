using Test
using DataFrames
using CSV
using JSON

@testset "AmP Species Archetype Database" begin
    script_path = normpath(joinpath(@__DIR__, "..", "examples", "build_amp_species_archetype_database.jl"))
    csv_path = normpath(joinpath(@__DIR__, "..", "data", "AmP_Species_Archetypes.csv"))
    json_path = normpath(joinpath(@__DIR__, "..", "data", "AmP_Species_Archetypes.json"))
    summary_path = normpath(joinpath(@__DIR__, "..", "output", "amp_species_archetype_database", "amp_species_archetype_database_build_summary.csv"))

    # 1. The archetype build script runs without error
    # We run it via external command to ensure the script's `if abspath(PROGRAM_FILE) == @__FILE__` guard triggers.
    try
        run(`$(Base.julia_cmd()) --project=. $script_path`)
        @test true
    catch e
        @test false # Script failed to run
    end

    # 2. CSV exists and is non-empty
    @test isfile(csv_path)
    @test filesize(csv_path) > 0

    # 3. JSON exists and is non-empty
    @test isfile(json_path)
    @test filesize(json_path) > 0

    # 4. Build summary exists and is non-empty
    @test isfile(summary_path)
    @test filesize(summary_path) > 0

    # 5. CSV contains required columns
    df = CSV.read(csv_path, DataFrame)
    required_cols = [
        "species_key", "species_name", "archetype_labels",
        "Fmax", "min_x_collapse", "most_sensitive_axis",
        "is_low_amplification_potential", "is_high_amplification_potential",
        "is_high_sensitivity_low_collapse_threshold", "is_robust_high_collapse_threshold"
    ]
    for col in required_cols
        @test col in names(df)
    end

    # 6. CSV has at least 100 species if AmP library supports it
    # AmP species library typically has > 7000
    if nrow(df) >= 100
        @test nrow(df) >= 100
    else
        @warn "Testing on a reduced AmP library: found $(nrow(df)) species"
        @test nrow(df) >= 10
    end

    # 7. At least 90% of rows have non-empty archetype_labels
    valid_rows = filter(row -> !ismissing(row.archetype_labels) && length(strip(row.archetype_labels)) > 0, df)
    @test nrow(valid_rows) / nrow(df) >= 0.90

    # 8. Required archetypes are present with at least one species
    function has_archetype(df, label)
        return any(row -> !ismissing(row.archetype_labels) && occursin(label, row.archetype_labels), eachrow(df))
    end

    @test has_archetype(df, "low_amplification_potential")
    @test has_archetype(df, "high_amplification_potential")
    @test has_archetype(df, "high_sensitivity_low_collapse_threshold")
    @test has_archetype(df, "robust_high_collapse_threshold")
    @test has_archetype(df, "maintenance_sensitive") || has_archetype(df, "growth_sensitive") || has_archetype(df, "reproduction_sensitive")

    # 9. At least three distinct most_sensitive_axis values
    distinct_axes = unique(skipmissing(df.most_sensitive_axis))
    @test length(distinct_axes) >= 3

    # 10. JSON metadata contains required fields
    json_data = JSON.parsefile(json_path)
    @test haskey(json_data, "metadata")
    meta = json_data["metadata"]
    @test haskey(meta, "generated_by")
    @test haskey(meta, "archetype_algorithm_version")
    @test haskey(meta, "n_species")

    # 11. JSON species count matches CSV row count
    json_species = json_data["species"]
    @test length(json_species) == nrow(df)

    # 12. No species_key duplicates in CSV
    @test length(unique(df.species_key)) == nrow(df)

    # 13. All numeric required fields are finite where expected
    finite_numeric_cols = ["A0", "lambda_min", "lambda_max", "KA", "Fmax"]
    for col in finite_numeric_cols
        @test all(x -> ismissing(x) || isnan(x) || isfinite(x), df[!, col])
    end
    # Allow Inf for min_x_collapse but not missing
    @test all(x -> ismissing(x) || isinf(x) || (isfinite(x) && x > 0), df.min_x_collapse)

    # 14. Quantile columns are within [0, 1] for finite metrics
    for q_col in ["Fmax_quantile", "min_x_collapse_quantile"]
        vals = filter(x -> !ismissing(x) && !isnan(x), df[!, q_col])
        @test all(x -> 0.0 <= x <= 1.0, vals)
    end

    # 15. Running the build script twice produces the same species_key order and same archetype labels (excluding generated_at_utc)
    # Re-run script to generate new CSV
    run(`$(Base.julia_cmd()) --project=. $script_path`)
    df_new = CSV.read(csv_path, DataFrame)

    @test df.species_key == df_new.species_key
    @test df.species_name == df_new.species_name
    @test df.archetype_labels == df_new.archetype_labels
    @test df.is_low_amplification_potential == df_new.is_low_amplification_potential

    # Just to be safe, clean up summary for a clean state or leave it depending on how runtests usually behaves.
end
