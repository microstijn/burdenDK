using Test
using DataFrames
using CSV
using TwoTimescaleResilience

@testset "Mixture Effect Overlap Demo" begin
    # 1. The example script runs without error
    script_path = normpath(joinpath(@__DIR__, "..", "examples", "mixture_effect_model_overlap_demo.jl"))

    # We execute it to ensure outputs are fresh and test it runs cleanly
    @test try
        include(script_path)
        true
    catch e
        println("Error running demo script: ", e)
        false
    end

    out_dir = normpath(joinpath(@__DIR__, "..", "output", "mixture_effect_model_overlap_demo"))
    summary_path = joinpath(out_dir, "mixture_effect_model_overlap_summary.csv")
    long_path = joinpath(out_dir, "mixture_effect_model_overlap_axis_long.csv")
    plot_path = joinpath(out_dir, "mixture_effect_model_overlap_comparison.png")

    # 2 & 3. Output CSV and PNG exist and are non-empty
    @test isfile(summary_path) && filesize(summary_path) > 0
    @test isfile(long_path) && filesize(long_path) > 0
    @test isfile(plot_path) && filesize(plot_path) > 0

    # Read the summary dataframe for testing
    df = CSV.read(summary_path, DataFrame)

    # 4. Summary CSV contains all three mixture_effect_model values
    models = unique(df.mixture_effect_model)
    @test "axis_toxic_unit_sum" in models
    @test "independent_action_axis_effects" in models
    @test "grouped_ca_then_ia_axis_effects" in models

    # Helper to get E_maintenance
    function get_E_m(scenario::String, model::String)
        return df[(df.scenario .== scenario) .& (df.mixture_effect_model .== model), :E_maintenance][1]
    end

    # Helper to get X_maintenance
    function get_X_m(scenario::String, model::String)
        return df[(df.scenario .== scenario) .& (df.mixture_effect_model .== model), :X_maintenance][1]
    end

    # 5. Single-contributor scenario
    @test get_E_m("single_contributor_control", "axis_toxic_unit_sum") ≈ 0.5
    @test get_E_m("single_contributor_control", "independent_action_axis_effects") ≈ 0.5
    @test get_E_m("single_contributor_control", "grouped_ca_then_ia_axis_effects") ≈ 0.5

    # 6. Two distinct effect groups scenario
    @test get_E_m("two_distinct_effect_groups", "axis_toxic_unit_sum") ≈ 2/3
    @test get_E_m("two_distinct_effect_groups", "independent_action_axis_effects") ≈ 0.75
    @test get_E_m("two_distinct_effect_groups", "grouped_ca_then_ia_axis_effects") ≈ 0.75

    # 7. Same effect group scenario
    @test get_E_m("two_same_effect_group", "axis_toxic_unit_sum") ≈ 2/3
    @test get_E_m("two_same_effect_group", "independent_action_axis_effects") ≈ 0.75
    @test get_E_m("two_same_effect_group", "grouped_ca_then_ia_axis_effects") ≈ 2/3

    # 8. Mixed grouping scenario
    val_tu = get_E_m("mixed_grouping", "axis_toxic_unit_sum")
    val_ia = get_E_m("mixed_grouping", "independent_action_axis_effects")
    val_grouped = get_E_m("mixed_grouping", "grouped_ca_then_ia_axis_effects")

    @test val_tu ≈ 0.75
    @test val_ia ≈ 0.875
    @test val_grouped ≈ 5/6
    @test val_tu < val_grouped < val_ia

    # 9. Multiple axes scenario
    @test get_E_m("multiple_axes", "axis_toxic_unit_sum") ≈ 2/3
    @test get_E_m("multiple_axes", "independent_action_axis_effects") ≈ 0.75
    @test get_E_m("multiple_axes", "grouped_ca_then_ia_axis_effects") ≈ 0.75

    # growth: all models ≈ 4/5
    for m in models
        val_g = df[(df.scenario .== "multiple_axes") .& (df.mixture_effect_model .== m), :E_growth][1]
        @test val_g ≈ 4/5
    end

    # 10. All E_axis values are bounded: 0 <= E_axis <= 1
    for ax in ["E_assimilation", "E_maintenance", "E_growth", "E_reproduction"]
        col = getproperty(df, Symbol(ax))
        @test all(0.0 .<= col .<= 1.0)
    end

    # 11. X_axis audit values are correct
    @test get_X_m("single_contributor_control", "axis_toxic_unit_sum") ≈ 1.0
    @test get_X_m("two_same_effect_group", "axis_toxic_unit_sum") ≈ 2.0
    @test get_X_m("mixed_grouping", "axis_toxic_unit_sum") ≈ 3.0

    # multiple_axes growth X = 4
    val_X_g = df[(df.scenario .== "multiple_axes") .& (df.mixture_effect_model .== "axis_toxic_unit_sum"), :X_growth][1]
    @test val_X_g ≈ 4.0
end
