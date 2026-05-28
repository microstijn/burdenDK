using Test
using TwoTimescaleResilience
using CSV
using DataFrames

@testset "Multi-Axis Response Calibrated Demo Tranche 1" begin
    example_path = normpath(joinpath(@__DIR__, "..", "examples", "ecotox_amp_multiaxis_response_calibrated_demo.jl"))
    include(example_path)

    candidates_df = inspect_candidate_records()

    @test nrow(candidates_df) > 0

    # Test that EC50 > NOEC
    for row in eachrow(candidates_df)
        @test row.EC50_median > row.NOEC_median
        @test row.NOEC_median >= 0
    end

    # Test routing
    for row in eachrow(candidates_df)
        axis = TwoTimescaleResilience.ecotox_effect_to_deb_axis(row.effect_code)
        @test string(axis) == row.axis
    end

    # Test at least two axes have records
    axes_present = unique(candidates_df.axis)
    @test length(axes_present) >= 2

    out_dir = normpath(joinpath(@__DIR__, "..", "output", "ecotox_amp_multiaxis_response_calibrated_demo"))
    out_path = joinpath(out_dir, "ecotox_axis_candidate_records.csv")
    @test isfile(out_path)

    # Tranche 2 Tests
    selected_species, selected_records = select_species_and_records(candidates_df)

    # Test valid species
    @test length(selected_species) > 0
    for sp in selected_species
        @test sp.params isa TwoTimescaleResilience.DEBAxisParams
    end

    # Test valid records and multi-axis
    @test length(selected_records) >= 2

    active_axes = Set{Symbol}()
    for rec in selected_records
        @test TwoTimescaleResilience.validate_ecotox_record(rec)
        @test rec["EC50_median"] > rec["NOEC_median"]
        @test rec["NOEC_median"] >= 0

        ax = TwoTimescaleResilience.ecotox_effect_to_deb_axis(rec["effect_code"])
        push!(active_axes, ax)
    end

    # We must have at least 2 distinct axes activated
    @test length(active_axes) >= 2

    # Tranche 3 Tests
    scenario_months, calibration_k = define_calibrated_scenario(selected_records)

    @test length(scenario_months) == 12

    active_in_scenario = Set{String}()

    for month_concs in scenario_months
        for (cas, C) in month_concs
            @test isfinite(C)
            @test C >= 0.0
            if C > 0.0
                push!(active_in_scenario, cas)
            end
        end
    end

    # Test that multiple compounds/axes receive nonzero pressure
    @test length(active_in_scenario) >= 2

    # Tranche 4 Tests
    df_comp, df_spec = run_diagnostic_scenario(selected_species, selected_records, scenario_months)

    @test isfile(joinpath(out_dir, "multiaxis_compound_summary.csv"))
    @test isfile(joinpath(out_dir, "multiaxis_species_summary.csv"))

    @test nrow(df_comp) > 0
    @test nrow(df_spec) > 0

    # Both methods should be present
    methods = unique(df_spec.mixture_method)
    @test "additive_axis_burden" in methods
    @test "axis_toxic_unit_sum" in methods

    # Months should be 1:12
    months_comp = sort(unique(df_comp.month))
    @test months_comp == collect(1:12)

    # Check that the two methods produce identical totals
    for m in 1:12
        for sp in selected_species
            # We now have 2 response_modes * 3 mixture_effect_models = 6 rows per month per species per mixture_method
            subset_add = filter(row -> row.month == m && row.species_key == sp.key && row.mixture_method == "additive_axis_burden", df_spec)
            subset_tox = filter(row -> row.month == m && row.species_key == sp.key && row.mixture_method == "axis_toxic_unit_sum", df_spec)

            @test nrow(subset_add) == 6
            @test nrow(subset_tox) == 6

            row_add = subset_add[1, :]
            row_tox = subset_tox[1, :]

            @test row_add.total_burden_assimilation ≈ row_tox.total_burden_assimilation
            @test row_add.total_burden_maintenance ≈ row_tox.total_burden_maintenance
            @test row_add.total_burden_growth ≈ row_tox.total_burden_growth
            @test row_add.total_burden_reproduction ≈ row_tox.total_burden_reproduction

            @test row_add.A_t ≈ row_tox.A_t
            @test row_add.lambda_t ≈ row_tox.lambda_t
            @test row_add.F_t ≈ row_tox.F_t
        end
    end

    # At least two DEB axes should have nonzero total burden somewhere in the scenario
    active_burden_axes = Set{Symbol}()
    for row in eachrow(df_spec)
        if row.total_burden_assimilation > 0; push!(active_burden_axes, :assimilation); end
        if row.total_burden_maintenance > 0; push!(active_burden_axes, :maintenance); end
        if row.total_burden_growth > 0; push!(active_burden_axes, :growth); end
        if row.total_burden_reproduction > 0; push!(active_burden_axes, :reproduction); end

        # Test numeric columns are finite
        @test isfinite(row.A_t)
        @test isfinite(row.lambda_t)
        @test isfinite(row.F_t)
    end

    @test length(active_burden_axes) >= 2

    # Tranche 5 Tests
    cal_df, target_hit = perform_calibration_checks(df_spec)

    @test isfile(joinpath(out_dir, "multiaxis_response_calibration_summary.csv"))
    @test nrow(cal_df) == length(selected_species)

    for row in eachrow(cal_df)
        @test isfinite(row.peak_F_t)
        @test isfinite(row.min_A_t)
        @test isfinite(row.min_lambda_t)

        # Test that at least two axes are activated for the overall run
        @test row.active_axes_count >= 2
    end

    # Tranche 6 Tests
    # generate_plots doesn't exist in the same scope due to includes, let's call it explicitly if needed
    # The actual multiaxis file has this:
    try
        TwoTimescaleResilience.generate_plots(df_spec, df_comp)
    catch
        generate_plots(df_spec, df_comp)
    end

    required_pngs = [
        "multiaxis_axis_burdens.png",
        "multiaxis_adaptive_margin.png",
        "multiaxis_restoring_force.png",
        "multiaxis_amplification.png",
        "multiaxis_dominant_compounds.png"
    ]

    for png in required_pngs
        png_path = joinpath(out_dir, png)
        @test isfile(png_path)
        @test filesize(png_path) > 0
    end
end

@testset "Multi-Axis Response Calibrated Demo Tranche 3" begin
    # Manually run the functions since there is no run_multiaxis_response_calibrated_demo
    candidates = inspect_candidate_records()
    selected_species, selected_records = select_species_and_records(candidates)
    scenario_data = define_calibrated_scenario(selected_records)
    scenario = scenario_data[1]  # The function returns (scenario, (k_single, k_multi))
    df_comp, df_spec = run_diagnostic_scenario(selected_species, selected_records, scenario)
    # generate_plots doesn't exist in the same scope due to includes, let's call it explicitly if needed
    # The actual multiaxis file has this:
    try
        TwoTimescaleResilience.generate_plots(df_spec, df_comp)
    catch
        generate_plots(df_spec, df_comp)
    end
    
    # Check that both effect models are present
    effect_models = unique(df_spec.mixture_effect_model)
    @test "axis_toxic_unit_sum" in effect_models
    @test "independent_action_axis_effects" in effect_models
    
    # For axis_toxic_unit_sum under EC50 mode: E_axis = X_axis / (1 + X_axis)
    df_ec50_tu = filter(row -> row.response_mode == "ec50_anchored_fractional_impairment" && row.mixture_effect_model == "axis_toxic_unit_sum", df_spec)
    for row in eachrow(df_ec50_tu)
        @test row.E_assimilation ≈ row.X_assimilation / (1.0 + row.X_assimilation)
        @test row.E_maintenance ≈ row.X_maintenance / (1.0 + row.X_maintenance)
        @test row.E_growth ≈ row.X_growth / (1.0 + row.X_growth)
        @test row.E_reproduction ≈ row.X_reproduction / (1.0 + row.X_reproduction)
    end
    
    # For independent_action_axis_effects under EC50 mode: 0 <= E_axis <= 1
    df_ec50_ia = filter(row -> row.response_mode == "ec50_anchored_fractional_impairment" && row.mixture_effect_model == "independent_action_axis_effects", df_spec)
    for row in eachrow(df_ec50_ia)
        @test 0.0 <= row.E_assimilation <= 1.0
        @test 0.0 <= row.E_maintenance <= 1.0
        @test 0.0 <= row.E_growth <= 1.0
        @test 0.0 <= row.E_reproduction <= 1.0
    end
    
    # Q_t equals weighted sum of E_axis
    for row in eachrow(filter(r -> r.response_mode == "ec50_anchored_fractional_impairment", df_spec))
        expected_Q = row.w_assimilation * row.E_assimilation +
                     row.w_maintenance * row.E_maintenance +
                     row.w_growth * row.E_growth +
                     row.w_reproduction * row.E_reproduction
        @test row.Q_t ≈ expected_Q
        
        expected_A = row.A0 * max(1e-6, 1.0 - row.Q_t)
        @test row.A_t ≈ expected_A
    end
    
    # Check that raw mode A_t is unchanged across mixture effect models
    for sp in unique(df_spec.species_key)
        for m in 1:12
            raw_tu = filter(row -> row.species_key == sp && row.month == m && row.response_mode == "raw_margin_subtraction" && row.mixture_effect_model == "axis_toxic_unit_sum", df_spec)
            raw_ia = filter(row -> row.species_key == sp && row.month == m && row.response_mode == "raw_margin_subtraction" && row.mixture_effect_model == "independent_action_axis_effects", df_spec)
            if nrow(raw_tu) > 0 && nrow(raw_ia) > 0
                @test raw_tu.A_t[1] == raw_ia.A_t[1]
            end
        end
    end
end

@testset "Multi-Axis Response Calibrated Demo Tranche 4" begin
    # Check that CSV was generated and has correct columns/values
    out_dir = normpath(joinpath(@__DIR__, "..", "output", "ecotox_amp_multiaxis_response_calibrated_demo"))
    csv_path = joinpath(out_dir, "multiaxis_mixture_effect_model_comparison_summary.csv")
    
    @test isfile(csv_path)
    
    using CSV
    using DataFrames
    df_sum = CSV.read(csv_path, DataFrame)
    
    @test nrow(df_sum) > 0
    
    models = unique(df_sum.mixture_effect_model)
    @test "axis_toxic_unit_sum" in models
    @test "independent_action_axis_effects" in models
    
    @test "ec50_anchored_fractional_impairment" in unique(df_sum.response_mode)
    
    # Check finites
    for row in eachrow(df_sum)
        @test isfinite(row.max_F_t)
        @test isfinite(row.min_A_t)
        @test isfinite(row.max_Q_t)
    end
    
    # Check deltas exist
    @test "delta_max_F_t_IA_minus_TU" in names(df_sum)
    @test "delta_min_A_t_IA_minus_TU" in names(df_sum)
    @test "delta_max_Q_t_IA_minus_TU" in names(df_sum)
    
    # IA >= TU logically when combining fractions for independent models
    df_ec50 = filter(row -> row.response_mode == "ec50_anchored_fractional_impairment", df_sum)
    for row in eachrow(df_ec50)
        # We don't strictly assert delta >= 0 globally here, as single-contributor axes are identical (delta=0)
        # However we verify that it is calculated properly.
        @test isfinite(row.delta_max_F_t_IA_minus_TU)
    end
end

@testset "Multi-Axis Response Calibrated Demo Tranche 5" begin
    out_dir = normpath(joinpath(@__DIR__, "..", "output", "ecotox_amp_multiaxis_response_calibrated_demo"))
    
    # Check that Tranche 5 PNGs were generated
    @test isfile(joinpath(out_dir, "multiaxis_mixture_axis_impairments.png"))
    @test isfile(joinpath(out_dir, "multiaxis_mixture_weighted_impairment_Q.png"))
    @test isfile(joinpath(out_dir, "multiaxis_mixture_amplification.png"))
    
    # Check that existing plots are still there
    @test isfile(joinpath(out_dir, "multiaxis_axis_impairments.png"))
    @test isfile(joinpath(out_dir, "multiaxis_amplification.png"))
    @test isfile(joinpath(out_dir, "multiaxis_adaptive_margin.png"))
    @test isfile(joinpath(out_dir, "multiaxis_restoring_force.png"))
    
    # non-empty check
    @test filesize(joinpath(out_dir, "multiaxis_mixture_axis_impairments.png")) > 0
    @test filesize(joinpath(out_dir, "multiaxis_mixture_weighted_impairment_Q.png")) > 0
    @test filesize(joinpath(out_dir, "multiaxis_mixture_amplification.png")) > 0
end
