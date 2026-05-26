using Test
using TwoTimescaleResilience

@testset "Response Modes Tests" begin
    @testset "EC50 Anchored Fractional Impairment" begin
        # X = 0 gives E = 0
        zero_input = (assimilation=0.0, maintenance=0.0, growth=0.0, reproduction=0.0)
        zero_res = ec50_anchored_fractional_impairment(zero_input)
        @test zero_res.assimilation == 0.0
        @test zero_res.maintenance == 0.0
        @test zero_res.growth == 0.0
        @test zero_res.reproduction == 0.0

        # X = 1 gives E = 0.5
        one_input = (assimilation=1.0, maintenance=1.0, growth=1.0, reproduction=1.0)
        one_res = ec50_anchored_fractional_impairment(one_input)
        @test one_res.assimilation == 0.5
        @test one_res.maintenance == 0.5
        @test one_res.growth == 0.5
        @test one_res.reproduction == 0.5

        # X = 999 gives E ≈ 0.999
        large_input = [999.0, 999.0, 999.0, 999.0]
        large_res = ec50_anchored_fractional_impairment(large_input)
        @test isapprox(large_res.assimilation, 0.999, atol=1e-5)

        # negative X throws ArgumentError
        @test_throws ArgumentError ec50_anchored_fractional_impairment([-1.0, 0.0, 0.0, 0.0])

        # NaN throws ArgumentError
        @test_throws ArgumentError ec50_anchored_fractional_impairment([NaN, 0.0, 0.0, 0.0])

        # Inf throws ArgumentError
        @test_throws ArgumentError ec50_anchored_fractional_impairment([Inf, 0.0, 0.0, 0.0])
    end
end

    @testset "Species-Specific Axis Weights" begin
        # Default valid DEBAxisParams
        params = DEBAxisParams(alpha_axes=(0.30, 0.35, 0.20, 0.15))
        
        # Test "auto" method
        res_auto = axis_weights_for_species(params)
        @test res_auto.axis_weight_method == "normalized_alpha_axes"
        @test isapprox(res_auto.w_assimilation, 0.30)
        @test isapprox(res_auto.w_maintenance, 0.35)
        @test isapprox(res_auto.w_growth, 0.20)
        @test isapprox(res_auto.w_reproduction, 0.15)
        sum_weights = res_auto.w_assimilation + res_auto.w_maintenance + res_auto.w_growth + res_auto.w_reproduction
        @test isapprox(sum_weights, 1.0)
        
        # Test valid method="normalized_alpha_axes"
        res_norm = axis_weights_for_species(params, method="normalized_alpha_axes")
        @test res_norm.axis_weight_method == "normalized_alpha_axes"
        @test isapprox(res_norm.w_assimilation + res_norm.w_maintenance + res_norm.w_growth + res_norm.w_reproduction, 1.0)
        
        # Fallback explicitly requested
        res_fallback = axis_weights_for_species(params, method="equal_weight_diagnostic_fallback")
        @test res_fallback.axis_weight_method == "equal_weight_diagnostic_fallback"
        @test res_fallback.w_assimilation == 0.25
        @test res_fallback.w_maintenance == 0.25
        @test res_fallback.w_growth == 0.25
        @test res_fallback.w_reproduction == 0.25
        
        # Auto method with invalid DEBAxisParams (negative alpha)
        invalid_params_1 = DEBAxisParams(alpha_axes=(-0.1, 0.5, 0.5, 0.1))
        res_invalid_1 = axis_weights_for_species(invalid_params_1)
        @test res_invalid_1.axis_weight_method == "equal_weight_diagnostic_fallback"
        
        # Auto method with invalid DEBAxisParams (sum <= 0)
        invalid_params_2 = DEBAxisParams(alpha_axes=(0.0, 0.0, 0.0, 0.0))
        res_invalid_2 = axis_weights_for_species(invalid_params_2)
        @test res_invalid_2.axis_weight_method == "equal_weight_diagnostic_fallback"
        
        # normalized_alpha_axes method with invalid DEBAxisParams throws ArgumentError
        @test_throws ArgumentError axis_weights_for_species(invalid_params_1, method="normalized_alpha_axes")
        @test_throws ArgumentError axis_weights_for_species(invalid_params_2, method="normalized_alpha_axes")
        
        # Unsupported method throws ArgumentError
        @test_throws ArgumentError axis_weights_for_species(params, method="unknown_method")
    end

    @testset "Runtime Response Function" begin
        params = DEBAxisParams(A0=1.0, alpha_axes=(0.30, 0.35, 0.20, 0.15))
        zero_input = (assimilation=0.0, maintenance=0.0, growth=0.0, reproduction=0.0)
        pressure_input = (assimilation=1.0, maintenance=0.5, growth=0.0, reproduction=0.0)
        
        # Test raw margin subtraction (zero pressure)
        res_raw_zero = compute_adaptive_margin_response(zero_input, params, response_mode="raw_margin_subtraction")
        @test res_raw_zero.A_t == deb_adaptive_margin(zero_input, params)
        @test res_raw_zero.lambda_t == restoring_force_from_margin(res_raw_zero.A_t, params)
        @test res_raw_zero.F_t == amplification_from_margin(res_raw_zero.A_t, params)
        
        # Test raw margin subtraction (with pressure)
        res_raw_pressure = compute_adaptive_margin_response(pressure_input, params, response_mode="raw_margin_subtraction")
        @test res_raw_pressure.A_t == deb_adaptive_margin(pressure_input, params)
        @test res_raw_pressure.lambda_t == restoring_force_from_margin(res_raw_pressure.A_t, params)
        @test res_raw_pressure.F_t == amplification_from_margin(res_raw_pressure.A_t, params)
        
        # Test EC50 anchored fractional impairment (zero pressure)
        res_ec50_zero = compute_adaptive_margin_response(zero_input, params, response_mode="ec50_anchored_fractional_impairment")
        @test res_ec50_zero.A_t == params.A0
        @test res_ec50_zero.Q_t == 0.0
        
        # Test EC50 anchored fractional impairment (with pressure)
        res_ec50_pressure = compute_adaptive_margin_response(pressure_input, params, response_mode="ec50_anchored_fractional_impairment")
        @test res_ec50_pressure.A_t < params.A0
        @test res_ec50_pressure.Q_t > 0.0 && res_ec50_pressure.Q_t < 1.0
        
        # Check audit fields
        @test res_ec50_pressure.E_assimilation == 0.5
        @test isapprox(res_ec50_pressure.E_maintenance, 0.5 / 1.5)
        
        # Check aliases
        @test res_ec50_pressure.A == res_ec50_pressure.A_t
        @test res_ec50_pressure.lambda == res_ec50_pressure.lambda_t
        @test res_ec50_pressure.amplification == res_ec50_pressure.F_t
        
        # Validation
        @test_throws ArgumentError compute_adaptive_margin_response(zero_input, params, response_mode="unknown_mode")
        @test_throws ArgumentError compute_adaptive_margin_response(zero_input, params, A_floor_fraction=0.0)
        @test_throws ArgumentError compute_adaptive_margin_response(zero_input, params, A_floor_fraction=1.5)
        @test_throws ArgumentError compute_adaptive_margin_response(zero_input, params, A_floor_fraction=-0.1)
    end

    @testset "Multiaxis Example Outputs" begin
        # Test outputs from tranche 4
        out_dir = normpath(joinpath(@__DIR__, "..", "output", "ecotox_amp_multiaxis_response_calibrated_demo"))
        
        # Test existence
        @test isfile(joinpath(out_dir, "multiaxis_species_summary.csv"))
        @test isfile(joinpath(out_dir, "multiaxis_compound_summary.csv"))
        
        # We can read and verify basic structure
        using CSV, DataFrames
        df_spec = CSV.read(joinpath(out_dir, "multiaxis_species_summary.csv"), DataFrame)
        df_comp = CSV.read(joinpath(out_dir, "multiaxis_compound_summary.csv"), DataFrame)
        
        # Dimensions
        @test nrow(df_spec) == 144
        @test nrow(df_comp) == 432
        
        # Check required columns
        @test "response_mode" in names(df_spec)
        @test "response_mode" in names(df_comp)
        
        # Check response modes
        @test "raw_margin_subtraction" in df_spec.response_mode
        @test "ec50_anchored_fractional_impairment" in df_spec.response_mode
        
        # Verify mixture methods
        @test "additive_axis_burden" in df_spec.mixture_method
        @test "axis_toxic_unit_sum" in df_spec.mixture_method
        
        # Test invariants for EC50 vs Raw
        raw_df = filter(row -> row.response_mode == "raw_margin_subtraction", df_spec)
        ec50_df = filter(row -> row.response_mode == "ec50_anchored_fractional_impairment", df_spec)
        
        # EC50 anchored fractional impairment should result in larger F_t when under pressure
        for i in 1:nrow(raw_df)
            if raw_df.X_maintenance[i] > 0 || raw_df.X_growth[i] > 0 || raw_df.X_reproduction[i] > 0
                @test ec50_df.F_t[i] >= raw_df.F_t[i]
            end
        end
    end

    @testset "Check exact columns in df_spec" begin
        out_dir = normpath(joinpath(@__DIR__, "..", "output", "ecotox_amp_multiaxis_response_calibrated_demo"))
        using CSV, DataFrames
        df_spec = CSV.read(joinpath(out_dir, "multiaxis_species_summary.csv"), DataFrame)
        
        required_cols = [
            "response_mode", "axis_weight_method", "axis_weight_scope",
            "X_assimilation", "X_maintenance", "X_growth", "X_reproduction",
            "E_assimilation", "E_maintenance", "E_growth", "E_reproduction",
            "w_assimilation", "w_maintenance", "w_growth", "w_reproduction",
            "Q_t", "A0", "A_t", "lambda0", "lambda_t", "F_t"
        ]
        
        for c in required_cols
            @test c in names(df_spec)
        end
        
        # Test X_axis equals total_burden_axis
        for i in 1:nrow(df_spec)
            @test isapprox(df_spec.X_assimilation[i], df_spec.total_burden_assimilation[i])
            @test isapprox(df_spec.X_maintenance[i], df_spec.total_burden_maintenance[i])
            @test isapprox(df_spec.X_growth[i], df_spec.total_burden_growth[i])
            @test isapprox(df_spec.X_reproduction[i], df_spec.total_burden_reproduction[i])
            
            # E_axis equals X_axis / (1 + X_axis)
            @test isapprox(df_spec.E_assimilation[i], df_spec.X_assimilation[i] / (1 + df_spec.X_assimilation[i]))
            
            # weights sum to 1
            w_sum = df_spec.w_assimilation[i] + df_spec.w_maintenance[i] + df_spec.w_growth[i] + df_spec.w_reproduction[i]
            @test isapprox(w_sum, 1.0)
            
            # Q_t is weighted sum
            q_expected = df_spec.w_assimilation[i]*df_spec.E_assimilation[i] +
                         df_spec.w_maintenance[i]*df_spec.E_maintenance[i] +
                         df_spec.w_growth[i]*df_spec.E_growth[i] +
                         df_spec.w_reproduction[i]*df_spec.E_reproduction[i]
            @test isapprox(df_spec.Q_t[i], q_expected)
        end
    end

    @testset "Response Mode Comparison Summary" begin
        out_dir = normpath(joinpath(@__DIR__, "..", "output", "ecotox_amp_multiaxis_response_calibrated_demo"))
        using CSV, DataFrames
        comp_df = CSV.read(joinpath(out_dir, "multiaxis_response_mode_comparison_summary.csv"), DataFrame)
        
        # We expect 2 response modes * 3 species = 6 rows
        @test nrow(comp_df) == 6
        
        # Test required columns
        req_cols = [
            "species_key", "response_mode", "max_F_t", "month_max_F_t",
            "min_A_t", "month_min_A_t", "max_Q_t", "month_max_Q_t",
            "activated_axes", "axis_weight_method", "axis_weight_scope",
            "delta_max_F_t_ec50_minus_raw", "delta_min_A_t_ec50_minus_raw"
        ]
        
        for c in req_cols
            @test c in names(comp_df)
        end
        
        # Test specific properties
        @test "raw_margin_subtraction" in comp_df.response_mode
        @test "ec50_anchored_fractional_impairment" in comp_df.response_mode
        
        for r in eachrow(comp_df)
            @test isfinite(r.max_F_t)
            @test isfinite(r.min_A_t)
            @test isfinite(r.max_Q_t)
            @test isfinite(r.delta_max_F_t_ec50_minus_raw)
            @test isfinite(r.delta_min_A_t_ec50_minus_raw)
            @test r.axis_weight_method == "normalized_alpha_axes" || r.axis_weight_method == "equal_weight_diagnostic_fallback"
            @test r.axis_weight_scope == "all_axes"
            
            if r.response_mode == "ec50_anchored_fractional_impairment"
                @test r.delta_max_F_t_ec50_minus_raw >= 0.0
                @test r.delta_min_A_t_ec50_minus_raw <= 0.0
            else
                @test isapprox(r.delta_max_F_t_ec50_minus_raw, 0.0)
                @test isapprox(r.delta_min_A_t_ec50_minus_raw, 0.0)
            end
        end
    end
