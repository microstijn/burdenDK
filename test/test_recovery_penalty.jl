using Test
using TwoTimescaleResilience

@testset "Recovery Penalty and Extended DEBAxisParams" begin
    # 1. With recovery penalty flags false, restoring_force_from_margin_and_axes equals restoring_force_from_margin.
    params_default = DEBAxisParams()
    A = 0.8
    axes_default = (assimilation=0.5, maintenance=0.2, growth=0.1, reproduction=0.0)
    @test restoring_force_from_margin_and_axes(A, axes_default, params_default) ≈ restoring_force_from_margin(A, params_default)

    # 2. With use_axis_recovery_penalty=true and positive axes, lambda is lower than or equal to the base lambda.
    params_penalty = DEBAxisParams(use_axis_recovery_penalty=true)
    lambda_base = restoring_force_from_margin(A, params_penalty)
    lambda_penalty = restoring_force_from_margin_and_axes(A, axes_default, params_penalty)
    @test lambda_penalty <= lambda_base

    # 3. Same A but different axis composition gives different lambda.
    axes1 = (assimilation=0.5, maintenance=0.0, growth=0.0, reproduction=0.0)
    axes2 = (assimilation=0.0, maintenance=0.5, growth=0.0, reproduction=0.0)
    params_diff_comp = DEBAxisParams(
        recovery_axes = (0.1, 0.8, 0.1, 0.05),
        use_axis_recovery_penalty = true
    )
    lambda1 = restoring_force_from_margin_and_axes(A, axes1, params_diff_comp)
    lambda2 = restoring_force_from_margin_and_axes(A, axes2, params_diff_comp)
    @test lambda2 < lambda1

    # 4. Positive Z with beta_Z > 0 increases lambda when use_buffer_recovery_factor=true.
    params_buffer = DEBAxisParams(
        use_buffer_recovery_factor=true,
        beta_Z=0.5
    )
    lambda_zero_Z = restoring_force_from_margin_and_axes(A, axes_default, params_buffer; Z=0.0)
    lambda_pos_Z = restoring_force_from_margin_and_axes(A, axes_default, params_buffer; Z=1.0)
    @test lambda_pos_Z > lambda_zero_Z

    # 5. Negative Z with beta_Z > 0 decreases lambda when use_buffer_recovery_factor=true.
    lambda_neg_Z = restoring_force_from_margin_and_axes(A, axes_default, params_buffer; Z=-1.0)
    @test lambda_neg_Z < lambda_zero_Z

    # 6. Lambda is always between lambda_min and lambda_max.
    params_clamp = DEBAxisParams(
        use_buffer_recovery_factor=true,
        beta_Z=10.0 # Extreme beta_Z to force lambda > lambda_max
    )
    lambda_extreme = restoring_force_from_margin_and_axes(A, axes_default, params_clamp; Z=10.0)
    @test lambda_extreme <= params_clamp.lambda_max

    params_clamp_low = DEBAxisParams(
        use_axis_recovery_penalty=true,
        recovery_axes=(10.0, 10.0, 10.0, 10.0) # Extreme penalty to force lambda < lambda_min
    )
    lambda_extreme_low = restoring_force_from_margin_and_axes(A, axes_default, params_clamp_low)
    @test lambda_extreme_low >= params_clamp_low.lambda_min

    # 7. Grid version preserves shape.
    Agrid = [0.8 0.5; 0.2 0.9]
    sA_grid = [0.1 0.0; 0.5 0.2]
    sM_grid = [0.2 0.1; 0.0 0.0]
    sG_grid = [0.0 0.0; 0.0 0.0]
    sR_grid = [0.0 0.0; 0.0 0.0]
    axes_grid = (assimilation=sA_grid, maintenance=sM_grid, growth=sG_grid, reproduction=sR_grid)

    lambdagrid = restoring_force_from_margin_and_axes_grid(Agrid, axes_grid, params_default)
    @test size(lambdagrid) == size(Agrid)

    # 8. Grid version propagates NaN.
    Agrid_nan = [0.8 NaN; 0.2 0.9]
    lambdagrid_nan = restoring_force_from_margin_and_axes_grid(Agrid_nan, axes_grid, params_default)
    @test isnan(lambdagrid_nan[1, 2])
    @test !isnan(lambdagrid_nan[1, 1])

    sA_grid_nan = [NaN 0.0; 0.5 0.2]
    axes_grid_nan = (assimilation=sA_grid_nan, maintenance=sM_grid, growth=sG_grid, reproduction=sR_grid)
    lambdagrid_nan2 = restoring_force_from_margin_and_axes_grid(Agrid, axes_grid_nan, params_default)
    @test isnan(lambdagrid_nan2[1, 1])

    Zgrid_nan = [1.0 1.0; 1.0 NaN]
    lambdagrid_nan3 = restoring_force_from_margin_and_axes_grid(Agrid, axes_grid, params_default; Zgrid=Zgrid_nan)
    @test isnan(lambdagrid_nan3[2, 2])
end
