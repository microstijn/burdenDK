using Test
using TwoTimescaleResilience

@testset "Reduced DEB Response Simulation" begin
    deb_params = DEBAxisParams(
        alpha_axes=(0.30, 0.35, 0.20, 0.15),
        use_buffer_recovery_factor=true,
        beta_Z=0.5
    )

    # event_cost_from_axes tests
    # Format A: NamedTuple of vectors
    axes_nt_vecs = (
        assimilation=[1.0, 0.0],
        maintenance=[0.0, 1.0],
        growth=[0.0, 0.0],
        reproduction=[0.0, 0.0]
    )
    cost_A = event_cost_from_axes(axes_nt_vecs, deb_params)
    @test cost_A[1] ≈ 0.30
    @test cost_A[2] ≈ 0.35

    # Format B: Vector of NamedTuples
    axes_vec_nt = [
        (assimilation=1.0, maintenance=0.0, growth=0.0, reproduction=0.0),
        (assimilation=0.0, maintenance=1.0, growth=0.0, reproduction=0.0)
    ]
    cost_B = event_cost_from_axes(axes_vec_nt, deb_params)
    @test cost_B == cost_A

    # simulate_reduced_deb_response tests
    t = [0.0, 1.0, 2.0, 3.0]
    A_bg = 0.8
    axes_bg = (assimilation=0.1, maintenance=0.1, growth=0.0, reproduction=0.0)

    # 1. Zero C_event and y0=0 gives y all zeros.
    C_event_zero = [0.0, 0.0, 0.0, 0.0]
    res_zero = simulate_reduced_deb_response(t, C_event_zero, A_bg, axes_bg, deb_params)
    @test all(res_zero.y .== 0.0)

    # 2. Positive C_event gives nonnegative y.
    C_event_pos = [0.1, 0.1, 0.1, 0.1]
    res_pos = simulate_reduced_deb_response(t, C_event_pos, A_bg, axes_bg, deb_params)
    @test all(res_pos.y .>= 0.0)
    @test res_pos.y[end] > 0.0

    # 3. Same event under lower lambda gives larger AUC.
    # To get lower lambda, we can increase stress axes (or we can use negative Z)
    res_normal = simulate_reduced_deb_response(t, C_event_pos, A_bg, axes_bg, deb_params; q=1.0)

    A_bg_low_margin = 0.2 # this lowers lambda
    res_low_lambda = simulate_reduced_deb_response(t, C_event_pos, A_bg_low_margin, axes_bg, deb_params; q=1.0)

    auc_normal = sum(res_normal.y) # rough proxy for auc
    auc_low = sum(res_low_lambda.y)
    @test auc_low > auc_normal

    # 4. For constant C_event and constant lambda over long time, final y approximately equals q*C/lambda.
    t_long = collect(0.0:1.0:100.0)
    C_event_const = fill(0.1, length(t_long))
    res_long = simulate_reduced_deb_response(t_long, C_event_const, A_bg, axes_bg, deb_params; q=2.0)

    final_y = res_long.y[end]
    expected_y = 2.0 * 0.1 / res_long.lambda[end]
    @test final_y ≈ expected_y atol=1e-3

    # 5. Z_timeseries changes lambda when buffer recovery factor is enabled.
    Z_timeseries = fill(1.0, length(t)) # Positive Z increases lambda since beta_Z > 0
    res_with_Z = simulate_reduced_deb_response(t, C_event_pos, A_bg, axes_bg, deb_params; Z_timeseries=Z_timeseries)
    @test all(res_with_Z.lambda .> res_normal.lambda)
end
