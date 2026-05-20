using Test
using TwoTimescaleResilience

@testset "Condition Buffer Z" begin
    params = ConditionBufferParams(
        rho_A=0.20,
        rho_M=0.20,
        rho_E=0.10,
        delta_Z=0.05,
        omega_Z=0.25,
        Z_min=-1.0,
        Z_max=1.0
    )

    # 1. Default zero stress gives dZdt = rho_A - rho_M = 0. (assuming Z=0, C_event=0)
    axes_zero = (assimilation=0.0, maintenance=0.0, growth=0.0, reproduction=0.0)
    dZdt_zero = condition_buffer_derivative(0.0, axes_zero, 0.0, params)
    @test dZdt_zero ≈ 0.0

    # 2. Higher maintenance lowers dZdt.
    axes_maint = (assimilation=0.0, maintenance=0.5, growth=0.0, reproduction=0.0)
    dZdt_maint = condition_buffer_derivative(0.0, axes_maint, 0.0, params)
    @test dZdt_maint < dZdt_zero

    # 3. Positive event cost lowers dZdt.
    dZdt_cost = condition_buffer_derivative(0.0, axes_zero, 1.0, params)
    @test dZdt_cost < dZdt_zero

    # 4. update clamps to Z_min/Z_max.
    dt_large = 100.0
    # Make dZdt positive and large
    axes_good = (assimilation=0.0, maintenance=-0.5, growth=0.0, reproduction=0.0) # Unrealistic but tests clamping
    Z_high = update_condition_buffer(0.0, axes_good, -1.0, dt_large, params)
    @test Z_high ≈ params.Z_max

    # Make dZdt negative and large
    axes_bad = (assimilation=1.0, maintenance=1.0, growth=0.0, reproduction=0.0)
    Z_low = update_condition_buffer(0.0, axes_bad, 10.0, dt_large, params)
    @test Z_low ≈ params.Z_min

    # 5. simulate returns same length as t.
    t = [0.0, 1.0, 2.0, 3.0]
    axes_timeseries = [axes_zero, axes_zero, axes_zero, axes_zero]
    C_event_timeseries = [0.0, 0.1, 0.2, 0.0]
    Z_sim = simulate_condition_buffer(t, axes_timeseries, C_event_timeseries, params)
    @test length(Z_sim) == length(t)

    # 6. simulate errors if t is not increasing.
    t_bad = [0.0, 1.0, 1.0, 3.0]
    @test_throws ArgumentError simulate_condition_buffer(t_bad, axes_timeseries, C_event_timeseries, params)

    t_bad2 = [0.0, 2.0, 1.0, 3.0]
    @test_throws ArgumentError simulate_condition_buffer(t_bad2, axes_timeseries, C_event_timeseries, params)

    # 7. positive Z increases adaptive margin.
    deb_params = DEBAxisParams()
    A_zero_Z = adaptive_margin_with_buffer(axes_zero, 0.0, deb_params, params)
    A_pos_Z = adaptive_margin_with_buffer(axes_zero, 0.5, deb_params, params)
    @test A_pos_Z > A_zero_Z

    # 8. negative Z decreases adaptive margin.
    A_neg_Z = adaptive_margin_with_buffer(axes_zero, -0.5, deb_params, params)
    @test A_neg_Z < A_zero_Z
end
