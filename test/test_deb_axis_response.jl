using TwoTimescaleResilience
using Test

@testset "DEB-Like Response Simulation" begin
    @testset "Test 6.1 -- zero pulse cost gives zero response" begin
        ntime = 10
        t = collect(range(0.0, step=0.1, length=ntime))

        pulse_axes = (
            assimilation = zeros(ntime),
            maintenance = zeros(ntime),
            growth = zeros(ntime),
            reproduction = zeros(ntime)
        )

        params = DEBAxisParams()
        A_bg = 1.0
        q = 1.0
        dt = 0.1

        y = simulate_deb_axis_response(t, pulse_axes, A_bg, params, q; y0=0.0, dt=dt)

        @test all(y .== 0.0)
    end

    @testset "Test 6.2 -- positive pulse cost gives nonnegative response" begin
        ntime = 10
        t = collect(range(0.0, step=0.1, length=ntime))

        pulse_axes = (
            assimilation = rand(ntime),
            maintenance = rand(ntime),
            growth = zeros(ntime),
            reproduction = zeros(ntime)
        )

        params = DEBAxisParams()
        A_bg = 1.0
        q = 1.0
        dt = 0.1

        y = simulate_deb_axis_response(t, pulse_axes, A_bg, params, q; y0=0.0, dt=dt)

        @test all(y .>= 0.0)
    end

    @testset "Test 6.3 -- stronger background margin depletion increases AUC" begin
        ntime = 50
        t = collect(range(0.0, step=0.1, length=ntime))

        # constant positive pulse
        pulse_axes = (
            assimilation = fill(0.5, ntime),
            maintenance = fill(0.0, ntime),
            growth = fill(0.0, ntime),
            reproduction = fill(0.0, ntime)
        )

        params = DEBAxisParams()
        q = 1.0
        dt = 0.1

        A_bg_high = 1.5
        A_bg_low = 0.2

        y_high_bg = simulate_deb_axis_response(t, pulse_axes, A_bg_high, params, q; y0=0.0, dt=dt)
        y_low_bg = simulate_deb_axis_response(t, pulse_axes, A_bg_low, params, q; y0=0.0, dt=dt)

        auc_high = sum(y_high_bg) * dt
        auc_low = sum(y_low_bg) * dt

        # lower margin A_bg_low => lower lambda => higher response y => higher AUC
        @test auc_low >= auc_high
    end

    @testset "Test 6.4 -- constant cost steady-state" begin
        ntime = 500
        dt = 0.5
        t = collect(range(0.0, step=dt, length=ntime))

        sA_const = 0.5
        pulse_axes = (
            assimilation = fill(sA_const, ntime),
            maintenance = zeros(ntime),
            growth = zeros(ntime),
            reproduction = zeros(ntime)
        )

        params = DEBAxisParams()
        A_bg = 1.0
        q = 2.0

        y = simulate_deb_axis_response(t, pulse_axes, A_bg, params, q; y0=0.0, dt=dt)

        alphaA = params.alpha_axes[1]
        cost_const = alphaA * sA_const
        A_total_const = A_bg - cost_const
        lambda_const = restoring_force_from_margin(A_total_const, params)

        y_star = q * cost_const / lambda_const

        @test isapprox(y[end], y_star; atol=1e-3)
    end
end
