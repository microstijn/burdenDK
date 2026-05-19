using Test
using TwoTimescaleResilience

# Simple helper trapezoid function for AUC testing
function simple_trapezoid_auc(t, y)
    dt = t[2] - t[1]
    return sum(y[1:end-1] .+ y[2:end]) / 2.0 * dt
end

@testset "Tranche 5: Response Dynamics" begin
    bg_params = BackgroundParams()
    bg_stressors_low = [BackgroundStressor("s1", 0.1, 1.0)]
    bg_stressors_high = [BackgroundStressor("s1", 1.0, 1.0)]

    # Test 5.1 -- zero pulse gives zero response
    p_zero = PulseStressor("zero", 0.0, 10.0, 20.0, 1.0, 1.0)
    sim_zero = simulate_two_timescale("zero", bg_stressors_low, [p_zero], bg_params)
    @test all(sim_zero.y .== 0.0)

    # Setup real pulse
    p_real = PulseStressor("p1", 10.0, 10.0, 30.0, 0.5, 1.0)

    sim_low = simulate_two_timescale("low", bg_stressors_low, [p_real], bg_params; tmax=200.0, dt=0.1)
    sim_high = simulate_two_timescale("high", bg_stressors_high, [p_real], bg_params; tmax=200.0, dt=0.1)

    # Test 5.2 -- stronger background gives lower lambda
    @test sim_high.B > sim_low.B
    @test sim_high.lambda <= sim_low.lambda

    # Test 5.3 -- same pulse gives larger AUC under high background
    auc_y_low = simple_trapezoid_auc(sim_low.t, sim_low.y)
    auc_y_high = simple_trapezoid_auc(sim_high.t, sim_high.y)
    @test auc_y_high >= auc_y_low

    # Test 5.4 -- analytical AUC identity
    # AUC_y ≈ q/lambda(B) * AUC_M
    q = 1.0
    auc_M_low = simple_trapezoid_auc(sim_low.t, sim_low.M)
    analytical_auc_y_low = (q / sim_low.lambda) * auc_M_low
    ratio = auc_y_low / analytical_auc_y_low
    @test 0.98 <= ratio <= 1.02

    # Test 5.5 -- amplification factor predicts AUC ratio
    sim_zero_bg = simulate_two_timescale("zero_bg", BackgroundStressor[], [p_real], bg_params; tmax=200.0, dt=0.1)
    auc_y_zero_bg = simple_trapezoid_auc(sim_zero_bg.t, sim_zero_bg.y)

    auc_ratio = auc_y_high / auc_y_zero_bg
    F_B = sim_high.amplification_factor

    @test isapprox(auc_ratio, F_B; rtol=1e-2)
end
