using Test
using TwoTimescaleResilience

@testset "Tranche 4: Acute Pulse Stressors and TKTD Burdens" begin
    # Test 4.1 -- rectangular pulse
    amp = 5.0
    t_start = 2.0
    t_end = 5.0
    @test rectangular_pulse(1.9, amp, t_start, t_end) == 0.0
    @test rectangular_pulse(2.0, amp, t_start, t_end) == amp
    @test rectangular_pulse(3.0, amp, t_start, t_end) == amp
    @test rectangular_pulse(5.0, amp, t_start, t_end) == amp
    @test rectangular_pulse(5.1, amp, t_start, t_end) == 0.0

    # Test 4.2 -- constant exposure steady state
    p_steady = PulseStressor("steady", 10.0, 0.0, 10000.0, 0.1, 1.0)
    t_vec = collect(0.0:0.1:1000.0)
    dt = 0.1
    P = pulse_exposure_matrix(t_vec, [p_steady])
    D = burden_matrix_exact_update(t_vec, P, [p_steady], dt)
    # After sufficiently long time, D(t) -> P0
    @test isapprox(D[end, 1], 10.0; atol=1e-3)

    # Test 4.3 -- analytical pulse AUC
    # dD/dt = k(P-D)
    # integral_0^infty D(t)dt = P0 * Delta
    p_auc = PulseStressor("auc_test", 5.0, 10.0, 20.0, 0.5, 1.0)
    Delta = 20.0 - 10.0
    P0 = 5.0
    expected_auc = P0 * Delta

    dt = 0.01
    t_vec_auc = collect(0.0:dt:100.0)
    P_auc = pulse_exposure_matrix(t_vec_auc, [p_auc])
    D_auc = burden_matrix_exact_update(t_vec_auc, P_auc, [p_auc], dt)

    # Trapezoidal rule for AUC
    numeric_auc = sum(D_auc[1:end-1, 1] .+ D_auc[2:end, 1]) / 2.0 * dt
    @test isapprox(numeric_auc, expected_auc; rtol=1e-2)

    # Test 4.4 -- multi-pulse burden additivity
    p1 = PulseStressor("p1", 1.0, 0.0, 10.0, 1.0, 0.5)
    p2 = PulseStressor("p2", 2.0, 0.0, 10.0, 1.0, 0.25)
    Drow = [2.0, 4.0]
    # M(t) = w1 D1(t) + w2 D2(t) = 0.5(2.0) + 0.25(4.0) = 1.0 + 1.0 = 2.0
    @test acute_mixture_burden(Drow, [p1, p2]) == 2.0

    # Test 4.5 -- pulse interaction
    interaction = zeros(2, 2)
    interaction[1, 2] = 0.5
    # M(t) = 2.0 + eta * D1 * D2 = 2.0 + 0.5(2.0)(4.0) = 2.0 + 4.0 = 6.0
    @test acute_mixture_burden(Drow, [p1, p2]; interaction=interaction) == 6.0
end
