using Test
using TwoTimescaleResilience

@testset "Tranche 6: Metrics Extraction" begin
    # Test 6.1 -- trapezoid AUC linear function
    t = collect(0.0:0.1:1.0)
    y = copy(t)
    @test isapprox(trapezoid_auc(t, y), 0.5; atol=1e-10)

    # Test 6.2 -- peak time known vector
    t2 = [0.0, 1.0, 2.0, 3.0]
    y2 = [0.0, 2.0, 5.0, 1.0]

    sim_fake = TwoTimescaleSimulationResult(
        name="fake",
        background_stressors=BackgroundStressor[],
        pulse_stressors=PulseStressor[],
        B=0.0, A=1.0, lambda=1.0, amplification_factor=1.0,
        t=t2, P=zeros(4,0), D=zeros(4,0), M=y2, y=y2
    )

    metrics = compute_metrics(sim_fake)
    @test metrics.ymax == 5.0
    @test metrics.t_ymax == 2.0
    @test metrics.Mmax == 5.0
    @test metrics.t_Mmax == 2.0

    # Test 6.3 -- recovery time
    # Synthetic vector: pulse ends at t=2.0
    # Peak is 10. Threshold for recovery (fraction=0.05) is 0.5.
    t3 = collect(0.0:1.0:10.0)
    y3 = [0.0, 10.0, 8.0, 4.0, 2.0, 1.0, 0.4, 0.2, 0.1, 0.0, 0.0]
    p_fake = PulseStressor("fake", 10.0, 0.0, 2.0, 1.0, 1.0)

    sim_recovery = TwoTimescaleSimulationResult(
        name="rec",
        background_stressors=BackgroundStressor[],
        pulse_stressors=[p_fake],
        B=0.0, A=1.0, lambda=1.0, amplification_factor=1.0,
        t=t3, P=zeros(11,1), D=zeros(11,1), M=y3, y=y3
    )

    # Pulse ends at t=2.0. Looking after t=2.0.
    # At t=6.0, y=0.4 which is < 0.5.
    # recovery time = 6.0 - 2.0 = 4.0
    @test recovery_time_after_last_pulse(sim_recovery; fraction=0.05) == 4.0

    # Test 6.4 -- residual response
    @test metrics.residual_y == y2[end]
end
