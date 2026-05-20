using Test
using TwoTimescaleResilience

@testset "DEB Amplification Pipeline Recovery Math" begin
    # Mock data
    layer1 = [0.1 0.2; 0.3 NaN]
    layer2 = [0.0 0.1; 0.2 0.0]
    layers = [layer1, layer2]

    mapping = DEBAxisMapping(
        W = [1.0 0.0; 0.0 1.0; 0.0 0.0; 0.0 0.0],
        clamp_nonnegative=true,
        clamp_unit=false
    )

    params_default = DEBAxisParams()

    # 1. Old behaviour unchanged when flags are false and no buffer is provided.
    res_old = deb_amplification_pipeline(layers, mapping, params_default)

    @test !haskey(res_old, :Z)
    @test res_old.A[1, 1] ≈ 1.0 - (0.3 * 0.1 + 0.35 * 0.0)
    @test isnan(res_old.amplification[2, 2])

    # 2. With recovery penalty enabled, lambda differs from simple restoring_force_from_margin_grid.
    params_penalty = DEBAxisParams(use_axis_recovery_penalty=true, recovery_axes=(0.1, 0.8, 0.1, 0.05))
    res_penalty = deb_amplification_pipeline(layers, mapping, params_penalty)

    # Simple lambda
    lambdagrid_simple = restoring_force_from_margin_grid(res_old.A, params_default)

    # In res_penalty, maintenance is > 0 in some cells, so lambda should be lower
    @test res_penalty.lambda[1, 2] < lambdagrid_simple[1, 2]

    # 3. Amplification increases when lambda decreases.
    @test res_penalty.amplification[1, 2] > res_old.amplification[1, 2]

    # 4. Positive buffer_grid can increase lambda if buffer recovery factor enabled.
    params_buffer = DEBAxisParams(use_buffer_recovery_factor=true, beta_Z=0.5)
    buffer_params = ConditionBufferParams()
    buffer_grid = [1.0 1.0; 1.0 1.0] # positive buffer

    res_buffer = deb_amplification_pipeline(layers, mapping, params_buffer; buffer_grid=buffer_grid, buffer_params=buffer_params)

    @test haskey(res_buffer, :Z)
    @test res_buffer.lambda[1, 1] > res_old.lambda[1, 1]

    # 5. NaN propagation still works.
    buffer_grid_nan = [1.0 1.0; NaN 1.0]
    res_nan = deb_amplification_pipeline(layers, mapping, params_buffer; buffer_grid=buffer_grid_nan, buffer_params=buffer_params)

    @test isnan(res_nan.Z[2, 1])
    @test isnan(res_nan.A[2, 1])
    @test isnan(res_nan.lambda[2, 1])
    @test isnan(res_nan.amplification[2, 1])

    # Layer NaN should also propagate
    @test isnan(res_nan.A[2, 2])
end
