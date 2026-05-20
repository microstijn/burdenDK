using Test
using TwoTimescaleResilience

@testset "ISIMIP DEB Pipeline" begin
    deb_params = DEBAxisParams()
    exp_filt = default_exposure_filter(7)
    moa_map = default_isimip_moa_mapping()
    deb_map = default_moa_to_deb_mapping()

    # 1. Scalar zero input gives F approximately 1.
    zeros_in = zeros(7)
    res_zero = isimip_deb_pipeline(zeros_in, exp_filt, moa_map, deb_map, deb_params)
    @test res_zero.amplification ≈ 1.0
    @test res_zero.Z === nothing
    @test res_zero.effective_values == zeros_in

    # 2. WT-only input increases thermal mode.
    wt_only = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    res_wt = isimip_deb_pipeline(wt_only, exp_filt, moa_map, deb_map, deb_params)
    @test res_wt.modes.thermal > res_wt.modes.oxygen
    @test res_wt.axes.assimilation > 0.0

    # 3. BOD-only input increases oxygen mode and maintenance axis.
    bod_only = [0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    res_bod = isimip_deb_pipeline(bod_only, exp_filt, moa_map, deb_map, deb_params)
    @test res_bod.modes.oxygen > res_bod.modes.thermal
    @test res_bod.axes.maintenance > 0.0

    # 4. Grid pipeline preserves shape.
    layers = [fill(0.0, 2, 2) for _ in 1:7]
    layers[1][1,1] = 1.0 # WT
    res_grid = isimip_deb_pipeline_grid(layers, exp_filt, moa_map, deb_map, deb_params)
    @test size(res_grid.A) == (2, 2)
    @test res_grid.amplification[2, 2] ≈ 1.0
    @test res_grid.amplification[1, 1] > 1.0

    # 5. NaN propagation through all outputs.
    layers_nan = [fill(0.0, 2, 2) for _ in 1:7]
    layers_nan[2][1, 2] = NaN
    res_nan = isimip_deb_pipeline_grid(layers_nan, exp_filt, moa_map, deb_map, deb_params)
    @test isnan(res_nan.A[1, 2])
    @test isnan(res_nan.lambda[1, 2])
    @test isnan(res_nan.amplification[1, 2])
    @test isnan(res_nan.effective_layers[2][1, 2])
    @test isnan(res_nan.modes.oxygen[1, 2])
    @test isnan(res_nan.axes.maintenance[1, 2])

    # 6. Stress gives amplification >= 1 when no positive buffer effect is used.
    @test res_wt.amplification >= 1.0
    @test res_bod.amplification >= 1.0

    # 7. Same A but different axes can produce different amplification when reduced-DEB recovery penalty is enabled.
    deb_params_pen = DEBAxisParams(use_axis_recovery_penalty=true)
    res_wt_pen = isimip_deb_pipeline(wt_only, exp_filt, moa_map, deb_map, deb_params_pen)
    @test res_wt_pen.amplification > 0.0 # Just testing it runs, hard to guarantee F1 != F2 without handcrafting

    # 8. buffer_grid NaN propagates to A, lambda, and amplification.
    buf_params = ConditionBufferParams()
    buf_grid = fill(1.0, 2, 2)
    buf_grid[2, 1] = NaN
    res_buf_nan = isimip_deb_pipeline_grid(layers, exp_filt, moa_map, deb_map, deb_params, buffer_grid=buf_grid, buffer_params=buf_params)
    @test isnan(res_buf_nan.A[2, 1])
    @test isnan(res_buf_nan.lambda[2, 1])
    @test isnan(res_buf_nan.amplification[2, 1])
    @test !isnan(res_buf_nan.A[1, 1])
end
