using Test
using TwoTimescaleResilience

@testset "ExposureFilters" begin
    filter = ExposureFilter(multipliers=[0.5, 0.2])

    # 1. Vector multiplication works.
    vals = [10.0, 5.0]
    res = apply_exposure_filter(vals, filter)
    @test res[1] ≈ 5.0
    @test res[2] ≈ 1.0

    # 2. Grid multiplication preserves shape.
    layers = [fill(10.0, 2, 2), fill(5.0, 2, 2)]
    res_grid = apply_exposure_filter_grid(layers, filter)
    @test size(res_grid[1]) == (2, 2)
    @test size(res_grid[2]) == (2, 2)
    @test res_grid[1][1, 1] ≈ 5.0
    @test res_grid[2][2, 2] ≈ 1.0

    # 3. NaN propagation works.
    layers_nan = [fill(10.0, 2, 2), fill(5.0, 2, 2)]
    layers_nan[1][1, 2] = NaN
    res_grid_nan = apply_exposure_filter_grid(layers_nan, filter)
    @test isnan(res_grid_nan[1][1, 2])
    @test !isnan(res_grid_nan[1][1, 1])

    # NaN in second layer should leave the first layer intact for that cell
    layers_nan2 = [fill(10.0, 2, 2), fill(5.0, 2, 2)]
    layers_nan2[2][1, 2] = NaN
    res_grid_nan2 = apply_exposure_filter_grid(layers_nan2, filter)
    @test res_grid_nan2[1][1, 2] ≈ 5.0
    @test isnan(res_grid_nan2[2][1, 2])

    # 4. Length mismatch throws an error.
    @test_throws ArgumentError apply_exposure_filter([1.0], filter)
    @test_throws ArgumentError apply_exposure_filter_grid([layers[1]], filter)
    layers_mismatch = [fill(10.0, 2, 2), zeros(3, 3)]
    @test_throws ArgumentError apply_exposure_filter_grid(layers_mismatch, filter)

    # 5. Aquatic filter leaves values unchanged.
    aq_filter = aquatic_exposure_filter(2)
    res_aq = apply_exposure_filter(vals, aq_filter)
    @test res_aq == vals

    # 6. Human filter has length 7 and values in [0,1].
    hum_filter = human_exposure_filter()
    @test length(hum_filter.multipliers) == 7
    @test all(0.0 .<= hum_filter.multipliers .<= 1.0)
end
