using Test
using TwoTimescaleResilience

@testset "Tranche 3: Raster/Matrix Background Operations" begin
    # Test 3.1 -- grid weighted sum
    L1 = fill(1.0, 2, 3)
    L2 = fill(2.0, 2, 3)
    weights = [0.25, 0.75]

    Bgrid = compute_background_index_grid([L1, L2], weights)
    @test size(Bgrid) == (2, 3)
    @test all(Bgrid .≈ 1.75)

    # Test 3.2 -- grid interaction
    interaction = zeros(2, 2)
    interaction[1, 2] = 2.0
    Bgrid_int = compute_background_index_grid([L1, L2], weights; interaction=interaction)
    @test all(Bgrid_int .≈ 5.75)

    # Test 3.3 -- grid shape preservation
    layers = [rand(7, 11), rand(7, 11)]
    w = [0.5, 0.5]
    Bgrid_shape = compute_background_index_grid(layers, w)
    @test size(Bgrid_shape) == (7, 11)

    params = BackgroundParams()
    Agrid = adaptive_margin_grid(Bgrid_shape, params)
    @test size(Agrid) == (7, 11)

    Rgrid = restoring_force_grid(Bgrid_shape, params)
    @test size(Rgrid) == (7, 11)

    Fgrid_shape = amplification_factor_grid(Bgrid_shape, params)
    @test size(Fgrid_shape) == (7, 11)

    # Test 3.4 -- amplification raster identity at zero background
    Bgrid_zero = zeros(5, 5)
    Fgrid_zero = amplification_factor_grid(Bgrid_zero, params)
    @test all(Fgrid_zero .== 1.0)

    # Test 3.5 -- amplification increases with burden
    B1 = fill(0.1, 4, 4)
    B2 = fill(0.8, 4, 4)
    F1 = amplification_factor_grid(B1, params)
    F2 = amplification_factor_grid(B2, params)
    @test all(F2 .>= F1)
end
