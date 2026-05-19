using Test
using TwoTimescaleResilience

@testset "Tranche 9: Synthetic Raster Demo" begin
    nrows, ncols = 20, 20
    layers = synthetic_background_layers(nrows, ncols)

    # Test 9.1 -- synthetic layers bounded
    for L in layers
        @test all(L .>= 0.0)
        @test all(L .<= 1.0)
    end

    # Run the demo locally
    output_dir = tempname()
    params = BackgroundParams()
    Bgrid, Agrid, lambdagrid, Fgrid = run_synthetic_raster_demo(params; output_dir=output_dir)

    # Test 9.2 -- amplification grid lower bound
    @test all(Fgrid .>= 1.0)

    # Test 9.3 -- high-stress cells amplify more
    max_idx = argmax(Bgrid)
    min_idx = argmin(Bgrid)
    @test Fgrid[max_idx] >= Fgrid[min_idx]

    # Test 9.4 -- files created
    expected_files = [
        "background_index_B.asc",
        "adaptive_margin_A.asc",
        "restoring_force_lambda.asc",
        "amplification_factor.asc",
        "background_index_B.png",
        "amplification_factor.png"
    ]

    for f in expected_files
        path = joinpath(output_dir, f)
        @test isfile(path)
        @test filesize(path) > 0
    end

    # Cleanup
    rm(output_dir, recursive=true, force=true)
end
