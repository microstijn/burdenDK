using Test
using TwoTimescaleResilience
using NCDatasets

@testset "Tranche 10: NetCDF Loading Stub and Interface Design" begin
    # Test 10.1 -- normalisation minmax
    x = [0.0 5.0; 10.0 15.0]
    x_norm = normalise_layer(x)
    @test isapprox(x_norm[1, 1], 0.0)
    @test isapprox(x_norm[1, 2], 1/3)
    @test isapprox(x_norm[2, 1], 2/3)
    @test isapprox(x_norm[2, 2], 1.0)

    # Test 10.2 -- clipping
    # Without minmax, default logic if we applied it.
    # The normalise_layer func does minmax then clips.
    # So if original was 0..15, norm is 0..1.
    # We will test clipping on an already normalized array by providing bounds.
    x_clip = normalise_layer(x; lower=0.2, upper=0.8)
    @test isapprox(x_clip[1, 1], 0.0) # < 0.2
    @test isapprox(x_clip[1, 2], 1/3) # 0.333
    @test isapprox(x_clip[2, 1], 2/3) # 0.666
    @test isapprox(x_clip[2, 2], 1.0) # > 0.8

    # Test 10.3 -- missing values
    x_nan = [0.0 NaN; 10.0 10.0]
    x_nan_norm = normalise_layer(x_nan)
    @test isnan(x_nan_norm[1, 2])
    @test isapprox(x_nan_norm[1, 1], 0.0)
    @test isapprox(x_nan_norm[2, 1], 1.0)

    # Test 10.4 -- NetCDF tiny roundtrip
    filename = tempname() * ".nc"

    # Create tiny NetCDF
    NCDataset(filename, "c") do ds
        defDim(ds, "lon", 2)
        defDim(ds, "lat", 2)
        v = defVar(ds, "test_var", Float64, ("lon", "lat"))
        v[:, :] = [1.0 2.0; 3.0 4.0]
    end

    # Read it back using the extension function
    layer_read = load_nc_layer(filename, "test_var")
    @test size(layer_read) == (2, 2)
    @test layer_read[1, 1] == 1.0
    @test layer_read[2, 2] == 4.0

    rm(filename, force=true)
end
