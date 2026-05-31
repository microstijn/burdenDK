using Test
using TwoTimescaleResilience
using NCDatasets

module TestDynQual
    # Include the example script to test its functions
    # The script has a `if abspath(PROGRAM_FILE) == @__FILE__` guard, so it won't run.
    Base.include(@__MODULE__, joinpath(dirname(@__DIR__), "examples", "dynqual_synthetic_isimip_pressure_demo.jl"))
end

@testset "DynQual Example Script Local Helpers" begin
    # We will test the local helpers without running the full script.
    
    mktempdir() do temp_dir
        # Create a tiny dummy NetCDF file
        dummy_nc = joinpath(temp_dir, "dummy.nc")
        NCDataset(dummy_nc, "c") do ds
            defDim(ds, "lon", 2)
            defDim(ds, "lat", 2)
            defDim(ds, "time", 3)
            
            lon_var = defVar(ds, "lon", Float32, ("lon",))
            lat_var = defVar(ds, "lat", Float32, ("lat",))
            time_var = defVar(ds, "time", Int32, ("time",))
            
            main_var = defVar(ds, "organic_monthlyAvg_1980_2019", Float32, ("lon", "lat", "time"))
            
            lon_var[:] = [10.0, 11.0]
            lat_var[:] = [20.0, 21.0]
            time_var[:] = [1, 2, 3]
            
            main_var[:, :, 1] = [1.0 2.0; 3.0 NaN]
            main_var[:, :, 2] = [4.0 5.0; NaN 6.0]
            main_var[:, :, 3] = [7.0 8.0; 9.0 10.0]
        end
        
        # Test dataset opening and variable detection
        ds = TestDynQual.open_dynqual_dataset(dummy_nc)
        @test ds !== nothing
        
        lon_v, lat_v, time_v = TestDynQual.detect_lon_lat_time_vars(ds)
        @test lon_v == "lon"
        @test lat_v == "lat"
        @test time_v == "time"
        
        main_v = TestDynQual.detect_main_variable(ds, "organic_monthlyAvg_1980_2019")
        @test main_v == "organic_monthlyAvg_1980_2019"
        
        lons = ds[lon_v][:]
        lats = ds[lat_v][:]
        
        # Test subsetting (no stride)
        lon_idx, lat_idx = TestDynQual.subset_indices(lons, lats, (lon=(9.0, 12.0), lat=(19.0, 22.0)), 1)
        @test lon_idx == [1, 2]
        @test lat_idx == [1, 2]
        
        # Test read_month_slice
        slice_1 = TestDynQual.read_month_slice(ds, main_v, lon_idx, lat_idx, 1)
        @test size(slice_1) == (2, 2)
        @test slice_1[1, 1] ≈ 1.0f0
        @test isnan(slice_1[2, 2])
        
        close(ds)
        
        # Test estimate_log_quantiles_sampled
        p02, p98, nt, nf, ns, frac = TestDynQual.estimate_log_quantiles_sampled(
            dummy_nc, "organic_monthlyAvg_1980_2019", lon_idx, lat_idx; time_stride=1, max_samples=100
        )
        
        @test nt == 12  # 3 slices of 4 cells
        @test nf == 10  # 10 valid non-NaN cells
        @test ns == 10
        @test 0.0 <= frac <= 0.2
        @test p98 > p02
        
        # Test robust_log_scale_slice
        scaled = TestDynQual.robust_log_scale_slice(slice_1, p02, p98)
        @test size(scaled) == (2, 2)
        @test 0.0 <= scaled[1, 1] <= 1.0
        @test scaled[2, 2] == 0.0f0 # NaN became 0
    end
end
