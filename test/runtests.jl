using Test
using TwoTimescaleResilience

@testset "TwoTimescaleResilience Tests" begin
    include("test_background.jl")
    include("test_multistressor.jl")
    include("test_grids.jl")
    include("test_pulses.jl")
    include("test_simulation.jl")
    include("test_metrics.jl")
    include("test_plotting.jl")
    include("test_ascii.jl")
    include("test_synthetic_raster.jl")
    include("test_netcdf.jl")
end
