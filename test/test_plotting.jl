using Test
using TwoTimescaleResilience

@testset "Tranche 7: CairoMakie Visualization" begin
    # Test 7.1 -- plot function creates file
    grid = [1.0 2.0; 3.0 4.0]
    filename = tempname() * ".png"

    plot_grid(grid; filename=filename)
    @test isfile(filename)
    @test filesize(filename) > 0

    # Test 7.2 -- plotted grid is not mutated
    grid_copy = copy(grid)
    plot_grid(grid; filename=tempname() * ".png")
    @test grid == grid_copy

    # Test 7.3 -- amplification grid plot
    F = [1.1 1.2; 1.5 2.0]
    filename_F = tempname() * ".png"
    plot_amplification_grid(F; filename=filename_F)
    @test isfile(filename_F)
    @test filesize(filename_F) > 0

    # Clean up test files
    rm(filename, force=true)
    rm(filename_F, force=true)
end
