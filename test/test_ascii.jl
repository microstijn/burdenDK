using Test
using TwoTimescaleResilience

@testset "Tranche 8: ASCII Grid Output" begin
    # Test 8.1 -- write and read identity
    grid = [1.0 2.0; 3.0 4.0]
    filename = tempname() * ".asc"

    write_ascii_grid(filename, grid)
    grid_read = read_ascii_grid(filename)

    @test grid == grid_read

    # Test 8.2 -- dimensions preserved
    @test size(grid_read) == (2, 2)

    # Test 8.3 -- nodata handling
    grid_with_nan = [1.0 NaN; 3.0 4.0]
    filename_nan = tempname() * ".asc"

    write_ascii_grid(filename_nan, grid_with_nan; nodata=-9999.0)
    grid_read_nan = read_ascii_grid(filename_nan)

    @test grid_read_nan[1, 1] == 1.0
    @test isnan(grid_read_nan[1, 2])
    @test grid_read_nan[2, 1] == 3.0
    @test grid_read_nan[2, 2] == 4.0

    # Clean up test files
    rm(filename, force=true)
    rm(filename_nan, force=true)
end
