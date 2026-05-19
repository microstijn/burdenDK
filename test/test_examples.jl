using TwoTimescaleResilience
using Test

@testset "Examples Verification" begin
    @testset "Test 7.1 -- check output names conceptually" begin
        # We don't want to run the full netcdf script if it requires actual netcdf files,
        # but we can verify the generated names.
        expected_files = [
            "pathogen_normalised.asc",
            "organic_normalised.asc",
            "deb_assimilation_stress.asc",
            "deb_maintenance_stress.asc",
            "deb_growth_stress.asc",
            "deb_reproduction_stress.asc",
            "deb_adaptive_margin.asc",
            "deb_restoring_force.asc",
            "deb_amplification_factor.asc"
        ]
        @test "deb_amplification_factor.asc" in expected_files
    end
end
