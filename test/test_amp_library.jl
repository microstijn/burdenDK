using Test
using TwoTimescaleResilience

@testset "AmP Library Translation - Tranche 1" begin
    library_path = joinpath(@__DIR__, "..", "data", "AmP_Species_Library.json")

    @testset "Loading" begin
        library = load_amp_species_library(library_path)
        @test haskey(library, "Abatus_cordatus")
    end

    @testset "Species key normalization" begin
        @test amp_species_key("Abatus cordatus") == "Abatus_cordatus"
        @test amp_species_key("Abatus_cordatus") == "Abatus_cordatus"
        @test amp_species_key("  Abatus cordatus  ") == "Abatus_cordatus"
        @test amp_species_key("Abatus   cordatus") == "Abatus_cordatus"
    end

    @testset "Record validation" begin
        library = load_amp_species_library(library_path)
        @test validate_amp_record(library["Abatus_cordatus"]) == true
    end

    @testset "Conversion to DEBAxisParams" begin
        library = load_amp_species_library(library_path)
        params = amp_species_deb_params(library, "Abatus_cordatus")

        @test params isa DEBAxisParams
        @test params.A0 ≈ 1539.9871108901957
        @test params.alpha_axes isa NTuple{4, Float64}
        @test params.alpha_axes[1] ≈ 0.0006493560841700461
        @test params.alpha_axes[2] ≈ 0.4250074376301585
        @test params.alpha_axes[3] ≈ 0.77712
        @test params.alpha_axes[4] ≈ 0.22287999999999997
        @test params.lambda_min ≈ 0.005783592166791565   # = k_M (somatic maintenance rate constant)
        @test params.lambda_max ≈ 0.011568702452292915
    end

    @testset "Flexible species lookup" begin
        library = load_amp_species_library(library_path)
        params_key = amp_species_deb_params(library, "Abatus_cordatus")
        params_name = amp_species_deb_params(library, "Abatus cordatus")

        @test params_key.A0 ≈ params_name.A0
        @test params_key.alpha_axes == params_name.alpha_axes
        @test params_key.lambda_min ≈ params_name.lambda_min
        @test params_key.lambda_max ≈ params_name.lambda_max
    end

    @testset "Restoring force sanity check" begin
        library = load_amp_species_library(library_path)
        params = amp_species_deb_params(library, "Abatus_cordatus")
        # Linear recovery curve: at the pristine margin A0 the restoring force is
        # exactly lambda_max (no half-saturation compression).
        @test restoring_force_from_margin(params.A0, params) ≈ params.lambda_max
        @test restoring_force_from_margin(0.0, params) ≈ params.lambda_min
    end

    @testset "Invalid record tests" begin
        # Missing A0
        invalid_record_1 = Dict(
            "alpha_axes" => [0.1, 0.2, 0.3, 0.4],
            "lambda_bounds" => Dict("KA" => 10.0, "lambda_min" => 0.1, "lambda_max" => 1.0)
        )
        @test_throws ArgumentError validate_amp_record(invalid_record_1)

        # alpha_axes length != 4
        invalid_record_2 = Dict(
            "A0" => 100.0,
            "alpha_axes" => [0.1, 0.2, 0.3],
            "lambda_bounds" => Dict("KA" => 10.0, "lambda_min" => 0.1, "lambda_max" => 1.0)
        )
        @test_throws ArgumentError validate_amp_record(invalid_record_2)

        # missing lambda_min (KA is no longer required or validated)
        invalid_record_3 = Dict(
            "A0" => 100.0,
            "alpha_axes" => [0.1, 0.2, 0.3, 0.4],
            "lambda_bounds" => Dict("lambda_max" => 1.0)
        )
        @test_throws ArgumentError validate_amp_record(invalid_record_3)

        # lambda_max < lambda_min
        invalid_record_4 = Dict(
            "A0" => 100.0,
            "alpha_axes" => [0.1, 0.2, 0.3, 0.4],
            "lambda_bounds" => Dict("KA" => 10.0, "lambda_min" => 0.5, "lambda_max" => 0.1)
        )
        @test_throws ArgumentError validate_amp_record(invalid_record_4)

        # Test KeyError for missing species
        library = load_amp_species_library(library_path)
        @test_throws KeyError amp_species_deb_params(library, "Not_A_Real_Species")
    end
end
