using TwoTimescaleResilience
using Test

@testset "DEB Amplification Pipeline" begin
    @testset "Test 4.1 -- pipeline consistency" begin
        L1 = fill(0.2, 2, 2)
        L2 = fill(0.4, 2, 2)
        layers = [L1, L2]

        mapping = default_pathogen_organic_deb_mapping(interaction_strength=0.25)
        params = DEBAxisParams(A0=1.0, alpha_axes=(1.0, 1.0, 1.0, 1.0))

        result = deb_amplification_pipeline(layers, mapping, params)

        # Manually compute
        sA = 0.2
        sM = 0.25
        sG = 0.1
        sR = 0.11

        A_expected = 1.0 - (1.0*sA + 1.0*sM + 1.0*sG + 1.0*sR)
        lambda_expected = restoring_force_from_margin(A_expected, params)
        F_expected = restoring_force_from_margin(params.A0, params) / lambda_expected

        @test all(result.axes.assimilation .≈ sA)
        @test all(result.axes.maintenance .≈ sM)
        @test all(result.axes.growth .≈ sG)
        @test all(result.axes.reproduction .≈ sR)
        @test all(result.A .≈ A_expected)
        @test all(result.lambda .≈ lambda_expected)
        @test all(result.amplification .≈ F_expected)
    end

    @testset "Test 4.2 -- amplification consistency" begin
        L1 = rand(3, 3)
        L2 = rand(3, 3)
        layers = [L1, L2]

        mapping = default_pathogen_organic_deb_mapping()
        params = DEBAxisParams(A0=2.0)

        result = deb_amplification_pipeline(layers, mapping, params)

        lambda_A0 = restoring_force_from_margin(params.A0, params)

        for r in 1:3
            for c in 1:3
                expected_F = lambda_A0 / result.lambda[r, c]
                @test isapprox(result.amplification[r, c], expected_F; atol=1e-8)
            end
        end
    end

    @testset "Test 4.3 -- no-stress baseline" begin
        L1 = zeros(2, 2)
        L2 = zeros(2, 2)
        layers = [L1, L2]

        mapping = default_pathogen_organic_deb_mapping()
        params = DEBAxisParams()

        result = deb_amplification_pipeline(layers, mapping, params)

        @test all(result.axes.assimilation .≈ 0.0)
        @test all(result.axes.maintenance .≈ 0.0)
        @test all(result.axes.growth .≈ 0.0)
        @test all(result.axes.reproduction .≈ 0.0)

        @test all(result.A .≈ params.A0)

        lambda_A0 = restoring_force_from_margin(params.A0, params)
        @test all(result.lambda .≈ lambda_A0)

        @test all(result.amplification .≈ 1.0)
    end

    @testset "Test 4.4 -- higher stress gives higher amplification" begin
        L1_low = fill(0.1, 2, 2)
        L2_low = fill(0.1, 2, 2)

        L1_high = fill(0.8, 2, 2)
        L2_high = fill(0.8, 2, 2)

        mapping = default_pathogen_organic_deb_mapping()
        params = DEBAxisParams()

        result_low = deb_amplification_pipeline([L1_low, L2_low], mapping, params)
        result_high = deb_amplification_pipeline([L1_high, L2_high], mapping, params)

        mean_amp_low = sum(result_low.amplification) / 4
        mean_amp_high = sum(result_high.amplification) / 4

        @test mean_amp_high >= mean_amp_low
    end
end
