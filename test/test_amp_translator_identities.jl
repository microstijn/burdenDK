using Test
using TwoTimescaleResilience

# ---------------------------------------------------------------------------
# Characterization (regression) tests pinning the AS-BUILT AmP -> capacity
# mapping produced by src/AmP_Translator.jl and consumed via src/amp_library.jl.
#
# The slow recovery floor is the DEB somatic maintenance rate constant
#   lambda_min = min(k_M, lambda_max),   k_M = [p_M]/[E_G]
# (clamped so reserve-rich species, g < 1, get no timescale separation). This
# replaced the earlier lambda_min = p_M/A0, which forced lambda_max/lambda_min
# == 1/kappa and made amplification a function of kappa alone (the kappa-collapse).
# See docs/notes/lambda_min_maintenance_rate.tex.
#
# Identities pinned:
#   - KA == 0.3 * A0
#   - lambda_min == min(k_M, lambda_max)            (k_M from auxiliary_metrics)
#   - unclamped species: lambda_max/lambda_min == g (energy investment ratio)
#                        and Fmax == 1 + (g - 1)/1.3
#   - clamped species (g <= 1): lambda_min == lambda_max, Fmax == 1 (resilient)
# ---------------------------------------------------------------------------

@testset "AmP_Translator as-built identities (characterization)" begin
    library_path = joinpath(@__DIR__, "..", "data", "AmP_Species_Library.json")
    library = load_amp_species_library(library_path)

    # ---- Library-wide identities across every shipped species ----
    @testset "Library-wide identities" begin
        checked = 0
        clamped = 0
        for (_species, record) in library
            haskey(record, "auxiliary_metrics") || continue
            aux = record["auxiliary_metrics"]
            (haskey(aux, "k_M") && haskey(aux, "g")) || continue
            local params
            try
                params = amp_record_to_deb_params(record)
            catch
                continue
            end
            checked += 1

            k_M = Float64(aux["k_M"])
            g   = Float64(aux["g"])

            # KA is the (still-present) 0.3*A0 shape constant.
            @test params.KA ≈ 0.3 * params.A0

            # Slow floor is the clamped maintenance rate constant.
            @test params.lambda_min ≈ min(k_M, params.lambda_max)

            if params.lambda_min < params.lambda_max - 1e-12
                # Unclamped: the timescale-separation ratio is the energy investment ratio g.
                @test params.lambda_max / params.lambda_min ≈ g rtol = 1e-6
                @test amplification_from_margin(0.0, params) ≈ 1 + (g - 1) / 1.3 rtol = 1e-6
            else
                # Clamped (g <= 1): no separation, no amplification.
                clamped += 1
                @test amplification_from_margin(0.0, params) ≈ 1.0
            end
        end
        @test checked > 0
        @test clamped > 0   # a substantial fraction of AmP species are reserve-rich (g <= 1)
    end

    # ---- Pinned exact values for one species (Abatus_cordatus, unclamped, g ~ 2) ----
    @testset "Pinned: Abatus_cordatus" begin
        params = amp_species_deb_params(library, "Abatus_cordatus")
        aux = library["Abatus_cordatus"]["auxiliary_metrics"]

        @test params.lambda_min ≈ 0.005783592166791565
        @test params.lambda_min ≈ Float64(aux["k_M"])
        @test params.KA ≈ 0.3 * params.A0
        @test Float64(aux["g"]) ≈ 2.0002624871647248
        @test params.lambda_max / params.lambda_min ≈ 2.0002624871647248

        # Fmax = 1 + (g - 1)/1.3, closed form and pinned numeric value.
        @test amplification_from_margin(0.0, params) ≈ 1 + (2.0002624871647248 - 1) / 1.3
        @test amplification_from_margin(0.0, params) ≈ 1.7694326824344038 atol = 1e-10
    end
end
