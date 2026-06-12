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
# The recovery curve is now LINEAR in margin (no KA half-saturation constant), so
# at the pristine margin A0 the restoring force is exactly lambda_max and
#   Fmax == lambda(A0)/lambda_min == lambda_max/lambda_min == g.
#
# Identities pinned:
#   - lambda_min == min(k_M, lambda_max)            (k_M from auxiliary_metrics)
#   - unclamped species: lambda_max/lambda_min == g (energy investment ratio)
#                        and Fmax == g
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

            # Slow floor is the clamped maintenance rate constant.
            @test params.lambda_min ≈ min(k_M, params.lambda_max)

            if params.lambda_min < params.lambda_max - 1e-12
                # Unclamped: the timescale-separation ratio is the energy investment ratio g.
                @test params.lambda_max / params.lambda_min ≈ g rtol = 1e-6
                # Linear recovery curve: Fmax = lambda(A0)/lambda_min = lambda_max/lambda_min = g.
                @test amplification_from_margin(0.0, params) ≈ g rtol = 1e-6
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
        @test Float64(aux["g"]) ≈ 2.0002624871647248
        @test params.lambda_max / params.lambda_min ≈ 2.0002624871647248

        # Linear recovery curve: Fmax = lambda_max/lambda_min = g (pinned numeric).
        @test amplification_from_margin(0.0, params) ≈ 2.0002624871647248
    end
end
