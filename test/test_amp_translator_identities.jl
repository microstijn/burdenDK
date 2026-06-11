using Test
using TwoTimescaleResilience

# ---------------------------------------------------------------------------
# Characterization (regression) tests pinning the AS-BUILT AmP -> capacity
# mapping produced by src/AmP_Translator.jl and consumed via src/amp_library.jl.
#
# These document the current behavior so the planned margin-equation refactor
# (see docs/claude/TwoTimescaleResilience_source_audit_2026-06-11.md) is
# auditable. They are EXPECTED to fail once the hidden `KA = 0.3*A0` knob is
# removed and the Fmax-kappa collapse is broken -- that failure is the signal
# that as-built behavior changed, and these assertions should be updated then.
#
# Identities pinned (derived in the audit note):
#   - kappa == alpha_axes[3]            (alpha_G = kap in AmP_Translator.jl)
#   - KA == 0.3 * A0                    (the hidden gain constant)
#   - lambda_max / lambda_min == 1/kappa
#   - Fmax := lambda(A0)/lambda(0) == 1 + (1/kappa - 1)/1.3
# ---------------------------------------------------------------------------

@testset "AmP_Translator as-built identities (characterization)" begin
    library_path = joinpath(@__DIR__, "..", "data", "AmP_Species_Library.json")
    library = load_amp_species_library(library_path)

    # ---- Library-wide identities across every shipped species ----
    @testset "Library-wide identities" begin
        checked = 0
        for (species_key, record) in library
            # Only consider records the adapter accepts.
            local params
            try
                params = amp_record_to_deb_params(record)
            catch
                continue
            end
            checked += 1

            kappa = params.alpha_axes[3]   # alpha_G == kap

            # KA is the hidden 0.3*A0 gain.
            @test params.KA ≈ 0.3 * params.A0

            # 1/kappa identity for the lambda bounds.
            @test params.lambda_max / params.lambda_min ≈ 1.0 / kappa

            # Fmax collapses to a pure function of kappa.
            Fmax_via_model = amplification_from_margin(0.0, params)          # lambda(A0)/lambda(0)
            Fmax_closed    = 1.0 + (1.0 / kappa - 1.0) / 1.3
            @test Fmax_via_model ≈ Fmax_closed

            # restoring force at A0 is lambda_min + (lambda_max-lambda_min)/1.3.
            lambda0_closed = params.lambda_min + (params.lambda_max - params.lambda_min) / 1.3
            @test restoring_force_from_margin(params.A0, params) ≈ lambda0_closed
        end
        @test checked > 0   # guard: the library actually yielded usable records
    end

    # ---- Pinned exact values for one species (Abatus_cordatus) ----
    @testset "Pinned: Abatus_cordatus" begin
        params = amp_species_deb_params(library, "Abatus_cordatus")
        kappa = params.alpha_axes[3]

        @test kappa ≈ 0.77712
        @test params.KA ≈ 0.3 * params.A0
        @test params.KA ≈ 461.99613326705867
        @test params.lambda_max / params.lambda_min ≈ 1.0 / 0.77712

        # Fmax for this species, closed form and pinned numeric value.
        @test amplification_from_margin(0.0, params) ≈ 1.0 + (1.0 / 0.77712 - 1.0) / 1.3
        @test amplification_from_margin(0.0, params) ≈ 1.22062 atol = 1e-4
    end
end
