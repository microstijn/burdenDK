using Test
using TwoTimescaleResilience

@testset "TwoTimescaleResilience Fast Tests" begin
    # -------------------------------------------------------------------------
    # Core mathematical invariants
    # -------------------------------------------------------------------------
    include("test_deb_axes.jl")
    include("test_deb_axis_response.jl")
    include("test_mixture_aggregation.jl")
    include("test_deb_axes_grid.jl")
    
    # -------------------------------------------------------------------------
    # Current development focus:
    # threshold-free vulnerability feature vectors + standardisation
    # -------------------------------------------------------------------------
    include("test_vulnerability_feature_vectors.jl")
    include("test_vulnerability_regime_clustering.jl")
    include("test_vulnerability_tranche_comparison.jl")

    include("test_vulnerability_regime_outputs.jl")

    # -------------------------------------------------------------------------
    # Lightweight empirical adapters / runtime basics
    # Keep these only if they are fast in practice.
    # -------------------------------------------------------------------------
    include("test_ecotox_library.jl")
    include("test_amp_library.jl")

    # -------------------------------------------------------------------------
    # Temporarily disabled: heavier / integration / IO / plotting / examples
    # Re-enable for full validation before merge or release.
    # -------------------------------------------------------------------------

    # include("test_deb_axes_grid.jl")
    # include("test_response_modes.jl")
    # include("test_default_mappings.jl")
    # include("test_deb_pipeline.jl")
    # include("test_pulse_deb_axes.jl")
    # include("test_examples.jl")
    # include("test_background.jl")
    # include("test_multistressor.jl")
    # include("test_grids.jl")
    # include("test_pulses.jl")
    # include("test_simulation.jl")
    # include("test_metrics.jl")
    # include("test_plotting.jl")
    # include("test_ascii.jl")
    # include("test_synthetic_raster.jl")
    # include("test_netcdf.jl")

    # include("test_recovery_penalty.jl")
    # include("test_condition_buffer.jl")
    # include("test_reduced_deb_response.jl")
    # include("test_mode_of_action.jl")
    # include("test_default_isimip_moa.jl")
    # include("test_exposure_filters.jl")
    # include("test_moa_deb_mapping.jl")
    # include("test_isimip_deb_pipeline.jl")
    # include("test_species_profiles.jl")
    # include("test_isimip_event_response.jl")
    # include("test_examples_isimip_moa.jl")
    # include("test_deb_pipeline_recovery_math.jl")

    # include("test_ecotox.jl")
    # include("test_amp_species_profile.jl")
    # include("test_amp_pipeline_integration.jl")
    # include("test_monthly_memory_demo.jl")

    # AmP species archetype database: likely data-building/heavy
    # include("test_amp_species_archetypes.jl")

    # Mixture-effect overlap demonstration example
    # include("test_mixture_effect_overlap_demo.jl")

    # Synthetic grid mixture demo
    # @testset "Synthetic Grid Mixture Demo" begin
    #     include("test_synthetic_grid_mixture_demo.jl")
    # end

    include("test_analytical_warm_start.jl")

    RUN_EXTENDED_TESTS = get(ENV, "TTR_RUN_EXTENDED_TESTS", "false") == "true"
    RUN_EXAMPLE_TESTS = get(ENV, "TTR_RUN_EXAMPLE_TESTS", "false") == "true"
    RUN_PLOTTING_TESTS = get(ENV, "TTR_RUN_PLOTTING_TESTS", "false") == "true"
    if RUN_EXTENDED_TESTS || RUN_EXAMPLE_TESTS
        include("test_archetype_compound_memory_10yr_grid_demo.jl")
    end
    if RUN_EXTENDED_TESTS || RUN_EXAMPLE_TESTS || RUN_PLOTTING_TESTS
        include("test_archetype_compound_memory_10yr_grid_plots.jl")
    end
end
# Multi-Tranche Grid Demo Test
include("test_archetype_compound_memory_multitranche_grid_demo.jl")
# Add DynQual helper tests if the flag is set
if get(ENV, "TTR_RUN_DYNQUAL_TESTS", "false") == "true" || get(ENV, "TTR_RUN_EXAMPLE_TESTS", "false") == "true"
    @info "Running DynQual example script tests..."
    include("test_dynqual_demo_helpers.jl")
else
    @info "Skipping DynQual helper tests. Set TTR_RUN_DYNQUAL_TESTS=true to run."
end
