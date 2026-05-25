module TwoTimescaleResilience

export load_nc_layer, normalise_layer, write_ascii_grid, read_ascii_grid, plot_scenario_comparison, plot_grid, plot_background_layers, plot_amplification_grid,
       trapezoid_auc, recovery_time_after_last_pulse, compute_metrics, TwoTimescaleSimulationResult, simulate_two_timescale, PulseStressor, rectangular_pulse, pulse_exposure_matrix, burden_matrix_exact_update, acute_mixture_burden,
       synthetic_background_layers, run_synthetic_raster_demo, compute_background_index_grid, adaptive_margin_grid, restoring_force_grid, amplification_factor_grid,
       positive_part, phi_background, adaptive_margin, restoring_force, amplification_factor, BackgroundParams, BackgroundStressor, background_index

export amp_species_key, load_amp_species_library, validate_amp_record, amp_record_to_deb_params, amp_species_deb_params, amp_species_profile

export load_ecotox_library, validate_ecotox_record, ecotox_active_stress
export ecotox_effect_to_deb_axis, deb_axis_index
export ecotox_record_to_deb_burden, ecotox_records_to_deb_burden, ecotox_records_to_deb_burden_stateful!
export ecotox_burden_to_response
export ecotox_filter_records, ecotox_records_for_taxon

export load_compound_memory_library, validate_compound_memory_record
export compound_retention, ecotox_default_retention

export EcotoxExposureState, get_internal_burden, set_internal_burden!, reset_internal_burdens!
export update_internal_burden!

include("deb_axes.jl")
include("condition_buffer.jl")
include("reduced_deb_response.jl")
include("mode_of_action.jl")
include("default_isimip_moa.jl")
include("exposure_filters.jl")
include("moa_deb_mapping.jl")
include("isimip_deb_pipeline.jl")
include("species_defaults.jl")
include("amp_library.jl")
include("isimip_event_response.jl")
include("ecotox_library.jl")
include("background.jl")
include("grids.jl")
include("pulses.jl")
include("simulation.jl")
include("metrics.jl")
include("plotting.jl")
include("ascii.jl")
include("netcdf.jl")

end # module
