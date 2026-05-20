export isimip_deb_pipeline
export isimip_deb_pipeline_grid

function isimip_deb_pipeline(values, exposure_filter, moa_mapping, moa_deb_mapping, deb_params;
                             buffer_params=nothing,
                             Z=nothing)

    # 1. Exposure filter
    effective_values = apply_exposure_filter(values, exposure_filter)

    # 2. Mode of action
    modes = mode_of_action(effective_values, moa_mapping)

    # 3. DEB axes
    axes = moa_to_deb_axes(modes, moa_deb_mapping)

    # 4. Adaptive margin
    if buffer_params !== nothing && Z !== nothing
        A = adaptive_margin_with_buffer(axes, Z, deb_params, buffer_params)
    else
        A = deb_adaptive_margin(axes, deb_params)
    end

    # 5. Restoring force
    if deb_params.use_axis_recovery_penalty || deb_params.use_buffer_recovery_factor
        lambda = restoring_force_from_margin_and_axes(A, axes, deb_params; Z=Z)
    else
        lambda = restoring_force_from_margin(A, deb_params)
    end

    # 6. Compute lambda_control from A0, zero axes, and Z=0
    zero_axes = (assimilation=0.0, maintenance=0.0, growth=0.0, reproduction=0.0)
    A_control = deb_params.A0
    if deb_params.use_axis_recovery_penalty || deb_params.use_buffer_recovery_factor
        lambda_control = restoring_force_from_margin_and_axes(A_control, zero_axes, deb_params; Z=0.0)
    else
        lambda_control = restoring_force_from_margin(A_control, deb_params)
    end

    # 7. Amplification factor
    F = lambda_control / lambda

    return (
        effective_values = effective_values,
        modes = modes,
        axes = axes,
        Z = Z,
        A = A,
        lambda = lambda,
        amplification = F
    )
end

function isimip_deb_pipeline_grid(layers, exposure_filter, moa_mapping, moa_deb_mapping, deb_params;
                                  buffer_grid=nothing,
                                  buffer_params=nothing)

    effective_layers = apply_exposure_filter_grid(layers, exposure_filter)
    modes = mode_of_action_grid(effective_layers, moa_mapping)
    axes = moa_to_deb_axes_grid(modes, moa_deb_mapping)

    dims = size(layers[1])

    # Pre-allocate output grids
    Agrid = Matrix{Float64}(undef, dims)
    lambdagrid = Matrix{Float64}(undef, dims)
    Fgrid = Matrix{Float64}(undef, dims)

    # Control restoring force (same for all cells)
    zero_axes = (assimilation=0.0, maintenance=0.0, growth=0.0, reproduction=0.0)
    A_control = deb_params.A0
    if deb_params.use_axis_recovery_penalty || deb_params.use_buffer_recovery_factor
        lambda_control = restoring_force_from_margin_and_axes(A_control, zero_axes, deb_params; Z=0.0)
    else
        lambda_control = restoring_force_from_margin(A_control, deb_params)
    end

    use_buffer = buffer_grid !== nothing && buffer_params !== nothing
    use_recovery = deb_params.use_axis_recovery_penalty || deb_params.use_buffer_recovery_factor

    for c in 1:dims[2]
        for r in 1:dims[1]
            has_nan = false
            for i in 1:length(layers)
                if isnan(layers[i][r, c])
                    has_nan = true
                    break
                end
            end

            if use_buffer && isnan(buffer_grid[r, c])
                has_nan = true
            end

            if has_nan
                Agrid[r, c] = NaN
                lambdagrid[r, c] = NaN
                Fgrid[r, c] = NaN
                continue
            end

            cell_axes = (
                assimilation = axes.assimilation[r, c],
                maintenance = axes.maintenance[r, c],
                growth = axes.growth[r, c],
                reproduction = axes.reproduction[r, c]
            )

            Z_cell = use_buffer ? buffer_grid[r, c] : nothing

            if use_buffer
                A_cell = adaptive_margin_with_buffer(cell_axes, Z_cell, deb_params, buffer_params)
            else
                A_cell = deb_adaptive_margin(cell_axes, deb_params)
            end

            if use_recovery
                lambda_cell = restoring_force_from_margin_and_axes(A_cell, cell_axes, deb_params; Z=Z_cell)
            else
                lambda_cell = restoring_force_from_margin(A_cell, deb_params)
            end

            Agrid[r, c] = A_cell
            lambdagrid[r, c] = lambda_cell
            Fgrid[r, c] = lambda_control / lambda_cell
        end
    end

    return (
        effective_layers = effective_layers,
        modes = modes,
        axes = axes,
        Z = buffer_grid,
        A = Agrid,
        lambda = lambdagrid,
        amplification = Fgrid
    )
end
