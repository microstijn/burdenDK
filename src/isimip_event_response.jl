export event_burdens_to_modes
export event_burdens_to_deb_axes
export simulate_isimip_deb_event_response

function event_burdens_to_modes(D::Matrix{Float64}, moa_mapping::ModeOfActionMapping)
    n_time, n_vars = size(D)

    if n_vars != size(moa_mapping.U, 2)
        throw(ArgumentError("Number of columns in D ($n_vars) must match columns in moa_mapping.U ($(size(moa_mapping.U, 2)))"))
    end

    n_modes = size(moa_mapping.U, 1)

    modes_series = [Vector{Float64}(undef, n_time) for _ in 1:n_modes]

    val_buf = zeros(Float64, n_vars)
    for t in 1:n_time
        for i in 1:n_vars
            val_buf[i] = D[t, i]
        end
        mv = mode_of_action_vector(val_buf, moa_mapping)
        for m in 1:n_modes
            modes_series[m][t] = mv[m]
        end
    end

    if n_modes == 8
        return (
            thermal = modes_series[1],
            oxygen = modes_series[2],
            osmotic = modes_series[3],
            immune = modes_series[4],
            eutrophication = modes_series[5],
            toxic = modes_series[6],
            feeding = modes_series[7],
            physical = modes_series[8]
        )
    else
        return NamedTuple{Tuple(Symbol("mode_$i") for i in 1:n_modes)}(Tuple(modes_series))
    end
end

function event_burdens_to_deb_axes(mode_timeseries, moa_deb_mapping::MoAToDEBMapping)
    return moa_to_deb_axes_timeseries(mode_timeseries, moa_deb_mapping)
end

function simulate_isimip_deb_event_response(t, D, moa_mapping, moa_deb_mapping,
                                            A_background, axes_background, deb_params;
                                            q=1.0,
                                            y0=0.0,
                                            Z_timeseries=nothing)

    # 1. Map event burdens to modes
    modes_event = event_burdens_to_modes(D, moa_mapping)

    # 2. Map modes to axes
    axes_event = event_burdens_to_deb_axes(modes_event, moa_deb_mapping)

    # 3. Compute event cost from axes (using the function from Prompt A)
    C_event = event_cost_from_axes(axes_event, deb_params)

    # 4. Simulate response
    response = simulate_reduced_deb_response(
        t,
        C_event,
        A_background,
        axes_background,
        deb_params;
        q=q,
        y0=y0,
        Z_timeseries=Z_timeseries
    )

    return (
        modes_event = modes_event,
        axes_event = axes_event,
        C_event = C_event,
        response = response
    )
end
