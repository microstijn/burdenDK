export ConditionBufferParams
export condition_buffer_derivative, update_condition_buffer, simulate_condition_buffer
export adaptive_margin_with_buffer

Base.@kwdef struct ConditionBufferParams
    rho_A::Float64 = 0.20
    rho_M::Float64 = 0.20
    rho_E::Float64 = 0.10
    delta_Z::Float64 = 0.05
    omega_Z::Float64 = 0.25
    Z_min::Float64 = -1.0
    Z_max::Float64 = 1.0
end

function condition_buffer_derivative(Z, axes, C_event, params::ConditionBufferParams)
    s = _deb_axes_to_vector(axes)
    sA = s[1]
    sM = s[2]

    dZdt =
        params.rho_A * (1.0 - sA) -
        params.rho_M * (1.0 + sM) -
        params.rho_E * C_event -
        params.delta_Z * Z

    return dZdt
end

function update_condition_buffer(Z, axes, C_event, dt, params::ConditionBufferParams)
    Z_next = Z + dt * condition_buffer_derivative(Z, axes, C_event, params)
    return clamp(Z_next, params.Z_min, params.Z_max)
end

function simulate_condition_buffer(t, axes_timeseries, C_event_timeseries, params::ConditionBufferParams; Z0=0.0)
    ntime = length(t)
    if length(axes_timeseries) != ntime
        throw(ArgumentError("length(axes_timeseries) must equal length(t)"))
    end
    if length(C_event_timeseries) != ntime
        throw(ArgumentError("length(C_event_timeseries) must equal length(t)"))
    end

    for i in 1:(ntime - 1)
        if t[i+1] <= t[i]
            throw(ArgumentError("t must be strictly increasing"))
        end
    end

    Z = zeros(Float64, ntime)

    if ntime > 0
        Z[1] = Z0
    end

    for i in 1:(ntime - 1)
        dt = t[i+1] - t[i]

        # Support vector of tuples or similar indexable structure for axes
        axes_i = axes_timeseries[i]
        C_event_i = C_event_timeseries[i]

        Z[i+1] = update_condition_buffer(Z[i], axes_i, C_event_i, dt, params)
    end

    return Z
end

function adaptive_margin_with_buffer(axes, Z, deb_params::DEBAxisParams, buffer_params::ConditionBufferParams)
    s = _deb_axes_to_vector(axes)
    alpha = collect(deb_params.alpha_axes)
    return deb_params.A0 + buffer_params.omega_Z * Z - sum(alpha .* s)
end
