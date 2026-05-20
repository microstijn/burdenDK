export event_cost_from_axes
export simulate_reduced_deb_response

function event_cost_from_axes(axes_timeseries, deb_params::DEBAxisParams)
    alpha = collect(deb_params.alpha_axes)

    if axes_timeseries isa NamedTuple
        # NamedTuple of vectors
        sA = axes_timeseries.assimilation
        sM = axes_timeseries.maintenance
        sG = axes_timeseries.growth
        sR = axes_timeseries.reproduction

        ntime = length(sA)
        C_event = zeros(Float64, ntime)

        for i in 1:ntime
            axes_i = [sA[i], sM[i], sG[i], sR[i]]
            C_event[i] = sum(alpha .* axes_i)
        end
        return C_event
    elseif axes_timeseries isa AbstractVector
        # Vector of NamedTuples
        ntime = length(axes_timeseries)
        C_event = zeros(Float64, ntime)

        for i in 1:ntime
            axes_i = _deb_axes_to_vector(axes_timeseries[i])
            C_event[i] = sum(alpha .* axes_i)
        end
        return C_event
    else
        throw(ArgumentError("axes_timeseries must be a NamedTuple of vectors or a Vector of NamedTuples"))
    end
end

function simulate_reduced_deb_response(t::Vector{Float64}, C_event::Vector{Float64}, A_background::Real, axes_background, deb_params::DEBAxisParams;
                                      q::Float64=1.0,
                                      y0::Float64=0.0,
                                      Z_timeseries=nothing)
    ntime = length(t)

    if length(C_event) != ntime
        throw(ArgumentError("length(C_event) must equal length(t)"))
    end

    if Z_timeseries !== nothing && length(Z_timeseries) != ntime
        throw(ArgumentError("length(Z_timeseries) must equal length(t)"))
    end

    for i in 1:(ntime - 1)
        if t[i+1] <= t[i]
            throw(ArgumentError("t must be strictly increasing"))
        end
    end

    y = zeros(Float64, ntime)
    lambdas = zeros(Float64, ntime)

    if ntime > 0
        y[1] = y0
    end

    for i in 1:ntime
        Z_i = Z_timeseries !== nothing ? Z_timeseries[i] : nothing
        lambdas[i] = restoring_force_from_margin_and_axes(A_background, axes_background, deb_params; Z=Z_i)
    end

    for i in 1:(ntime - 1)
        dt = t[i+1] - t[i]
        lambda_i = lambdas[i]

        if lambda_i > 0
            y[i+1] = y[i] * exp(-lambda_i * dt) + (q * C_event[i] / lambda_i) * (1.0 - exp(-lambda_i * dt))
        else
            y[i+1] = y[i] + q * C_event[i] * dt
        end
    end

    return (
        t = t,
        y = y,
        lambda = lambdas,
        C_event = C_event
    )
end
