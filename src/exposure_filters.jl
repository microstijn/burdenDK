export ExposureFilter
export apply_exposure_filter, apply_exposure_filter_grid
export aquatic_exposure_filter, human_exposure_filter, default_exposure_filter

Base.@kwdef struct ExposureFilter
    multipliers::Vector{Float64}
    name::String = "generic"
    description::String = ""
end

function apply_exposure_filter(values::AbstractVector{<:Real}, filter::ExposureFilter)
    if length(values) != length(filter.multipliers)
        throw(ArgumentError("Length of values ($(length(values))) must equal length of multipliers ($(length(filter.multipliers)))"))
    end
    return filter.multipliers .* values
end

function apply_exposure_filter_grid(layers::Vector{Matrix{Float64}}, filter::ExposureFilter)
    n_vars = length(layers)
    if n_vars == 0
        throw(ArgumentError("At least one layer must be provided"))
    end
    if n_vars != length(filter.multipliers)
        throw(ArgumentError("Number of layers ($n_vars) must equal length of multipliers ($(length(filter.multipliers)))"))
    end

    dims = size(layers[1])
    for i in 2:n_vars
        if size(layers[i]) != dims
            throw(ArgumentError("All layers must have the same size"))
        end
    end

    out_layers = [Matrix{Float64}(undef, dims) for _ in 1:n_vars]

    for c in 1:dims[2]
        for r in 1:dims[1]
            for i in 1:n_vars
                v = layers[i][r, c]
                if isnan(v)
                    out_layers[i][r, c] = NaN
                else
                    out_layers[i][r, c] = v * filter.multipliers[i]
                end
            end
        end
    end

    return out_layers
end

function aquatic_exposure_filter(nvars::Int)
    return ExposureFilter(
        multipliers = fill(1.0, nvars),
        name = "aquatic",
        description = "Direct aquatic exposure, multipliers are 1.0"
    )
end

function default_exposure_filter(nvars::Int)
    return aquatic_exposure_filter(nvars)
end

function human_exposure_filter()
    return ExposureFilter(
        multipliers = [0.10, 0.25, 0.50, 0.60, 0.35, 0.40, 0.30],
        name = "human",
        description = "This is a placeholder contact/use filter, not calibrated human health exposure."
    )
end
