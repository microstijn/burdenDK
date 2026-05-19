"""
    load_nc_layer(path::String, variable::String; time_index=nothing, bbox=nothing)

Loads a NetCDF layer and returns `(layer, lon, lat)`.
`bbox` can be `(minlon, maxlon, minlat, maxlat)`.
"""
function load_nc_layer end

using Statistics

function normalise_layer(layer::Matrix{Float64}; method=:minmax, lower=nothing, upper=nothing, p=0.01)
    nrows, ncols = size(layer)
    out = zeros(Float64, nrows, ncols)

    if method == :minmax || method == :robust
        valid_vals = filter(x -> !isnan(x), layer)
        
        if isempty(valid_vals)
            min_val = 0.0
            max_val = 1.0
        elseif method == :robust
            min_val = quantile(valid_vals, p)
            max_val = quantile(valid_vals, 1.0 - p)
        else
            min_val = minimum(valid_vals)
            max_val = maximum(valid_vals)
        end
        
        range_val = max_val - min_val
        if range_val == 0.0
            range_val = 1.0
        end

        for i in 1:nrows
            for j in 1:ncols
                val = layer[i, j]
                if isnan(val)
                    out[i, j] = NaN
                else
                    norm_val = (val - min_val) / range_val

                    # Clamp to [0, 1] for robust method or if specifically requested
                    if method == :robust
                        norm_val = clamp(norm_val, 0.0, 1.0)
                    end

                    if lower !== nothing && norm_val < lower
                        norm_val = 0.0
                    end
                    if upper !== nothing && norm_val > upper
                        norm_val = 1.0
                    end

                    out[i, j] = norm_val
                end
            end
        end
    else
        error("Unknown normalisation method: $method")
    end

    return out
end
