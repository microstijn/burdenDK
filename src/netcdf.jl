function load_nc_layer end

function normalise_layer(layer::Matrix{Float64}; method=:minmax, lower=nothing, upper=nothing)
    nrows, ncols = size(layer)
    out = zeros(Float64, nrows, ncols)

    if method == :minmax
        valid_vals = filter(x -> !isnan(x), layer)
        min_val = isempty(valid_vals) ? 0.0 : minimum(valid_vals)
        max_val = isempty(valid_vals) ? 1.0 : maximum(valid_vals)
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
