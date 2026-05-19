module TwoTimescaleResilienceNCDatasetsExt

using TwoTimescaleResilience
using NCDatasets

function TwoTimescaleResilience.load_nc_layer(path::String, variable::String; time_index=nothing)
    NCDataset(path, "r") do ds
        var_data = ds[variable]
        if time_index !== nothing
            # Assuming format is (lon, lat, time) or similar
            layer = var_data[:, :, time_index]
        else
            layer = var_data[:, :]
        end

        # Convert to standard Matrix{Float64} handling missing values
        nrows, ncols = size(layer)
        out = zeros(Float64, nrows, ncols)
        for i in 1:nrows
            for j in 1:ncols
                val = layer[i, j]
                out[i, j] = ismissing(val) ? NaN : Float64(val)
            end
        end
        return out
    end
end

end
