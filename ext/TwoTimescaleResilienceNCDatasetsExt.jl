module TwoTimescaleResilienceNCDatasetsExt

using TwoTimescaleResilience
using NCDatasets

function TwoTimescaleResilience.load_nc_layer(path::String, variable::String; time_index=nothing, bbox=nothing)
    NCDataset(path, "r") do ds
        lon_data = haskey(ds, "lon") ? ds["lon"][:] : (haskey(ds, "longitude") ? ds["longitude"][:] : collect(1:size(ds[variable], 1)))
        lat_data = haskey(ds, "lat") ? ds["lat"][:] : (haskey(ds, "latitude") ? ds["latitude"][:] : collect(1:size(ds[variable], 2)))
        
        var_data = ds[variable]
        
        if bbox !== nothing
            minlon, maxlon, minlat, maxlat = bbox
            lon_idx = findall(x -> minlon <= x <= maxlon, lon_data)
            lat_idx = findall(x -> minlat <= x <= maxlat, lat_data)
            
            lon_out = lon_data[lon_idx]
            lat_out = lat_data[lat_idx]
            
            if time_index !== nothing
                layer = var_data[lon_idx, lat_idx, time_index]
            else
                layer = var_data[lon_idx, lat_idx]
            end
        else
            lon_out = lon_data
            lat_out = lat_data
            if time_index !== nothing
                layer = var_data[:, :, time_index]
            else
                layer = var_data[:, :]
            end
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
        return out, lon_out, lat_out
    end
end

end
