# ==============================================================================
# Beamer raster export add-on for nc_monthly_longterm_isimip_moa_deb_inspection.jl
#
# Paste this block into the script after monthly_summary(...) helper and before
# the main workflow. Then call export_example_month_rasters(...) after the example
# month pipeline, and export_selected_time_rasters(...) after the long-term loop.
# ============================================================================== 

# Toggle raster export for slides / Beamer figures.
export_beamer_rasters = true
beamer_output_dir = joinpath(output_dir, "beamer_rasters")
mkpath(beamer_output_dir)

# ASCII export is useful for reproducibility; PNG export is useful for Beamer.
export_ascii_rasters = true
export_png_rasters = true

# Optional: keep consistent colour ranges across selected years for key outputs.
# For most normalised input/mode/axis rasters, 0..1 is appropriate.
input_colorrange = (0.0, 1.0)
mode_colorrange = (0.0, 1.0)
axis_colorrange = (0.0, 1.0)
A_colorrange = nothing
lambda_colorrange = nothing
F_colorrange = nothing

function safe_name(s::AbstractString)
    return replace(s, " " => "_", "/" => "_", "-" => "_", ":" => "_", "." => "p")
end

function export_ascii_grid_if_requested(grid, filename)
    if export_ascii_rasters
        try
            write_ascii_grid(filename, grid)
        catch err
            @warn "Could not write ASCII grid $filename" exception=(err, catch_backtrace())
        end
    end
end

function export_png_grid_if_requested(grid, lon, lat; title, filename, colormap=:viridis, colorrange=nothing)
    if export_png_rasters
        try
            plot_geo_grid(grid, lon, lat; title=title, filename=filename, colormap=colormap, colorrange=colorrange)
        catch err
            @warn "Could not write PNG grid $filename" exception=(err, catch_backtrace())
        end
    end
end

function export_beamer_grid(grid, lon, lat; key, group, name, title, colormap=:viridis, colorrange=nothing)
    if !export_beamer_rasters
        return nothing
    end

    group_dir = joinpath(beamer_output_dir, safe_name(group))
    mkpath(group_dir)

    base = "$(safe_name(key))__$(safe_name(group))__$(safe_name(name))"
    png_file = joinpath(group_dir, base * ".png")
    asc_file = joinpath(group_dir, base * ".asc")

    export_png_grid_if_requested(grid, lon, lat; title=title, filename=png_file, colormap=colormap, colorrange=colorrange)
    export_ascii_grid_if_requested(grid, asc_file)

    return (png=png_file, asc=asc_file)
end

function export_example_month_rasters(; key, lon, lat,
    raw_layers, norm_layers, effective_layers, modes, axes, Agrid, lambdagrid, Fgrid)

    println("Exporting beamer rasters for example month: $key")

    # Raw variables. Colour ranges are automatic because raw real NetCDF values may not be 0..1.
    for (name, grid) in pairs(raw_layers)
        export_beamer_grid(grid, lon, lat;
            key=key, group="01_raw_inputs", name=String(name),
            title="Raw $(String(name)) ($key)", colormap=:viridis, colorrange=nothing)
    end

    # Normalised variables.
    for (name, grid) in pairs(norm_layers)
        export_beamer_grid(grid, lon, lat;
            key=key, group="02_normalised_inputs", name=String(name),
            title="Normalised $(String(name)) ($key)", colormap=:viridis, colorrange=input_colorrange)
    end

    # Effective exposure layers are a vector in canonical order.
    effective_names = (:WT, :BOD, :TDS, :FC, :Nutrient, :Chemical, :Plastic)
    for (idx, grid) in enumerate(effective_layers)
        name = String(effective_names[idx])
        export_beamer_grid(grid, lon, lat;
            key=key, group="03_effective_exposure", name=name,
            title="Effective exposure $(name) ($key)", colormap=:viridis, colorrange=input_colorrange)
    end

    # Mode-of-action rasters.
    for name in (:thermal, :oxygen, :osmotic, :immune, :eutrophication, :toxic, :feeding, :physical)
        grid = getproperty(modes, name)
        export_beamer_grid(grid, lon, lat;
            key=key, group="04_modes", name=String(name),
            title="MoA $(String(name)) ($key)", colormap=:magma, colorrange=mode_colorrange)
    end

    # DEB axis rasters.
    for name in (:assimilation, :maintenance, :growth, :reproduction)
        grid = getproperty(axes, name)
        export_beamer_grid(grid, lon, lat;
            key=key, group="05_deb_axes", name=String(name),
            title="DEB axis $(String(name)) ($key)", colormap=:plasma, colorrange=axis_colorrange)
    end

    # Final response-operator rasters.
    export_beamer_grid(Agrid, lon, lat;
        key=key, group="06_response_operator", name="adaptive_margin_A",
        title="Adaptive margin A ($key)", colormap=:viridis, colorrange=A_colorrange)

    export_beamer_grid(lambdagrid, lon, lat;
        key=key, group="06_response_operator", name="restoring_force_lambda",
        title="Restoring force lambda ($key)", colormap=:viridis, colorrange=lambda_colorrange)

    export_beamer_grid(Fgrid, lon, lat;
        key=key, group="06_response_operator", name="amplification_F",
        title="Amplification factor F ($key)", colormap=:inferno, colorrange=F_colorrange)

    return nothing
end

function export_selected_time_rasters(; lon, lat,
    raw_layers_by_time, norm_layers_by_time, effective_layers_by_time,
    modes_by_time, axes_by_time, A_by_time, lambda_by_time, amplification_by_time)

    if !export_beamer_rasters
        return nothing
    end

    println("Exporting beamer rasters for all selected times...")

    for key in sort(collect(keys(amplification_by_time)))
        export_example_month_rasters(
            key=key,
            lon=lon,
            lat=lat,
            raw_layers=raw_layers_by_time[key],
            norm_layers=norm_layers_by_time[key],
            effective_layers=effective_layers_by_time[key],
            modes=modes_by_time[key],
            axes=axes_by_time[key],
            Agrid=A_by_time[key],
            lambdagrid=lambda_by_time[key],
            Fgrid=amplification_by_time[key]
        )
    end

    # First-last difference rasters for Beamer.
    sorted_keys = sort(collect(keys(amplification_by_time)))
    if length(sorted_keys) >= 2
        first_key = first(sorted_keys)
        last_key = last(sorted_keys)

        Fdiff = amplification_by_time[last_key] .- amplification_by_time[first_key]
        fv = finite_values(Fdiff)
        max_abs = isempty(fv) ? 1.0 : maximum(abs.(fv))

        export_beamer_grid(Fdiff, lon, lat;
            key="$(first_key)_to_$(last_key)",
            group="07_differences",
            name="amplification_F_difference",
            title="Amplification difference F: $(last_key) minus $(first_key)",
            colormap=:RdBu_11,
            colorrange=(-max_abs, max_abs))

        Adiff = A_by_time[last_key] .- A_by_time[first_key]
        fvA = finite_values(Adiff)
        max_abs_A = isempty(fvA) ? 1.0 : maximum(abs.(fvA))

        export_beamer_grid(Adiff, lon, lat;
            key="$(first_key)_to_$(last_key)",
            group="07_differences",
            name="adaptive_margin_A_difference",
            title="Adaptive margin difference A: $(last_key) minus $(first_key)",
            colormap=:RdBu_11,
            colorrange=(-max_abs_A, max_abs_A))
    end

    return nothing
end

# ==============================================================================
# Call locations
# ==============================================================================
# After the example-month pipeline and after WT_norm...Fgrid are defined, call:
#
example_key = string(example_year, "-", lpad(example_month, 2, "0"))
export_example_month_rasters(
    key=example_key,
    lon=lon,
    lat=lat,
    raw_layers=(WT=WT_raw, BOD=BOD_raw, TDS=TDS_raw, FC=FC_raw,
                Nutrient=Nutrient_raw, Chemical=Chemical_raw, Plastic=Plastic_raw),
    norm_layers=(WT=WT_norm, BOD=BOD_norm, TDS=TDS_norm, FC=FC_norm,
                 Nutrient=Nutrient_norm, Chemical=Chemical_norm, Plastic=Plastic_norm),
    effective_layers=effective_layers,
    modes=modes,
    axes=axes,
    Agrid=Agrid,
    lambdagrid=lambdagrid,
    Fgrid=Fgrid
)

export_selected_time_rasters(
    lon=lon,
    lat=lat,
    raw_layers_by_time=raw_layers_by_time,
    norm_layers_by_time=norm_layers_by_time,
    effective_layers_by_time=effective_layers_by_time,
    modes_by_time=modes_by_time,
    axes_by_time=axes_by_time,
    A_by_time=A_by_time,
    lambda_by_time=lambda_by_time,
    amplification_by_time=amplification_by_time
)
