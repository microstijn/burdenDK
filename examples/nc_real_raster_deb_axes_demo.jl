# nc_real_raster_deb_axes_demo.jl
using CairoMakie
using TwoTimescaleResilience
# Note: To process actual NetCDF files, make sure NCDatasets is loaded in your environment
# and uncomment the following line:
# using NCDatasets 

"""
    run_nc_real_raster_deb_axes_demo(; 
        pathogen_nc_path="pathogen.nc", 
        organic_nc_path="organic.nc",
        output_dir="output",
        month_index=1,
        bbox=nothing
    )

Runs the DEB-like axes amplification pipeline on real NetCDF datasets.
The script assumes dimensions are (longitude, latitude, time).
"""
function run_nc_real_raster_deb_axes_demo(; 
    pathogen_nc_path="pathogen.nc", 
    organic_nc_path="organic.nc",
    output_dir="output_deb_demo",
    month_index=1,
    bbox=nothing,
    use_dummy_data=false
)
    mkpath(output_dir)
    
    if use_dummy_data
        println("Using dummy synthetic data for demonstration purposes...")
        lon = collect(range(-10.0, 10.0, length=50))
        lat = collect(range(-10.0, 10.0, length=40))
        pathogen_raw = rand(50, 40) .* 100.0
        organic_raw = rand(50, 40) .* 100.0
    else
        println("Loading NetCDF layers from month index ", month_index)
        # These functions come from ext/TwoTimescaleResilienceNCDatasetsExt.jl
        pathogen_raw, lon, lat = load_nc_layer(pathogen_nc_path, "pathogen", time_index=month_index, bbox=bbox)
        organic_raw, _, _ = load_nc_layer(organic_nc_path, "organic", time_index=month_index, bbox=bbox)
    end
    
    println("Normalising layers...")
    pathogen_norm = normalise_layer(pathogen_raw; method=:robust)
    organic_norm = normalise_layer(organic_raw; method=:robust)
    
    layers = [pathogen_norm, organic_norm]
    
    println("Configuring DEB axes mapping...")
    mapping = default_pathogen_organic_deb_mapping(interaction_strength=0.25)
    params = DEBAxisParams()
    
    println("Running DEB amplification pipeline...")
    result = deb_amplification_pipeline(layers, mapping, params)
    
    println("Exporting ASCII grids...")
    # Convert from (longitude, latitude) to (latitude, longitude) for ASCII format
    # The read_ascii_grid/write_ascii_grid works with matrices in shape (nrows, ncols) -> (lat, lon) typically,
    # but the exact orientation depends on the export convention. We use permutedims to export in (lat, lon) format.
    
    write_ascii_grid(joinpath(output_dir, "pathogen_normalised.asc"), permutedims(pathogen_norm))
    write_ascii_grid(joinpath(output_dir, "organic_normalised.asc"), permutedims(organic_norm))
    
    write_ascii_grid(joinpath(output_dir, "deb_assimilation_stress.asc"), permutedims(result.axes.assimilation))
    write_ascii_grid(joinpath(output_dir, "deb_maintenance_stress.asc"), permutedims(result.axes.maintenance))
    write_ascii_grid(joinpath(output_dir, "deb_growth_stress.asc"), permutedims(result.axes.growth))
    write_ascii_grid(joinpath(output_dir, "deb_reproduction_stress.asc"), permutedims(result.axes.reproduction))
    
    write_ascii_grid(joinpath(output_dir, "deb_adaptive_margin.asc"), permutedims(result.A))
    write_ascii_grid(joinpath(output_dir, "deb_restoring_force.asc"), permutedims(result.lambda))
    write_ascii_grid(joinpath(output_dir, "deb_amplification_factor.asc"), permutedims(result.amplification))
    
    println("Plotting PNG grids...")
    fig = Figure(size=(1200, 1000))
    
    # Pathogen
    ax1 = Axis(fig[1, 1], title="Pathogen (Norm)", aspect=DataAspect())
    hm1 = heatmap!(ax1, lon, lat, pathogen_norm, colormap=:Reds)
    Colorbar(fig[1, 2], hm1)
    
    # Organic
    ax2 = Axis(fig[1, 3], title="Organic (Norm)", aspect=DataAspect())
    hm2 = heatmap!(ax2, lon, lat, organic_norm, colormap=:Greens)
    Colorbar(fig[1, 4], hm2)
    
    # Assimilation Stress
    ax3 = Axis(fig[2, 1], title="Assimilation Stress", aspect=DataAspect())
    hm3 = heatmap!(ax3, lon, lat, result.axes.assimilation, colormap=:viridis)
    Colorbar(fig[2, 2], hm3)
    
    # Maintenance Stress
    ax4 = Axis(fig[2, 3], title="Maintenance Stress", aspect=DataAspect())
    hm4 = heatmap!(ax4, lon, lat, result.axes.maintenance, colormap=:viridis)
    Colorbar(fig[2, 4], hm4)
    
    # Adaptive Margin
    ax5 = Axis(fig[3, 1], title="Adaptive Margin (A_DEB)", aspect=DataAspect())
    hm5 = heatmap!(ax5, lon, lat, result.A, colormap=:plasma)
    Colorbar(fig[3, 2], hm5)
    
    # Amplification Factor
    ax6 = Axis(fig[3, 3], title="Amplification Factor (F)", aspect=DataAspect())
    hm6 = heatmap!(ax6, lon, lat, result.amplification, colormap=:inferno)
    Colorbar(fig[3, 4], hm6)
    
    save(joinpath(output_dir, "deb_pipeline_results.png"), fig)
    println("Done. Output saved to: ", output_dir)
end

# To run with dummy data (since real NetCDF files are not provided):
# run_nc_real_raster_deb_axes_demo(use_dummy_data=true)
