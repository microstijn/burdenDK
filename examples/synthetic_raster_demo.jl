using TwoTimescaleResilience

println("Running synthetic raster demo...")
params = BackgroundParams()
run_synthetic_raster_demo(params; output_dir=".")
println("Demo complete. Check the current directory for .asc and .png files.")
