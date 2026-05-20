using TwoTimescaleResilience

println("--- Species Profile Comparison 3x3 Demo ---")

# 1. Create seven synthetic 3x3 layers
layers = [fill(0.0, 3, 3) for _ in 1:7]
layers[1] .= [20.0 22.0 24.0; 21.0 25.0 26.0; 19.0 20.0 21.0] # WT
layers[2] .= [2.0 5.0 8.0; 3.0 6.0 9.0; 1.0 2.0 3.0] # BOD
layers[3] .= [100.0 150.0 200.0; 110.0 160.0 210.0; 90.0 100.0 110.0] # TDS
layers[4] .= [10.0 50.0 100.0; 20.0 60.0 200.0; 5.0 10.0 20.0] # FC
layers[5] .= [0.5 1.0 2.0; 0.6 1.2 2.5; 0.4 0.8 1.5] # Nutrient
layers[6] .= [0.01 0.05 0.1; 0.02 0.06 0.15; 0.005 0.01 0.02] # Chemical
layers[7] .= [0.1 0.5 1.0; 0.2 0.6 1.2; 0.05 0.1 0.2] # Plastic

# 2. Profiles
profiles = [
    aquatic_invertebrate_profile(),
    fish_profile(),
    bivalve_profile(),
    human_profile()
]

# 3. Run and compare
for p in profiles
    println("\nRunning profile: $(p.name)")
    res = isimip_deb_pipeline_grid(layers, p.exposure_filter, p.moa_mapping, p.moa_deb_mapping, p.deb_params)

    println("Mean Amplification: ", sum(res.amplification)/9)
    println("Max Amplification:  ", maximum(res.amplification))
    println("Fgrid raster:")
    display(res.amplification)
end
