using TwoTimescaleResilience

println("--- ISIMIP MoA DEB 3x3 Demo ---")

# 1. Create seven synthetic 3x3 layers
layers = [fill(0.0, 3, 3) for _ in 1:7]
# WT
layers[1] .= [20.0 22.0 24.0;
              21.0 25.0 26.0;
              19.0 20.0 21.0]
# BOD
layers[2] .= [2.0 5.0 8.0;
              3.0 6.0 9.0;
              1.0 2.0 3.0]
# TDS
layers[3] .= [100.0 150.0 200.0;
              110.0 160.0 210.0;
              90.0 100.0 110.0]
# FC
layers[4] .= [10.0 50.0 100.0;
              20.0 60.0 200.0;
              5.0 10.0 20.0]
# Nutrient
layers[5] .= [0.5 1.0 2.0;
              0.6 1.2 2.5;
              0.4 0.8 1.5]
# Chemical
layers[6] .= [0.01 0.05 0.1;
              0.02 0.06 0.15;
              0.005 0.01 0.02]
# Plastic
layers[7] .= [0.1 0.5 1.0;
              0.2 0.6 1.2;
              0.05 0.1 0.2]

# 2. Use fish profile
profile = fish_profile()

# 3. Run grid pipeline
pipeline_res = isimip_deb_pipeline_grid(layers, profile.exposure_filter, profile.moa_mapping, profile.moa_deb_mapping, profile.deb_params)

# 4. Expose intermediate variables
effective_layers = pipeline_res.effective_layers
modes = pipeline_res.modes
axes = pipeline_res.axes
Agrid = pipeline_res.A
lambdagrid = pipeline_res.lambda
Fgrid = pipeline_res.amplification

# 5. Print results
println("\nFish Profile Amplification Factor (Fgrid):")
display(Fgrid)

println("\nThermal Mode:")
display(modes.thermal)

println("\nOxygen Mode:")
display(modes.oxygen)

println("\nAssimilation Axis:")
display(axes.assimilation)

println("\nMaintenance Axis:")
display(axes.maintenance)

println("\nAdaptive Margin (Agrid):")
display(Agrid)

println("\nRestoring Force (lambda):")
display(lambdagrid)
