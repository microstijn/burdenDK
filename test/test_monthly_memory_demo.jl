using Test
using TwoTimescaleResilience
using CSV
using DataFrames

# Include the example script to test its helper functions
include("../examples/ecotox_amp_multispecies_multicompound_monthly_memory_demo.jl")

@testset "Monthly Memory Demo Tranche 1" begin
    # Test 12-month scenario generator
    concs = generate_scenario_concentrations(1.0)
    @test length(concs) == 12

    # Test recovery to zero after pulse
    @test concs[6] == 0.0 # month 6 is after pulse
    @test concs[11] == 0.0 # month 11 is after moderate pulse
end

@testset "Monthly Memory Demo Tranche 2" begin
    amp = load_amp_species_library()

    # Check that species exist
    @test haskey(amp, "Abatus_cordatus")
    @test haskey(amp, "Podarcis_muralis")
    @test haskey(amp, "Thalia_democratica")

    params_ac = amp_species_deb_params(amp, "Abatus cordatus")
    @test typeof(params_ac) == TwoTimescaleResilience.DEBAxisParams

    ecotox = load_ecotox_library()
    memory = load_compound_memory_library()

    nacl_cas = "7647-14-5"
    cd_cas = "7440-43-9"
    hg_cas = "7439-97-6"

    rho_nacl = compound_retention(nacl_cas; memory_library=memory)
    rho_cd = compound_retention(cd_cas; memory_library=memory)
    rho_hg = compound_retention(hg_cas; memory_library=memory)

    K_nacl = compound_bioaccumulation_factor(nacl_cas; memory_library=memory)
    K_cd = compound_bioaccumulation_factor(cd_cas; memory_library=memory)
    K_hg = compound_bioaccumulation_factor(hg_cas; memory_library=memory)

    @test rho_cd > rho_nacl
    @test rho_hg > rho_nacl
    @test K_cd >= K_nacl
    @test K_hg >= K_nacl

    rec_nacl = get_ecotox_record(ecotox, nacl_cas)
    rec_cd = get_ecotox_record(ecotox, cd_cas)
    rec_hg = get_ecotox_record(ecotox, hg_cas)

    @test !isempty(rec_nacl)
    @test !isempty(rec_cd)
    @test !isempty(rec_hg)
end

@testset "Monthly Memory Demo Tranche 3" begin
    memory = load_compound_memory_library()

    # Controlled toy test with known rho and K
    state = EcotoxExposureState()
    cd_cas = "7440-43-9"
    rho_cd = compound_retention(cd_cas; memory_library=memory)
    K_cd = compound_bioaccumulation_factor(cd_cas; memory_library=memory)

    ecotox = load_ecotox_library()
    rec_cd = get_ecotox_record(ecotox, cd_cas)

    C1 = 10.0
    C2 = 0.0

    r1 = Dict{String, Any}()
    for (k,v) in rec_cd
        r1[k] = v
    end
    r1["NOEC_median"] = Float64(rec_cd["NOEC_median"])
    r1["EC50_median"] = Float64(rec_cd["EC50_median"])

    stateful_concs_1 = Dict(rec_cd["cas_norm"] => C1)
    ecotox_records_to_deb_burden_stateful!(state, stateful_concs_1, [rec_cd]; memory_library=memory)
    B1 = get_internal_burden(state, cd_cas)

    expected_B1 = rho_cd * 0.0 + (1 - rho_cd) * K_cd * C1
    @test isapprox(B1, expected_B1; atol=1e-5)

    stateful_concs_2 = Dict(rec_cd["cas_norm"] => C2)
    ecotox_records_to_deb_burden_stateful!(state, stateful_concs_2, [rec_cd]; memory_library=memory)
    B2 = get_internal_burden(state, cd_cas)

    expected_B2 = rho_cd * expected_B1 + (1 - rho_cd) * K_cd * C2
    @test isapprox(B2, expected_B2; atol=1e-5)

    # Test rho=0, K=1 (NaCl)
    state_nacl = EcotoxExposureState()
    nacl_cas = "7647-14-5"
    rec_nacl = get_ecotox_record(ecotox, nacl_cas)
    r2 = Dict{String, Any}()
    for (k,v) in rec_nacl
        r2[k] = v
    end
    r2["NOEC_median"] = Float64(rec_nacl["NOEC_median"])
    r2["EC50_median"] = Float64(rec_nacl["EC50_median"])

    C_nacl = 5.0
    stateful_concs_nacl = Dict(rec_nacl["cas_norm"] => C_nacl)
    ecotox_records_to_deb_burden_stateful!(state_nacl, stateful_concs_nacl, [rec_nacl]; memory_library=memory)
    B_nacl = get_internal_burden(state_nacl, nacl_cas)
    @test isapprox(B_nacl, C_nacl; atol=1e-5)

    # Test persistence of Cd
    @test B2 > 0.0 # B2 is after C2=0.0
end

@testset "Monthly Memory Demo Tranche 4" begin
    # Test aggregation logic without mixture toxicity
    amp = load_amp_species_library()
    params = amp_species_deb_params(amp, "Abatus cordatus")

    # 1. Zero burden -> amplification 1.0, baseline A
    zero_burden = (assimilation=0.0, maintenance=0.0, growth=0.0, reproduction=0.0)
    zero_resp = ecotox_burden_to_response(zero_burden, params)

    @test zero_resp.A == params.A0
    @test isapprox(zero_resp.amplification, 1.0; atol=1e-5)

    # 2. Some burden -> reduced A, lambda > lambda0, F > 1.0
    some_burden = (assimilation=0.0, maintenance=10.0, growth=0.0, reproduction=0.0)
    some_resp = ecotox_burden_to_response(some_burden, params)

    @test some_resp.A < zero_resp.A
    @test some_resp.lambda < zero_resp.lambda
    @test some_resp.amplification > 1.0

    # Check bounded lambda
    @test some_resp.lambda >= params.lambda_min
    @test some_resp.lambda <= params.lambda_max

    # Aggregation test wrapper
    function aggregate_burdens(burdens)
        return (
            assimilation = sum(b.assimilation for b in burdens),
            maintenance = sum(b.maintenance for b in burdens),
            growth = sum(b.growth for b in burdens),
            reproduction = sum(b.reproduction for b in burdens)
        )
    end

    # Multi-compound summation
    b1 = (assimilation=0.0, maintenance=5.0, growth=0.0, reproduction=0.0)
    b2 = (assimilation=0.0, maintenance=2.0, growth=0.0, reproduction=0.0)
    agg = aggregate_burdens([b1, b2])

    @test agg.maintenance == 7.0
    agg_resp = ecotox_burden_to_response(agg, params)
    @test agg_resp.amplification > 1.0
end

@testset "Monthly Memory Demo Tranche 5" begin
    # Test that the script creates the expected output CSVs
    output_dir = normpath(joinpath(@__DIR__, "..", "output", "ecotox_amp_multispecies_multicompound_monthly_memory_demo"))

    # Run the main function to generate the CSVs
    main()

    @test isdir(output_dir)

    compound_csv = joinpath(output_dir, "monthly_compound_summary.csv")
    species_csv = joinpath(output_dir, "monthly_species_summary.csv")

    @test isfile(compound_csv)
    @test isfile(species_csv)

    df_comp = CSV.read(compound_csv, DataFrame)
    df_spec = CSV.read(species_csv, DataFrame)

    @test nrow(df_comp) == 108 # 3 species * 3 compounds * 12 months
    @test nrow(df_spec) == 36  # 3 species * 12 months

    # Check compound columns
    expected_comp_cols = ["species_name", "month", "cas_norm", "chemical_name", "C_t", "B_t", "x_t", "burden_assimilation", "burden_maintenance", "burden_growth", "burden_reproduction"]
    for col in expected_comp_cols
        @test col in names(df_comp)
    end

    # Check species columns
    expected_spec_cols = ["species_key", "species_name", "month", "A_t", "lambda_t", "lambda0", "F_t"]
    for col in expected_spec_cols
        @test col in names(df_spec)
    end

    # Test months 1-12
    @test sort(unique(df_comp.month)) == collect(1:12)
    @test sort(unique(df_spec.month)) == collect(1:12)
end

@testset "Monthly Memory Demo Tranche 6" begin
    # Assuming main() was already run in Tranche 5, check for PNG plots
    output_dir = normpath(joinpath(@__DIR__, "..", "output", "ecotox_amp_multispecies_multicompound_monthly_memory_demo"))

    plots = [
        "monthly_concentrations.png",
        "monthly_internal_burdens.png",
        "monthly_active_stress.png",
        "monthly_axis_burdens.png",
        "monthly_adaptive_margin.png",
        "monthly_restoring_force.png",
        "monthly_amplification.png"
    ]

    for p in plots
        path = joinpath(output_dir, p)
        @test isfile(path)
        if isfile(path)
            # Check file is non-empty
            @test filesize(path) > 100
        end
    end
end

@testset "Monthly Memory Demo Tranche 7" begin
    # The end-to-end examples test
    # Re-run it once safely to make sure it functions as a standalone command
    try
        main()
        @test true
    catch e
        @test false
    end
end
