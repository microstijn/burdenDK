using Test
using TwoTimescaleResilience

@testset "ECOTOX Runtime Adapter: Tranche 1" begin
    lib_path = joinpath(@__DIR__, "..", "data", "ECOTOX_Toxicity_Library.json")

    @test isfile(lib_path)

    lib = load_ecotox_library(lib_path)
    @test length(lib) > 0

    r = first(lib)
    @test validate_ecotox_record(r) == true

    for key in ["cas", "cas_norm", "cas_hyphenated", "taxon_class", "effect_code", "NOEC_median", "EC50_median", "n_NOEC", "n_EC50"]
        @test haskey(r, key)
    end

    bad = deepcopy(r)
    bad["EC50_median"] = bad["NOEC_median"]
    @test_throws ArgumentError validate_ecotox_record(bad)

    bad2 = deepcopy(r)
    delete!(bad2, "effect_code")
    @test_throws ArgumentError validate_ecotox_record(bad2)

    bad3 = deepcopy(r)
    delete!(bad3, "NOEC_median")
    @test_throws ArgumentError validate_ecotox_record(bad3)

    bad4 = deepcopy(r)
    bad4["NOEC_median"] = -1.0
    @test_throws ArgumentError validate_ecotox_record(bad4)

    bad5 = deepcopy(r)
    bad5["cas_norm"] = ""
    @test_throws ArgumentError validate_ecotox_record(bad5)

    @test ecotox_active_stress(1.0, 2.0, 4.5) == 0.0
    @test ecotox_active_stress(2.0, 2.0, 4.5) == 0.0
    @test ecotox_active_stress(4.5, 2.0, 4.5) ≈ 1.0
    @test ecotox_active_stress(12.0, 2.0, 4.5) ≈ 4.0

    @test_throws ArgumentError ecotox_active_stress(3.0, 2.0, 2.0)
    @test_throws ArgumentError ecotox_active_stress(3.0, -1.0, 2.0)
    @test_throws ArgumentError ecotox_active_stress(NaN, 2.0, 4.5)
    @test_throws ArgumentError ecotox_active_stress(3.0, NaN, 4.5)
    @test_throws ArgumentError ecotox_active_stress(3.0, 2.0, NaN)

    NOEC = r["NOEC_median"]
    EC50 = r["EC50_median"]

    @test ecotox_active_stress(NOEC, NOEC, EC50) == 0.0
    @test ecotox_active_stress(EC50, NOEC, EC50) ≈ 1.0

    @test ecotox_effect_to_deb_axis("MOR") == :maintenance
    @test ecotox_effect_to_deb_axis("mor") == :maintenance
    @test ecotox_effect_to_deb_axis(" GRO ") == :growth
    @test ecotox_effect_to_deb_axis("REP") == :reproduction
    @test ecotox_effect_to_deb_axis("FEC") == :reproduction
    @test ecotox_effect_to_deb_axis("DVP") == :growth
    @test ecotox_effect_to_deb_axis("BCM") == :assimilation
    @test ecotox_effect_to_deb_axis("BEH") == :maintenance
    @test ecotox_effect_to_deb_axis("POP") == :maintenance
    @test ecotox_effect_to_deb_axis("ENZ") == :maintenance
    @test ecotox_effect_to_deb_axis("FDB") == :maintenance

    @test deb_axis_index(:assimilation) == 1
    @test deb_axis_index(:maintenance) == 2
    @test deb_axis_index(:growth) == 3
    @test deb_axis_index(:reproduction) == 4

    @test_throws ArgumentError ecotox_effect_to_deb_axis("UNKNOWN")
    @test_throws ArgumentError ecotox_effect_to_deb_axis("")
    @test_throws ArgumentError deb_axis_index(:unknown)

    effects = unique(string(rec["effect_code"]) for rec in lib)
    for eff in effects
        try
            ecotox_effect_to_deb_axis(eff)
        catch
            # Expected that some effect codes may be unsupported
        end
    end

    record = Dict(
        "cas" => "50000",
        "cas_norm" => "50000",
        "cas_hyphenated" => "50-00-0",
        "taxon_class" => "Actinopterygii",
        "effect_code" => "MOR",
        "NOEC_median" => 2.0,
        "EC50_median" => 4.5,
        "n_NOEC" => 1,
        "n_EC50" => 1
    )

    burden = ecotox_record_to_deb_burden(12.0, record)
    @test burden.assimilation == 0.0
    @test burden.maintenance ≈ 4.0
    @test burden.growth == 0.0
    @test burden.reproduction == 0.0

    burden0 = ecotox_record_to_deb_burden(1.0, record)
    @test burden0.maintenance == 0.0

    record_gro = deepcopy(record)
    record_gro["effect_code"] = "GRO"
    burden_gro = ecotox_record_to_deb_burden(4.5, record_gro)
    @test burden_gro.growth ≈ 1.0
    @test burden_gro.maintenance == 0.0

    record_rep = deepcopy(record)
    record_rep["effect_code"] = "REP"
    burden_rep = ecotox_record_to_deb_burden(4.5, record_rep)
    @test burden_rep.reproduction ≈ 1.0

    record_bcm = deepcopy(record)
    record_bcm["effect_code"] = "BCM"
    burden_bcm = ecotox_record_to_deb_burden(4.5, record_bcm)
    @test burden_bcm.assimilation ≈ 1.0

    bad_record = deepcopy(record)
    bad_record["EC50_median"] = 2.0
    @test_throws ArgumentError ecotox_record_to_deb_burden(4.5, bad_record)

    record_mor = Dict("cas_norm" => "50000", "cas" => "50000", "cas_hyphenated" => "50-00-0", "taxon_class" => "Actinopterygii", "effect_code" => "MOR", "NOEC_median" => 2.0, "EC50_median" => 4.5, "n_NOEC" => 1, "n_EC50" => 1)
    record_gro_m = Dict("cas_norm" => "50000", "cas" => "50000", "cas_hyphenated" => "50-00-0", "taxon_class" => "Actinopterygii", "effect_code" => "GRO", "NOEC_median" => 2.0, "EC50_median" => 4.5, "n_NOEC" => 1, "n_EC50" => 1)
    record_rep_m = Dict("cas_norm" => "7440439", "cas" => "7440439", "cas_hyphenated" => "7440-43-9", "taxon_class" => "Actinopterygii", "effect_code" => "REP", "NOEC_median" => 1.0, "EC50_median" => 3.0, "n_NOEC" => 1, "n_EC50" => 1)

    concentrations = Dict(
        "50-00-0" => 4.5,
        "7440-43-9" => 3.0
    )

    burden_multi = ecotox_records_to_deb_burden(concentrations, [record_mor, record_gro_m, record_rep_m])
    @test burden_multi.maintenance ≈ 1.0
    @test burden_multi.growth ≈ 1.0
    @test burden_multi.reproduction ≈ 1.0
    @test burden_multi.assimilation == 0.0

    burden_missing = ecotox_records_to_deb_burden(Dict("50-00-0" => 4.5), [record_rep_m])
    @test burden_missing.assimilation == 0.0
    @test burden_missing.maintenance == 0.0
    @test burden_missing.growth == 0.0
    @test burden_missing.reproduction == 0.0

    burden_empty = ecotox_records_to_deb_burden(concentrations, [])
    @test burden_empty.assimilation == 0.0
    @test burden_empty.maintenance == 0.0
    @test burden_empty.growth == 0.0
    @test burden_empty.reproduction == 0.0

    concentrations2 = Dict("50000" => 4.5)
    burden2 = ecotox_records_to_deb_burden(concentrations2, [record_mor])
    @test burden2.maintenance ≈ 1.0

    burden_zero = (
        assimilation = 0.0,
        maintenance = 0.0,
        growth = 0.0,
        reproduction = 0.0
    )

    params = DEBAxisParams(
        A0 = 100.0,
        alpha_axes = (1.0, 2.0, 3.0, 4.0),
        lambda_min = 0.1,
        lambda_max = 1.0,
        KA = 50.0
    )

    response_zero = ecotox_burden_to_response(burden_zero, params)
    @test response_zero.A ≈ params.A0
    @test response_zero.amplification ≈ 1.0

    burden_nonzero = (
        assimilation = 1.0,
        maintenance = 2.0,
        growth = 0.0,
        reproduction = 0.0
    )

    response = ecotox_burden_to_response(burden_nonzero, params)
    @test response.A ≈ 95.0
    @test response.lambda <= restoring_force_from_margin(params.A0, params)
    @test response.amplification >= 1.0

    rec_synthetic = Dict(
        "cas" => "123",
        "cas_norm" => "123",
        "cas_hyphenated" => "12-3",
        "taxon_class" => "Actinopterygii",
        "effect_code" => "MOR",
        "NOEC_median" => 2.0,
        "EC50_median" => 4.5,
        "n_NOEC" => 1,
        "n_EC50" => 1
    )
    b_syn = ecotox_record_to_deb_burden(12.0, rec_synthetic)
    resp_syn = ecotox_burden_to_response(b_syn, params)

    @test b_syn.maintenance ≈ 4.0
    @test resp_syn.A < params.A0
    @test resp_syn.amplification >= 1.0

    # Tranche 7
    lib = load_ecotox_library(lib_path)
    if length(lib) > 0
        r_first = first(lib)
        cas_h = r_first["cas_hyphenated"]
        cas_n = r_first["cas_norm"]

        a = ecotox_filter_records(lib; cas=cas_h)
        b = ecotox_filter_records(lib; cas=cas_n)

        @test length(a) == length(b)
        @test length(a) > 0

        tax = r_first["taxon_class"]
        t = ecotox_filter_records(lib; taxon_class=tax)
        @test length(t) > 0
        @test all(lowercase(strip(string(x["taxon_class"]))) == lowercase(strip(string(tax))) for x in t)

        eff = r_first["effect_code"]
        e = ecotox_filter_records(lib; effect_code=lowercase(eff))
        @test length(e) > 0
        @test all(uppercase(strip(string(x["effect_code"]))) == uppercase(strip(string(eff))) for x in e)

        c = ecotox_filter_records(lib; cas=cas_h, taxon_class=tax, effect_code=eff)
        @test length(c) > 0

        c_t = ecotox_records_for_taxon(lib, tax; cas=cas_h, effect_code=eff)
        @test length(c) == length(c_t)

        none = ecotox_filter_records(lib; cas="000-00-0")
        @test isempty(none)

        # test vector inputs
        vec_test = ecotox_filter_records(lib; cas=[cas_h, "000-00-0"])
        @test length(vec_test) == length(a)

        vec_test_tax = ecotox_filter_records(lib; taxon_class=[tax, "UnknownTaxon"])
        @test length(vec_test_tax) == length(t)
    end
end

@testset "ECOTOX Memory Library: Tranche 1" begin
    using CSV
    memory_path = joinpath(@__DIR__, "..", "data", "Compound_Memory_Library.csv")
    @test isfile(memory_path)

    df = CSV.File(memory_path)

    cols = propertynames(df)
    @test :cas_norm in cols
    @test :cas_hyphenated in cols
    @test :chemical_name in cols
    @test :memory_class in cols
    @test :retention_rho_monthly in cols
    @test :basis in cols
    @test :confidence in cols
    @test :notes in cols

    cas_norms = string.(df.cas_norm)
    @test "7440439" in cas_norms
    @test "7647145" in cas_norms
    @test "7439976" in cas_norms
    @test "9002884" in cas_norms

    @test all(0.0 .<= df.retention_rho_monthly .< 1.0)
end

@testset "ECOTOX Memory Library: Tranche 2" begin
    memory = load_compound_memory_library()
    @test length(memory) > 0

    @test validate_compound_memory_record(first(memory)) == true
    @test all(validate_compound_memory_record(r) for r in memory)

    bad1 = deepcopy(first(memory))
    bad1 = typeof(bad1)(k => k == "cas_norm" ? "" : v for (k,v) in pairs(bad1))
    @test_throws ArgumentError validate_compound_memory_record(bad1)

    bad2 = deepcopy(first(memory))
    bad2 = typeof(bad2)(k => k == "chemical_name" ? "" : v for (k,v) in pairs(bad2))
    @test_throws ArgumentError validate_compound_memory_record(bad2)

    bad3 = deepcopy(first(memory))
    bad3 = typeof(bad3)(k => k == "retention_rho_monthly" ? -0.1 : v for (k,v) in pairs(bad3))
    @test_throws ArgumentError validate_compound_memory_record(bad3)

    bad4 = deepcopy(first(memory))
    bad4 = typeof(bad4)(k => k == "retention_rho_monthly" ? 1.0 : v for (k,v) in pairs(bad4))
    @test_throws ArgumentError validate_compound_memory_record(bad4)

    bad5 = deepcopy(first(memory))
    bad5 = typeof(bad5)(k => k == "retention_rho_monthly" ? NaN : v for (k,v) in pairs(bad5))
    @test_throws ArgumentError validate_compound_memory_record(bad5)

    bad6 = deepcopy(first(memory))
    delete!(bad6, "cas_norm")
    @test_throws ArgumentError validate_compound_memory_record(bad6)
end

@testset "ECOTOX Memory Library: Tranche 3" begin
    memory = load_compound_memory_library()

    @test compound_retention("7647-14-5"; memory_library=memory) == 0.0
    @test compound_retention("7647145"; memory_library=memory) == 0.0

    @test compound_retention("7440-43-9"; memory_library=memory) ≈ 0.90
    @test compound_retention("7440439"; memory_library=memory) ≈ 0.90

    @test compound_retention("7439-97-6"; memory_library=memory) ≈ 0.95
    @test compound_retention("9002-88-4"; memory_library=memory) ≈ 0.95

    @test compound_retention("000-00-0"; memory_library=memory) == 0.0
    @test compound_retention("---"; memory_library=memory) == 0.0

    @test ecotox_default_retention("7440-43-9"; memory_library=memory) ≈ compound_retention("7440-43-9"; memory_library=memory)
end

@testset "ECOTOX Memory Library: Tranche 4" begin
    state = EcotoxExposureState()

    @test get_internal_burden(state, "7440-43-9") == 0.0

    set_internal_burden!(state, "7440-43-9", 12.5)

    @test get_internal_burden(state, "7440-43-9") ≈ 12.5
    @test get_internal_burden(state, "7440439") ≈ 12.5

    @test_throws ArgumentError set_internal_burden!(state, "7440-43-9", -1.0)
    @test_throws ArgumentError set_internal_burden!(state, "7440-43-9", NaN)

    reset_internal_burdens!(state)
    @test get_internal_burden(state, "7440-43-9") == 0.0
end

@testset "ECOTOX Memory Library: Tranche 5" begin
    state = EcotoxExposureState()
    B = update_internal_burden!(state, "7647-14-5", 10.0; retention=0.0)
    @test B ≈ 10.0

    B2 = update_internal_burden!(state, "7647-14-5", 0.0; retention=0.0)
    @test B2 ≈ 0.0

    state = EcotoxExposureState()

    B1 = update_internal_burden!(state, "7440-43-9", 10.0; retention=0.9)
    @test B1 ≈ 1.0

    B2 = update_internal_burden!(state, "7440-43-9", 10.0; retention=0.9)
    @test B2 ≈ 1.9

    B3 = update_internal_burden!(state, "7440-43-9", 0.0; retention=0.9)
    @test B3 ≈ 1.71

    memory = load_compound_memory_library()
    state = EcotoxExposureState()
    B = update_internal_burden!(state, "7440-43-9", 10.0; memory_library=memory)
    @test B ≈ 1.0

    @test_throws ArgumentError update_internal_burden!(state, "7440-43-9", -1.0)
    @test_throws ArgumentError update_internal_burden!(state, "7440-43-9", NaN)
    @test_throws ArgumentError update_internal_burden!(state, "7440-43-9", 10.0; retention=-0.1)
    @test_throws ArgumentError update_internal_burden!(state, "7440-43-9", 10.0; retention=1.0)
end

@testset "ECOTOX Memory Library: Tranche 6" begin
    record = Dict(
        "cas" => "7440439",
        "cas_norm" => "7440439",
        "cas_hyphenated" => "7440-43-9",
        "taxon_class" => "Actinopterygii",
        "effect_code" => "MOR",
        "NOEC_median" => 2.0,
        "EC50_median" => 4.0,
        "n_NOEC" => 1,
        "n_EC50" => 1
    )

    state = EcotoxExposureState()

    # Month 1
    concentrations = Dict("7440-43-9" => 10.0)
    burden1 = ecotox_records_to_deb_burden_stateful!(state, concentrations, [record]; retention=0.9)

    @test get_internal_burden(state, "7440-43-9") ≈ 1.0
    @test burden1.maintenance == 0.0

    # Month 2
    burden2 = ecotox_records_to_deb_burden_stateful!(state, concentrations, [record]; retention=0.9)
    @test get_internal_burden(state, "7440-43-9") ≈ 1.9
    @test burden2.maintenance == 0.0

    # Month 3
    burden3 = ecotox_records_to_deb_burden_stateful!(state, concentrations, [record]; retention=0.9)
    @test get_internal_burden(state, "7440-43-9") ≈ 2.71
    @test burden3.maintenance ≈ (2.71 - 2.0) / (4.0 - 2.0)

    # Decay test
    concentrations_zero = Dict("7440-43-9" => 0.0)
    B_before = get_internal_burden(state, "7440-43-9")
    burden_decay = ecotox_records_to_deb_burden_stateful!(state, concentrations_zero, [record]; retention=0.9)
    B_after = get_internal_burden(state, "7440-43-9")
    @test B_after < B_before
    @test B_after > 0.0

    # Non-accumulative equivalence
    state_salt = EcotoxExposureState()
    record_salt = deepcopy(record)
    record_salt["cas"] = "7647145"
    record_salt["cas_norm"] = "7647145"
    record_salt["cas_hyphenated"] = "7647-14-5"

    burden_stateful = ecotox_records_to_deb_burden_stateful!(
        state_salt,
        Dict("7647-14-5" => 10.0),
        [record_salt];
        retention=0.0
    )

    burden_stateless = ecotox_records_to_deb_burden(
        Dict("7647-14-5" => 10.0),
        [record_salt]
    )

    @test burden_stateful.maintenance ≈ burden_stateless.maintenance

    # Data-driven retention
    memory = load_compound_memory_library()
    state_data = EcotoxExposureState()
    burden_data = ecotox_records_to_deb_burden_stateful!(
        state_data,
        Dict("7440-43-9" => 10.0),
        [record];
        memory_library=memory
    )
    @test get_internal_burden(state_data, "7440-43-9") ≈ 1.0
end
