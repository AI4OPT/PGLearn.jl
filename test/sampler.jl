using Distributions
using LinearAlgebra
using TOML

function test_glocal()
    d = PGLearn.Glocal(
        Uniform(0.0, 1.0),
        Distributions.MvNormal(zeros(4), Diagonal(ones(4)))
    )

    @test eltype(d) == Float64
    @test length(d) == 4

    @test size(rand(d)) == (4,)
    @test size(rand(d, 1)) == (4, 1)
    @test size(rand(d, 2)) == (4, 2)

    return nothing
end

function test_ScaledLogNormal()
    d = ScaledLogNormal(0.8, 1.2, 0.05 .* ones(3))

    @test length(d) == 3

    @test isa(d, PGLearn.Glocal)
    @test d.d_α == Uniform(0.8, 1.2)
    @test isa(d.d_η, Distributions.MvLogNormal)

    # Sanity checks
    @test_throws ErrorException ScaledLogNormal(0.8, 0.7, ones(3))   # l > u
    @test_throws ErrorException ScaledLogNormal(0.8, 1.2, -ones(3))  # σ < 0

    return nothing
end

function test_ScaledUniform()
    d = ScaledUniform(0.8, 1.2, 0.05 .* ones(5))

    @test length(d) == 5

    @test isa(d, PGLearn.Glocal)
    @test d.d_α == Uniform(0.8, 1.2)
    @test isa(d.d_η, Distributions.Product)

    # Sanity checks
    @test_throws ErrorException ScaledUniform(0.8, 0.7, ones(3))   # l > u
    @test_throws ErrorException ScaledUniform(0.8, 1.2, -ones(3))  # σ < 0

    return nothing
end

function test_LoadScaler()
    data = PGLearn.OPFData(make_basic_network(pglib("pglib_opf_case14_ieee")))

    # ScaledLogNormal
    options = Dict(
        "noise_type" => "ScaledLogNormal",
        "l" => 0.8,
        "u" => 1.2,
        "sigma" => 0.05,
    )
    ls = LoadScaler(data, options)
    @test ls.d.d_α == Uniform(0.8, 1.2)
    @test isa(ls.d.d_η, MvLogNormal)
    @test length(ls.d.d_η) == 2*data.L
    @test ls.pd_ref == data.pd
    @test ls.qd_ref == data.qd

    # ScaledUniform
    options = Dict(
        "noise_type" => "ScaledUniform",
        "l" => 0.7,
        "u" => 1.5,
        "sigma" => 0.05,
    )
    ls = LoadScaler(data, options)
    @test ls.d.d_α == Uniform(0.7, 1.5)
    @test isa(ls.d.d_η, Distributions.Product)
    @test length(ls.d.d_η) == 2*data.L
    @test ls.pd_ref == data.pd
    @test ls.qd_ref == data.qd

    # Different noise level for buses
    options = Dict(
        "noise_type" => "ScaledUniform",
        "l" => 0.8,
        "u" => 1.2,
        "sigma" => [i/100 for i in 1:data.L],
    )
    ls = LoadScaler(data, options)
    @test ls.d.d_α == Uniform(0.8, 1.2)
    @test isa(ls.d.d_η, Distributions.Product)
    @test length(ls.d.d_η) == 2*data.L
    # check local distributions
    @test Distributions.mean(ls.d.d_η) ≈ ones(2*data.L)
    V = Distributions.var(ls.d.d_η)
    @test V[1:data.L] ≈ [(2*i/100)^2 / 12 for i in 1:data.L]
    @test V[(data.L+1):(2*data.L)] ≈ [(2*i/100)^2 / 12 for i in 1:data.L]
    @test ls.pd_ref == data.pd
    @test ls.qd_ref == data.qd

    return nothing
end

function test_LoadScaler_sanity_checks()
    data = PGLearn.OPFData(make_basic_network(pglib("pglib_opf_case14_ieee")))

    # Invalid noise type
    options = Dict()
    @test_throws ErrorException LoadScaler(data, options)  # missing key
    options["noise_type"] = "InvalidNoiseType"
    @test_throws ErrorException LoadScaler(data, options)  # key exists but bad value

    # Missing or invalid global parameters
    options["noise_type"] = "ScaledLogNormal"
    options["l"] = 0.8
    options["u"] = 1.2
    for v in [Inf, NaN, missing, nothing, 1+im]
        options["l"] = v
        @test_throws ErrorException LoadScaler(data, options)  # bad `l`
        options["l"] = 0.8

        options["u"] = v
        @test_throws ErrorException LoadScaler(data, options)  # bad `u`
        options["u"] = 1.2
    end
    options["l"] = 1.3
    @test_throws ErrorException LoadScaler(data, options)  # l > u
    options["l"] = 0.8

    # Invalid sigma values
    for σ in [Inf, -1, im, ["1", "2"], ones(2, 2), "0.05", [0.05, Inf]]
        options["sigma"] = σ
        @test_throws ErrorException LoadScaler(data, options)
    end

    return nothing
end

function test_sampler()
    data = PGLearn.OPFData(make_basic_network(pglib("pglib_opf_case14_ieee")))
    _data = deepcopy(data)  # keep a deepcopy nearby
    sampler_config = Dict(
        "load" => Dict(
            "noise_type" => "ScaledLogNormal",
            "l" => 0.8,
            "u" => 1.2,
            "sigma" => 0.05,        
        )
    )
    
    opf_sampler  = SimpleOPFSampler(data, sampler_config)
    data1 = rand(MersenneTwister(42), opf_sampler)

    # No side-effect checks
    @test data == _data   # initial data should not have been modified
    @test data !== data1  # new data should be a different dictionary

    # Same RNG and seed should give the same data
    data2 = rand(MersenneTwister(42), opf_sampler)
    @test data2 == data1

    return nothing
end

function test_nminus1_sampler()
    data = PGLearn.OPFData(make_basic_network(pglib("pglib_opf_case14_ieee")))
    sampler_config = Dict{String,Any}(
        "load" => Dict(
            "noise_type" => "ScaledLogNormal",
            "l" => 0.8,
            "u" => 1.2,
            "sigma" => 0.05,        
        ),
        "status" => Dict(
            "type"=> "NMinus1",
        )
    )

    opf_sampler = SimpleOPFSampler(data, sampler_config)

    data1 = rand(MersenneTwister(42), opf_sampler)

    # Exactly one generator or branch should be disabled
    G, E = data.G, data.E
    @test sum(data.gen_status) + sum(data.branch_status) == (G + E)  # original data
    @test sum(data1.gen_status) + sum(data1.branch_status) == (G + E - 1)  # N-1

    # Same RNG and seed should give the same data
    data2 = rand(MersenneTwister(42), opf_sampler)
    @test data2 == data1

    # Unsupporting config should error
    sampler_config["status"]["type"] = "error"
    @test_throws ErrorException SimpleOPFSampler(data, sampler_config)

    return nothing
end

function test_inplace_sampler()
    data = PGLearn.OPFData(make_basic_network(pglib("pglib_opf_case14_ieee")))
    sampler_config = Dict(
        "load" => Dict(
            "noise_type" => "ScaledLogNormal",
            "l" => 0.8,
            "u" => 1.2,
            "sigma" => 0.05,        
        )
    )

    rng = MersenneTwister(42)
    opf_sampler  = SimpleOPFSampler(data, sampler_config)
    rand!(rng, opf_sampler, data)

    return nothing
end

function test_sampler_script()
    sampler_script = joinpath(@__DIR__, "..", "exp", "sampler.jl")
    temp_dir = mktempdir()
    config = Dict(
        "pglib_case" => "pglib_opf_case14_ieee",
        "export_dir" => temp_dir,
        "floating_point_type" => "Float64",
        "sampler" => Dict(
            "compute_clique_decomposition" => true,
            "load" => Dict(
                "noise_type" => "ScaledLogNormal",
                "l" => 0.6,
                "u" => 0.8,
                "sigma" => 0.05,
            ),
            "reserve" => Dict( # tiny reserve requirement
                "type" => "E2ELR",
                "l" => 0.0,
                "u" => 0.1,
                "factor" => 5.0,
            )
        ),
        "OPF" => Dict(
            "DCOPF" => Dict(
                "type" => "DCOPF",
                "solver" => Dict(
                    "name" => "Clarabel",
                )
            ),
            "ACOPF" => Dict(
                "type" => "ACOPF",
                "solver" => Dict(
                    "name" => "Ipopt",
                    "attributes" => Dict(
                        "tol" => 1e-6,
                    )
                )
            ),
            "ED" => Dict(
                "type" => "EconomicDispatch",
                "solver" => Dict(
                    "name" => "Clarabel",
                )
            ),
            "SOCOPF" => Dict(
                "type" => "SOCOPF",
                "solver" => Dict(
                    "name" => "Clarabel",
                )
            ),
            "SOCOPF128" => Dict(
                "type" => "SOCOPF",
                "solver" => Dict(
                    "name" => "Clarabel128",
                    "attributes" => Dict(
                        "max_iter" => 2000,
                        "max_step_fraction" => 0.995,
                        "equilibrate_enable" => true,
                        "tol_gap_abs" => 1e-12,
                        "tol_gap_rel" => 1e-12,
                        "tol_feas" => 1e-12,
                        "tol_infeas_rel" => 1e-12,
                        "tol_ktratio" => 1e-10,
                        "reduced_tol_gap_abs" => 1e-8,
                        "reduced_tol_gap_rel" => 1e-8,
                        "reduced_tol_feas" => 1e-8,
                        "reduced_tol_infeas_abs" => 1e-8,
                        "reduced_tol_infeas_rel" => 1e-8,
                        "reduced_tol_ktratio" => 1e-7,
                        "static_regularization_enable" => false,
                        "dynamic_regularization_enable" => true,
                        "dynamic_regularization_eps" => 1e-28,
                        "dynamic_regularization_delta" => 1e-14,
                        "iterative_refinement_reltol" => 1e-18,
                        "iterative_refinement_abstol" => 1e-18,
                    )
                )
            ),
            "SDPOPF" => Dict(
                "type" => "SDPOPF",
                "solver" => Dict(
                    "name" => "Clarabel",
                    "attributes" => Dict(
                        "static_regularization_constant" => 1e-7,
                    )
                )
            ),
            "SparseSDPOPF" => Dict(
                "type" => "SparseSDPOPF",
                "solver" => Dict(
                    "name" => "Clarabel",
                    "attributes" => Dict(
                        "static_regularization_constant" => 1e-7,
                        "chordal_decomposition_enable" => false
                    )
                )
            )
        )
    )

    config_file = joinpath(temp_dir, "config.toml")
    open(config_file, "w") do io
        TOML.print(io, config)
    end

    case_file, case_name = PGLearn._get_case_info(config)
    smin, smax = 1, 4
    
    # setup env for exp
    @test success(run(setenv(`$(joinpath(Sys.BINDIR, "julia")) --project=. -e "using Pkg; Pkg.instantiate()"`, dir=joinpath(@__DIR__, ".."))))
    
    # run sampler script
    proc = run(setenv(`$(joinpath(Sys.BINDIR, "julia")) --project=. $sampler_script $config_file $smin $smax`, dir=joinpath(@__DIR__, "..")))

    @test success(proc)

    OPFs = collect(keys(config["OPF"]))

    h5_dir = joinpath(@__DIR__, "..", config["export_dir"], "res_h5")

    @test isdir(h5_dir)

    input_file_path = joinpath(h5_dir, "$(case_name)_input_s$smin-s$smax.h5")
    @test isfile(input_file_path)
    # Check that input data file is structured as expected
    h5open(input_file_path, "r") do h5
        @test haskey(h5, "seed")
        @test size(h5["seed"]) == (4,)
        @test eltype(h5["seed"]) == Int

        @test haskey(h5, "pd")
        @test size(h5["pd"]) == (11, 4)
        @test eltype(h5["pd"]) == Float64
        @test haskey(h5, "qd")
        @test size(h5["pd"]) == (11, 4)
        @test eltype(h5["qd"]) == Float64

        @test haskey(h5, "branch_status")
        @test size(h5["branch_status"]) == (20, 4)
        @test eltype(h5["branch_status"]) == Bool
        @test haskey(h5, "gen_status")
        @test size(h5["gen_status"]) == (5, 4)
        @test eltype(h5["gen_status"]) == Bool

        @test haskey(h5, "reserve_requirement")
        @test size(h5["reserve_requirement"]) == (4,)
        @test eltype(h5["reserve_requirement"]) == Float64
        @test haskey(h5, "rmin")
        @test size(h5["rmin"]) == (5,4)
        @test eltype(h5["rmin"]) == Float64
        @test haskey(h5, "rmax")
        @test size(h5["rmax"]) == (5,4)
        @test eltype(h5["rmax"]) == Float64
    end

    h5_paths = [
        joinpath(h5_dir, "$(case_name)_$(opf)_s$smin-s$smax.h5")
        for opf in OPFs
    ]
    @test all(isfile.(h5_paths))
    
    for h5_path in h5_paths
        h5open(h5_path, "r") do h5
            @test haskey(h5, "meta")
            @test haskey(h5, "primal")
            @test haskey(h5, "dual")
            n_seed = length(h5["meta"]["seed"])
            @test n_seed == smax - smin + 1
            for i in 1:n_seed
                @test h5["meta"]["termination_status"][i] ∈ ["OPTIMAL", "ALMOST_OPTIMAL", "LOCALLY_SOLVED", "ALMOST_LOCALLY_SOLVED"]
                @test h5["meta"]["primal_status"][i] ∈ ["FEASIBLE_POINT", "NEARLY_FEASIBLE_POINT"]
                @test h5["meta"]["dual_status"][i] ∈ ["FEASIBLE_POINT", "NEARLY_FEASIBLE_POINT"]
                @test h5["meta"]["seed"][i] == smin + i - 1
            end
        end
    end
end

@testset "Sampler" begin
    @testset test_glocal()
    @testset test_ScaledLogNormal()
    @testset test_ScaledUniform()
    @testset test_LoadScaler()
    @testset test_LoadScaler_sanity_checks()
    @testset test_sampler()
    @testset test_nminus1_sampler()
    @testset test_inplace_sampler()
    @testset test_sampler_script()
end