using PGLearn: OPFData

include("utils.jl")
include("ptdf.jl")

function test_opf_pm(OPF::Type{<:PGLearn.AbstractFormulation}, casename::String)
    network = make_basic_network(pglib(casename))
    @testset "Full data" begin test_opf_pm(OPF, network) end

    network_drop = deepcopy(network)

    is_bridge = PGLearn.bridges(PGLearn.OPFData(network_drop))
    non_bridge = [network_drop["branch"]["$e"] for (e, b) in enumerate(is_bridge) if !b]
    drop_branch = argmin(branch->branch["rate_a"], non_bridge)["index"]
    network_drop["branch"]["$drop_branch"]["br_status"] = 0

    drop_gen = argmin(gen->gen[2]["pmax"], network_drop["gen"])[1]
    network_drop["gen"]["$drop_gen"]["gen_status"] = 0

    if OPF == PGLearn.EconomicDispatch
        @test_throws ErrorException test_opf_pm(OPF, network_drop) # ED does not yet support branch status
    else
        @testset "Branch/Gen Status" begin test_opf_pm(OPF, network_drop) end
    end
end

"""
    test_opf_pm(OPF, data)

Build & solve opf using PGLearn, and compare against PowerModels implementation.
"""
function test_opf_pm(::Type{OPF}, data::Dict) where{OPF <: PGLearn.AbstractFormulation}
    error("""`test_opf_pm($(OPF), data)` not implemented.
    You must implement a function with the following signature:
        function test_opf_pm(::Type{OPF}, data::Dict) where{OPF <: $(OPF)}
            # unit tests ...
            return nothing
        end
    """)

    return nothing
end

include("acp.jl")
include("dcp.jl")
include("socwr.jl")
include("sdpwrm.jl")
include("sparse_sdpwrm.jl")
include("ed.jl")

# other tests
include("quad_obj.jl")

@testset "OPF" begin
    @testset "$(OPF)" for OPF in PGLearn.SUPPORTED_OPF_MODELS
        if OPF in [PGLearn.SDPOPF, PGLearn.SparseSDPOPF]
            cases = PGLIB_CASES_SDP
        else
            cases = PGLIB_CASES
        end
        @testset "$(casename)" for casename in cases
            test_opf_pm(OPF, casename)
        end

        @testset "QuadObj" begin test_quad_obj_warn(OPF) end
    end

    # SOCOPF
    @testset _test_socwr_DualSolFormat()
    @testset _test_socwr_DualFeasibility()

    # SDPOPF
    @testset _test_sdpwrm_DualSolFormat(PGLearn.SDPOPF)
    @testset _test_sdpwrm_DualFeasibility(PGLearn.SDPOPF)

    # SparseSDPOPF
    @testset _test_sdpwrm_DualSolFormat(PGLearn.SparseSDPOPF)
    @testset _test_sdpwrm_DualFeasibility(PGLearn.SparseSDPOPF)
end

@testset "OPFData" begin
    test_voltage_phasor_bounds_scalar()
    @testset "$casename" for casename in PGLIB_CASES
        network = make_basic_network(pglib(casename))
        data = OPFData(network)
        test_opfdata(data, network)
        test_voltage_phasor_bounds(data, network)
        @testset "PTDF" begin
            test_ptdf(network, data)
        end
    end
end
