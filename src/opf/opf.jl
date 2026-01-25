using SparseArrays

abstract type AbstractFormulation end


mutable struct OPFData
    case::String
    base_mva::Float64

    N::Int  # number of buses
    E::Int  # number of branches
    G::Int  # number of generators
    L::Int  # number of loads
    
    A::SparseMatrixCSC{Float64,Int}  # from/to branch incidence matrix
    Ag::SparseMatrixCSC{Float64,Int}  # generator incidence matrix

    # Bus data
    vnom::Vector{Float64}
    vmin::Vector{Float64}
    vmax::Vector{Float64}
    gs::Vector{Float64}  # shunt
    bs::Vector{Float64}  # shunt
    bus_arcs_fr::Vector{Vector{Int}}  # indices of branches exiting the bus
    bus_arcs_to::Vector{Vector{Int}}  # indices of branches entering the bus
    bus_gens::Vector{Vector{Int}}  # indices of generators at the bus
    bus_loads::Vector{Vector{Int}}
    ref_bus::Int  # index of slack bus

    # Load data
    pd::Vector{Float64}  # active power demand
    qd::Vector{Float64}  # reactive power demand

    # Generator data
    pgmin::Vector{Float64}
    pgmax::Vector{Float64}
    qgmin::Vector{Float64}
    qgmax::Vector{Float64}
    c0::Vector{Float64}  # constant cost
    c1::Vector{Float64}  # linear cost
    c2::Vector{Float64}  # quadratic cost
    gen_status::Vector{Bool}  # generator status
    rmin::Vector{Float64}  # reserve lower bound
    rmax::Vector{Float64}  # reserve upper bound

    reserve_requirement::Float64  # minimum reserve requirement

    # Branch data
    bus_fr::Vector{Int}  # from bus
    bus_to::Vector{Int}  # to bus
    g::Vector{Float64}
    b::Vector{Float64}
    gff::Vector{Float64}
    gft::Vector{Float64}
    gtf::Vector{Float64}
    gtt::Vector{Float64}
    bff::Vector{Float64}
    bft::Vector{Float64}
    btf::Vector{Float64}
    btt::Vector{Float64}
    smax::Vector{Float64}  # max apparent power flow
    dvamin::Vector{Float64}  # angle difference min
    dvamax::Vector{Float64}  # angle difference max
    branch_status::Vector{Bool}  # branch status

    # Clique decomposition data for SparseSDPOPF
    clique_decomposition::Union{Nothing,Vector{Vector{Int}}}
end

"""
    OPFData(network::Dict{String,Any}; compute_clique_decomposition::Bool=false)

Convert a PowerModels data dictionary to `OPFData` structure.

The PowerModels data dictionary must be in basic format.
"""
function OPFData(network::Dict{String,Any}; compute_clique_decomposition::Bool=false)
    @assert network["basic_network"] "Network data is not in basic format."
    @assert network["per_unit"] == true "Network data is not per-unit scaled."

    N = length(network["bus"])
    E = length(network["branch"])
    G = length(network["gen"])
    L = length(network["load"])

    # Bus data
    vnom = [network["bus"]["$i"]["base_kv"] for i in 1:N]
    vmin = [network["bus"]["$i"]["vmin"] for i in 1:N]
    vmax = [network["bus"]["$i"]["vmax"] for i in 1:N]

    # Aggregate shunts at the bus level
    gs = zeros(Float64, N)
    bs = zeros(Float64, N)
    for (_, shunt) in network["shunt"]
        shunt["status"] == 1 || continue  # skip inactive buses
        i = shunt["shunt_bus"]
        gs[i] += shunt["gs"]
        bs[i] += shunt["bs"]
    end

    # Reference bus
    ref_buses = [i for i in 1:N if network["bus"]["$i"]["bus_type"] == 3]
    @assert length(ref_buses) == 1 "There must be exactly one reference bus"
    ref_bus = ref_buses[1]

    # Load data
    pd = zeros(Float64, L)
    qd = zeros(Float64, L)
    bus_loads = [Int[] for _ in 1:N]
    for l in 1:L
        load = network["load"]["$l"]
        i = load["load_bus"]
        push!(bus_loads[i], l)
        pd[l] = (load["status"] == 1) * load["pd"]
        qd[l] = (load["status"] == 1) * load["qd"]
    end
    sort!.(bus_loads)

    # Generator data
    pgmin = zeros(Float64, G)
    pgmax = zeros(Float64, G)
    qgmin = zeros(Float64, G)
    qgmax = zeros(Float64, G)
    c0 = zeros(Float64, G)
    c1 = zeros(Float64, G)
    c2 = zeros(Float64, G)
    gen_status = zeros(Bool, G)
    bus_gens = [Int[] for _ in 1:N]
    Ag_i = zeros(Int, G)
    Ag_j = zeros(Int, G)
    Ag_v = zeros(Float64, G)
    rmin = zeros(Float64, G)
    rmax = zeros(Float64, G)
    for g in 1:G
        gen = network["gen"]["$g"]

        i = gen["gen_bus"]
        push!(bus_gens[i], g)

        pgmin[g] = gen["pmin"]
        pgmax[g] = gen["pmax"]
        qgmin[g] = gen["qmin"]
        qgmax[g] = gen["qmax"]
        # ⚠️ cost data assumes quadratic cost everywhere
        c0[g] = gen["cost"][3]
        c1[g] = gen["cost"][2]
        c2[g] = gen["cost"][1]

        gen_status[g] = (gen["gen_status"] == 1)

        # Generator incidence matrix
        Ag_i[g] = gen["gen_bus"]
        Ag_j[g] = g
        Ag_v[g] = gen["gen_status"]

        # Reserve bounds
        rmin[g] = get(gen, "rmin", 0.0)
        rmax[g] = get(gen, "rmax", 0.0)
    end
    # sort everything again
    sort!.(bus_gens)

    Ag = sparse(Ag_i, Ag_j, Ag_v, N, G)

    # Minimum reserve requirement
    reserve_requirement = get(network, "reserve_requirement", 0.0)

    # Branch data
    bus_fr = zeros(Int, E)
    bus_to = zeros(Int, E)
    branch_g = zeros(Float64, E)
    branch_b = zeros(Float64, E)
    gff = zeros(Float64, E)
    gft = zeros(Float64, E)
    gtf = zeros(Float64, E)
    gtt = zeros(Float64, E)
    bff = zeros(Float64, E)
    bft = zeros(Float64, E)
    btf = zeros(Float64, E)
    btt = zeros(Float64, E)
    smax = zeros(Float64, E)
    dvamin = zeros(Float64, E)
    dvamax = zeros(Float64, E)
    branch_status = zeros(Bool, E)
    bus_arcs_fr = [Int[] for _ in 1:N]
    bus_arcs_to = [Int[] for _ in 1:N]
    A_i = zeros(Int, 2*E)
    A_j = zeros(Int, 2*E)
    A_v = zeros(Float64, 2*E)
    for e in 1:E
        branch = network["branch"]["$e"]
        i::Int = branch["f_bus"]
        j::Int = branch["t_bus"]

        bus_fr[e] = i
        bus_to[e] = j
        push!(bus_arcs_fr[i], e)
        push!(bus_arcs_to[j], e)

        z::ComplexF64 = branch["br_r"] + im * branch["br_x"]
        y = inv(z)  # compute branch admittance
        isfinite(y) || error("Branch $e has zero impedance")
        g, b = real(y), imag(y)

        branch_g[e] = g
        branch_b[e] = b

        # Compute tap ratio
        τ::Float64 = get(branch, "tap", 1.0)
        ϕ::Float64 = get(branch, "shift", 0.0)
        tr = τ * cos(ϕ)
        ti = τ * sin(ϕ)
        ttm = abs2(τ)

        g_fr::Float64 = branch["g_fr"]
        b_fr::Float64 = branch["b_fr"]
        g_to::Float64 = branch["g_to"]
        b_to::Float64 = branch["b_to"]

        # The formula below match the PM implementation,
        #  given Ohm's in the form:
        #    gff * wf + gft * wr + bft * wi == pf
        #   -bff * wf - bft * wr + gft * wi == qf
        #    gtt * wt + gtf * wr - btf * wi == pt
        #   -btt * wt - btf * wr - gtf * wi == qt
        # Note: OPOMO and Gurobi optimods use the same form of Ohm's law
        gff[e] = (g+g_fr) / ttm
        gft[e] = (-g*tr+b*ti) / ttm
        gtf[e] = (-g*tr-b*ti) / ttm
        gtt[e] = (g+g_to)
        bff[e] = (b+b_fr) / ttm
        bft[e] = (-b*tr-g*ti) / ttm
        btf[e] = (-b*tr+g*ti) / ttm
        btt[e] = (b+b_to)

        # Angle deviation constraints
        dvamin[e] = branch["angmin"]
        dvamax[e] = branch["angmax"]

        # Thermal limits
        smax[e] = branch["rate_a"]

        branch_status[e] = branch["br_status"] == 1

        # Branch incidence matrix
        A_i[e] = e
        A_j[e] = i
        A_v[e] = 1

        A_i[E+e] = e
        A_j[E+e] = j
        A_v[E+e] = -1
    end
    sort!.(bus_arcs_fr)
    sort!.(bus_arcs_to)

    A = sparse(A_i, A_j, A_v, E, N)

    if compute_clique_decomposition
        clique_decomposition = _get_clique_decomposition(network)
    else 
        clique_decomposition = nothing
    end

    return OPFData(
        network["name"], network["baseMVA"],
        N, E, G, L, A, Ag,
        vnom, vmin, vmax, gs, bs,
        bus_arcs_fr, bus_arcs_to, bus_gens, bus_loads, ref_bus,
        pd, qd,
        pgmin, pgmax,
        qgmin, qgmax,
        c0, c1, c2,
        gen_status,
        rmin, rmax, reserve_requirement,
        bus_fr, bus_to,
        branch_g, branch_b,
        gff, gft, gtf, gtt,
        bff, bft, btf, btt,
        smax, dvamin, dvamax,
        branch_status,
        clique_decomposition,
    )
end

"""
    _get_clique_decomposition(network::Dict{String,Any})

Compute the clique decomposition of a PowerModels data dictionary.
The output is a `Vector{Vector{Int}}` containing all the cliques in the chordal completion of network.
"""
function _get_clique_decomposition(network::Dict{String,Any})
    return instantiate_model(network, SparseSDPWRMPowerModel, PM.build_opf).ext[:SDconstraintDecomposition].decomp
end

function to_dict(data::OPFData)
    d = Dict{String,Any}()
    for field in fieldnames(OPFData)
        v = getfield(data, field)
        if isa(v, SparseMatrixCSC)
            I, J, V = findnz(v)
            M, N = size(v)
            d[string(field)] = Dict("I" => I, "J" => J, "V" => V, "M" => M, "N" => N)
        elseif isa(v, AbstractArray)
            d[string(field)] = copy(v)
        else
            d[string(field)] = v
        end
    end
    return d
end

# use == instead of === on fields. see julia#4648
Base.:(==)(a::OPFData, b::OPFData) = all(getfield(a, field) == getfield(b, field) for field in fieldnames(OPFData))


mutable struct OPFModel{OPF <: AbstractFormulation}
    data::OPFData
    model::JuMP.GenericModel
end

include("utils.jl")
include("ptdf.jl")

include("acp.jl")      # ACPPowerModel
include("dcp.jl")      # DCPPowerModel
include("ed.jl")       # EconomicDispatch
include("socwr.jl")    # SOCOPF & SOCOPFQuad
include("sdpwrm.jl")   # SDPWRMPowerModel
include("sparse_sdpwrm.jl")

# Contains a list of all supported OPF models
const SUPPORTED_OPF_MODELS = Type{<:AbstractFormulation}[
    ACOPF,
    DCOPF,
    EconomicDispatch,
    SOCOPFQuad,
    SOCOPF,
    SDPOPF,
    SparseSDPOPF,
]

# A name --> type dictionary
# Used for passing the OPF type as a string (e.g. through config file)
const OPF2TYPE = Dict{String,Type{<:AbstractFormulation}}(
    "ACOPF" => ACOPF,
    "DCOPF" => DCOPF,
    "EconomicDispatch" => EconomicDispatch,
    "SOCOPFQuad" => SOCOPFQuad,
    "SOCOPF" => SOCOPF,
    "SDPOPF" => SDPOPF,
    "SparseSDPOPF" => SparseSDPOPF,
)

function build_opf(OPF::Type{<:AbstractFormulation}, network::Dict, optimizer; compute_clique_decomposition::Bool=false, kwargs...)
    return build_opf(OPF, OPFData(network; compute_clique_decomposition), optimizer; kwargs...)
end

function solve!(opf::OPFModel{<:AbstractFormulation})
    maybe_set_ad(opf)
    optimize!(opf.model)
end
function maybe_set_ad(opf::OPFModel{<:AbstractFormulation})
    !isnothing(nonlinear_model(opf.model)) && set_attribute(
        opf.model,
        MOI.AutomaticDifferentiationBackend(),
        MOI.Nonlinear.SymbolicMode(),
    )
end