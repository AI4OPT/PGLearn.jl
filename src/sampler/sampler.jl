# common distributions
include("glocal.jl")

abstract type AbstractOPFSampler end

function Random.rand(::AbstractRNG, ::AbstractOPFSampler)
    error("`rand` function not implemented for $(typeof(s)).")
end

struct SimpleOPFSampler{LS,RS,SS} <: AbstractOPFSampler
    data::OPFData
    load_sampler::LS
    reserve_sampler::RS
    status_sampler::SS
end

function SimpleOPFSampler(data::OPFData, config::Dict)
    data = deepcopy(data)

    # Instantiate load sampler
    load_sampler = LoadScaler(data, config["load"])

    # Instantiate reserve sampler
    get!(config, "reserve", Dict())
    reserve_sampler = ReserveScaler(data, config["reserve"])

    # Instantiate status sampler
    get!(config, "status", Dict())
    status_sampler = StatusSampler(data, config["status"])

    return SimpleOPFSampler(data, load_sampler, reserve_sampler, status_sampler)
end

function Random.rand(rng::AbstractRNG, opf_sampler::SimpleOPFSampler)
    data = deepcopy(opf_sampler.data)
    rand!(rng, opf_sampler, data)
end

function Random.rand(seed::Int, opf_sampler::SimpleOPFSampler)
    rng = MersenneTwister(seed)
    data = deepcopy(opf_sampler.data)
    rand!(rng, opf_sampler, data)
end

"""
    rand!(rng::AbstractRNG, s::SimpleOPFSampler, data::Dict)

Sample one new OPF instance and modify `data` in-place.

`data` must be a `Dict` in PowerModels format, representing the same network
    (i.e., same grid components with same indexing) as the one used to create `s`.
"""
function Random.rand!(rng::AbstractRNG, s::SimpleOPFSampler, data::OPFData)
    pd, qd = rand(rng, s.load_sampler)
    _set_loads!(data, pd, qd)

    MRR, rmin, rmax = rand(rng, s.reserve_sampler)
    _set_reserve!(data, MRR, rmin, rmax)

    br_status, gen_status = rand(rng, s.status_sampler)
    _set_status!(data, br_status, gen_status)

    return data
end

function _set_loads!(data, pd, qd)
    L = data.L
    length(pd) == L || throw(DimensionMismatch())
    length(qd) == L || throw(DimensionMismatch())
    
    data.pd .= pd
    data.qd .= qd

    return nothing
end

function _set_reserve!(data, MRR, rmin, rmax)
    G = data.G
    length(rmin) == G || throw(DimensionMismatch())
    length(rmax) == G || throw(DimensionMismatch())

    data.rmin .= rmin
    data.rmax .= rmax
    data.reserve_requirement = MRR

    return nothing
end

function _set_status!(data, br_status, gen_status)
    length(br_status) == data.E || throw(DimensionMismatch())
    length(gen_status) == data.G || throw(DimensionMismatch())

    data.branch_status .= br_status
    data.gen_status .= gen_status

    return nothing
end

include("load.jl")
include("reserve.jl")
include("status.jl")