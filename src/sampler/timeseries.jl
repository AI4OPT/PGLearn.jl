"""
    TimeSeriesOPFSampler

An OPFSampler that grabs input data from an H5 file.
"""
struct TimeSeriesOPFSampler{T} <: AbstractOPFSampler
    data::OPFData
    pd::Matrix{T}
    qd::Matrix{T}
    MRR::Vector{T}
    rmin::Matrix{T}
    rmax::Matrix{T}
    br_status::Matrix{Bool}
    gen_status::Matrix{Bool}
end

function TimeSeriesOPFSampler(data::OPFData, config::Dict)
    data = deepcopy(data)

    h5_path = get(config, "h5_path", "")
    h5 = load_h5(h5_path, "/")

    pd = get(h5, "pd", nothing)
    qd = get(h5, "qd", nothing)

    MRR = get(h5, "MRR", nothing)
    rmin = get(h5, "rmin", nothing)
    rmax = get(h5, "rmax", nothing)

    br_status = get(h5, "br_status", nothing)
    gen_status = get(h5, "gen_status", nothing)

    return TimeSeriesOPFSampler(data, pd, qd, MRR, rmin, rmax, br_status, gen_status)
end


function Random.rand(seed::Int, opf_sampler::TimeSeriesOPFSampler)
    data = deepcopy(opf_sampler.data)
    rand!(seed, opf_sampler, data)
end

function Random.rand!(seed::Int, opf_sampler::TimeSeriesOPFSampler, data::OPFData)
    pd = notnothing(pd) ? opf_sampler.pd[:, seed] : opf_sampler.data.pd
    qd = notnothing(qd) ? opf_sampler.qd[:, seed] : opf_sampler.data.qd
    _set_loads!(data, pd, qd)

    MRR = notnothing(opf_sampler.MRR) ? opf_sampler.MRR[seed] : opf_sampler.data.reserve_requirement
    rmin = notnothing(opf_sampler.rmin) ? opf_sampler.rmin[:, seed] : opf_sampler.data.rmin
    rmax = notnothing(opf_sampler.rmax) ? opf_sampler.rmax[:, seed] : opf_sampler.data.rmax
    _set_reserve!(data, MRR, rmin, rmax)

    br_status = notnothing(opf_sampler.br_status) ? opf_sampler.br_status[:, seed] : opf_sampler.data.branch_status
    gen_status = notnothing(opf_sampler.gen_status) ? opf_sampler.gen_status[:, seed] : opf_sampler.data.gen_status
    _set_status!(data, br_status, gen_status)

    return data
end
