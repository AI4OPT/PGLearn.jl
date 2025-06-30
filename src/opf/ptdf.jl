using SparseArrays
using PowerModels
using LinearAlgebra

abstract type AbstractPTDF end


struct LazyPTDF{FT} <: AbstractPTDF
    N::Int  # number of buses
    E::Int  # number of branches
    islack::Int  # Index of slack bus

    A::SparseMatrixCSC{Float64,Int}  # incidence matrix
    b::Vector{Float64}  # branch susceptances
    BA::SparseMatrixCSC{Float64,Int}  # B*A
    AtBA::SparseMatrixCSC{Float64,Int}  # AᵀBA

    F::FT  # Factorization of AᵀBA

    # TODO: cache
end

function LazyPTDF(data::OPFData)
    N, E, A, b, ref_idx = data.N, data.E, data.A, data.b, data.ref_bus

    all(data.branch_status) || error("LazyPTDF does not support disabled branches.")

    B = Diagonal(b)
    BA = B * A
    S = AtBA = A' * BA

    # Adjust row & column corresponding to slack bus
    S[ref_idx, :] .= 0.0
    S[:, ref_idx] .= 0.0
    S[ref_idx, ref_idx] = 1.0

    F = ldlt(S)

    return LazyPTDF{typeof(F)}(N, E, ref_idx, A, b, BA, AtBA, F)
end

"""
    compute_flow!(pf, pg, Φ::LazyPTDF)

Compute power flow `pf = Φ*pg` lazily, without forming the PTDF matrix.

Namely, `pf` is computed as `pf = BA * (F \\ pg)`, where `F` is an LDLᵀ factorization of `AᵀBA`.
"""
function compute_flow!(pf, pg, Φ::LazyPTDF)
    θ = Φ.F \ pg
    θ[Φ.islack, :] .= 0  # slack voltage angle is zero
    mul!(pf, Φ.BA, θ)
    return pf
end

"""
    ptdf_row(Φ::LazyPTDF, e::Int)

Return the `e`-th row of (lazy) PTDF matrix `Φ`.
"""
function ptdf_row(Φ::LazyPTDF, e::Int)
    1 <= e <= Φ.E || error("Invalid row index: $e (must be between 1 and $E)")

    u = Φ.BA[e, :]
    y = Φ.F \ u
    y[Φ.islack] = 0.0

    return y
end


struct FullPTDF{T} <: AbstractPTDF
    N::Int  # number of buses
    E::Int  # number of branches

    matrix::Matrix{T}  # PTDF matrix
end

function FullPTDF(data::OPFData)
    lazy = LazyPTDF(data)
    return FullPTDF(lazy)
end

function FullPTDF(lazy::LazyPTDF)
    (lazy.E > 4096) && @warn "FullPTDF: large PTDF matrix (E = $(lazy.E)). Consider using LazyPTDF instead."
    T = eltype(lazy.F)
    PTDF = zeros(T, lazy.E, lazy.N)
    for e in 1:lazy.E
        PTDF[e, :] .= ptdf_row(lazy, e)
    end

    return FullPTDF(lazy.N, lazy.E, PTDF)
end

"""
    compute_flow!(pf, pg, Φ::FullPTDF)

Compute power flow `pf = Φ*pg` given PTDF matrix `Φ` and nodal injections `pg`.
"""
function compute_flow!(pf, pg, Φ::FullPTDF)
    mul!(pf, Φ.matrix, pg)
    return pf
end

"""
    ptdf_row(Φ::FullPTDF, e::Int)

Return the `e`-th row of PTDF matrix `Φ`.
"""
function ptdf_row(Φ::FullPTDF, e::Int)
    return Φ.matrix[e, :]
end