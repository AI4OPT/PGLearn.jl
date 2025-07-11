using Random
using Test

using HDF5
using Quadmath

using Clarabel
using Ipopt
using HiGHS
using JuMP

using PowerModels
const PM = PowerModels
PM.silence()
using PGLib

using PGLearn

const IPOPT_SOLVER = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "linear_solver" => "mumps", "print_level" => 1, "tol" => 1e-6)
const CLRBL_SOLVER = JuMP.optimizer_with_attributes(Clarabel.Optimizer, "verbose" => true)
const CLRBL_SOLVER_SDP = JuMP.optimizer_with_attributes(Clarabel.Optimizer, "verbose" => true, "static_regularization_constant" => 1e-7, "chordal_decomposition_enable" => false)
const CLRBL128_SOLVER = JuMP.optimizer_with_attributes(Clarabel.Optimizer{Float128}, "verbose" => true)
const HIGHS_SOLVER = JuMP.optimizer_with_attributes(HiGHS.Optimizer, "output_flag" => true)

const OPT_SOLVERS = Dict(
    PGLearn.ACOPF             => IPOPT_SOLVER,
    PGLearn.SOCOPFQuad        => IPOPT_SOLVER,
    PGLearn.SOCOPF            => CLRBL_SOLVER,
    PGLearn.SDPOPF            => CLRBL_SOLVER_SDP,
    PGLearn.SparseSDPOPF      => CLRBL_SOLVER_SDP,
    PGLearn.DCOPF             => HIGHS_SOLVER,
    PGLearn.EconomicDispatch  => HIGHS_SOLVER,
)


const PGLIB_CASES = ["14_ieee", "30_ieee", "57_ieee", "89_pegase", "118_ieee"]
const PGLIB_CASES_SDP = ["14_ieee"]
@testset "PGLearn" begin
    include("utils/io.jl")
    include("utils/float.jl")    
    include("bridges.jl")
    include("sampler.jl")
    include("opf/opf.jl")
    include("postprocess.jl")
end