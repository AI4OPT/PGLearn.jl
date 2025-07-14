using LinearAlgebra
using SparseArrays

function test_opf_pm(::Type{PGLearn.SparseSDPOPF}, data::Dict)
    OPF = PGLearn.SparseSDPOPF

    data["basic_network"] || error("Input data must be in basic format to test")
    N = length(data["bus"])
    E = length(data["branch"])
    G = length(data["gen"])

    # Solve OPF with PowerModels
    solver = OPT_SOLVERS[OPF]
    res_pm = PM.solve_opf(data, PM.SparseSDPWRMPowerModel, solver)

    # Build and solve OPF with PGLearn
    solver = OPT_SOLVERS[OPF]
    opf = PGLearn.build_opf(OPF, data, solver; compute_clique_decomposition=true)
    PGLearn.solve!(opf)
    res = PGLearn.extract_result(opf)

    # Check that the right problem was indeed solved
    @test res["meta"]["formulation"] == string(OPF)
    @test res["meta"]["termination_status"] ∈ ["OPTIMAL", "ALMOST_OPTIMAL"]
    @test res["meta"]["primal_status"] ∈ ["FEASIBLE_POINT", "NEARLY_FEASIBLE_POINT"]
    @test res["meta"]["dual_status"] ∈ ["FEASIBLE_POINT", "NEARLY_FEASIBLE_POINT"]
    # Check objective value against PowerModels
    @test isapprox(res["meta"]["primal_objective_value"], res_pm["objective"], atol=1e-4, rtol=1e-4)
    @test isapprox(res["meta"]["primal_objective_value"], res["meta"]["dual_objective_value"], rtol=1e-6)

    # Force PM solution into our model, and check that the solution is feasible
    # TODO: use JuMP.primal_feasibility_report instead
    #    (would require extracting a variable => value Dict)
    sol_pm = res_pm["solution"]
    var2val_pm = Dict(
        :pg => Float64[
            get(get(sol_pm["gen"], "$g", Dict()), "pg", 0) for g in 1:G
        ],
        :qg => Float64[
            get(get(sol_pm["gen"], "$g", Dict()), "qg", 0) for g in 1:G
        ],
        :w => Float64[sol_pm["bus"]["$i"]["w"] for i in 1:N]
    )
    model = opf.model
    @constraint(model, model[:pg] .== var2val_pm[:pg])
    @constraint(model, model[:qg] .== var2val_pm[:qg])
    @constraint(model, model[:w] .== var2val_pm[:w])

    optimize!(model)
    @test termination_status(model) ∈ [OPTIMAL, ALMOST_OPTIMAL]
    @test primal_status(model) ∈ [FEASIBLE_POINT, NEARLY_FEASIBLE_POINT]
    # Also check that we get the same objective value as PowerModels
    @test isapprox(objective_value(opf.model), res_pm["objective"], atol=1e-3, rtol=1e-3)

    return nothing
end