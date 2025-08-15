using SparseArrays

struct EconomicDispatch <: AbstractFormulation end

const THERMAL_PENALTY = 150000.0
const MAX_PTDF_ITERATIONS = 128
const MAX_PTDF_PER_ITERATION = 8
const ITERATIVE_PTDF_TOL = 1e-6

function build_opf(::Type{EconomicDispatch}, data::OPFData, optimizer;
    T=Float64,
    soft_thermal_limit::Bool=false,
    thermal_penalty=THERMAL_PENALTY,
    iterative_ptdf_tol=ITERATIVE_PTDF_TOL,
    max_ptdf_iterations=MAX_PTDF_ITERATIONS,
    max_ptdf_per_iteration=MAX_PTDF_PER_ITERATION,
)
    thermal_penalty >= 0.0 || error("$OPF option transmission_penalty must be non-negative")
    max_ptdf_iterations > 0 || error("$OPF option max_ptdf_iterations must be a positive integer")
    max_ptdf_per_iteration > 0 || error("$OPF option max_ptdf_per_iteration must be a positive integer")

    # Grab some data
    E, G = data.E, data.G
    pgmin, pgmax = data.pgmin, data.pgmax
    c0, c1, c2 = data.c0, data.c1, data.c2
    smax = data.smax
    branch_status, gen_status = data.branch_status, data.gen_status

    all(branch_status) || error("EconomicDispatch does not support disabled branches.")

    PD = sum(data.pd) # total demand

    # reserves data
    MRR = data.reserve_requirement
    rmin = data.rmin
    rmax = data.rmax

    model = JuMP.GenericModel{T}(optimizer)
    model.ext[:opf_model] = EconomicDispatch
    model.ext[:solve_metadata] = Dict{Symbol,Any}(
        :iterative_ptdf_tol => iterative_ptdf_tol,
        :max_ptdf_iterations => max_ptdf_iterations,
        :max_ptdf_per_iteration => max_ptdf_per_iteration,
    )

    #
    #   I. Variables
    #

    # Active dispatch
    @variable(model, pg[g in 1:G])
    # Active reserves
    @variable(model, r[g in 1:G])

    # Active branch flow
    @variable(model, pf[e in 1:E])

    # Active branch flow slack
    @variable(model, δf[1:E] >= 0)

    # 
    #   II. Constraints
    #

    # Generation bounds (zero if generator is off)
    set_lower_bound.(pg, gen_status .* pgmin)
    # ⚠️ we do not set upper bounds on pg because they are redundant
    #       with the constraint pg + r <= pgmax

    # Reserve bounds (zero if generator is off)
    set_lower_bound.(r, gen_status .* rmin)
    set_upper_bound.(r, gen_status .* rmax)

    # Flow and flow slack bounds
    active_smax = branch_status .* smax

    # Thermal limits
    # Note that `soft_thermal_limit` is a boolean flag such that 
    #   * if `true`, the penalty variables `δf` appears as expected
    #   * if `false`, the penalty variables `δf` do not appear in the constraint
    @constraint(model, pf_lower_bound[e in 1:E], pf[e] + soft_thermal_limit * δf[e] >= -active_smax[e])
    @constraint(model, pf_upper_bound[e in 1:E], pf[e] - soft_thermal_limit * δf[e] <= active_smax[e])

    # Maximum generator output
    @constraint(model, gen_max_output[g in 1:G], pg[g] + r[g] <= gen_status[g] * pgmax[g])

    # Total reserve requirement
    @constraint(model,
        global_reserve_requirement,
        sum(gen_status[g] * r[g] for g in 1:G) >= MRR
    )

    @constraint(model,
        global_power_balance,
        sum(gen_status[g] * pg[g] for g in 1:G) == PD
    )

    model.ext[:PTDF] = LazyPTDF(data)

    model[:ptdf_flow] = Vector{ConstraintRef}(undef, E)
    model.ext[:tracked_branches] = zeros(Bool, E)
    model.ext[:ptdf_iterations] = 0
    
    #
    #   III. Objective
    #
    l, u = extrema(c2[g] for g in 1:G if gen_status[g])
    (l == u == 0.0) || @warn "Data $(data.case) has quadratic cost terms; those terms are being ignored"

    @objective(model, Min, 
        sum(c1[g] * pg[g] + c0[g] for g in 1:G if gen_status[g])
        + thermal_penalty * sum(δf)
    )

    return OPFModel{EconomicDispatch}(data, model)
end


function solve!(opf::OPFModel{EconomicDispatch}) 
    maybe_set_ad(opf)
    model = opf.model

    data = opf.data

    # Grab some data
    N, E, G, L = data.N, data.E, data.G, data.L
    Ag, smax = data.Ag, data.smax
    tol = model.ext[:solve_metadata][:iterative_ptdf_tol]

    # Get bus-wise pg VariableRef
    pg_bus = Ag * model[:pg]
    
    # Compute load-induced power flows
    pd_nodal = zeros(N)
    for (i, loads) in enumerate(data.bus_loads)
        for l in loads
            pd_nodal[i] += data.pd[l]
        end
    end
    ptdfb = zeros(E)
    compute_flow!(ptdfb, pd_nodal, model.ext[:PTDF]) # ptdfb = PTDF * pd

    # Initialize lazy pf buffer
    pf_ = zeros(E)

    # Initialize metadata
    solved = false
    niter = 0
    solve_time = 0.0
    model.ext[:termination_info] = Dict{Symbol,Any}(
        :termination_status => nothing,
        :primal_status => nothing,
        :dual_status => nothing,
        :solve_time => nothing,
    )
    st = nothing

    # Begin lazy PTDF loop
    t0 = time()
    while !solved && niter < model.ext[:solve_metadata][:max_ptdf_iterations]
        # Solve model
        optimize!(opf.model)
        
        # Exit if not solved optimally
        st = termination_status(model)
        st ∈ (MOI.OPTIMAL, MOI.LOCALLY_SOLVED) || break

        # Get pg and the corresponding pf_
        pg_ = value.(model[:pg])
        p_ = Ag * pg_ - pd_nodal
        compute_flow!(pf_, p_, model.ext[:PTDF])
        
        # Check pf_ for violations
        n_violated = 0
        for e in 1:E
            # Skip check if already tracked
            if (model.ext[:tracked_branches][e]) || (-smax[e] - tol <= pf_[e] <= smax[e] + tol)
                continue
            end
            
            # If violated, add the corresponding constraint
            n_violated += 1
            if n_violated <= model.ext[:solve_metadata][:max_ptdf_per_iteration]
                row = ptdf_row(model.ext[:PTDF], e)
                model[:ptdf_flow][e] = @constraint(
                    model,
                    dot(row, pg_bus) - model[:pf][e]
                    == ptdfb[e]
                )
                model.ext[:tracked_branches][e] = true
            end
        end
        solved = n_violated == 0
        niter += 1
    end
    solve_time = time() - t0

    if niter == model.ext[:solve_metadata][:max_ptdf_iterations]
        model.ext[:termination_info] = Dict{Symbol,Any}(
            :primal_status => MOI.INFEASIBLE_POINT,
            :dual_status => MOI.FEASIBLE_POINT,
        )
        st = MOI.ITERATION_LIMIT
    elseif st == MOI.TIME_LIMIT
        model.ext[:termination_info] = Dict{Symbol,Any}(
            :primal_status => MOI.INFEASIBLE_POINT,
            :dual_status => MOI.FEASIBLE_POINT,
        )
    elseif st ∈ (MOI.INFEASIBLE, MOI.LOCALLY_INFEASIBLE)
        model.ext[:termination_info] = Dict{Symbol,Any}(
            :primal_status => MOI.UNKNOWN_RESULT_STATUS,
            :dual_status => MOI.INFEASIBILITY_CERTIFICATE,
        )
    elseif st == MOI.DUAL_INFEASIBLE
        model.ext[:termination_info] = Dict{Symbol,Any}(
            :primal_status => MOI.INFEASIBILITY_CERTIFICATE,
            :dual_status => MOI.UNKNOWN_RESULT_STATUS,
        )
    elseif st ∈ (MOI.OPTIMAL, MOI.LOCALLY_SOLVED)
        model.ext[:termination_info] = Dict{Symbol,Any}(
            :primal_status => MOI.FEASIBLE_POINT,
            :dual_status => MOI.FEASIBLE_POINT,
        )
    else
        model.ext[:termination_info] = Dict{Symbol,Any}(
            :primal_status => MOI.UNKNOWN_RESULT_STATUS,
            :dual_status => MOI.UNKNOWN_RESULT_STATUS,
        )
    end

    model.ext[:termination_info][:solve_time] = solve_time
    model.ext[:termination_info][:ptdf_iterations] = niter
    model.ext[:termination_info][:termination_status] = st

    if has_values(model)
        # save the final pf based on pg
        pg_ = value.(model[:pg])
        p_ = Ag * pg_ - pd_nodal
        compute_flow!(pf_, p_, model.ext[:PTDF])
        model.ext[:ptdf_pf] = pf_
    end

    return
end

function extract_primal(opf::OPFModel{EconomicDispatch}) 
    model = opf.model
    T = JuMP.value_type(typeof(model))

    data = opf.data

    E, G = data.E, data.G

    primal_solution = Dict{String,Any}(
        "pg" => zeros(T, G),
        "pf" => zeros(T, E),
        "δf" => zeros(T, E),
        "r"  => zeros(T, G),
    )

    if has_values(model)
        # generator
        primal_solution["pg"] = value.(model[:pg])
        primal_solution["r"]  = value.(model[:r])

        # branch
        primal_solution["pf"] = value.(model.ext[:ptdf_pf])
        primal_solution["δf"] = value.(model[:δf])
    end

    return primal_solution
end

function extract_dual(opf::OPFModel{EconomicDispatch}) 
    model = opf.model
    T = JuMP.value_type(typeof(model))

    data = opf.data

    E, G = data.E, data.G

    dual_solution = Dict{String,Any}(
        # global
        "power_balance" => zero(T),
        "reserve_requirement" => zero(T),
        # generator
        # branch
        "ptdf_flow" => zeros(T, E),
        # Variable lower/upper bound
        "pg"     => zeros(T, G),
        "r"      => zeros(T, G),
        "pf"     => zeros(T, E),
        "δf"     => zeros(T, E),
    )


    if has_duals(model)
        # global
        dual_solution["power_balance"] = dual(model[:global_power_balance])
        dual_solution["reserve_requirement"] = dual(model[:global_reserve_requirement])

        # branch
        dual_solution["ptdf_flow"] = [
            # if a constraint wasn't added to the model, its dual is zero
            isassigned(model[:ptdf_flow], e) ? dual(model[:ptdf_flow][e]) : zero(T)
            for e in 1:E
        ]

        # Duals of variable lower/upper bounds
        # We store λ = λₗ + λᵤ, where λₗ, λᵤ are the dual variables associated to
        #   lower and upper bounds, respectively.
        # Recall that, in JuMP's convention, we have λₗ ≥ 0, λᵤ ≤ 0, hence
        #   λₗ = max(λ, 0) and λᵤ = min(λ, 0).
        # generator
        dual_solution["pg"] = dual.(LowerBoundRef.(model[:pg])) + dual.(model[:gen_max_output])
        dual_solution["r"] = dual.(LowerBoundRef.(model[:r])) + dual.(UpperBoundRef.(model[:r]))
        # branch
        dual_solution["pf"] = dual.(model[:pf_lower_bound]) + dual.(model[:pf_upper_bound])
        dual_solution["δf"] = dual.(LowerBoundRef.(model[:δf]))
    end

    return dual_solution
end
