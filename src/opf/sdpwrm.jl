struct SDPOPF <: AbstractFormulation end

"""
    build_opf(SDPOPF, data, optimizer)

Build an SDPOPF model.
"""
function build_opf(::Type{SDPOPF}, data::OPFData, optimizer;
    T=Float64,
)
    # Grab some data
    N, E, G = data.N, data.E, data.G
    vmin, vmax = data.vmin, data.vmax
    i0 = data.ref_bus
    gs, bs = data.gs, data.bs
    pd, qd = data.pd, data.qd
    bus_arcs_fr, bus_arcs_to = data.bus_arcs_fr, data.bus_arcs_to
    bus_gens = data.bus_gens
    pgmin, pgmax = data.pgmin, data.pgmax
    qgmin, qgmax = data.qgmin, data.qgmax
    c0, c1, c2 = data.c0, data.c1, data.c2
    gen_status = data.gen_status
    bus_fr, bus_to = data.bus_fr, data.bus_to
    gff, gft, gtf, gtt = data.gff, data.gft, data.gtf, data.gtt
    bff, bft, btf, btt = data.bff, data.bft, data.btf, data.btt
    dvamin, dvamax, smax = data.dvamin, data.dvamax, data.smax
    branch_status = data.branch_status

    wr_min, wr_max, wi_min, wi_max = compute_voltage_phasor_bounds(data)

    model = JuMP.GenericModel{T}(optimizer)
    model.ext[:opf_model] = SDPOPF

    #
    #   I. Variables
    #
    #=
        Some generators and branches may be inactive, as indicated by `branch_status` and `gen_status`.
        Primal variables associated to inactive components are still declared, and handled as follows:
        * lower/upper bounds are set to zero
        * constraint coefficients are set to zero
    =#

    # voltage magnitude and product
    @variable(model, WR[1:N, 1:N], Symmetric)
    @variable(model, WI[1:N, 1:N])

    # Active and reactive dispatch
    @variable(model, pg[g in 1:G])
    @variable(model, qg[g in 1:G])
    
    # Directional branch flows
    @variable(model, pf[e in 1:E])
    @variable(model, qf[e in 1:E])
    @variable(model, pt[e in 1:E])
    @variable(model, qt[e in 1:E])

    # 
    #   II. Constraints
    #

    # Voltage magnitude bounds
    set_lower_bound.([WR[i, i] for i in 1:N], vmin.^2)
    set_upper_bound.([WR[i, i] for i in 1:N], vmax.^2)

    # Voltage product bounds
    # Only apply bounds on (i, j) entries and not (j, i)
    # If a bus pair has parallel branches, `set_lower_bound` and `set_upper_bound` are called
    # with the same bounds for multiple times (calls subsequent to the first are no-op).
    set_lower_bound.([WR[i,j] for (i, j) in zip(bus_fr, bus_to)], wr_min)
    set_upper_bound.([WR[i,j] for (i, j) in zip(bus_fr, bus_to)], wr_max)
    set_lower_bound.([WI[i,j] for (i, j) in zip(bus_fr, bus_to)], wi_min)
    set_upper_bound.([WI[i,j] for (i, j) in zip(bus_fr, bus_to)], wi_max)

    # Active generation bounds (both zero if generator is off)
    set_lower_bound.(pg, gen_status .* pgmin)
    set_upper_bound.(pg, gen_status .* pgmax)

    # Reactive generation bounds (both zero if generator is off)
    set_lower_bound.(qg, gen_status .* qgmin)
    set_upper_bound.(qg, gen_status .* qgmax)

    # Active flow bounds (both zero if branch is off)
    set_lower_bound.(pf, branch_status .* -smax)
    set_upper_bound.(pf, branch_status .* smax)
    set_lower_bound.(pt, branch_status .* -smax)
    set_upper_bound.(pt, branch_status .* smax)

    # Reactive flow bounds (both zero if branch is off)
    set_lower_bound.(qf, branch_status .* -smax)
    set_upper_bound.(qf, branch_status .* smax)
    set_lower_bound.(qt, branch_status .* -smax)
    set_upper_bound.(qt, branch_status .* smax)

    # Nodal power balance
    @constraint(model, 
        kcl_p[i in 1:N],
        sum(gen_status[g] * pg[g] for g in bus_gens[i])
        - sum(branch_status[e] * pf[e] for e in bus_arcs_fr[i])
        - sum(branch_status[e] * pt[e] for e in bus_arcs_to[i])
        - gs[i] * WR[i, i]
        ==
        sum(pd[l] for l in data.bus_loads[i])
    )
    @constraint(model,
        kcl_q[i in 1:N],
        sum(gen_status[g] * qg[g] for g in bus_gens[i]) 
        - sum(branch_status[e] * qf[e] for e in bus_arcs_fr[i])
        - sum(branch_status[e] * qt[e] for e in bus_arcs_to[i])
        + bs[i] * WR[i, i]
        ==
        sum(qd[l] for l in data.bus_loads[i])
    )

    # Branch power flow physics and limit constraints
    # If e_1 and e_2 are parallel branches that connect from bus i to j, then
    # wf[e_1] and wf[e_2] will refer to the same entry in W.
    # Similarly for wt, wr and wi.
    # ⚠️ This is different from SOCOPF, where these are defined over the branches distinctively.
    @expression(model, wf[e in 1:E], WR[bus_fr[e], bus_fr[e]])
    @expression(model, wt[e in 1:E], WR[bus_to[e], bus_to[e]])
    @expression(model, wr[e in 1:E], WR[bus_fr[e], bus_to[e]])
    @expression(model, wi[e in 1:E], WI[bus_fr[e], bus_to[e]])

    # Ohm's law
    @constraint(model, ohm_pf[e in 1:E],
        branch_status[e] * ( gff[e] * wf[e] + gft[e] * wr[e] + bft[e] * wi[e]) - pf[e] == 0
    )
    @constraint(model, ohm_qf[e in 1:E],
        branch_status[e] * (-bff[e] * wf[e] - bft[e] * wr[e] + gft[e] * wi[e]) - qf[e] == 0
    )
    @constraint(model, ohm_pt[e in 1:E],
        branch_status[e] * ( gtt[e] * wt[e] + gtf[e] * wr[e] - btf[e] * wi[e]) - pt[e] == 0
    )
    @constraint(model, ohm_qt[e in 1:E],
        branch_status[e] * (-btt[e] * wt[e] - btf[e] * wr[e] - gtf[e] * wi[e]) - qt[e] == 0
    )

    # Thermal limit
    @constraint(model, sm_fr[e in 1:E], [smax[e], pf[e], qf[e]] in SecondOrderCone())
    @constraint(model, sm_to[e in 1:E], [smax[e], pt[e], qt[e]] in SecondOrderCone())

    # Voltage angle difference limit
    @constraint(model, va_diff_lb[e in 1:E], branch_status[e] * wi[e] - branch_status[e] * tan(dvamin[e]) * wr[e] >= 0)
    @constraint(model, va_diff_ub[e in 1:E], branch_status[e] * wi[e] - branch_status[e] * tan(dvamax[e]) * wr[e] <= 0)
    
    # PSD constraint
    @constraint(model, S, [WR WI; -WI WR] in PSDCone())
    
    #
    #   III. Objective
    #
    l, u = extrema(c2)
    (l == u == 0.0) || @warn "Data $(data.case) has quadratic cost terms; those terms are being ignored"
    @objective(model,
        Min,
        sum(c1[g] * pg[g] + c0[g] for g in 1:G if gen_status[g])
    )

    return OPFModel{SDPOPF}(data, model)
end

function extract_primal(opf::OPFModel{SDPOPF})
    model = opf.model
    T = JuMP.value_type(typeof(model))

    data = opf.data

    N, E, G = data.N, data.E, data.G

    primal_solution = Dict{String,Any}(
        # bus
        "w" => zeros(T, N),
        # generator
        "pg" => zeros(T, G),
        "qg" => zeros(T, G),
        # branch
        "wr" => zeros(T, E),
        "wi" => zeros(T, E),
        "pf" => zeros(T, E),
        "qf" => zeros(T, E),
        "pt" => zeros(T, E),
        "qt" => zeros(T, E),
    )
    if has_values(model)
        # bus
        primal_solution["w"] = value.([model[:WR][i, i] for i in 1:N])  # diagonal of W

        # generator
        primal_solution["pg"] = value.(model[:pg])
        primal_solution["qg"] = value.(model[:qg])

        # branch
        # W is dense, so extract only the off-diagonal entries of W that correspond to branches
        # to save space. Other off-diagonal entries of W only appear in the PSD constraint and
        # not in any other constraint.
        # These entries can be recovered by solving a PSD matrix completion problem.
        # If there are multiple branches from bus i to j, the same (i, j) entry of W is extracted
        # for each of the branches.
        primal_solution["wr"] = value.(model[:wr])
        primal_solution["wi"] = value.(model[:wi])
        primal_solution["pf"] = value.(model[:pf])
        primal_solution["qf"] = value.(model[:qf])
        primal_solution["pt"] = value.(model[:pt])
        primal_solution["qt"] = value.(model[:qt])
    end

    return primal_solution
end

function extract_dual(opf::OPFModel{SDPOPF})
    model = opf.model
    T = JuMP.value_type(typeof(model))

    data = opf.data

    N, E, G = data.N, data.E, data.G
    bus_fr, bus_to = data.bus_fr, data.bus_to

    dual_solution = Dict{String,Any}(
        # bus
        "kcl_p"      => zeros(T, N),
        "kcl_q"      => zeros(T, N),
        "s"          => zeros(T, N),
        # generator
        # N/A
        # branch
        "ohm_pf"     => zeros(T, E),
        "ohm_pt"     => zeros(T, E),
        "ohm_qf"     => zeros(T, E),
        "ohm_qt"     => zeros(T, E),
        "va_diff"    => zeros(T, E),
        "sm_fr"      => zeros(T, E, 3),
        "sm_to"      => zeros(T, E, 3),
        "sr"         => zeros(T, E),
        "si"         => zeros(T, E),
        # variables lower/upper bounds
        # bus
        "w"          => zeros(T, N),
        # generator
        "pg"         => zeros(T, G),
        "qg"         => zeros(T, G),
        # branch
        "wr"         => zeros(T, E),
        "wi"         => zeros(T, E),
        "pf"         => zeros(T, E),
        "qf"         => zeros(T, E),
        "pt"         => zeros(T, E),
        "qt"         => zeros(T, E),
    )

    if has_duals(model)
        S = dual.(model[:S])  # 2N * 2N, with four N * N blocks
        S = (S + S') / 2  # ensure symmetry

        # Bus-level constraints
        dual_solution["kcl_p"] = dual.(model[:kcl_p])
        dual_solution["kcl_q"] = dual.(model[:kcl_q])
        # Diagonal mean of S[1,1] and S[2,2] blocks
        dual_solution["s"] = [(S[i, i] + S[i + N, i + N]) / 2 for i in 1:N]

        # Generator-level constraints
        # N/A

        # Branch-level constraints
        dual_solution["ohm_pf"] = dual.(model[:ohm_pf])
        dual_solution["ohm_pt"] = dual.(model[:ohm_pt])
        dual_solution["ohm_qf"] = dual.(model[:ohm_qf])
        dual_solution["ohm_qt"] = dual.(model[:ohm_qt])
        dual_solution["va_diff"] = dual.(model[:va_diff_lb]) + dual.(model[:va_diff_ub])  # same as bound constraints
        dual_solution["sm_fr"] = dual.(model[:sm_fr])
        dual_solution["sm_to"] = dual.(model[:sm_to])
        # Extract only the off-diagonal entries of S that correspond to branches, since entries
        # that don't are 0.
        # If there are multiple branches from bus i to j, the same (i, j) entry of S is extracted for
        # each of the branches.
        # Mean of S[1,1] and S[2,2] blocks
        dual_solution["sr"] = [(S[bus_fr[e], bus_to[e]] + S[bus_fr[e] + N, bus_to[e] + N]) / 2 for e in 1:E]
        # Mean of upper triangle and lower triangle of S[1,2] block
        dual_solution["si"] = [(S[bus_fr[e], bus_to[e] + N] - S[bus_to[e], bus_fr[e] + N]) / 2 for e in 1:E]

        # For conic constraints, JuMP will return Vector{Vector{T}}
        # reshape duals of conic constraints into matrix shape
        dual_solution["sm_fr"] = mapreduce(permutedims, vcat, dual_solution["sm_fr"])
        dual_solution["sm_to"] = mapreduce(permutedims, vcat, dual_solution["sm_to"])

        # Duals of variable lower/upper bounds
        # We store λ = λₗ + λᵤ, where λₗ, λᵤ are the dual variables associated to
        #   lower and upper bounds, respectively.
        # Recall that, in JuMP's convention, we have λₗ ≥ 0, λᵤ ≤ 0, hence
        #   λₗ = max(λ, 0) and λᵤ = min(λ, 0).

        # bus
        dual_solution["w"]  = dual.(LowerBoundRef.([model[:WR][i, i] for i in 1:N])) + dual.(UpperBoundRef.([model[:WR][i, i] for i in 1:N]))
        # generator
        dual_solution["pg"] = dual.(LowerBoundRef.(model[:pg])) + dual.(UpperBoundRef.(model[:pg]))
        dual_solution["qg"] = dual.(LowerBoundRef.(model[:qg])) + dual.(UpperBoundRef.(model[:qg]))
        # branch
        # If there are multiple branches from bus i to j, the same (i, j) entries of μ_WR and μ_WI are extracted for
        # each of the branches.
        dual_solution["wr"] = dual.(LowerBoundRef.([model[:WR][i, j] for (i, j) in zip(bus_fr, bus_to)])) + dual.(UpperBoundRef.([model[:WR][i, j] for (i, j) in zip(bus_fr, bus_to)]))
        dual_solution["wi"] = (dual.(LowerBoundRef.([model[:WI][i, j] for (i, j) in zip(bus_fr, bus_to)])) + dual.(UpperBoundRef.([model[:WI][i, j] for (i, j) in zip(bus_fr, bus_to)]))) .* 2
        dual_solution["pf"] = dual.(LowerBoundRef.(model[:pf])) + dual.(UpperBoundRef.(model[:pf]))
        dual_solution["qf"] = dual.(LowerBoundRef.(model[:qf])) + dual.(UpperBoundRef.(model[:qf]))
        dual_solution["pt"] = dual.(LowerBoundRef.(model[:pt])) + dual.(UpperBoundRef.(model[:pt]))
        dual_solution["qt"] = dual.(LowerBoundRef.(model[:qt])) + dual.(UpperBoundRef.(model[:qt]))
    end

    return dual_solution
end