struct ACOPFParam <: AbstractFormulation end

"""
    build_opf(ACOPFParam, data, optimizer)

Build a parametric ACOPF model where pd (active load) and qd (reactive load)
are MOI.Parameter variables.
Returns the model and the parameter variables (pd, qd) for differentiation.
"""
function build_opf(::Type{ACOPFParam}, data::OPFData, optimizer;
    T=Float64,
)
    model = JuMP.GenericModel{T}(optimizer)
    return build_opf!(ACOPFParam, model, data)
end

function build_opf!(::Type{ACOPFParam}, model::JuMP.AbstractModel, data::OPFData)
    N, E, G, L = data.N, data.E, data.G, data.L
    vmin, vmax = data.vmin, data.vmax
    i0 = data.ref_bus
    gs, bs = data.gs, data.bs
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

    model.ext[:opf_model] = ACOPFParam

    #
    #   0. Parameters (pd and qd as MOI.Parameter)
    #
    @variable(model, pd[l in 1:L] in MOI.Parameter(data.pd[l]))
    @variable(model, qd[l in 1:L] in MOI.Parameter(data.qd[l]))

    #
    #   I. Variables
    #

    # nodal voltage
    @variable(model, vm[1:N], start=1.0)
    @variable(model, va[1:N])

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
    set_lower_bound.(vm, vmin)
    set_upper_bound.(vm, vmax)

    # Active generation bounds
    set_lower_bound.(pg, gen_status .* pgmin)
    set_upper_bound.(pg, gen_status .* pgmax)

    # Reactive generation bounds
    set_lower_bound.(qg, gen_status .* qgmin)
    set_upper_bound.(qg, gen_status .* qgmax)

    # Active flow bounds
    set_lower_bound.(pf, branch_status .* -smax)
    set_upper_bound.(pf, branch_status .* smax)
    set_lower_bound.(pt, branch_status .* -smax)
    set_upper_bound.(pt, branch_status .* smax)

    # Reactive flow bounds
    set_lower_bound.(qf, branch_status .* -smax)
    set_upper_bound.(qf, branch_status .* smax)
    set_lower_bound.(qt, branch_status .* -smax)
    set_upper_bound.(qt, branch_status .* smax)

    # Slack bus
    @constraint(model, slack_bus, va[i0] == 0.0)

    # Nodal power balance with parametric pd
    @constraint(model,
        kcl_p[i in 1:N],
        sum(gen_status[g] * pg[g] for g in bus_gens[i])
        - sum(branch_status[e] * pf[e] for e in bus_arcs_fr[i])
        - sum(branch_status[e] * pt[e] for e in bus_arcs_to[i])
        - gs[i] * vm[i]^2
        ==
        sum(pd[l] for l in data.bus_loads[i])
    )
    @constraint(model,
        kcl_q[i in 1:N],
        sum(gen_status[g] * qg[g] for g in bus_gens[i])
        - sum(branch_status[e] * qf[e] for e in bus_arcs_fr[i])
        - sum(branch_status[e] * qt[e] for e in bus_arcs_to[i])
        + bs[i] * vm[i]^2
        ==
        sum(qd[l] for l in data.bus_loads[i])
    )

    # Ohm's law
    @expression(model, wf[e in 1:E], vm[bus_fr[e]]^2)
    @expression(model, wt[e in 1:E], vm[bus_to[e]]^2)
    @expression(model, wr[e in 1:E], vm[bus_fr[e]] * vm[bus_to[e]] * cos(va[bus_fr[e]] - va[bus_to[e]]))
    @expression(model, wi[e in 1:E], vm[bus_fr[e]] * vm[bus_to[e]] * sin(va[bus_fr[e]] - va[bus_to[e]]))

    @constraint(model,
        ohm_pf[e in 1:E],
        branch_status[e] * ( gff[e] * wf[e] + gft[e] * wr[e] + bft[e] * wi[e]) - pf[e] == 0
    )
    @constraint(model,
        ohm_qf[e in 1:E],
        branch_status[e] * (-bff[e] * wf[e] - bft[e] * wr[e] + gft[e] * wi[e]) - qf[e] == 0
    )
    @constraint(model,
        ohm_pt[e in 1:E],
        branch_status[e] * ( gtt[e] * wt[e] + gtf[e] * wr[e] - btf[e] * wi[e]) - pt[e] == 0
    )
    @constraint(model,
        ohm_qt[e in 1:E],
        branch_status[e] * (-btt[e] * wt[e] - btf[e] * wr[e] - gtf[e] * wi[e]) - qt[e] == 0
    )

    # Thermal limit
    @constraint(model, sm_fr[e in 1:E], pf[e]^2 + qf[e]^2 <= smax[e]^2)
    @constraint(model, sm_to[e in 1:E], pt[e]^2 + qt[e]^2 <= smax[e]^2)

    # Voltage angle difference limit
    @constraint(model,
        va_diff[e in 1:E],
        dvamin[e] <= branch_status[e] * (va[bus_fr[e]] - va[bus_to[e]]) <= dvamax[e]
    )

    #
    #   III. Objective
    #
    l, u = extrema(c2)
    (l == u == 0.0) || @warn "Data $(data.case) has quadratic cost terms; those terms are being ignored"
    @objective(model,
        Min,
        sum(c1[g] * pg[g] + c0[g] for g in 1:G if gen_status[g])
    )

    return OPFModel{ACOPFParam}(data, model), pd, qd
end
