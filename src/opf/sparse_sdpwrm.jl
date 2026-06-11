# Part of this implementation is modified from PowerModels.jl (https://github.com/lanl-ansi/PowerModels.jl/blob/master/src/form/wrm.jl).

struct SparseSDPOPF <: AbstractFormulation end

"""
    build_opf(SparseSDPOPF, data, optimizer)

Build an SparseSDPOPF model.
"""
function build_opf(::Type{SparseSDPOPF}, data::OPFData, optimizer;
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
    model.ext[:opf_model] = SparseSDPOPF

    #
    #   I. Variables
    #
    #=
        Some generators and branches may be inactive, as indicated by `branch_status` and `gen_status`.
        Primal variables associated to inactive components are still declared, and handled as follows:
        * lower/upper bounds are set to zero
        * constraint coefficients are set to zero
    =#
    groups = data.clique_decomposition
    isnothing(groups) && @error "Clique decomposition data required for solving SparseSDPOPF. Please set compute_clique_decomposition to true when creating the OPFData or in the sampler configuration."
    # WR and WI variables of each group (clique)
    voltage_product_groups = Vector{Dict{Symbol, Matrix}}(undef, length(groups))
    w_map = Dict()
    wr_map = Dict()
    wi_map = Dict()
    visited_buses = []
    visited_directed_buspairs = []
    # Voltage magnitude and product
    for (gidx, group) in enumerate(groups)
        n = length(group)
        voltage_product_groups[gidx] = Dict()
        WR_g = model[Symbol("WR_$(gidx)")] = voltage_product_groups[gidx][:WR] = @variable(
            model,
            [1:n, 1:n],
            Symmetric;
            base_name="WR_$(gidx)"
        )
        WI_g = model[Symbol("WI_$(gidx)")] = voltage_product_groups[gidx][:WI] = @variable(
            model,
            [1:n, 1:n];
            base_name="WI_$(gidx)"
        )
        
        for i in 1:n
            i_bus = group[i]
            # Voltage magnitude bounds
            # note that bounds are added for all the linked variables and not just one of them
            set_lower_bound.(WR_g[i, i], vmin[i_bus]^2)
            set_upper_bound.(WR_g[i, i], vmax[i_bus]^2)
            # if i_bus is currently unmatched, link it to WR_g[i, i]
            if !(i_bus in visited_buses)
                push!(visited_buses, i_bus)
                w_map[i_bus] = WR_g[i, i]
            end
        end

        # Match the entries of WR_g and of WI_g to branches
        # Iterate over both directions of each bus pair in `group` since a branch may exist in either direction
        offdiag_indices = [(i, j) for i in 1:n, j in 1:n if i != j]
        for (i, j) in offdiag_indices
            i_bus, j_bus = group[i], group[j]
            # `e` is only for indexing the computed bounds values
            e = findfirst(e -> bus_fr[e] == i_bus && bus_to[e] == j_bus, 1:E)
            # if (i_bus, j_bus) is a valid directed bus pair in the system
            if !isnothing(e)
                # Only apply bounds on (i, j) entries and not (j, i).
                # Currently, these variables are not fixed to 0 when all branches between the bus pair are off.
                # The outer loop over the cliques will apply bounds to all variables linked to the directed
                # bus pair (i_bus, j_bus).
                set_upper_bound(WR_g[i, j], wr_max[e])
                set_lower_bound(WR_g[i, j], wr_min[e])
                set_upper_bound(WI_g[i, j], wi_max[e])
                set_lower_bound(WI_g[i, j], wi_min[e])
                # If the directed bus pair is currently unmatched, match it to WR_g[i, j] and WI_g[i, j]
                # for `wr` and `wi` expressions used in constraints.
                if !((i_bus, j_bus) in visited_directed_buspairs)
                    push!(visited_directed_buspairs, (i_bus, j_bus))
                    wr_map[(i_bus, j_bus)] = WR_g[i, j]
                    wi_map[(i_bus, j_bus)] = WI_g[i, j]
                end
            end
        end
        # TODO: for n == 2 we can use SOC constraint
        @constraint(model, [WR_g WI_g; -WI_g WR_g] in PSDCone(); base_name="S_$(gidx)")
    end

    # linking constraints
    # All variables representing the same variable in the dense formulation are required to be equal
    pair_matrix(group) = [(i, j) for i in group, j in group]
    overlapping_pairs = _get_overlapping_pairs(groups)
    for (i, j) in overlapping_pairs
        gi, gj = groups[i], groups[j]  # with indices in the dense formulation
        var_i, var_j = voltage_product_groups[i], voltage_product_groups[j]

        Gi, Gj = pair_matrix(gi), pair_matrix(gj)
        # Get the indices of entries in Gi and Gj that represent the same variables in the dense formulation
        overlap_i, overlap_j = _overlap_indices(Gi, Gj)
        indices = zip(overlap_i, overlap_j)
        for (idx_i, idx_j) in indices
            # Dual variables of the linking constraints cancel out when computing S
            # so we don't need to name them
            JuMP.@constraint(model, var_i[:WR][idx_i] == var_j[:WR][idx_j])
            JuMP.@constraint(model, var_i[:WI][idx_i] == var_j[:WI][idx_j])
        end
    end

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

    @expression(model, w[i in 1:N], w_map[i])

    # Nodal power balance
    @constraint(model, 
        kcl_p[i in 1:N],
        sum(gen_status[g] * pg[g] for g in bus_gens[i])
        - sum(branch_status[e] * pf[e] for e in bus_arcs_fr[i])
        - sum(branch_status[e] * pt[e] for e in bus_arcs_to[i])
        - gs[i] * w[i]
        ==
        sum(pd[l] for l in data.bus_loads[i])
    )
    @constraint(model,
        kcl_q[i in 1:N],
        sum(gen_status[g] * qg[g] for g in bus_gens[i]) 
        - sum(branch_status[e] * qf[e] for e in bus_arcs_fr[i])
        - sum(branch_status[e] * qt[e] for e in bus_arcs_to[i])
        + bs[i] * w[i]
        ==
        sum(qd[l] for l in data.bus_loads[i])
    )

    # Branch power flow physics and limit constraints
    # If e_1 and e_2 are parallel branches that connect from bus i to j, then
    # wf[e_1] and wf[e_2] will refer to the same entry in the same linked variable.
    # Similarly for wt, wr and wi.
    # ⚠️ This is different from SOCOPF, where these are defined over the branches distinctively.
    @expression(model, wf[e in 1:E], w_map[bus_fr[e]])
    @expression(model, wt[e in 1:E], w_map[bus_to[e]])
    @expression(model, wr[e in 1:E], wr_map[(bus_fr[e], bus_to[e])])
    @expression(model, wi[e in 1:E], wi_map[(bus_fr[e], bus_to[e])])

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
    
    #
    #   III. Objective
    #
    l, u = extrema(c2)
    (l == u == 0.0) || @warn "Data $(data.case) has quadratic cost terms; those terms are being ignored"
    @objective(model,
        Min,
        sum(c1[g] * pg[g] + c0[g] for g in 1:G if gen_status[g])
    )

    return OPFModel{SparseSDPOPF}(data, model)
end

function extract_primal(opf::OPFModel{SparseSDPOPF})
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
        primal_solution["w"] = value.([model[:w][i] for i in 1:N])  # diagonal of W

        # generator
        primal_solution["pg"] = value.(model[:pg])
        primal_solution["qg"] = value.(model[:qg])

        # branch
        # W is dense, so extract only the off-diagonal entries of W that correspond to branches
        # of the original network to save space. Other off-diagonal entries of W only appear
        # in the PSD constraints and not in any other constraint.
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

function extract_dual(opf::OPFModel{SparseSDPOPF})
    model = opf.model
    T = JuMP.value_type(typeof(model))

    data = opf.data

    N, E, G = data.N, data.E, data.G
    bus_fr, bus_to = data.bus_fr, data.bus_to

    dual_solution = Dict{String,Any}(
        # bus
        "kcl_p"      => zeros(T, N),
        "kcl_q"      => zeros(T, N),
        "s"          => zeros(T, N),  # diagonal of the upper-left block of S
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
        "sr"         => zeros(T, E),  # sparse entries of the upper left block of S
        "si"         => zeros(T, E),  # sparse entries of the upper right block of S
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
        groups = opf.data.clique_decomposition
        # Construct S (which has the sparsity pattern of the original network) by summing S_g of
        # every group g in the clique decomposition. S is PSD and satisfies the "sum" of all dual
        # constraints on S_g.
        # In this summed constraint, the dual variables associated with the linking constraints all
        # cancel out.
        # Compute mu_w too.
        for (gidx, group) in enumerate(groups)
            n = length(group)
            S_g = dual.(constraint_by_name(model, "S_$(gidx)"))  # 2n * 2n, with four n * n blocks
            S_g = (S_g + S_g') / 2  # ensure symmetry
            WR_g = model[Symbol("WR_$(gidx)")]
            WI_g = model[Symbol("WI_$(gidx)")]
            for i in 1:n
                i_bus = group[i]
                # Diagonal mean of S_g[1,1] and S_g[2,2] blocks
                dual_solution["s"][i_bus] += (S_g[i, i] + S_g[i + n, i + n]) / 2
                # For mu_w, we also need to sum the values multiple times for the same bus, since bounds
                # have been imposed on all linked w variables
                dual_solution["w"][i_bus] += dual(LowerBoundRef(WR_g[i, i])) + dual(UpperBoundRef(WR_g[i, i]))
            end
            # Extract only the off-diagonal entries of S that correspond to branches, since entries
            # that don't (including the ones that correspond to edges added in the chordal extension) are 0
            # (due to dual constraints).
            offdiag_indices = [(i, j) for i in 1:n, j in 1:n if i != j]
            for (i, j) in offdiag_indices
                i_bus, j_bus = group[i], group[j]
                # If there are multiple branches from bus i to j, the same (i, j) entries of S, μ_WR and μ_WI are
                # extracted for each of the branches.
                es = findall(e -> bus_fr[e] == i_bus && bus_to[e] == j_bus, 1:E)
                for e in es
                    # Mean of S_g[1,1] and S_g[2,2] blocks
                    dual_solution["sr"][e] += (S_g[i, j] + S_g[i + n, j + n]) / 2
                    # Mean of upper triangle and lower triangle of S_g[1,2] block
                    dual_solution["si"][e] += (S_g[i, j + n] - S_g[j, i + n]) / 2
                    dual_solution["wr"][e] += dual(LowerBoundRef(WR_g[i, j])) + dual(UpperBoundRef(WR_g[i, j]))
                    dual_solution["wi"][e] += (dual(LowerBoundRef(WI_g[i, j])) + dual(UpperBoundRef(WI_g[i, j])))
                end
            end
        end

        # Bus-level constraints
        dual_solution["kcl_p"] = dual.(model[:kcl_p])
        dual_solution["kcl_q"] = dual.(model[:kcl_q])

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
        # dual_solution["w"] already computed
        # generator
        dual_solution["pg"] = dual.(LowerBoundRef.(model[:pg])) + dual.(UpperBoundRef.(model[:pg]))
        dual_solution["qg"] = dual.(LowerBoundRef.(model[:qg])) + dual.(UpperBoundRef.(model[:qg]))
        # branch
        dual_solution["pf"] = dual.(LowerBoundRef.(model[:pf])) + dual.(UpperBoundRef.(model[:pf]))
        dual_solution["qf"] = dual.(LowerBoundRef.(model[:qf])) + dual.(UpperBoundRef.(model[:qf]))
        dual_solution["pt"] = dual.(LowerBoundRef.(model[:pt])) + dual.(UpperBoundRef.(model[:pt]))
        dual_solution["qt"] = dual.(LowerBoundRef.(model[:qt])) + dual.(UpperBoundRef.(model[:qt]))
    end

    return dual_solution
end

"""
    _get_overlapping_pairs(groups)

Get the indices of pairs of cliques in groups that overlap.

Refer to https://github.com/lanl-ansi/PowerModels.jl/blob/be6af59202a6868b20a41214cb341b883d62e5f0/src/form/wrm.jl#L273-L274
"""
function _get_overlapping_pairs(groups::Vector{Vector{Int}})
    # Compute the clique tree
    tree = PM._prim(PM._overlap_graph(groups))
    # Get the Cartesian indices of the adjacency matrix of the clique tree, which are exactly the indices of pairs
    # of cliques that overlap
    overlapping_pairs = [Tuple(CartesianIndices(tree)[i]) for i in (LinearIndices(tree))[findall(x->x!=0, tree)]]
    return overlapping_pairs
end

"""
    idx_a, idx_b = _overlap_indices(A, B)
Given two arrays (sizes need not match) that share some values, return:

- linear index of shared values in A
- linear index of shared values in B

Thus, A[idx_a] == B[idx_b].

Refer to https://github.com/lanl-ansi/PowerModels.jl/blob/be6af59202a6868b20a41214cb341b883d62e5f0/src/form/wrm.jl#L469
"""
function _overlap_indices(A::Array, B::Array, symmetric=true)
    return PM._overlap_indices(A, B, symmetric)
end
