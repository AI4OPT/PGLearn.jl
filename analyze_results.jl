using HDF5

# Load the HDF5 file
h5file = "/home/runner/work/PGLearn.jl/PGLearn.jl/data/14_ieee/res_h5/pglib_opf_case14_ieee_SOCOPF128_s1-s64.h5"

h5open(h5file, "r") do file
    # Read termination status
    termination_status = read(file, "meta/termination_status")
    primal_status = read(file, "meta/primal_status")
    dual_status = read(file, "meta/dual_status")
    primal_objective = read(file, "meta/primal_objective_value")
    dual_objective = read(file, "meta/dual_objective_value")
    
    println("="^80)
    println("Analysis of solver statuses for 64 samples")
    println("="^80)
    
    # Count termination statuses
    println("\n1. Termination Status Distribution:")
    status_counts = Dict{String, Int}()
    for status in termination_status
        status_counts[status] = get(status_counts, status, 0) + 1
    end
    for (status, count) in sort(collect(status_counts), by=x->x[2], rev=true)
        println("   - $status: $count/64 ($(round(100*count/64, digits=1))%)")
    end
    
    # Count primal statuses
    println("\n2. Primal Status Distribution:")
    primal_counts = Dict{String, Int}()
    for status in primal_status
        primal_counts[status] = get(primal_counts, status, 0) + 1
    end
    for (status, count) in sort(collect(primal_counts), by=x->x[2], rev=true)
        println("   - $status: $count/64 ($(round(100*count/64, digits=1))%)")
    end
    
    # Count dual statuses
    println("\n3. Dual Status Distribution:")
    dual_counts = Dict{String, Int}()
    for status in dual_status
        dual_counts[status] = get(dual_counts, status, 0) + 1
    end
    for (status, count) in sort(collect(dual_counts), by=x->x[2], rev=true)
        println("   - $status: $count/64 ($(round(100*count/64, digits=1))%)")
    end
    
    # Check duality gap for OPTIMAL solutions
    println("\n4. Duality Gap Analysis:")
    optimal_indices = findall(x -> x == "OPTIMAL", termination_status)
    almost_optimal_indices = findall(x -> x == "ALMOST_OPTIMAL", termination_status)
    
    println("   OPTIMAL solutions: $(length(optimal_indices))")
    if length(optimal_indices) > 0
        gaps = abs.(primal_objective[optimal_indices] .- dual_objective[optimal_indices])
        println("   - Max gap: $(maximum(gaps))")
        println("   - Mean gap: $(sum(gaps)/length(gaps))")
        println("   - All gaps <= 1e-8: $(all(gaps .<= 1e-8))")
    end
    
    println("\n   ALMOST_OPTIMAL solutions: $(length(almost_optimal_indices))")
    if length(almost_optimal_indices) > 0
        gaps = abs.(primal_objective[almost_optimal_indices] .- dual_objective[almost_optimal_indices])
        println("   - Max gap: $(maximum(gaps))")
        println("   - Mean gap: $(sum(gaps)/length(gaps))")
        println("   - All gaps <= 1e-8: $(all(gaps .<= 1e-8))")
    end
    
    println("\n" * "="^80)
    println("RESULTS SUMMARY:")
    println("="^80)
    
    # Check if conditions are met
    optimal_count = get(status_counts, "OPTIMAL", 0)
    feasible_primal = get(primal_counts, "FEASIBLE_POINT", 0)
    feasible_dual = get(dual_counts, "FEASIBLE_POINT", 0)
    
    println("✓ Mostly OPTIMAL termination status: $(optimal_count >= 32) ($(optimal_count)/64)")
    println("✓ Mostly FEASIBLE_POINT primal status: $(feasible_primal >= 32) ($(feasible_primal)/64)")
    println("✓ Mostly FEASIBLE_POINT dual status: $(feasible_dual >= 32) ($(feasible_dual)/64)")
    
    if length(optimal_indices) > 0
        gaps = abs.(primal_objective[optimal_indices] .- dual_objective[optimal_indices])
        gaps_ok = all(gaps .<= 1e-8)
        println("✓ Duality gap <= 1e-8 for OPTIMAL: $gaps_ok")
    end
    
    println("="^80)
end
