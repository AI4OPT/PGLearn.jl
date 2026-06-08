using PGLearn
using HDF5
using JSON
using TOML
using Base.Threads
using Random

main(fconfig::AbstractString) = main(TOML.parsefile(fconfig))

function main(config::Dict)
    export_dir = pop!(config, "export_dir")
    case_file, case_name = PGLearn._get_case_info(config)
    all_h5_files = filter(endswith(".h5"), readdir(joinpath(export_dir, "res_h5"), join=true))

    slurm_config = pop!(config, "slurm")

    # Process each dataset
    OPFs = sort(collect(keys(config["OPF"])))
    for dataset_name in ["input"; OPFs]
        rgx = Regex("$(case_name)_$(dataset_name)_s\\d+-s\\d+.h5")
        res_files = filter(s -> startswith(basename(s), rgx), all_h5_files)
        nfiles = length(res_files)
        Ds = Vector{Dict{String,Any}}(undef, nfiles)
        @threads for i in 1:nfiles
            Ds[i] = h5read(res_files[i], "/")
        end
        
        # Merge minibatches, sort, and export to disk
        D = PGLearn._merge_h5(Ds)
        # _dedupe_and_sort_h5! expects the input dictionary to have a D["meta"]["seed"],
        #     so we artificially create this for input data...
        if dataset_name == "input"
            D["meta"] = Dict("seed" => copy(D["seed"]))
        end
        PGLearn._dedupe_and_sort_h5!(D)
        # ... and delete it here
        if dataset_name == "input"
            delete!(D, "meta")
        end

        # Save dataset to disk
        if dataset_name == "input"
            PGLearn.save_h5(joinpath(export_dir, "input.h5"), D)
        else
            mkpath(joinpath(export_dir, dataset_name))
            PGLearn.save_h5(joinpath(export_dir, dataset_name, "meta.h5"), D["meta"])
            PGLearn.save_h5(joinpath(export_dir, dataset_name, "primal.h5"), D["primal"])
            PGLearn.save_h5(joinpath(export_dir, dataset_name, "dual.h5"), D["dual"])
        end

        GC.gc()
    end

    # Cleanup temp h5 files
    rm(joinpath(export_dir, "res_h5"), recursive=true)

    if "--no-split" in ARGS || get(slurm_config, "no_split", false)
        return
    end

    rng = MersenneTwister(42)

    feasible = nothing
    n_samples = nothing

    # find inputs that are feasible for all OPFs
    for opf in OPFs
        h5open(joinpath(export_dir, opf, "meta.h5"), "r") do f
            if isnothing(feasible)
                n_samples = length(read(f["termination_status"]))
                feasible = trues(n_samples)
            end
            feasible .&= (read(f["termination_status"]) .== "OPTIMAL") .| (read(f["termination_status"]) .== "LOCALLY_SOLVED") .| (read(f["termination_status"]) .== "SLOW_PROGRESS")
            feasible .&= (read(f["primal_status"]) .== "FEASIBLE_POINT")
            feasible .&= (read(f["dual_status"]) .== "FEASIBLE_POINT")
        end
    end

    n_feasible = sum(feasible)
    @info "Number of feasible samples: $n_feasible / $n_samples"

    feasible_idx = findall(feasible)
    infeasible_idx = sort(findall(.!feasible))

    shuffle!(rng, feasible_idx)

    n_train = Int(floor(get(slurm_config, "train_ratio", 0.8) * n_feasible))
    train_idx = sort(feasible_idx[1:n_train])
    test_idx = sort(feasible_idx[n_train+1:end])

    @info "Number of train samples: $(length(train_idx))"
    @info "Number of test samples: $(length(test_idx))"

    mkpath(joinpath(export_dir, "train"))
    mkpath(joinpath(export_dir, "test"))
    mkpath(joinpath(export_dir, "infeasible"))

    # g[k] = f[k][idx]  ∀ k ∈ keys(f)
    function _copy_to_new_h5(g, f, idx)
        for k in keys(f)
            val = read(f[k])
            g[k] = collect(selectdim(val, ndims(val), idx))
        end
    end

    # Split consolidated H5s into train/test/infeasible sets
    all_h5_paths = ["input"; [joinpath(opf, file) for opf in OPFs for file in ["primal", "dual", "meta"]]]
    for (idx, split) in zip([train_idx, test_idx, infeasible_idx], ["train", "test", "infeasible"])
        for p in all_h5_paths
            src_path = joinpath(export_dir, p * ".h5")
            h5open(src_path, "r") do f
                dst_path = joinpath(export_dir, split, p * ".h5")
                mkpath(dirname(dst_path))
                h5open(dst_path, "w") do g
                    _copy_to_new_h5(g, f, idx)
                end  # close dst
            end  # close src
        end
    end
    # Cleanup workspace
    rm(joinpath(export_dir, "input.h5"))
    for opf in OPFs
        rm(joinpath(export_dir, opf), recursive=true)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    fconfig = ARGS[1]
    # TODO: check ARGS and that fconfig is valid path
    main(fconfig)
    exit(0)
end