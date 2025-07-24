using TOML
using PowerModels
PowerModels.silence()
using PGLib

using PGLearn

include(joinpath(@__DIR__, "..", "exp", "sampler.jl"))

function main_ref(config; export_ref=false)
    case_file, case_name = PGLearn._get_case_info(config)
    export_dir = pop!(config, "export_dir")
    pop!(config, "slurm")

    data = PGLearn.OPFData(make_basic_network(PowerModels.parse_file(case_file));
        compute_clique_decomposition=get(config["sampler"], "compute_clique_decomposition", false)
    )

    d = main(data, config)
    d["data"] = PGLearn.to_dict(data)
    d["config"] = config

    if export_ref
        mkpath(export_dir)
        PGLearn.save_json(joinpath(export_dir, "case.json"), d)
    end
    return d
end

if abspath(PROGRAM_FILE) == @__FILE__
    config_file = ARGS[1]
    config = TOML.parsefile(config_file)
    main_ref(config; export_ref=true)
    exit(0)
end
