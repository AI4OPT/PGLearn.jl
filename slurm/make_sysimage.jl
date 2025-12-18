using PackageCompiler

# first run `julia --project=. -t1 --trace-compile=app/precompile.jl exp/sampler.jl <path/to/config.toml> 1 1`
create_sysimage(
    [
        "Clarabel",
        "HSL_jll",
        "HiGHS",
        "Ipopt",
        "JuMP",
        "LinearAlgebra",
        "MKL",
        "Mosek",
        "MosekTools",
        "PGLearn",
        "PGLib",
        "Pkg",
        "PowerModels",
        "Quadmath",
        "Random",
        "TOML",
    ];
    sysimage_path="app/julia.so",
    precompile_statements_file="app/precompile.jl"
);