using PackageCompiler

# first run `julia --project=exp -t1 --trace-compile=app/precompile.jl exp/sampler.jl <path/to/config.toml> 1 1`
create_sysimage(
    [
        "CodecBzip2",
        "CodecZlib",
        "Distributions",
        "Graphs",
        "HDF5",
        "JSON",
        "JuMP",
        "LinearAlgebra",
        "PGLib",
        "PowerModels",
        "Random",
        "SparseArrays",
        "Clarabel",
        "HSL_jll",
        "HiGHS",
        "Ipopt",
        "MKL",
        "Mosek",
        "MosekTools",
        "Quadmath",
        "TOML",
        "PGLearn",
    ];
    sysimage_path="app/julia.so",
    precompile_statements_file="app/precompile.jl"
);