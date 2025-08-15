Copyright Georgia Tech 2022-2024

[![Build][build-img]][build-url]
[![codecov][codecov-img]][codecov-url]
[![docs-dev][docs-dev-img]][docs-dev-url]
 
[build-img]: https://github.com/ai4opt/PGLearn.jl/actions/workflows/ci.yml/badge.svg?branch=main
[build-url]: https://github.com/ai4opt/PGLearn.jl/actions?query=workflow%3ACI
[codecov-img]: https://codecov.io/gh/AI4OPT/PGLearn.jl/graph/badge.svg
[codecov-url]: https://codecov.io/gh/AI4OPT/PGLearn.jl
[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://ai4opt.github.io/PGLearn.jl/dev/

# PGLearn
Instance generator for various OPF problems.

## Installation instructions

This repository is a non-registered Julia package.

* Option 1: install as a Julia package. You can use it, but not modify the code
    ```julia
    using Pkg
    Pkg.add(url="git@github.com:AI4OPT/PGLearn.jl.git")
    ```

* Option 2: clone the repository. Use this if you want to change the code
    ```bash
    git clone git@github.com:AI4OPT/PGLearn.jl.git
    ```
    To use the package after cloning the repo
    ```bash
    $ cd PGLearn
    $ julia --project=.
    julia> using PGLearn
    ```

    If you are modifying the source code, it is recommened to use the package [`Revise.jl`](https://github.com/timholy/Revise.jl)
    so that you can use the changes without having to start Julia.
    Make sure you load `Revise` before loading `PGLearn` in your julia session.
    ```julia
    using Revise
    using PGLearn
    ```

### Using HSL solvers (ma27, ma57)

Please refer to [`Ipopt.jl` installation instructions](https://github.com/jump-dev/Ipopt.jl?tab=readme-ov-file#linear-solvers)
    for how to install non-default linear solvers, e.g., HSL, Pardiso, etc...
Note that a specific license may be required.

To use HSL linear solvers when solving OPF instances, set the parameter "linear_solver" to "ma27" or "ma57" in the config file.
The recommended solver for Ipopt is `ma27`.
```toml
solver.name = "Ipopt"
solver.attributes.linear_solver = "ma27"
```

## Quick start

### Building and solving an OPF problem

```julia
using PGLib, PowerModels
using PGLearn
using Ipopt

# Load data from PGLib case
pm_network = PowerModels.make_basic_network(pglib("14_ieee"))
data = PGLearn.OPFData(pm_network)

# Build OPF optimization model
# To switch formulations, replace ACOPF with, e.g., DCOPF or SOCOPF
opf = PGLearn.build_opf(PGLearn.ACOPF, data, Ipopt.Optimizer)

# Solve and extract result
PGLearn.solve!(opf)
res = PGLearn.extract_result(opf)
```


### Generating random OPF instances

```julia
using Random 
using PGLib, PowerModels
using PGLearn

# Load data from PGLib case
pm_network = PowerModels.make_basic_network(pglib("14_ieee"))
data = PGLearn.OPFData(pm_network)

# Load scaler using global scaling + uncorrelated LogNormal noise
config = Dict(
    "load" => Dict(
        "noise_type" => "ScaledUniform",
        "l" => 0.8,
        "u" => 1.2,
        "sigma" => 0.10,       
    )
)
opf_sampler = SimpleOPFSampler(data, config)

# Generate a new instance
rng = MersenneTwister(42)
new_data = rand(rng, opf_sampler)

data.pd[1] # old 
0.217

new_data.pd[1]  # new
0.21480423013573954  # ⚠️ exact value depends on Julia version)
```

To generate multiple instances, run the above code in a loop
```julia
dataset = [
    PGLearn.rand(rng, opf_sampler)
    for i in 1:100
]
```

## Generating datasets

A script for generating multiple ACOPF instances is given in [`exp/sampler.jl`](exp/sampler.jl).
Make sure to install the required dependencies first:
```bash
# cd path/to/PGLearn.jl
julia --project=exp -e "using Pkg; Pkg.instantiate()"
# NOTE: if using Julia version ≤ 1.10, run instead:
# julia --project=exp -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
```

It is called from the command-line as follows:
```bash
julia --project=exp exp/sampler.jl <path/to/config.toml> <seed_min> <seed_max>
```
where
* `<path/to/config.toml>` is a path to a valid configuration file in TOML format (see [`exp/config.toml`](exp/config.toml) for an example)
* `<seed_min>` and `<seed_max>` are minimum and maximum values for the random seed. Must be integer values.
    The script will generate instances for values `smin, smin+1, ..., smax-1, smax`.

## Acknowledgements

This material is based upon work supported by the National Science Foundation under Grant No. 2112533 NSF AI Institute for Advances in Optimization ([AI4OPT](https://www.ai4opt.org/)). 
Any opinions, findings, and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the National Science Foundation.
