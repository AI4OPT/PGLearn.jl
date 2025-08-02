#!/bin/bash
#SBATCH --job-name=sysimage_OPF               # Job name
#SBATCH --account={{:charge_account}}         # charge account
#SBATCH --nodes=1                             # Use one node
#SBATCH --ntasks=1                            # Run a single task
#SBATCH --cpus-per-task=1                     # Give 1 CPU
#SBATCH --mem-per-cpu={{:sysimage_memory}}    # Memory per processor
#SBATCH --time=01:00:00                       # Time limit hrs:min:sec
#SBATCH -o {{{:logs_dir}}}/sysimage.out       # Combined output and error messages file
#SBATCH -q{{:queue}}

. {{{:env_path}}}

cd {{{:opfgenerator_dir}}}

mkdir app

julia --project=exp -t1 --trace-compile=app/precompile.jl {{{:sampler_script}}} {{{:config_file}}} 1 1

julia --project=slurm slurm/make_sysimage.jl