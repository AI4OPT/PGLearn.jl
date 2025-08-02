#!/bin/bash
#SBATCH --job-name=ref_OPF                    # Job name
#SBATCH --account={{:charge_account}}         # charge account
#SBATCH --nodes=1                             # Use one node
#SBATCH --ntasks=1                            # Run a single task
#SBATCH --cpus-per-task=1                     # Give 1 CPU
#SBATCH --mem-per-cpu={{:ref_memory}}         # Memory per processor
#SBATCH --time=01:00:00                       # Time limit hrs:min:sec
#SBATCH -o {{{:logs_dir}}}/ref.out            # Combined output and error messages file
#SBATCH -q{{:queue}}

. {{{:env_path}}}

cd {{{:opfgenerator_dir}}}

julia --project=exp -t1 {{{:opfgenerator_dir}}}/slurm/make_ref.jl {{{:config_file}}}