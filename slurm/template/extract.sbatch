#!/bin/bash
#SBATCH --job-name=extract_OPF                # Job name
#SBATCH --account={{:charge_account}}         # charge account
#SBATCH --nodes=1                             # Use one node
#SBATCH --ntasks=1                            # Run a single task
#SBATCH --cpus-per-task={{:cpus_per_task}}    # Give CPUs to each task
#SBATCH --mem={{:extract_memory}}             # Total memory
#SBATCH --time=08:00:00                       # Time limit hrs:min:sec
#SBATCH -o {{{:logs_dir}}}/extract.out        # Combined output and error messages file
#SBATCH -q{{:queue}}

. {{{:env_path}}}

cd {{{:opfgenerator_dir}}}

julia --project=slurm -t{{:cpus_per_task}} slurm/merge.jl {{{:config_file}}}