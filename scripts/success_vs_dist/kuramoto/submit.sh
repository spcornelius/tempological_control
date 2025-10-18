#!/bin/bash
#SBATCH --job-name=kuramoto_dist
#SBATCH --ntasks=200
#SBATCH --mem-per-cpu=2GB
#SBATCH --time=05:00:00
#SBATCH --output=kuramoto_dist.log
#SBATCH --error=kuramoto_dist.err

script_path=$(scontrol show job $SLURM_JOBID | awk -F= '/Command=/{print $2}')
cd $(dirname $script_path)

julia --project --sysimage=$HOME/research/temporal_sync/$JULIA_SYSIMAGE_NAME kuramoto.jl