#!/bin/bash
#SBATCH --job-name=sim_regen
#SBATCH --output=logs/sim_regen_%j.out
#SBATCH --error=logs/sim_regen_%j.err
#SBATCH --partition=cpuq
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=03:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=alessia.mapelli@fht.org

# Regenerates the simulation results without refitting (see
# Sim_regenerate_results.R). Single-threaded: it reads 480 collections of
# node files and writes the figures. Submit from Simulation_analysis/:
#   mkdir -p logs && sbatch slurm_regenerate_results.sbatch

cd "$SLURM_SUBMIT_DIR"
export LC_ALL=C.UTF-8
source /center/healthds/singularity_functions
Rscript --vanilla Sim_regenerate_results.R
