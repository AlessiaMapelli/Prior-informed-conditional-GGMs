#!/bin/bash
#SBATCH --job-name=sim2_analyse
#SBATCH --output=logs/analyse_%j.out
#SBATCH --error=logs/analyse_%j.err
#SBATCH --partition=cpuq
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00

# Summary tables and figures of the simulation study (analyse_results.R).
# P_VALUES restricts the analysis to some blocks, separated by ":" because
# --export already splits on commas, e.g.
#   sbatch --export=ALL,P_VALUES=50:100 slurm/analyse.sbatch
# Without it all settings found in results/ are summarised.
cd "$SLURM_SUBMIT_DIR"
export LC_ALL=C.UTF-8
source /center/healthds/singularity_functions
Rscript --vanilla analyse_results.R results ${P_VALUES}
