#!/bin/bash
#SBATCH --job-name=stabtime
#SBATCH --partition=cpuq
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G
#SBATCH --time=06:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alessia.mapelli@fht.org

# Timing test of stability_bic in the P/n >= 0.5 cells: n = 1000, p = 500,
# rep 1, Noisy prior, nodes 1-16 (one node per core). Writes nodes_timing.rds.
# A bash script is needed: `sbatch --wrap` runs under /bin/sh, where sourcing
# /center/healthds/singularity_functions fails at line 43 (function `aws-cli`).
# Submit from the Simulation_analysis_v2 folder, one job per q:
#   for Q in 1 3; do sbatch --output=logs/stabtime_q${Q}_%j.out slurm/stability_timing.sbatch $Q; done

Q=${1:?usage: sbatch slurm/stability_timing.sbatch <q>}
cd "$SLURM_SUBMIT_DIR"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export LC_ALL=C.UTF-8
source /center/healthds/singularity_functions

echo "Stability timing: n=1000 p=500 q=$Q on $SLURM_CPUS_PER_TASK cores ($(hostname))"
time Rscript --vanilla run_task.R --n 1000 --p 500 --q "$Q" --rep 1 \
  --ncores "$SLURM_CPUS_PER_TASK" --nodes 1-16 --priors noisy \
  --selection stability_bic --out results_stability_timing
