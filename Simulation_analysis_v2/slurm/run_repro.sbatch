#!/bin/bash
#SBATCH --job-name=repro
#SBATCH --output=logs/repro_%A_%a.out
#SBATCH --error=logs/repro_%A_%a.err
#SBATCH --partition=cpuq
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=3G
#SBATCH --time=02:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alessia.mapelli@fht.org

# Reproducibility experiment (README.md, "Reproducibility experiment"): one
# array task per row of slurm/settings_repro.csv (n, p, q, dataset), each
# dataset drawn from the fixed true network of (p, q) and fitted with
# bic_relaxed, ebic_relaxed and stability_bic under the noisy and no prior.
#   tasks   1- 50  n = 10000
#   tasks  51-100  n = 5000
#   tasks 101-150  n = 1000
# Submit from the Simulation_analysis_v2 folder:
#   sbatch --array=1-150%30 slurm/run_repro.sbatch
# Finished tasks (networks.rds present) are skipped when resubmitted.

cd "$SLURM_SUBMIT_DIR"
ROW=$(sed -n "${SLURM_ARRAY_TASK_ID}p" slurm/settings_repro.csv)
IFS=, read -r N P Q REP <<< "$ROW"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export LC_ALL=C.UTF-8
source /center/healthds/singularity_functions

echo "Repro task $SLURM_ARRAY_TASK_ID: n=$N p=$P q=$Q dataset=$REP on $SLURM_CPUS_PER_TASK cores ($(hostname))"
Rscript --vanilla run_repro_task.R --n "$N" --p "$P" --q "$Q" --rep "$REP" \
  --ncores "$SLURM_CPUS_PER_TASK" --out results_reproducibility
echo "Repro task $SLURM_ARRAY_TASK_ID completed"
