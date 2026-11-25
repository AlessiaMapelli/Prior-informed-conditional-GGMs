#!/bin/bash
#SBATCH --job-name=ablation
#SBATCH --output=logs/ablation_%A_%a.out
#SBATCH --error=logs/ablation_%A_%a.err
#SBATCH --partition=cpuq
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=3G
#SBATCH --time=01:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alessia.mapelli@fht.org

# Ablation of the selective prior (README.md, "Ablation: prior on all terms"):
# the main-study datasets of p = 100, q = 3 (replications 1-100) refitted with
# the prior also weighting the covariate-specific terms, under the noisy and
# perfect priors (bic_relaxed). The prCGGM arms (no prior; prior on the
# baseline only) are the main-study fits of the same datasets in results/.
#   tasks   1-100  n = 20000
#   tasks 101-200  n = 10000
#   tasks 201-300  n = 5000
#   tasks 301-400  n = 1000
# Submit from the Simulation_analysis_v2 folder:
#   sbatch --array=1-400%50 slurm/run_ablation.sbatch
# Finished tasks are skipped when resubmitted.

cd "$SLURM_SUBMIT_DIR"
ROW=$(sed -n "${SLURM_ARRAY_TASK_ID}p" slurm/settings_ablation.csv)
IFS=, read -r N P Q REP <<< "$ROW"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export LC_ALL=C.UTF-8
source /center/healthds/singularity_functions

echo "Ablation task $SLURM_ARRAY_TASK_ID: n=$N p=$P q=$Q rep=$REP on $SLURM_CPUS_PER_TASK cores ($(hostname))"
Rscript --vanilla run_task.R --n "$N" --p "$P" --q "$Q" --rep "$REP" \
  --ncores "$SLURM_CPUS_PER_TASK" --priors noisy,perfect --selection bic_relaxed \
  --prior_on_covariates TRUE --out results_ablation
echo "Ablation task $SLURM_ARRAY_TASK_ID completed"
