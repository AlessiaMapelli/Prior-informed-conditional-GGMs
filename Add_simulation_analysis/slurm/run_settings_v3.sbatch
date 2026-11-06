#!/bin/bash
#SBATCH --job-name=addsim3
#SBATCH --output=logs/addsim3_%A_%a.out
#SBATCH --error=logs/addsim3_%A_%a.err
#SBATCH --partition=cpuq
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=3G
#SBATCH --time=02:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alessia.mapelli@fht.org

# One array task per row of slurm/settings_v3.csv (sample size x replication).
# Within a task the p node regressions run in parallel on the allocated
# cores. Submit from the Add_simulation_analysis folder:
#   sbatch --array=1-400%20 slurm/run_settings_v3.sbatch
# (%20 = at most 20 settings at the same time; rows are ordered from the
# largest n, so the longest tasks start first.)
# A task that hits the time limit can be resubmitted with the same array
# index: the nodes already fitted are skipped.
# Version 3: 100 replications per sample size (400 tasks), OLS refit only,
# baseline target at the mean of the continuous covariates and the reference
# level of the binary one. Resources from the v2 run (16 cores): 8-10 min per
# setting at n = 1000, 12 min at 5000, 19 min at 10000, 28-35 min at 20000;
# at most 16 GB at n = 20000. Node files and simulated data are deleted after
# the evaluation (results_v3/<setting>/rep*/ keeps evaluation, times, truth).

cd "$SLURM_SUBMIT_DIR"

ROW=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" slurm/settings_v3.csv)
N=$(echo "$ROW" | cut -d, -f1)
REP=$(echo "$ROW" | cut -d, -f2)

# One thread per worker: the parallelism is over nodes
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export LC_ALL=C.UTF-8   # avoids the LC_CTYPE warning on the nodes

source /center/healthds/singularity_functions

echo "Task $SLURM_ARRAY_TASK_ID: n=$N rep=$REP on $SLURM_CPUS_PER_TASK cores ($(hostname))"
Rscript --vanilla run_setting.R --n "$N" --rep "$REP" --p 100 --prior noisy \
  --ncores "$SLURM_CPUS_PER_TASK" --stab_B 50 --out results_v3
echo "Task $SLURM_ARRAY_TASK_ID completed"
