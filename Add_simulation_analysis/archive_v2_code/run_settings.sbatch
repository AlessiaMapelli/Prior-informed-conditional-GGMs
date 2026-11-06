#!/bin/bash
#SBATCH --job-name=addsim
#SBATCH --output=logs/addsim_%A_%a.out
#SBATCH --error=logs/addsim_%A_%a.err
#SBATCH --partition=cpuq
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=3G
#SBATCH --time=04:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alessia.mapelli@fht.org

# One array task per row of slurm/settings.csv (sample size x replication).
# Within a task the p node regressions run in parallel on the allocated
# cores. Submit from the Add_simulation_analysis folder:
#   sbatch --array=1-20%10 slurm/run_settings.sbatch
# (%10 = at most 10 settings at the same time; rows are ordered from the
# largest n, so the longest tasks start first.)
# A task that hits the time limit can be resubmitted with the same array
# index: the nodes already fitted are skipped.
# Resources from the timing test (cnode01, 8 cores): about 215 s per node
# per core at n = 1000-5000 and 425 s at n = 20000; about 2 GB per worker
# at n = 20000. With 16 cores a setting takes 30-65 minutes.

cd "$SLURM_SUBMIT_DIR"

ROW=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" slurm/settings.csv)
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
  --ncores "$SLURM_CPUS_PER_TASK" --stab_B 50 --out results_v2
echo "Task $SLURM_ARRAY_TASK_ID completed"
