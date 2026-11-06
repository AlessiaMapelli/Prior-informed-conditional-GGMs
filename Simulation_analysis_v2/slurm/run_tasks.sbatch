#!/bin/bash
#SBATCH --job-name=sim2
#SBATCH --output=logs/sim2_%A_%a.out
#SBATCH --error=logs/sim2_%A_%a.err
#SBATCH --partition=cpuq
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G
#SBATCH --time=04:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alessia.mapelli@fht.org

# Simulation study: one array task per row of slurm/settings.csv (one dataset
# (n, p, q, replication) fitted under the three priors). The rows are ordered
# by p, so that each block can be submitted with its own resources:
#   tasks    1- 800  p = 500
#   tasks  801-1600  p = 250
#   tasks 1601-2400  p = 100
#   tasks 2401-3200  p = 50
# The array index is added to OFFSET (default 0), so that blocks larger than
# the cluster's maximum array size can be submitted in pieces, e.g.
#   sbatch --array=1-800%20 --export=ALL,OFFSET=0 slurm/run_tasks.sbatch
# Resources (--time, --mem-per-cpu) are set on the command line per block,
# from the timing test (see README.md, "Submission").
# Finished tasks are skipped when resubmitted.

cd "$SLURM_SUBMIT_DIR"
TASK=$(( ${OFFSET:-0} + SLURM_ARRAY_TASK_ID ))
ROW=$(sed -n "$((TASK + 1))p" slurm/settings.csv)
IFS=, read -r N P Q REP <<< "$ROW"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export LC_ALL=C.UTF-8
source /center/healthds/singularity_functions

echo "Task $TASK: n=$N p=$P q=$Q rep=$REP on $SLURM_CPUS_PER_TASK cores ($(hostname))"
Rscript --vanilla run_task.R --n "$N" --p "$P" --q "$Q" --rep "$REP" \
  --ncores "$SLURM_CPUS_PER_TASK" --out results
echo "Task $TASK completed"
