#!/bin/bash
#SBATCH --job-name=seltest
#SBATCH --output=logs/seltest_%A_%a.out
#SBATCH --error=logs/seltest_%A_%a.err
#SBATCH --partition=cpuq
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G
#SBATCH --time=06:00:00

# Selection-rule test (selection_test/run_selection_test.R): p = 100, q = 3,
# n in {1000, 5000, 10000, 20000}, replications 1-5, three priors per task.
#   tasks  1- 5  n = 20000      tasks 11-15  n = 5000
#   tasks  6-10  n = 10000      tasks 16-20  n = 1000
# Submit from the Simulation_analysis_v2 folder:
#   sbatch --array=1-20 selection_test/run_selection_test.sbatch
cd "$SLURM_SUBMIT_DIR"
NS=(20000 10000 5000 1000)
IDX=$(( SLURM_ARRAY_TASK_ID - 1 ))
N=${NS[$(( IDX / 5 ))]}
REP=$(( IDX % 5 + 1 ))

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export LC_ALL=C.UTF-8
source /center/healthds/singularity_functions

echo "Selection test: n=$N p=100 q=3 rep=$REP on $SLURM_CPUS_PER_TASK cores ($(hostname))"
Rscript --vanilla selection_test/run_selection_test.R --n "$N" --p 100 --q 3 --rep "$REP" \
  --ncores "$SLURM_CPUS_PER_TASK"
echo "Task $SLURM_ARRAY_TASK_ID completed"
