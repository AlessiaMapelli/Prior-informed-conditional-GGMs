#!/bin/bash
#SBATCH --job-name=seltest2
#SBATCH --output=logs/seltest2_%A_%a.out
#SBATCH --error=logs/seltest2_%A_%a.err
#SBATCH --partition=cpuq
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G
#SBATCH --time=02:00:00

# Part 2 of the sensitivity analysis (selection_test/run_selection_test.R
# --stab 0): BIC (penalised and relaxed) and EBIC relaxed, current vs extended
# grid, p = 100, q = 3, three priors, 100 replications per n (the datasets of
# the main study). Without stability selection.
#   tasks   1-100  n = 20000      tasks 201-300  n = 5000
#   tasks 101-200  n = 10000      tasks 301-400  n = 1000
# Submit from the Simulation_analysis_v2 folder:
#   sbatch --array=1-400%20 selection_test/run_part2.sbatch
# Finished tasks are skipped when resubmitted.
cd "$SLURM_SUBMIT_DIR"
NS=(20000 10000 5000 1000)
IDX=$(( SLURM_ARRAY_TASK_ID - 1 ))
N=${NS[$(( IDX / 100 ))]}
REP=$(( IDX % 100 + 1 ))

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export LC_ALL=C.UTF-8
source /center/healthds/singularity_functions

echo "Part 2: n=$N p=100 q=3 rep=$REP on $SLURM_CPUS_PER_TASK cores ($(hostname))"
Rscript --vanilla selection_test/run_selection_test.R --n "$N" --p 100 --q 3 --rep "$REP" \
  --ncores "$SLURM_CPUS_PER_TASK" --stab 0 --out selection_test/results_100
echo "Task $SLURM_ARRAY_TASK_ID completed"
