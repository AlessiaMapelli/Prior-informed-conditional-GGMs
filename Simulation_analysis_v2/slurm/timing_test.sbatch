#!/bin/bash
#SBATCH --job-name=sim2_timing
#SBATCH --output=logs/timing_%A_%a.out
#SBATCH --error=logs/timing_%A_%a.err
#SBATCH --partition=cpuq
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=4G
#SBATCH --time=03:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=alessia.mapelli@fht.org

# Timing test of the simulation study: one task per row of
# slurm/settings_timing.csv (p = 500, 250 at n = 1000 and 20000, and p = 100
# at n = 20000, q = 3), all three priors, nodes 1-16 only (one node per core),
# so each prior takes about the time of one node.
# Submit from the Simulation_analysis_v2 folder:
#   sbatch --array=1-5 slurm/timing_test.sbatch
# Afterwards:
#   sacct -j <jobid> --format=JobID,Elapsed,TotalCPU,MaxRSS,ReqMem,AllocCPUS,State
#   grep -h "Data generated\|Prior" logs/timing_<jobid>_*.out

cd "$SLURM_SUBMIT_DIR"
ROW=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" slurm/settings_timing.csv)
IFS=, read -r N P Q REP <<< "$ROW"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export LC_ALL=C.UTF-8
source /center/healthds/singularity_functions

echo "Timing task $SLURM_ARRAY_TASK_ID: n=$N p=$P q=$Q on $SLURM_CPUS_PER_TASK cores ($(hostname))"
time Rscript --vanilla run_task.R --n "$N" --p "$P" --q "$Q" --rep "$REP" \
  --ncores "$SLURM_CPUS_PER_TASK" --nodes 1-16 --out results_timing
