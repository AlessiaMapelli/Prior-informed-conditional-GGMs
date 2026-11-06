#!/bin/bash
#SBATCH --job-name=addsim_test
#SBATCH --output=logs/addsim_test_%j.out
#SBATCH --error=logs/addsim_test_%j.err
#SBATCH --partition=cpuq
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=2G
#SBATCH --time=04:00:00

# Timing test: fits one node per core (nodes 1 to cpus-per-task) for each
# sample size, with all selection rules and 50 stability subsamples.
# Submit from the Add_simulation_analysis folder:
#   sbatch slurm/test_timing.sbatch            # n = 1000, 5000, 20000
#   sbatch slurm/test_timing.sbatch 20000      # a single sample size
# Afterwards, report the log and the resources used:
#   sacct -j <jobid> --format=JobID,Elapsed,TotalCPU,MaxRSS,ReqMem,AllocCPUS
# Results go to results_test/ and do not interfere with the full run.

cd "$SLURM_SUBMIT_DIR"
SIZES=${1:-"1000 5000 20000"}
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

source /center/healthds/singularity_functions

for N in $SIZES; do
  echo "=== n=$N, $SLURM_CPUS_PER_TASK cores ($(hostname))"
  time Rscript --vanilla run_setting.R --n "$N" --rep 1 --p 100 --prior noisy \
    --ncores "$SLURM_CPUS_PER_TASK" --nodes "1-$SLURM_CPUS_PER_TASK" --stab_B 50 \
    --out results_test
done
