#!/bin/bash
#SBATCH --job-name=t2d_perm
#SBATCH --output=logs/perm_%A_%a.out
#SBATCH --error=logs/perm_%A_%a.err
#SBATCH --partition=cpuq
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=08:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=<your email>

# Permutation null for the T2D-specific network (T2D_permutation_null.R).
# One array task = one permutation (mean step + 366 node refits; ~3 h for the
# 5-year fit, n = 24,446; somewhat more for the 10-year fit, n = 32,439;
# task t -> permutation t - 1, permutation 0 = real labels). Submit from the
# project folder after
#   Rscript --vanilla ../T2D_application_code/T2D_permutation_null.R prepare 100
# with (B = 100 -> 101 tasks):
#   sbatch --array=1-101%50 ../T2D_application_code/slurm_permutation_null.sbatch
# Finished chunks are skipped when resubmitted.
# The fit to test is set by T2D_RES (default results_bic_relaxed_10y), e.g.
#   sbatch --array=1-101%50 --export=ALL,T2D_RES=results_bic_relaxed ../T2D_application_code/slurm_permutation_null.sbatch

cd "$SLURM_SUBMIT_DIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export LC_ALL=C.UTF-8
# Optional personal library holding prCGGM >= 0.2.0 (export PRCGGM_LIB=... before sbatch)
PRCGGM_LIB=${PRCGGM_LIB:-}
T2D_CODE_DIR=${T2D_CODE_DIR:-../T2D_application_code}
export T2D_CODE_DIR SINGULARITYENV_T2D_CODE_DIR=$T2D_CODE_DIR APPTAINERENV_T2D_CODE_DIR=$T2D_CODE_DIR
export R_LIBS=$PRCGGM_LIB SINGULARITYENV_R_LIBS=$PRCGGM_LIB APPTAINERENV_R_LIBS=$PRCGGM_LIB
export T2D_RES=${T2D_RES:-results_bic_relaxed_10y}
export SINGULARITYENV_T2D_RES=$T2D_RES APPTAINERENV_T2D_RES=$T2D_RES
# Cluster-specific R set-up (e.g. a script that makes Rscript available in a
# container): export R_SETUP=<file> before sbatch
if [ -n "${R_SETUP:-}" ]; then source "$R_SETUP"; fi

echo "task $SLURM_ARRAY_TASK_ID on $(hostname), $(date), results $T2D_RES"
Rscript --vanilla "$T2D_CODE_DIR/T2D_permutation_null.R" node "$SLURM_ARRAY_TASK_ID"
echo "done $(date)"
