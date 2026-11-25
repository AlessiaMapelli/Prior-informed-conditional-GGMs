#!/bin/bash
#SBATCH --job-name=prcggm_T2D
#SBATCH --output=logs/prcggm_node_%A_%a.out
#SBATCH --error=logs/prcggm_node_%A_%a.err
#SBATCH --partition=cpuq
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=<your email>

# prCGGM (>= 0.2.0) node-wise estimation for the T2D case study: one array task
# per protein. Submitted by Network_estimation_case_study_T2D_UKB.R from the
# project folder:
#   sbatch --array=1-366 ../T2D_application_code/slurm_prcggm_node.sbatch <input.rds> <output_dir>
# Resource guess for n = 24,446 and P = 1,460 predictors (16 hyperparameter
# combinations, no CV): check the first tasks with
#   sacct -j <jobid> --format=JobID,Elapsed,MaxRSS
# and resubmit unfinished nodes with a larger --time/--mem if needed.

INPUT_FILE=$1
OUTPUT_DIR=$2
NODE_INDEX=$SLURM_ARRAY_TASK_ID

cd "$SLURM_SUBMIT_DIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export LC_ALL=C.UTF-8
# Rscript --vanilla ignores ~/.Renviron, so the library is set here
# Optional personal library holding prCGGM >= 0.2.0 (export PRCGGM_LIB=... before sbatch)
PRCGGM_LIB=${PRCGGM_LIB:-}
T2D_CODE_DIR=${T2D_CODE_DIR:-../T2D_application_code}
export T2D_CODE_DIR SINGULARITYENV_T2D_CODE_DIR=$T2D_CODE_DIR APPTAINERENV_T2D_CODE_DIR=$T2D_CODE_DIR
export R_LIBS=$PRCGGM_LIB
# also passed explicitly into the container, in case it starts with a clean
# environment
export SINGULARITYENV_R_LIBS=$PRCGGM_LIB APPTAINERENV_R_LIBS=$PRCGGM_LIB
# Cluster-specific R set-up (e.g. a script that makes Rscript available in a
# container): export R_SETUP=<file> before sbatch
if [ -n "${R_SETUP:-}" ]; then source "$R_SETUP"; fi

echo "prCGGM: node $NODE_INDEX on $(hostname), $(date)"
# The library is also set inside R, so the job does not depend on R_LIBS
Rscript --vanilla -e "if (nzchar('$PRCGGM_LIB')) .libPaths(c('$PRCGGM_LIB', .libPaths())); cat('prCGGM', as.character(packageVersion('prCGGM')), '\n'); prCGGM::prcggm_node(input = '$INPUT_FILE', node = $NODE_INDEX, output_dir = '$OUTPUT_DIR')"
echo "prCGGM: node $NODE_INDEX completed, $(date)"
