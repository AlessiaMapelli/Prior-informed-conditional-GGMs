#!/bin/bash
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=aalessia.mapelli@fht.org
#SBATCH --job-name=job_name
#SBATCH --output=logs/ggReg_node_%A_%a.out
#SBATCH --error=logs/ggReg_node_%A_%a.err
#SBATCH --time=01:00:00
#SBATCH --mem=2GB
#SBATCH --partition=cpuq

# SLURM array script for GGReg node-wise estimation
# Each array task processes one node i
# Usage: sbatch --array=1-p slurm_ggReg_node.sbatch input_data_path output_path name_output

# Create logs directory if it doesn't exist
mkdir -p logs

# Parse command line arguments
INPUT_DATA_PATH=$1
OUTPUT_PATH=$2
NAME_OUTPUT=$3

# Get the array task ID (this is the node index i)
NODE_INDEX=$SLURM_ARRAY_TASK_ID

# Load required modules (adjust based on your cluster setup)
source /center/healthds/singularity_functions

echo "Starting job for node $NODE_INDEX"
echo "Input data: $INPUT_DATA_PATH"
echo "Output path: $OUTPUT_PATH"

# Run the R script for this specific node
Rscript --vanilla ggReg_node_job.R $INPUT_DATA_PATH $NODE_INDEX $OUTPUT_PATH $NAME_OUTPUT

echo "Job completed for node $NODE_INDEX"