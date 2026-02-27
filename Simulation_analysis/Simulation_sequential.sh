#!/bin/bash

JOBS_LIMIT=100

# Check that an argument was provided
if [ -z "$1" ]; then
  echo "Usage: $0 <path_to_csv_file>"
  exit 1
fi

csv_file_conf="$1"

# Check if CSV file exists
if [[ ! -f "$csv_file_conf" ]]; then
    echo "Error: CSV file not found: $csv_file_conf"
    exit 1
fi

#Option 1: Process the entire CSV file (excluding header)
total_lines=$(tail -n +2 "$csv_file_conf" | wc -l)
echo "Processing all CSV data rows (total: $total_lines)"
line_num=1

# Option 2: Define the line range to process (data lines, not including header)
# start_line=257
# end_line=260
# echo "Processing lines $start_line to $end_line from CSV file (total data rows: $total_lines)"
# line_num=$start_line

tail -n +2 "$csv_file_conf" | while IFS=',' read -r p slurm_script_path input_data_path output_path name_output; do
#sed -n "$((start_line+1)),$((end_line+1))p" "$csv_file_conf" | while IFS=',' read -r p slurm_script_path input_data_path output_path name_output; do
    # Wait if job limit is reached
    while [ "$(squeue -u $USER | wc -l)" -ge "${JOBS_LIMIT}" ]; do
        echo "Jobs limit reached, sleeping for 4 minutes..."
        sleep 240
    done
    
    echo "Processing: Line ${line_num} of ${total_lines}"
    #echo "Processing: Line ${line_num} of ${end_line} (range: ${start_line}-${end_line})"
    
    # Remove any potential whitespace/quotes from variables
    p=$(echo "$p" | tr -d ' "')
    slurm_script_path=$(echo "$slurm_script_path" | tr -d ' "')
    input_data_path=$(echo "$input_data_path" | tr -d ' "')
    output_path=$(echo "$output_path" | tr -d ' "')
    name_output=$(echo "$name_output" | tr -d ' "')
    
    # Validate that required fields are not empty
    if [[ -z "$p" || -z "$slurm_script_path" || -z "$input_data_path" || -z "$output_path" || -z "$name_output" ]]; then
        echo "Warning: Skipping line ${line_num} - missing required fields"
        ((line_num++))
        continue
    fi
    
    # Submit the job
    echo "Submitting job with parameters:"
    echo "  Array size: $p"
    echo "  Script: $slurm_script_path"
    echo "  Input: $input_data_path"
    echo "  Output: $output_path"
    echo "  Name: $name_output"
    
    RES=$(sbatch --array=1-${p} "${slurm_script_path}" "${input_data_path}" "${output_path}" "${name_output}")
    
    if [[ $? -eq 0 ]]; then
        echo "Successfully submitted job: ${RES}"
    else
        echo "Error: Failed to submit job for line ${line_num}"
    fi
    
    echo "---"
    ((line_num++))
done