#!/bin/bash

# Script to run SLiM simulation in parallel
# Usage: ./run_slim_parallel.sh <slim_script> [num_cpus] [num_iterations]

# Default values
DEFAULT_CPUS=4
DEFAULT_ITERATIONS=100

# Check if script name is provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 <slim_script> [num_cpus] [num_iterations]"
    echo ""
    echo "Available SLiM scripts:"
    echo "  breakpoint_instants_additiven15.slim"
    echo "  breakpoint_instants_additiven40.slim"  
    echo "  breakpoint_instants_recessiven40.slim"
    echo ""
    echo "Examples:"
    echo "  $0 breakpoint_instants_additiven40.slim"
    echo "  $0 breakpoint_instants_additiven15.slim 8 50"
    echo "  $0 breakpoint_instants_recessiven40.slim 2 200"
    exit 1
fi

# Parse arguments
SLIM_SCRIPT=$1
NUM_CPUS=${2:-$DEFAULT_CPUS}
NUM_ITERATIONS=${3:-$DEFAULT_ITERATIONS}

# Validate SLiM script exists
if [ ! -f "$SLIM_SCRIPT" ]; then
    echo "ERROR: SLiM script '$SLIM_SCRIPT' not found!"
    echo "Make sure the file exists in the current directory."
    exit 1
fi

# Extract base name for output directory (remove .slim extension)
BASE_NAME=$(basename "$SLIM_SCRIPT" .slim)
OUTPUT_DIR="${BASE_NAME}_outputs"

# Create main output directory
mkdir -p "$OUTPUT_DIR"

echo "========================================"
echo "SLiM Parallel Simulation Runner"
echo "========================================"
echo "Script: $SLIM_SCRIPT"
echo "CPUs: $NUM_CPUS"
echo "Iterations: $NUM_ITERATIONS"
echo "Output directory: $OUTPUT_DIR"
echo ""
echo "Starting simulations..."
echo ""

# Function to run a single simulation
run_simulation() {
    local i=$1
    local script=$2
    local output_dir=$3
    
    echo "Starting iteration $i..."
    
    # Create directory for this iteration
    mkdir -p "$output_dir/iteration_$i"
    
    # Run SLiM with iteration number as argument
    # The -s flag sets a different random seed for each iteration
    # The -d flag defines the ITER constant
    slim -s $i -d ITER=$i "$script" > "$output_dir/iteration_$i/slim_output.log" 2>&1
    
    # Check if the simulation completed successfully
    if [ $? -eq 0 ]; then
        echo "Iteration $i completed successfully"
    else
        echo "ERROR: Iteration $i failed! Check $output_dir/iteration_$i/slim_output.log"
        return 1
    fi
}

# Export the function so it can be used by xargs
export -f run_simulation

# Run simulations in parallel
seq 1 $NUM_ITERATIONS | xargs -n 1 -P $NUM_CPUS -I {} bash -c "run_simulation {} '$SLIM_SCRIPT' '$OUTPUT_DIR'"

echo ""
echo "========================================"
echo "All iterations completed!"
echo "========================================"
echo ""
echo "Output structure:"
echo "$OUTPUT_DIR/"
echo "├── iteration_1/"
echo "│   ├── slim_output.log"
echo "│   ├── shift1_before_sample.vcf"
echo "│   ├── shift1_after_sample.vcf"
echo "│   ├── shift2_before_sample.vcf"
echo "│   ├── shift2_after_sample.vcf"
echo "│   ├── shift3_before_sample.vcf"
echo "│   ├── shift3_after_sample.vcf"
echo "│   └── last_sample.vcf"
echo "├── iteration_2/"
echo "│   └── ..."
echo "└── iteration_$NUM_ITERATIONS/"
echo "    └── ..."

# Optional: Show directory sizes
echo ""
echo "Directory sizes:"
du -sh "$OUTPUT_DIR"/iteration_*

# Show summary
echo ""
echo "========================================"
echo "Summary"
echo "========================================"
echo "Script: $SLIM_SCRIPT"
echo "Iterations completed: $NUM_ITERATIONS"
echo "Parallel processes used: $NUM_CPUS"
echo "Output saved to: $OUTPUT_DIR/"
echo "Each iteration contains 7 VCF files (before/after each selection shift, plus final sample)"
