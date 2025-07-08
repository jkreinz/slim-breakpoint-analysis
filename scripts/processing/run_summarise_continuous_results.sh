#!/bin/bash

# Script to summarize results from SLiM continuous simulation outputs
# Usage: ./run_summarise_results.sh <output_directory> <file_prefix> [num_iterations]

# Default values
DEFAULT_ITERATIONS=100

# Check if required arguments are provided
if [ $# -lt 2 ]; then
    echo "Usage: $0 <output_directory> <file_prefix> [num_iterations]"
    echo ""
    echo "Examples:"
    echo "  $0 breakpoint_continuous_additive_n200_outputs m2_allsamples_nomulti_iter"
    echo "  $0 breakpoint_continuous_additive_n70_outputs m2_allsamples_n70_nomulti_iter 50"
    echo "  $0 breakpoint_continuous_recessive_n200_outputs m2_allsamples_nomulti_recessive_iter 100"
    echo ""
    echo "This will create a summary CSV file in the specified output directory."
    exit 1
fi

# Parse arguments
OUTPUT_DIR=$1
FILE_PREFIX=$2
NUM_ITERATIONS=${3:-$DEFAULT_ITERATIONS}

# Validate output directory exists
if [ ! -d "$OUTPUT_DIR" ]; then
    echo "ERROR: Output directory '$OUTPUT_DIR' not found!"
    exit 1
fi

# Extract simulation type from output directory name for summary file naming
SIM_TYPE=$(basename "$OUTPUT_DIR" | sed 's/_outputs$//')
SUMMARY_FILE="${OUTPUT_DIR}/${SIM_TYPE}_summary.csv"

echo "========================================"
echo "Summarizing Results"
echo "========================================"
echo "Output directory: $OUTPUT_DIR"
echo "File prefix: $FILE_PREFIX"
echo "Iterations: $NUM_ITERATIONS"
echo "Summary file: $SUMMARY_FILE"
echo ""

# Create header for summary file
echo "iteration,best_aic_bp,best_bic_bp,breakpoint_1,breakpoint_1_se,breakpoint_2,breakpoint_2_se,breakpoint_3,breakpoint_3_se,breakpoint_4,breakpoint_4_se,s_period_1,s_period_1_se,s_period_2,s_period_2_se,s_period_3,s_period_3_se,s_period_4,s_period_4_se,s_period_5,s_period_5_se" > "$SUMMARY_FILE"

# Counter for processed iterations
processed=0
failed=0

# Process each iteration
for iter_dir in "$OUTPUT_DIR"/iteration_{1..$NUM_ITERATIONS}; do
    if [ -d "$iter_dir" ]; then
        iter_num=$(basename "$iter_dir" | sed 's/iteration_//')
        echo "Processing iteration $iter_num..."
        
        cd "$iter_dir"
        
        # Run R script with file prefix as argument
        if Rscript "$(dirname "$0")/summarise_continuous_results.R" "$iter_num" "$FILE_PREFIX" 2>/dev/null; then
            # Append results to summary file
            if [ -f "summary_${SIM_TYPE}_iter${iter_num}.csv" ]; then
                cat "summary_${SIM_TYPE}_iter${iter_num}.csv" >> "../$SUMMARY_FILE"
                ((processed++))
            else
                echo "  Warning: Summary file not created for iteration $iter_num"
                ((failed++))
            fi
        else
            echo "  Error: R script failed for iteration $iter_num"
            ((failed++))
        fi
        
        cd - > /dev/null
    else
        echo "Warning: Directory $iter_dir not found"
    fi
done

echo ""
echo "========================================"
echo "Summary Complete"
echo "========================================"
echo "Processed: $processed iterations"
echo "Failed: $failed iterations"
echo "Results saved to: $SUMMARY_FILE"
