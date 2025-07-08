#!/bin/bash

# Script to process VCF files from SLiM simulation outputs
# Usage: ./process_vcf_files.sh <output_directory> [num_iterations]

# Default values
DEFAULT_ITERATIONS=100

# Check if output directory is provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 <output_directory> [num_iterations]"
    echo ""
    echo "Examples:"
    echo "  $0 breakpoint_instants_additiven40_outputs"
    echo "  $0 breakpoint_continuous_additive_n70_outputs 50"
    echo ""
    echo "This script will process all VCF files in the specified directory structure,"
    echo "filtering for m2 mutations only and excluding multiallelic variants."
    exit 1
fi

# Parse arguments
OUTPUT_DIR=$1
NUM_ITERATIONS=${2:-$DEFAULT_ITERATIONS}

# Validate output directory exists
if [ ! -d "$OUTPUT_DIR" ]; then
    echo "ERROR: Output directory '$OUTPUT_DIR' not found!"
    echo "Make sure you've run the simulation script first."
    exit 1
fi

echo "========================================"
echo "VCF Processing Script"
echo "========================================"
echo "Processing directory: $OUTPUT_DIR"
echo "Expected iterations: $NUM_ITERATIONS"
echo ""
echo "Filtering criteria:"
echo "  - MT=2 mutations only"
echo "  - Exclude multiallelic variants (no comma in ALT)"
echo ""

# Counter for processed files
total_processed=0
total_iterations=0

# Process each iteration directory
for iter_dir in "$OUTPUT_DIR"/iteration_{1..$NUM_ITERATIONS}; do
    if [ -d "$iter_dir" ]; then
        iter_num=$(basename "$iter_dir" | sed 's/iteration_//')
        echo "Processing iteration $iter_num"
        cd "$iter_dir"
        
        # Counter for this iteration
        files_processed=0
        
        # Process each VCF file (look for various naming patterns)
        for vcf in *.vcf; do
            if [ -f "$vcf" ]; then
                # Skip already processed files
                if [[ "$vcf" == m2_* ]]; then
                    continue
                fi
                
                output_vcf="m2_${vcf}"
                
                # Keep header and filter for MT=2, exclude multiallelic (no comma in ALT)
                cat <(grep "^#" "$vcf") <(grep -v "^#" "$vcf" | grep 'MT=2' | awk -F'\t' '$5 !~ /,/') > "$output_vcf"
                
                echo "  Filtered $vcf -> $output_vcf"
                ((files_processed++))
                ((total_processed++))
            fi
        done
        
        if [ $files_processed -eq 0 ]; then
            echo "  WARNING: No VCF files found in $iter_dir"
        else
            echo "  Processed $files_processed VCF files in iteration $iter_num"
        fi
        
        cd - > /dev/null
        ((total_iterations++))
    else
        echo "WARNING: Directory $iter_dir not found"
    fi
done

echo ""
echo "========================================"
echo "Processing Summary"
echo "========================================"
echo "Iterations processed: $total_iterations"
echo "Total VCF files processed: $total_processed"
echo ""

# Show some statistics about the processed files
echo "Sample of processed files:"
find "$OUTPUT_DIR" -name "m2_*.vcf" | head -10

echo ""
echo "Directory sizes after processing:"
du -sh "$OUTPUT_DIR"/iteration_* | head -5
if [ $total_iterations -gt 5 ]; then
    echo "... (showing first 5 of $total_iterations iterations)"
fi
