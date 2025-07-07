#!/bin/bash

# Script to run SLiM simulation 10 times with organized output

# Create main output directory
mkdir -p instants_simulation_outputs

echo "Starting 10 iterations of the SLiM simulation..."
echo "Each iteration will be saved in simulation_outputs_n70/iteration_X/"
echo ""

# Run 10 iterations
for i in {1..2}; do
    echo "Starting iteration $i of 10..."
    
    # Create directory for this iteration
    mkdir -p instants_simulation_outputs/iteration_$i
    
    # Run SLiM with iteration number as argument
    # The -s flag sets a different random seed for each iteration
    # The -d flag defines the ITER constant
    slim -s $i -d ITER=$i breakpoint_instants.slim > instants_simulation_outputs/iteration_$i/slim_output.log 2>&1
 
    # Check if the simulation completed successfully
    if [ $? -eq 0 ]; then
        echo "Iteration $i completed successfully"
    else
        echo "ERROR: Iteration $i failed! Check instants_simulation_outputs/iteration_$i/slim_output.log"
    fi
    
    echo ""
done

echo "All iterations completed!"
echo ""
echo "Output structure:"
echo "simulation_outputs/"
echo "├── iteration_1/"
echo "│   ├── slim_output.log"
echo "│   ├── t000_sample.vcf"
echo "│   ├── t004_sample.vcf"
echo "│   └── ... (65 VCF files total)"
echo "├── iteration_2/"
echo "│   └── ..."
echo "└── iteration_10/"
echo "    └── ..."

# Optional: Show directory sizes
echo ""
echo "Directory sizes:"
du -sh simulation_outputs_n70/iteration_*
