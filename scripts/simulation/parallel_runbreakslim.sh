#!/bin/bash

# Script to run SLiM simulation 10 times in parallel using 4 CPUs

# Create main output directory
mkdir -p simulation_outputs_recessive

echo "Starting 10 iterations of the SLiM simulation in parallel (4 CPUs)..."
echo "Each iteration will be saved in simulation_outputs_recessive/iteration_X/"
echo ""

# Function to run a single simulation
run_simulation() {
    local i=$1
    echo "Starting iteration $i..."
    
    # Create directory for this iteration
    mkdir -p simulation_outputs_recessive/iteration_$i
    
    # Run SLiM with iteration number as argument
    # The -s flag sets a different random seed for each iteration
    # The -d flag defines the ITER constant
    slim -s $i -d ITER=$i breakpoint_k3_h010.slim > simulation_outputs_recessive/iteration_$i/slim_output.log 2>&1
    
    # Check if the simulation completed successfully
    if [ $? -eq 0 ]; then
        echo "Iteration $i completed successfully"
    else
        echo "ERROR: Iteration $i failed! Check simulation_outputs_recessive/iteration_$i/slim_output.log"
    fi
}

# Export the function so it can be used by xargs
export -f run_simulation

# Run simulations in parallel using 4 processes
# seq 1 10 generates numbers 1 through 10
# xargs -n 1 takes one number at a time
# -P 4 runs up to 4 processes in parallel
# -I {} replaces {} with the input number
#seq 1 100 | xargs -n 1 -P 5 -I {} bash -c 'run_simulation {}'
echo 55 3 84 81 59 15 90 37 | xargs -n 1 -P 5 -I {} bash -c 'run_simulation {}'

echo ""
echo "All iterations completed!"
echo ""
echo "Output structure:"
echo "simulation_outputs_recessive/"
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
du -sh simulation_outputs_recessive/iteration_*

# Show total time and summary
echo ""
echo "Summary: 10 simulations completed using 4 parallel processes"
