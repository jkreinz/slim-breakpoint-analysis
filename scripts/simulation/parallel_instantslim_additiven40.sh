#!/bin/bash
# Script to run SLiM simulation in parallel using 4 CPUs
# Create main output directory
mkdir -p instants_simulation_outputs
echo "Starting iterations of the SLiM simulation in parallel (4 CPUs)..."
echo "Each iteration will be saved in instants_simulation_outputs/iteration_X/"
echo ""

# Function to run a single simulation
run_simulation() {
    local i=$1
    echo "Starting iteration $i..."
    # Create directory for this iteration
    mkdir -p instants_simulation_outputs/iteration_$i
    # Run SLiM with iteration number as argument
    # The -s flag sets a different random seed for each iteration
    # The -d flag defines the ITER constant
    slim -s $i -d ITER=$i breakpoint_instants_additiven40.slim > instants_simulation_outputs/iteration_$i/slim_output.log 2>&1
    # Check if the simulation completed successfully
    if [ $? -eq 0 ]; then
        echo "Iteration $i completed successfully"
    else
        echo "ERROR: Iteration $i failed! Check instants_simulation_outputs/iteration_$i/slim_output.log"
    fi
}

# Export the function so it can be used by xargs
export -f run_simulation

# Run simulations in parallel using 4 processes
# You can use either approach:

# Option 1: Run specific iterations (like your original)
#echo 55 3 84 81 59 15 90 37 | xargs -n 1 -P 4 -I {} bash -c 'run_simulation {}'

# Option 2: Run sequential iterations (uncomment if preferred)
 seq 1 100 | xargs -n 1 -P 2 -I {} bash -c 'run_simulation {}'

echo ""
echo "All iterations completed!"
echo ""
echo "Output structure:"
echo "instants_simulation_outputs/"
echo "├── iteration_55/"
echo "│   ├── slim_output.log"
echo "│   ├── shift1_before_sample.vcf"
echo "│   ├── shift1_after_sample.vcf"
echo "│   ├── shift2_before_sample.vcf"
echo "│   ├── shift2_after_sample.vcf"
echo "│   ├── shift3_before_sample.vcf"
echo "│   └── shift3_after_sample.vcf"
echo "├── iteration_3/"
echo "│   └── ..."
echo "└── iteration_37/"
echo "    └── ..."

# Optional: Show directory sizes
echo ""
echo "Directory sizes:"
du -sh instants_simulation_outputs/iteration_*

# Show total time and summary
echo ""
echo "Summary: Simulations completed using 4 parallel processes"
echo "Each iteration contains 6 VCF files (before/after each selection shift)"
