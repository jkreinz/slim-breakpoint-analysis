# SLiM Breakpoint Detection Analysis

Population genetics simulations using SLiM to evaluate breakpoint detection methods for temporal genomic data.

## Overview

This project implements and evaluates methods for detecting breakpoints in selection pressure using temporal sampling of genomic data. Despite varying sample size, scheme, and dominance, all simulations take the same form - after a burn in of 2000 generations where all mutations are evolving neturally, the population experiences shifts in selection at t=2000 to s=0.05, at t=2050 to s=-0.025, and at t=2100 to s=0.10 and continues to evolve until t=2141. The downstream analysis includes two distinct pipelines:

1. **Instantaneous Selection Pipeline**: Evaluates instanteous selection coefficient changes at discrete time points
2. **Continuous Selection Pipeline**: Analyzes selection dynamics over continuous time series with logistic regression

### SLIM parameter space
- **Multiple Sampling Strategies**: Varying sample sizes (n=15 x 6 timepoints or n=40 x 6 timepoints for instantaneous; n=70 or n=200 for continuous sampling through time, in a fashion similar to herbarium/museum collection availability)
- **Genetic Architectures**: Additive (h=0.5) and recessive (h=0.15) dominance models

## Pipeline Structure

### 1. Instantaneous Selection Pipeline
Analyzes selection coefficient changes at specific breakpoints with discrete sampling, and estimates selection coefficients around known shifts in selection.

**Simulation Types:**
- `breakpoint_instants_additiven15.slim` - Additive model, n=15 samples
- `breakpoint_instants_additiven40.slim` - Additive model, n=40 samples  
- `breakpoint_instants_recessiven40.slim` - Recessive model, n=40 samples

**Output Structure:**
```
breakpoint_instants_[type]_outputs/
├── iteration_1/
│   ├── shift1_before_sample.vcf
│   ├── shift1_after_sample.vcf
│   ├── shift2_before_sample.vcf
│   ├── shift2_after_sample.vcf
│   ├── shift3_before_sample.vcf
│   ├── shift3_after_sample.vcf
│   └── last_sample_n15_h5.vcf
└── iteration_N/
    └── ...
```

### 2. Continuous Selection Pipeline
Analyzes selection dynamics through continuous temporal sampling, along with power to infer shifts (breakpoints) in selection, and the accuracy, and bias of selection coefficients and timing of selective shifts.


**Simulation Types:**
- `breakpoint_continuous_additive_n200.slim` - Additive model, n=200 population
- `breakpoint_continuous_additive_n70.slim` - Additive model, n=70 population
- `breakpoint_continuous_recessive_n200.slim` - Recessive model, n=200 population

**Output Structure:**
```
breakpoint_continuous_[type]_outputs/
├── iteration_1/
│   ├── t128_sample_h50_n200.vcf
│   ├── t130_sample_h50_n200.vcf
│   ├── t132_sample_h50_n200.vcf
│   └── ...
└── iteration_N/
    └── ...
```

## Quick Start

### Prerequisites
- SLiM 4.0+
- R 4.0+ with required packages:
  ```r
  install.packages(c("vcfR", "dplyr", "ggplot2", "purrr", "tidyr", 
                     "viridis", "patchwork", "segmented", "data.table"))
  ```
- vcftools
- bgzip/tabix

### Running Instantaneous Selection Analysis

```bash
# 1. Run parallel simulations
./run_slim_parallel.sh breakpoint_instants_additiven40.slim <num cpu> <num iterations> #Runs instantaneous selection simulations in parallel
./run_slim_parallel.sh breakpoint_instants_additiven15.slim <num cpu> <num iterations> #Runs instantaneous selection simulations in parallel
./run_slim_parallel.sh breakpoint_instants_recessive40.slim <num cpu> <num iterations> #Runs instantaneous selection simulations in parallel

# 2. Process VCF files (filter for m2 mutations)
./process_vcf_files.sh breakpoint_instants_additiven40_outputs <num iterations>
./process_vcf_files.sh breakpoint_instants_additiven15_outputs <num iterations>
./process_vcf_files.sh breakpoint_instants_recessive40_outputs <num iterations>

# 3. Run analysis and visualization
Rscript breakpoint_instantaneousS_analysis.R
```

### Running Continuous Selection Analysis

```bash
# 1. Run parallel simulations
./run_slim_continuous_parallel.sh breakpoint_continuous_additive_n200.slim <num cpu> <num iterations> #Runs continuous selection simulations in parallel
./run_slim_continuous_parallel.sh breakpoint_continuous_additive_n70.slim <num cpu> <num iterations> #Runs continuous selection simulations in parallel
./run_slim_continuous_parallel.sh breakpoint_continuous_recessive_n200.slim <num cpu> <num iterations> #Runs continuous selection simulations in parallel

# 2. Process VCF files (filter for m2 mutations)  
./process_vcf_files.sh breakpoint_continuous_additive_n200_outputs <num iterations>
./process_vcf_files.sh breakpoint_continuous_recessive_n200_outputs <num iterations>
./process_vcf_files.sh breakpoint_continuous_additive_n70_outputs <num iterations>

# 3. Generate summary statistics (breakpoint model statistics for 1-5 breakpoints)
./run_summarise_results.sh breakpoint_continuous_additive_n200_outputs m2_allsamples_nomulti_iter <num iterations>
./run_summarise_results.sh breakpoint_continuous_recessive_n200_outputs m2_allsamples_nomulti_iter <num iterations>
./run_summarise_results.sh breakpoint_continuous_additive_n70_outputs m2_allsamples_nomulti_iter <num iterations>

# 4. Run comprehensive analysis and produce figures
Rscript breakpointslimulations_continuous_analysis_updatednames.R
```


## File Organization

```
project/
├── README.md
├── docs/
├── results/
└── scripts/
    ├── analysis_and_visualization/
    │   ├── breakpoint_instantaneousS_analysis.R           # Instantaneous analysis
    │   └── breakpointslimulations_continuous_analysis_updatednames.R  # Continuous analysis
    ├── processing/
    │   ├── process_slim_outout.sh                         # VCF processing (both pipelines)
    │   ├── run_summarise_continuous_results.sh           # Continuous summarization
    │   └── summarise_continuous_results.R                # Summary statistics
    └── simulation/
        ├── breakpoint_instants_additiven15.slim          # Instantaneous SLiM script
        ├── breakpoint_instants_additiven40.slim          # Instantaneous SLiM script
        ├── breakpoint_instants_recessiven40.slim         # Instantaneous SLiM script
        ├── breakpoint_continuous_additive_n200.slim      # Continuous SLiM script
        ├── breakpoint_continuous_additive_n70.slim       # Continuous SLiM script
        ├── breakpoint_continuous_recessive_n200.slim     # Continuous SLiM script
        ├── parallel_continuous_slim.sh                   # Continuous simulations
        └── parallel_instantS_slim.sh                     # Instantaneous simulations

# Output directories (created after running simulations):
outputs/
├── breakpoint_instants_*_outputs/                        # Instantaneous results
└── breakpoint_continuous_*_outputs/                      # Continuous results
