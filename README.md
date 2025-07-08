# SLiM Breakpoint Detection Analysis

Population genetics simulations using SLiM to evaluate breakpoint detection methods for temporal genomic data.

## Overview

This project implements and evaluates methods for detecting breakpoints in selection pressure using temporal sampling of genomic data. The analysis includes:

- SLiM simulations with varying selection regimes
- Multiple sampling strategies (n=15, n=40, n=200)
- Additive and recessive genetic architectures
- Statistical analysis of detection power and accuracy

## Key Features

- **Simulation Framework**: Parallel SLiM simulations with configurable parameters
- **Data Processing**: Automated VCF processing and filtering pipelines
- **Visualization**: Comprehensive R scripts for trajectory and bias analysis
- **Statistical Analysis**: Breakpoint detection accuracy and selection coefficient estimation

## Quick Start

### Prerequisites
- SLiM 4.0+
- R 4.0+ with required packages (vcfR, dplyr, ggplot2, purrr, tidyr)
- vcftools
- bgzip/tabix

### Running Simulations
```bash
# Run additive simulations (n=15)
bash scripts/simulation/parallel_instantslim_additiven15.sh

# Process VCF outputs
bash scripts/processing/process_instants_slim.sh

# Generate visualizations
Rscript scripts/visualization/instanteouss_visualize_final.R
