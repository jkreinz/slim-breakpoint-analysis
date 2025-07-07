#!/bin/bash

# Define the sample columns to keep for 70-sample subset
# Based on the reduced sampling scheme: early sparse, then mostly 1 per timepoint, few duplicates in recent periods
# Column numbers (1-indexed, accounting for first 9 VCF columns): 
# 10=i0, 12=t008_sample_i0, 15=t016_sample_i0, etc.

for iter_dir in simulation_outputs/iteration_{1..100}; do
    if [ -d "$iter_dir" ]; then
        iter_num=$(basename "$iter_dir" | sed 's/iteration_//')
        cd "$iter_dir"
        
        # Check if the original m2_allsamples file exists
        if [ -f "m2_allsamples_iter${iter_num}.vcf" ]; then
            echo "Processing iteration ${iter_num}..."
            
            # Create list of individuals to keep (70 total)
            cat > individuals_n70.txt << 'EOF'
t004_sample_i0
t008_sample_i0
t016_sample_i0
t024_sample_i0
t032_sample_i0
t040_sample_i0
t048_sample_i0
t056_sample_i0
t060_sample_i0
t062_sample_i0
t064_sample_i0
t066_sample_i0
t068_sample_i0
t070_sample_i0
t072_sample_i0
t074_sample_i0
t076_sample_i0
t078_sample_i0
t080_sample_i0
t082_sample_i0
t084_sample_i0
t086_sample_i0
t088_sample_i0
t090_sample_i0
t092_sample_i0
t094_sample_i0
t096_sample_i0
t098_sample_i0
t100_sample_i0
t102_sample_i0
t104_sample_i0
t106_sample_i0
t108_sample_i0
t110_sample_i0
t112_sample_i0
t114_sample_i0
t116_sample_i0
t118_sample_i0
t118_sample_i1
t120_sample_i0
t122_sample_i0
t122_sample_i1
t124_sample_i0
t126_sample_i0
t126_sample_i1
t128_sample_i0
t130_sample_i0
t130_sample_i1
t132_sample_i0
t134_sample_i0
t134_sample_i1
t136_sample_i0
t138_sample_i0
t138_sample_i1
t140_sample_i0
t141_sample_i0
t141_sample_i1
t139_sample_i0
t137_sample_i0
t137_sample_i1
t137_sample_i2
t137_sample_i3
t139_sample_i1
t139_sample_i2
t139_sample_i3
t141_sample_i2
t141_sample_i3
t140_sample_i1
t140_sample_i2
t140_sample_i3
EOF
            
            # Subsample the VCF to 70 individuals using vcftools
            vcftools --vcf m2_allsamples_iter${iter_num}.vcf \
                     --keep individuals_n70.txt \
                     --recode --recode-INFO-all \
                     --out m2_allsamples_n70_iter${iter_num}
            
            # Rename the output file
            mv m2_allsamples_n70_iter${iter_num}.recode.vcf m2_allsamples_n70_iter${iter_num}.vcf
            
            # Apply the same processing as the original script
            sed 's/\t\./\t0|0/g' m2_allsamples_n70_iter${iter_num}.vcf > m2_allsamples_n70_homoref_iter${iter_num}.vcf
            
            # Extract selection coefficients (same as original)
            grep -v "^#" m2_allsamples_n70_homoref_iter${iter_num}.vcf | \
                awk -F'\t' '{print $8}' | \
                awk -F';' '{for(i=1;i<=NF;i++) if($i~/^S=/) print substr($i,3)}' > selection_coefficients_n70_iter${iter_num}.txt
            
            # Create 012 format files
            vcftools --vcf m2_allsamples_n70_homoref_iter${iter_num}.vcf --012 --out m2_allsamples_n70_nomulti_iter${iter_num}
            
            # Clean up temporary file
            rm individuals_n70.txt
            
            echo "Completed iteration ${iter_num} - subsampled to 70 individuals"
        else
            echo "File m2_allsamples_iter${iter_num}.vcf not found in $iter_dir"
        fi
        
        cd - > /dev/null
    fi
done

echo "Subsampling complete. All files have '_n70' appended to distinguish from original 200-sample files."
