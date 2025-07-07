#!/bin/bash

for iter_dir in simulation_outputs_recessive/iteration_{55,3,84,81,59,15,90,37}; do
    if [ -d "$iter_dir" ]; then
        iter_num=$(basename "$iter_dir" | sed 's/iteration_//')
        cd "$iter_dir"
        
        for vcf in t*_sample.vcf; do
            if [ -f "$vcf" ]; then
                bgzip -f "$vcf" && tabix -p vcf "${vcf}.gz"
            fi
        done
        
        vcf-merge t*_sample.vcf.gz > allsamples_iter${iter_num}.vcf
        
        cat <(grep "^#" allsamples_iter${iter_num}.vcf) <(grep -v "^#" allsamples_iter${iter_num}.vcf | grep 'MT=2') | \
            grep -v "MULTIALLELIC" > m2_allsamples_iter${iter_num}.vcf
        
        sed 's/\t\./\t0|0/g' m2_allsamples_iter${iter_num}.vcf > m2_allsamples_homoref_iter${iter_num}.vcf
        
        grep -v "^#" m2_allsamples_homoref_iter${iter_num}.vcf | \
            awk -F'\t' '{print $8}' | \
            awk -F';' '{for(i=1;i<=NF;i++) if($i~/^S=/) print substr($i,3)}' > selection_coefficients_iter${iter_num}.txt
        
        vcftools --vcf m2_allsamples_homoref_iter${iter_num}.vcf --012 --out m2_allsamples_nomulti_recessive_iter${iter_num}
        
        cd - > /dev/null
    fi
done
