for iter_dir in instants_simulation_outputs/iteration_{1..100}; do
    if [ -d "$iter_dir" ]; then
        iter_num=$(basename "$iter_dir" | sed 's/iteration_//')
        echo "Processing iteration $iter_num"
        cd "$iter_dir"
        
        # Process each VCF file to filter for m2 mutations only
        for vcf in *_sample_n15_h5.vcf; do
            if [ -f "$vcf" ]; then
                output_vcf="m2_${vcf}"
                # Keep header and filter for MT=2, exclude multiallelic (no comma in ALT)
                cat <(grep "^#" "$vcf") <(grep -v "^#" "$vcf" | grep 'MT=2' | awk -F'\t' '$5 !~ /,/') > "$output_vcf"
                echo "  Filtered $vcf -> $output_vcf"
            fi
        done
        
        cd - > /dev/null
    fi
done

