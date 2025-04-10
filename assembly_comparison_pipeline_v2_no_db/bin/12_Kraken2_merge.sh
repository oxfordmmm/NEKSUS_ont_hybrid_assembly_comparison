#!/bin/bash

# Define the output file
current_dir=$(pwd)
current_prefix=${current_dir%/bin}
output_file=${current_prefix}/bin/summaries/12_kraken_summary.tsv

# Clear the output file if it exists
> "$output_file"

echo -e "assembler\tsample\tcontig\tlength\tcontig_species\tmain_chromosome_species" > "$output_file"

# Loop through each subdirectory (assembler) in the current/directory
for assembler_dir in ${current_prefix}/kraken2/*/; do
    # Check if it is a directory
    if [ -d "$assembler_dir" ]; then
        assembler=$(basename "$assembler_dir")

        # Loop through each {sample}_k2_output.tsv file in the assembler directory
        for tsv_file in "$assembler_dir"*"_k2_output.tsv"; do
            # Extract the sample name from the file name (before the _k2_output.tsv)
            sample=$(basename "$tsv_file" "_k2_output.tsv")

            # Read the corresponding main_chromosome_species.txt file for the sample
            main_chromosome_species_file="$assembler_dir${sample}_chromosome_species.txt"
            
            if [ -f "$main_chromosome_species_file" ]; then
                # Read the contents of the main_chromosome_species.txt file
                main_chromosome_species=$(<"$main_chromosome_species_file")
            else
                main_chromosome_species="N/A"  # If the file is not found, set it to "N/A"
            fi

            # Process each file line-by-line
            while IFS=$'\t' read -r completeness contig contig_species length _; do
                # Append the required data to the output file, including the new column
                echo -e "$assembler\t$sample\t$contig\t$length\t$contig_species\t$main_chromosome_species" >> "$output_file"
            done < "$tsv_file"
        done
    fi
done

echo "Concatenated TSV file has been created: $output_file"
