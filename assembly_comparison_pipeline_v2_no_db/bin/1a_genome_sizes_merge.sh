#!/bin/bash

# Define the parent directory containing the assembler directories as the current directory containing the bin/ dir
current_dir=$(pwd)
# Remove '/bin' from the end of the directory path
current_prefix=${current_dir%/bin}
parent_dir="$current_prefix/genome_sizes"

# Define the output file
output_file="$current_dir/summaries/genome_sizes_merged_sup.tsv"

echo -e "sample\tgenome_size" > "$output_file"

# Change to the parent directory
cd "$parent_dir" || exit 

        # Loop through each text file in the assembler directory
        for txt_file in "$parent_dir"/*.txt; do
            # Check if the text file exists
            if [[ -f "$txt_file" ]]; then
                # Get the basename of the text file without the extension
                sample=$(basename "$txt_file" _genome_size.txt)

                # Read the TSV file line by line, skipping the header
                genome_size=$(cat "$txt_file")

                    # Prepend the assembler name and append to the output file
                    echo -e "$sample\t$genome_size" >> "$output_file"
            fi
        done

echo "Merging complete. Output written to $output_file"
