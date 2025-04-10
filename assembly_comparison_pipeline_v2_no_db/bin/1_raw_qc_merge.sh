#!/bin/bash

# Define the parent directory containing the assembler directories as the current directory containing the bin/ dir
current_dir=$(pwd)
# Remove '/bin' from the end of the directory path
current_prefix=${current_dir%/bin}
parent_dir="$current_prefix/raw_QC_CSVs"

# Define the output file
output_file="$current_dir/summaries/raw_QC_merged_sup.tsv"

# Initialize a variable to track whether the header has been written
header_written=false

# Change to the parent directory
cd "$parent_dir" || exit 

# Loop through each directory in the parent directory
for assembler in */; do
    # Remove the trailing slash
    assembler=${assembler%/}
    
    # Check if the directory is not empty
    if [[ -d "$assembler" ]]; then
        # Loop through each TSV file in the assembler directory
        for tsv_file in "$assembler"/*.tsv; do
            # Check if the TSV file exists
            if [[ -f "$tsv_file" ]]; then
                # Get the basename of the TSV file without the extension
                file_basename=$(basename "$tsv_file" .tsv)
                # Get the header from the TSV file
                header=$(head -n 1 "$tsv_file")
                if [[ "$header_written" == false ]]; then
                    echo -e "assembler\tsample\t$header" > "$output_file"
                    header_written=true
                fi
                
                # Read the TSV file line by line, skipping the header
                tail -n +2 "$tsv_file" | while IFS= read -r line; do
                    # Prepend the assembler name and append to the output file
                    echo -e "$assembler\t$file_basename\t$line" >> "$output_file"
                done
            fi
        done
    fi
done

echo "Merging complete. Output written to $output_file"
