#!/bin/bash

current_dir=$(pwd)
# Remove '/bin' from the end of the directory path
current_prefix=${current_dir%/bin}



# Define the output file
output_file="$current_dir/summaries/amrfinderplus_merged_sup.tsv"

# Initialize a variable to keep track of whether the header has been written
header_written=false

# Create or clear the output file
> "$output_file"

# Loop through all the .tsv files in the AMRFinderPlus subdirectories
for file in $current_prefix/AMRFinderPlus/*/*.tsv; do
    # Extract the assembler and sample names from the file path
    assembler=$(basename "$(dirname "$file")")
    sample=$(basename "$file" .tsv)
    
    # Read the file and process its content
    while IFS= read -r line || [ -n "$line" ]; do
        if [ "$header_written" = false ]; then
            # Write the header to the output file
            echo -e "assembler\tsample\t$line" >> "$output_file"
            header_written=true
        else
            # Skip the header line in subsequent files
            if [[ "$line" != $(head -n 1 "$file") ]]; then
                # Write the content to the output file with added columns
                echo -e "$assembler\t$sample\t$line" >> "$output_file"
            fi
        fi
    done < "$file"
done

echo "Combined TSV file has been created: $output_file"
