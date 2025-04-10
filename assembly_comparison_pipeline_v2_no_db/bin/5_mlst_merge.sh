#!/bin/bash


current_dir=$(pwd)
# Remove '/bin' from the end of the directory path
current_prefix=${current_dir%/bin}
parent_dir="$current_prefix/MLST/"

# Define the output file
output_file="$current_dir/summaries/MLST_merged.tsv"

# Initialize a variable to track whether the header has been written
header="assembler\tsample\tspecies\tmlst\tgene1\tgene2\tgene3\tgene4\tgene5\tgene5_2\tgene6"

# Change to the parent directory
cd "$parent_dir" || exit

# Write the custom header to the output file
echo -e "$header" > "$output_file"

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

                # Read the single line in the TSV file and remove the first column
                line=$(cat "$tsv_file" | cut -f2-)

                # Prepend the assembler name and file basename, then append to the output file
                echo -e "$assembler\t$file_basename\t$line" >> "$output_file"
            fi
        done
    fi
done

echo "Merging complete. Output written to $output_file"
