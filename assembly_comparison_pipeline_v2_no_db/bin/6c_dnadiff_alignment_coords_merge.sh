#!/bin/bash
# Define the parent directory containing the assembler directories as the current directory containing the bin
current_dir=$(pwd)

# Remove '/bin' from the end of the directory path
current_prefix=${current_dir%/bin}
parent_dir="$current_prefix/dnadiff"

output_file=$current_dir/summaries/dnadiff_unique_coords.tsv

# Output TSV file
echo -e "assembler\tsample\tautocycler_ref_contig\tqry_contig" > $output_file

# Find all .1coords files matching the pattern
find "$current_prefix/dnadiff" -type f -name "*.1coords" | while read file; do
    # Extract assembler and sample names from the path
    assembler=$(echo "$file" | awk -F'/' '{print $(NF-1)}')
    sample=$(echo "$file" | awk -F'/' '{print $(NF-0)}' | sed 's/.1coords$//')
    
    # Extract unique column 12 and 13 pairs
    awk '{print $12, $13}' "$file" | sort -u | while read col12 col13; do
        echo -e "$assembler\t$sample\t$col12\t$col13" >> $output_file
    done
done
