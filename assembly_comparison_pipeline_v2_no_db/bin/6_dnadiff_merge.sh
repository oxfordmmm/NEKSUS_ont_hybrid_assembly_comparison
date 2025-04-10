#!/bin/bash

current_dir=$(pwd)
# Remove '/bin' from the end of the directory path
current_prefix=${current_dir%/bin}

input_path="$current_prefix/dnadiff/*/*.report"

# Define the output TSV file
output_file="$current_dir/summaries/dnadiff_merged_sup.tsv"

# Create the header row
echo -e "sample\tassembler\tTotalSNPs\tInsertions\tTotalIndels\tUnalignedBases" > "$output_file"

# Loop through all .report files in subdirectories of dnadiff
for report_file in $input_path; do
    # Extract sample name from the file basename
    sample=$(basename "$report_file" .report)
    # Extract assembler name from the directory name
    assembler=$(dirname "$report_file" | xargs basename)
    # Extract values from the report file
    total_snps=$(grep "^TotalSNPs" "$report_file" | awk '{print $3}')
    insertions=$(grep "^Insertions" "$report_file" | awk '{print $3}')
    total_indels=$(grep "^TotalIndels" "$report_file" | awk '{print $3}')
    unaligned_bases=$(grep "^UnalignedBases" "$report_file" | awk '{print $3}')
    # Append the values to the output file
    echo -e "$sample\t$assembler\t$total_snps\t$insertions\t$total_indels\t$unaligned_bases" >> "$output_file"
done

echo "Merged report files into $output_file"  
