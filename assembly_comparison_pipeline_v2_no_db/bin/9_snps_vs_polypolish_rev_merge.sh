#!/bin/bash

#Define input file
input_path="/gpfs3/well/bag/users/nch313/GNBSI_pilot4/snps_vs_polypolish_rev/*/*.snps"
# Define the output file
output_file="/gpfs3/well/bag/users/nch313/GNBSI_pilot4/bin/summaries/snps_vs_polypolish_rev_merged_sup.tsv"

# Write the header to the output file
echo -e "sample\tassembler\tposition_ref\tbase_ref\tbase_qry\tposition_qry\tbuff_ref\tbuff_qry\tdist_ref\tdist_qry\tfrm_ref\tfrm_qry\tcontig_ref\tcontig_qry" > "$output_file"

# Iterate through all .snps files in the directory structure
for snps_file in $input_path; do
    # Extract the assembler and sample from the file path
    assembler=$(basename "$(dirname "$snps_file")")
    sample=$(basename "$snps_file" .snps)

    # Read the .snps file and append its content to the output file with additional columns
    awk -v sample="$sample" -v assembler="$assembler" '
        BEGIN { FS = "\t"; OFS = "\t" }
        {
            print sample, assembler, $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12
        }
    ' "$snps_file" >> "$output_file"
done

echo "Merging completed. Output saved to $output_file"