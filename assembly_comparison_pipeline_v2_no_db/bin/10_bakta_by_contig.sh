#!/bin/bash

current_dir=$(pwd)
# Remove '/bin' from the end of the directory path
current_prefix=${current_dir%/bin}

# Define the output file
output_file="$current_dir/summaries/bakta_by_contig.tsv"

# Print the header to the output file
echo -e "assembler\tsample\tcontig_no\tCDS\trRNA\ttRNA\ttmRNA\tncRNA\tncRNA-region\tsorf\toriC\toriT\toriV\tpseudogene\thypothetical\tsignal_peptide\tCRISPR_array" > "$output_file"


# Loop through each .tsv file in the directory structure
for file in $current_prefix/bakta/*/*_bakta/*.tsv; do
    assembler=$(basename $(dirname $(dirname "$file")))
    sample=$(basename "$file" _bakta/${sample}.tsv | sed 's/_bakta//')
    echo "Processing assembler: $assembler, sample: $sample..."

    # Read the file and extract the information
    awk -v assembler="$assembler" -v sample="$sample" -v OFS="\t" '
    BEGIN { FS=OFS="\t" }
    NR > 6 {
        contig_no[$1]++;
        if ($2 == "cds") { CDS[$1]++ }
        else if ($2 == "rRNA") { rRNA[$1]++ }
        else if ($2 == "tRNA") { tRNA[$1]++ }
        else if ($2 == "tmRNA") { tmRNA[$1]++ }
        else if ($2 == "ncRNA") { ncRNA[$1]++ }
        else if ($2 == "ncRNA-region") { ncRNA_region[$1]++ }
        else if ($2 == "sorf") { sorf[$1]++ }
        else if ($2 == "oriC") { oriC[$1]++ }
        else if ($2 == "oriT") { oriT[$1]++ }
        else if ($2 == "oriV") { oriV[$1]++ }

        
        if ($8 ~ /^pseudogene/) { pseudogene[$1]++ }
        if ($8 ~ /^hypothetical/) { hypothetical[$1]++ }
        if ($8 ~ /^signal peptide/) { signal_peptide[$1]++ }
        if ($8 ~ /^CRISPR array/) { CRISPR_array[$1]++ }
    }
    END {
        for (contig in contig_no) {
            print assembler, sample, contig, 
                (CDS[contig] ? CDS[contig] : 0), 
                (rRNA[contig] ? rRNA[contig] : 0), 
                (tRNA[contig] ? tRNA[contig] : 0), 
                (tmRNA[contig] ? tmRNA[contig] : 0), 
                (ncRNA[contig] ? ncRNA[contig] : 0), 
                (ncRNA_region[contig] ? ncRNA_region[contig] : 0),
                (sorf[contig] ? sorf[contig] : 0), 
                (oriC[contig] ? oriC[contig] : 0), 
                (oriT[contig] ? oriT[contig] : 0), 
                (oriV[contig] ? oriV[contig] : 0),
                (pseudogene[contig] ? pseudogene[contig] : 0),
                (hypothetical[contig] ? hypothetical[contig] : 0),
                (signal_peptide[contig] ? signal_peptide[contig] : 0),
                (CRISPR_array[contig] ? CRISPR_array[contig] : 0);
        }
    }
    ' "$file" >> "$output_file"
done

echo "Summary file created: $output_file"

#Now clean up the output file:

# Define the input and output files
input_csv="$current_dir/summaries/bakta_by_contig.tsv"
output_csv="$current_dir/summaries/cleaned_bakta_by_contig.tsv"

# Print the header to the output file
head -n 1 "$input_csv" > "$output_csv"

# Process the file, skipping rows with 'hypothetical' in the second column
# and removing '.tsv' from the end of sample names in the second column
awk -F"\t" -v OFS="\t" '
NR > 1 {
    if ($2 !~ /hypothetical/) {
        gsub(/\.tsv$/, "", $2)
        print $0
    }
}
' "$input_csv" >> "$output_csv"

echo "Cleaned summary file created: $output_csv"
