#!/bin/bash

# Define the parent directory containing the assembler directories as the current directory containing the bin
current_dir=$(pwd)

# Remove '/bin' from the end of the directory path
current_prefix=${current_dir%/bin}
parent_dir="$current_prefix/pypolca"

output_file=$current_dir/summaries/pypolca_reports_merge.tsv


# Output TSV file
echo -e "assembler\tsample\tsubstitution_errors\tindel_errors\tassembly_size\tconsensus_quality_before_polishing\tconsensus_qv_before_polishing" > $output_file

# Find all pypolca.report files matching the pattern
find $parent_dir -type f -path "${parent_dir}/*/*/pypolca_out/pypolca.report" | while read report; do
    # Extract assembler and sample names from the path
    assembler=$(echo "$report" | awk -F'/' '{print $7}')
    sample=$(echo "$report" | awk -F'/' '{print $8}')
    
    # Extract relevant stats from the report file
    substitution_errors=$(grep "Substitution Errors Found" "$report" | awk '{print $NF}')
    indel_errors=$(grep "Insertion/Deletion Errors Found" "$report" | awk '{print $NF}')
    assembly_size=$(grep "Assembly Size" "$report" | awk '{print $NF}')
    consensus_quality=$(grep "Consensus Quality Before Polishing" "$report" | awk '{print $NF}')
    consensus_qv=$(grep "Consensus QV Before Polishing" "$report" | awk '{print $NF}')
    
    # Append to the summary file
    echo -e "$assembler\t$sample\t$substitution_errors\t$indel_errors\t$assembly_size\t$consensus_quality\t$consensus_qv" >> $output_file
done

