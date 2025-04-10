#!/bin/bash

# Define the parent directory containing the assembler directories as the current directory containin
current_dir=$(pwd)
# Remove '/bin' from the end of the directory path
current_prefix=${current_dir%/bin}
parent_dir="$current_prefix/assemblies"

# Move Hybracter assemblies
# Define the base source and target directories
BASE_SRC_DIR="$parent_dir/autocycler_consensus"
BASE_TGT_DIR="$parent_dir/autocycler_consensus"

# Loop through each sample directory
for sample_dir in "$BASE_SRC_DIR"/*; do
    if [[ -d "$sample_dir" ]]; then
        sample=$(basename "$sample_dir")

        # Handle both complete and incomplete directories
        for completeness in complete incomplete; do
            src_file="$BASE_SRC_DIR/${sample}/autocycler_out/consensus_assembly.fasta"
            gfa_file="$BASE_SRC_DIR/${sample}/autocycler_out/consensus_assembly.gfa"
            if [[ -f $src_file ]]; then
                tgt_file="$BASE_TGT_DIR/${sample}.fasta"
                tgt_gfa_file="$BASE_TGT_DIR/${sample}.gfa"
                # Copy and rename the file
                cp "$src_file" "$tgt_file"
                echo "Copied and renamed $src_file to $tgt_file"
                cp "$gfa_file" "$tgt_gfa_file"
                echo "Copied and renamed $gfa_file to $tgt_gfa_file"
                break
            fi
        done

        if [[ ! -f "$BASE_TGT_DIR/${sample}.fasta" ]]; then
            echo "File not found for sample: $sample"
        fi
    fi
done


# Define the base source and target directories
BASE_SRC_DIR="$parent_dir/flye_pypolca"
BASE_TGT_DIR="$parent_dir/flye_pypolca"

# Loop through each sample directory
for sample_dir in "$BASE_SRC_DIR"/*; do
    if [[ -d "$sample_dir" ]]; then
        sample=$(basename "$sample_dir")

        # Handle both complete and incomplete directories
        for completeness in complete incomplete; do
            src_file="$BASE_SRC_DIR/${sample}/pypolca_out/pypolca_corrected.fasta"
            if [[ -f $src_file ]]; then
                tgt_file="$BASE_TGT_DIR/${sample}.fasta"

                # Copy and rename the file
                cp "$src_file" "$tgt_file"
                echo "Copied and renamed $src_file to $tgt_file"
                break
            fi
        done

        if [[ ! -f "$BASE_TGT_DIR/${sample}.fasta" ]]; then
            echo "File not found for sample: $sample"
        fi
    fi
done


# Define the base source and target directories
BASE_SRC_DIR="$parent_dir/autocycler_consensus_pypolca"
BASE_TGT_DIR="$parent_dir/autocycler_consensus_pypolca"

# Loop through each sample directory
for sample_dir in "$BASE_SRC_DIR"/*; do
    if [[ -d "$sample_dir" ]]; then
        sample=$(basename "$sample_dir")

        # Handle both complete and incomplete directories
        for completeness in complete incomplete; do
            src_file="$BASE_SRC_DIR/${sample}/pypolca_out/pypolca_corrected.fasta"
            if [[ -f $src_file ]]; then
                tgt_file="$BASE_TGT_DIR/${sample}.fasta"

                # Copy and rename the file
                cp "$src_file" "$tgt_file"
                echo "Copied and renamed $src_file to $tgt_file"
                break
            fi
        done

        if [[ ! -f "$BASE_TGT_DIR/${sample}.fasta" ]]; then
            echo "File not found for sample: $sample"
        fi
    fi
done
