#!/bin/bash

# Define the parent directory containing the assembler directories as the current directory containin
current_dir=$(pwd)
# Remove '/bin' from the end of the directory path
current_prefix=${current_dir%/bin}
parent_dir="$current_prefix/assemblies"

# Move Hybracter assemblies
# Define the base source and target directories
BASE_SRC_DIR="$parent_dir/hybracter_hybrid"
BASE_TGT_DIR="$parent_dir/hybracter_hybrid"

# Loop through each sample directory
for sample_dir in "$BASE_SRC_DIR"/*; do
    if [[ -d "$sample_dir" ]]; then
        sample=$(basename "$sample_dir")

        # Handle both complete and incomplete directories
        for completeness in complete incomplete; do
            src_file="$BASE_SRC_DIR/${sample}/hybracter_hybrid_out/FINAL_OUTPUT/${completeness}/${sample}_final.fasta"
            gfa_file="$BASE_SRC_DIR/${sample}/hybracter_hybrid_out/processing/assemblies/${sample}/assembly_graph.gfa"
            plasmid_file="$BASE_SRC_DIR/${sample}/hybracter_hybrid_out/processing/plassembler/${sample}/plassembler_plasmids.gfa"
            if [[ -f $src_file ]]; then
                tgt_file="$BASE_TGT_DIR/${sample}.fasta"
                tgt_gfa_file="$BASE_TGT_DIR/${sample}.gfa"
                mkdir -p "$BASE_TGT_DIR/plasmids"
                tgt_plasmid_file="$BASE_TGT_DIR/plasmids/${sample}_plasmids.gfa"
                # Copy and rename the file
                cp "$src_file" "$tgt_file"
                echo "Copied and renamed $src_file to $tgt_file"
                cp "$gfa_file" "$tgt_gfa_file"
                echo "Copied and renamed $gfa_file to $tgt_gfa_file"
                cp "$plasmid_file" "$tgt_plasmid_file"
                echo "Copied and renamed $plasmid_file to $tgt_plasmid_file"
                break
            fi
        done

        if [[ ! -f "$BASE_TGT_DIR/${sample}.fasta" ]]; then
            echo "File not found for sample: $sample"
        fi
    fi
done

BASE_SRC_DIR="$parent_dir/hybracter_long"
BASE_TGT_DIR="$parent_dir/hybracter_long"


# Loop through each sample directory
for number in 01 02 03 04; do 
    BASE_TGT_DIR="$parent_dir/hybracter_long_$number"
    mkdir -p "$BASE_TGT_DIR"
    for sample_dir in "$BASE_SRC_DIR"/*; do
        if [[ -d "$sample_dir" ]]; then
            sample=$(basename "$sample_dir")

            # Handle both complete and incomplete directories
            for completeness in complete incomplete; do
                src_file="$BASE_SRC_DIR/${sample}/${number}/hybracter_long_out/FINAL_OUTPUT/${completeness}/${sample}_${number}_final.fasta"
                echo "Checking: $src_file"
                gfa_file="$BASE_SRC_DIR/${sample}/${number}/hybracter_long_out/processing/assemblies/${sample}_${number}/assembly_graph.gfa"
                echo "Checking: $gfa_file"
                plasmid_file="$BASE_SRC_DIR/${sample}/${number}/hybracter_long_out/processing/plassembler/${sample}_${number}/plassembler_plasmids.gfa"
                echo "Checking: $plasmid_file"
                if [[ -f $src_file ]]; then
                    tgt_file="$BASE_TGT_DIR/${sample}.fasta"
                    tgt_gfa_file="$BASE_TGT_DIR/${sample}.gfa"
                    mkdir -p "$BASE_TGT_DIR/plasmids"
                    tgt_plasmid_file="$BASE_TGT_DIR/plasmids/${sample}_plasmids.gfa"
                    # Copy and rename the file
                    cp "$src_file" "$tgt_file"
                    echo "Copied and renamed $src_file to $tgt_file"
                    cp "$gfa_file" "$tgt_gfa_file"
                    echo "Copied and renamed $gfa_file to $tgt_gfa_file"                    
                    cp "$plasmid_file" "$tgt_plasmid_file"
                    echo "Copied and renamed $plasmid_file to $tgt_plasmid_file"
                    break
                fi
            done

            if [[ ! -f "$BASE_TGT_DIR/${sample}.fasta" ]]; then
                echo "File not found for sample: $sample"
            fi
        fi
    done
done
