#!/bin/bash

# Determine the current prefix as the directory of this script without /bin at the end
current_prefix=$(dirname "$(realpath "$0")")
current_prefix=${current_prefix%/bin}

# Output file
output_file="${current_prefix}/bin/summaries/bakta_cds_summary.tsv"
echo -e "assembler\tsample\tcontig_id\tcds_no\tcds_total_length\tcds_mean_length\tcoding_density" > "$output_file"

# Iterate over all matching TSV files
for tsv_file in "$current_prefix"/bakta/*/*_bakta/*.tsv; do
    assembler=$(basename "$(dirname "$(dirname "$tsv_file")")")
    sample=$(basename "${tsv_file}" .tsv)
    sample=${sample%_bakta}  # Remove '_bakta' suffix if present
    
    # Get the coding density from the corresponding sample.txt file
    sample_txt="$(dirname "$tsv_file")/$sample.txt"
    if [[ -f "$sample_txt" ]]; then
        coding_density=$(grep -i "^coding density:" "$sample_txt" | awk '{print $3}')
    else
        coding_density="NA"
    fi

    # Process the TSV file
    awk -v assembler="$assembler" -v sample="$sample" -v coding_density="$coding_density" 'BEGIN {FS="\t"; OFS="\t"}
    !/^#/ {
        contig_id = $1;
        if ($2 == "cds") {
            cds_count[contig_id]++;
            cds_length[contig_id] += ($4 > $3) ? ($4 - $3) : ($3 - $4);
        }
    }
    END {
        for (id in cds_count) {
            mean_length = cds_count[id] > 0 ? cds_length[id] / cds_count[id] : 0;
            print assembler, sample, id, cds_count[id], cds_length[id], mean_length, coding_density;
        }
    }' "$tsv_file" >> "$output_file"
done
