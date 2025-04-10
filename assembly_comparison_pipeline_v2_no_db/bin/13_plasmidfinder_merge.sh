#!/bin/bash


current_dir=$(pwd)
# Remove '/bin' from the end of the directory path
current_prefix=${current_dir%/bin}

# Output TSV file
output_tsv="$current_dir/summaries/13_plasmidfinder_summary.tsv"

echo -e "Assembler\tSample\tContig_name\tPlasmid\tAccession\tCoverage\tIdentity\tPosition_in_ref\tPosition_in_contig\tTemplate_length" > "$output_tsv"

# Base directory containing the plasmidfinder output
base_dir="$current_prefix/PlasmidFinder"


# Loop through all assembler subdirectories
for assembler_dir in "$base_dir"/*; do
  # Ensure it's a directory
  [ -d "$assembler_dir" ] || continue
    assembler=$(basename "$assembler_dir" | sed 's/^plasmidfinder_output_//')

  # Loop through all {sample}_plasmidfinder.txt files in the subdirectory
  for txt_file in "$assembler_dir"/*_plasmidfinder.txt; do
    # Ensure the file exists
    [ -f "$txt_file" ] || continue

    # Extract the sample name (filename without "_plasmidfinder.txt")
    sample=$(basename "$txt_file" "_plasmidfinder.txt")

      grep -E "'contig_name':|'plasmid':|'accession':|'coverage':|'identity':|'position_in_ref':|'positions_in_contig':|'template_length':" "$txt_file" | awk -v assembler="$assembler" -v sample="$sample" '
      BEGIN {
        contig_count = 0
        plasmid_count = 0
        accession_count = 0
        coverage_count = 0
        identity_count = 0
        position_ref_count = 0
        position_contig_count = 0
        template_length_count = 0
      }
      /'"'"'contig_name'"'"'/ {
        gsub(/'"'"'/, "", $2) # Remove single quotes
        gsub(/,$/, "", $2)    # Remove the comma at the end
        contig_name[contig_count++] = $2
      }
      /'"'"'plasmid'"'"'/ {
        gsub(/'"'"'/, "", $2) # Remove single quotes
        gsub(/,$/, "", $2)    # Remove the comma at the end
        plasmid_name[plasmid_count++] = $2
      }
      /'"'"'accession'"'"'/ {
        gsub(/'"'"'/, "", $2) # Remove single quotes
        gsub(/,$/, "", $2)    # Remove the comma at the end
        accession[accession_count++] = $2
      }
      /'"'"'coverage'"'"'/ {
        gsub(/,$/, "", $2)    # Remove the comma at the end (no single quotes)
        coverage[coverage_count++] = $2
      }
      /'"'"'identity'"'"'/ {
        gsub(/,$/, "", $2)    # Remove the comma at the end (no single quotes)
        identity[identity_count++] = $2
      }
      /'"'"'position_in_ref'"'"'/ {
        gsub(/'"'"'/, "", $2) # Remove single quotes
        gsub(/,$/, "", $2)    # Remove the comma at the end
        position_ref[position_ref_count++] = $2
      }
      /'"'"'positions_in_contig'"'"'/ {
        gsub(/'"'"'/, "", $2) # Remove single quotes
        gsub(/,$/, "", $2)    # Remove the comma at the end
        position_contig[position_contig_count++] = $2
      }
      /'"'"'template_length'"'"'/ {
        gsub(/'"'"'/, "", $2) # Remove single quotes
        gsub(/,$/, "", $2)    # Remove the comma at the end
        gsub(/}+$/, "", $2)   # Remove trailing } (1 or more)
        template_length[template_length_count++] = $2
      }
      END {
        # Determine the maximum number of findings
        max_count = contig_count
        if (plasmid_count > max_count) max_count = plasmid_count
        if (accession_count > max_count) max_count = accession_count
        if (coverage_count > max_count) max_count = coverage_count
        if (identity_count > max_count) max_count = identity_count
        if (position_ref_count > max_count) max_count = position_ref_count
        if (position_contig_count > max_count) max_count = position_contig_count
        if (template_length_count > max_count) max_count = template_length_count

        # Print each finding to separate lines
        for (i = 0; i < max_count; i++) {
          printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
                 assembler,
                 sample,
                 contig_name[i],
                 plasmid_name[i],
                 accession[i],
                 coverage[i],
                 identity[i],
                 position_ref[i],
                 position_contig[i],
                 template_length[i])
        }
      }
    ' >> "$output_tsv"
  done
done

echo "Summary saved to $output_tsv"
