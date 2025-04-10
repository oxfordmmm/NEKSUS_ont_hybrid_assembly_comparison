#!/usr/bin/env python3
import os
import re
import pandas as pd
from Bio import SeqIO


# Get current directory and adjust paths
current_dir = os.getcwd()
if current_dir.endswith("/bin"):
    current_dir = current_dir[:-4]

# Set the assembly directory
assembly_dir = os.path.join(current_dir, "assemblies")

# Output file path
output_file = os.path.join(current_dir, "summaries", "contigs_summary.tsv")



# Open output file for writing
	
	
with open(output_file, "w") as out_file:
 # Write header
	out_file.write("sample\tassembler\tcontig_name\tlength\tcircular\tdepth\n")
	    
	    # Iterate over directories in assembly_dir
	for subdir in os.listdir(assembly_dir):
	    sub_dir_path = os.path.join(assembly_dir, subdir)
	    if not os.path.isdir(sub_dir_path):
	        continue
	        
	        # Extract assembler name from the subdirectory name
	    assembler = subdir
	        
	        # Iterate over files in the current subdirectory
	    for file in os.listdir(sub_dir_path):
	            # Check if the file is a FASTA file
	        if file.endswith(".fasta"):
	                # Extract sample name from the file name
	            sample = os.path.splitext(file)[0]
	                
	                # Parse the FASTA file
	            fasta_file = os.path.join(sub_dir_path, file)
	            for record in SeqIO.parse(fasta_file, "fasta"):
	                    # Parse the FASTA header to extract circularity information and depth
	                header = record.description
	                circular_match = re.search(r"circular=([A-Za-z]+)", header)
                        suggest_circular_match = re.search(r"suggestCircular=(yes|no)", header)

                        if circular_match:
                        circular = circular_match.group(1)
                    elif suggest_circular_match:
                        circular = "true" if suggest_circular_match.group(1) == "yes" else "false"
                    else:
                        circular = "NA"

	                # Extract depth from `depth=` or `plasmid_copy_number_long=`
                    depth_match = re.search(r"depth=(\d+)", header)
                    plasmid_copy_match = re.search(r"plasmid_copy_number_long=(\d+\.\d+)x?", header)

                    if depth_match:
                        depth = depth_match.group(1)
                    elif plasmid_copy_match:
                        depth = plasmid_copy_match.group(1)
                    else:
                        depth = "NA"

	                    # Write contig summary to output file
	                out_file.write(f"{sample}\t{assembler}\t{record.id}\t{len(record.seq)}\t{circular}\t{depth}\n")
					#Or modify to the following if using Python version <3.6 (i.e. default on some older computer clusters) as F-string introduced >=v3.6
					#out_file.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(sample, assembler, record.id, len(record.seq), circular, depth))

#Add Flye contig metadata from .txt. files

df = pd.read_csv(output_file, sep='\t')

# Ensure the new column 'copy_number' exists
if 'copy_number' not in df.columns:
    df['copy_number'] = 'NA'

# Function to process each .txt file
def process_txt_file(file_path, sample_name):
    # Read the .txt file into a DataFrame
    txt_df = pd.read_csv(file_path, sep='\t')

    # Iterate through each row of the .txt DataFrame
    for _, row in txt_df.iterrows():
        # Extract necessary values
        contig_name = row['#seq_name']
        depth = row['cov.']
        circular = row['circ.'] == 'Y'  # Convert 'Y' to boolean True
        copy_number = row.get('mult.', 'NA')

        # Update the corresponding row in the output DataFrame
        mask = (df['sample'] == sample_name) & (df['contig_name'] == contig_name) & (df['assembler'] == "FLYE")
        df.loc[mask, 'depth'] = depth
        df.loc[mask, 'circular'] = circular
        df.loc[mask, 'copy_number'] = copy_number

# Loop through flye_01 to flye_04 directories
for i in range(1, 5):
    input_directory = os.path.join(assembly_dir, f"flye_0{i}")

    if os.path.isdir(input_directory):
        for file_name in os.listdir(input_directory):
            if file_name.endswith('.txt'):
                file_path = os.path.join(input_directory, file_name)
                sample_name = file_name.replace('.txt', '')
                process_txt_file(file_path, sample_name)
    else:
        print(f"Warning: {input_directory} does not exist.")

# Save the updated DataFrame back to the TSV file
df.to_csv(output_file, sep='\t', index=False)

print(f"Updated contigs summary saved to {output_file}")
