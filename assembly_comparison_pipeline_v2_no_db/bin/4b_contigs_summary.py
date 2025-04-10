import os
import re
import pathlib
from Bio import SeqIO


# Get current directory and adjust paths
current_dir = os.getcwd()
if current_dir.endswith("/bin"):
    current_dir = current_dir[:-4]

# Set assembly directory & output path
assembly_dir = os.path.join(current_dir, "assemblies")
output_dir = os.path.join(current_dir, "bin", "summaries")

# Ensure the summaries directory exists
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

output_file = os.path.join(output_dir, "contigs_summary.tsv")


def parse_gfa(gfa_file):
    """
    Parse GFA file to extract contig information: length, circularity, and depth.
    """
    contigs = {}
    circular, not_circular = set(), set()

    print(f"Parsing GFA file: {gfa_file}")  # Debugging print

    with open(gfa_file, "rt") as gfa:
        for line in gfa:
            parts = line.strip().split("\t")

            # Segment (S) line: Extract contig name and sequence length
            if parts[0] == "S":
                contig_name = parts[1]
                sequence = parts[2]
                contigs[contig_name] = {
                    "length": len(sequence),
                    "circular_gfa": "false",
                    "depth_gfa": "NA",
                    "circular_fasta": "NA",
                    "depth_fasta": "NA"
                }

                # Extract depth information (searching both DP and dp)
                depth_match_float = re.search(r"DP:f:([\d.]+)", line)  # Floating-point depth
                depth_match_int = re.search(r"(?:DP|dp):i:(\d+)", line)
                copy_number_match = re.search(r"plasmid_copy_number_long=([\d.]+)x?", line)

                if depth_match_float:
                    contigs[contig_name]["depth_gfa"] = depth_match_float.group(1)
                elif depth_match_int:
                    contigs[contig_name]["depth_gfa"] = depth_match_int.group(1)
                elif copy_number_match:
                    contigs[contig_name]["depth_gfa"] = copy_number_match.group(1)


                # Link (L) line: Detect circular contigs
            elif parts[0] == "L":
                seg_1, strand_1 = parts[1], parts[2]
                seg_2, strand_2 = parts[3], parts[4]
                cigar = parts[5]

                if seg_1 == seg_2 and strand_1 == strand_2 and cigar == "0M":
                    circular.add(seg_1)
                else:
                    not_circular.add(seg_1)
                    not_circular.add(seg_2)

    # Identify truly circular contigs
    circular -= not_circular
    for contig in circular:
        if contig in contigs:
            contigs[contig]["circular_gfa"] = "true"

    return contigs


def parse_fasta(fasta_file, contigs):
    """
    Parse FASTA file to extract circularity and depth information.
    """

    print(f"Parsing FASTA file: {fasta_file}")  # Debugging print

    with open(fasta_file, "rt") as fasta:
        for record in SeqIO.parse(fasta, "fasta"):
            header = record.description
            contig_name = record.id
            seq_length = len(record.seq)

            # If the contig is already in the dictionary, update it; otherwise, add it
            if contig_name not in contigs:
                contigs[contig_name] = {
                    "length": seq_length,
                    "circular_gfa": "NA",
                    "depth_gfa": "NA",
                    "circular_fasta": "NA",
                    "depth_fasta": "NA"
                }
            else:
                contigs[contig_name]["length"] = seq_length  # Ensure length is captured

            # Extract circularity from circular= or suggestCircular=
            circular_match = re.search(r"circular=([A-Za-z]+)", header)
            suggest_circular_match = re.search(r"suggestCircular=(yes|no)", header)

            if circular_match:
                contigs[contig_name]["circular_fasta"] = circular_match.group(1)
            elif suggest_circular_match:
                contigs[contig_name]["circular_fasta"] = "true" if suggest_circular_match.group(1) == "yes" else "false"
            else:
                contigs[contig_name]["circular_fasta"] = "NA"

            # Extract depth from `depth=` or `plasmid_copy_number_long=`
            depth_match = re.search(r"depth=(\d+)", header)
            plasmid_copy_match = re.search(r"plasmid_copy_number_long=(\d+\.\d+)x?", header)

            if depth_match:
                contigs[contig_name]["depth_fasta"] = depth_match.group(1)
            elif plasmid_copy_match:
                contigs[contig_name]["depth_fasta"] = plasmid_copy_match.group(1)


def main():
    # Ensure arguments are passed correctly
    contigs_summary = []

    # Iterate through each assembler directory
    for assembler in os.listdir(assembly_dir):
        assembler_path = os.path.join(assembly_dir, assembler)
        if not os.path.isdir(assembler_path):
            continue

        # List all files in the assembler directory
        for file in os.listdir(assembler_path):
            # Check for .gfa and .fasta files in the current directory
            if file.endswith(".gfa"):
                sample = os.path.splitext(file)[0]
                gfa_file = os.path.join(assembler_path, file)

                # Parse GFA file
                contigs = {}
                print(f"Found GFA file: {gfa_file}")  # Debugging print
                contigs.update(parse_gfa(gfa_file))

                # Check for corresponding FASTA file
                fasta_file = os.path.join(assembler_path, f"{sample}.fasta")
                if os.path.isfile(fasta_file):
                    print(f"Found FASTA file: {fasta_file}")  # Debugging print
                    parse_fasta(fasta_file, contigs)

                # Collect contig information for the TSV
                for contig_name, info in contigs.items():
                    contigs_summary.append([
                        sample, assembler, contig_name, info['length'],
                        info['circular_gfa'], info['depth_gfa'], 
                        info['circular_fasta'], info['depth_fasta']
                    ])

    # Write the collected information to the output file
    if contigs_summary:
        with open(output_file, "w") as out_file:
            out_file.write("sample\tassembler\tcontig_name\tlength\tcircular_gfa\tdepth_gfa\tcircular_fasta\tdepth_fasta\n")
            for row in contigs_summary:
                out_file.write("\t".join(map(str, row)) + "\n")

        print(f"Contigs summary saved to {output_file}")
    else:
        print("No contigs data found.")


if __name__ == "__main__":
    main()
