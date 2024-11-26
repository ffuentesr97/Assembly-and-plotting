from Bio import SeqIO

# Name of the input FASTA file
input_fasta = "HH103_new_definitivo.fasta"

# Read the input FASTA file and split into replicons
replicons = SeqIO.parse(input_fasta, "fasta")

# Counter to name the output files
count = 1

# Iterate over each replicon and write it to a separate FASTA file
for replicon in replicons:
    output_fasta = f"replicon_{count}.fasta"
    with open(output_fasta, "w") as output_handle:
        SeqIO.write(replicon, output_handle, "fasta")
    count += 1

print(f"{count-1} individual FASTA files have been created.")
