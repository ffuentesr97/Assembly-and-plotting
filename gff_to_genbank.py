import sys
import os
from Bio import SeqIO
from BCBio import GFF

# Ensure the script is called with the correct number of arguments
if len(sys.argv) != 3:
    print("Usage: python genbank.py <fasta_file> <gff_file>")
    sys.exit(1)

# Get the input file names from the command line arguments
fasta_file = sys.argv[1]
gff_file = sys.argv[2]

# Extract the base name of the fasta file to use for the output file
base_name = os.path.splitext(os.path.basename(fasta_file))[0]
genbank_file = f"{base_name}.gbk"

# Read the sequence from the FASTA file
with open(fasta_file, "r") as fasta_handle:
    fasta_record = SeqIO.read(fasta_handle, "fasta")

# Add the molecule_type annotation
fasta_record.annotations["molecule_type"] = "DNA"

# Read the annotations from the GFF3 file
with open(gff_file, "r") as gff_handle:
    annotations = list(GFF.parse(gff_handle, base_dict={fasta_record.id: fasta_record}))

# Print the identifiers found in the annotations
print("Identifiers found in the GFF3 file:")
for annotation in annotations:
    print(annotation.id)

# Ensure there is exactly one entry in the annotations
if len(annotations) != 1:
    raise ValueError("The GFF3 file must contain annotations for exactly one replicon.")

# Get the annotated record
annotated_record = annotations[0]

# Write the annotated record in GenBank format
with open(genbank_file, "w") as genbank_handle:
    SeqIO.write(annotated_record, genbank_handle, "genbank")

print(f"The GenBank file has been successfully created: {genbank_file}")
