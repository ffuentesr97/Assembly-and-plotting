from Bio import SeqIO

# List of replicons to remove
replicons_to_remove = ['chr', 'plasmid_a1', 'plasmid_a2', 'plasmid_b', 'plasmid_c', 'plasmid_e']

# Load the input FASTA file
input_fasta = 'HH103_new.fasta'
output_fasta = 'HH103_new_filtered.fasta'

# Filter the sequences
with open(output_fasta, 'w') as output_handle:
    for record in SeqIO.parse(input_fasta, 'fasta'):
        if record.id not in replicons_to_remove:
            SeqIO.write(record, output_handle, 'fasta')

print(f'The file {output_fasta} has been created without the specified replicons.')
