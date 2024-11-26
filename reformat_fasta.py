import sys

def reformat_fasta(input_fasta, output_fasta, line_length=80):
    """
    Reformats a FASTA file so that sequences are written with a specified line length.

    Parameters:
        input_fasta (str): Path to the input FASTA file.
        output_fasta (str): Path to the reformatted output FASTA file.
        line_length (int): Maximum length of each line in the sequence. Default is 80 characters.
    """
    with open(input_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
        sequence = ''
        header = ''
        for line in infile:
            line = line.strip()
            if line.startswith('>'):
                if sequence:
                    # Write the sequence in chunks of size `line_length`
                    for i in range(0, len(sequence), line_length):
                        outfile.write(sequence[i:i+line_length] + '\n')
                header = line
                outfile.write(header + '\n')  # Write the header
                sequence = ''
            else:
                # Concatenate the sequence
                sequence += line
        # Write the last sequence if it exists
        if sequence:
            for i in range(0, len(sequence), line_length):
                outfile.write(sequence[i:i+line_length] + '\n')

if __name__ == "__main__":
    # Verify if the correct number of arguments is provided
    if len(sys.argv) != 3:
        print("Usage: python reformat_fasta.py <input_fasta> <output_fasta>")
        sys.exit(1)

    # Get the file names from the command line arguments
    input_fasta = sys.argv[1]
    output_fasta = sys.argv[2]

    # Call the function with the provided files
    reformat_fasta(input_fasta, output_fasta)
