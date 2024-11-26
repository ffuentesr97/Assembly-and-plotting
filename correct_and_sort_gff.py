def correct_and_sort_gff3(input_file, output_file):
    try:
        # Read the GFF3 file and store the lines in a list
        with open(input_file, 'r') as f:
            lines = f.readlines()

        # Process each line to correct start and end positions if necessary
        corrected_lines = []
        for line in lines:
            fields = line.strip().split('\t')
            start = int(fields[3])
            end = int(fields[4])

            # Check and correct if the start position is greater than the end position
            if start > end:
                # Swap start and end
                fields[3], fields[4] = str(end), str(start)
                corrected_line = '\t'.join(fields) + '\n'
                corrected_lines.append(corrected_line)
            else:
                corrected_lines.append(line)

        # Sort the corrected GFF3 lines by replicon and start position
        sorted_lines = sorted(corrected_lines, key=lambda x: (get_replicon(x), get_start(x)))

        # Write the sorted lines to the output file
        with open(output_file, 'w') as f:
            f.writelines(sorted_lines)

    except Exception as e:
        print(f"Error processing the file: {str(e)}")

def get_replicon(line):
    try:
        # Get the replicon from a GFF3 file line
        fields = line.strip().split('\t')
        replicon = fields[0]
        return replicon
    except IndexError:
        return ""

def get_start(line):
    try:
        # Get the start position of an annotation in a GFF3 file line
        fields = line.strip().split('\t')
        start = int(fields[3])
        return start
    except (IndexError, ValueError):
        return float('inf')  # Return infinite value to handle errors

# Example usage
input_file = "HH103_new_def_def_def.gff"
output_file = "HH103_new_def_def_def_sorted.gff"

correct_and_sort_gff3(input_file, output_file)
