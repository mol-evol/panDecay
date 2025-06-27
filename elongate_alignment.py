#!/usr/bin/env python3
"""
This script takes a NEXUS alignment file and elongates it by repeating
each sequence 10 times. The original alignment length is preserved in each
repeat, resulting in a total length 10 times the original.

Usage:
    python elongate_alignment.py input_file.nex output_file.nex
"""

import sys
import re

def parse_nexus_file(filename):
    """Parse a NEXUS file and extract the alignment information."""
    with open(filename, 'r') as f:
        content = f.read()

    # Extract dimensions
    dimensions_match = re.search(r'dimensions\s+ntax=(\d+)\s+nchar=(\d+)', content)
    if dimensions_match:
        ntax = int(dimensions_match.group(1))
        nchar = int(dimensions_match.group(2))
    else:
        raise ValueError("Could not find dimensions in NEXUS file")

    # Check format settings
    format_match = re.search(r'format\s+([^;]+);', content)
    format_settings = format_match.group(1) if format_match else ""
    matchchar_match = re.search(r'matchchar\s*=\s*([^\s]+)', format_settings)
    matchchar = matchchar_match.group(1) if matchchar_match else None

    # Extract matrix
    matrix_match = re.search(r'matrix\s+([\s\S]+?);', content)
    if not matrix_match:
        raise ValueError("Could not find matrix in NEXUS file")

    matrix_text = matrix_match.group(1)

    # Parse sequences from matrix
    sequences = {}
    taxon_sequence_blocks = {}

    for line in matrix_text.split('\n'):
        line = line.strip()
        if line and not line.startswith('['):  # Skip comments
            if re.match(r'^[A-Za-z0-9_]+\s+', line):  # Line starts with taxon name
                parts = re.split(r'\s+', line, 1)
                if len(parts) > 1:
                    taxon = parts[0]
                    sequence = parts[1].replace(" ", "")

                    if taxon not in taxon_sequence_blocks:
                        taxon_sequence_blocks[taxon] = []

                    taxon_sequence_blocks[taxon].append(sequence)

    # Combine sequence blocks for each taxon
    for taxon, blocks in taxon_sequence_blocks.items():
        combined_sequence = ''.join(blocks)

        # If a matchchar is used, resolve it by replacing with the first sequence
        if matchchar and taxon != list(taxon_sequence_blocks.keys())[0]:
            first_taxon = list(taxon_sequence_blocks.keys())[0]
            first_sequence = ''.join(taxon_sequence_blocks[first_taxon])

            # Replace matchchar with corresponding character from first sequence
            resolved_sequence = ""
            for i, char in enumerate(combined_sequence):
                if char == matchchar:
                    if i < len(first_sequence):
                        resolved_sequence += first_sequence[i]
                    else:
                        resolved_sequence += 'N'  # Use N if outside first sequence
                else:
                    resolved_sequence += char

            sequences[taxon] = resolved_sequence
        else:
            sequences[taxon] = combined_sequence

    return {
        'content': content,
        'dimensions': {
            'ntax': ntax,
            'nchar': nchar
        },
        'matchchar': matchchar,
        'sequences': sequences
    }

def elongate_alignment(nexus_data):
    """Elongate each sequence in the alignment by repeating it 10 times."""
    elongated_sequences = {}
    for taxon, sequence in nexus_data['sequences'].items():
        elongated_sequences[taxon] = sequence * 10

    return elongated_sequences

def create_elongated_nexus(nexus_data, elongated_sequences, output_file):
    """Create a new NEXUS file with the elongated sequences."""
    with open(output_file, 'w') as f:
        # Extract and write the header
        header_match = re.search(r'(#NEXUS\s+begin\s+data;[\s\S]+?format\s+[^;]+);', nexus_data['content'])
        if header_match:
            # Modify the header to update dimensions
            header = header_match.group(1)
            header = re.sub(
                r'dimensions\s+ntax=\d+\s+nchar=\d+',
                f'dimensions ntax={nexus_data["dimensions"]["ntax"]} nchar={nexus_data["dimensions"]["nchar"]*10}',
                header
            )

            # Ensure matchchar is removed to avoid issues
            if nexus_data['matchchar']:
                header = re.sub(
                    r'matchchar=[^\s]+',
                    'matchchar=none',
                    header
                )

            f.write(header + ";\n")
        else:
            # Create a basic header
            f.write("#NEXUS\n\nbegin data;\n")
            f.write(f"\tdimensions ntax={nexus_data['dimensions']['ntax']} nchar={nexus_data['dimensions']['nchar']*10};\n")
            f.write("\tformat missing=? gap=- datatype=dna;\n")

        # Option settings
        f.write("\toptions gapmode=missing;\n")

        # Write the matrix with interleaved format
        f.write("\tmatrix\n")

        # Write sequences in blocks of 70 characters
        line_length = 70
        max_taxon_length = max(len(taxon) for taxon in elongated_sequences.keys())
        num_bases = len(next(iter(elongated_sequences.values())))

        for i in range(0, num_bases, line_length):
            for taxon in sorted(elongated_sequences.keys()):
                sequence_block = elongated_sequences[taxon][i:i+line_length]
                padding = ' ' * (max_taxon_length - len(taxon) + 4)
                f.write(f"{taxon}{padding}{sequence_block}\n")
            # Add a blank line between blocks (except after the last block)
            if i + line_length < num_bases:
                f.write("\n")

        # Close the matrix and block
        f.write("\t;\nendblock;\n")

def main():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} input_file.nex output_file.nex")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    try:
        nexus_data = parse_nexus_file(input_file)
        elongated_sequences = elongate_alignment(nexus_data)
        create_elongated_nexus(nexus_data, elongated_sequences, output_file)
        print(f"Successfully created elongated alignment: {output_file}")
        print(f"Original length: {nexus_data['dimensions']['nchar']} bp")
        print(f"New length: {nexus_data['dimensions']['nchar'] * 10} bp")
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
