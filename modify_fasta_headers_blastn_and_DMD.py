import sys

def print_help():
    # Print help message
    print("""
    Usage: python modify_fasta_headers.py <blast_results_file> <fasta_file> <output_file> <append_so>

    Description:
        This script modifies the headers of a FASTA file based on tabular BLAST/DIAMOND results.

    Arguments:
        <blast_results_file>: Path to the tabular BLAST/DIAMOND results file.
        <fasta_file>: Path to the FASTA file.
        <output_file>: Path to the output file where modified FASTA sequences will be saved.
        <append_so>: String to append in the header.
    """)

def modify_fasta_headers(blast_results_file, fasta_file, output_file, append_so):
    # Load the BLAST results
    blast_results = {}

    # Read the tabular BLAST results file
    with open(blast_results_file) as f:
        for line in f:
            fields = line.strip().split('\t')
            query_id = fields[0]  # Extract query ID from the first column
            query_len = fields[1]  # Extract query length
            column_17 = fields[16]  # Extract column 17 (adjust index to account for zero-based indexing) that is the subject description
            if query_id not in blast_results:
                blast_results[query_id] = (query_len, column_17)  # Store subject ID and column 17 in a dictionary

    # Modify the FASTA headers based on BLAST results
    with open(fasta_file) as f:
        with open(output_file, "w") as out_f:
            skip_sequence = False  # Flag to skip sequence
            for line in f:
                if line.startswith(">"):  # If line is a header line
                    header = line.strip()[1:]  # Remove ">" from the header
                    header_parts = header.split()
                    contig_id = header_parts[0]  # Extract contig ID from the header
                    if contig_id in blast_results:  # If contig ID is present in the BLAST results
                        skip_sequence = False  # Reset flag to skip sequence
                        # Construct new header with contig ID, subject ID, and column 17 separated by tabs
                        new_header = ">{}\t{}\t{}\t{}\n".format(append_so, contig_id, blast_results[contig_id][0], blast_results[contig_id][1])
                        out_f.write(new_header)  # Write the new header to the output file
                    else:
                        skip_sequence = True  # Set flag to skip sequence
                else:
                    if not skip_sequence:  # If flag to skip sequence is not set
                        out_f.write(line)  # Write sequence lines as they are

    print("Header modification completed!")  # Print completion message

if __name__ == "__main__":
    if len(sys.argv) != 5 or sys.argv[1] in {"-h", "--help"}:  # If incorrect number of arguments or help option specified
        print_help()  # Print help message
        sys.exit(1)  # Exit the script
    
    # Extract command-line arguments
    blast_results_file = sys.argv[1]
    fasta_file = sys.argv[2]
    output_file = sys.argv[3]
    append_so = sys.argv[4]

    # Call the function to modify FASTA headers based on BLAST results
    modify_fasta_headers(blast_results_file, fasta_file, output_file, append_so)
