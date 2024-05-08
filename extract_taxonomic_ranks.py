import sys

def read_rankedlineage_file(filename):
    """
    Read the rankedlineage file and store taxonomic information for each taxid.
    
    Parameters:
        filename (str): The name of the rankedlineage file.
        
    Returns:
        dict: A dictionary mapping taxids to their taxonomic lineage.
    """
    data = {}
    with open(filename, "r") as file:
        for line in file:
            fields = line.strip().split("|")
            taxid = fields[0].strip()
            lineage = '|'.join(field.strip() for field in fields[1:])  # Join all fields after taxid with "|"
            data[taxid] = lineage
    return data

def main():
    """
    Main function to read taxids from input file and generate taxonomic lineage output.
    """
    if len(sys.argv) < 3 or sys.argv[1] in {"-h", "--help"}:
        print("Usage: python script.py <input_file> <output_file>")
        print("Reads a list of taxids from input_file and generates taxonomic lineage information in output_file.")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]

    rankedlineage_data = read_rankedlineage_file("rankedlineage.dmp")

    header = "TaxID\tTax_name\tSpecies\tGenus\tFamily\tOrder\tClass\tPhylum\tKingdom\tSuperkingdom\n"

    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        outfile.write(header)
        for line in infile:
            taxid = line.strip()
            lineage = rankedlineage_data.get(taxid, "NA")  # Default to "NA" if taxid not found
            lineage = lineage.rstrip("|")  # Remove the last "|"
            lineage = lineage.replace("|", "\t")  # Replace remaining "|" with tabs
            outfile.write(f"{taxid}\t{lineage}\n")

if __name__ == "__main__":
    main()
