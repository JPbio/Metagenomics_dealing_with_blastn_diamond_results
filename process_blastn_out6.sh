#!/bin/bash

# Function to display help message
print_help() {
    echo "Usage: $0 -f <input.fasta> -t <input_blastn.tab>"
    echo
    echo "Description:"
    echo "    This script processes BLASTN results and FASTA sequences to generate various outputs."
    echo
    echo "Arguments:"
    echo "    -f, --fasta       Path to the input FASTA file."
    echo "    -t, --blastn-tab  Path to the input BLASTN tabular file."
    echo "    -h, --help        Display this help message."
    exit 1
}

# Check if any arguments are provided
if [ "$#" -eq 0 ]; then
    print_help
fi

# Parse command-line arguments
while getopts ":f:t:h" opt; do
    case ${opt} in
        f )
            fasta_file=$OPTARG
            ;;
        t )
            blastn_tab_file=$OPTARG
            ;;
        h )
            print_help
            ;;
        \? )
            echo "Invalid option: $OPTARG" 1>&2
            exit 1
            ;;
        : )
            echo "Invalid option: $OPTARG requires an argument" 1>&2
            exit 1
            ;;
    esac
done
shift $((OPTIND -1))

# Check if both fasta and blastn files are provided
if [ -z "$fasta_file" ] || [ -z "$blastn_tab_file" ]; then
    echo "Both input FASTA and BLASTN files are required." 1>&2
    print_help
fi

# Filter all viral hits
grep -i "virus" "$blastn_tab_file" > viral_input_blastn.tab

# Filter all nonviral hits
grep -v -i "virus" "$blastn_tab_file" > nonviral_input_blastn.tab

# Filter the top best hits for each contig based on the first occurrence
awk '!seen[$1]++' "$blastn_tab_file" > BestHits_input_blastn.tab

# Filter the top viral hits
grep -i "virus" BestHits_input_blastn.tab > viral_BestHits_input_blastn.tab

# Filter the top nonviral hits
grep -v -i "virus" BestHits_input_blastn.tab > nonviral_BestHits_input_blastn.tab

# Add header to the tabular blastn results
header="qseqid\tqlen\tqstart\tqend\tlength\tsseqid\tslen\tsstart\tsend\tpident\tqcovs\tevalue\tbitscore\tframes\tsframe\tsstrand\tstitle\tstaxid\tssciname\tsblastname\tsskingdom"
for file in *blastn*.tab; do
    # Prepend the header to the file
    echo -e "$header" | cat - "$file" > temp && mv temp "$file"   
done < <(yes "yes")

# Hit-based strand correction

# Get the contigs to be inverted
grep -A1 -f <(awk -F'\t' '$15 < 0' BestHits_input_blastn.tab | cut -f1)  "$fasta_file" | grep -v ^-- > Neg_strand_input.fasta

# Filter the contigs that are already strand correct
grep -A1 -f <(comm -13 <(awk -F'\t' '$15 < 0' BestHits_input_blastn.tab | cut -f1 | sort) <(cut -f1 BestHits_input_blastn.tab | sort) | grep -v ^"qseqid") "$fasta_file" | grep -v ^-- > Pos_strand_input.fasta

# Get reverse-complement of negative strand hits
fastx_reverse_complement -i Neg_strand_input.fasta -o revcomp_Neg_strand_input.fasta
rm -f Neg_strand_input.fasta

# Merge the corrected contigs for the final fasta of blast hits
cat Pos_strand_input.fasta revcomp_Neg_strand_input.fasta > contigs_gt200_correct_strand.fasta
rm -f Pos_strand_input.fasta
rm -f revcomp_Neg_strand_input.fasta

# Filter the unaligned sequences for the next steps
grep -v -f <(cut -f1 BestHits_input_blastn.tab | sort | uniq) <(grep ">" "$fasta_file" | cut -f1 -d " ") > aux_unaligned_blastn.txt
grep -A1 -f aux_unaligned_blastn.txt "$fasta_file" | grep -v ^-- > contigs_unaligned_blastn.fasta
rm -f aux_unaligned_blastn.txt

# Add information of the best hits to the headers of viral and nonvial blasnt contigs
python modify_fasta_headers_blastn.py viral_BestHits_input_blastn.tab contigs_gt200_correct_strand.fasta contigs_viral_BestHits_blastn.fasta
python modify_fasta_headers_blastn.py nonviral_BestHits_input_blastn.tab contigs_gt200_correct_strand.fasta contigs_nonviral_BestHits_blastn.fasta
