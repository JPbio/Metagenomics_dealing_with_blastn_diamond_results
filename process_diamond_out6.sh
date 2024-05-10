#!/bin/bash

# Function to display help message
print_help() {
    echo "Usage: $0 -f <input.fasta> -t <input_diamond.tab>"
    echo
    echo "Description:"
    echo "    This script processes DIAMOND results and FASTA sequences to generate various outputs."
    echo
    echo "Arguments:"
    echo "    -f, --fasta       Path to the input FASTA file."
    echo "    -t, --diamond-tab  Path to the input DIAMOND tabular file."
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
            diamond_tab_file=$OPTARG
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

# Check if both fasta and diamond files are provided
if [ -z "$fasta_file" ] || [ -z "$diamond_tab_file" ]; then
    echo "Both input FASTA and DIAMOND files are required." 1>&2
    print_help
fi

#Get list od Taxids
cut -f 18 "$diamond_tab_file" | cut -f 1 -d ";" | sort | uniq > uniq_taxids_diamond.txt #Here I'm getting only the first Taxid separated by ";". Not sure if it could have implications in the interpretability,need think more about it and maybe deal with it.

#Get Ranks taxons information for each TaxID
python extract_taxonomic_ranks.py uniq_taxids_diamond.txt ranks_taxids_diamond.tab
rm -f uniq_taxids_diamond.txt

#keeping only one taxid per hit in the tabular result
diamond_tab_file=$(awk -F'\t' '{split($18, arr, ";"); $18 = arr[1]; for (i = 1; i <= NF; i++) printf "%s%s", $i, (i == NF ? "\n" : "\t")}' <<< "$diamond_tab_file")


#Appending the taxonomy information to the tabular blast result
awk -v diamond_tab="$diamond_tab_file" '
    BEGIN{FS=OFS="\t"}
    FNR==NR{test_out[$1]=$0; next}
    {
        if ($18 in test_out) {
            split(test_out[$18], fields, "\t")
            printf "%s%s%s\n", $0, OFS, substr(test_out[$18], length(fields[1]) + 2)
        } else {
            print $0
        }
    }
' ranks_taxids_diamond.tab "$diamond_tab_file" > diamond_gt200_AllHits_outfmt6_TAX.tab


# Filter all viral hits
awk -F'\t' '$27 == "Viruses"' diamond_gt200_AllHits_outfmt6_TAX.tab > viral_AllHits_diamond.tab

# Filter all nonviral hits
awk -F'\t' '$27 != "Viruses"' diamond_gt200_AllHits_outfmt6_TAX.tab > nonviral_AllHits_diamond.tab

# Filter the top best hits for each contig based on the first occurrence
awk '!seen[$1]++' diamond_gt200_AllHits_outfmt6_TAX.tab  > BestHits_diamond.tab

# Filter the top viral hits
awk -F'\t' '$27 == "Viruses"' BestHits_diamond.tab > viral_BestHits_diamond.tab

# Filter the top nonviral hits
awk -F'\t' '$27 != "Viruses"' BestHits_diamond.tab > nonviral_BestHits_diamond.tab

# Add header to the tabular blastn results
header="qseqid\tqlen\tqstart\tqend\tqcovhsp\tlength\tsseqid\tslen\tsstart\tsend\tscovhsp\tpident\tevalue\tbitscore\tqstrand\tqframe\tstitle\tstaxid\tTax_name\tSpecies\tGenus\tFamily\tOrder\tClass\tPhylum\tKingdom\tSuperkingdom"
for file in *diamond*.tab; do
    # Prepend the header to the file
    echo -e "$header" | cat - "$file" > temp && mv temp "$file"
done < <(yes "yes")

#Linearizing diamond fasta outputs
fasta_formatter -i diamond_aligned.fasta  -o diamond_aligned_linear.fasta
rm -f diamond_aligned.fasta 
fasta_formatter -i diamond_unaligned.fasta  -o contigs_nohits_FINAL_diamond.fasta #Final contigs without similarity: "dark matter"
rm -f diamond_unaligned.fasta



## Hit-based strand correction
# Get the contigs to be inverted
grep -A1 -f <(awk -F'\t' '$16 < 0' BestHits_diamond.tab | cut -f1)  diamond_aligned_linear.fasta | grep -v ^-- > Neg_strand_input.fasta

# Filter the contigs that are already strand correct
grep -A1 -f <(comm -13 <(awk -F'\t' '$16 < 0' BestHits_diamond.tab | cut -f1 | sort) <(cut -f1 BestHits_diamond.tab | sort) | grep -v ^"qseqid") diamond_aligned_linear.fasta | grep -v ^-- > Pos_strand_input.fasta

# Get reverse-complement of negative strand hits
fastx_reverse_complement -i Neg_strand_input.fasta -o revcomp_Neg_strand_input.fasta
rm -f Neg_strand_input.fasta

# Merge the corrected contigs for the final fasta of blast hits
cat Pos_strand_input.fasta revcomp_Neg_strand_input.fasta > contigs_gt200_diamond_correct_strand.fasta
rm -f Pos_strand_input.fasta
rm -f revcomp_Neg_strand_input.fasta

# Add information of the best hits to the headers of viral and nonvial blasnt contigs
python modify_fasta_headers_blastn_and_DMD.py viral_BestHits_diamond.tab contigs_gt200_diamond_correct_strand.fasta contigs_viral_BestHits_diamond.fasta bX_
python modify_fasta_headers_blastn_and_DMD.py nonviral_BestHits_diamond.tab contigs_gt200_diamond_correct_strand.fasta contigs_nonviral_BestHits_diamond.fasta bX_

# add the prefix "bX" to all contigs IDs
for file in *diamond*tab; do
    if [[ ! $file =~ ^rank ]]; then
        sed -i '/^qseqid/! s/^/bX_/' "$file"
    fi
done

#### END
