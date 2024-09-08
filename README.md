# Here I'll provide scripts and explanations about how I process tabular blastn and DIAMOND results for my viral metagenomic analysis.I'll show command lines and scripts to download and format databases, filter the results, add alignments hit info to fasta headers and add taxonony information to the results;

#Formatting the nr database to get taxonomy info from DIAMOND hits

  1 Download the nr database
  
    wget https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
    pigz -d -p <number of threads> nr_cluster.fasta.gz #you can just use gunzip. but without paralelazing it will take a while
   
  2 Download the Taxonomy data from NCBI (you'll need both)

    wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip
    unzip taxdmp.zip

    wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz
    gunzip prot.accession2taxid.FULL.gz

  3 Download the Taxonomy database from which we will extract the Ranks informatios

    wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
    tar -xzvf new_taxdump.tar.gz
      
  4 Download DIAMOND executable file (check for more recent versions on their github before)

    wget http://github.com/bbuchfink/diamond/releases/download/v2.1.9/diamond-linux64.tar.gz
    tar -xzvf diamond-linux64.tar.gz

  5 Formatting the database using the Taxonomy DB files

    ./diamond makedb --in nr.gz --taxonnodes nodes.dmp --taxonnames names.dmp --taxonmap prot.accession2taxid.FULL -d nr_tax

  6- Running DIAMOND against tth entire NR indexed with taxonomy info and generating a tabular result
    ./diamond blastx -q teste.fasta -d nr_tax.dmnd -k 10 -p 31 -e 0.001 -f 6 qseqid qlen qstart qend qcovhsp length sseqid slen sstart send scovhsp pident evalue bitscore qstrand qframe stitle staxids -c 1 -b 20 --very-sensitive -o teste_out --al dmnd_hits.fasta --alfmt fasta 2> dmnd_log
    






