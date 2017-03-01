species=$1

################################################
# Feb 28, 2017
# Use Blast to map genes to reference genomes
# Note that paths are specific to Pollard Sever
################################################# 


genome_id=`cat ~/BenNanditaProject/metaphlan2_to_patric/genome_info.txt  | grep -w $species | awk -F  "\t" '$3 == 1' | cut -f1`

dir='/pollard/shattuck0/ngarud/'
# convert B. uniformis PATRIC representative genome to uppercase
file=`ls /pollard/shattuck0/snayfach/databases/PATRIC/genomes/${genome_id}/${genome_id}.PATRIC.ffn*`

if file --mime-type "$file" | grep -q gzip$; then  

    zcat /pollard/shattuck0/snayfach/databases/PATRIC/genomes/${genome_id}/${genome_id}.PATRIC.ffn.gz | awk '{ if ($0 !~ />/) {print toupper($0)} else {print $0} }'  > ~/BenNanditaProject/metaphlan2_to_patric/blast/tmp.fa

else
    cat /pollard/shattuck0/snayfach/databases/PATRIC/genomes/${genome_id}/${genome_id}.PATRIC.ffn | awk '{ if ($0 !~ />/) {print toupper($0)} else {print $0} }'  > ~/BenNanditaProject/metaphlan2_to_patric/blast/tmp.fa

fi

# make the datbase using the B. uniformis representative genome
makeblastdb -in ~/BenNanditaProject/metaphlan2_to_patric/blast/tmp.fa -out ~/BenNanditaProject/metaphlan2_to_patric/blast/${genome_id}_db -dbtype nucl

# Metaphlan2 genes:
# grep the B. uniformis genes from markers_info.txt
genus=`echo $species | cut -f1 -d'_'`
species_name=`echo $species | cut -f2 -d'_'`

cat ~/BenNanditaProject/metaphlan2_to_patric/markers_info.txt | grep $genus | grep $species_name | cut -f1 > ~/BenNanditaProject/metaphlan2_to_patric/blast/${species}_metaphlan2_genes.txt

# now match these genes to what is in markers.fasta
genes=`cat ~/BenNanditaProject/metaphlan2_to_patric/blast/${species}_metaphlan2_genes.txt`
samtools faidx ~/BenNanditaProject/metaphlan2_to_patric/blast/markers.fasta $genes > ~/BenNanditaProject/metaphlan2_to_patric/blast/${species}_metaphlan2.fa


# run blast
blastn -db ~/BenNanditaProject/metaphlan2_to_patric/blast/${genome_id}_db -query ~/BenNanditaProject/metaphlan2_to_patric/blast/${species}_metaphlan2.fa -outfmt 6 -out ~/BenNanditaProject/metaphlan2_to_patric/blast/${species}_blast_metaphlan2_PATRIC.txt


# list of matching genes:
cat ~/BenNanditaProject/metaphlan2_to_patric/blast/${species}_blast_metaphlan2_PATRIC.txt | cut -f2 | cut -f2 -d'|' | sort | uniq > ~/BenNanditaProject/ben_nandita_hmp_data/metaphlan2_genes/${species}_metaphlan2_genes_mapped.txt

