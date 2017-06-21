species=$1



############################################################################
# Blast within-host gene changes against other genomes. Is there a match?  #
############################################################################



#################################################
# run blast against the genes that are changing #
#################################################


file=~/ben_nandita_hmp_analysis/${species}_within_host_gene_changes.txt

# need to get the sequences of the genes that are changing. Grep the gene names and reads from this file. 
rm ~/tmp_intermediate_files/${species}_within_host_gene_changes.fasta 
while read gene; do
    samtools faidx ~/ben_nandita_hmp_data/midas_db/pan_genomes/${species}/centroids.ffn.gz $gene >> ~/tmp_intermediate_files/${species}_within_host_gene_changes.fasta
done < $file


# Run blast
blastn -db ~/tmp_intermediate_files/all_species_pan_genomes_genes_db  -query ~/tmp_intermediate_files/${species}_within_host_gene_changes.fasta -outfmt 6 -out ~/tmp_intermediate_files/${species}_all_within_host_gene_changes_blast.txt 

# check which species' genomes that the genes are mapping to. Are they random or closely related genomes?
rm ~/tmp_intermediate_files/${species}_all_within_host_changes_matching_genome.txt 
while read line; do
    gene=`echo $line | cut -f1 -d' '`
    first_no=`echo $line | cut -f2 -d' '| cut -f1 -d'.'`
    second_no=`echo $line | cut -f2 -d' '| cut -f2 -d'.'`
    percent=`echo $line | cut -f3 -d' '`
    bit_score=`echo $line | cut -f12 -d' '`
    matching_genome=`cat ~/ben_nandita_hmp_data/midas_db/genome_info.txt | grep -w ${first_no}.${second_no} | cut -f6`
    echo -e "$gene\t$matching_genome\t$percent\t$bit_score" >> ~/tmp_intermediate_files/${species}_all_within_host_changes_matching_genome.txt
done < ~/tmp_intermediate_files/${species}_all_within_host_gene_changes_blast.txt



#####################################################################################
# repeat for a random list of genes from the pangenome of the species of interest.  #
#####################################################################################
number_genes=`wc -l $file | cut -f1 -d' '`
zcat ~/ben_nandita_hmp_data/midas_db/pan_genomes/${species}/centroids.ffn.gz | grep '>' | shuf -n $number_genes | cut -f2 -d'>' > ~/tmp_intermediate_files/${species}_random_gene_set.txt

# need to get the sequences of the random genes

rm ~/tmp_intermediate_files/${species}_random_gene_set.fasta 
while read gene; do
    samtools faidx ~/ben_nandita_hmp_data/midas_db/pan_genomes/${species}/centroids.ffn.gz $gene >> ~/tmp_intermediate_files/${species}_random_gene_set.fasta
done < ~/tmp_intermediate_files/${species}_random_gene_set.txt


# Run blast on the random gene set 
blastn -db ~/tmp_intermediate_files/all_species_pan_genomes_genes_db  -query ~/tmp_intermediate_files/${species}_random_gene_set.fasta -outfmt 6 -out ~/tmp_intermediate_files/${species}_all_random_gene_set_blast.txt 

# check which species' genomes that the genes are mapping to. Are they random or closely related genomes?
rm ~/tmp_intermediate_files/${species}_all_random_gene_set_matching_genome.txt
while read line; do
    gene=`echo $line | cut -f1 -d' '`
    first_no=`echo $line | cut -f2 -d' '| cut -f1 -d'.'`
    second_no=`echo $line | cut -f2 -d' '| cut -f2 -d'.'`
    percent=`echo $line | cut -f3 -d' '`
    bit_score=`echo $line | cut -f12 -d' '`
    matching_genome=`cat ~/ben_nandita_hmp_data/midas_db/genome_info.txt | grep -w ${first_no}.${second_no} | cut -f6`
    echo -e "$gene\t$matching_genome\t$percent\t$bit_score" >> ~/tmp_intermediate_files/${species}_all_random_gene_set_matching_genome.txt
done < ~/tmp_intermediate_files/${species}_all_random_gene_set_blast.txt


# Python script: 
# quantify the number of species to which the gene can map to with 100% accuracy. 

python ~/ben_nandita_hmp_scripts/find_orthologs.py $species


