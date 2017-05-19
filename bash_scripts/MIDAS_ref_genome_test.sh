
########################
# May 5, 2017
# Rerun midas using the 6 differnt reference genomes.
# use bowtie mapping threshold of 94%
########################                                                              
# Goal: replace the following two files genome.features.gz  genome.fna.gz in /pollard/home/ngarud/midas_db/rep_genomes/Bacteroides_uniformis_57318
                                                                                      # First, identify which are the 7 different genomes

cat /pollard/home/ngarud/midas_db/genome_info.txt | grep Bacteroides_uniformis| cut -f1
1235787.3
1339348.3
411479.10
457393.3
585543.3
997889.3
997890.3

# create a new midas DB for each reference genome

for genome_id in 1235787.3 1339348.3 457393.3 585543.3 997889.3 997890.3; do
    nohup nice cp -R ~/ben_nandita_hmp_data/midas_db ~/BenNanditaProject/MIDAS_ref_genome_test/midas_db_new_ref_genomes/midas_db_${genome_id} &
done


# create genome_features files
python ~/ben_nandita_hmp_scripts/create_genome_features_file.py Bacteroides_uniformis_57318

# copy the new genomes in PATRIC  and the new genome_features files   
for genome_id in 1235787.3 1339348.3 457393.3 585543.3 997889.3 997890.3; do
    nohup nice cp /pollard/shattuck0/snayfach/databases/PATRIC/genomes/${genome_id}/${genome_id}.fna.gz ~/BenNanditaProject/MIDAS_ref_genome_test/midas_db_new_ref_genomes/midas_db_${genome_id}/rep_genomes/Bacteroides_uniformis_57318/genome.fna.gz &

    nohup nice cp ~/BenNanditaProject/MIDAS_ref_genome_test/genome_features_files/${genome_id}_features.gz ~/BenNanditaProject/MIDAS_ref_genome_test/midas_db_new_ref_genomes/midas_db_${genome_id}/rep_genomes/Bacteroides_uniformis_57318/genome.features.gz &
done

# run MIDAS
qsub ~/ben_nandita_hmp_scripts/qsub_scripts/qsub_script_HMP_MIDAS_snps_diff_ref_genomes

# create bai file. 
for genome in 1235787.3 1339348.3 457393.3 585543.3 700037453 997889.3 997890.3; do 
    cd /pollard/home/ngarud/BenNanditaProject/MIDAS_ref_genome_test/MIDAS_output/${genome}/700037453/snps/temp
    nohup nice samtools index genomes.bam  genomes.bam.bai &
done

