import os
import sys
import parse_midas_data

subject_sample_map = parse_midas_data.parse_subject_sample_map()

for subject in subject_sample_map.keys(): # loop over subjects (hosts)
    for sample in subject_sample_map[subject].keys(): # loop over samples
        print "Accessions for sample", sample , "and subject", subject
        os.system('rm /netapp/home/ngarud/shattuck/BenNanditaProject/joined_fastq_files_hmp/' + sample + '_1.fastq.gz')
        os.system('rm /netapp/home/ngarud/shattuck/BenNanditaProject/joined_fastq_files_hmp/' + sample + '_2.fastq.gz')
        
        for accession in subject_sample_map[subject][sample]:
            print accession
            os.system('cat /pollard/shattuck0/snayfach/metagenomes/HMP/fastq/' + accession + '_1.fastq.gz >> /netapp/home/ngarud/shattuck/BenNanditaProject/MIDAS_intermediate_files_hmp/joined_fastq_files_hmp_combine_tech_reps/' + sample + '_1.fastq.gz')
            os.system('cat /pollard/shattuck0/snayfach/metagenomes/HMP/fastq/' + accession + '_2.fastq.gz >> /netapp/home/ngarud/shattuck/BenNanditaProject/MIDAS_intermediate_files_hmp/joined_fastq_files_hmp_combine_tech_reps/' + sample + '_2.fastq.gz')

