import os
import sys
import parse_midas_data
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("subject_id", help="name of subject to process")
args = parser.parse_args()
subject = args.subject_id

if subject.split('-')[0] == 'FAT_DON_11':
    subject='FAT_DON_11'

# dictionary with all the fastq file names for the sample of interest
files={}
inFN="/netapp/home/ngarud/shattuck/BenNanditaProject/MIDAS_intermediate_files_Li_2016_FMT/Li_2016_FMT_paired_accessions.txt"
inFile=open(inFN, 'r')

for line in inFile:
    items=line.strip().split()
    accession=items[0]
    subject_id=items[1]
    subject_id_split=subject_id.split('-')
    if subject_id_split[0]=='FAT_DON_11':
        subject_id=subject_id_split[0]
    if subject_id not in files:
        files[subject_id]=[]
    files[subject_id].append(accession)


os.system('rm /netapp/home/ngarud/shattuck/BenNanditaProject/MIDAS_intermediate_files_Li_2016_FMT/joined_fastq_files_hmp_combine_tech_reps/' + subject + '_1.fastq.gz')
os.system('rm /netapp/home/ngarud/shattuck/BenNanditaProject/MIDAS_intermediate_files_Li_2016_FMT/joined_fastq_files_hmp_combine_tech_reps/' + subject + '_2.fastq.gz')

# iterate through the accessions identified for subject of interest and concatenate them. 
for accession in files[subject]:
    os.system('cat /netapp/home/ngarud/shattuck/metagenomic_fastq_files/Li_2016_FMT/fastq_files/' + accession + '_1.fastq.gz >> /netapp/home/ngarud/shattuck/BenNanditaProject/MIDAS_intermediate_files_Li_2016_FMT/joined_fastq_files_hmp_combine_tech_reps/' + subject + '_1.fastq.gz')

    os.system('cat /netapp/home/ngarud/shattuck/metagenomic_fastq_files/Li_2016_FMT/fastq_files/' + accession + '_2.fastq.gz >> /netapp/home/ngarud/shattuck/BenNanditaProject/MIDAS_intermediate_files_Li_2016_FMT/joined_fastq_files_hmp_combine_tech_reps/' + subject + '_2.fastq.gz')


