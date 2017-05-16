import os
import sys
import parse_midas_data
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("subject_id", help="name of subject to process")
args = parser.parse_args()
subject = args.subject_id

subject_sample_time_map = parse_midas_data.parse_subject_sample_time_map()

for visno in subject_sample_time_map[subject].keys(): # loop over samples
    print "Subject:", subject, "Samples for visno ", visno
    if len(subject_sample_time_map[subject][visno]) >1:
        new_sample_name=subject_sample_time_map[subject][visno][0][0] + 'c'
        print new_sample_name
        os.system('rm /netapp/home/ngarud/shattuck/BenNanditaProject/MIDAS_intermediate_files_hmp/joined_fastq_files_hmp_combine_sample_reps/' + new_sample_name + '_1.fastq.gz')
        os.system('rm /netapp/home/ngarud/shattuck/BenNanditaProject/MIDAS_intermediate_files_hmp/joined_fastq_files_hmp_combine_sample_reps/' + new_sample_name + '_2.fastq.gz')
        for i in range(0,len(subject_sample_time_map[subject][visno])):
            sample = subject_sample_time_map[subject][visno][i][0]
            print sample
            os.system('cat /netapp/home/ngarud/shattuck/BenNanditaProject/MIDAS_intermediate_files_hmp/joined_fastq_files_hmp_combine_tech_reps/' + sample + '_1.fastq.gz >> /netapp/home/ngarud/shattuck/BenNanditaProject/MIDAS_intermediate_files_hmp/joined_fastq_files_hmp_combine_sample_reps/' + new_sample_name + '_1.fastq.gz')
            os.system('cat /netapp/home/ngarud/shattuck/BenNanditaProject/MIDAS_intermediate_files_hmp/joined_fastq_files_hmp_combine_tech_reps/' + sample + '_2.fastq.gz >> /netapp/home/ngarud/shattuck/BenNanditaProject/MIDAS_intermediate_files_hmp/joined_fastq_files_hmp_combine_sample_reps/' + new_sample_name + '_2.fastq.gz' )

