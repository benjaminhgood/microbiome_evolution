import os
import sys
import parse_midas_data

subject_sample_time_map = parse_midas_data.parse_subject_sample_time_map()

outFile=open("/netapp/home/ngarud/shattuck/BenNanditaProject/MIDAS_intermediate_files_hmp/sample_paths.txt",'w')
for subject in subject_sample_time_map.keys(): # loop over subjects (hosts)
    for visno in subject_sample_time_map[subject].keys(): # loop over samples
        if len(subject_sample_time_map[subject][visno]) ==1:
            sample_name=subject_sample_time_map[subject][visno][0][0]
            outFile.write('/netapp/home/ngarud/shattuck/BenNanditaProject/MIDAS_intermediate_files_hmp/MIDAS_1.2.2_output/' + sample_name + '\n')
