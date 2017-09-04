import os
import sys
import parse_midas_data

#outFile=open('/pollard/home/ngarud/BenNanditaProject/MIDAS_intermediate_files_hmp/subject_ids_hmp.txt','w')

outFile=open('/pollard/home/ngarud/BenNanditaProject/MIDAS_intermediate_files_hmp/sample_ids_unique_hosts_hmp.txt','w')

subject_sample_time_map = parse_midas_data.parse_subject_sample_time_map() 
for subject in subject_sample_time_map.keys(): # loop over subjects (hosts)
#    outFile.write(subject +'\n')
    times=subject_sample_time_map[subject].keys()
    if len(subject_sample_time_map[subject][times[0]]) >1:
        outFile.write(subject_sample_time_map[subject][times[0]][0][0] + 'c' + '\n')
    else:
        outFile.write(subject_sample_time_map[subject][times[0]][0][0] +'\n')
