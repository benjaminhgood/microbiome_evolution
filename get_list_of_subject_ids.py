import os
import sys
import parse_midas_data

outFile=open('/pollard/home/ngarud/BenNanditaProject/MIDAS_intermediate_files_hmp/subject_ids_hmp.txt','w')

subject_sample_time_map = parse_midas_data.parse_subject_sample_time_map() 
for subject in subject_sample_time_map.keys(): # loop over subjects (hosts)
    outFile.write(subject +'\n')
