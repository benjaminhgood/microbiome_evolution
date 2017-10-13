import parse_midas_data
import parse_patric
import sys
import numpy
from numpy.random import normal
import diversity_utils
import stats_utils
import os

subject_sample_time_map = parse_midas_data.parse_subject_sample_time_map()
#covert subject_sample_time_map so that visno replicates are gone. 
subject_sample_time_map = parse_midas_data.collapse_visno_reps_subject_sample_time_map(subject_sample_time_map)

# get the species list for each time pair
for subject_id in subject_sample_time_map:
    species_lists=[]
    for timept in subject_sample_time_map[subject_id]:
        sample_id=subject_sample_time_map[subject_id][timept][0][0]
        last_char=sample_id[-1]
        if last_char=='c':
            dir='MIDAS_1.2.2_samples_combined_output'
            fastq_dir='joined_fastq_files_hmp_combine_sample_reps'
        else:
            dir='MIDAS_1.2.2_output'
            fastq_dir='joined_fastq_files_hmp_combine_tech_reps'
        inFN=os.path.expanduser('~/BenNanditaProject/MIDAS_intermediate_files_hmp/' + dir + '/%s/species/species_profile.txt') %sample_id
        species_list=parse_midas_data.parse_intermediate_species_file(sample_id, inFN)
        species_lists.append(species_list)

        # for each sample, I want to rename the snps and gene directories so that I don't overwrite the old db output.
        #move_command = 'mv ' + os.path.expanduser('~/BenNanditaProject/MIDAS_intermediate_files_hmp/' + dir + '/%s/snps/ ~/BenNanditaProject/MIDAS_intermediate_files_hmp/' + dir + '/%s/snps_new_db') % (sample_id, sample_id)  
        #os.system(move_command)
        #move_command = 'mv ' + os.path.expanduser('~/BenNanditaProject/MIDAS_intermediate_files_hmp/' + dir + '/%s/genes/ ~/BenNanditaProject/MIDAS_intermediate_files_hmp/' + dir + '/%s/genes_new_db') % (sample_id, sample_id)  
        #os.system(move_command)

 
    # find the union of the different time points:
    species_union=species_lists[0]
    if len(species_lists)>1:
        for timept in range(1, len(species_lists)):
            species_union=species_lists[timept].union(species_union)

    # additionally, add the crassphage genome to each species_union
    #species_union.add("crassphage")

    # now that we have the list of species to run each sample with, I want to create a file with the species union in it.  
    for timept in subject_sample_time_map[subject_id]:
        sample_id=subject_sample_time_map[subject_id][timept][0][0]
        outFN=os.path.expanduser('~/BenNanditaProject/MIDAS_intermediate_files_hmp/MIDAS_species_lists_time_pairs/%s_species_union.txt') % sample_id 
        outFile=open(outFN,'w')
        outFile.write('species_id\tcount_reads\tcoverage\trelative_abundance\n') # header
        for species in species_union:
            outFile.write(species + '\n')



    # to assess what is actually changing in terms of species composition from one time point to another:
    if len(species_lists)>1:
        time_1_diff=len(species_lists[0]-species_lists[1]) 
        time_2_diff=len(species_lists[1]-species_lists[0]) 
        intersection=len(species_lists[0] & species_lists[1])
        print str(time_1_diff) + '\t' + str(time_2_diff) + '\t' +str(intersection) + '\t' + str(len(species_union)) 
