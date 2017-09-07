import parse_midas_data
import parse_patric
import sys
import numpy
from numpy.random import normal
import diversity_utils
import stats_utils
import os
import shutil

# iterate through the genome_metadata file from patric to identify HMP genomes.
HMP_genomes=parse_patric.get_HMP_reference_genomes()

# get a list of genome_ids:
genome_ids=parse_midas_data.genome_ids_dictionary()

# cluster each genome in HMP_genomes according to their species ID
# aggregate the representive genome flag to see if at least one HMP representative genome is being used as the reference genome. 
genome_names={}
genome_names_hmp_ids={}
for genome in HMP_genomes:
    if genome in genome_ids:
        species_name=genome_ids[genome][0]
        representative= genome_ids[genome][1]
        if species_name not in genome_names:
            genome_names[species_name] = [representative]
            genome_names_hmp_ids[species_name]=[genome]
        else:
            genome_names[species_name].append(representative)
            genome_names_hmp_ids[species_name].append(genome)

# iterate through the genome_names dict to see which species are not using an HMP genome when there is one available. If there are none using HMP, then pick the one with the most bps. Store this genome choice in a dictopnary
genome_choice={}

for species_name in genome_names:
    if '1' not in genome_names[species_name]:
        # iterate through and find the genome with the most bps
        genome_max_bps=''
        maxbps=0
        for genome_id in genome_names_hmp_ids[species_name]:
            if HMP_genomes[genome_id][1] > maxbps:
                maxbps=HMP_genomes[genome_id][1]
                genome_max_bps=genome_id
            genome_choice[species_name]=genome_max_bps


#write new genome features file for each of the genomes in genome_choice:

for species in genome_choice:
    genome_id=genome_choice[species]
    #os.system('mkdir /pollard/home/ngarud/BenNanditaProject/MIDAS_HMP_ref_genome_swap/%s' % species)
    #outFile='/pollard/home/ngarud/BenNanditaProject/MIDAS_HMP_ref_genome_swap/%s/%s_features.gz' % (species, genome_id) # note that paths are updated (NRG 09/06/07)
    outFile='/pollard/shattuck0/ngarud/midas_db_HMP_refs/rep_genomes/%s/genome.features.gz' % (species)
    parse_patric.new_genome_features_file(genome_id, outFile)
    

# may need to update the paths below. (NRG 09/06/17)

# swap out all the fasta files and genome features files in the existing MIDAS db
for species in genome_choice:
    genome_id=genome_choice[species]
    #copy the fasta file
    shutil.copyfile('/pollard/shattuck0/snayfach/databases/PATRIC/genomes/%s/%s.fna.gz' %(genome_id,genome_id), '/pollard/home/ngarud/midas_db_HMP_refs/rep_genomes/%s/genome.fna.gz' %species)
    shutil.copyfile('/pollard/home/ngarud/BenNanditaProject/MIDAS_HMP_ref_genome_swap/%s/%s_features.gz' %(species, genome_id), '/pollard/home/ngarud/midas_db_HMP_refs/rep_genomes/%s/genome.features.gz' %species)

#update the genome_info.txt file with an indicator of whether or not the genome is the rep. genome

# read in the species_info.txt file 
inFN='/pollard/home/ngarud/midas_db_HMP_refs/genome_info.txt'
inFile=open(inFN,'r')

# outfile for the new species_info file
outFN='/pollard/home/ngarud/midas_db_HMP_refs/genome_info_new.txt'
outFile=open(outFN,'w')

header=inFile.readline()
outFile.write(header)

# write in the crassphage genome:
outFile.write("NC_024711\tcrassphage\t1\t97000\t1\tcrassphage\n")

for line in inFile:
    items=line.strip().split('\t')
    genome_id = items[0]
    genome_name=items[1]
    rep_genome=items[2]
    length=items[3]
    contigs=items[4]
    species_id=items[5] 
    if species_id in genome_choice:
        if genome_id == genome_choice[species_id]:
            rep_genome='1'
        else:
            rep_genome='0'
        outFile.write(genome_id + '\t' + genome_name + '\t' + rep_genome + '\t' + length + '\t' + contigs + '\t' + species_id +'\n')
    else:
        outFile.write(line)
    

