import matplotlib  
matplotlib.use('Agg') 
import parse_midas_data
import parse_patric
import sys
import numpy
from numpy.random import normal
import diversity_utils
import stats_utils
import pylab
import os
import os.path
import gzip
import pandas

################################################################################
#
# Standard header to read in argument information
#
################################################################################
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("species_name", help="name of species to process")
args = parser.parse_args()

species_name = args.species_name

################################################################################

# get a list of marker genes
marker_genes = parse_midas_data.load_marker_genes(species_name)

# load the centroid gene map (maps from the centroid to the ref genome)
centroid_gene_map=parse_midas_data.load_centroid_gene_map(species_name)

# invert the centroid gene map:
centroid_gene_map_inverted={}
for key in centroid_gene_map:
    value=centroid_gene_map[key]
    centroid_gene_map_inverted[value]=key

# load the complete centroid gene map
complete_centroid_gene_map=parse_midas_data.load_complete_centroid_gene_map(species_name)

 
# iterate through the marker genes and initialize dictionaries to store the data in. 

all_marker_genes={} # store the hits to marker genes in this dictionary
all_genomes={} # store the hits to genome (aggregated accross marker genes) in this dictionary. 
all_genomes_all_samples={} 
for gene in marker_genes:
    #map refernece genome marker gene to the centroid gene
    centroid_gene=centroid_gene_map_inverted[gene]
    # find all the genes in the 95% centroid
    all_genes_in_centroid=complete_centroid_gene_map[centroid_gene]
    for gene_name in all_genes_in_centroid:
        all_marker_genes[gene_name] = 0
        first_no=gene_name.split('.')[0]
        second_no=gene_name.split('.')[1]
        genome=first_no + '.' + second_no
        if genome not in all_genomes:
            all_genomes[genome]=0
            all_genomes_all_samples[genome]=[]

#outFile to store all the results for each sample:

outFN_blast=('%s/%s_blast_parse.txt' % (os.path.expanduser("~/tmp_intermediate_files"),species_name))
outFile_blast=open(outFN_blast,'w')

            
# get a list of samples
inFN_samples=(os.path.expanduser("~/BenNanditaProject/MIDAS_intermediate_files_hmp/List_of_non_sample_reps.txt"))
inFile_samples=open(inFN_samples,'r')

samples=[]
for line in inFile_samples:
    samples.append(line.strip())
              
# read in the blast output and find the best hit to marker genes:
for sample in samples:

    # clear the all_genomes dict:
    for genome in all_genomes:
        all_genomes[genome]=0
    inFN='%s/%s_%s_1_marker_genes_blast.txt' % (os.path.expanduser("~/tmp_intermediate_files"), species_name, sample)
    if os.path.isfile(inFN): 
 #       print sample
        inFile=open(inFN,'r')
    
        query_prev=''
        for line in inFile:
            items=line.strip().split('\t')
            query=items[0]
            marker_gene=items[1]
            bit_score=float(items[11])
            first_no=marker_gene.split('.')[0]
            second_no=marker_gene.split('.')[1]
            genome=first_no + '.' + second_no

            if query !=query_prev:
                max_bit_score=bit_score
                all_marker_genes[marker_gene] +=1
                all_genomes[genome]+=1
            elif bit_score==max_bit_score:
                all_marker_genes[marker_gene] +=1
                all_genomes[genome] +=1
            query_prev=query

        # add the results to a full dictionary
        for genome in all_genomes:
            all_genomes_all_samples[genome].append(all_genomes[genome])

df=pandas.DataFrame(all_genomes_all_samples)
df.boxplot()
pylab.savefig('%s/%s_reads_mapping.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight',figsize=(18,10)) 

print(df)
#print all_genomes_all_samples

print(df.sum(axis=0))
