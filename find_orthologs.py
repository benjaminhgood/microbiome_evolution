import matplotlib  
matplotlib.use('Agg') 
import pylab
import sys
import numpy
import stats_utils
import matplotlib.colors as colors
import matplotlib.cm as cmx
from math import log10,ceil
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy.random import randint
import pandas
import os
import parse_midas_data
import pickle
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


################
# definitions  #
################


######################################
# read in ortholog information       #
######################################      
def read_ortholog_information(infile):
    gene_statistics={}
    matching_genome_list={}
    num_genes_with_match=0 # I am storing the number of genes with a match to any species other than the species of interest. THe match must be 100%
    for line in infile:
        items=line.strip().split()
        gene=items[0]
        matching_genome=items[1]
        percent_id=float(items[2])
        if gene not in gene_statistics:
            gene_statistics[gene]=[0,0,0] # first is for number of genomes matching, second is for number of genomes matching with 100%, # thurs is for whether the speceis matches at 100%.
        if percent_id == 100.00:
            if matching_genome==species_name:
                gene_statistics[gene][2]+=1
            else:
                gene_statistics[gene][1]+=1
                num_genes_with_match +=1
                if matching_genome not in matching_genome_list:
                    matching_genome_list[matching_genome]=1
                else:
                    matching_genome_list[matching_genome]+=1
        else:
            gene_statistics[gene][0]+=1
    return gene_statistics, matching_genome_list, num_genes_with_match

#####################
# main part of code #
#####################

# genes changing within hosts
infile=open(os.path.expanduser("~/tmp_intermediate_files/" + species_name + "_all_within_host_changes_matching_genome.txt"),'r') 

gene_statistics_within_host, matching_genome_list_within_host, num_genes_with_match_within_host = read_ortholog_information(infile)

# repeat with random genes
infile=open(os.path.expanduser("~/tmp_intermediate_files/" + species_name + "_all_random_gene_set_matching_genome.txt"),'r') 

gene_statistics_random, matching_genome_list_random, num_genes_with_match_random = read_ortholog_information(infile)

###########################################################################
# Plot: histogram showing number of orthologs
###########################################################################  

label_size = 10
matplotlib.rcParams['ytick.labelsize'] = label_size
width = 1.0 

df=pandas.DataFrame(gene_statistics_within_host)
df=df.transpose()
num_genes_with_match_within_host=(df[1]!= 0).sum() 
df=df.rename(columns = {0:'< 100%, non-ref match', 1: '100%, non-ref match', 2: '100%, ref match'})
ax =df.plot(kind='barh', stacked=True, figsize=(6,10), title=species_name, width=width)
ax.set_xlabel("Number of orthologs")
ax.set_ylabel("Gene")

pylab.savefig('%s/%s_within_host_gene_orthologs.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight') 


############################################
# Histogram of matching species statistics #
############################################

df=pandas.DataFrame(matching_genome_list_within_host, index=[0])
df=pandas.melt(df)

label_size = 8
matplotlib.rcParams['ytick.labelsize'] = label_size 
pos = numpy.arange(len(df.index))
width = 1.0 

ax =df.plot(kind='barh', stacked=True, figsize=(6,10), title=species_name, width=width)
ax.set_xlabel("Number of genes")
ax.set_ylabel("Species")

ax.set_yticks(pos + (width / 2))
ax.set_yticklabels(df['variable'].tolist())

pylab.savefig('%s/%s_within_host_gene_orthologs_species_distribution.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight') 

###########################################################################
# Plot: histogram showing number of orthologs
###########################################################################  

label_size = 10
matplotlib.rcParams['ytick.labelsize'] = label_size

df=pandas.DataFrame(gene_statistics_random)
df=df.transpose()
num_genes_with_match_random=(df[1]!= 0).sum()
df=df.rename(columns = {0:'< 100%, non-ref match', 1: '100%, non-ref match', 2: '100%, ref match'})
ax =df.plot(kind='barh', stacked=True, figsize=(6,10), title=species_name, width=width)
ax.set_xlabel("Number of orthologs")
ax.set_ylabel("Gene")


pylab.savefig('%s/%s_random_gene_orthologs.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight') 



############################################
# Histogram of matching species statistics #
############################################

df=pandas.DataFrame(matching_genome_list_random, index=[0])
df=pandas.melt(df)

label_size = 8
matplotlib.rcParams['ytick.labelsize'] = label_size 
pos = numpy.arange(len(df.index))
width = 1.0 

ax =df.plot(kind='barh', stacked=True, figsize=(6,10), title=species_name, width=width)
ax.set_xlabel("Number of genes")
ax.set_ylabel("Species")

ax.set_yticks(pos + (width / 2))
ax.set_yticklabels(df['variable'].tolist())

pylab.savefig('%s/%s_random_gene_orthologs_species_distribution.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight') 



##################
# save the number of genes with match in a pickle. Make a plot later looking accross species

saved_data=dict(
num_genes_with_match_within_host=num_genes_with_match_within_host,
num_genes_with_match_random=num_genes_with_match_random)

saved_data_file=os.path.expanduser('~/tmp_intermediate_files/ortholog_%s.dat' % species_name)
with open(saved_data_file, 'wb') as outfile:
    pickle.dump(saved_data, outfile, protocol=pickle.HIGHEST_PROTOCOL)
