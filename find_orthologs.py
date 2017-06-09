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

######################################
# read in ortholog information       #
######################################      
infile_gene_changes=open(os.path.expanduser("~/tmp_intermediate_files/" + species_name + "_all_within_host_changes_matching_genome.txt"),'r') 

gene_statistics={}
for line in infile_gene_changes:
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
    else:
        gene_statistics[gene][0]+=1




###########################################################################
# Plot: histogram showing number of orthologs
###########################################################################  

# try with dataframe
pylab.subplot(2, 2, 1)     
label_size = 10
matplotlib.rcParams['ytick.labelsize'] = label_size
width = 1.0 

df=pandas.DataFrame(gene_statistics)
df=df.transpose()
df=df.rename(columns = {0:'< 100%, non-ref match', 1: '100%, non-ref match', 2: '100%, ref match'})
ax =df.plot(kind='barh', stacked=True, figsize=(6,10), title=species_name, width=width)
ax.set_xlabel("Number of orthologs")
ax.set_ylabel("Gene")

pylab.savefig('%s/%s_within_host_gene_orthologs.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight') 


#############################
# repeat with random genes  #
#############################

infile_random_genes=open(os.path.expanduser("~/tmp_intermediate_files/" + species_name + "_all_random_gene_set_matching_genome.txt"),'r') 


gene_statistics={}
for line in infile_random_genes:
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
    else:
        gene_statistics[gene][0]+=1






###########################################################################
# Plot: histogram showing number of orthologs
###########################################################################  

# try with dataframe
pylab.subplot(2, 2, 2)     
label_size = 10
matplotlib.rcParams['ytick.labelsize'] = label_size

df=pandas.DataFrame(gene_statistics)
df=df.transpose()
df=df.rename(columns = {0:'< 100%, non-ref match', 1: '100%, non-ref match', 2: '100%, ref match'})
ax =df.plot(kind='barh', stacked=True, figsize=(6,10), title=species_name, width=width)
ax.set_xlabel("Number of orthologs")
ax.set_ylabel("Gene")


pylab.savefig('%s/%s_random_gene_orthologs.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight') 
