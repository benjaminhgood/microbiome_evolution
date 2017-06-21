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

species_names = parse_midas_data.parse_good_species_list() 
# load species info
files={}

for species_name in  species_names:
    if (os.path.exists(os.path.expanduser('~/tmp_intermediate_files/ortholog_%s.dat' % species_name))):
        files[species_name]=pickle.load(open(os.path.expanduser('~/tmp_intermediate_files/ortholog_%s.dat' % species_name),'rb'))
    else:
        print species_name


# iterate through the species, add the counts to vectors
num_genes_with_match_within_host=[]
num_genes_with_match_random=[]

for species_name in files:
    num_genes_with_match_within_host.append(files[species_name]['num_genes_with_match_within_host'])
    num_genes_with_match_random.append(files[species_name]['num_genes_with_match_random'])


#plot


pylab.figure()
pylab.xlabel('number of within-host genes with a 100% match')
pylab.ylabel('number of random genes with a 100% match')
pylab.title('Number of genes with a 100% match to an off-target species')

pylab.plot(num_genes_with_match_within_host, num_genes_with_match_random, "ro")
pylab.plot([0,800],[0,800] , 'k-')
pylab.savefig('%s/within_vs_random_gene_orthologs.png' % (parse_midas_data.analysis_directory),bbox_inches='tight') 


