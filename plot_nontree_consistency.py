import matplotlib  
matplotlib.use('Agg') 
import sample_utils
import config
import parse_midas_data
import os.path
import pylab
import sys
import numpy
from numpy.random import choice

import species_phylogeny_utils
import diversity_utils
import gene_diversity_utils
import calculate_substitution_rates
import calculate_snv_distances
import figure_utils

import stats_utils
import matplotlib.colors as colors
import matplotlib.cm as cmx
from math import log10,ceil,log,exp
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy.random import randint
from mpl_toolkits.axes_grid.inset_locator import inset_axes

from numpy.random import randint, binomial, choice
from scipy.stats import poisson as poisson_distribution
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster.hierarchy import cophenet
from scipy.cluster.hierarchy import fcluster

from scipy.stats import gaussian_kde

mpl.rcParams['font.size'] = 5
mpl.rcParams['axes.labelpad'] = 2
mpl.rcParams['lines.linewidth'] = 0.5
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'


################################################################################
#
# Standard header to read in argument information
#
################################################################################
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--debug", help="Loads only a subset of SNPs for speed", action="store_true")
parser.add_argument("--chunk-size", type=int, help="max number of records to load", default=1000000000)

args = parser.parse_args()

debug = args.debug
chunk_size = args.chunk_size

################################################################################

low_divergence_threshold = 2e-04
min_sample_size = 33 # 46 gives at least 1000 pairs, 33 gives at least 500 (actually 528)
allowed_variant_types = set(['4D'])

#####
#
# Sset up figure
#
#####

pylab.figure(1,figsize=(3.42,2))
pylab.xlabel('Divergence, $d$')
pylab.ylabel('Phylogenetic inconsistency')
pylab.xlim([1e-04,1e-01])
pylab.ylim([0,1.05])

#####
#
# Do calculation
#
#####

good_species_list = parse_midas_data.parse_good_species_list()
if debug:
    good_species_list = good_species_list[0:3]
    
sys.stderr.write("Loading sample metadata...\n")
subject_sample_map = sample_utils.parse_subject_sample_map()
sample_continent_map = sample_utils.parse_sample_continent_map()
sys.stderr.write("Done!\n")

for species_name in good_species_list:

    sys.stderr.write("Loading haploid samples...\n")
    # Only plot samples above a certain depth threshold that are "haploids"
    snp_samples = diversity_utils.calculate_haploid_samples(species_name, debug=debug)
    
    if len(snp_samples) < min_sample_size:
        sys.stderr.write("Not enough haploid samples!\n")
        continue
        
    sys.stderr.write("Calculating unique samples...\n")
    # Only consider one sample per person
    snp_samples = snp_samples[sample_utils.calculate_unique_samples(subject_sample_map, sample_list=snp_samples)]


    if len(snp_samples) < min_sample_size:
        sys.stderr.write("Not enough unique samples!\n")
        continue


    # Load SNP information for species_name
    sys.stderr.write("Loading SNPs for %s...\n" % species_name)
    sys.stderr.write("(core genes only...)\n")
    snv_distance_map = calculate_snv_distances.load_snv_distance_map(species_name)
    
    ds = numpy.logspace(-4,-1.5,50) # 15 points are plotted
    total_snps = numpy.zeros_like(ds)
    inconsistent_snps = numpy.zeros_like(ds)
      
    for location_tuple in snv_distance_map:
        var_type, derived_allele_counts, ancestral_allele_counts, between_d, within_d1, within_d2 = snv_distance_map[location_tuple]
        
        if var_type in allowed_variant_types:
            
            within_d = min([within_d1, within_d2])
            good_idxs = (ds>=between_d)
            inconsistent_idxs = good_idxs*(within_d>=2*ds)
            
            total_snps[good_idxs] += 1
            inconsistent_snps[inconsistent_idxs] += 1
    
    fraction_inconsistent = inconsistent_snps*1.0/(total_snps+(total_snps==0))
    
    print ds
    print total_snps
    print inconsistent_snps
    
    
    if species_name=="Bacteroides_vulgatus_57955":
        color = 'b'
        linewidth=1
        zorder=2
        alpha=1
    else:
        color = 'r'
        alpha =0.3
        linewidth=0.5
        zorder = 1
    pylab.semilogx(ds[total_snps>0], fraction_inconsistent[total_snps>0],'-',color=color,linewidth=linewidth,zorder=zorder,alpha=alpha)

pylab.plot([1],[-1],'b-',linewidth=1, alpha=1,label='B. vulgatus')
pylab.plot([1],[-1],'r-',linewidth=0.5, alpha=0.3,label='Other species')
pylab.legend(loc='upper right',frameon=False,fontsize=5,numpoints=1,handlelength=1)
pylab.savefig('non_tree_consistency.pdf',bbox_inches='tight')