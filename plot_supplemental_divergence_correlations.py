import matplotlib  
matplotlib.use('Agg') 
import config
import parse_midas_data
import parse_HMP_data
import os.path
import pylab
import sys
import numpy

import diversity_utils
import gene_diversity_utils
import calculate_substitution_rates

import stats_utils
import matplotlib.colors as colors
import matplotlib.cm as cmx
from math import log10,ceil
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy.random import randint

from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster.hierarchy import cophenet
from scipy.cluster.hierarchy import fcluster

from scipy.stats import gaussian_kde

mpl.rcParams['font.size'] = 6
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

good_species_list = ['Bacteroides_vulgatus_57955', 'Bacteroides_uniformis_57318', 'Alistipes_putredinis_61533']

####################################################
#
# Set up Figure (3 panels, arranged in 1x3 grid)
#
####################################################

pylab.figure(1,figsize=(7,1.5))
fig = pylab.gcf()
# make three panels panels
outer_grid  = gridspec.GridSpec(1,3,width_ratios=[1,1,1],wspace=0.1)

#######
#
# SNP divergence vs Gene divergence in B. vulgatus
#
#######
gene_axis = plt.Subplot(fig, outer_grid[0])
fig.add_subplot(gene_axis)

gene_axis.set_ylabel('SNP divergence\n %s' % (good_species_list[0]))
gene_axis.set_xlabel('Gene divergence\n %s' % (good_species_list[0]))

gene_axis.set_ylim([1e-06,1e-01])
#gene_axis.set_xlim([1e-02,1])

gene_axis.spines['top'].set_visible(False)
gene_axis.spines['right'].set_visible(False)
gene_axis.get_xaxis().tick_bottom()
gene_axis.get_yaxis().tick_left()

#######
#
# SNP divergence (B vulgatus) vs SNP divergence (A putredinis)
#
#######
species_axis_1 = plt.Subplot(fig, outer_grid[1])
fig.add_subplot(species_axis_1)

species_axis_1.set_xlabel('SNP divergence\n %s' % (good_species_list[1]))

species_axis_1.set_ylim([1e-06,1e-01])
species_axis_1.set_xlim([1e-06,1e-01])

species_axis_1.spines['top'].set_visible(False)
species_axis_1.spines['right'].set_visible(False)
species_axis_1.get_xaxis().tick_bottom()
species_axis_1.get_yaxis().tick_left()



#######
#
# SNP divergence (B vulgatus) vs SNP divergence (A putredinis)
#
#######
species_axis_2 = plt.Subplot(fig, outer_grid[2])
fig.add_subplot(species_axis_2)

species_axis_2.set_xlabel('SNP divergence\n %s' % (good_species_list[2]))

species_axis_2.set_ylim([1e-06,1e-01])
species_axis_2.set_xlim([1e-06,1e-01])

species_axis_2.spines['top'].set_visible(False)
species_axis_2.spines['right'].set_visible(False)
species_axis_2.get_xaxis().tick_bottom()
species_axis_2.get_yaxis().tick_left()


########
#
# Now do calculation and plot figures
#
########

sys.stderr.write("Loading sample metadata...\n")
subject_sample_map = parse_HMP_data.parse_subject_sample_map()
sample_order_map = parse_HMP_data.parse_sample_order_map()
sys.stderr.write("Done!\n")

snp_divergence_map = {species_name: {} for species_name in good_species_list}
gene_divergence_map = {species_name: {} for species_name in good_species_list}
for species_name in good_species_list:

    sys.stderr.write("Loading haploid samples...\n")
    snp_samples = diversity_utils.calculate_haploid_samples(species_name, debug=debug)
        
    sys.stderr.write("Calculating unique samples...\n")
    # Only consider one sample per person
    snp_samples = snp_samples[parse_midas_data.calculate_unique_samples(subject_sample_map, sample_list=snp_samples)]
        
    sys.stderr.write("Loading pre-computed substitution rates for %s...\n" % species_name)
    substitution_rate_map = calculate_substitution_rates.load_substitution_rate_map(species_name)
    sys.stderr.write("Calculating snp matrix...\n")
    dummy_samples, snp_difference_matrix, snp_opportunity_matrix = calculate_substitution_rates.calculate_matrices_from_substitution_rate_map(substitution_rate_map, 'core', allowed_samples=snp_samples)
    snp_samples = dummy_samples
    sys.stderr.write("Done!\n")
    
    sys.stderr.write("Calculating gene matrix...\n")
    gene_samples, gene_difference_matrix, gene_opportunity_matrix = calculate_substitution_rates.calculate_matrices_from_substitution_rate_map(substitution_rate_map, 'genes', allowed_samples=snp_samples)
    snp_samples = gene_samples
    sys.stderr.write("Done!\n")
    
    # Focus on the subset of samples that have sufficient gene depth and snp depth
    desired_samples = gene_samples
    
    # Figure out which pairs of indices in desired_samples belong to diff subjects
    desired_same_sample_idxs, desired_same_subject_idxs, desired_diff_subject_idxs =  parse_midas_data.calculate_subject_pairs( subject_sample_map, desired_samples)

    # Turn these into indices for snp and gene matrices
    snp_sample_idx_map = parse_midas_data.calculate_sample_idx_map(desired_samples, snp_samples)
    gene_sample_idx_map = parse_midas_data.calculate_sample_idx_map(desired_samples, gene_samples)
  
    same_subject_snp_idxs = parse_midas_data.apply_sample_index_map_to_indices(snp_sample_idx_map, desired_same_subject_idxs)  
    same_subject_gene_idxs = parse_midas_data.apply_sample_index_map_to_indices(gene_sample_idx_map, desired_same_subject_idxs)  

    diff_subject_snp_idxs = parse_midas_data.apply_sample_index_map_to_indices(snp_sample_idx_map, desired_diff_subject_idxs)  
    diff_subject_gene_idxs = parse_midas_data.apply_sample_index_map_to_indices(gene_sample_idx_map, desired_diff_subject_idxs)  
    
    for sample_pair_idx in xrange(0,len(diff_subject_snp_idxs[0])):
    
        snp_i = diff_subject_snp_idxs[0][sample_pair_idx]
        snp_j = diff_subject_snp_idxs[1][sample_pair_idx]
        
        gene_i = diff_subject_gene_idxs[0][sample_pair_idx]
        gene_j = diff_subject_gene_idxs[1][sample_pair_idx]
        
        sample_i = desired_samples[gene_i]
        sample_j = desired_samples[gene_j]
        
        # This will serve as a key in snp_divergence_map
        sample_pair = frozenset([sample_i,sample_j])
        
        # Focus on pairs of samples with sufficient coverage 
        if snp_opportunity_matrix[snp_i,snp_j]>0:
            snp_d = snp_difference_matrix[snp_i,snp_j]*1.0/snp_opportunity_matrix[snp_i,snp_j]
            snp_divergence_map[species_name][sample_pair] = snp_d
            
        if gene_opportunity_matrix[gene_i, gene_j]>0:
            gene_d = gene_difference_matrix[gene_i, gene_j]*1.0/gene_opportunity_matrix[gene_i, gene_j]
            gene_divergence_map[species_name][sample_pair] = gene_d


#################
#
# Plot figures!
#
#################

# First calculate SNP vs gene divergence in B. vulgatus
species_name = good_species_list[0]
snp_divergences = []
gene_divergences = []
# Loop over sample pairs that are in both snp_divergence_map and gene_divergence_map
for sample_pair in (set(snp_divergence_map[species_name].keys()) & set(gene_divergence_map[species_name].keys()) ):
    
    snp_divergences.append( snp_divergence_map[species_name][sample_pair] )
    gene_divergences.append( gene_divergence_map[species_name][sample_pair] )

snp_divergences = numpy.array(snp_divergences)
gene_divergences = numpy.array(gene_divergences)

# Null expectation (medians line up)
median_ratio = numpy.median(snp_divergences)/numpy.median(gene_divergences)
gene_axis.loglog([1e-02,1],[1e-02*median_ratio,1*median_ratio],'k-',linewidth=0.25)
    
gene_axis.loglog(gene_divergences, snp_divergences, 'r.', markersize=2,alpha=0.5,markeredgewidth=0, rasterized=True)

# Then SNP divergence between two species
species_1 = good_species_list[0]
species_2 = good_species_list[1]
snp_divergences_1 = []
snp_divergences_2 = []

# Loop over sample pairs that are in both snp_divergence_map and gene_divergence_map
for sample_pair in (set(snp_divergence_map[species_1].keys()) & set(snp_divergence_map[species_2].keys()) ):
    
    snp_divergences_1.append( snp_divergence_map[species_1][sample_pair] )
    snp_divergences_2.append( snp_divergence_map[species_2][sample_pair] )

snp_divergences_1 = numpy.array(snp_divergences_1)
snp_divergences_2 = numpy.array(snp_divergences_2)

# Null expectation (medians line up)
median_ratio = numpy.median(snp_divergences_1)/numpy.median(snp_divergences_2)
species_axis_1.loglog([1e-06,1e-01],[1e-06*median_ratio,1e-01*median_ratio],'k-',linewidth=0.25)

# Observed values
species_axis_1.loglog(snp_divergences_2, snp_divergences_1, 'r.', markersize=2,alpha=0.5,markeredgewidth=0, rasterized=True)

# Then SNP divergence between other two species
species_1 = good_species_list[0]
species_2 = good_species_list[2]
snp_divergences_1 = []
snp_divergences_2 = []

# Loop over sample pairs that are in both snp_divergence_map and gene_divergence_map
for sample_pair in (set(snp_divergence_map[species_1].keys()) & set(snp_divergence_map[species_2].keys()) ):
    
    snp_divergences_1.append( snp_divergence_map[species_1][sample_pair] )
    snp_divergences_2.append( snp_divergence_map[species_2][sample_pair] )

snp_divergences_1 = numpy.array(snp_divergences_1)
snp_divergences_2 = numpy.array(snp_divergences_2)

# Null expectation (medians line up)
median_ratio = numpy.median(snp_divergences_1)/numpy.median(snp_divergences_2)
species_axis_2.loglog([1e-06,1e-01],[1e-06*median_ratio,1e-01*median_ratio],'k-',linewidth=0.25)

species_axis_2.loglog(snp_divergences_2, snp_divergences_1, 'r.', markersize=2,alpha=0.5,markeredgewidth=0,rasterized=True)

# Since y-axes are shared, do not duplicate ticklables
species_axis_1.set_yticklabels([])
species_axis_2.set_yticklabels([])


sys.stderr.write("Saving figure...\t")
fig.savefig('%s/supplemental_divergence_correlations.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight',dpi=600)
sys.stderr.write("Done!\n")

 