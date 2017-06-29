# SFS plots, comparing syn and non syn 
# (a) regular SFS (percentage of population)
# (b) Singletons, doubletons, and pi-weighted

# Summary across species?
# dN/dS vs dS plot? 
# core vs variable genes? 

# Within-snp gene changes

import matplotlib  
matplotlib.use('Agg') 
import config
import parse_midas_data
###
#
# For today while the new data processes
#
import os
#parse_midas_data.data_directory = os.path.expanduser("~/ben_nandita_hmp_data_062517/")
#########################################
import parse_HMP_data


import pylab
import sys
import numpy

import diversity_utils
import gene_diversity_utils

import stats_utils
import matplotlib.colors as colors
import matplotlib.cm as cmx
from math import log10,ceil
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy.random import randint

mpl.rcParams['font.size'] = 6
mpl.rcParams['lines.linewidth'] = 0.5
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'

species_name = "Bacteroides_vulgatus_57955"

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

min_coverage = config.min_median_coverage
alpha = 0.5 # Confidence interval range for rate estimates
low_divergence_threshold = 5e-04
min_change = 0.8
allowed_variant_types = set(['1D','2D','3D','4D'])
max_clade_d = 1e-02

# Load subject and sample metadata
sys.stderr.write("Loading sample metadata...\n")
subject_sample_map = parse_HMP_data.parse_subject_sample_map()
sample_country_map = parse_HMP_data.parse_sample_country_map()
sys.stderr.write("Done!\n")
       
# Only plot samples above a certain depth threshold that are "haploids"
snp_samples = diversity_utils.calculate_haploid_samples(species_name, debug=debug)
snp_samples = snp_samples[ diversity_utils.parse_midas_data.calculate_unique_samples(subject_sample_map, snp_samples)]


####################################################
#
# Set up Figure (2 panels, arranged horizontally)
#
####################################################
# This figure spreads them all out

pylab.figure(1,figsize=(5,2))
fig = pylab.gcf()
# make three panels panels
outer_grid  = gridspec.GridSpec(1,2, width_ratios=[2, 1], wspace=0.25)

##############################################################################
#
# Panel (a). Pooled SFS as function of minor allele freq
#
##############################################################################

sfs_axis = plt.Subplot(fig, outer_grid[0])
fig.add_subplot(sfs_axis)

sfs_axis.set_xlabel('Minor allele freq in largest clade, $f$')
sfs_axis.set_ylabel('Scaled fraction of sites, $f(1-f) P(f)$')
sfs_axis.set_xlim([0,0.5])
sfs_axis.spines['top'].set_visible(False)
sfs_axis.spines['right'].set_visible(False)
sfs_axis.get_xaxis().tick_bottom()
sfs_axis.get_yaxis().tick_left()

##############################################################################
#
# Panel (b). Comparison of 1D and 4D for minor allele *count* 
#
##############################################################################

count_axis = plt.Subplot(fig, outer_grid[1])
fig.add_subplot(count_axis)

count_axis.set_ylabel('Fraction nonsynonymous (%)')
count_axis.set_xlim([-0.5,3.5])
count_axis.set_ylim([0,100])
count_axis.set_xticks([0,1,2,3])
count_axis.set_xticklabels(['1','2','3','4+'])
count_axis.spines['top'].set_visible(False)
count_axis.spines['right'].set_visible(False)
count_axis.get_xaxis().tick_bottom()
count_axis.get_yaxis().tick_left()

###############
#
# Do calculation!
#
###############

sys.stderr.write("Loading core genes...\n")
core_genes = parse_midas_data.load_core_genes(species_name)
sys.stderr.write("Done! Core genome consists of %d genes\n" % len(core_genes))


# Analyze SNPs, looping over chunk sizes. 
# Clunky, but necessary to limit memory usage on cluster

# Load SNP information for species_name
sys.stderr.write("Loading %d samples for %s...\n" % (len(snp_samples), species_name))

snp_difference_matrix = numpy.array([])
snp_opportunity_matrix = numpy.array([])

final_line_number = 0
while final_line_number >= 0:
    
    sys.stderr.write("Loading chunk starting @ %d...\n" % final_line_number)
    dummy_samples, allele_counts_map, passed_sites_map, final_line_number = parse_midas_data.parse_snps(species_name, debug=debug, allowed_variant_types=allowed_variant_types, allowed_genes=core_genes, allowed_samples=snp_samples,chunk_size=chunk_size,initial_line_number=final_line_number)
    sys.stderr.write("Done! Loaded %d genes\n" % len(allele_counts_map.keys()))
    
    # Calculate fixation matrix
    sys.stderr.write("Calculating matrix of snp differences...\n")
    chunk_snp_difference_matrix, chunk_snp_opportunity_matrix = diversity_utils.calculate_fixation_matrix(allele_counts_map, passed_sites_map, allowed_variant_types=allowed_variant_types, allowed_genes=core_genes, min_change=min_change)    
    sys.stderr.write("Done!\n")
    
    if snp_difference_matrix.shape[0]==0:
        snp_difference_matrix = numpy.zeros_like(chunk_snp_difference_matrix)*1.0
        snp_opportunity_matrix = numpy.zeros_like(snp_difference_matrix)*1.0
    
    snp_difference_matrix += chunk_snp_difference_matrix
    snp_opportunity_matrix += chunk_snp_opportunity_matrix
    
     
substitution_rate = snp_difference_matrix*1.0/snp_opportunity_matrix 

coarse_grained_idxs, coarse_grained_cluster_list = diversity_utils.cluster_samples(substitution_rate, min_d=low_divergence_threshold, max_ds=[max_clade_d])
coarse_grained_cluster_idxss = coarse_grained_cluster_list[0] 
coarse_grained_cluster_sizes = numpy.array([cluster_idxs.sum() for cluster_idxs in coarse_grained_cluster_idxss])

print "Top level:", len(coarse_grained_cluster_idxss), coarse_grained_cluster_sizes

# only focus on the members of the largest clade
remapped_cluster_idxss = [cluster_idxs[coarse_grained_idxs] for cluster_idxs in coarse_grained_cluster_idxss]
largest_clade_idxs = remapped_cluster_idxss[0]
largest_clade_size = largest_clade_idxs.sum()

coarse_grained_samples = snp_samples[coarse_grained_idxs]
largest_clade_samples = set(coarse_grained_samples[largest_clade_idxs])

if len(coarse_grained_samples)<3:
    sys.stderr.write("Too few samples! Quitting...\n")
    sys.exit(1)
    
sys.stderr.write("Continuing with %d samples...\n" % len(coarse_grained_samples))

# Load SNP information for species_name
sys.stderr.write("Re-loading %s...\n" % species_name)

snp_difference_matrix = numpy.array([])
snp_opportunity_matrix = numpy.array([])
    
   
maf_bins = numpy.arange(1,largest_clade_size+1)*1.0/largest_clade_size
maf_bins -= (maf_bins[1]-maf_bins[0])/2
maf_bins[0]=-0.1
maf_bins[-1] = 1.1
mafs = numpy.arange(1,largest_clade_size)*1.0/largest_clade_size
    
count_bins = numpy.arange(1,largest_clade_size+1)-0.5
count_locations = numpy.arange(1,largest_clade_size)    
    
synonymous_sfs = numpy.zeros_like(mafs)
nonsynonymous_sfs = numpy.zeros_like(mafs)

synonymous_count_sfs = numpy.zeros_like(count_locations)
nonsynonymous_count_sfs = numpy.zeros_like(count_locations)

synonymous_pi_weighted_counts = 0
nonsynonymous_pi_weighted_counts = 0

    
final_line_number = 0
while final_line_number >= 0:
    
    sys.stderr.write("Loading chunk starting @ %d...\n" % final_line_number)
    dummy_samples, allele_counts_map, passed_sites_map, final_line_number = parse_midas_data.parse_snps(species_name, debug=debug, allowed_variant_types=allowed_variant_types, allowed_samples=coarse_grained_samples,allowed_genes=core_genes, chunk_size=chunk_size,initial_line_number=final_line_number)
    sys.stderr.write("Done! Loaded %d genes\n" % len(allele_counts_map.keys()))
    
    # Calculate fixation matrix
    sys.stderr.write("Calculating matrix of snp differences...\n")
    chunk_snp_difference_matrix, chunk_snp_opportunity_matrix = diversity_utils.calculate_fixation_matrix(allele_counts_map, passed_sites_map, allowed_genes=core_genes, min_change=min_change)    
    sys.stderr.write("Done!\n")
    
    if snp_difference_matrix.shape[0]==0:
        snp_difference_matrix = numpy.zeros_like(chunk_snp_difference_matrix)*1.0
        snp_opportunity_matrix = numpy.zeros_like(snp_difference_matrix)*1.0
    
    snp_difference_matrix += chunk_snp_difference_matrix
    snp_opportunity_matrix += chunk_snp_opportunity_matrix
    
        
    sys.stderr.write("Calculating the SFS...\n")
    chunk_synonymous_freqs = diversity_utils.calculate_pooled_freqs(allele_counts_map, passed_sites_map, allowed_sample_idxs=largest_clade_idxs, allowed_variant_types = set(['4D']), allowed_genes=core_genes)
    chunk_nonsynonymous_freqs = diversity_utils.calculate_pooled_freqs(allele_counts_map, passed_sites_map, allowed_sample_idxs=largest_clade_idxs, allowed_variant_types = set(['1D']), allowed_genes=core_genes)
        
    chunk_synonymous_sfs, dummy = numpy.histogram(chunk_synonymous_freqs, bins=maf_bins) 
    synonymous_sfs += chunk_synonymous_sfs
        
    chunk_nonsynonymous_sfs, dummy = numpy.histogram(chunk_nonsynonymous_freqs, bins=maf_bins) 
    nonsynonymous_sfs += chunk_nonsynonymous_sfs
    
    sys.stderr.write("Calculating count SFS...\n")
    chunk_synonymous_counts, chunk_synonymous_weights = diversity_utils.calculate_pooled_counts(allele_counts_map, passed_sites_map, allowed_sample_idxs=largest_clade_idxs, allowed_variant_types = set(['4D']), allowed_genes=core_genes,pi_min_k=4)
    chunk_nonsynonymous_counts, chunk_nonsynonymous_weights = diversity_utils.calculate_pooled_counts(allele_counts_map, passed_sites_map, allowed_sample_idxs=largest_clade_idxs, allowed_variant_types = set(['1D']), allowed_genes=core_genes,pi_min_k=4)
        
    chunk_synonymous_count_sfs, dummy = numpy.histogram(chunk_synonymous_counts, bins=count_bins) 
    synonymous_count_sfs += chunk_synonymous_count_sfs
        
    chunk_nonsynonymous_count_sfs, dummy = numpy.histogram(chunk_nonsynonymous_counts, bins=count_bins) 
    nonsynonymous_count_sfs += chunk_nonsynonymous_count_sfs

    synonymous_pi_weighted_counts += chunk_synonymous_weights
    nonsynonymous_pi_weighted_counts += chunk_nonsynonymous_weights

sys.stderr.write("%d singletons!\n" % (synonymous_count_sfs[0]+nonsynonymous_count_sfs[0]))
sys.stderr.write("%d doubletons!\n" % (synonymous_count_sfs[1]+nonsynonymous_count_sfs[1]))
sys.stderr.write("%d tripletons!\n" % (synonymous_count_sfs[2]+nonsynonymous_count_sfs[2]))
sys.stderr.write("%d rest-tons!\n" % (synonymous_count_sfs[3:]+nonsynonymous_count_sfs[3:]).sum())

sys.stderr.write("%g pi-weighted!\n" % (nonsynonymous_pi_weighted_counts+synonymous_pi_weighted_counts))


print count_locations
print synonymous_count_sfs
print nonsynonymous_count_sfs

###############
#
# Plot figures!
#
###############

sfs_axis.plot(mafs, synonymous_sfs*mafs*(1-mafs)/(synonymous_sfs*mafs*(1-mafs)).sum(), 'b.-',label='4D')
sfs_axis.plot(mafs, nonsynonymous_sfs*mafs*(1-mafs)/(nonsynonymous_sfs*mafs*(1-mafs)).sum(),'r.-',label='1D')

singleton_fraction = nonsynonymous_count_sfs[0]*1.0/(synonymous_count_sfs[0]+nonsynonymous_count_sfs[0])
doubleton_fraction =  nonsynonymous_count_sfs[1]*1.0/(synonymous_count_sfs[1]+nonsynonymous_count_sfs[1])
tripleton_fraction =  nonsynonymous_count_sfs[2]*1.0/(synonymous_count_sfs[2]+nonsynonymous_count_sfs[2])
restton_fraction =   nonsynonymous_count_sfs[3:].sum()*1.0/((synonymous_count_sfs[3:]+nonsynonymous_count_sfs[3:]).sum())


pi_weighted_fraction = nonsynonymous_pi_weighted_counts*1.0/(nonsynonymous_pi_weighted_counts+synonymous_pi_weighted_counts)

count_axis.bar([-0.3],[singleton_fraction*100],width=0.6,color='r',linewidth=0)
count_axis.bar([1-0.3],[doubleton_fraction*100],width=0.6,color='r',linewidth=0)
count_axis.bar([2-0.3],[tripleton_fraction*100],width=0.6,color='r',linewidth=0)
#count_axis.bar([3-0.3],[restton_fraction*100],width=0.6,color='r',linewidth=0)
count_axis.bar([3-0.3],[pi_weighted_fraction*100],width=0.6,color='r',linewidth=0)


sys.stderr.write("Saving figure...\t")
fig.savefig('%s/figure_5.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')
sys.stderr.write("Done!\n")
