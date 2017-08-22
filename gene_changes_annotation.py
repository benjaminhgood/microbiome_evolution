# Within-snp gene changes

import matplotlib  
matplotlib.use('Agg') 
import config
import parse_midas_data
import os
import parse_HMP_data


import pylab
import sys
import numpy

import diversity_utils
import gene_diversity_utils
import sfs_utils
import calculate_substitution_rates
import calculate_temporal_changes

import stats_utils
import matplotlib.colors as colors
import matplotlib.cm as cmx
from math import log10,ceil
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy.random import randint, choice

mpl.rcParams['font.size'] = 5
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
parser.add_argument('--other-species', type=str, help='Run the script for a different species')

args = parser.parse_args()

debug = args.debug
chunk_size = args.chunk_size
other_species = args.other_species

if other_species:
    species_name = other_species
    other_species_str = "_%s" % species_name
else:
    other_species_str = ""
    
################################################################################

# parameters:

min_coverage = config.min_median_coverage
clade_divergence_threshold = 1e-02
modification_divergence_threshold = 1e-03 #the threshold for deciding when something is a modification vs a replacement. Like most other things, it is an arbitrary choice. 

# Load subject and sample metadata
sys.stderr.write("Loading sample metadata...\n")
subject_sample_map = parse_HMP_data.parse_subject_sample_map()
sample_country_map = parse_HMP_data.parse_sample_country_map()
sample_order_map = parse_HMP_data.parse_sample_order_map()
sys.stderr.write("Done!\n")
       
# Only plot samples above a certain depth threshold that are "haploids"
snp_samples = diversity_utils.calculate_haploid_samples(species_name, debug=debug)



####################
# Analyze the data #
####################

reference_genes = parse_midas_data.load_reference_genes(species_name)

# Analyze SNPs, looping over chunk sizes. 
# Clunky, but necessary to limit memory usage on cluster

import sfs_utils
sys.stderr.write("Loading SFSs for %s...\t" % species_name)
samples, sfs_map = parse_midas_data.parse_within_sample_sfs(species_name, allowed_variant_types=set(['1D','2D','3D','4D'])) 
sys.stderr.write("Done!\n")


sys.stderr.write("Loading pre-computed substitution rates for %s...\n" % species_name)
substitution_rate_map = calculate_substitution_rates.load_substitution_rate_map(species_name)
sys.stderr.write("Calculating matrix...\n")
dummy_samples, snp_difference_matrix, snp_opportunity_matrix = calculate_substitution_rates.calculate_matrices_from_substitution_rate_map(substitution_rate_map, 'all', allowed_samples=snp_samples)
snp_samples = dummy_samples
sys.stderr.write("Done!\n")

sys.stderr.write("Loading pre-computed temporal changes for %s...\n" % species_name)
temporal_change_map = calculate_temporal_changes.load_temporal_change_map(species_name)
sys.stderr.write("Done!\n")

snp_substitution_rate = snp_difference_matrix*1.0/(snp_opportunity_matrix+(snp_opportunity_matrix==0))
sys.stderr.write("Done!\n")   

# Load gene coverage information for species_name
sys.stderr.write("Loading pangenome data for %s...\n" % species_name)
gene_samples, gene_names, gene_presence_matrix, gene_depth_matrix, marker_coverages, gene_reads_matrix = parse_midas_data.parse_pangenome_data(species_name,allowed_samples=snp_samples)
sys.stderr.write("Done!\n")

sys.stderr.write("Loaded gene info for %d samples\n" % len(gene_samples))

gene_copynum_matrix = gene_depth_matrix*1.0/(marker_coverages+(marker_coverages==0))

clipped_gene_copynum_matrix = numpy.clip(gene_depth_matrix,0.1,1e09)/(marker_coverages+0.1*(marker_coverages==0))

low_copynum_matrix = (gene_copynum_matrix<=3)
good_copynum_matrix = (gene_copynum_matrix>=0.5)*(gene_copynum_matrix<=3) # why isn't this till 2? NRG

prevalence_idxs = (parse_midas_data.calculate_unique_samples(subject_sample_map, gene_samples))*(marker_coverages>=min_coverage)
    
prevalences = gene_diversity_utils.calculate_fractional_gene_prevalences(gene_depth_matrix[:,prevalence_idxs], marker_coverages[prevalence_idxs])

pangenome_prevalences = numpy.array(prevalences,copy=True)
pangenome_prevalences.sort()
        
# Calculate matrix of number of genes that differ
sys.stderr.write("Calculating matrix of gene differences...\n")
gene_gain_matrix, gene_loss_matrix, gene_opportunity_matrix = gene_diversity_utils.calculate_coverage_based_gene_hamming_matrix_gain_loss(gene_reads_matrix, gene_depth_matrix, marker_coverages)

gene_difference_matrix = gene_gain_matrix + gene_loss_matrix

# Now need to make the gene samples and snp samples match up
desired_samples = gene_samples

num_haploids = len(desired_samples)
     
# Calculate which pairs of idxs belong to the same sample, which to the same subject
# and which to different subjects
desired_same_sample_idxs, desired_same_subject_idxs, desired_diff_subject_idxs = parse_midas_data.calculate_ordered_subject_pairs(sample_order_map, desired_samples)

sys.stderr.write("%d temporal samples\n" % len(desired_same_subject_idxs[0]))

# get the idx for samples being considered for snp changes to match with gene changes.

snp_sample_idx_map = parse_midas_data.calculate_sample_idx_map(desired_samples, snp_samples)
gene_sample_idx_map = parse_midas_data.calculate_sample_idx_map(desired_samples, gene_samples)

# indexes for time pairs  
same_subject_snp_idxs = parse_midas_data.apply_sample_index_map_to_indices(snp_sample_idx_map, desired_same_subject_idxs)  
same_subject_gene_idxs = parse_midas_data.apply_sample_index_map_to_indices(gene_sample_idx_map, desired_same_subject_idxs)  

# Calculate subset of "modification timepoints" 
modification_pair_idxs = set([])

for sample_pair_idx in xrange(0,len(same_subject_snp_idxs[0])):
#    
    snp_i = same_subject_snp_idxs[0][sample_pair_idx]
    snp_j = same_subject_snp_idxs[1][sample_pair_idx]

    if snp_substitution_rate[snp_i, snp_j] < modification_divergence_threshold:
        modification_pair_idxs.add( sample_pair_idx ) 


# In this loop, pull out the gene idxes for genes that are changing. Use this to look up the gene IDs. 

for sample_pair_idx in xrange(0,len(same_subject_snp_idxs[0])):    
    snp_i = same_subject_snp_idxs[0][sample_pair_idx]
    snp_j = same_subject_snp_idxs[1][sample_pair_idx]
    #
    sample_i = snp_samples[snp_i]
    sample_j = snp_samples[snp_j]
    if marker_coverages[i]>min_coverage and marker_coverages[j]>min_coverage:
        # iterate through modificaiton pair idxs -- otherwise we get confused with replacements. 
        if sample_pair_idx in modification_pair_idxs:
            # Calculate set of genes that are present in at least one sample
            present_gene_idxs = []
            present_gene_idxs.extend( numpy.nonzero( (gene_copynum_matrix[:,i]>0.5)*(gene_copynum_matrix[:,i]<2))[0] )
            present_gene_idxs.extend( numpy.nonzero( (gene_copynum_matrix[:,j]>0.5)*(gene_copynum_matrix[:,j]<2))[0] )
            #
            # calculate the indexes of genes that are changing. 
            pair_specific_gene_idxs = gene_diversity_utils.calculate_gene_differences_between_idxs(i, j, gene_reads_matrix, gene_depth_matrix, marker_coverages)
            if len(pair_specific_gene_idxs)>0:
                print pair_specific_gene_idxs
