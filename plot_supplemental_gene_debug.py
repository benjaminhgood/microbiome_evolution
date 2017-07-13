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

debug=True
    
################################################################################

min_coverage = config.min_median_coverage
alpha = 0.5 # Confidence interval range for rate estimates
low_pi_threshold = 1e-03
clade_divergence_threshold = 1e-02
modification_divergence_threshold = 1e-03
min_change = 0.8
include_high_copynum = False
#include_high_copynum = True

# Load subject and sample metadata
sys.stderr.write("Loading sample metadata...\n")
subject_sample_map = parse_HMP_data.parse_subject_sample_map()
sample_country_map = parse_HMP_data.parse_sample_country_map()
sample_order_map = parse_HMP_data.parse_sample_order_map()
sys.stderr.write("Done!\n")
       
# Only plot samples above a certain depth threshold that are "haploids"
snp_samples = diversity_utils.calculate_haploid_samples(species_name, debug=debug)
# Only use the subset from North America 
# The only temporal samples are from here, best not contaminate the between-subject
# comparisons with out of sample effects
#snp_samples = snp_samples[parse_HMP_data.calculate_country_samples(sample_country_map, sample_list=snp_samples, allowed_countries=set(["United States"]))]

####################################################
#
# Set up Figure (4 panels, arranged in 2x2 grid)
#
####################################################

pylab.figure(1,figsize=(3.42,2))




# Load SNP information for species_name
sys.stderr.write("Loading SNPs for %s...\n" % species_name)
sys.stderr.write("(not just core genes...)\n")
pi_matrix_syn = numpy.array([])
avg_pi_matrix_syn = numpy.array([])
snp_difference_matrix = numpy.array([])
snp_difference_matrix_mutation = numpy.array([])
snp_difference_matrix_reversion = numpy.array([])
snp_opportunity_matrix = numpy.array([])
final_line_number = 0

while final_line_number >= 0:
    
    sys.stderr.write("Loading chunk starting @ %d...\n" % final_line_number)
    dummy_samples, allele_counts_map, passed_sites_map, final_line_number = parse_midas_data.parse_snps(species_name, debug=debug, allowed_samples=snp_samples,chunk_size=chunk_size,initial_line_number=final_line_number)
    sys.stderr.write("Done! Loaded %d genes\n" % len(allele_counts_map.keys()))
    
    print len(dummy_samples), "dummy samples!"
    
    # Calculate fixation matrix
    sys.stderr.write("Calculating matrix of snp differences...\n")
    chunk_snp_difference_matrix_mutation, chunk_snp_difference_matrix_reversion, chunk_snp_opportunity_matrix =     diversity_utils.calculate_fixation_matrix_mutation_reversion(allele_counts_map, passed_sites_map, min_change=min_change)   # 
    sys.stderr.write("Done!\n")
    
    if snp_difference_matrix.shape[0]==0:
        snp_difference_matrix = numpy.zeros_like(chunk_snp_difference_matrix_mutation)*1.0
        snp_difference_matrix_mutation = numpy.zeros_like(snp_difference_matrix)*1.0
        snp_difference_matrix_reversion = numpy.zeros_like(snp_difference_matrix)*1.0
        snp_opportunity_matrix = numpy.zeros_like(snp_difference_matrix)*1.0
    
    snp_difference_matrix += chunk_snp_difference_matrix_mutation+chunk_snp_difference_matrix_reversion
    snp_difference_matrix_mutation += chunk_snp_difference_matrix_mutation
    snp_difference_matrix_reversion += chunk_snp_difference_matrix_reversion
    
    
    snp_opportunity_matrix += chunk_snp_opportunity_matrix

    snp_samples = dummy_samples

snp_substitution_rate = snp_difference_matrix*1.0/(snp_opportunity_matrix+(snp_opportunity_matrix==0))
sys.stderr.write("Done!\n")   

# Load gene coverage information for species_name
sys.stderr.write("Loading pangenome data for %s...\n" % species_name)
gene_samples, gene_names, gene_presence_matrix, gene_depth_matrix, marker_coverages, gene_reads_matrix = parse_midas_data.parse_pangenome_data(species_name,allowed_samples=snp_samples)
sys.stderr.write("Done!\n")

sys.stderr.write("Loaded gene info for %d samples\n" % len(gene_samples))

gene_copynum_matrix = gene_depth_matrix*1.0/(marker_coverages+(marker_coverages==0))

prevalence_idxs = (parse_midas_data.calculate_unique_samples(subject_sample_map, gene_samples))*(marker_coverages>=min_coverage)
    
prevalences = gene_diversity_utils.calculate_fractional_gene_prevalences(gene_depth_matrix[:,prevalence_idxs], marker_coverages[prevalence_idxs])

pangenome_prevalences = numpy.array(prevalences,copy=True)
pangenome_prevalences.sort()



# Now need to make the gene samples and snp samples match up
desired_samples = gene_samples


num_haploids = len(desired_samples)
     
# Calculate which pairs of idxs belong to the same sample, which to the same subject
# and which to different subjects
#desired_same_sample_idxs, desired_same_subject_idxs, desired_diff_subject_idxs = parse_midas_data.calculate_subject_pairs(subject_sample_map, desired_samples)

# Calculate which pairs of idxs belong to the same sample, which to the same subject
# and which to different subjects
desired_same_sample_idxs, desired_same_subject_idxs, desired_diff_subject_idxs = parse_midas_data.calculate_ordered_subject_pairs(sample_order_map, desired_samples)

sys.stderr.write("%d temporal samples\n" % len(desired_same_subject_idxs[0]))

snp_sample_idx_map = parse_midas_data.calculate_sample_idx_map(desired_samples, snp_samples)
gene_sample_idx_map = parse_midas_data.calculate_sample_idx_map(desired_samples, gene_samples)
  

same_sample_snp_idxs  = parse_midas_data.apply_sample_index_map_to_indices(snp_sample_idx_map,  desired_same_sample_idxs)  
same_sample_gene_idxs = parse_midas_data.apply_sample_index_map_to_indices(gene_sample_idx_map, desired_same_sample_idxs)  


same_subject_snp_idxs = parse_midas_data.apply_sample_index_map_to_indices(snp_sample_idx_map, desired_same_subject_idxs)  
same_subject_gene_idxs = parse_midas_data.apply_sample_index_map_to_indices(gene_sample_idx_map, desired_same_subject_idxs)  

diff_subject_snp_idxs = parse_midas_data.apply_sample_index_map_to_indices(snp_sample_idx_map, desired_diff_subject_idxs)  
diff_subject_gene_idxs = parse_midas_data.apply_sample_index_map_to_indices(gene_sample_idx_map, desired_diff_subject_idxs)  


absent_thresholds = numpy.array([0, 0.01,0.05,0.1,0.25,0.49])

fold_changes = numpy.array([2,3,4,5,6,7,8,9,10,15,20])

total_changes = []
for absent_threshold in absent_thresholds:
    # Calculate matrix of number of genes that differ
    sys.stderr.write("Calculating matrix of gene differences for threshold %g...\n" % absent_threshold)
    gene_gain_matrix, gene_loss_matrix, gene_opportunity_matrix = gene_diversity_utils.calculate_coverage_based_gene_hamming_matrix_gain_loss(gene_depth_matrix, marker_coverages, absent_threshold=absent_threshold)
    

    gene_difference_matrix = gene_gain_matrix + gene_loss_matrix
    

    changes = 0

    for sample_pair_idx in xrange(0,len(same_subject_snp_idxs[0])):
        
        snp_i = same_subject_snp_idxs[0][sample_pair_idx]
        snp_j = same_subject_snp_idxs[1][sample_pair_idx]
    
             
        i = same_subject_gene_idxs[0][sample_pair_idx]
        j = same_subject_gene_idxs[1][sample_pair_idx]
    
        if marker_coverages[i]<min_coverage or marker_coverages[j]<min_coverage:
            # can't look at gene changes!
            continue
        
        else:
    
            if snp_substitution_rate[snp_i,snp_j] < modification_divergence_threshold:
                
                changes += gene_difference_matrix[i,j]
                
            else:
                continue
        
    total_changes.append(changes)
    
total_changes = numpy.array(total_changes)
total_changes = numpy.clip(total_changes, 3e-01,1e09)

pylab.semilogy(absent_thresholds, numpy.ones_like(absent_thresholds)*1,'k:')
pylab.semilogy(absent_thresholds, total_changes,'b.-')

pylab.xlabel('Copynum threshold for absence')
pylab.ylabel('Total number of gene modifications')

sys.stderr.write("Saving figure...\t")
pylab.savefig('%s/supplemental_gene_foldchange%s.pdf' % (parse_midas_data.analysis_directory, other_species_str),bbox_inches='tight')
pylab.savefig('%s/supplemental_gene_foldchange%s.png' % (parse_midas_data.analysis_directory, other_species_str),bbox_inches='tight',dpi=300)
sys.stderr.write("Done!\n")

    
