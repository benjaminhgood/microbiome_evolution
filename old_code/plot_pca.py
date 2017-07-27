import matplotlib  
matplotlib.use('Agg') 
import parse_midas_data
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


################################################################################
#
# Standard header to read in argument information
#
################################################################################
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("species_name", help="name of species to process")
parser.add_argument("--debug", help="Loads only a subset of SNPs for speed", action="store_true")
parser.add_argument("--chunk-size", type=int, help="max number of records to load", default=1000000000)
args = parser.parse_args()

species_name = args.species_name
debug = args.debug
chunk_size = args.chunk_size
################################################################################

min_change = 0.8
min_coverage = 20
alpha = 0.5 # Confidence interval range for rate estimates

# Set up figure
fig = plt.figure(figsize=(10, 4))

# Set up grids to hold figure panels
outer_grid = gridspec.GridSpec(1, 2, width_ratios=[1,1],wspace=0.1)

pca_axis = plt.Subplot(fig, outer_grid[0])
fig.add_subplot(pca_axis)
fig.suptitle(species_name)

low_divergence_pca_axis = plt.Subplot(fig, outer_grid[1])
fig.add_subplot(low_divergence_pca_axis)




# Load subject and sample metadata
sys.stderr.write("Loading HMP metadata...\n")
subject_sample_map = parse_midas_data.parse_subject_sample_map()
sys.stderr.write("Done!\n")
    
# Load genomic coverage distributions
sample_coverage_histograms, samples = parse_midas_data.parse_coverage_distribution(species_name)
median_coverages = numpy.array([stats_utils.calculate_nonzero_median_from_histogram(sample_coverage_histogram) for sample_coverage_histogram in sample_coverage_histograms])
sample_coverage_map = {samples[i]: median_coverages[i] for i in xrange(0,len(samples))}

# Load pi information for species_name
sys.stderr.write("Loading within-sample diversity for %s...\n" % species_name)
samples, total_pis, total_pi_opportunities = parse_midas_data.parse_within_sample_pi(species_name, debug)
sys.stderr.write("Done!\n")
pis = total_pis/total_pi_opportunities

median_coverages = numpy.array([sample_coverage_map[samples[i]] for i in xrange(0,len(samples))])

# Only plot samples above a certain depth threshold that are "haploids"
snp_samples = samples[(median_coverages>=min_coverage)*(pis<=1e-03)]
# Restrict to single timepoint single timepoints per person
unique_subject_idxs = parse_midas_data.calculate_unique_samples(subject_sample_map, snp_samples)
snp_samples = snp_samples[unique_subject_idxs]

# Analyze SNPs, looping over chunk sizes. 
# Clunky, but necessary to limit memory usage on cluster

# Load SNP information for species_name
sys.stderr.write("Loading %s...\n" % species_name)

genotype_matrix = numpy.array([])
passed_sites_matrix = numpy.array([])
snp_difference_matrix = numpy.array([])
snp_opportunity_matrix = numpy.array([])

final_line_number = 0
while final_line_number >= 0:
    
    sys.stderr.write("Loading chunk starting @ %d...\n" % final_line_number)
    dummy_samples, allele_counts_map, passed_sites_map, final_line_number = parse_midas_data.parse_snps(species_name, debug=debug, allowed_samples=snp_samples,chunk_size=chunk_size,initial_line_number=final_line_number)
    sys.stderr.write("Done! Loaded %d genes\n" % len(allele_counts_map.keys()))
    
    # Calculate fixation matrix
    sys.stderr.write("Calculating matrix of snp differences...\n")
    chunk_snp_difference_matrix, chunk_snp_opportunity_matrix = diversity_utils.calculate_fixation_matrix(allele_counts_map, passed_sites_map, min_change=min_change)    
    sys.stderr.write("Done!\n")
    
    if snp_difference_matrix.shape[0]==0:
        snp_difference_matrix = numpy.zeros_like(chunk_snp_difference_matrix)*1.0
        snp_opportunity_matrix = numpy.zeros_like(snp_difference_matrix)*1.0
    
    snp_difference_matrix += chunk_snp_difference_matrix
    snp_opportunity_matrix += chunk_snp_opportunity_matrix
    
    sys.stderr.write("Calculating genotype matrix for chunk...\n")
    for gene_name in allele_counts_map.keys():
        for variant_type in allele_counts_map[gene_name].keys():
    
            if len(allele_counts_map[gene_name][variant_type]['alleles']) ==0:
                continue
            
            chunk_genotype_matrix, chunk_passed_sites_matrix = diversity_utils.calculate_consensus_genotypes(allele_counts_map[gene_name][variant_type]['alleles'])
    
            if chunk_genotype_matrix.shape[0]<1:
                continue
    
            if genotype_matrix.shape[0]==0:
                genotype_matrix = numpy.zeros((0,len(dummy_samples)))
                passed_sites_matrix = numpy.zeros((0,len(dummy_samples)))
                
        
            genotype_matrix = numpy.vstack([genotype_matrix, chunk_genotype_matrix])
            passed_sites_matrix = numpy.vstack([passed_sites_matrix, chunk_passed_sites_matrix])
    sys.stderr.write("Done!\n")
    
substitution_rate = snp_difference_matrix*1.0/snp_opportunity_matrix 

sys.stderr.write("Calculating PCA coordinates...\n")
pca_coords, variances = diversity_utils.calculate_pca_coordinates(genotype_matrix, passed_sites_matrix)
sys.stderr.write("Done!\n")

pca_axis.set_ylabel('PC2 (%0.1f)' % (100*variances[1]))
pca_axis.set_xlabel('PC1 (%0.1f)' % (100*variances[0]))



# Calculate which pairs of idxs belong to the same sample, which to the same subject
# and which to different subjects
same_sample_idxs, same_subject_idxs, diff_subject_idxs = parse_midas_data.calculate_subject_pairs(subject_sample_map, snp_samples)

# Find pairs of samples from different hosts that have very low divergence
low_divergence_diff_sample_idx_idxs = (substitution_rate[diff_subject_idxs]<1e-04)
low_divergence_sample_idxs = list(set(diff_subject_idxs[0][low_divergence_diff_sample_idx_idxs]) |  set(diff_subject_idxs[1][low_divergence_diff_sample_idx_idxs]))
low_divergence_samples = snp_samples[low_divergence_sample_idxs]

pca_axis.plot(pca_coords[0], pca_coords[1],'b.',alpha=0.5,markersize=3)
pca_axis.plot(pca_coords[0][low_divergence_sample_idxs], pca_coords[1][low_divergence_sample_idxs],'r.',alpha=0.5)


if len(low_divergence_samples)>2:
    # Now redo everything again with low_divergence_samples

    sys.stderr.write("Continuing with %d low-divergence samples...\n" % len(low_divergence_samples))

    # Load SNP information for species_name
    sys.stderr.write("Re-loading %s...\n" % species_name)

    genotype_matrix = numpy.array([])
    passed_sites_matrix = numpy.array([])
    snp_difference_matrix = numpy.array([])
    snp_opportunity_matrix = numpy.array([])

    final_line_number = 0
    while final_line_number >= 0:
    
        sys.stderr.write("Loading chunk starting @ %d...\n" % final_line_number)
        dummy_samples, allele_counts_map, passed_sites_map, final_line_number = parse_midas_data.parse_snps(species_name, debug=debug, allowed_samples=low_divergence_samples,chunk_size=chunk_size,initial_line_number=final_line_number)
        sys.stderr.write("Done! Loaded %d genes\n" % len(allele_counts_map.keys()))
    
        # Calculate fixation matrix
        sys.stderr.write("Calculating matrix of snp differences...\n")
        chunk_snp_difference_matrix, chunk_snp_opportunity_matrix = diversity_utils.calculate_fixation_matrix(allele_counts_map, passed_sites_map, min_change=min_change)    
        sys.stderr.write("Done!\n")
    
        if snp_difference_matrix.shape[0]==0:
            snp_difference_matrix = numpy.zeros_like(chunk_snp_difference_matrix)*1.0
            snp_opportunity_matrix = numpy.zeros_like(snp_difference_matrix)*1.0
    
        snp_difference_matrix += chunk_snp_difference_matrix
        snp_opportunity_matrix += chunk_snp_opportunity_matrix
    
        sys.stderr.write("Calculating genotype matrix for chunk...\n")
        for gene_name in allele_counts_map.keys():
            for variant_type in allele_counts_map[gene_name].keys():
    
                if len(allele_counts_map[gene_name][variant_type]['alleles']) ==0:
                    continue
            
                chunk_genotype_matrix, chunk_passed_sites_matrix = diversity_utils.calculate_consensus_genotypes(allele_counts_map[gene_name][variant_type]['alleles'])
    
                if chunk_genotype_matrix.shape[0]<1:
                    continue
    
                if genotype_matrix.shape[0]==0:
                    genotype_matrix = numpy.zeros((0,len(dummy_samples)))
                    passed_sites_matrix = numpy.zeros((0,len(dummy_samples)))
                
        
                genotype_matrix = numpy.vstack([genotype_matrix, chunk_genotype_matrix])
                passed_sites_matrix = numpy.vstack([passed_sites_matrix, chunk_passed_sites_matrix])
        sys.stderr.write("Done!\n")
    
    substitution_rate = snp_difference_matrix*1.0/snp_opportunity_matrix 

    sys.stderr.write("Calculating PCA coordinates...\n")
    pca_coords, variances = diversity_utils.calculate_pca_coordinates(genotype_matrix, passed_sites_matrix)
    sys.stderr.write("Done!")

    low_divergence_pca_axis.set_ylabel('low-divergence PC2 (%0.1f)' % (100*variances[1]))
    low_divergence_pca_axis.set_xlabel('low-divergence PC1 (%0.1f)' % (100*variances[0]))
    low_divergence_pca_axis.yaxis.set_label_position("right")
    
    low_divergence_pca_axis.plot(pca_coords[0], pca_coords[1],'r.',alpha=0.5)

pca_axis.set_xticklabels([])
pca_axis.set_yticklabels([])
low_divergence_pca_axis.set_xticklabels([])
low_divergence_pca_axis.set_yticklabels([])

fig.savefig('%s/%s_pca.pdf' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight')
fig.savefig('%s/%s_pca.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight',dpi=300)

    
