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
parser.add_argument("species_name_1", help="name of first species to process")
parser.add_argument("species_name_2", help="name of first species to process")
parser.add_argument("--debug", help="Loads only a subset of SNPs for speed", action="store_true")
parser.add_argument("--chunk-size", type=int, help="max number of records to load", default=1000000000)
args = parser.parse_args()

species_name_1 = args.species_name_1
species_name_2 = args.species_name_2
debug = args.debug
chunk_size = args.chunk_size
################################################################################

min_change = 0.8
min_coverage = 20
alpha = 0.5 # Confidence interval range for rate estimates


# Load subject and sample metadata
sys.stderr.write("Loading HMP metadata...\n")
subject_sample_map = parse_midas_data.parse_subject_sample_map()
sys.stderr.write("Done!\n")

sys.stderr.write("First processing %s...\n" % species_name_1)   
    
# Load genomic coverage distributions
sample_coverage_histograms, samples = parse_midas_data.parse_coverage_distribution(species_name_1)
median_coverages = numpy.array([stats_utils.calculate_nonzero_median_from_histogram(sample_coverage_histogram) for sample_coverage_histogram in sample_coverage_histograms])
sample_coverage_map = {samples[i]: median_coverages[i] for i in xrange(0,len(samples))}

# Load pi information for first species
sys.stderr.write("Loading within-sample diversity for %s...\n" % species_name_1)
samples, total_pis, total_pi_opportunities = parse_midas_data.parse_within_sample_pi(species_name_1, debug)
sys.stderr.write("Done!\n")
pis = total_pis/total_pi_opportunities

median_coverages = numpy.array([sample_coverage_map[samples[i]] for i in xrange(0,len(samples))])

# Only plot samples above a certain depth threshold that are "haploids"
first_snp_samples = samples[(median_coverages>=min_coverage)*(pis<=1e-03)]

# Analyze SNPs, looping over chunk sizes. 
# Clunky, but necessary to limit memory usage on cluster

# Load SNP information for species_name
sys.stderr.write("Loading %s...\n" % species_name_1)

first_snp_difference_matrix = numpy.array([])
first_snp_opportunity_matrix = numpy.array([])
final_line_number = 0
while final_line_number >= 0:
    
    sys.stderr.write("Loading chunk starting @ %d...\n" % final_line_number)
    dummy_samples, allele_counts_map, passed_sites_map, final_line_number = parse_midas_data.parse_snps(species_name_1, debug=debug, allowed_samples=first_snp_samples,chunk_size=chunk_size,initial_line_number=final_line_number)
    sys.stderr.write("Done! Loaded %d genes\n" % len(allele_counts_map.keys()))
    
    # Calculate fixation matrix
    sys.stderr.write("Calculating matrix of snp differences...\n")
    chunk_snp_difference_matrix, chunk_snp_opportunity_matrix =     diversity_utils.calculate_fixation_matrix(allele_counts_map, passed_sites_map, min_change=min_change)    
    sys.stderr.write("Done!\n")
    
    if first_snp_difference_matrix.shape[0]==0:
        first_snp_difference_matrix = numpy.zeros_like(chunk_snp_difference_matrix)*1.0
        first_snp_opportunity_matrix = numpy.zeros_like(first_snp_difference_matrix)*1.0
    
    first_snp_difference_matrix += chunk_snp_difference_matrix
    first_snp_opportunity_matrix += chunk_snp_opportunity_matrix

sys.stderr.write("Done!\n")   

sys.stderr.write("Now processing %s...\n" % species_name_2)

# Load genomic coverage distributions for second species
sample_coverage_histograms, samples = parse_midas_data.parse_coverage_distribution(species_name_2)
median_coverages = numpy.array([stats_utils.calculate_nonzero_median_from_histogram(sample_coverage_histogram) for sample_coverage_histogram in sample_coverage_histograms])
sample_coverage_map = {samples[i]: median_coverages[i] for i in xrange(0,len(samples))}

# Load pi information for second species
sys.stderr.write("Loading within-sample diversity for %s...\n" % species_name_2)
samples, total_pis, total_pi_opportunities = parse_midas_data.parse_within_sample_pi(species_name_2, debug)
sys.stderr.write("Done!\n")
pis = total_pis/total_pi_opportunities

median_coverages = numpy.array([sample_coverage_map[samples[i]] for i in xrange(0,len(samples))])

# Only plot samples above a certain depth threshold that are "haploids"
second_snp_samples = set(samples[(median_coverages>=min_coverage)*(pis<=1e-03)])

joint_sample_idxs = []
for i in xrange(0,len(first_snp_samples)):
    if first_snp_samples[i] in second_snp_samples:
        joint_sample_idxs.append(i)
joint_sample_idxs = numpy.array(joint_sample_idxs)        
joint_snp_samples = first_snp_samples[joint_sample_idxs]

first_snp_difference_matrix = first_snp_difference_matrix[numpy.ix_(joint_sample_idxs,joint_sample_idxs)]

first_snp_opportunity_matrix = first_snp_opportunity_matrix[numpy.ix_(joint_sample_idxs,joint_sample_idxs)]


# Analyze SNPs, looping over chunk sizes. 
# Clunky, but necessary to limit memory usage on cluster

# Load SNP information for species_name
sys.stderr.write("Loading %s...\n" % species_name_2)

second_snp_difference_matrix = numpy.array([])
second_snp_opportunity_matrix = numpy.array([])
final_line_number = 0
while final_line_number >= 0:
    
    sys.stderr.write("Loading chunk starting @ %d...\n" % final_line_number)
    dummy_samples, allele_counts_map, passed_sites_map, final_line_number = parse_midas_data.parse_snps(species_name_2, debug=debug, allowed_samples=joint_snp_samples,chunk_size=chunk_size,initial_line_number=final_line_number)
    sys.stderr.write("Done! Loaded %d genes\n" % len(allele_counts_map.keys()))
    
    # Calculate fixation matrix
    sys.stderr.write("Calculating matrix of snp differences...\n")
    chunk_snp_difference_matrix, chunk_snp_opportunity_matrix =     diversity_utils.calculate_fixation_matrix(allele_counts_map, passed_sites_map, min_change=min_change)    
    sys.stderr.write("Done!\n")
    
    if second_snp_difference_matrix.shape[0]==0:
        second_snp_difference_matrix = numpy.zeros_like(chunk_snp_difference_matrix)*1.0
        second_snp_opportunity_matrix = numpy.zeros_like(second_snp_difference_matrix)*1.0
    
    second_snp_difference_matrix += chunk_snp_difference_matrix
    second_snp_opportunity_matrix += chunk_snp_opportunity_matrix

sys.stderr.write("Done!\n")   

# Now calculate rates
first_divergence_matrix = first_snp_difference_matrix /  first_snp_opportunity_matrix
second_divergence_matrix = second_snp_difference_matrix /  second_snp_opportunity_matrix
     
# Calculate which pairs of idxs belong to the same sample, which to the same subject
# and which to different subjects
same_sample_idxs, same_subject_idxs, diff_subject_idxs = parse_midas_data.calculate_subject_pairs(subject_sample_map, joint_snp_samples)

# Set up figure
fig = plt.figure(figsize=(3, 2.5))

# Set up grids to hold figure panels
outer_grid = gridspec.GridSpec(1, 1)

###################
#
# SNP Panel
#
###################

divergence_axis = plt.Subplot(fig, outer_grid[0])
fig.add_subplot(divergence_axis)

divergence_axis.set_ylabel('%s divergence' % species_name_2)
divergence_axis.set_xlabel('%s divergence' % species_name_1)
divergence_axis.set_xlim([1e-07,9e-02])
divergence_axis.set_ylim([1e-07,9e-02])

divergence_axis.loglog([1e-09,1e-09],[1,1],'g.',label='Within host')
divergence_axis.loglog([1e-09,1e-09],[1,1],'r.',label='Between host')

pylab.plot(first_divergence_matrix[same_subject_idxs], second_divergence_matrix[same_subject_idxs],'g.',alpha=0.5,markersize=3)
pylab.plot(first_divergence_matrix[diff_subject_idxs], second_divergence_matrix[diff_subject_idxs],'r.',alpha=0.5,markersize=3)


divergence_axis.legend(loc='upper right',frameon=False)

fig.savefig('%s/%s_%s_twospecies_substitutions.pdf' % (parse_midas_data.analysis_directory,species_name_1, species_name_2),bbox_inches='tight')
fig.savefig('%s/%s_%s_twospecies_substitutions.png' % (parse_midas_data.analysis_directory,species_name_1, species_name_2),bbox_inches='tight',dpi=300)

    
