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

from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster.hierarchy import cophenet
from scipy.cluster.hierarchy import fcluster

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
sys.stderr.write("Loading %d samples for %s...\n" % (len(snp_samples), species_name))

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
    
     
substitution_rate = snp_difference_matrix*1.0/snp_opportunity_matrix 


cluster_idxss = diversity_utils.cluster_samples(substitution_rate, min_d=0, max_d=1e09)
print len(cluster_idxss), [cluster_idxs[0].sum() for cluster_idxs in cluster_idxss]

cluster_idxss = diversity_utils.cluster_samples(substitution_rate, min_d=5e-04, max_d=1e-02)
#cluster_idxss = diversity_utils.cluster_samples(substitution_rate, min_d=0, max_d=1e-04)
#cluster_idxss = diversity_utils.cluster_samples(substitution_rate, min_d=0, max_d=1e-04)

print len(cluster_idxss), [cluster_idxs[0].sum() for cluster_idxs in cluster_idxss]


# calculate compressed distance matrix suitable for agglomerative clustering
Y = []
for i in xrange(0,substitution_rate.shape[0]):
    for j in xrange(i+1,substitution_rate.shape[1]):
        Y.append(substitution_rate[i,j]) 
Y = numpy.array(Y) 
    
Z = linkage(Y, method='average')        
    
c, coph_dists = cophenet(Z, Y)
sys.stderr.write("cophenetic correlation: %g\n" % c)

# Make a dendrogram
pylab.figure(1, figsize=(15, 5))
pylab.title('Full hierarchical clustering dendrogram for %s' % species_name)
pylab.xlabel('sample index')
pylab.ylabel('distance')
dendrogram(
    Z,
    leaf_rotation=90.,  # rotates the x axis labels
    leaf_font_size=8.,  # font size for the x axis labels
)
pylab.savefig('%s/%s_full_dendrogram.pdf' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight')
pylab.savefig('%s/%s_full_dendrogram.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight',dpi=300)

# Make a zoomed dendrogram
pylab.figure(2, figsize=(15, 5))
pylab.title('Zoomed hierarchical clustering dendrogram for %s' % species_name)
pylab.xlabel('sample index')
pylab.ylabel('distance')
dendrogram(
    Z,
    leaf_rotation=90.,  # rotates the x axis labels
    leaf_font_size=8.,  # font size for the x axis labels
)
pylab.ylim([0,5e-04])
pylab.savefig('%s/%s_zoomed_dendrogram.pdf' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight')
pylab.savefig('%s/%s_zoomed_dendrogram.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight',dpi=300)


low_divergence_idxs = numpy.logical_or(cluster_idxss[0][0],cluster_idxss[0][1])

low_divergence_cluster_idxss = [(cluster_idxs[0][low_divergence_idxs], cluster_idxs[1][low_divergence_idxs]) for cluster_idxs in cluster_idxss]

low_divergence_samples = snp_samples[low_divergence_idxs]

if len(low_divergence_samples)>2:

    sys.stderr.write("Continuing with %d samples...\n" % len(low_divergence_samples))

    # Load SNP information for species_name
    sys.stderr.write("Re-loading %s...\n" % species_name)

    snp_difference_matrix = numpy.array([])
    snp_opportunity_matrix = numpy.array([])
    
    total_polymorphic_sites = 0
    total_inconsistent_sites = 0

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
    
        sys.stderr.write("Calculating phylogenetic consistency...\n")
        chunk_inconsistent_sites, chunk_polymorphic_sites = diversity_utils.calculate_phylogenetic_consistency(allele_counts_map, passed_sites_map, low_divergence_cluster_idxss)    
        sys.stderr.write("Done!\n")
    
        total_inconsistent_sites += chunk_inconsistent_sites
        total_polymorphic_sites += chunk_polymorphic_sites
    
    substitution_rate = snp_difference_matrix*1.0/snp_opportunity_matrix 

    sys.stderr.write('%d SNPs, %d inconsistent (%g)\n' % (total_polymorphic_sites, total_inconsistent_sites, total_inconsistent_sites*1.0/total_polymorphic_sites))
    

    # calculate compressed distance matrix suitable for agglomerative clustering
    Y = []
    for i in xrange(0,substitution_rate.shape[0]):
        for j in xrange(i+1,substitution_rate.shape[1]):
            Y.append(substitution_rate[i,j]) 
    Y = numpy.array(Y) 
    
    Z = linkage(Y, method='average')        
    
    c, coph_dists = cophenet(Z, Y)

    sys.stderr.write("cophenetic correlation: %g\n" % c)

    pylab.figure(3, figsize=(15, 5))
    pylab.title('Low-divergence Dendrogram for %s' % species_name)
    pylab.xlabel('sample index')
    pylab.ylabel('distance')
    dendrogram(
        Z,
        leaf_rotation=90.,  # rotates the x axis labels
        leaf_font_size=8.,  # font size for the x axis labels
    )

    pylab.savefig('%s/%s_low_divergence_samples.pdf' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight')
    pylab.savefig('%s/%s_low_divergence_samples.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight',dpi=300)

else:
    
    sys.stderr.write("No low divergence samples!\n")
    
