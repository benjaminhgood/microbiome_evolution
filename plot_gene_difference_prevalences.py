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

mpl.rcParams['font.size'] = 8
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
    
# Load core gene set
sys.stderr.write("Loading core genes...\n")
core_genes = parse_midas_data.load_core_genes(species_name)
sys.stderr.write("Done! Core genome consists of %d genes\n" % len(core_genes))
    
# Load genomic coverage distributions
sample_coverage_histograms, samples = parse_midas_data.parse_coverage_distribution(species_name)
median_coverages = numpy.array([stats_utils.calculate_nonzero_median_from_histogram(sample_coverage_histogram) for sample_coverage_histogram in sample_coverage_histograms])
sample_coverage_map = {samples[i]: median_coverages[i] for i in xrange(0,len(samples))}

# Load pi information for species_name
sys.stderr.write("Loading within-sample diversity for %s...\n" % species_name)
samples, total_pis, total_pi_opportunities = parse_midas_data.parse_within_sample_pi(species_name, allowed_genes=core_genes, debug=debug)
sys.stderr.write("Done!\n")
pis = total_pis/total_pi_opportunities

median_coverages = numpy.array([sample_coverage_map[samples[i]] for i in xrange(0,len(samples))])

# Only plot samples above a certain depth threshold that are "haploids"
snp_samples = samples[(median_coverages>=min_coverage)*(pis<=1e-03)]

# Load gene coverage information for species_name
sys.stderr.write("Loading pangenome data for %s...\n" % species_name)
gene_samples, gene_names, gene_presence_matrix, gene_depth_matrix, marker_coverages, gene_reads_matrix = parse_midas_data.parse_pangenome_data(species_name,allowed_samples=snp_samples)
sys.stderr.write("Done!\n")
        
# Calculate matrix of number of genes that differ
sys.stderr.write("Calculating matrix of gene differences...\n")
gene_difference_matrix, gene_opportunity_matrix = gene_diversity_utils.calculate_coverage_based_gene_hamming_matrix(gene_depth_matrix, marker_coverages, min_log2_fold_change=4)

# Now need to make the gene samples and snp samples match up
desired_samples = gene_samples[marker_coverages>min_coverage]   

prevalence_idxs = (parse_midas_data.calculate_unique_samples(subject_sample_map, gene_samples))*(marker_coverages>=min_coverage)
    
prevalences = gene_diversity_utils.calculate_fractional_gene_prevalences(gene_depth_matrix[:,prevalence_idxs], marker_coverages[prevalence_idxs])

pangenome_prevalences = numpy.array(prevalences,copy=True)
pangenome_prevalences.sort()
     
# Calculate which pairs of idxs belong to the same sample, which to the same subject
# and which to different subjects
desired_same_sample_idxs, desired_same_subject_idxs, desired_diff_subject_idxs = parse_midas_data.calculate_subject_pairs(subject_sample_map, desired_samples)

snp_sample_idx_map = parse_midas_data.calculate_sample_idx_map(desired_samples, snp_samples)
gene_sample_idx_map = parse_midas_data.calculate_sample_idx_map(desired_samples, gene_samples)
    
same_sample_snp_idxs  = parse_midas_data.apply_sample_index_map_to_indices(snp_sample_idx_map,  desired_same_sample_idxs)  

same_sample_gene_idxs = parse_midas_data.apply_sample_index_map_to_indices(gene_sample_idx_map, desired_same_sample_idxs)  

same_subject_snp_idxs = parse_midas_data.apply_sample_index_map_to_indices(snp_sample_idx_map, desired_same_subject_idxs)  

same_subject_gene_idxs = parse_midas_data.apply_sample_index_map_to_indices(gene_sample_idx_map, desired_same_subject_idxs)  

diff_subject_snp_idxs = parse_midas_data.apply_sample_index_map_to_indices(snp_sample_idx_map, desired_diff_subject_idxs)  

diff_subject_gene_idxs = parse_midas_data.apply_sample_index_map_to_indices(gene_sample_idx_map, desired_diff_subject_idxs)  

within_host_gene_idx_map = {}
for sample_pair_idx in xrange(0,len(same_subject_snp_idxs[0])):
    
    i = same_subject_gene_idxs[0][sample_pair_idx]
    j = same_subject_gene_idxs[1][sample_pair_idx]
    gene_differences = gene_diversity_utils.calculate_gene_differences_between(i, j, gene_depth_matrix, marker_coverages)

    for gene_idx, depth_tuple_1, depth_tuple_2 in gene_differences:
        if gene_idx not in within_host_gene_idx_map:
            within_host_gene_idx_map[gene_idx]=0
            
        within_host_gene_idx_map[gene_idx]+=1


between_host_gene_idx_map = {}
for sample_pair_idx in xrange(0,len(diff_subject_snp_idxs[0])):
        
    i = diff_subject_gene_idxs[0][sample_pair_idx]
    j = diff_subject_gene_idxs[1][sample_pair_idx]
    gene_differences = gene_diversity_utils.calculate_gene_differences_between(i, j, gene_depth_matrix, marker_coverages)

    for gene_idx, depth_tuple_1, depth_tuple_2 in gene_differences:
        if gene_idx not in between_host_gene_idx_map:
            between_host_gene_idx_map[gene_idx]=0
            
        between_host_gene_idx_map[gene_idx]+=1

within_host_gene_sfs = []
within_host_gene_prevalences = []
for gene_idx in within_host_gene_idx_map.keys():
    within_host_gene_sfs.append(within_host_gene_idx_map[gene_idx])
    for i in xrange(0, within_host_gene_idx_map[gene_idx]):
        within_host_gene_prevalences.append(prevalences[gene_idx])
within_host_gene_sfs.sort()
within_host_gene_sfs = numpy.array(within_host_gene_sfs)
within_host_gene_prevalences.sort()
within_host_gene_prevalences = numpy.array(within_host_gene_prevalences)

print within_host_gene_sfs.mean(), within_host_gene_sfs.std(), within_host_gene_sfs.max()

between_host_gene_sfs = []
between_host_gene_prevalences = []
for gene_idx in between_host_gene_idx_map.keys():
    between_host_gene_sfs.append(between_host_gene_idx_map[gene_idx])
    for i in xrange(0, between_host_gene_idx_map[gene_idx]):
        between_host_gene_prevalences.append(prevalences[gene_idx])
between_host_gene_sfs.sort()
between_host_gene_sfs = numpy.array(between_host_gene_sfs)
between_host_gene_prevalences.sort()
between_host_gene_prevalences = numpy.array(between_host_gene_prevalences)

    
# Done calculating... now plot figure!

# Set up figure
fig = plt.figure(figsize=(3.42,2))

# Set up grids to hold figure panels
outer_grid = gridspec.GridSpec(1, 1)

###################
#
# Prevalence
#
###################

prevalence_axis = plt.Subplot(fig, outer_grid[0])
fig.add_subplot(prevalence_axis)
fig.suptitle(species_name)

prevalence_axis.set_ylabel('Fraction genes $\geq p$')
prevalence_axis.set_xlabel('Prevalence of gene, $p$')
prevalence_axis.set_xlim([0,1])
prevalence_axis.set_ylim([0,1])

xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(pangenome_prevalences)
prevalence_axis.step(xs,ns*1.0/ns[0],'k-',label='Total pan genome')

xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(between_host_gene_prevalences)
prevalence_axis.step(xs,ns*1.0/ns[0],'r-',label='Between host differences')

xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(within_host_gene_prevalences)
prevalence_axis.step(xs,ns*1.0/ns[0],'g-',label='Within host differences')


prevalence_axis.legend(loc='upper right',frameon=False,fontsize=6)

fig.savefig('%s/%s_gene_differences_prevalences.pdf' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight')
fig.savefig('%s/%s_gene_differences_prevalences.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight',dpi=300)
