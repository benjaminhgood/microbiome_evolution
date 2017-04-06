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
min_coverage = 40
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

# Only plot samples above a certain depth threshold
# Don't worry about haploids and diploids this time
snp_samples = samples[(median_coverages>=min_coverage)]

# Analyze SNPs, looping over chunk sizes. 
# Clunky, but necessary to limit memory usage on cluster

# Load SNP information for species_name
sys.stderr.write("Loading %s...\n" % species_name)

pi_matrix_syn = numpy.array([])
avg_pi_matrix_syn = numpy.array([])

snp_difference_matrix = numpy.array([])
snp_opportunity_matrix = numpy.array([])

polymorphism_difference_matrix = numpy.array([])
polymorphism_opportunity_matrix = numpy.array([])

final_line_number = 0

while final_line_number >= 0:
    
    sys.stderr.write("Loading chunk starting @ %d...\n" % final_line_number)
    dummy_samples, allele_counts_map, passed_sites_map, final_line_number = parse_midas_data.parse_snps(species_name, debug=debug, allowed_samples=snp_samples,chunk_size=chunk_size,initial_line_number=final_line_number)
    sys.stderr.write("Done! Loaded %d genes\n" % len(allele_counts_map.keys()))
    
    # Calculate fixation matrix
    sys.stderr.write("Calculating matrix of snp differences...\n")
    chunk_snp_difference_matrix, chunk_snp_opportunity_matrix =     diversity_utils.calculate_fixation_matrix(allele_counts_map, passed_sites_map, min_change=min_change)    
    sys.stderr.write("Done!\n")
    
    sys.stderr.write("Calculating matrix of polymorphism differences...\n")
    chunk_polymorphism_difference_matrix, chunk_polymorphism_opportunity_matrix =     diversity_utils.calculate_new_snp_matrix(allele_counts_map, passed_sites_map, min_freq=0.05, max_freq=0.10)    
    sys.stderr.write("Done!\n")
    
    if snp_difference_matrix.shape[0]==0:
        snp_difference_matrix = numpy.zeros_like(chunk_snp_difference_matrix)*1.0
        snp_opportunity_matrix = numpy.zeros_like(snp_difference_matrix)*1.0
        polymorphism_difference_matrix = numpy.zeros_like(snp_difference_matrix)*1.0
        polymorphism_opportunity_matrix = numpy.zeros_like(snp_difference_matrix)*1.0

    
    snp_difference_matrix += chunk_snp_difference_matrix
    snp_opportunity_matrix += chunk_snp_opportunity_matrix
    polymorphism_difference_matrix += chunk_polymorphism_difference_matrix
    polymorphism_opportunity_matrix += chunk_polymorphism_opportunity_matrix

sys.stderr.write("Done!\n")   

# Load gene coverage information for species_name
sys.stderr.write("Loading pangenome data for %s...\n" % species_name)
gene_samples, gene_names, gene_presence_matrix, gene_depth_matrix, marker_coverages, gene_reads_matrix = parse_midas_data.parse_pangenome_data(species_name,allowed_samples=snp_samples)
sys.stderr.write("Done!\n")
        
# Calculate matrix of number of genes that differ
sys.stderr.write("Calculating matrix of gene differences...\n")
gene_difference_matrix, gene_opportunity_matrix = gene_diversity_utils.calculate_coverage_based_gene_hamming_matrix(gene_depth_matrix, marker_coverages, min_log2_fold_change=4)

# Now need to make the gene samples and snp samples match up
desired_samples = gene_samples[marker_coverages>min_coverage]   
     
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

same_subject_pi_plowers = []
same_subject_pi_pmids = []
same_subject_pi_puppers = []

same_subject_snp_plowers = []
same_subject_snp_pmids = []
same_subject_snp_puppers = []

same_subject_polymorphism_plowers = []
same_subject_polymorphism_pmids = []
same_subject_polymorphism_puppers = []

same_subject_gene_plowers = []
same_subject_gene_pmids = []
same_subject_gene_puppers = []

for sample_pair_idx in xrange(0,len(same_subject_snp_idxs[0])):
    
    i = same_subject_snp_idxs[0][sample_pair_idx]
    j = same_subject_snp_idxs[1][sample_pair_idx]

    i2 = desired_same_subject_idxs[0][sample_pair_idx]
    j2 = desired_same_subject_idxs[1][sample_pair_idx] 
 
    if pis[i2]>pis[j2]:
        pi_idx = i2
    else:
        pi_idx = j2

    # Calculate max pi between two samples
    plower,pupper = stats_utils.calculate_poisson_rate_interval(total_pis[pi_idx], total_pi_opportunities[pi_idx],alpha)
    pmid = total_pis[pi_idx]*1.0/total_pi_opportunities[pi_idx]
    
    same_subject_pi_plowers.append(plower)
    same_subject_pi_pmids.append(pmid)
    same_subject_pi_puppers.append(pupper)
    
    # Calculate SNP changes
    plower,pupper = stats_utils.calculate_poisson_rate_interval(snp_difference_matrix[i,j], snp_opportunity_matrix[i,j],alpha)
    
    pmid = snp_difference_matrix[i,j]*1.0/snp_opportunity_matrix[i,j]
    
    same_subject_snp_plowers.append(plower)
    same_subject_snp_pmids.append(pmid)
    same_subject_snp_puppers.append(pupper)
     
    # Calculate polymorphism changes
    plower,pupper = stats_utils.calculate_poisson_rate_interval(polymorphism_difference_matrix[i,j], polymorphism_opportunity_matrix[i,j],alpha)
    
    pmid = polymorphism_difference_matrix[i,j]*1.0/polymorphism_opportunity_matrix[i,j]
    
    same_subject_polymorphism_plowers.append(plower)
    same_subject_polymorphism_pmids.append(pmid)
    same_subject_polymorphism_puppers.append(pupper)
    
    i = same_subject_gene_idxs[0][sample_pair_idx]
    j = same_subject_gene_idxs[1][sample_pair_idx]
    
    
    plower,pupper = stats_utils.calculate_poisson_rate_interval(gene_difference_matrix[i,j], gene_opportunity_matrix[i,j])
    pmid = gene_difference_matrix[i,j]*1.0/gene_opportunity_matrix[i,j]
    
    same_subject_gene_plowers.append(plower)
    same_subject_gene_pmids.append(pmid)
    same_subject_gene_puppers.append(pupper)

# clip lower bounds 
same_subject_gene_plowers = numpy.clip(same_subject_gene_plowers,1e-09,1e09)
same_subject_snp_plowers = numpy.clip(same_subject_snp_plowers,1e-09,1e09)
same_subject_pi_plowers = numpy.clip(same_subject_pi_plowers,1e-09,1e09)
same_subject_polymorphism_plowers = numpy.clip(same_subject_polymorphism_plowers,1e-09,1e09)


diff_subject_pi_plowers = []
diff_subject_pi_pmids = []
diff_subject_pi_puppers = []

diff_subject_snp_plowers = []
diff_subject_snp_pmids = []
diff_subject_snp_puppers = []

diff_subject_polymorphism_plowers = []
diff_subject_polymorphism_pmids = []
diff_subject_polymorphism_puppers = []

diff_subject_gene_plowers = []
diff_subject_gene_pmids = []
diff_subject_gene_puppers = []

for sample_pair_idx in xrange(0,len(diff_subject_snp_idxs[0])):
    
    i = diff_subject_snp_idxs[0][sample_pair_idx]
    j = diff_subject_snp_idxs[1][sample_pair_idx]

    i2 = desired_diff_subject_idxs[0][sample_pair_idx]
    j2 = desired_diff_subject_idxs[1][sample_pair_idx] 
 
    if pis[i2]>pis[j2]:
        pi_idx = i2
    else:
        pi_idx = j2

    # Calculate max pi between two samples
    plower,pupper = stats_utils.calculate_poisson_rate_interval(total_pis[pi_idx], total_pi_opportunities[pi_idx],alpha)
    pmid = total_pis[pi_idx]*1.0/total_pi_opportunities[pi_idx]
    
    diff_subject_pi_plowers.append(plower)
    diff_subject_pi_pmids.append(pmid)
    diff_subject_pi_puppers.append(pupper)
    
    # Calculate SNP changes
    plower,pupper = stats_utils.calculate_poisson_rate_interval(snp_difference_matrix[i,j], snp_opportunity_matrix[i,j],alpha)
    
    pmid = snp_difference_matrix[i,j]*1.0/snp_opportunity_matrix[i,j]
    
    diff_subject_snp_plowers.append(plower)
    diff_subject_snp_pmids.append(pmid)
    diff_subject_snp_puppers.append(pupper)
     
    # Calculate polymorphism changes
    plower,pupper = stats_utils.calculate_poisson_rate_interval(polymorphism_difference_matrix[i,j], polymorphism_opportunity_matrix[i,j],alpha)
    
    pmid = polymorphism_difference_matrix[i,j]*1.0/polymorphism_opportunity_matrix[i,j]
    
    diff_subject_polymorphism_plowers.append(plower)
    diff_subject_polymorphism_pmids.append(pmid)
    diff_subject_polymorphism_puppers.append(pupper)
    
    i = diff_subject_gene_idxs[0][sample_pair_idx]
    j = diff_subject_gene_idxs[1][sample_pair_idx]
    
    plower,pupper = stats_utils.calculate_poisson_rate_interval(gene_difference_matrix[i,j], gene_opportunity_matrix[i,j])
    pmid = gene_difference_matrix[i,j]*1.0/gene_opportunity_matrix[i,j]
    
    diff_subject_gene_plowers.append(plower)
    diff_subject_gene_pmids.append(pmid)
    diff_subject_gene_puppers.append(pupper)

# clip lower bounds 
diff_subject_gene_plowers = numpy.clip(diff_subject_gene_plowers,1e-09,1e09)
diff_subject_snp_plowers = numpy.clip(diff_subject_snp_plowers,1e-09,1e09)
diff_subject_pi_plowers = numpy.clip(diff_subject_pi_plowers,1e-09,1e09)
diff_subject_polymorphism_plowers = numpy.clip(diff_subject_polymorphism_plowers,1e-09,1e09)

# Set up figure
fig = plt.figure(figsize=(3.42, 4))

# Set up grids to hold figure panels
outer_grid = gridspec.GridSpec(2, 1, height_ratios=[1,1], hspace=0.05)

###################
#
# SNP Panel
#
###################

snp_axis = plt.Subplot(fig, outer_grid[0])
fig.add_subplot(snp_axis)
fig.suptitle(species_name)

snp_axis.set_ylabel('Substitution rate')
snp_axis.set_ylim([1e-07,1e-01])

#snp_axis.semilogx([1e-09,1e-09],[1,1],'g.',label='Within host')
#snp_axis.semilogx([1e-09,1e-09],[1,1],'r.',label='Between host')

###################
#
# New polymorphism panel
#
###################

polymorphism_axis = plt.Subplot(fig, outer_grid[1])
fig.add_subplot(polymorphism_axis)

polymorphism_axis.set_ylabel('New polymorphism rate')
polymorphism_axis.set_ylim([1e-07,1e-01])
polymorphism_axis.set_xlabel('Max within-host diversity')
polymorphism_axis.set_ylim([1e-07,1e-01])

snp_axis.loglog( numpy.clip(diff_subject_pi_pmids, 1e-07, 1e09), numpy.clip(diff_subject_snp_pmids, 1e-07, 1e09), 'r.',markersize=2)

snp_axis.plot( numpy.clip(same_subject_pi_pmids, 1e-07, 1e09), numpy.clip(same_subject_snp_pmids, 1e-07, 1e09), 'g.',markersize=3)

polymorphism_axis.loglog( numpy.clip(diff_subject_pi_pmids, 1e-07, 1e09), numpy.clip(diff_subject_polymorphism_pmids, 1e-07, 1e09), 'r.',markersize=2)

polymorphism_axis.plot( numpy.clip(same_subject_pi_pmids, 1e-07, 1e09), numpy.clip(same_subject_polymorphism_pmids, 1e-07, 1e09), 'g.',markersize=3)

#snp_axis.legend(loc='lower left',frameon=False)

fig.savefig('%s/%s_fixations_snps_vs_pi.pdf' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight')
fig.savefig('%s/%s_fixations_snps_vs_pi.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight',dpi=300)

    
