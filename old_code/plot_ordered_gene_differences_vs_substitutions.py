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

mpl.rcParams['font.size'] = 7
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
low_pi_threshold = 1e-03
low_divergence_threshold = 1e-03

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

# Analyze SNPs, looping over chunk sizes. 
# Clunky, but necessary to limit memory usage on cluster

# Load SNP information for species_name
sys.stderr.write("Loading SNPs for %s...\n" % species_name)
sys.stderr.write("(not just core genes...)\n")
pi_matrix_syn = numpy.array([])
avg_pi_matrix_syn = numpy.array([])
snp_difference_matrix = numpy.array([])
snp_opportunity_matrix = numpy.array([])
final_line_number = 0

while final_line_number >= 0:
    
    sys.stderr.write("Loading chunk starting @ %d...\n" % final_line_number)
    dummy_samples, allele_counts_map, passed_sites_map, final_line_number = parse_midas_data.parse_snps(species_name, debug=debug, allowed_samples=snp_samples,chunk_size=chunk_size,initial_line_number=final_line_number)
    sys.stderr.write("Done! Loaded %d genes\n" % len(allele_counts_map.keys()))
    
    # Calculate fixation matrix
    sys.stderr.write("Calculating matrix of snp differences...\n")
    chunk_snp_difference_matrix, chunk_snp_opportunity_matrix =     diversity_utils.calculate_fixation_matrix(allele_counts_map, passed_sites_map, min_change=min_change)   # 
    sys.stderr.write("Done!\n")
    
    if snp_difference_matrix.shape[0]==0:
        snp_difference_matrix = numpy.zeros_like(chunk_snp_difference_matrix)*1.0
        snp_opportunity_matrix = numpy.zeros_like(snp_difference_matrix)*1.0
    
    snp_difference_matrix += chunk_snp_difference_matrix
    snp_opportunity_matrix += chunk_snp_opportunity_matrix


snp_substitution_rate = snp_difference_matrix*1.0/(snp_opportunity_matrix+(snp_opportunity_matrix==0))
sys.stderr.write("Done!\n")   

# Load gene coverage information for species_name
sys.stderr.write("Loading pangenome data for %s...\n" % species_name)
gene_samples, gene_names, gene_presence_matrix, gene_depth_matrix, marker_coverages, gene_reads_matrix = parse_midas_data.parse_pangenome_data(species_name,allowed_samples=snp_samples)
sys.stderr.write("Done!\n")

prevalence_idxs = (parse_midas_data.calculate_unique_samples(subject_sample_map, gene_samples))*(marker_coverages>=min_coverage)
    
prevalences = gene_diversity_utils.calculate_fractional_gene_prevalences(gene_depth_matrix[:,prevalence_idxs], marker_coverages[prevalence_idxs])

pangenome_prevalences = numpy.array(prevalences,copy=True)
pangenome_prevalences.sort()
        
# Calculate matrix of number of genes that differ
sys.stderr.write("Calculating matrix of gene differences...\n")
gene_difference_matrix, gene_opportunity_matrix = gene_diversity_utils.calculate_coverage_based_gene_hamming_matrix(gene_depth_matrix, marker_coverages, min_log2_fold_change=4)

# Now need to make the gene samples and snp samples match up
desired_samples = gene_samples[marker_coverages>min_coverage]      

num_haploids = len(desired_samples)
     
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

typical_same_subject_snp_opportunities = numpy.median(snp_opportunity_matrix[same_subject_snp_idxs])
typical_diff_subject_snp_opportunities = numpy.median(snp_opportunity_matrix[diff_subject_snp_idxs])


typical_same_subject_gene_opportunities = numpy.median(gene_opportunity_matrix[same_subject_gene_idxs])
typical_diff_subject_gene_opportunities = numpy.median(gene_opportunity_matrix[diff_subject_gene_idxs])

#print typical_same_subject_snp_opportunities, typical_diff_subject_snp_opportunities
#print typical_same_subject_gene_opportunities, typical_diff_subject_gene_opportunities

Lsnps = typical_diff_subject_snp_opportunities
Lgenes = typical_diff_subject_gene_opportunities

same_subject_snp_plowers = []
same_subject_snp_puppers = []
same_subject_gene_plowers = []
same_subject_gene_puppers = []
within_host_gene_idx_map = {}
for sample_pair_idx in xrange(0,len(same_subject_snp_idxs[0])):
#    
    snp_i = same_subject_snp_idxs[0][sample_pair_idx]
    snp_j = same_subject_snp_idxs[1][sample_pair_idx]
#    
    plower,pupper = stats_utils.calculate_poisson_rate_interval(snp_difference_matrix[snp_i, snp_j], snp_opportunity_matrix[snp_i, snp_j],alpha)
#    
    same_subject_snp_plowers.append(plower)
    same_subject_snp_puppers.append(pupper)
#    
    #snp_differences = diversity_utils.calculate_snp_differences_between(i,j,allele_counts_map, passed_sites_map, min_change=min_change)
#
    i = same_subject_gene_idxs[0][sample_pair_idx]
    j = same_subject_gene_idxs[1][sample_pair_idx]
    gene_differences = gene_diversity_utils.calculate_gene_differences_between(i, j, gene_depth_matrix, marker_coverages)

    if snp_substitution_rate[snp_i,snp_j] < low_divergence_threshold:
        for gene_idx, depth_tuple_1, depth_tuple_2 in gene_differences:
            if gene_idx not in within_host_gene_idx_map:
                within_host_gene_idx_map[gene_idx]=0
#            
            within_host_gene_idx_map[gene_idx]+=1
#
    plower,pupper = stats_utils.calculate_poisson_rate_interval(gene_difference_matrix[i,j], gene_opportunity_matrix[i,j])
#    
    same_subject_gene_plowers.append(plower)
    same_subject_gene_puppers.append(pupper)
#

# clip lower bounds 
same_subject_gene_plowers = numpy.clip(same_subject_gene_plowers,1e-06,1e09)
same_subject_snp_plowers = numpy.clip(same_subject_snp_plowers,1e-08,1e09)
# Sort all four lists by ascending lower bound on SNP changes, then gene changes
same_subject_snp_plowers, same_subject_gene_plowers, same_subject_snp_puppers,  same_subject_gene_puppers = (numpy.array(x) for x in zip(*sorted(zip(same_subject_snp_plowers, same_subject_gene_plowers, same_subject_snp_puppers, same_subject_gene_puppers))))


diff_subject_snp_plowers = []
diff_subject_snp_puppers = []
diff_subject_gene_plowers = []
diff_subject_gene_puppers = []
between_host_gene_idx_map = {}
low_divergence_between_host_gene_idx_map = {}
for sample_pair_idx in xrange(0,len(diff_subject_snp_idxs[0])):
    
    snp_i = diff_subject_snp_idxs[0][sample_pair_idx]
    snp_j = diff_subject_snp_idxs[1][sample_pair_idx]
    
    plower,pupper = stats_utils.calculate_poisson_rate_interval(snp_difference_matrix[snp_i,snp_j], snp_opportunity_matrix[snp_i, snp_j])
    
    diff_subject_snp_plowers.append(plower)
    diff_subject_snp_puppers.append(pupper)
        
    i = diff_subject_gene_idxs[0][sample_pair_idx]
    j = diff_subject_gene_idxs[1][sample_pair_idx]
    gene_differences = gene_diversity_utils.calculate_gene_differences_between(i, j, gene_depth_matrix, marker_coverages)

    for gene_idx, depth_tuple_1, depth_tuple_2 in gene_differences:
        if gene_idx not in between_host_gene_idx_map:
            between_host_gene_idx_map[gene_idx]=0
            
        between_host_gene_idx_map[gene_idx]+=1
    
    if snp_substitution_rate[snp_i,snp_j] < low_divergence_threshold:
        for gene_idx, depth_tuple_1, depth_tuple_2 in gene_differences:
            if gene_idx not in low_divergence_between_host_gene_idx_map:
                low_divergence_between_host_gene_idx_map[gene_idx]=0
            
            low_divergence_between_host_gene_idx_map[gene_idx]+=1
            
    plower,pupper = stats_utils.calculate_poisson_rate_interval(gene_difference_matrix[i,j], gene_opportunity_matrix[i,j],alpha)
    
    diff_subject_gene_plowers.append(plower)
    diff_subject_gene_puppers.append(pupper)

# clip lower bounds 
diff_subject_gene_plowers = numpy.clip(diff_subject_gene_plowers,1e-06,1e09)
diff_subject_snp_plowers = numpy.clip(diff_subject_snp_plowers,1e-08,1e09)
# Sort all four lists by ascending lower bound on SNP changes, then gene changes
diff_subject_snp_plowers, diff_subject_gene_plowers, diff_subject_snp_puppers,  diff_subject_gene_puppers = (numpy.array(x) for x in zip(*sorted(zip(diff_subject_snp_plowers, diff_subject_gene_plowers, diff_subject_snp_puppers, diff_subject_gene_puppers))))

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

low_divergence_between_host_gene_sfs = []
low_divergence_between_host_gene_prevalences = []
for gene_idx in low_divergence_between_host_gene_idx_map.keys():
    low_divergence_between_host_gene_sfs.append(low_divergence_between_host_gene_idx_map[gene_idx])
    for i in xrange(0, low_divergence_between_host_gene_idx_map[gene_idx]):
        low_divergence_between_host_gene_prevalences.append(prevalences[gene_idx])

low_divergence_between_host_gene_sfs.sort()
low_divergence_between_host_gene_sfs = numpy.array(low_divergence_between_host_gene_sfs)
low_divergence_between_host_gene_prevalences.sort()
low_divergence_between_host_gene_prevalences = numpy.array(low_divergence_between_host_gene_prevalences)


# Done calculating... now plot figure!

# Set up figure
fig = plt.figure(figsize=(7, 3))

# Set up grids to hold figure panels
outer_grid = gridspec.GridSpec(2, 1, height_ratios=[1,1], hspace=0.05)

###################
#
# SNP Panel
#
###################

snp_axis = plt.Subplot(fig, outer_grid[0])
fig.add_subplot(snp_axis)

snp_axis.set_title("%s %s (%s)" % tuple(species_name.split("_")),fontsize=7)

snp_axis.set_ylabel('SNP substitution rate')
snp_axis.set_ylim([1e-07,1e-01])

snp_axis.semilogy([1e-09,1e-09],[1,1],'g-',label='Within host')
snp_axis.semilogy([1e-09,1e-09],[1,1],'r-',label='Between host')

###################
#
# Gene Panel
#
###################

gene_axis = plt.Subplot(fig, outer_grid[1])
fig.add_subplot(gene_axis)

gene_axis.set_ylabel('Gene gain/loss rate')
gene_axis.set_ylim([1e-05,0.9])

gene_axis.set_xlabel('Sample pairs ($n=%d$)' % num_haploids)

scale_factor = 100.0/len(same_subject_snp_plowers)

y = 0
for snp_plower, snp_pupper, gene_plower, gene_pupper in zip(same_subject_snp_plowers, same_subject_snp_puppers, same_subject_gene_plowers, same_subject_gene_puppers):

    y-=2
    
    snp_axis.semilogy([y,y], [snp_plower,snp_pupper],'g-',linewidth=0.25)
    snp_axis.semilogy([y], [snp_plower],'g.',markersize=1.5)
        
    gene_axis.semilogy([y,y], [gene_plower,gene_pupper], 'g-',linewidth=0.25)
    gene_axis.semilogy([y],[gene_plower],  'g.',markersize=1.5)

y-=5
snp_axis.semilogy([y,y],[1e-09,1e09],'-',linewidth=0.25,color='0.7')
gene_axis.semilogy([y,y],[1e-09,1e09],'-',linewidth=0.25,color='0.7')

y-=4

    
    
for snp_plower, snp_pupper, gene_plower, gene_pupper in zip(diff_subject_snp_plowers, diff_subject_snp_puppers, diff_subject_gene_plowers, diff_subject_gene_puppers)[0:100]:

    y-=1
    
    snp_axis.semilogy([y,y],[snp_plower,snp_pupper],'r-',linewidth=0.35)
    gene_axis.semilogy([y,y],[gene_plower,gene_pupper],'r-',linewidth=0.35)

y-=5
snp_axis.semilogy([y,y],[1e-09,1e09],'-',linewidth=0.25,color='0.7')
gene_axis.semilogy([y,y],[1e-09,1e09],'-',linewidth=0.25,color='0.7')
y-=4

if len(diff_subject_snp_plowers)<=500:
    idxs = numpy.arange(0,len(diff_subject_snp_plowers))
else:
    idxs = randint(0,len(diff_subject_snp_plowers),500)
    idxs.sort()

for idx in idxs:
    
    snp_plower = diff_subject_snp_plowers[idx]
    snp_pupper = diff_subject_snp_puppers[idx]
    gene_plower = diff_subject_gene_plowers[idx]
    gene_pupper = diff_subject_gene_puppers[idx]
    
    y-=1
    
    snp_axis.semilogy([y,y],[snp_plower,snp_pupper],'r-',linewidth=0.35)
    gene_axis.semilogy([y,y],[gene_plower,gene_pupper],'r-',linewidth=0.35)
        
y-=3
snp_axis.set_xlim([y-1,0])
gene_axis.set_xlim([y-1,0])    

snp_axis.set_xticks([])
gene_axis.set_xticks([])

#snp_axis.legend(loc='upper right',frameon=False)

#labels = snp_axis.get_yticklabels()
#print labels[0].get_text()
#labels[0].set_text('0')
#snp_axis.set_yticklabels(labels)

fig.savefig('%s/%s_ordered_gene_differences_vs_substitutions.pdf' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight')
fig.savefig('%s/%s_ordered_gene_differences_vs_substitutions.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight',dpi=300)

###################
#
# Prevalence
#
###################

# Set up figure
prevalence_fig = plt.figure(figsize=(3.42,2))

# Set up grids to hold figure panels
outer_grid = gridspec.GridSpec(1, 1)


prevalence_axis = plt.Subplot(prevalence_fig, outer_grid[0])
prevalence_fig.add_subplot(prevalence_axis)
prevalence_fig.suptitle(species_name,fontsize=7)

prevalence_axis.set_ylabel('Fraction genes $\geq p$')
prevalence_axis.set_xlabel('Prevalence of gene, $p$')
prevalence_axis.set_xlim([0,1])
prevalence_axis.set_ylim([0,1])

xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(pangenome_prevalences)
prevalence_axis.step(xs,ns*1.0/ns[0],'k-',label='Total pan genome')

xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(between_host_gene_prevalences)
prevalence_axis.step(xs,ns*1.0/ns[0],'r-',label='Between host differences')

if len(low_divergence_between_host_gene_prevalences) > 0:
    print low_divergence_between_host_gene_prevalences
    print low_divergence_between_host_gene_prevalences.mean()
    print len(low_divergence_between_host_gene_prevalences), len(between_host_gene_prevalences)
    xs, ns =     stats_utils.calculate_unnormalized_survival_from_vector( low_divergence_between_host_gene_prevalences)
    prevalence_axis.step(xs,ns*1.0/ns[0],'r-',label=('d<%g' % low_divergence_threshold), alpha=0.5)

xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(within_host_gene_prevalences)
prevalence_axis.step(xs,ns*1.0/ns[0],'g-',label='Within host differences')



prevalence_axis.legend(loc='upper right',frameon=False,fontsize=6)

prevalence_fig.savefig('%s/%s_gene_differences_prevalences.pdf' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight')
prevalence_fig.savefig('%s/%s_gene_differences_prevalences.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight',dpi=300)
    
    
