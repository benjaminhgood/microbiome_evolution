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
       
# Only plot samples above a certain depth threshold that are involved in timecourse
snp_samples = diversity_utils.calculate_temporal_samples(species_name)

# The subset of samples that are haploid
haploid_samples = set(diversity_utils.calculate_haploid_samples(species_name))

# Only use the subset from North America 
# The only temporal samples are from here, best not contaminate the between-subject
# comparisons with out of sample effects
#snp_samples = snp_samples[parse_HMP_data.calculate_country_samples(sample_country_map, sample_list=snp_samples, allowed_countries=set(["United States"]))]

####################################################
#
# Set up Figure (4 panels, arranged in 2x2 grid)
#
####################################################

pylab.figure(1,figsize=(5,2))
fig = pylab.gcf()
# make three panels panels
outer_grid  = gridspec.GridSpec(1,2, width_ratios=[2,1], wspace=0.25)

differences_grid = gridspec.GridSpecFromSubplotSpec(2, 2, height_ratios=[1,1],
                subplot_spec=outer_grid[0], hspace=0.05, width_ratios=[1,1], wspace=0.025)
                
gene_grid = gridspec.GridSpecFromSubplotSpec(2, 1, height_ratios=[1,1],
                subplot_spec=outer_grid[1], hspace=0.45)

###################
#
# SNP change panel
#
###################

snp_axis = plt.Subplot(fig, differences_grid[0,0])
fig.add_subplot(snp_axis)

#snp_axis.set_title("%s %s (%s)" % tuple(species_name.split("_")),fontsize=7)

snp_axis.set_ylabel('SNP changes')
snp_axis.set_ylim([2e-01,1e05])

snp_axis.semilogy([1e-09,1e-09],[1e08,1e08],'g-',label='Within host')
snp_axis.semilogy([1e-09,1e-09],[1e08,1e08],'r-',label='Between host')

snp_axis.spines['top'].set_visible(False)
snp_axis.spines['right'].set_visible(False)
snp_axis.get_xaxis().tick_bottom()
snp_axis.get_yaxis().tick_left()


within_snp_axis = plt.Subplot(fig, differences_grid[0,1])
fig.add_subplot(within_snp_axis)

#snp_axis.set_title("%s %s (%s)" % tuple(species_name.split("_")),fontsize=7)

within_snp_axis.set_ylim([2e-01,1e05])

within_snp_axis.spines['top'].set_visible(False)
within_snp_axis.spines['right'].set_visible(False)
within_snp_axis.spines['left'].set_visible(False)

within_snp_axis.get_xaxis().tick_bottom()
within_snp_axis.get_yaxis().tick_left()



###################
#
# Gene change panel
#
###################

gene_axis = plt.Subplot(fig, differences_grid[1,0])
fig.add_subplot(gene_axis)

gene_axis.set_ylabel('Gene changes')
gene_axis.set_ylim([2e-01,1e04])

gene_axis.set_xlabel('Between-host')

gene_axis.spines['top'].set_visible(False)
gene_axis.spines['right'].set_visible(False)
gene_axis.get_xaxis().tick_bottom()
gene_axis.get_yaxis().tick_left()


within_gene_axis = plt.Subplot(fig, differences_grid[1,1])
fig.add_subplot(within_gene_axis)

within_gene_axis.set_ylim([2e-01,1e04])

within_gene_axis.spines['top'].set_visible(False)
within_gene_axis.spines['right'].set_visible(False)
within_gene_axis.spines['left'].set_visible(False)
within_gene_axis.get_xaxis().tick_bottom()
within_gene_axis.get_yaxis().tick_left()

within_gene_axis.set_xlabel('Within-host')

##############################################################################
#
# Gene change prevalence panel
#
##############################################################################

prevalence_axis = plt.Subplot(fig, gene_grid[0])
fig.add_subplot(prevalence_axis)

prevalence_axis.set_ylabel('Fraction gene changes ')
prevalence_axis.set_xlabel('Gene prevalence, $p$')
prevalence_axis.set_xlim([0,1])
#prevalence_axis.set_ylim([0,1.1])

prevalence_axis.spines['top'].set_visible(False)
prevalence_axis.spines['right'].set_visible(False)
prevalence_axis.get_xaxis().tick_bottom()
prevalence_axis.get_yaxis().tick_left()


##############################################################################
#
# Gene change parallelism panel
#
##############################################################################

parallelism_axis = plt.Subplot(fig, gene_grid[1])
fig.add_subplot(parallelism_axis)

parallelism_axis.set_ylabel('Fraction gene changes')
parallelism_axis.set_xlabel('Gene multiplicity, $m$')
parallelism_axis.set_xlim([0.5,3.5])
parallelism_axis.set_ylim([0,1.05])

parallelism_axis.set_xticks([1,2,3])

parallelism_axis.spines['top'].set_visible(False)
parallelism_axis.spines['right'].set_visible(False)
parallelism_axis.get_xaxis().tick_bottom()
parallelism_axis.get_yaxis().tick_left()


# Analyze SNPs, looping over chunk sizes. 
# Clunky, but necessary to limit memory usage on cluster

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
        
# Calculate matrix of number of genes that differ
sys.stderr.write("Calculating matrix of gene differences...\n")
gene_gain_matrix, gene_loss_matrix, gene_opportunity_matrix = gene_diversity_utils.calculate_coverage_based_gene_hamming_matrix_gain_loss(gene_depth_matrix, marker_coverages, min_log2_fold_change=4, include_high_copynum=include_high_copynum)

gene_difference_matrix = gene_gain_matrix + gene_loss_matrix

# Now need to make the gene samples and snp samples match up
desired_samples = gene_samples


num_haploids = len(desired_samples)
     
# Calculate which pairs of idxs belong to the same sample, which to the same subject
# and which to different subjects
#desired_same_sample_idxs, desired_same_subject_idxs, desired_diff_subject_idxs = parse_midas_data.calculate_subject_pairs(subject_sample_map, desired_samples)

# Calculate which pairs of idxs belong to the same sample, which to the same subject
# and which to different subjects
desired_same_sample_idxs, desired_same_subject_idxs, desired_diff_subject_idxs = parse_midas_data.calculate_ordered_subject_pairs(sample_order_map, desired_samples)


snp_sample_idx_map = parse_midas_data.calculate_sample_idx_map(desired_samples, snp_samples)
gene_sample_idx_map = parse_midas_data.calculate_sample_idx_map(desired_samples, gene_samples)
  

same_sample_snp_idxs  = parse_midas_data.apply_sample_index_map_to_indices(snp_sample_idx_map,  desired_same_sample_idxs)  
same_sample_gene_idxs = parse_midas_data.apply_sample_index_map_to_indices(gene_sample_idx_map, desired_same_sample_idxs)  


same_subject_snp_idxs = parse_midas_data.apply_sample_index_map_to_indices(snp_sample_idx_map, desired_same_subject_idxs)  
same_subject_gene_idxs = parse_midas_data.apply_sample_index_map_to_indices(gene_sample_idx_map, desired_same_subject_idxs)  

diff_subject_snp_idxs = parse_midas_data.apply_sample_index_map_to_indices(snp_sample_idx_map, desired_diff_subject_idxs)  
diff_subject_gene_idxs = parse_midas_data.apply_sample_index_map_to_indices(gene_sample_idx_map, desired_diff_subject_idxs)  



same_subject_ploidy_changes = []

same_subject_snp_changes = []
same_subject_snp_mutations = []
same_subject_snp_reversions = []

same_subject_gene_changes = []
same_subject_gene_gains = []
same_subject_gene_losses = []

num_temporal_samples = 0


for sample_pair_idx in xrange(0,len(same_subject_snp_idxs[0])):
#    
    snp_i = same_subject_snp_idxs[0][sample_pair_idx]
    snp_j = same_subject_snp_idxs[1][sample_pair_idx]
    
   
    if snp_samples[snp_i] in haploid_samples:
    
        if snp_samples[snp_j] in haploid_samples:
            ploidy_change = 'haploid->haploid'
        else:
            ploidy_change = 'haploid->diploid'
            
    else:
        if snp_samples[snp_j] in haploid_samples:
            ploidy_change = 'diploid->haploid'
        else:
            ploidy_change = 'diploid->diploid'
    
    num_temporal_samples += 1
    
    
    
    snp_changes = snp_difference_matrix[snp_i,snp_j]
    snp_mutations = snp_difference_matrix_mutation[snp_i,snp_j]
    snp_reversions = snp_difference_matrix_reversion[snp_i,snp_j]

    same_subject_ploidy_changes.append(ploidy_change)
    same_subject_snp_changes.append(snp_changes)
    same_subject_snp_mutations.append(snp_mutations)
    same_subject_snp_reversions.append(snp_reversions)
    
    i = same_subject_gene_idxs[0][sample_pair_idx]
    j = same_subject_gene_idxs[1][sample_pair_idx]
    
    if marker_coverages[i]<min_coverage or marker_coverages[j]<min_coverage:
        # can't look at gene changes!
        num_gene_changes = -1
        num_gene_gains = -1
        num_gene_losses = -1
        
    else:
    
        num_gene_changes = gene_difference_matrix[i,j]
        num_gene_gains = gene_gain_matrix[i,j]
        num_gene_losses = gene_loss_matrix[i,j]
        
    same_subject_gene_changes.append(num_gene_changes)
    same_subject_gene_gains.append(num_gene_gains)
    same_subject_gene_losses.append(num_gene_losses)
  
    if snp_changes > 0.5:
        
        sys.stderr.write("%s->%s: %d SNP changes (%d/%d), %d gene changes (%d/%d), %s\n" % (snp_samples[snp_i], snp_samples[snp_j], snp_changes, snp_mutations, snp_reversions, num_gene_changes, num_gene_gains, num_gene_losses, ploidy_change))
        
    
sys.stderr.write("%d total temporal samples\n" % num_temporal_samples)

sys.exit(0)

# clip lower bounds 
same_subject_gene_changes = numpy.clip(same_subject_gene_changes,0.3,1e09)
same_subject_snp_changes = numpy.clip(same_subject_snp_changes,0.3,1e09)

# Sort all lists by ascending lower bound on SNP changes, then gene changes
same_subject_snp_changes, same_subject_gene_changes, same_subject_snp_reversions, same_subject_snp_mutations, same_subject_gene_gains, same_subject_gene_losses, same_subject_ploidy_changes = (numpy.array(x) for x in zip(*sorted(zip(same_subject_snp_changes, same_subject_gene_changes, same_subject_snp_reversions, same_subject_snp_mutations, same_subject_gene_gains, same_subject_gene_losses, same_subject_ploidy_changes))))


diff_subject_snp_plowers = []
diff_subject_snp_puppers = []
diff_subject_gene_plowers = []
diff_subject_gene_puppers = []
between_host_gene_idx_map = {}
between_host_gene_idxs = []
low_divergence_between_host_gene_idxs = []
low_divergence_between_host_gene_idx_map = {}
for sample_pair_idx in xrange(0,len(diff_subject_snp_idxs[0])):
    
    snp_i = diff_subject_snp_idxs[0][sample_pair_idx]
    snp_j = diff_subject_snp_idxs[1][sample_pair_idx]
    
    if not ((snp_samples[snp_i] in haploid_samples) and (snp_samples[snp_j] in haploid_samples)):
        # both have to be haploids
        continue
        
    
    plower = snp_difference_matrix[snp_i,snp_j]
    pupper = plower*1.1
 
    #plower,pupper = stats_utils.calculate_poisson_rate_interval(snp_difference_matrix[snp_i,snp_j], snp_opportunity_matrix[snp_i, snp_j])
    
    diff_subject_snp_plowers.append(plower)
    diff_subject_snp_puppers.append(pupper)
        
    i = diff_subject_gene_idxs[0][sample_pair_idx]
    j = diff_subject_gene_idxs[1][sample_pair_idx]
    gene_differences = gene_diversity_utils.calculate_gene_differences_between(i, j, gene_depth_matrix, marker_coverages,include_high_copynum=False)

    if snp_substitution_rate[snp_i,snp_j] < clade_divergence_threshold:
    
        for gene_idx, depth_tuple_1, depth_tuple_2 in gene_differences:
            if gene_idx not in between_host_gene_idx_map:
                between_host_gene_idx_map[gene_idx]=0
            
            
            between_host_gene_idxs.append(gene_idx)
            
            between_host_gene_idx_map[gene_idx]+=1
    
    if snp_substitution_rate[snp_i,snp_j] < modification_divergence_threshold:
        # A modification, not a replacement!
        for gene_idx, depth_tuple_1, depth_tuple_2 in gene_differences:
            if gene_idx not in low_divergence_between_host_gene_idx_map:
                low_divergence_between_host_gene_idx_map[gene_idx]=0
            
            low_divergence_between_host_gene_idxs.append(gene_idx)
            low_divergence_between_host_gene_idx_map[gene_idx]+=1
    
    if marker_coverages[i]<min_coverage or marker_coverages[j]<min_coverage:
        # can't look at gene changes!
        plower = -1
        pupper = -1
    else:
        plower = gene_difference_matrix[i,j]
        pupper = plower*1.1
            
    #plower,pupper = stats_utils.calculate_poisson_rate_interval(gene_difference_matrix[i,j], gene_opportunity_matrix[i,j],alpha)
    
    diff_subject_gene_plowers.append(plower)
    diff_subject_gene_puppers.append(pupper)

# clip lower bounds 
diff_subject_gene_plowers = numpy.clip(diff_subject_gene_plowers,0.3,1e09)
diff_subject_snp_plowers = numpy.clip(diff_subject_snp_plowers,0.3,1e09)
# Sort all four lists by ascending lower bound on SNP changes, then gene changes
diff_subject_snp_plowers, diff_subject_gene_plowers, diff_subject_snp_puppers,  diff_subject_gene_puppers = (numpy.array(x) for x in zip(*sorted(zip(diff_subject_snp_plowers, diff_subject_gene_plowers, diff_subject_snp_puppers, diff_subject_gene_puppers))))



within_host_gene_sfs = []
within_host_gene_multiplicities = []
within_host_gene_prevalences = []
for gene_idx in within_host_gene_idx_map.keys():
    within_host_gene_sfs.append(within_host_gene_idx_map[gene_idx])
    for i in xrange(0, within_host_gene_idx_map[gene_idx]):
        within_host_gene_prevalences.append(prevalences[gene_idx])
        within_host_gene_multiplicities.append(within_host_gene_idx_map[gene_idx])

within_host_null_gene_prevalences = []
for gene_idx in within_host_null_gene_idxs:
    within_host_null_gene_prevalences.append(prevalences[gene_idx])    

within_host_gene_sfs.sort()
within_host_gene_sfs = numpy.array(within_host_gene_sfs)
within_host_gene_prevalences.sort()
within_host_gene_prevalences = numpy.array(within_host_gene_prevalences)
within_host_null_gene_prevalences.sort()
within_host_null_gene_prevalences = numpy.array(within_host_null_gene_prevalences)

within_host_gene_multiplicities.sort()
within_host_gene_multiplicities = numpy.array(within_host_gene_multiplicities)


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

# Bootstrap between-host multiplicities
between_host_gene_multiplicities = []
num_bootstraps = 1
for i in xrange(0,num_bootstraps):
    
    bootstrapped_gene_idxs = choice(between_host_gene_idxs,len(within_host_gene_idxs),replace=False)
    
    # Create idx map
    bootstrapped_gene_idx_map = {}
    for gene_idx in bootstrapped_gene_idxs:
        if gene_idx not in bootstrapped_gene_idx_map:
            bootstrapped_gene_idx_map[gene_idx]=0
            
        bootstrapped_gene_idx_map[gene_idx]+=1
        
    for gene_idx in bootstrapped_gene_idxs:
        between_host_gene_multiplicities.append(bootstrapped_gene_idx_map[gene_idx])
        
between_host_gene_multiplicities.sort()
between_host_gene_multiplicities = numpy.array(between_host_gene_multiplicities)
    

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


y = 0
for snp_changes, snp_mutations, snp_reversions, gene_changes, gene_gains, gene_losses, ploidy_change in zip(same_subject_snp_changes, same_subject_snp_mutations, same_subject_snp_reversions, same_subject_gene_changes, same_subject_gene_gains, same_subject_gene_losses, same_subject_ploidy_changes):

    y-=2
    
    #within_snp_axis.semilogy([y,y], [snp_plower,snp_pupper],'g-',linewidth=0.25)
    within_snp_axis.semilogy([y], [snp_changes],'g.',markersize=1.5)
        
    #within_gene_axis.semilogy([y,y], [gene_plower,gene_pupper], 'g-',linewidth=0.25)
    within_gene_axis.semilogy([y],[gene_changes],  'g.',markersize=1.5)

    print "Mutations=%g, Reversions=%g, Gains=%g, Losses=%g, %s" % (snp_mutations, snp_reversions, gene_gains, gene_losses, ploidy_change)

y-=4

within_snp_axis.semilogy([y,y],[1e-09,1e09],'-',linewidth=0.25,color='k')
within_gene_axis.semilogy([y,y],[1e-09,1e09],'-',linewidth=0.25,color='k')


within_snp_axis.set_xlim([y-0.2,0])
within_gene_axis.set_xlim([y-0.2,0])    

within_snp_axis.set_xticks([])
within_gene_axis.set_xticks([])

within_snp_axis.set_yticks([])
within_snp_axis.minorticks_off()

within_gene_axis.set_yticks([])
within_gene_axis.minorticks_off()

#within_snp_axis.set_yticklabels([])
#within_gene_axis.set_yticklabels([])




y=0    
    
for snp_plower, snp_pupper, gene_plower, gene_pupper in zip(diff_subject_snp_plowers, diff_subject_snp_puppers, diff_subject_gene_plowers, diff_subject_gene_puppers)[0:50]:

    y-=1
    
    snp_axis.semilogy([y],[snp_plower],'r.',linewidth=0.35,markersize=1.5)
    gene_axis.semilogy([y],[gene_plower],'r.',linewidth=0.35,markersize=1.5)

    
    #snp_axis.semilogy([y,y],[snp_plower,snp_pupper],'r-',linewidth=0.35)
    #gene_axis.semilogy([y,y],[gene_plower,gene_pupper],'r-',linewidth=0.35)

y-=4

snp_axis.set_xlim([y-1,0])
gene_axis.set_xlim([y-1,0])    

snp_axis.set_xticks([])
gene_axis.set_xticks([])

#snp_axis.legend(loc='upper right',frameon=False)

#labels = snp_axis.get_yticklabels()
#print labels[0].get_text()
#labels[0].set_text('0')
#snp_axis.set_yticklabels(labels)

prevalence_bins = numpy.linspace(0,1,11)
prevalence_locations = prevalence_bins[:-1]+(prevalence_bins[1]-prevalence_bins[0])/2

#h = numpy.histogram(within_host_null_gene_prevalences,bins=prevalence_bins)[0]
#prevalence_axis.plot(prevalence_locations, h*1.0/h.sum(),'k-',label='Random')

h = numpy.histogram(between_host_gene_prevalences,bins=prevalence_bins)[0]
prevalence_axis.plot(prevalence_locations, h*1.0/h.sum(),'r.-',label='Between-host',markersize=3)

if len(low_divergence_between_host_gene_prevalences) > 0:
    print low_divergence_between_host_gene_prevalences
    print low_divergence_between_host_gene_prevalences.mean()
    print len(low_divergence_between_host_gene_prevalences), len(between_host_gene_prevalences)
    
    h = numpy.histogram(low_divergence_between_host_gene_prevalences,bins=prevalence_bins)[0]
    prevalence_axis.plot(prevalence_locations, h*1.0/h.sum(),'r.-',label=('d<%g' % modification_divergence_threshold), alpha=0.5,markersize=3)

h = numpy.histogram(within_host_gene_prevalences,bins=prevalence_bins)[0]
prevalence_axis.plot(prevalence_locations, h*1.0/h.sum(),'g.-',label='Within-host',markersize=3)

print len(within_host_gene_prevalences), "within-host changes"

prevalence_axis.legend(loc='upper right',frameon=False,fontsize=4)

multiplicity_bins = numpy.arange(0,5)+0.5
multiplicity_locs = numpy.arange(1,5)

between_host_multiplicity_histogram = numpy.histogram(between_host_gene_multiplicities,bins=multiplicity_bins)[0]

within_host_multiplicity_histogram = numpy.histogram(within_host_gene_multiplicities,bins=multiplicity_bins)[0]


#xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(between_host_gene_multiplicities)
#parallelism_axis.step(xs+0.5,ns*1.0/ns[0],'r-',label='Between-host')

#xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(within_host_gene_multiplicities)
#parallelism_axis.step(xs+0.5,ns*1.0/ns[0],'g-',label='Within-host')


parallelism_axis.bar(multiplicity_locs, between_host_multiplicity_histogram*1.0/between_host_multiplicity_histogram.sum(), width=0.3,color='r',linewidth=0)

parallelism_axis.bar(multiplicity_locs-0.3, within_host_multiplicity_histogram*1.0/within_host_multiplicity_histogram.sum(), width=0.3,color='g',linewidth=0)

#prevalence_axis.set_ylim([0,0.6])

sys.stderr.write("Saving figure...\t")
fig.savefig('%s/supplemental_polyploid_changes%s.pdf' % (parse_midas_data.analysis_directory, other_species_str),bbox_inches='tight')
sys.stderr.write("Done!\n")

    
