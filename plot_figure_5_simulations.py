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
import core_gene_utils

import stats_utils
import matplotlib.colors as colors
import matplotlib.cm as cmx
from math import log10,ceil
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy.random import randint, choice

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
clade_divergence_threshold = 1e-02 # specifically tuned for B. vulgatus
modification_divergence_threshold = 3e-04 #the threshold for deciding when something is a modification vs a replacement. Like most other things, it is an arbitrary choice. in this case, specifically tuned for B. vulgatus. 

# Load subject and sample metadata
sys.stderr.write("Loading sample metadata...\n")
subject_sample_map = parse_HMP_data.parse_subject_sample_map()
sample_country_map = parse_HMP_data.parse_sample_country_map()
sample_order_map = parse_HMP_data.parse_sample_order_map()
sys.stderr.write("Done!\n")
       
# Only plot samples above a certain depth threshold that are "haploids"
snp_samples = diversity_utils.calculate_haploid_samples(species_name, debug=debug)

####################################################
#
# Set up Figure (6 panels, arranged in 3x2 grid)
#
####################################################

#pylab.figure(1,figsize=(3.43,3.1))
pylab.figure(1,figsize=(2.5,4.2))

fig = pylab.gcf()


# make 2 panels
#outer_grid  = gridspec.GridSpec(1,2, width_ratios=[1,1], wspace=0.4)

#differences_grid = gridspec.GridSpecFromSubplotSpec(3, 1, height_ratios=[1,1,1], subplot_spec=outer_grid[0], hspace=0.1)
                
#gene_grid = gridspec.GridSpecFromSubplotSpec(3, 1, height_ratios=[1,1,1], subplot_spec=outer_grid[1], hspace=0.5)

## Supp figure
pylab.figure(2,figsize=(3.42,0.9))
fig2 = pylab.gcf()

## Supp figure
pylab.figure(3,figsize=(1.7,0.9))
fig3 = pylab.gcf()


differences_grid = gridspec.GridSpec(3, 1, height_ratios=[1,1,1], hspace=0.1)
                
gene_grid = gridspec.GridSpec(1,2, width_ratios=[1,1], wspace=0.4)            
other_gene_grid = gridspec.GridSpec(1,1)




###################
#
# SNP change panel
#
###################

snp_axis = plt.Subplot(fig, differences_grid[0])
fig.add_subplot(snp_axis)

#snp_axis.set_title("%s %s (%s)" % tuple(species_name.split("_")),fontsize=7)

snp_axis.set_ylabel('SNV changes',labelpad=2)
snp_axis.set_ylim([1.2e-01,1e05])

snp_axis.semilogy([1e-09,1e-09],[1e08,1e08],'b.',markersize=1.5,label='Within-host')
snp_axis.semilogy([1e-09,1e-09],[1e08,1e08],'r.',linewidth=0.35,markersize=1.5,label='Between-host (lowest)')


snp_axis.spines['top'].set_visible(False)
snp_axis.spines['right'].set_visible(False)
snp_axis.get_xaxis().tick_bottom()
snp_axis.get_yaxis().tick_left()

#snp_axis.set_title('SNV changes',fontsize=5)

within_snp_axis = snp_axis 

###################
#
# Gene change panel
#
###################

gene_loss_axis = plt.Subplot(fig, differences_grid[1])
fig.add_subplot(gene_loss_axis)

gene_loss_axis.set_ylabel('Gene losses',labelpad=2)
#gene_loss_axis.set_title('Gene losses', fontsize=5,y=0.9)
gene_loss_axis.set_ylim([1.2e-01,1e04])

gene_loss_axis.spines['top'].set_visible(False)
gene_loss_axis.spines['right'].set_visible(False)
gene_loss_axis.get_xaxis().tick_bottom()
gene_loss_axis.get_yaxis().tick_left()

gene_gain_axis = plt.Subplot(fig, differences_grid[2])
fig.add_subplot(gene_gain_axis)

gene_gain_axis.set_ylabel('Gene gains',labelpad=2)
#gene_gain_axis.set_title('Gene gains', fontsize=5,y=0.9)
gene_gain_axis.set_ylim([1.2e-01,1e04])

gene_gain_axis.set_xlabel('Sample pairs')

gene_gain_axis.spines['top'].set_visible(False)
gene_gain_axis.spines['right'].set_visible(False)
gene_gain_axis.get_xaxis().tick_bottom()
gene_gain_axis.get_yaxis().tick_left()

snp_axis.fill_between([-1e06,1e06],[1e-01,1e-01],[0.6,0.6],color='0.8')
gene_loss_axis.fill_between([-1e06,1e06],[1e-01,1e-01],[0.6,0.6],color='0.8')
gene_gain_axis.fill_between([-1e06,1e06],[1e-01,1e-01],[0.6,0.6],color='0.8')

##############################################################################
#
# Gene change prevalence panel
#
##############################################################################

prevalence_axis = plt.Subplot(fig2, gene_grid[0])
fig2.add_subplot(prevalence_axis)

prevalence_axis.set_ylabel('Fraction genes $\leq p$',labelpad=2)
prevalence_axis.set_xlabel('Prevalence, $p$',labelpad=2)
prevalence_axis.set_xlim([0,1.05])
#prevalence_axis.set_ylim([0,1.1])

prevalence_axis.spines['top'].set_visible(False)
prevalence_axis.spines['right'].set_visible(False)
prevalence_axis.get_xaxis().tick_bottom()
prevalence_axis.get_yaxis().tick_left()

##############################################################################
#
# Gene linkage panel
#
##############################################################################

linkage_axis = plt.Subplot(fig2, gene_grid[1])
fig2.add_subplot(linkage_axis)

linkage_axis.set_ylabel('Fraction of genes $\geq c$',labelpad=2)
linkage_axis.set_xlabel('Neighboring fold change, $c$',labelpad=2)
linkage_axis.set_xlim([1,10])
linkage_axis.set_ylim([0,1])
#prevalence_axis.set_ylim([0,1.1])

linkage_axis.spines['top'].set_visible(False)
linkage_axis.spines['right'].set_visible(False)
linkage_axis.get_xaxis().tick_bottom()
linkage_axis.get_yaxis().tick_left()


##############################################################################
#
# Gene change multiplicity panel
#
##############################################################################

#multiplicity_axis = plt.Subplot(fig, gene_grid[1])
#fig.add_subplot(multiplicity_axis)
#multiplicity_axis = plt.Subplot(supplemental_fig, supplemental_outer_grid[0])
#supplemental_fig.add_subplot(multiplicity_axis)


#multiplicity_axis.set_ylabel('Fraction gene changes',labelpad=2)
#multiplicity_axis.set_xlabel('Gene multiplicity, $m$',labelpad=2)
#multiplicity_axis.set_xlim([0.5,3.5])
#multiplicity_axis.set_ylim([0,1.05])

#multiplicity_axis.set_xticks([1,2,3])

#multiplicity_axis.spines['top'].set_visible(False)
#multiplicity_axis.spines['right'].set_visible(False)
#multiplicity_axis.get_xaxis().tick_bottom()
#multiplicity_axis.get_yaxis().tick_left()

##############################################################################
#
# Gene change parallelism panel
#
##############################################################################

parallelism_axis = plt.Subplot(fig3, other_gene_grid[0])
fig3.add_subplot(parallelism_axis)
#parallelism_axis = plt.Subplot(supplemental_fig, supplemental_outer_grid[0])
#supplemental_fig.add_subplot(parallelism_axis)


parallelism_axis.set_ylabel('Fraction genes $\geq c$',labelpad=2)
parallelism_axis.set_xlabel('Largest parallel fold-change, $c$',labelpad=2)
#parallelism_axis.set_xlim([0,5])
parallelism_axis.set_ylim([0,1.1])

parallelism_axis.spines['top'].set_visible(False)
parallelism_axis.spines['right'].set_visible(False)
parallelism_axis.get_xaxis().tick_bottom()
parallelism_axis.get_yaxis().tick_left()


#####################
# Analyze the data  #
#####################

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

sys.stderr.write("Loading whitelisted genes...\n")
non_shared_genes = core_gene_utils.parse_non_shared_reference_genes(species_name)
shared_pangenome_genes = core_gene_utils.parse_shared_genes(species_name)
sys.stderr.write("Done! %d shared genes and %d non-shared genes\n" % (len(shared_pangenome_genes), len(non_shared_genes)))


# Load gene coverage information for species_name
sys.stderr.write("Loading pangenome data for %s...\n" % species_name)
gene_samples, gene_names, gene_presence_matrix, gene_depth_matrix, marker_coverages, gene_reads_matrix = parse_midas_data.parse_pangenome_data(species_name,allowed_samples=snp_samples,disallowed_genes=shared_pangenome_genes)
gene_names = list(gene_names)
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

# indexes for different subject pairs
diff_subject_snp_idxs = parse_midas_data.apply_sample_index_map_to_indices(snp_sample_idx_map, desired_diff_subject_idxs)  
diff_subject_gene_idxs = parse_midas_data.apply_sample_index_map_to_indices(gene_sample_idx_map, desired_diff_subject_idxs)  

# Calculate median between-host differences
#modification_divergence_threshold = numpy.median(snp_substitution_rate[diff_subject_snp_idxs])/4.0

median_between_host_nucleotide_changes = numpy.median(snp_difference_matrix[diff_subject_snp_idxs])
median_between_host_gene_gains = numpy.median(gene_gain_matrix[diff_subject_gene_idxs])
median_between_host_gene_losses = numpy.median(gene_loss_matrix[diff_subject_gene_idxs])


# Calculate subset of "modification timepoints" 
modification_pair_idxs = set([])

for sample_pair_idx in xrange(0,len(same_subject_snp_idxs[0])):
#    
    snp_i = same_subject_snp_idxs[0][sample_pair_idx]
    snp_j = same_subject_snp_idxs[1][sample_pair_idx]
    #
    if snp_substitution_rate[snp_i, snp_j] < modification_divergence_threshold:
        modification_pair_idxs.add( sample_pair_idx ) 



between_host_gene_idxs = [] # indexes of genes that changed between hosts
low_divergence_between_host_gene_idxs = [] # indexes of genes that changed between particularly low divergence hosts

# Store the total amount of SNP and gene changes in these arrays. 
diff_subject_snp_changes = []
diff_subject_gene_changes = []
diff_subject_gene_gains = []
diff_subject_gene_losses = []

for sample_pair_idx in xrange(0,len(diff_subject_snp_idxs[0])):
    #
    snp_i = diff_subject_snp_idxs[0][sample_pair_idx]
    snp_j = diff_subject_snp_idxs[1][sample_pair_idx]
    #
    diff_subject_snp_changes.append( snp_difference_matrix[snp_i, snp_j] )
    #
    # Now do genes
    i = diff_subject_gene_idxs[0][sample_pair_idx]
    j = diff_subject_gene_idxs[1][sample_pair_idx]
    #
    if (marker_coverages[i]<min_coverage) or (marker_coverages[j]<min_coverage):
        diff_subject_gene_changes.append( -1 )
    else:
        diff_subject_gene_changes.append( gene_difference_matrix[i,j] )
        diff_subject_gene_gains.append( gene_gain_matrix[i,j] )
        diff_subject_gene_losses.append( gene_loss_matrix[i,j] )
        #
        #
        # why are we conditionng on the substitution rate being less than this threshold? NRG
        if snp_substitution_rate[snp_i, snp_j] < clade_divergence_threshold:
            #
            # Now actually calculate genes that differ! 
            gene_idxs = gene_diversity_utils.calculate_gene_differences_between_idxs(i,j, gene_reads_matrix, gene_depth_matrix, marker_coverages)
            #
            between_host_gene_idxs.extend(gene_idxs)
            #
            if snp_substitution_rate[snp_i, snp_j] < modification_divergence_threshold:
                low_divergence_between_host_gene_idxs.extend(gene_idxs)

diff_subject_snp_changes = numpy.array(diff_subject_snp_changes)
diff_subject_gene_changes = numpy.array(diff_subject_gene_changes)
diff_subject_gene_gains = numpy.array(diff_subject_gene_gains)
diff_subject_gene_losses = numpy.array(diff_subject_gene_losses)

# What does this sort do? NRG
diff_subject_snp_changes, diff_subject_gene_changes, diff_subject_gene_gains, diff_subject_gene_losses =  (numpy.array(x) for x in zip(*sorted(zip(diff_subject_snp_changes, diff_subject_gene_changes, diff_subject_gene_gains, diff_subject_gene_losses))))


#these arrays will store the number of SNP and gene changes for different same-subject sample pairs. 

same_subject_snp_changes = []
same_subject_snp_mutations = []
same_subject_snp_reversions = []

same_subject_gene_changes = []
same_subject_gene_gains = []
same_subject_gene_losses = []

within_host_gene_idxs = [] # the indexes of genes that actually changed between samples
within_host_null_gene_idxs = [] # a null distribution of gene indexes. chosen randomly from genes "present" in genome
# NRG: why are these null genes chosen?

# NRG what does this store?
within_host_next_fold_changes = []
within_host_null_next_fold_changes = []
within_host_between_next_fold_changes = []

within_host_neighbor_fold_changes = []
within_host_null_neighbor_fold_changes = []
within_host_between_neighbor_fold_changes = []



# sizes of estimated blocks of gene changes within hosts
within_host_blockss = []

total_modification_error_rate = 0
total_modifications = 0

# Calculate the gene changes in each host
for sample_pair_idx in xrange(0,len(same_subject_snp_idxs[0])):
#    
    snp_i = same_subject_snp_idxs[0][sample_pair_idx]
    snp_j = same_subject_snp_idxs[1][sample_pair_idx]
    #
    sample_i = snp_samples[snp_i]
    sample_j = snp_samples[snp_j]
    #    
    L, perr, mutations, reversions = calculate_temporal_changes.calculate_mutations_reversions_from_temporal_change_map(temporal_change_map, sample_i, sample_j) 
    perr = L*perr   
    #
    num_mutations = len(mutations)
    num_reversions = len(reversions)
    num_snp_changes = num_mutations+num_reversions
    #
    if perr>1: # perr is probably a bad name, because it's the false positive rate summed across the entire genome (i.e., the expected number of errors across the genome, which can definitely be >1 if per site error rates are not tiny).  
        num_mutations = 0
        num_reversions = 0
        num_snp_changes = -1
    #
    same_subject_snp_changes.append(num_snp_changes)
    same_subject_snp_mutations.append(num_mutations)
    same_subject_snp_reversions.append(num_reversions)
    #
    #
    #Including replacement events in within-host evolution statistics messes everything up, so we want to only focus on changes that are not obvious replacements (i.e., modifications):
    if snp_substitution_rate[snp_i, snp_j] < modification_divergence_threshold:
        if perr<=1:
            total_modification_error_rate += perr # why are we calculating this? NRG
            total_modifications += num_snp_changes
    #       
    if num_snp_changes>-1:
        print sample_i, sample_j, num_snp_changes, num_mutations, num_reversions, perr
    #
    i = same_subject_gene_idxs[0][sample_pair_idx]
    j = same_subject_gene_idxs[1][sample_pair_idx]
    #
    if marker_coverages[i]<min_coverage or marker_coverages[j]<min_coverage:
        # can't look at gene changes!
        #
        same_subject_gene_changes.append(-1)
        same_subject_gene_gains.append(-1)
        same_subject_gene_losses.append(-1)
        gene_perr = 1
    else:
        # where is the output for this used? NRG
        gene_L, gene_perr, gains, losses = calculate_temporal_changes.calculate_gains_losses_from_temporal_change_map(temporal_change_map, sample_i, sample_j)
        
        if gene_L>0:
            gene_perr = gene_L*gene_perr
            #
        #
        print sample_i, sample_j, gene_difference_matrix[i,j], gene_gain_matrix[i,j], gene_loss_matrix[i,j], gene_perr
        #
        same_subject_gene_changes.append(gene_difference_matrix[i,j])
        same_subject_gene_gains.append(gene_gain_matrix[i,j])
        same_subject_gene_losses.append(gene_loss_matrix[i,j])
        #
        # Only include samples that are deemed to be modifications 
        # (don't include replacements)
        if sample_pair_idx in modification_pair_idxs and (gene_perr>-0.5) and (gene_perr<1):
            #
            all_gene_changes = gains+losses
            #
            gene_blocks = gene_diversity_utils.merge_nearby_gene_differences(all_gene_changes)
            #
            within_host_blockss.append(gene_blocks)
            #
            #
            # Calculate set of genes that are present in at least one sample
            present_gene_idxs = []
            present_gene_idxs.extend( numpy.nonzero( (gene_copynum_matrix[:,i]>0.5)*(gene_copynum_matrix[:,i]<2))[0] )
            present_gene_idxs.extend( numpy.nonzero( (gene_copynum_matrix[:,j]>0.5)*(gene_copynum_matrix[:,j]<2))[0] )
            #
            pair_specific_gene_idxs = gene_diversity_utils.calculate_gene_differences_between_idxs(i, j, gene_reads_matrix, gene_depth_matrix, marker_coverages)
            #
            if len(pair_specific_gene_idxs)==0:
                continue
            #
            # Ultimately, gene changes within and between hosts are compared with randomly chosen genes (hence: null)
            pair_specific_null_gene_idxs = choice(present_gene_idxs, len(pair_specific_gene_idxs)*10 )
            pair_specific_between_gene_idxs = choice(between_host_gene_idxs, len(pair_specific_gene_idxs)*10 )
            #
            within_host_gene_idxs.extend(pair_specific_gene_idxs)
            within_host_null_gene_idxs.extend(pair_specific_null_gene_idxs)                       #
            #
            #####
            # Calculate copynum fold change of neighboring genes
            neighbor_fold_changes = []
            null_neighbor_fold_changes = []
            between_neighbor_fold_changes = []
            #
            for fold_changes, gene_idxs in zip([neighbor_fold_changes, null_neighbor_fold_changes, between_neighbor_fold_changes], [within_host_gene_idxs, within_host_null_gene_idxs,pair_specific_between_gene_idxs]):
                #
                neighboring_gene_idxs = []
                for gene_idx in gene_idxs:
                    neighboring_gene_idxs.extend( gene_diversity_utils.get_nearby_gene_idxs(gene_names, gene_idx) )
                #
                neighboring_gene_idxs = numpy.array(neighboring_gene_idxs)
                #
                # calculate log-fold change
                logfoldchanges = numpy.fabs( numpy.log2(clipped_gene_copynum_matrix[neighboring_gene_idxs,j] / clipped_gene_copynum_matrix[neighboring_gene_idxs,i] ) )
                # check if the copy num of either i or j is eithin the accepted range of 0.5 or 3. 
                good_idxs = numpy.logical_or( good_copynum_matrix[neighboring_gene_idxs, i], good_copynum_matrix[neighboring_gene_idxs, j] ) 
                good_idxs *= numpy.logical_and( low_copynum_matrix[neighboring_gene_idxs, i], low_copynum_matrix[neighboring_gene_idxs, j] ) 
                #
                bad_idxs = numpy.logical_not( good_idxs ) 
                logfoldchanges[bad_idxs] = -1
                #
                fold_changes.extend(logfoldchanges[logfoldchanges>-0.5])
                #   
            within_host_neighbor_fold_changes.extend(neighbor_fold_changes)
            within_host_null_neighbor_fold_changes.extend(null_neighbor_fold_changes)
            within_host_between_neighbor_fold_changes.extend(between_neighbor_fold_changes)
            #
            other_fold_changes = []
            null_other_fold_changes = []
            between_other_fold_changes = []
            #
            #How parallel are gene changes across hosts? 
            #This next section is calculating the next-most extreme log-fold change, so the inner for loop is trying to figure out what the next-most extreme event is. 
            for other_sample_pair_idx in xrange(0,len(same_subject_snp_idxs[0])):
                #
                other_i = same_subject_gene_idxs[0][other_sample_pair_idx]
                other_j = same_subject_gene_idxs[1][other_sample_pair_idx]
                #
                # Make sure we don't count the same thing twice! 
                if other_sample_pair_idx == sample_pair_idx:
                    continue
                #
                # Make sure it is not a replacement!    
                if other_sample_pair_idx not in modification_pair_idxs:
                    continue
                # dont' inculde stuff where coverage is low.    
                if (marker_coverages[other_i]<min_coverage) or (marker_coverages[other_j]<min_coverage):
                    continue     
                #
                # calculate log-fold change
                logfoldchanges = numpy.fabs( numpy.log2(clipped_gene_copynum_matrix[pair_specific_gene_idxs,other_j] / clipped_gene_copynum_matrix[pair_specific_gene_idxs,other_i] ) )
                # check if the copy num of either i or j is eithin the accepted range of 0.5 or 3. 
                good_idxs = numpy.logical_or( good_copynum_matrix[pair_specific_gene_idxs, other_i], good_copynum_matrix[pair_specific_gene_idxs, other_j] ) 
                good_idxs *= numpy.logical_and( low_copynum_matrix[pair_specific_gene_idxs, other_i], low_copynum_matrix[pair_specific_gene_idxs, other_j] ) 
                #
                bad_idxs = numpy.logical_not( good_idxs ) 
                logfoldchanges[bad_idxs] = -1
                #
                other_fold_changes.append(logfoldchanges)
                #
                # calculate log-fold change
                logfoldchanges = numpy.fabs( numpy.log2(clipped_gene_copynum_matrix[pair_specific_null_gene_idxs,other_j] / clipped_gene_copynum_matrix[pair_specific_null_gene_idxs,other_i] ) )
                #
                # only include genes that are at low copynum at both timepoints # Why?? NRG
                # and have a "good" copynum at at least one point
                good_idxs = numpy.logical_or( good_copynum_matrix[pair_specific_null_gene_idxs, other_i], good_copynum_matrix[pair_specific_null_gene_idxs, other_j] ) 
                good_idxs *= numpy.logical_and( low_copynum_matrix[pair_specific_null_gene_idxs, other_i], low_copynum_matrix[pair_specific_null_gene_idxs, other_j] ) 
                bad_idxs = numpy.logical_not( good_idxs ) 
                #
                logfoldchanges[bad_idxs] = -1
                #
                null_other_fold_changes.append( logfoldchanges )
                #
                # calculate log-fold change
                logfoldchanges = numpy.fabs( numpy.log2(clipped_gene_copynum_matrix[pair_specific_between_gene_idxs,other_j] / clipped_gene_copynum_matrix[pair_specific_between_gene_idxs,other_i] ) )
                good_idxs = numpy.logical_or( good_copynum_matrix[pair_specific_between_gene_idxs, other_i], good_copynum_matrix[pair_specific_between_gene_idxs, other_j] )                # 
                good_idxs *= numpy.logical_and( low_copynum_matrix[pair_specific_between_gene_idxs, other_i], low_copynum_matrix[pair_specific_between_gene_idxs, other_j] )                #
                # 
                bad_idxs = numpy.logical_not( good_idxs ) 
                logfoldchanges[bad_idxs] = -1
                #
                between_other_fold_changes.append( logfoldchanges )
                #
            other_fold_changes = numpy.array(other_fold_changes)
            null_other_fold_changes = numpy.array(null_other_fold_changes)
            between_other_fold_changes = numpy.array(between_other_fold_changes)
            #
            #
            # Pull out the largest fold change computed across all pairs of hosts other than the i and j of interest at the top of the for loop
            for gene_idx in xrange(0,other_fold_changes.shape[1]):
                fold_changes = other_fold_changes[:,gene_idx]
                fold_changes = fold_changes[fold_changes>-0.5]
                if len(fold_changes)>0:
                    #print "Observed biggest change: %g, median change %g" % (fold_changes.max(), numpy.median(fold_changes))
                    #print fold_changes
                    within_host_next_fold_changes.append( (fold_changes).max() )
                    #
            for gene_idx in xrange(0,null_other_fold_changes.shape[1]):
                fold_changes = null_other_fold_changes[:,gene_idx]
                fold_changes = fold_changes[fold_changes>-0.5]
                if len(fold_changes)>0:
                    within_host_null_next_fold_changes.append( (fold_changes).max() )
                #
            for gene_idx in xrange(0, between_other_fold_changes.shape[1]):
                fold_changes = between_other_fold_changes[:,gene_idx]
                fold_changes = fold_changes[fold_changes>-0.5]
                if len(fold_changes)>0:
                    within_host_between_next_fold_changes.append( (fold_changes).max() )


within_host_next_fold_changes = numpy.array(within_host_next_fold_changes)
within_host_between_next_fold_changes = numpy.array(within_host_between_next_fold_changes)
within_host_null_next_fold_changes = numpy.array(within_host_null_next_fold_changes)

within_host_next_fold_changes = numpy.power(2, within_host_next_fold_changes)
within_host_between_next_fold_changes = numpy.power(2, within_host_between_next_fold_changes)
within_host_null_next_fold_changes = numpy.power(2, within_host_null_next_fold_changes)

within_host_neighbor_fold_changes = numpy.power(2, within_host_neighbor_fold_changes)
within_host_null_neighbor_fold_changes = numpy.power(2, within_host_null_neighbor_fold_changes)

within_host_between_neighbor_fold_changes = numpy.power(2, within_host_between_neighbor_fold_changes)


print "%g modifications, %g expected" % (total_modifications, total_modification_error_rate)


# Sort all lists by ascending lower bound on SNP changes, then gene changes
same_subject_snp_changes, same_subject_gene_changes, same_subject_snp_reversions, same_subject_snp_mutations, same_subject_gene_gains, same_subject_gene_losses = (numpy.array(x) for x in zip(*sorted(zip(same_subject_snp_changes, same_subject_gene_changes, same_subject_snp_reversions, same_subject_snp_mutations, same_subject_gene_gains, same_subject_gene_losses))))

# Calculate distribution of num blocks
within_host_num_blocks = []
within_host_block_sizes = []
for within_host_blocks in within_host_blockss:
    within_host_num_blocks.append(len(within_host_blocks))
    for block in within_host_blocks:
        within_host_block_sizes.append(len(block))
within_host_num_blocks = numpy.array(within_host_num_blocks)
within_host_num_blocks.sort()

within_host_block_sizes = numpy.array(within_host_block_sizes)
within_host_block_sizes.sort()


# Construct site frequency spectra

within_host_gene_idx_counts = {}
for gene_idx in within_host_gene_idxs:
    if gene_idx not in within_host_gene_idx_counts:
        within_host_gene_idx_counts[gene_idx] = 0
    within_host_gene_idx_counts[gene_idx] += 1

within_host_gene_sfs = list(within_host_gene_idx_counts.values())

within_host_gene_prevalences = []
within_host_gene_multiplicities = []    
for gene_idx in within_host_gene_idxs:    
    within_host_gene_multiplicities.append( within_host_gene_idx_counts[gene_idx] )
    within_host_gene_prevalences.append(prevalences[gene_idx])

        
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

# Calculate counts of between host gene changes
between_host_gene_idx_counts = {}
for gene_idx in between_host_gene_idxs:
    if gene_idx not in between_host_gene_idx_counts:
        between_host_gene_idx_counts[gene_idx] = 0
    between_host_gene_idx_counts[gene_idx] += 1


between_host_gene_sfs = list(between_host_gene_idx_counts.values())

between_host_gene_prevalences = []
for gene_idx in between_host_gene_idxs:
    between_host_gene_prevalences.append(prevalences[gene_idx])
    
between_host_gene_prevalences.sort()
between_host_gene_prevalences = numpy.array(between_host_gene_prevalences)

# Bootstrap between-host multiplicities
between_host_gene_multiplicities = []
between_host_gene_sfs = []
num_bootstraps = 1
for i in xrange(0,num_bootstraps):
#    
    bootstrapped_gene_idxs = choice(between_host_gene_idxs,len(within_host_gene_idxs),replace=False)
#    
    # Create idx map
    bootstrapped_gene_idx_map = {}
    for gene_idx in bootstrapped_gene_idxs:
        if gene_idx not in bootstrapped_gene_idx_map:
            bootstrapped_gene_idx_map[gene_idx]=0
#            
        bootstrapped_gene_idx_map[gene_idx]+=1
#        
    for gene_idx in bootstrapped_gene_idxs:
        between_host_gene_multiplicities.append( bootstrapped_gene_idx_map[gene_idx] )
#    
    between_host_gene_sfs.extend( bootstrapped_gene_idx_map.values() )    
    
        
between_host_gene_multiplicities.sort()
between_host_gene_multiplicities = numpy.array(between_host_gene_multiplicities)

low_divergence_between_host_gene_prevalences = []
for gene_idx in low_divergence_between_host_gene_idxs:
    low_divergence_between_host_gene_prevalences.append(prevalences[gene_idx])

low_divergence_between_host_gene_prevalences.sort()
low_divergence_between_host_gene_prevalences = numpy.array(low_divergence_between_host_gene_prevalences)

# Done calculating... now plot figure!

# First plot within-host changes
y = 0
for snp_changes, snp_mutations, snp_reversions, gene_changes, gene_gains, gene_losses in zip(same_subject_snp_changes, same_subject_snp_mutations, same_subject_snp_reversions, same_subject_gene_changes, same_subject_gene_gains, same_subject_gene_losses):

    if snp_changes>-0.5 and snp_changes<0.5:
        snp_changes = 0.3
    
    if gene_changes>-0.5 and gene_changes<0.5:
        gene_changes = 0.3
    
    if gene_changes>-0.5 and gene_losses<0.5:
        gene_losses = 0.3
    
    if gene_changes>-0.5 and gene_gains<0.5:
        gene_gains = 0.3
    

    y-=2
    
    within_snp_axis.semilogy([y], [snp_changes],'b.',markersize=3,zorder=1)
        
    gene_loss_axis.semilogy([y],[gene_losses],  'b.',markersize=3,zorder=1)
    gene_gain_axis.semilogy([y],[gene_gains],'b.',markersize=3,zorder=1)
    
    print "Mutations=%g, Reversions=%g, Gains=%g, Losses=%g" % (snp_mutations, snp_reversions, gene_gains, gene_losses)

ymin = y-3
ymax = 3

snp_axis.set_xlim([ymin,ymax])
gene_loss_axis.set_xlim([ymin,ymax])    
gene_gain_axis.set_xlim([ymin,ymax])    
snp_axis.set_xticks([])
gene_loss_axis.set_xticks([])
gene_gain_axis.set_xticks([])

# Plot typical between-host values
line, = snp_axis.semilogy([ymin,ymax],[median_between_host_nucleotide_changes, median_between_host_nucleotide_changes],'r-',zorder=0,linewidth=0.25,alpha=0.5)
#line.set_dashes((0.5,0.5))
line, = gene_loss_axis.semilogy([ymin,ymax],[median_between_host_gene_losses, median_between_host_gene_losses],'r-',zorder=0,alpha=0.5)
#line.set_dashes((0.5,0.5))
line, = gene_gain_axis.semilogy([ymin,ymax],[median_between_host_gene_gains, median_between_host_gene_gains],'r-',zorder=0,alpha=0.5)
#line.set_dashes((0.5,0.5))

num_to_plot = 50
dy = (ymax-ymin)*1.0/num_to_plot

y=ymax    
    
for snp_changes, gene_changes, gene_gains, gene_losses in zip(diff_subject_snp_changes, diff_subject_gene_changes, diff_subject_gene_gains, diff_subject_gene_losses)[0:num_to_plot]:

    
    if snp_changes>-0.5 and snp_changes<0.5:
        snp_changes = 0.3
    
    if gene_changes>-0.5 and gene_changes<0.5:
        gene_changes = 0.3
    

    y-=dy
    
    snp_axis.semilogy([y],[snp_changes],'r.',linewidth=0.35,markersize=1.5,zorder=0)
    gene_loss_axis.semilogy([y], [gene_losses],'r.',linewidth=0.35,markersize=1.5,zorder=0)
    gene_gain_axis.semilogy([y], [gene_gains],'r.',linewidth=0.35,markersize=1.5,zorder=0)



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
#prevalence_axis.plot(prevalence_locations, h*1.0/h.sum(),'r.-',label='Between-host',markersize=3)

if len(low_divergence_between_host_gene_prevalences) > 0:
    print low_divergence_between_host_gene_prevalences
    print low_divergence_between_host_gene_prevalences.mean()
    print len(low_divergence_between_host_gene_prevalences), len(between_host_gene_prevalences)
    
    h = numpy.histogram(low_divergence_between_host_gene_prevalences,bins=prevalence_bins)[0]
    #prevalence_axis.plot(prevalence_locations, h*1.0/h.sum(),'r.-',label=('d<%g' % modification_divergence_threshold), alpha=0.5,markersize=3)

h = numpy.histogram(within_host_gene_prevalences,bins=prevalence_bins)[0]
#prevalence_axis.plot(prevalence_locations, h*1.0/h.sum(),'b.-',label='Within-host',markersize=3)

print len(within_host_gene_prevalences), "within-host changes"

# CDF version

xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(within_host_gene_prevalences)
prevalence_axis.step(xs,1-ns*1.0/ns[0],'b-',label='Within-host',zorder=2)

xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(between_host_gene_prevalences)
prevalence_axis.step(xs,1-ns*1.0/ns[0],'r-',label='Between-host',zorder=1)

xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(within_host_null_gene_prevalences)
prevalence_axis.step(xs,1-ns*1.0/ns[0],'k-',label='Random',zorder=0)

prevalence_axis.set_ylim([0,1.1])
prevalence_axis.set_xlim([0,1.05])



multiplicity_bins = numpy.arange(0,5)+0.5
multiplicity_locs = numpy.arange(1,5)

between_host_multiplicity_histogram = numpy.histogram(between_host_gene_multiplicities,bins=multiplicity_bins)[0]

within_host_multiplicity_histogram = numpy.histogram(within_host_gene_multiplicities,bins=multiplicity_bins)[0]



#multiplicity_axis.bar(multiplicity_locs, between_host_multiplicity_histogram*1.0/between_host_multiplicity_histogram.sum(), width=0.3,color='r',linewidth=0)

#multiplicity_axis.bar(multiplicity_locs-0.3, within_host_multiplicity_histogram*1.0/within_host_multiplicity_histogram.sum(), width=0.3,color='g',linewidth=0)


#prevalence_axis.set_ylim([0,0.6])

xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(within_host_next_fold_changes)
parallelism_axis.step(xs,ns*1.0/ns[0],'b-',label='Within-host',zorder=2)

xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(within_host_between_next_fold_changes)
parallelism_axis.step(xs,ns*1.0/ns[0],'r-',label='Between-host (lowest)',zorder=1)

xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(within_host_null_next_fold_changes)
parallelism_axis.step(xs,ns*1.0/ns[0],'k-',label='Random',zorder=0)


#parallelism_axis.legend(loc='upper right',frameon=False,fontsize=4)

#snp_axis.legend(loc='upper right',frameon=False,fontsize=4, numpoints=1)

snp_axis.legend(loc=(0.05,0.9),frameon=False,fontsize=5, ncol=2, numpoints=1,handlelength=1)



parallelism_axis.semilogx([1],[-1],'k.')
parallelism_axis.set_xlim([1,10])

# Plot block size distribution
#xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(within_host_block_sizes)
#linkage_axis.step(xs,ns*1.0/ns[0],'b-',zorder=0)


xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(within_host_neighbor_fold_changes)
linkage_axis.step(xs,ns*1.0/ns[0],'b-',label='Within-host',zorder=2)

xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(within_host_between_neighbor_fold_changes)
linkage_axis.step(xs,ns*1.0/ns[0],'r-',label='Between-host',zorder=1)


xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(within_host_null_neighbor_fold_changes)
linkage_axis.step(xs,ns*1.0/ns[0],'k-',label='Random',zorder=0)

linkage_axis.semilogx([1],[-1],'k.')

linkage_axis.legend(loc='upper right',frameon=False,fontsize=4)


sys.stderr.write("Saving figure...\t")
fig.savefig('%s/figure_5%s.pdf' % (parse_midas_data.analysis_directory, other_species_str),bbox_inches='tight')
fig2.savefig('%s/supplemental_figure_5%s.pdf' % (parse_midas_data.analysis_directory, other_species_str),bbox_inches='tight')
fig3.savefig('%s/supplemental_figure_5c%s.pdf' % (parse_midas_data.analysis_directory, other_species_str),bbox_inches='tight')

sys.stderr.write("Done!\n")

    
