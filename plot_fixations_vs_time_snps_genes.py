import matplotlib  
matplotlib.use('Agg') 
import parse_midas_data
import pylab
import sys
import numpy
import diversity_utils
import gene_diversity_utils
import stats_utils
import os


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

#############################################################################

min_change = 0.8
# Minimum median coverage of sample to look at
min_coverage = 20

#################
# Load metadata #
#################

# Load subject and sample metadata
sys.stderr.write("Loading HMP metadata...\n")
subject_sample_map = parse_midas_data.parse_subject_sample_map()
sys.stderr.write("Done!\n")
    
# Load time metadata
subject_sample_time_map = parse_midas_data.parse_subject_sample_time_map()

######################
# Load coverage data #
######################

# Load genomic coverage distributions
sample_coverage_histograms, samples = parse_midas_data.parse_coverage_distribution(species_name)
median_coverages = numpy.array([stats_utils.calculate_median_from_histogram(sample_coverage_histogram) for sample_coverage_histogram in sample_coverage_histograms])
sample_coverage_map = {samples[i]: median_coverages[i] for i in xrange(0,len(samples))}
    
   

###############################################################
# Compute Pi within patients to figure out which are haploid  #
###############################################################

# Load pi information for species_name
sys.stderr.write("Loading within-sample diversity for %s...\n" % species_name)
samples, total_pis, total_pi_opportunities = parse_midas_data.parse_within_sample_pi(species_name, debug)
sys.stderr.write("Done!\n")
pis = total_pis/total_pi_opportunities

median_coverages = numpy.array([sample_coverage_map[samples[i]] for i in xrange(0,len(samples))])

###############################################################
# Indexes for SNP samples that have low piS and high coverage #
###############################################################

# Only plot samples above a certain depth threshold that are "haploids"
snp_samples = samples[(median_coverages>=min_coverage)*(pis<=1e-03)]

# get the time info for the snp_samples
time_pair_idxs_snps_only, visno_snps_only, day_snps_only = parse_midas_data.calculate_time_pairs(subject_sample_time_map, snp_samples)

### different patient idx: 
# to compare results to time_pair idxs, we want different patient pair idxs. This helps us to contextualize if we are seeing events within patients that resemble replacements or modifications. 

# Calculate which pairs of idxs belong to the same sample, which to the same subject
# and which to different subjects
snp_same_sample_idxs, snp_same_subject_idxs, snp_diff_subject_idxs = parse_midas_data.calculate_subject_pairs(subject_sample_map, snp_samples)

############################
# Fixation matrix for SNPs #
############################

# Analyze SNPs, looping over chunk sizes. 
# Clunky, but necessary to limit memory usage on cluster

# Load SNP information for species_name
sys.stderr.write("Loading %s...\n" % species_name)

snp_difference_matrix = numpy.array([])
snp_opportunity_matrix = numpy.array([])
final_line_number = 0

while final_line_number >= 0:
    
    sys.stderr.write("Loading chunk starting @ %d...\n" % final_line_number)
    dummy_samples, allele_counts_map, passed_sites_map, final_line_number = parse_midas_data.parse_snps(species_name, debug=debug, allowed_samples=snp_samples,chunk_size=chunk_size,initial_line_number=final_line_number)
    sys.stderr.write("Done! Loaded %d genes\n" % len(allele_counts_map.keys()))
    
    # Calculate fixation matrix
    sys.stderr.write("Calculating matrix of snp differences...\n")
    chunk_snp_difference_matrix, chunk_snp_opportunity_matrix =     diversity_utils.calculate_fixation_matrix(allele_counts_map, passed_sites_map, min_change=min_change)    
    sys.stderr.write("Done!\n")
    
    if snp_difference_matrix.shape[0]==0:
        snp_difference_matrix = numpy.zeros_like(chunk_snp_difference_matrix)*1.0
        snp_opportunity_matrix = numpy.zeros_like(snp_difference_matrix)*1.0
    
    snp_difference_matrix += chunk_snp_difference_matrix
    snp_opportunity_matrix += chunk_snp_opportunity_matrix

sys.stderr.write("Done!\n")   

    
##################################################################    
# Compute gene gains vs losses for different pairs of individuals
##################################################################

# Load gene coverage information for species_name
sys.stderr.write("Loading pangenome data for %s...\n" % species_name)
gene_samples, gene_names, gene_presence_matrix, gene_depth_matrix, marker_coverages, gene_reads_matrix = parse_midas_data.parse_pangenome_data(species_name,allowed_samples=snp_samples)
sys.stderr.write("Done!\n")

gene_hamming_matrix_gain, gene_hamming_matrix_loss, num_opportunities = gene_diversity_utils.calculate_coverage_based_gene_hamming_matrix_gain_loss(gene_depth_matrix, marker_coverages, min_log2_fold_change=4)

######################################################################################################
# Compute indexes genes and SNPs so that everyone has low piS, high coverage for both SNPs and Genes #
###################################################################################################### 

# get the intersection of gene_samples and snp_samples satisfying piS being low and coverage being high both both. 
desired_samples = numpy.array(list(set(snp_samples) & set(gene_samples[marker_coverages>min_coverage])))   
# figure out what the indexes are for the desired samples in the snp_samples and gene_samples 
snp_sample_idx_map = parse_midas_data.calculate_sample_idx_map(desired_samples, snp_samples)
gene_sample_idx_map = parse_midas_data.calculate_sample_idx_map(desired_samples, gene_samples)

### time pair idxs where patients may have up to 3 time points:

# Calculate which pairs of idxs belong to different time points for both high cov and low piS
# time_pair_idxs is comprised of 2 arrays. The first array corresponds to indecies for the 1st visno. The second array corresponds to indecies for the second or 3d visno
time_pair_idxs, visno_snps_genes, day_snps_genes = parse_midas_data.calculate_time_pairs(subject_sample_time_map, desired_samples)

# since the time_pair_idx are in terms of desired samples ordering,I need to convert the desired_samples idxs to the snp_sample and gene_sample orders. 
time_pair_snp_idxs=parse_midas_data.apply_sample_index_map_to_indices(snp_sample_idx_map,  time_pair_idxs) #use these idxs to get the relevant fields from total_fixation_matrix
time_pair_gene_idxs=parse_midas_data.apply_sample_index_map_to_indices(gene_sample_idx_map,  time_pair_idxs) # use these idxs to get the relevant fields from gene_hamming_matrix


#### time pair idxs where patients can have exactly 1 time point (so that points plotted are iid)

time_pair_idxs_unique, visno_snps_genes_unique, day_snps_genes_unique = parse_midas_data.calculate_unique_time_pairs(subject_sample_time_map, desired_samples)
time_pair_snp_idxs_unique=parse_midas_data.apply_sample_index_map_to_indices(snp_sample_idx_map,  time_pair_idxs_unique) 
time_pair_gene_idxs_unique=parse_midas_data.apply_sample_index_map_to_indices(gene_sample_idx_map,  time_pair_idxs_unique) 


### different patient idx: 
# to compare results to time_pair idxs, we want different patient pair idxs. This helps us to contextualize if we are seeing events within patients that resemble replacements or modifications. 

# Calculate which pairs of idxs belong to the same sample, which to the same subject
# and which to different subjects
desired_same_sample_idxs, desired_same_subject_idxs, desired_diff_subject_idxs = parse_midas_data.calculate_subject_pairs(subject_sample_map, desired_samples)

# convert desired_idxs to the snp_sample_idx and gene_sample_idx values
# Ignore desired_same_subject_idxs because this is captured by the time_pair_idxs above. Not sure yet what to do with desired_same_sample_idxs.

diff_subject_snp_idxs = parse_midas_data.apply_sample_index_map_to_indices(snp_sample_idx_map, desired_diff_subject_idxs)  
diff_subject_gene_idxs = parse_midas_data.apply_sample_index_map_to_indices(gene_sample_idx_map, desired_diff_subject_idxs)  
 


##########
# Plot:  #
##########

# plot fixations vs days
# note that this is for patients that only satisfy coverage and piS requirements for SNPs -- gene constraints are not included here. 

# compute fraction of snp differences and clip for log scale
fraction_snp_difference=snp_difference_matrix/snp_opportunity_matrix
fraction_snp_difference=numpy.clip(fraction_snp_difference,1e-13,1)
diff_subject_control_time_pt=numpy.repeat(-10,len(snp_diff_subject_idxs[0]))

pylab.figure()
pylab.xlabel('Days')
pylab.ylabel('Fraction "fixations"')
pylab.title(species_name)
pylab.ylim(1e-13,1)  
pylab.semilogy(day_snps_only, fraction_snp_difference[time_pair_idxs_snps_only], 'go')
pylab.semilogy(diff_subject_control_time_pt, fraction_snp_difference[snp_diff_subject_idxs], 'ro')

pylab.legend(['Between time points','Different subjects'],'upper right',prop={'size':6})
    
pylab.savefig('%s/%s_fixation_vs_days_%0.1f.png' % (parse_midas_data.analysis_directory, species_name, min_change),bbox_inches='tight', dpi=300)



################################################### 
# plot gene gains and losses vs fixations of SNPs #
###################################################

# clip for log scale
gene_hamming_matrix_gain = numpy.clip(gene_hamming_matrix_gain,1e-13,1e09)
gene_hamming_matrix_loss = numpy.clip(gene_hamming_matrix_loss,1e-13,1e09)

pylab.figure() 
pylab.xlabel('Num substitutions')
pylab.ylabel('Num gene differences')
pylab.ylim([1e-14,1e04])
pylab.xlim([1e-14,1e05])
pylab.title(species_name)

pylab.loglog(fraction_snp_difference[diff_subject_snp_idxs], gene_hamming_matrix_loss[diff_subject_gene_idxs] + gene_hamming_matrix_gain[diff_subject_gene_idxs],'ro')
pylab.loglog(fraction_snp_difference[time_pair_snp_idxs], gene_hamming_matrix_gain[time_pair_gene_idxs],'yo')
pylab.loglog(fraction_snp_difference[time_pair_snp_idxs], gene_hamming_matrix_loss[time_pair_gene_idxs],'bo')

pylab.legend(['diff subjects differences', 'gains','losses'],'upper right',prop={'size':6})

pylab.savefig('%s/%s_gene_gain_loss_vs_substitutions.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight',dpi=300)

### redo plot with unique time pairs (so that every point is iid)
pylab.figure() 
pylab.xlabel('Num substitutions')
pylab.ylabel('Num gene differences')
pylab.ylim([1e-14,1e04])
pylab.xlim([1e-14,1e05])
pylab.title(species_name)

pylab.loglog(fraction_snp_difference[diff_subject_snp_idxs], gene_hamming_matrix_loss[diff_subject_gene_idxs] + gene_hamming_matrix_gain[diff_subject_gene_idxs],'ro')
pylab.loglog(fraction_snp_difference[time_pair_snp_idxs_unique], gene_hamming_matrix_gain[time_pair_gene_idxs_unique],'yo')
pylab.loglog(fraction_snp_difference[time_pair_snp_idxs_unique], gene_hamming_matrix_loss[time_pair_gene_idxs_unique],'bo')

pylab.legend(['diff subjects differences', 'gains','losses'],'upper right',prop={'size':6})

pylab.savefig('%s/%s_gene_gain_loss_vs_substitutions_unique.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight',dpi=300)



#######################################
# Plot gene gains and losses vs days  #
#######################################

pylab.figure() 
pylab.xlabel('Days')
pylab.ylabel('Num gene differences')
pylab.ylim([1e-14,1e04])
pylab.title(species_name)

diff_subject_control_time_pt=numpy.repeat(-10,len(diff_subject_snp_idxs[0]))

pylab.semilogy(diff_subject_control_time_pt, gene_hamming_matrix_loss[diff_subject_gene_idxs] + gene_hamming_matrix_gain[diff_subject_gene_idxs],'ro')
pylab.semilogy(day_snps_genes, gene_hamming_matrix_gain[time_pair_gene_idxs],'yo')
jitter_day_snps_genes_unique = list(numpy.asarray(day_snps_genes) + 3)
pylab.semilogy(jitter_day_snps_genes, gene_hamming_matrix_loss[time_pair_gene_idxs],'bo')

pylab.legend(['diff subjects differences', 'gains','losses'],'upper right',prop={'size':6})

pylab.savefig('%s/%s_gene_gain_loss_vs_days.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight',dpi=300)


### redo plot with unique time pairs (so that every point is iid)

pylab.figure() 
pylab.xlabel('Days')
pylab.ylabel('Num gene differences')
pylab.ylim([1e-14,1e04])
pylab.title(species_name)

diff_subject_control_time_pt=numpy.repeat(-10,len(diff_subject_snp_idxs[0]))

pylab.semilogy(diff_subject_control_time_pt, gene_hamming_matrix_loss[diff_subject_gene_idxs] + gene_hamming_matrix_gain[diff_subject_gene_idxs],'ro')
pylab.semilogy(day_snps_genes_unique, gene_hamming_matrix_gain[time_pair_gene_idxs_unique],'yo')
jitter_day_snps_genes_unique = list(numpy.asarray(day_snps_genes_unique) + 3)
pylab.semilogy(jitter_day_snps_genes_unique, gene_hamming_matrix_loss[time_pair_gene_idxs_unique],'bo')

pylab.legend(['diff subjects differences', 'gains','losses'],'upper right',prop={'size':6})

pylab.savefig('%s/%s_gene_gain_loss_vs_days_unique.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight',dpi=300)



##############################################################
# Plot gene gain vs loss: are there more gains than losses ? #
#############################################################

pylab.figure() 
pylab.xlabel('Number of gene gains')
pylab.ylabel('Number of gene losses')
pylab.xlim([1e-14,1e05])
pylab.ylim([1e-14,1e05])
pylab.title(species_name)

pylab.loglog(gene_hamming_matrix_gain[time_pair_gene_idxs_unique], gene_hamming_matrix_loss[time_pair_gene_idxs_unique], 'go')
diff_subject_control=(gene_hamming_matrix_loss[diff_subject_gene_idxs] + gene_hamming_matrix_gain[diff_subject_gene_idxs])/2
pylab.loglog(diff_subject_control, diff_subject_control,'ro')

pylab.plot([1e-14,1e05], [1e-14,1e05], ls="--", c=".3")

pylab.legend(['diff subjects', 'within subject'],'upper right',prop={'size':6})

pylab.savefig('%s/%s_gene_gain_vs_loss_unique.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight',dpi=300)
