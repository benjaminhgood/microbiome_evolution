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
clipped_pis = (total_pis+1)/(total_pi_opportunities+1)


median_coverages = numpy.array([sample_coverage_map[samples[i]] for i in xrange(0,len(samples))])

###############################################################
# Indexes for SNP samples that have low piS and high coverage #
###############################################################

# Only plot samples above a certain depth threshold that are "haploids"
high_cov_samples = samples[(median_coverages>=min_coverage)]
high_cov_pis     = clipped_pis[(median_coverages>=min_coverage)]

# get the time info for the snp_samples
time_pair_idxs, visno, day = parse_midas_data.calculate_time_pairs(subject_sample_time_map, high_cov_samples)


#### time pair idxs where patients can have exactly 1 time point (so that points plotted are iid)
time_pair_idxs_unique, visno_snps_genes_unique, day_snps_genes_unique = parse_midas_data.calculate_unique_time_pairs(subject_sample_time_map, high_cov_samples)


### different patient idx: 
# to compare results to time_pair idxs, we want different patient pair idxs. This helps us to contextualize if we are seeing events within patients that resemble replacements or modifications. 

# Calculate which pairs of idxs belong to the same sample, which to the same subject
# and which to different subjects
snp_same_sample_idxs, snp_same_subject_idxs, snp_diff_subject_idxs = parse_midas_data.calculate_subject_pairs(subject_sample_map, high_cov_samples)





##########
# Plot:  #
##########


pylab.figure() 
pylab.xlabel('First time point')
pylab.ylabel('Subsequent time point')
pylab.xlim([1e-5,1e-1])
pylab.ylim([1e-5,1e-1])
pylab.title(species_name)

pylab.loglog(high_cov_pis[time_pair_idxs[0]],high_cov_pis[time_pair_idxs[1]], 'go',alpha=0.5)
pylab.plot([1e-5,1e-1], [1e-5,1e-1], ls="--", c=".3")


pylab.savefig('%s/%s_pis_ordered_by_time.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight',dpi=300)
