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



# Load subject and sample metadata
sys.stderr.write("Loading HMP metadata...\n")
subject_sample_map = parse_midas_data.parse_subject_sample_map()
sys.stderr.write("Done!\n")
    

# Load time data
subject_sample_time_map = parse_midas_data.parse_subject_sample_time_map()

# Load genomic coverage distributions
sample_coverage_histograms, samples = parse_midas_data.parse_coverage_distribution(species_name)
median_coverages = numpy.array([stats_utils.calculate_median_from_histogram(sample_coverage_histogram) for sample_coverage_histogram in sample_coverage_histograms])
    
sample_coverage_map = {samples[i]: median_coverages[i] for i in xrange(0,len(samples))}
    
# Load SNP information for species_name
sys.stderr.write("Loading %s...\n" % species_name)

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




# Calculate which pairs of idxs belong to different time points for both high cov and low piS
time_pair_idxs, visno, day = parse_midas_data.calculate_time_pairs(subject_sample_time_map, snp_samples)

    
    
# write to an intermediate file so that I can plot all species' time series plots on one plot (plot_fixations_vs_time_multispecies.py)

#numpy.savez(os.path.expanduser('~/tmp_intermediate_files/time_fix_%s_%s.npz' %(species_name, min_change)), day=day, total_fixation_matrix=total_fixation_matrix, fixation_matrix_non=fixation_matrix_non, fixation_matrix_syn=fixation_matrix_syn, time_pair_idxs=time_pair_idxs, visno=visno, high_coverage_low_pi_time_pair_idxs=high_coverage_low_pi_time_pair_idxs,high_coverage_low_pi_visno=high_coverage_low_pi_visno, high_coverage_low_pi_day=high_coverage_low_pi_day, fixation_opportunities_non=fixation_opportunities_non, fixation_opportunities_syn=fixation_opportunities_syn)






##########
# Plot:  #
##########

#plot fixations vs days

pylab.figure()
pylab.xlabel('days')
pylab.ylabel('Num "fixations"')
pylab.title(species_name)
  
pylab.semilogy(day, (snp_difference_matrix/snp_opportunity_matrix)[time_pair_idxs], 'ro')
    
pylab.savefig('%s/%s_fixation_vs_days_%0.1f.pdf' % (parse_midas_data.analysis_directory, species_name, min_change),bbox_inches='tight')
pylab.savefig('%s/%s_fixation_vs_days_%0.1f.png' % (parse_midas_data.analysis_directory, species_name, min_change),bbox_inches='tight', dpi=300)


