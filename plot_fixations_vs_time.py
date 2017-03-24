import matplotlib  
matplotlib.use('Agg') 
import parse_midas_data
import pylab
import sys
import numpy
import diversity_utils
import stats_utils
import os

min_change = 0.8


################################################################################
#
# Standard header to read in argument information
#
################################################################################

if len(sys.argv)>1:
    if len(sys.argv)>2:
        debug=True # debug does nothing in this script
        species_name=sys.argv[2]
    else:
        debug=False
        species_name=sys.argv[1]
else:
    sys.stderr.write("Usage: python plot_fixations_vs_time.py [debug] species_name")
#############################################################################

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
samples, allele_counts_map, passed_sites_map = parse_midas_data.parse_snps(species_name, debug)
sys.stderr.write("Done!\n")
    
median_coverages = numpy.array([sample_coverage_map[samples[i]] for i in xrange(0,len(samples))])
    
# Calculate fixation matrix for synonymous variants
fixation_matrix_syn, persite_fixation_matrix_syn = diversity_utils.calculate_fixation_matrix(allele_counts_map, passed_sites_map, variant_type='4D', min_change=min_change)
    
    
# Calculate fixation matrix for non synonynous variants
fixation_matrix_non, persite_fixation_matrix_non = diversity_utils.calculate_fixation_matrix(allele_counts_map, passed_sites_map, variant_type='1D', min_change=min_change)
    

fixation_matrix_non = numpy.clip(fixation_matrix_non, 1e-01, 1e09)    
fixation_matrix_syn = numpy.clip(fixation_matrix_syn, 1e-01, 1e09)    
total_fixation_matrix = fixation_matrix_syn + fixation_matrix_non
total_fixation_matrix = numpy.clip(total_fixation_matrix, 1e-01, 1e09)

# Calculate which pairs of idxs belong to the same sample, which to the same subject
# and which to different subjects
same_sample_idxs, same_subject_idxs, diff_subject_idxs = parse_midas_data.calculate_subject_pairs(subject_sample_map, samples)


time_pair_idxs, visno, day = parse_midas_data.calculate_time_pairs(subject_sample_time_map, samples)
    
    


# write to an intermediate file so that I can plot all species' time series plots on one plot (plot_fixations_vs_time_multispecies.py)

numpy.savez(os.path.expanduser('~/tmp_intermediate_files/time_fix_%s_%s.npz' %(species_name, min_change)), day, total_fixation_matrix[time_pair_idxs], fixation_matrix_non[time_pair_idxs], fixation_matrix_syn[time_pair_idxs], species_name, min_change)


#plot fixations vs days

pylab.figure()
pylab.xlabel('days')
pylab.ylabel('Num "fixations"')
pylab.title(species_name)
  
pylab.semilogy(day, total_fixation_matrix[time_pair_idxs], 'ro')
    
pylab.savefig('%s/%s_fixation_vs_days_%0.1f.pdf' % (parse_midas_data.analysis_directory, species_name, min_change),bbox_inches='tight')
pylab.savefig('%s/%s_fixation_vs_days_%0.1f.png' % (parse_midas_data.analysis_directory, species_name, min_change),bbox_inches='tight', dpi=300)


# plot syn vs nonsyn fixations vs days:
pylab.figure()
pylab.xlabel('days')
pylab.ylabel('Num "fixations"')
pylab.title(species_name)
  
pylab.semilogy(day, fixation_matrix_non[time_pair_idxs], 'ro')
pylab.semilogy(day, fixation_matrix_syn[time_pair_idxs], 'bo')
    
pylab.savefig('%s/%s_dN_and_dS_vs_days_%0.1f.pdf' % (parse_midas_data.analysis_directory, species_name, min_change),bbox_inches='tight')
pylab.savefig('%s/%s_dN_and_dS_vs_days_%0.1f.png' % (parse_midas_data.analysis_directory, species_name, min_change),bbox_inches='tight', dpi=300)


# plot dN/dS vs days:
pylab.figure()
pylab.xlabel('days')
pylab.ylabel('Dn/Ds')
pylab.title(species_name)
  
pylab.semilogy(day, (fixation_matrix_non/fixation_matrix_syn)[time_pair_idxs], 'ro')
    
pylab.savefig('%s/%s_dNdS_vs_days_%0.1f.pdf' % (parse_midas_data.analysis_directory, species_name, min_change),bbox_inches='tight')
pylab.savefig('%s/%s_dNdS_vs_days_%0.1f.png' % (parse_midas_data.analysis_directory, species_name, min_change),bbox_inches='tight', dpi=300)



#plot fixations vs visnos

pylab.figure()
pylab.xlabel('Visit Number')
pylab.ylabel('Num "fixations"')
pylab.title(species_name)
  
pylab.semilogy(visno, total_fixation_matrix[time_pair_idxs], 'ro')
    
pylab.savefig('%s/%s_fixation_vs_visnos_%0.1f.pdf' % (parse_midas_data.analysis_directory, species_name, min_change),bbox_inches='tight')
pylab.savefig('%s/%s_fixation_vs_visnos_%0.1f.png' % (parse_midas_data.analysis_directory, species_name, min_change),bbox_inches='tight', dpi=300)


