import matplotlib  
matplotlib.use('Agg') 
import parse_midas_data
import pylab
import sys
import numpy
import diversity_utils
import stats_utils
species_name=sys.argv[1]

min_change = 0.8


# Load subject and sample metadata
sys.stderr.write("Loading HMP metadata...\n")
subject_sample_map = parse_midas_data.parse_subject_sample_map()
sys.stderr.write("Done!\n")
    
# Load time data
subject_sample_time_map = parse_midas_data.parse_subject_sample_time_map()


# Load genomic coverage distributions
sample_coverage_histograms, samples = parse_midas_data.parse_coverage_distribution(species_name, combination_type="sample")
median_coverages = numpy.array([stats_utils.calculate_median_from_histogram(sample_coverage_histogram) for sample_coverage_histogram in sample_coverage_histograms])
    
sample_coverage_map = {samples[i]: median_coverages[i] for i in xrange(0,len(samples))}
    
# Load SNP information for species_name
sys.stderr.write("Loading %s...\n" % species_name)
samples, allele_counts_map, passed_sites_map = parse_midas_data.parse_snps(species_name, combination_type="sample", debug=False)
sys.stderr.write("Done!\n")
    
median_coverages = numpy.array([sample_coverage_map[samples[i]] for i in xrange(0,len(samples))])
    
# Calculate fixation matrix for synonymous variants
fixation_matrix_syn = diversity_utils.calculate_fixation_matrix(allele_counts_map, passed_sites_map, variant_type='4D', min_change=min_change)
    
    
# Calculate fixation matrix for non synonynous variants
fixation_matrix_non = diversity_utils.calculate_fixation_matrix(allele_counts_map, passed_sites_map, variant_type='1D', min_change=min_change)
    
    
total_fixation_matrix = fixation_matrix_syn + fixation_matrix_non
total_fixation_matrix = numpy.clip(total_fixation_matrix, 1e-01, 1e09)

# Calculate which pairs of idxs belong to the same sample, which to the same subject
# and which to different subjects
same_sample_idxs, same_subject_idxs, diff_subject_idxs = parse_midas_data.calculate_subject_pairs(subject_sample_map, samples)


time_pair_idxs, visno, day = parse_midas_data.calculate_time_pairs(subject_sample_time_map, samples)
    
    
total_fixation_matrix[time_pair_idxs]


#plot fixations vs days

pylab.figure()
pylab.xlabel('days')
pylab.ylabel('Num "fixations"')
pylab.title(species_name)
  
pylab.plot(day, total_fixation_matrix[time_pair_idxs], 'ro')
    
pylab.savefig('%s/%s_fixation_vs_days_%0.1f.pdf' % (parse_midas_data.analysis_directory, species_name, min_change),bbox_inches='tight')
pylab.savefig('%s/%s_fixation_vs_days_%0.1f.png' % (parse_midas_data.analysis_directory, species_name, min_change),bbox_inches='tight', dpi=300)



#plot fixations vs visnos

pylab.figure()
pylab.xlabel('Visit Number')
pylab.ylabel('Num "fixations"')
pylab.title(species_name)
  
pylab.plot(visno, total_fixation_matrix[time_pair_idxs], 'ro')
    
pylab.savefig('%s/%s_fixation_vs_visnos_%0.1f.pdf' % (parse_midas_data.analysis_directory, species_name, min_change),bbox_inches='tight')
pylab.savefig('%s/%s_fixation_vs_visnos_%0.1f.png' % (parse_midas_data.analysis_directory, species_name, min_change),bbox_inches='tight', dpi=300)

