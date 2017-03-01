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


# load list of metaphlan2 genes
metaphlan2_genes = parse_midas_data.load_metaphlan2_genes(species_name)

# Load subject and sample metadata
sys.stderr.write("Loading HMP metadata...\n")
subject_sample_map = parse_midas_data.parse_subject_sample_map()
sys.stderr.write("Done!\n")
    
# Load genomic coverage distributions
sample_coverage_histograms, samples = parse_midas_data.parse_coverage_distribution(species_name, combination_type="sample")
median_coverages = numpy.array([stats_utils.calculate_median_from_histogram(sample_coverage_histogram) for sample_coverage_histogram in sample_coverage_histograms])
    
sample_coverage_map = {samples[i]: median_coverages[i] for i in xrange(0,len(samples))}
    
# Load SNP information for species_name
sys.stderr.write("Loading %s...\n" % species_name)
samples, allele_counts_map, passed_sites_map = parse_midas_data.parse_snps(species_name, combination_type="sample", debug=False)
sys.stderr.write("Done!\n")
    
median_coverages = numpy.array([sample_coverage_map[samples[i]] for i in xrange(0,len(samples))])
    
# Calculate full matrix of synonymous pairwise differences
sys.stderr.write("Calculate synonymous pi matrix...\n")
pi_matrix_syn, avg_pi_matrix_syn = diversity_utils.calculate_pi_matrix(allele_counts_map, passed_sites_map, variant_type='4D', allowed_genes=metaphlan2_genes)

# Calculate fixation matrix
fixation_matrix_syn = diversity_utils.calculate_fixation_matrix(allele_counts_map, passed_sites_map, variant_type='4D', min_change=min_change, allowed_genes=metaphlan2_genes)
    
sys.stderr.write("Done!\n")
    
# Calculate full matrix of nonsynonymous pairwise differences
sys.stderr.write("Calculate nonsynonymous pi matrix...\n")
# Calculate allele count matrices
pi_matrix_non, avg_pi_matrix_non = diversity_utils.calculate_pi_matrix(allele_counts_map, passed_sites_map, variant_type='1D',allowed_genes=metaphlan2_genes)
# Calculate fixation matrix
fixation_matrix_non = diversity_utils.calculate_fixation_matrix(allele_counts_map, passed_sites_map, variant_type='1D', min_change=min_change,allowed_genes=metaphlan2_genes)
    
    
sys.stderr.write("Done!\n")
    
# Calculate which pairs of idxs belong to the same sample, which to the same subject
# and which to different subjects
same_sample_idxs, same_subject_idxs, diff_subject_idxs = parse_midas_data.calculate_subject_pairs(subject_sample_map, samples)
    
    
fst_matrix_syn = numpy.clip(pi_matrix_syn-avg_pi_matrix_syn,0,1)
pis = numpy.diag(pi_matrix_syn)
desired_idxs = numpy.nonzero(pis<1e-03)[0]
    
print diversity_utils.phylip_distance_matrix_str(fst_matrix_syn[numpy.ix_(desired_idxs, desired_idxs)], [samples[idx] for idx in desired_idxs])
    
fst_matrix_syn = numpy.clip(fst_matrix_syn,1e-06,1)
    
total_fixation_matrix = fixation_matrix_syn + fixation_matrix_non
total_fixation_matrix = numpy.clip(total_fixation_matrix, 1e-01, 1e09)
    
# Done calculating... now plot figure!
pylab.figure()
pylab.xlabel('$\\pi_s$')
pylab.ylabel('$\\pi_n/\\pi_s$')
pylab.ylim([0,1.1])
pylab.title(species_name)

pylab.semilogx(pi_matrix_syn[diff_subject_idxs], (pi_matrix_non/pi_matrix_syn)[diff_subject_idxs],'r.')

pylab.semilogx(pi_matrix_syn[same_subject_idxs], (pi_matrix_non/pi_matrix_syn)[same_subject_idxs],'g.')

pylab.semilogx(pi_matrix_syn[same_sample_idxs], (pi_matrix_non/pi_matrix_syn)[same_sample_idxs],'b.')

pylab.savefig('%s/%s_pNpS_vs_pi_metaphlan2_genes.pdf' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight')
    
pylab.figure()
pylab.xlabel('$F_s$')
pylab.ylabel('$\\pi_n/\\pi_s$')
pylab.ylim([0,1.1])
pylab.title(species_name)
  
pylab.semilogx(fst_matrix_syn[diff_subject_idxs], (pi_matrix_non/pi_matrix_syn)[diff_subject_idxs],'r.')

pylab.semilogx(fst_matrix_syn[same_subject_idxs], (pi_matrix_non/pi_matrix_syn)[same_subject_idxs],'g.')

pylab.semilogx(fst_matrix_syn[same_sample_idxs], (pi_matrix_non/pi_matrix_syn)[same_sample_idxs],'b.')
    
#pylab.savefig('%s_pNpS_vs_fst.pdf' % species_name,bbox_inches='tight')
    
    
pylab.figure()
pylab.xlabel('$\\pi_s$')
pylab.ylabel('$F_s$')
pylab.title(species_name)
  
pylab.loglog(avg_pi_matrix_syn[diff_subject_idxs], fst_matrix_syn[diff_subject_idxs],'r.', label='Different subjects')
    
pylab.loglog(avg_pi_matrix_syn[same_subject_idxs], fst_matrix_syn[same_subject_idxs], 'g.', label='Same subject, different timepoints')
    
#pylab.loglog(avg_pi_matrix_syn[same_sample_idxs], fst_matrix_syn[same_sample_idxs], 'b.')
pylab.loglog([1e-06,1e-01],[1e-06,1e-01],'k:')
pylab.xlim([1e-06,1e-01])
pylab.ylim([1e-06,1e-01])
pylab.savefig('%s/%s_fst_vs_pi_metaphlan2_genes.pdf' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight')
    
    
    
pylab.figure()
pylab.xlabel('Within-sample $\\pi_s$ (avg)')
pylab.ylabel('Between-sample $\\pi_s$')
pylab.title(species_name)
  
pylab.loglog(avg_pi_matrix_syn[diff_subject_idxs], pi_matrix_syn[diff_subject_idxs],'r.', label='Different subjects')
    
pylab.loglog(avg_pi_matrix_syn[same_subject_idxs], pi_matrix_syn[same_subject_idxs], 'g.', label='Same subject, different timepoints')
    
#pylab.loglog(avg_pi_matrix_syn[same_sample_idxs], pi_matrix_syn[same_sample_idxs], 'b.')
pylab.loglog([1e-06,1e-01],[1e-06,1e-01],'k:')
pylab.xlim([1e-06,1e-01])
pylab.ylim([1e-06,1e-01])
pylab.legend(loc='lower right',frameon=False)
pylab.savefig('%s/%s_pi_between_vs_within_metaphlan2_genes.pdf' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight')
pylab.savefig('%s/%s_pi_between_vs_within_metaphlan2_genes.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight', dpi=300)
    
pylab.figure()
pylab.xlabel('Within-sample $\\pi_s$ (avg)')
pylab.ylabel('Num "fixations"')
pylab.title(species_name)
  
pylab.loglog(avg_pi_matrix_syn[diff_subject_idxs], total_fixation_matrix[diff_subject_idxs],'r.', label='Different subjects')
    
    
pylab.loglog(avg_pi_matrix_syn[same_subject_idxs], total_fixation_matrix[same_subject_idxs], 'g.', label='Same subject, different timepoints')
    
#pylab.loglog(avg_pi_matrix_syn[same_sample_idxs], pi_matrix_syn[same_sample_idxs], 'b.')
pylab.xlim([1e-06,1e-01])
pylab.ylim([3e-02,1e04])
pylab.semilogx([1e-06,1e-01],[3e-01,3e-01],'k:')
    
pylab.savefig('%s/%s_fixation_vs_pi_%0.1f_metaphlan2_genes.pdf' % (parse_midas_data.analysis_directory, species_name, min_change),bbox_inches='tight')
pylab.savefig('%s/%s_fixation_vs_pi_%0.1f_metaphlan2_genes.png' % (parse_midas_data.analysis_directory, species_name, min_change),bbox_inches='tight', dpi=300)

pylab.figure()
pylab.xlabel('Median coverage')
pylab.ylabel('Within-sample $\\pi_s$')
print median_coverages.shape
print numpy.diag(pi_matrix_syn).shape
pylab.loglog(median_coverages, numpy.diag(pi_matrix_syn),'k.')
pylab.savefig('%s/%s_pi_vs_coverage_metaphlan2_genes.pdf' % (parse_midas_data.analysis_directory, species_name), bbox_inches='tight')
    

    
