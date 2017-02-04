import matplotlib  
matplotlib.use('Agg') 
import parse_midas_data
import pylab
import sys
import numpy
import diversity_utils
import stats_utils
species_name=sys.argv[1]

# Load subject and sample metadata
sys.stderr.write("Loading HMP metadata...\n")
subject_sample_map = parse_midas_data.parse_subject_sample_map()
sys.stderr.write("Done!\n")
    
# Load SNP information for species_name
sys.stderr.write("Loading %s...\n" % species_name)
samples, allele_counts_map, passed_sites_map = parse_midas_data.parse_snps(species_name, combination_type="sample", debug=False)
sys.stderr.write("Done!\n")
    
# Calculate full matrix of synonymous pairwise differences
sys.stderr.write("Calculate synonymous pi matrix...\n")
pi_matrix_syn, avg_pi_matrix_syn = diversity_utils.calculate_pi_matrix(allele_counts_map, passed_sites_map, variant_type='4D')
    
# Calculate which pairs of idxs belong to the same sample, which to the same subject
# and which to different subjects
same_sample_idxs, same_subject_idxs, diff_subject_idxs = parse_midas_data.calculate_subject_pairs(subject_sample_map, samples)
    
pi_within = pi_matrix_syn[same_sample_idxs]
pi_within.sort()
    
pi_timepoints = pi_matrix_syn[same_subject_idxs]
pi_timepoints.sort()
    
pi_between = pi_matrix_syn[diff_subject_idxs]
pi_between.sort()
    
pylab.figure(figsize=(5,3))
pylab.xlabel('Synonymous diversity, $\\pi_s$')
pylab.ylabel('Survival function')
pylab.gca().spines['top'].set_visible(False)
pylab.gca().spines['right'].set_visible(False)
pylab.gca().get_xaxis().tick_bottom()
pylab.gca().get_yaxis().tick_left()

pylab.step(pi_between, 1-numpy.arange(0,len(pi_between))*1.0/len(pi_between),color='r',where='post',label='Between people')
pylab.step(pi_within, 1-numpy.arange(0,len(pi_within))*1.0/len(pi_within),color='b',where='post',label='Within people')
pylab.step(pi_timepoints, 1-numpy.arange(0,len(pi_timepoints))*1.0/len(pi_timepoints),color='g',where='post',label='Across time')
pylab.legend(loc='lower left',frameon=False,fontsize=9)
pylab.semilogx([1e-02],[-1])
pylab.ylim([0,1])
pylab.xlim([1e-06,1e-01])
pylab.savefig('%s/%s_pi_distribution.pdf' % (parse_midas_data.analysis_directory,species_name), bbox_inches='tight')
#pylab.show()
    
    
