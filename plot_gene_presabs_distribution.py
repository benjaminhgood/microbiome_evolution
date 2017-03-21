import matplotlib  
matplotlib.use('Agg') 
import parse_midas_data
import pylab
import sys
import numpy
import diversity_utils
import stats_utils

########################################################################################
#
# Standard header to read in argument information
#
########################################################################################
if len(sys.argv)>1:
    if len(sys.argv)>2:
        debug=True # debug does nothing in this script
        species_name=sys.argv[2]
    else:
        debug=False
        species_name=sys.argv[1]
else:
    sys.stderr.write("Usage: python plot_pNpS_vs_pi.py [debug] species_name")
########################################################################################


# Load subject and sample metadata
sys.stderr.write("Loading HMP metadata...\n")
subject_sample_map = parse_midas_data.parse_subject_sample_map()
sys.stderr.write("Done!\n")
    
# Load gene presence/absence information for species_name
sys.stderr.write("Loading %s...\n" % species_name)
samples, gene_names, gene_presence_matrix = parse_midas_data.parse_gene_presences(species_name)
sys.stderr.write("Done!\n")
    
# Calculate matrix of number of genes that differ
sys.stderr.write("Calculate gene hamming matrix...\n")
gene_hamming_matrix = diversity_utils.calculate_gene_hamming_matrix(gene_presence_matrix)
    

# Calculate fraction of shared genes
sys.stderr.write("Calculate gene sharing matrix...\n")
gene_sharing_matrix = diversity_utils.calculate_gene_sharing_matrix(gene_presence_matrix)
    
    
# Calculate which pairs of idxs belong to the same sample, which to the same subject
# and which to different subjects
same_sample_idxs, same_subject_idxs, diff_subject_idxs = parse_midas_data.calculate_subject_pairs(subject_sample_map, samples)
    
   
hamming_timepoints = gene_hamming_matrix[same_subject_idxs]
hamming_timepoints.sort()
hamming_timepoints_dns, hamming_timepoints_survivals = stats_utils.calculate_unnormalized_survival_from_vector(hamming_timepoints, min_x=0.1, max_x=1e05)
hamming_timepoints_survivals /= hamming_timepoints_survivals[0]    
    
hamming_between = gene_hamming_matrix[diff_subject_idxs]
hamming_between.sort()
hamming_between_dns, hamming_between_survivals = stats_utils.calculate_unnormalized_survival_from_vector(hamming_between, min_x=0.1, max_x=1e05)
hamming_between_survivals /= hamming_between_survivals[0]

gene_counts = gene_presence_matrix.sum(axis=0)
gene_counts.sort()
gene_count_ns, gene_count_survivals = stats_utils.calculate_unnormalized_survival_from_vector(gene_counts, min_x=0.1, max_x=1e05)
gene_count_survivals /= gene_count_survivals[0]


print "Median gene count=", numpy.median(gene_presence_matrix.sum(axis=0))    
print "Median timepoints=", numpy.median(hamming_timepoints)
print "Median between=", numpy.median(hamming_between)

pylab.figure(1,figsize=(5,3))
pylab.title(species_name,fontsize=11)
pylab.xlabel('Number of genes',fontsize=11)
pylab.ylabel('Survival function',fontsize=11)
pylab.gca().spines['top'].set_visible(False)
pylab.gca().spines['right'].set_visible(False)
pylab.gca().get_xaxis().tick_bottom()
pylab.gca().get_yaxis().tick_left()


pylab.step(hamming_between_dns, hamming_between_survivals,color='r',label='Differ between people')
pylab.step(hamming_timepoints_dns, hamming_timepoints_survivals,color='g',label='Differ across time')
pylab.step(gene_count_ns, gene_count_survivals,color='b',label='Present per person')

pylab.legend(loc='lower left',frameon=False,fontsize=9)
pylab.semilogx([1e-02],[-1])
pylab.ylim([0,1.05])
pylab.xlim([1,1e04])
pylab.savefig('%s/%s_gene_hamming_distribution.pdf' % (parse_midas_data.analysis_directory,species_name), bbox_inches='tight')
#pylab.show()

pylab.figure(2,figsize=(5,3))
pylab.title(species_name,fontsize=11)
pylab.xlabel('Fraction shared genes',fontsize=11)
pylab.ylabel('Survival function',fontsize=11)
pylab.gca().spines['top'].set_visible(False)
pylab.gca().spines['right'].set_visible(False)
pylab.gca().get_xaxis().tick_bottom()
pylab.gca().get_yaxis().tick_left()

sharing_timepoints = gene_sharing_matrix[same_subject_idxs]
sharing_timepoints.sort()
sharing_timepoints_dns, sharing_timepoints_survivals = stats_utils.calculate_unnormalized_survival_from_vector(sharing_timepoints, min_x=0.1, max_x=1e05)
sharing_timepoints_survivals /= sharing_timepoints_survivals[0]    
    
sharing_between = gene_sharing_matrix[diff_subject_idxs]
sharing_between.sort()
sharing_between_dns, sharing_between_survivals = stats_utils.calculate_unnormalized_survival_from_vector(sharing_between, min_x=0.1, max_x=1e05)
sharing_between_survivals /= sharing_between_survivals[0]

pylab.step(sharing_between_dns, sharing_between_survivals,color='r',label='Between people')
pylab.step(sharing_timepoints_dns, sharing_timepoints_survivals,color='g',label='Across time')

pylab.legend(loc='lower left',frameon=False,fontsize=9)
#pylab.semilogx([1e-02],[-1])
pylab.ylim([0,1.05])
pylab.xlim([0,1.05])
pylab.savefig('%s/%s_gene_sharing_distribution.pdf' % (parse_midas_data.analysis_directory,species_name), bbox_inches='tight')
#pylab.show()
    
    
