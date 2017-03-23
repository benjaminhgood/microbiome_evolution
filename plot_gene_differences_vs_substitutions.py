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
    sys.stderr.write("Usage: python command.py [debug] species_name")
########################################################################################

min_change = 0.8
min_coverage = 20


# Load subject and sample metadata
sys.stderr.write("Loading HMP metadata...\n")
subject_sample_map = parse_midas_data.parse_subject_sample_map()
sys.stderr.write("Done!\n")
    
# Load genomic coverage distributions
sample_coverage_histograms, samples = parse_midas_data.parse_coverage_distribution(species_name)
median_coverages = numpy.array([stats_utils.calculate_nonzero_median_from_histogram(sample_coverage_histogram) for sample_coverage_histogram in sample_coverage_histograms])
sample_coverage_map = {samples[i]: median_coverages[i] for i in xrange(0,len(samples))}
    
  

# Load SNP information for species_name
sys.stderr.write("Loading %s...\n" % species_name)
snp_samples, allele_counts_map, passed_sites_map = parse_midas_data.parse_snps(species_name, debug)
sys.stderr.write("Done!\n")
 
median_coverages = numpy.array([sample_coverage_map[snp_samples[i]] for i in xrange(0,len(snp_samples))])
  
# Calculate full matrix of synonymous pairwise differences
sys.stderr.write("Calculate synonymous pi matrix...\n")
pi_matrix_syn, avg_pi_matrix_syn = diversity_utils.calculate_pi_matrix(allele_counts_map, passed_sites_map, variant_type='4D')

pi_matrix_syn = numpy.clip(pi_matrix_syn,1e-06,1)
avg_pi_matrix_syn = numpy.clip(avg_pi_matrix_syn,1e-06,1)
pis = numpy.diag(pi_matrix_syn)

# Calculate fixation matrix
fixation_matrix_syn, persite_fixation_matrix_syn = diversity_utils.calculate_fixation_matrix(allele_counts_map, passed_sites_map, variant_type='4D', min_change=min_change)    
sys.stderr.write("Done!\n")
    
# Calculate full matrix of nonsynonymous pairwise differences
sys.stderr.write("Calculate nonsynonymous pi matrix...\n")
# Calculate allele count matrices
pi_matrix_non, avg_pi_matrix_non = diversity_utils.calculate_pi_matrix(allele_counts_map, passed_sites_map, variant_type='1D')
# Calculate fixation matrix
fixation_matrix_non, persite_fixation_matrix_non = diversity_utils.calculate_fixation_matrix(allele_counts_map, passed_sites_map, variant_type='1D', min_change=min_change)
sys.stderr.write("Done!\n")

# Calculate total fixation matrix
total_fixation_matrix = fixation_matrix_syn + fixation_matrix_non

# Load gene presence/absence information for species_name
sys.stderr.write("Loading %s...\n" % species_name)
gene_samples, gene_names, gene_presence_matrix, gene_depth_matrix, marker_coverages, gene_reads_matrix = parse_midas_data.parse_pangenome_data(species_name)
sys.stderr.write("Done!\n")
    
total_num_genes = gene_presence_matrix.sum(axis=0)    
avg_num_genes = (total_num_genes[:,None]+total_num_genes[None,:])/2.0
    
# Calculate matrix of number of genes that differ
sys.stderr.write("Calculate gene hamming matrix...\n")
gene_hamming_matrix = diversity_utils.calculate_coverage_based_gene_hamming_matrix(gene_depth_matrix, marker_coverages, min_log2_fold_change=4)

# Now need to make the gene samples and snp samples match up
desired_samples = list(set(snp_samples[(median_coverages>=min_coverage)*(pis<1e-03)]) & set(gene_samples[marker_coverages>min_coverage]))   
     
# Calculate which pairs of idxs belong to the same sample, which to the same subject
# and which to different subjects
desired_same_sample_idxs, desired_same_subject_idxs, desired_diff_subject_idxs = parse_midas_data.calculate_subject_pairs(subject_sample_map, desired_samples)

snp_sample_idx_map = parse_midas_data.calculate_sample_idx_map(desired_samples, snp_samples)
gene_sample_idx_map = parse_midas_data.calculate_sample_idx_map(desired_samples, gene_samples)
    

same_sample_snp_idxs  = parse_midas_data.apply_sample_index_map_to_indices(snp_sample_idx_map,  desired_same_sample_idxs)  
same_sample_gene_idxs = parse_midas_data.apply_sample_index_map_to_indices(gene_sample_idx_map, desired_same_sample_idxs)  


same_subject_snp_idxs = parse_midas_data.apply_sample_index_map_to_indices(snp_sample_idx_map, desired_same_subject_idxs)  
same_subject_gene_idxs = parse_midas_data.apply_sample_index_map_to_indices(gene_sample_idx_map, desired_same_subject_idxs)  

diff_subject_snp_idxs = parse_midas_data.apply_sample_index_map_to_indices(snp_sample_idx_map, desired_diff_subject_idxs)  
diff_subject_gene_idxs = parse_midas_data.apply_sample_index_map_to_indices(gene_sample_idx_map, desired_diff_subject_idxs)  


# Done calculating... now plot figure!

# clip for log scale
gene_hamming_matrix = numpy.clip(gene_hamming_matrix,0.5,1e09)
total_fixation_matrix = numpy.clip(total_fixation_matrix,0.5,1e09)

pylab.figure(1)
pylab.xlabel('Num substitutions')
pylab.ylabel('Num gene differences')
pylab.ylim([1e-01,1e04])
pylab.xlim([1e-01,1e05])
pylab.title(species_name)

pylab.loglog(total_fixation_matrix[diff_subject_snp_idxs], gene_hamming_matrix[diff_subject_gene_idxs],'r.')
pylab.loglog(total_fixation_matrix[same_subject_snp_idxs], gene_hamming_matrix[same_subject_gene_idxs],'g.')
pylab.plot([1e-01,1e6],[1,1],'k:')
pylab.plot([1,1],[1e-01,1e04],'k:')


pylab.savefig('%s/%s_gene_differences_vs_substitutions.pdf' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight')
pylab.savefig('%s/%s_gene_differences_vs_substitutions.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight',dpi=300)

    
