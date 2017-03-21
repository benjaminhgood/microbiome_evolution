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


# Load subject and sample metadata
sys.stderr.write("Loading HMP metadata...\n")
subject_sample_map = parse_midas_data.parse_subject_sample_map()
sys.stderr.write("Done!\n")
    
    
# Load SNP information for species_name
sys.stderr.write("Loading %s...\n" % species_name)
snp_samples, allele_counts_map, passed_sites_map = parse_midas_data.parse_snps(species_name, debug)
sys.stderr.write("Done!\n")
    
# Calculate full matrix of synonymous pairwise differences
sys.stderr.write("Calculate synonymous pi matrix...\n")
pi_matrix_syn, avg_pi_matrix_syn = diversity_utils.calculate_pi_matrix(allele_counts_map, passed_sites_map, variant_type='4D')

pi_matrix_syn = numpy.clip(pi_matrix_syn,1e-06,1)
avg_pi_matrix_syn = numpy.clip(avg_pi_matrix_syn,1e-06,1)


# Load gene presence/absence information for species_name
sys.stderr.write("Loading %s...\n" % species_name)
gene_samples, gene_names, gene_presence_matrix, gene_depth_matrix, marker_coverages, gene_reads_matrix = parse_midas_data.parse_pangenome_data(species_name)
sys.stderr.write("Done!\n")
    
total_num_genes = gene_presence_matrix.sum(axis=0)    
avg_num_genes = (total_num_genes[:,None]+total_num_genes[None,:])/2.0
    
# Calculate matrix of number of genes that differ
sys.stderr.write("Calculate gene hamming matrix...\n")
gene_hamming_matrix = diversity_utils.calculate_coverage_based_gene_hamming_matrix(gene_depth_matrix, marker_coverages, min_log2_fold_change=4)
gene_hamming_matrix = numpy.clip(gene_hamming_matrix,0.5,1e05)

# Now need to make the gene samples and snp samples match up
desired_samples = list(set(snp_samples) & set(gene_samples) & set(gene_samples[marker_coverages>40]))
     
     
     
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
pylab.figure(1)
pylab.xlabel('Within-sample $\\pi_s$ (avg)')
pylab.ylabel('Num gene differences')
pylab.ylim([1e-01,1e04])
pylab.xlim([1e-06,1e-01])
pylab.title(species_name)

pylab.loglog(avg_pi_matrix_syn[diff_subject_snp_idxs], gene_hamming_matrix[diff_subject_gene_idxs],'r.')
pylab.loglog(avg_pi_matrix_syn[same_sample_snp_idxs], gene_hamming_matrix[same_sample_gene_idxs],'b.')
pylab.loglog(avg_pi_matrix_syn[same_subject_snp_idxs], gene_hamming_matrix[same_subject_gene_idxs],'g.')
pylab.plot([1e-06,1e-01],[1,1],'k:')


pylab.savefig('%s/%s_gene_differences_vs_pi.pdf' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight')
pylab.savefig('%s/%s_gene_differences_vs_pi.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight',dpi=300)

    
# Done calculating... now plot figure!
pylab.figure(2)
pylab.xlabel('Within-sample $\\pi_s$ (avg)')
pylab.ylabel('Total num genes')
#pylab.ylim([1,1e04])
pylab.xlim([1e-06,1e-01])
pylab.title(species_name)

pylab.loglog(avg_pi_matrix_syn[same_sample_snp_idxs], avg_num_genes[same_sample_gene_idxs],'b.')


pylab.savefig('%s/%s_gene_number_vs_pi.pdf' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight')
pylab.savefig('%s/%s_gene_number_vs_pi.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight',dpi=300)
 

    
