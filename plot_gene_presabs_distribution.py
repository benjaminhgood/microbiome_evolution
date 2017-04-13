import matplotlib  
matplotlib.use('Agg') 
import parse_midas_data
import pylab
import sys
import numpy
import diversity_utils
import gene_diversity_utils
import stats_utils

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
    sys.stderr.write("Usage: python command.py [debug] species_name")
################################################################################

# Load subject and sample metadata
sys.stderr.write("Loading HMP metadata...\n")
subject_sample_map = parse_midas_data.parse_subject_sample_map()
sys.stderr.write("Done!\n")
   
# Load gene presence/absence information for species_name
sys.stderr.write("Loading pangenome data for %s...\n" % species_name)
samples, gene_names, gene_presence_matrix, gene_depth_matrix, marker_coverages, gene_reads_matrix = parse_midas_data.parse_pangenome_data(species_name)
sys.stderr.write("Loaded %d genes across %d samples\n" % gene_depth_matrix.shape)
sys.stderr.write("Done!\n")

min_marker_coverage = 20
high_coverage_samples = samples[marker_coverages>=min_marker_coverage]
sys.stderr.write("Focusing on %d high coverage samples...\n" % len(high_coverage_samples))

# Load metaphlan2 genes
sys.stderr.write("Loading metaphlan2 genes...\n")
metaphlan2_genes = set(parse_midas_data.load_metaphlan2_genes(species_name))   
metaphlan2_gene_idxs = numpy.array([gene_name in metaphlan2_genes for gene_name in gene_names])
sys.stderr.write("Done! (%d genes)\n" % len(metaphlan2_genes))

# Load reference genes
sys.stderr.write("Loading reference genes...\n")
reference_genes = set(parse_midas_data.load_reference_genes(species_name))   
reference_gene_idxs = numpy.array([gene_name in reference_genes for gene_name in gene_names])  
sys.stderr.write("Done! (%d genes)\n" % len(reference_genes))

print reference_genes[0:10]
print gene_names[0:10]
  
# Calculate matrix of number of genes that differ
sys.stderr.write("Calculate gene hamming matrix...\n")
# Either: for all genes in pan-genome
gene_hamming_matrix, num_opportunities = gene_diversity_utils.calculate_coverage_based_gene_hamming_matrix(gene_depth_matrix, marker_coverages, min_log2_fold_change=4)
#
# Or: just the subset from the MIDAS reference genome
#gene_hamming_matrix = diversity_utils.calculate_coverage_based_gene_hamming_matrix(gene_depth_matrix[reference_gene_idxs,:], marker_coverages, min_log2_fold_change=4)
#

sample_idx_map = parse_midas_data.calculate_sample_idx_map(high_coverage_samples, samples)
        
# Calculate which pairs of idxs belong to the same sample, which to the same subject
# and which to different subjects
high_coverage_same_sample_idxs, high_coverage_same_subject_idxs, high_coverage_diff_subject_idxs = parse_midas_data.calculate_subject_pairs(subject_sample_map, high_coverage_samples)

same_sample_idxs = parse_midas_data.apply_sample_index_map_to_indices(sample_idx_map, high_coverage_same_sample_idxs)  
same_subject_idxs = parse_midas_data.apply_sample_index_map_to_indices(sample_idx_map, high_coverage_same_subject_idxs)  
diff_subject_idxs = parse_midas_data.apply_sample_index_map_to_indices(sample_idx_map, high_coverage_diff_subject_idxs)  
     
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
print "Median timepoints=", numpy.median(hamming_timepoints[hamming_timepoints>0.8])
print "Median between=", numpy.median(hamming_between[hamming_between>0.8])

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

pylab.legend(loc='upper right',frameon=False,fontsize=9)
pylab.semilogx([1e-02],[-1])
pylab.ylim([0,1.05])
pylab.xlim([1,1e04])
pylab.savefig('%s/%s_gene_hamming_distribution.pdf' % (parse_midas_data.analysis_directory,species_name), bbox_inches='tight')
pylab.savefig('%s/%s_gene_hamming_distribution.png' % (parse_midas_data.analysis_directory,species_name), bbox_inches='tight',dpi=300)



pylab.figure(3,figsize=(5,3))
pylab.title(species_name,fontsize=11)
pylab.xlabel('Number of samples with gene',fontsize=11)
pylab.ylabel('Fraction of genes',fontsize=11)
pylab.gca().spines['top'].set_visible(False)
pylab.gca().spines['right'].set_visible(False)
pylab.gca().get_xaxis().tick_bottom()
pylab.gca().get_yaxis().tick_left()
  
# Plot gene prevalence SFS
samples, gene_names, gene_presence_matrix, gene_depth_matrix, marker_coverages, gene_reads_matrix

gene_copy_numbers = (gene_depth_matrix/marker_coverages)[:,marker_coverages>=min_marker_coverage]
gene_presence_calls = gene_copy_numbers>0.1
   
gene_prevalence_map = {gene_names[i]: gene_presence_calls[i,:].sum() for i in xrange(0,len(gene_names))}
   
gene_prevalences = (gene_copy_numbers>0.1).sum(axis=1)
stringent_gene_prevalences = (gene_copy_numbers>0.5).sum(axis=1)
lax_gene_prevalences = (gene_copy_numbers>0.01).sum(axis=1)

prevalence_bins = numpy.arange(0,gene_copy_numbers.shape[1]+2)-0.5
prevalences = prevalence_bins[1:]+(prevalence_bins[1]-prevalence_bins[0])/2    
prevalences /= prevalences[-1]

prevalence_counts = numpy.histogram(gene_prevalences, prevalence_bins)[0]*1.0
prevalence_counts /= prevalence_counts.sum()

stringent_prevalence_counts = numpy.histogram(stringent_gene_prevalences, prevalence_bins)[0]*1.0
stringent_prevalence_counts /= stringent_prevalence_counts.sum()

lax_prevalence_counts = numpy.histogram(lax_gene_prevalences, prevalence_bins)[0]*1.0
lax_prevalence_counts /= lax_prevalence_counts.sum()


metaphlan2_prevalence_counts = numpy.histogram(gene_prevalences[metaphlan2_gene_idxs], prevalence_bins)[0]*1.0
metaphlan2_prevalence_counts /= metaphlan2_prevalence_counts.sum()



pylab.semilogy(prevalences, stringent_prevalence_counts, 'k.-',label='All genes (CN>0.5)')
pylab.semilogy(prevalences, prevalence_counts, 'b.-',label='All genes (CN>0.1)')
pylab.semilogy(prevalences, lax_prevalence_counts, 'g.-',label='All genes (CN>0.01)')

   
pylab.legend(loc='upper left',frameon=False,fontsize=8)
#pylab.semilogx([1e-02],[-1])
pylab.ylim([0,1.05])
pylab.xlim([-0.05, 1.05])
pylab.savefig('%s/%s_gene_presence_sfs.pdf' % (parse_midas_data.analysis_directory,species_name), bbox_inches='tight')
pylab.savefig('%s/%s_gene_presence_sfs.png' % (parse_midas_data.analysis_directory,species_name), bbox_inches='tight',dpi=300)
    
# Plot gene prevalence SFS for MIDAS reference genome
pylab.figure(4,figsize=(5,3))
pylab.title(species_name,fontsize=11)
pylab.xlabel('Present in fewer than $n$ samples',fontsize=11)
pylab.ylabel('Genes in reference genome',fontsize=11)
pylab.gca().spines['top'].set_visible(False)
pylab.gca().spines['right'].set_visible(False)
pylab.gca().get_xaxis().tick_bottom()
pylab.gca().get_yaxis().tick_left()

reference_gene_prevalences = (gene_copy_numbers[reference_gene_idxs,:]>0.1).sum(axis=1)
stringent_reference_gene_prevalences = (gene_copy_numbers[reference_gene_idxs,:]>0.5).sum(axis=1)
lax_reference_gene_prevalences = (gene_copy_numbers[reference_gene_idxs,:]>0.01).sum(axis=1)

prevalence_bins = numpy.arange(0,gene_copy_numbers.shape[1]+2)-0.5

reference_prevalence_counts = numpy.histogram(reference_gene_prevalences, prevalence_bins)[0].cumsum()
stringent_reference_prevalence_counts = numpy.histogram(stringent_reference_gene_prevalences, prevalence_bins)[0].cumsum()
lax_reference_prevalence_counts = numpy.histogram(lax_reference_gene_prevalences, prevalence_bins)[0].cumsum()

prevalences = numpy.arange(0,gene_copy_numbers.shape[1]+1)

print reference_prevalence_counts[-1], numpy.histogram(stringent_reference_gene_prevalences, prevalence_bins)[0].sum() 

pylab.semilogy(prevalences, stringent_reference_prevalence_counts, 'k.-',label='CN>0.5')
pylab.semilogy(prevalences, reference_prevalence_counts, 'b.-',label='CN>0.1')
pylab.semilogy(prevalences, lax_reference_prevalence_counts, 'g.-',label='CN>0.01')
pylab.semilogy([10,10],[1,1e04],'k:')
pylab.xlim([-1,prevalences[-1]+1])
pylab.legend(loc='lower right',frameon=False,fontsize=8)
pylab.savefig('%s/%s_reference_gene_presence_sfs.pdf' % (parse_midas_data.analysis_directory,species_name), bbox_inches='tight')
pylab.savefig('%s/%s_reference_gene_presence_sfs.png' % (parse_midas_data.analysis_directory,species_name), bbox_inches='tight',dpi=300)
 
sys.exit(0) 
within_subject_gene_change_prevalences = []
# Calculate gene content differences
for i,j in zip(same_subject_idxs[0],same_subject_idxs[1]):

     gene_differences = gene_diversity_utils.calculate_gene_differences_between(i, j, gene_names, gene_depth_matrix, marker_coverages, min_log2_fold_change=4)

     if len(gene_differences)>0:
         print "Differences between", i, j
         for idx in xrange(0,len(gene_differences)):
             print gene_differences[idx], gene_prevalence_map[gene_differences[idx][0]]
             