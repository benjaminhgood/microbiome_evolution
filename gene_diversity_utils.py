import numpy

# For each gene in gene_depth_matrix, calculates # of samples 
# in which it is "present". Returns vector of prevalences
def calculate_gene_prevalences(gene_depth_matrix, marker_coverages, min_copynum=0.5):

    gene_copynum_matrix = gene_depth_matrix * 1.0 / numpy.clip(marker_coverages,1,1e09)
    
    return (gene_copynum_matrix>=min_copynum).sum(axis=1)

# For each sample in gene_depth_matrix, calculates # of genes
# that are "present". Returns vector of gene numbers
def calculate_gene_numbers(gene_depth_matrix, marker_coverages, min_copynum=0.5):

    gene_copynum_matrix = gene_depth_matrix * 1.0 / numpy.clip(marker_coverages,1,1e09)
    
    return (gene_copynum_matrix>=min_copynum).sum(axis=0)
    

# Calculates the number of gene differences between pairs of samples
# Returns: matrix of # of gene differences
#          matrix of # of opportunities 
#          --> we should chat about what the denominator should be!
def calculate_coverage_based_gene_hamming_matrix(gene_depth_matrix, marker_coverages, min_log2_fold_change=4):

    marker_coverages = numpy.clip(marker_coverages,1,1e09)
    gene_copynum_matrix = numpy.clip(gene_depth_matrix,1,1e09)*1.0/marker_coverages[None,:]

    # copynum is between 0.5 and 2
    is_good_copynum = (gene_copynum_matrix>0.5)*(gene_copynum_matrix<2)
  
    # now size is about to get a lot bigger
    num_genes = gene_copynum_matrix.shape[0]
    num_samples = gene_copynum_matrix.shape[1]
    
    gene_hamming_matrix = numpy.zeros((num_samples, num_samples))
    num_opportunities = numpy.zeros((num_samples, num_samples))
    
    # chunk it up into groups of 1000 genes
    chunk_size = 1000
    num_chunks = long(num_genes/chunk_size)+1
    for i in xrange(0, num_chunks):
        
        lower_gene_idx = i*chunk_size
        upper_gene_idx = min([(i+1)*chunk_size, num_genes])
    
        sub_gene_copynum_matrix = gene_copynum_matrix[lower_gene_idx:upper_gene_idx,:] 
        sub_is_good_copynum = is_good_copynum[lower_gene_idx:upper_gene_idx,:] 
    
        is_good_opportunity = numpy.logical_or(sub_is_good_copynum[:,:,None],sub_is_good_copynum[:,None,:])

        fold_change_matrix = numpy.fabs( numpy.log2(sub_gene_copynum_matrix[:,:,None]/sub_gene_copynum_matrix[:,None,:]) ) * is_good_opportunity
        gene_hamming_matrix += (fold_change_matrix>=min_log2_fold_change).sum(axis=0)
        num_opportunities += is_good_opportunity.sum(axis=0)
    
    return gene_hamming_matrix, num_opportunities




def calculate_coverage_based_gene_hamming_matrix_gain_loss(gene_depth_matrix, marker_coverages, min_log2_fold_change=4):
    # in this definition, we keep track of whether there was a gene 'gain' or 'loss' by removing numpy.fabs. This info will be used to plot the gain/loss events between two successive time pints. 

    marker_coverages = numpy.clip(marker_coverages,1,1e09)
    gene_copynum_matrix = numpy.clip(gene_depth_matrix,1,1e09)*1.0/marker_coverages[None,:]

    # copynum is between 0.5 and 2
    is_good_copynum = (gene_copynum_matrix>0.5)*(gene_copynum_matrix<2)
  
    # now size is about to get a lot bigger
    num_genes = gene_copynum_matrix.shape[0]
    num_samples = gene_copynum_matrix.shape[1]
    
    gene_hamming_matrix_gain = numpy.zeros((num_samples, num_samples))
    gene_hamming_matrix_loss = numpy.zeros((num_samples, num_samples))
    num_opportunities = numpy.zeros((num_samples, num_samples))
    
    # chunk it up into groups of 1000 genes
    chunk_size = 1000
    num_chunks = long(num_genes/chunk_size)+1
    for i in xrange(0, num_chunks):
        
        lower_gene_idx = i*chunk_size
        upper_gene_idx = min([(i+1)*chunk_size, num_genes])
    
        sub_gene_copynum_matrix = gene_copynum_matrix[lower_gene_idx:upper_gene_idx,:] 
        sub_is_good_copynum = is_good_copynum[lower_gene_idx:upper_gene_idx,:] 
    
        is_good_opportunity = numpy.logical_or(sub_is_good_copynum[:,:,None],sub_is_good_copynum[:,None,:])

        fold_change_matrix = ( numpy.log2(sub_gene_copynum_matrix[:,:,None]/sub_gene_copynum_matrix[:,None,:]) ) * is_good_opportunity # numpy.fabs removed

        gene_hamming_matrix_gain += (fold_change_matrix>=min_log2_fold_change).sum(axis=0)
        gene_hamming_matrix_loss += (fold_change_matrix<=((-1)*min_log2_fold_change)).sum(axis=0)
        num_opportunities += is_good_opportunity.sum(axis=0)
    
    return gene_hamming_matrix_gain, gene_hamming_matrix_loss, num_opportunities




# Calculate polarized gene copynum changes from i to j that exceed threshold 
# Returns list of differences. Each difference is a tuple of form 
#
# (gene_name, (cov_i, marker_cov_i), (cov_j, marker_cov_j))
#
def calculate_gene_differences_between(i, j, gene_depth_matrix, marker_coverages, min_log2_fold_change=4):

    # Look at these two samples
    gene_depth_matrix = gene_depth_matrix[:,[i,j]]
    marker_coverages = marker_coverages[[i,j]]
    
    marker_coverages = numpy.clip(marker_coverages,1,1e09)
    gene_copynum_matrix = numpy.clip(gene_depth_matrix,1,1e09)*1.0/marker_coverages[None,:]

    gene_differences = []

    # copynum is between 0.5 and 2
    is_good_copynum = (gene_copynum_matrix>0.5)*(gene_copynum_matrix<2)
    
    fold_changes = numpy.fabs(numpy.log2(gene_copynum_matrix[:,1]/gene_copynum_matrix[:,0])) * numpy.logical_or(is_good_copynum[:,0],is_good_copynum[:,1])

    changed_genes = numpy.nonzero(fold_changes>=min_log2_fold_change)[0]

    if len(changed_genes>0):
    
        for gene_idx in changed_genes:
            
            gene_differences.append( (gene_idx, (gene_depth_matrix[gene_idx,0], marker_coverages[0]), (gene_depth_matrix[gene_idx,1], marker_coverages[1])) )
            
    return gene_differences
