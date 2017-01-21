import numpy

def calculate_rsquared(allele_counts_1, allele_counts_2):
    
    # allele counts = 1 x samples x alleles vector
    
    depths_1 = allele_counts_1.sum(axis=2)
    freqs_1 = allele_counts_1[:,:,0]*1.0/(depths_1+(depths_1==0))
    depths_2 = allele_counts_2.sum(axis=2)
    freqs_2 = allele_counts_2[:,:,0]*1.0/(depths_2+(depths_2==0))
    
    # consensus approximation
    freqs_1 = numpy.around(freqs_1)
    freqs_2 = numpy.around(freqs_2)
    
    
    joint_passed_sites = (depths_1>0)[None,:,:]*(depths_2>0)[:,None,:]
    # sites x sites x samples matrix
    
    joint_freqs = freqs_1[None,:,:]*freqs_2[:,None,:]
    # sites x sites x samples_matrix
    
    total_joint_passed_sites = joint_passed_sites.sum(axis=2)
    total_joint_passed_sites = total_joint_passed_sites+(total_joint_passed_sites==0)
    
    joint_pooled_freqs = (joint_freqs*joint_passed_sites).sum(axis=2)/total_joint_passed_sites
   
    joint_pooled_freqs *= (joint_pooled_freqs>1e-10)
    
    marginal_pooled_freqs_1 = (freqs_1[None,:,:]*joint_passed_sites).sum(axis=2)/total_joint_passed_sites
    marginal_pooled_freqs_1 *= (marginal_pooled_freqs_1>1e-10)
    
    marginal_pooled_freqs_2 = (freqs_2[:,None,:]*joint_passed_sites).sum(axis=2)/total_joint_passed_sites 
    marginal_pooled_freqs_2 *= (marginal_pooled_freqs_2>1e-10)
       
    rsquared_numerators = numpy.square(joint_pooled_freqs-marginal_pooled_freqs_1*marginal_pooled_freqs_2)
    
    
    rsquared_denominators = marginal_pooled_freqs_1*(1-marginal_pooled_freqs_1)*marginal_pooled_freqs_2*(1-marginal_pooled_freqs_2)
    
    rsquareds = rsquared_numerators/(rsquared_denominators+(rsquared_denominators==0))
    
    return rsquared_numerators, rsquared_denominators

def calculate_sample_freqs(allele_counts_map, passed_sites_map, variant_type='4D', allowed_genes=None):

    if allowed_genes == None:
        allowed_genes = set(passed_sites_map.keys())
     
    sample_freqs = [[] for i in xrange(0,allele_counts_map[allele_counts_map.keys()[0]][variant_type]['alleles'].shape[1])]
    
    passed_sites = numpy.zeros(passed_sites_map[passed_sites_map.keys()[0]][variant_type]['sites'].shape[0])*1.0
    
    for gene_name in allowed_genes:
    
        allele_counts = allele_counts_map[gene_name][variant_type]['alleles']

        if len(allele_counts)==0:
            continue
            
        depths = allele_counts.sum(axis=2)
        freqs = allele_counts[:,:,0]/(depths+(depths==0))
        freqs = numpy.fmin(freqs,1-freqs)
        for sample_idx in xrange(0,freqs.shape[1]):
            gene_freqs = freqs[:,sample_idx]
            sample_freqs[sample_idx].extend( gene_freqs[gene_freqs>0])
            
        passed_sites += numpy.diagonal(passed_sites_map[gene_name][variant_type]['sites'])
        
    
    return sample_freqs, passed_sites

        
def calculate_pooled_freqs(allele_counts_map, passed_sites_map,  variant_type='4D', allowed_genes=None):

    if allowed_genes == None:
        allowed_genes = set(passed_sites_map.keys())
     
    pooled_freqs = []
    
    for gene_name in allowed_genes:
    
        allele_counts = allele_counts_map[gene_name][variant_type]['alleles']

        if len(allele_counts)==0:
            continue
            
        depths = allele_counts.sum(axis=2)
        freqs = allele_counts/(depths+(depths==0))[:,:,None]
        gene_pooled_freqs = freqs[:,:,0].sum(axis=1)/(depths>0).sum(axis=1)
        pooled_freqs.extend(gene_pooled_freqs)
    
    pooled_freqs = numpy.array(pooled_freqs)
    return pooled_freqs

     

def calculate_fixation_matrix(allele_counts_map, passed_sites_map, variant_type='4D', allowed_genes=None, min_freq=0, min_change=0.8):

    if allowed_genes == None:
        allowed_genes = set(passed_sites_map.keys())
        
    fixation_matrix = numpy.zeros_like(passed_sites_map[passed_sites_map.keys()[0]][variant_type]['sites'])*1.0
    
    for gene_name in allowed_genes:
    
        allele_counts = allele_counts_map[gene_name][variant_type]['alleles']

        if len(allele_counts)==0:
            continue

        depths = allele_counts.sum(axis=2)
        alt_freqs = allele_counts[:,:,0]/(depths+(depths==0))
        alt_freqs[alt_freqs<min_freq] = 0.0
        alt_freqs[alt_freqs>(1-min_freq)] = 1.0
        passed_depths = (depths>0)[:,:,None]*(depths>0)[:,None,:]
    
        delta_freq = numpy.fabs(alt_freqs[:,:,None]-alt_freqs[:,None,:])
        delta_freq[passed_depths==0] = 0
        delta_freq[delta_freq<min_change] = 0
        
        fixation_matrix += delta_freq.sum(axis=0)
        
    return fixation_matrix
    
def calculate_pi_matrix(allele_counts_map, passed_sites_map, variant_type='4D', allowed_genes=None):

    if allowed_genes == None:
        allowed_genes = set(passed_sites_map.keys())
        
    pi_matrix = numpy.zeros_like(passed_sites_map[passed_sites_map.keys()[0]][variant_type]['sites'])*1.0
    avg_pi_matrix = numpy.zeros_like(pi_matrix)
    passed_sites = numpy.zeros_like(pi_matrix)
    
    for gene_name in allowed_genes:
        
        #print passed_sites_map[gene_name][variant_type].shape, passed_sites.shape
        #print gene_name, variant_type
        
        passed_sites += passed_sites_map[gene_name][variant_type]['sites']
           
        allele_counts = allele_counts_map[gene_name][variant_type]['alleles']

        if len(allele_counts)==0:
            continue
         

        depths = allele_counts.sum(axis=2)
        freqs = allele_counts/(depths+(depths==0))[:,:,None]
        self_freqs = (allele_counts-1)/(depths-1+2*(depths==0))[:,:,None]
        self_pis = ((depths>0)-(freqs*self_freqs).sum(axis=2))
             
        I,J = depths.shape
    
        # pi between sample j and sample l
        gene_pi_matrix = numpy.einsum('ij,il',(depths>0)*1.0,(depths>0)*1.0)-numpy.einsum('ijk,ilk',freqs,freqs)
    
        # average of pi within sample j and within sample i
        gene_avg_pi_matrix = (numpy.einsum('ij,il',self_pis,(depths>0)*1.0)+numpy.einsum('ij,il',(depths>0)*1.0,self_pis))/2
    
        diagonal_idxs = numpy.diag_indices(J)
        gene_pi_matrix[diagonal_idxs] = gene_avg_pi_matrix[diagonal_idxs]
    
        pi_matrix += gene_pi_matrix
        avg_pi_matrix += gene_avg_pi_matrix
        
    pi_matrix = pi_matrix/(passed_sites+(passed_sites==0))
    avg_pi_matrix = avg_pi_matrix/(passed_sites+(passed_sites==0))
    
    return pi_matrix, avg_pi_matrix

    
def phylip_distance_matrix_str(matrix, samples):
    
    lines = [str(len(samples))]
    for i in xrange(0,len(samples)):
        lines.append( "\t".join([samples[i]]+["%g" % x for x in matrix[i,:]]))
    
    return "\n".join(lines)
    
import numpy
from scipy.special import gammaln as loggamma

def fold_sfs(fs):
    n = len(fs)+1
    folded_fs = (fs + fs[::-1])[0:(n-1)/2]
    if (n-1) % 2 != 0:
        folded_fs[-1] *= 0.5
    return folded_fs


def estimate_sfs_naive_binning(allele_counts, target_depth=10):

    depths = allele_counts.sum(axis=1)
    
    allele_counts = allele_counts[depths>0]
    depths = depths[depths>0]
    
    freqs = allele_counts[:,0]/depths
    
    bins = (numpy.arange(0,target_depth+2)-0.5)/target_depth
    
    counts,dummy = numpy.histogram(freqs,bins)
    
    return counts

def estimate_sfs_downsampling(allele_counts, target_depth=10):
    
    depths = allele_counts.sum(axis=1)
    
    allele_counts = allele_counts[depths>0]
    depths = depths[depths>0]
    
    Dmin = min([depths.min(),target_depth]) # this is what we have to downsample to
    # if you don't like it, send us an allele_counts matrix
    # that has been thresholded to a higher min value
    
    count_density = numpy.zeros(Dmin+1)*1.0

    
    A = numpy.outer(allele_counts[:,0], numpy.ones(Dmin+1))
    D = numpy.outer(depths, numpy.ones(Dmin+1))
    ks = numpy.outer(numpy.ones_like(depths), numpy.arange(0,Dmin+1))
    
    count_density = numpy.exp(loggamma(A+1)-loggamma(A-ks+1)-loggamma(ks+1) + loggamma(D-A+1)-loggamma(D-A-(Dmin-ks)+1)-loggamma(Dmin-ks+1) + loggamma(D-Dmin+1) + loggamma(Dmin+1) - loggamma(D+1)).sum(axis=0)
    
    return count_density
