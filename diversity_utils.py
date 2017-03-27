import numpy

def calculate_rsquared_condition_freq(allele_counts_1, allele_counts_2, low_freq, high_freq):
    # Note: should actually be sigma_squared! 
    # sigma_squared= E[X]/E[Y], where X=(p_ab-pa*pb)^2 and Y=(pa*(1-pa)*pb*(1-pb))
    # rsquared=E[X/Y]
    # see McVean 2002 for more notes on the difference. 

    # allele counts = 1 x samples x alleles vector
    
    depths_1 = allele_counts_1.sum(axis=2)
    freqs_1 = allele_counts_1[:,:,0]*1.0/(depths_1+(depths_1==0))
    depths_2 = allele_counts_2.sum(axis=2)
    freqs_2 = allele_counts_2[:,:,0]*1.0/(depths_2+(depths_2==0))

    
    # consensus approximation
    freqs_1 = numpy.around(freqs_1)
    freqs_2 = numpy.around(freqs_2)

    # condition on allele frequency in the pooled population:
    pooled_freqs_1=freqs_1[:,:].sum(axis=1)/len(freqs_1[0])
    pooled_freqs_2=freqs_2[:,:].sum(axis=1)/len(freqs_2[0])

    # check if any freqs >0.5, if so, fold:
    pooled_freqs_1=numpy.where(pooled_freqs_1 > 0.5, 1-pooled_freqs_1, pooled_freqs_1)
    pooled_freqs_2=numpy.where(pooled_freqs_2 > 0.5, 1-pooled_freqs_2, pooled_freqs_2) 

    # this asks which pairs of sites have depths >0 at BOTH sites as well as which paris of sites both have pooled frequencies within the low_freq and high_freq ranges. 
    # None here takes the product of the elements in the two vectors and returns a matrix. 

    
    passed_sites_1=(depths_1>0)*(pooled_freqs_1 >= low_freq)[:,None]*(pooled_freqs_1 <=high_freq)[:,None]
    passed_sites_2=(depths_2>0)*(pooled_freqs_2 >= low_freq)[:,None]*(pooled_freqs_2 <= high_freq)[:,None]
    joint_passed_sites=passed_sites_1[None,:,:]*passed_sites_2[:,None,:]
    # sites x sites x samples matrix
    
    joint_freqs = freqs_1[None,:,:]*freqs_2[:,None,:]
    # sites x sites x samples_matrix
    
    # this tells us what the denominator is for the computation below for joint_pooled_freqs
    total_joint_passed_sites = joint_passed_sites.sum(axis=2)
    # add 1 to denominator if some pair is 0. 
    total_joint_passed_sites = total_joint_passed_sites+(total_joint_passed_sites==0)
    
    # compute p_ab
    joint_pooled_freqs = (joint_freqs*joint_passed_sites).sum(axis=2)/total_joint_passed_sites   
    # floting point issue
    joint_pooled_freqs *= (joint_pooled_freqs>1e-10)
    
    # compute p_a
    marginal_pooled_freqs_1 = (freqs_1[None,:,:]*joint_passed_sites).sum(axis=2)/total_joint_passed_sites
    marginal_pooled_freqs_1 *= (marginal_pooled_freqs_1>1e-10)

    # compute p_b
    marginal_pooled_freqs_2 = (freqs_2[:,None,:]*joint_passed_sites).sum(axis=2)/total_joint_passed_sites 
    marginal_pooled_freqs_2 *= (marginal_pooled_freqs_2>1e-10)
       
    # (p_ab-p_a*p_b)^2
    rsquared_numerators = numpy.square(joint_pooled_freqs-marginal_pooled_freqs_1*marginal_pooled_freqs_2)
    
    # (p_a*(1-p_a)*pb*(1-p_b))
    rsquared_denominators = marginal_pooled_freqs_1*(1-marginal_pooled_freqs_1)*marginal_pooled_freqs_2*(1-marginal_pooled_freqs_2)


    rsquareds = rsquared_numerators/(rsquared_denominators+(rsquared_denominators==0))
    
    return rsquared_numerators, rsquared_denominators


#####################################################################

def calculate_unbiased_sigmasquared(allele_counts_1, allele_counts_2):
    
    # Note: should actually be sigma_squared! 
    # sigma_squared= E[X]/E[Y], where X=(p_ab-pa*pb)^2 and Y=(pa*(1-pa)*pb*(1-pb))
    # rsquared=E[X/Y]
    # see McVean 2002 for more notes on the difference. 

    # allele counts = 1 x samples x alleles vector
    
    depths_1 = allele_counts_1.sum(axis=2)
    freqs_1 = allele_counts_1[:,:,0]*1.0/(depths_1+(depths_1==0))
    depths_2 = allele_counts_2.sum(axis=2)
    freqs_2 = allele_counts_2[:,:,0]*1.0/(depths_2+(depths_2==0))
    
    # consensus approximation
    freqs_1 = numpy.around(freqs_1)
    freqs_2 = numpy.around(freqs_2)
    

    # this asks which pairs of sites have depths >0 at BOTH sites
    # None here takes the product of the elements in the two vectors and returns a matrix. 
    joint_passed_sites = (depths_1>0)[None,:,:]*(depths_2>0)[:,None,:]
    # sites x sites x samples matrix
    
    joint_freqs = freqs_1[None,:,:]*freqs_2[:,None,:]
    # sites x sites x samples_matrix
    
    # this tells us what the denominator is for the computation below for joint_pooled_freqs
    total_joint_passed_sites = joint_passed_sites.sum(axis=2)
    # add 1 to denominator if some pair is 0. 
    total_joint_passed_sites = total_joint_passed_sites+(total_joint_passed_sites==0)
    
    # compute p_ab
    joint_pooled_freqs = (joint_freqs*joint_passed_sites).sum(axis=2)/total_joint_passed_sites   
    # floting point issue
    joint_pooled_freqs *= (joint_pooled_freqs>1e-10)
    
    # compute p_a
    marginal_pooled_freqs_1 = (freqs_1[None,:,:]*joint_passed_sites).sum(axis=2)/total_joint_passed_sites
    marginal_pooled_freqs_1 *= (marginal_pooled_freqs_1>1e-10)

    # compute p_b
    marginal_pooled_freqs_2 = (freqs_2[:,None,:]*joint_passed_sites).sum(axis=2)/total_joint_passed_sites 
    marginal_pooled_freqs_2 *= (marginal_pooled_freqs_2>1e-10)
       
    # (p_ab-p_a*p_b)^2
    rsquared_numerators = numpy.square(joint_pooled_freqs-marginal_pooled_freqs_1*marginal_pooled_freqs_2)
    
    # (p_a*(1-p_a)*pb*(1-p_b))
    rsquared_denominators = marginal_pooled_freqs_1*(1-marginal_pooled_freqs_1)*marginal_pooled_freqs_2*(1-marginal_pooled_freqs_2)
    
    rsquareds = rsquared_numerators/(rsquared_denominators+(rsquared_denominators==0))
    
    return rsquared_numerators, rsquared_denominators


#####################################################################

def calculate_rsquared(allele_counts_1, allele_counts_2):
    # Note: should actually be sigma_squared! 
    # sigma_squared= E[X]/E[Y], where X=(p_ab-pa*pb)^2 and Y=(pa*(1-pa)*pb*(1-pb))
    # rsquared=E[X/Y]
    # see McVean 2002 for more notes on the difference. 

    # allele counts = 1 x samples x alleles vector
    
    depths_1 = allele_counts_1.sum(axis=2)
    freqs_1 = allele_counts_1[:,:,0]*1.0/(depths_1+(depths_1==0))
    depths_2 = allele_counts_2.sum(axis=2)
    freqs_2 = allele_counts_2[:,:,0]*1.0/(depths_2+(depths_2==0))

    
    # consensus approximation
    freqs_1 = numpy.around(freqs_1)
    freqs_2 = numpy.around(freqs_2)
    

    # this asks which pairs of sites have depths >0 at BOTH sites
    # None here takes the product of the elements in the two vectors and returns a matrix. 
    joint_passed_sites = (depths_1>0)[None,:,:]*(depths_2>0)[:,None,:]
    # sites x sites x samples matrix
    
    joint_freqs = freqs_1[None,:,:]*freqs_2[:,None,:]
    # sites x sites x samples_matrix
    
    # this tells us what the denominator is for the computation below for joint_pooled_freqs
    total_joint_passed_sites = joint_passed_sites.sum(axis=2)
    # add 1 to denominator if some pair is 0. 
    total_joint_passed_sites = total_joint_passed_sites+(total_joint_passed_sites==0)
    
    # compute p_ab
    joint_pooled_freqs = (joint_freqs*joint_passed_sites).sum(axis=2)/total_joint_passed_sites   
    # floting point issue
    joint_pooled_freqs *= (joint_pooled_freqs>1e-10)
    
    # compute p_a
    marginal_pooled_freqs_1 = (freqs_1[None,:,:]*joint_passed_sites).sum(axis=2)/total_joint_passed_sites
    marginal_pooled_freqs_1 *= (marginal_pooled_freqs_1>1e-10)

    # compute p_b
    marginal_pooled_freqs_2 = (freqs_2[:,None,:]*joint_passed_sites).sum(axis=2)/total_joint_passed_sites 
    marginal_pooled_freqs_2 *= (marginal_pooled_freqs_2>1e-10)
       
    # (p_ab-p_a*p_b)^2
    rsquared_numerators = numpy.square(joint_pooled_freqs-marginal_pooled_freqs_1*marginal_pooled_freqs_2)
    
    # (p_a*(1-p_a)*pb*(1-p_b))
    rsquared_denominators = marginal_pooled_freqs_1*(1-marginal_pooled_freqs_1)*marginal_pooled_freqs_2*(1-marginal_pooled_freqs_2)
    
    rsquareds = rsquared_numerators/(rsquared_denominators+(rsquared_denominators==0))
    
    return rsquared_numerators, rsquared_denominators





##################################
def generate_haplotype(allele_counts_4D, allele_counts_1D, location_dictionary, species_name):

    freqs={}

    depths_4D = allele_counts_4D.sum(axis=2)
    freqs['4D'] = allele_counts_4D[:,:,0]*1.0/(depths_4D+(depths_4D==0))

    depths_1D = allele_counts_1D.sum(axis=2)
    freqs['1D'] = allele_counts_1D[:,:,0]*1.0/(depths_1D+(depths_1D==0))

    #explanation of numpy commands above:
    # allele_counts_1.sum(axis=2) this returns a sum over all sites alt + ref counts. 
    #(depths_1+(depths_1==0) this is done because if depths_1==0, then we've have a division error. addition of 1 when depths_1==0. 
    #allele_counts_1[:,:,0] means that the alt allele is grabbed. Multiply by 1.0 to convert to float
    
    # consensus approximation
    consensus={}
    consensus['4D'] = numpy.around(freqs['4D'])
    consensus['1D'] = numpy.around(freqs['1D'])
    

    locations=location_dictionary.keys()
    locations=sorted(locations)
   
    #s_consensus='' # store the haplotypes in a string for printing out later
    #s_annotation=''
    outFile_consensus=open('tmp_consensus_%s.txt' % species_name ,'w')
    outFile_anno=open('tmp_anno_%s.txt' % species_name ,'w')

    for loc in range(0, len(locations)):
        location=str(int(locations[loc])) 
        index=location_dictionary[locations[loc]][0]
        variant_type=location_dictionary[locations[loc]][1]
        alleles=consensus[variant_type][index].tolist()
        annotation=freqs[variant_type][index].tolist()

        for person in range(0, len(alleles)):
            alleles[person]=str(int(alleles[person]))
            if annotation[person] ==0:
                annotation[person]=str(0) # no difference from ref
            elif annotation[person] ==1:
                if variant_type=='4D':
                    annotation[person]=str(1) # fixed syn diff from ref
                else:
                    annotation[person]=str(2) # fixed nonsyn diff from ref
            else: 
                if variant_type=='4D':
                    annotation[person]=str(3) # polymorphic syn within host
                else:
                    annotation[person]=str(4) # polymorphic nonsyn within host
        s_consensus = location + ',' + ','.join(alleles) +'\n' 
        s_annotation = location + ',' + ','.join(annotation) + '\n'
        outFile_consensus.write(s_consensus)
        outFile_anno.write(s_annotation)

    return [s_consensus, s_annotation]

####################################

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



def calculate_fixation_matrix(allele_counts_map, passed_sites_map, allowed_variant_types=set([]), allowed_genes=set([]), min_freq=0, min_change=0.8):

    total_genes = set(passed_sites_map.keys())

    if len(allowed_genes)==0:
        allowed_genes = set(passed_sites_map.keys())
    
    allowed_genes = (allowed_genes & total_genes)     
    
    if len(allowed_variant_types)==0:
        allowed_variant_types = set(['1D','2D','3D','4D'])    
                    
    fixation_matrix = numpy.zeros_like(passed_sites_map.values()[0].values()[0]['sites'])*1.0  
    passed_sites = numpy.zeros_like(fixation_matrix)*1.0
    
    for gene_name in allowed_genes:
        
        for variant_type in passed_sites_map[gene_name].keys():
             
            if variant_type not in allowed_variant_types:
                continue
        
            passed_sites += passed_sites_map[gene_name][variant_type]['sites']
   
            allele_counts = allele_counts_map[gene_name][variant_type]['alleles']                        
            if len(allele_counts)==0:
                continue
            

            depths = allele_counts.sum(axis=2)
            alt_freqs = allele_counts[:,:,0]/(depths+(depths==0))
            alt_freqs[alt_freqs<min_freq] = 0.0
            alt_freqs[alt_freqs>=(1-min_freq)] = 1.0
            passed_depths = (depths>0)[:,:,None]*(depths>0)[:,None,:]
    
            delta_freq = numpy.fabs(alt_freqs[:,:,None]-alt_freqs[:,None,:])
            delta_freq[passed_depths==0] = 0
            delta_freq[delta_freq<min_change] = 0
        
            fixation_matrix += delta_freq.sum(axis=0)
        
    return fixation_matrix, passed_sites  

   
def calculate_pi_matrix(allele_counts_map, passed_sites_map, variant_type='4D', allowed_genes=None):

    if allowed_genes == None:
        allowed_genes = set(passed_sites_map.keys())
        
    pi_matrix = numpy.zeros_like(passed_sites_map[passed_sites_map.keys()[0]][variant_type]['sites'])*1.0
    avg_pi_matrix = numpy.zeros_like(pi_matrix)
    passed_sites = numpy.zeros_like(pi_matrix)
    
    for gene_name in allowed_genes:
        
        if gene_name in passed_sites_map.keys():
            #print passed_sites_map[gene_name][variant_type].shape, passed_sites.shape
            #print gene_name, variant_type
        
            passed_sites += passed_sites_map[gene_name][variant_type]['sites']
           
            allele_counts = allele_counts_map[gene_name][variant_type]['alleles']

            if len(allele_counts)==0:
                continue
         

            depths = allele_counts.sum(axis=2)
            freqs = allele_counts/(depths+(depths<0.1))[:,:,None]
            self_freqs = (allele_counts-1)/(depths-1+2*(depths<1.1))[:,:,None]
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
    
    

# Calculate polarized SNP changes from i to j that exceed threshold 
# Returns list of differences. Each difference is a tuple of form 
#
# (gene_name, (contig, location), (alt_i, depth_i), (alt_j, depth_j))
#
def calculate_snp_differences_between(i,j,allele_counts_map, passed_sites_map, allowed_variant_types=set([]), allowed_genes=set([]), min_freq=0, min_change=0.8):

    if len(allowed_genes)==0:
        allowed_genes = set(passed_sites_map.keys())
        
    if len(allowed_variant_types)==0:
        allowed_variant_types = set(['1D','2D','3D','4D'])    
    
    snp_changes = []
        
    for gene_name in allowed_genes:
        
        if gene_name not in allele_counts_map.keys():
            continue
            
        for variant_type in allele_counts_map[gene_name].keys():
            
            if variant_type not in allowed_variant_types:
                continue

            allele_counts = allele_counts_map[gene_name][variant_type]['alleles']
                        
            if len(allele_counts)==0:
                continue

            allele_counts = allele_counts[:,[i,j],:]
            depths = allele_counts.sum(axis=2)
            alt_freqs = allele_counts[:,:,0]/(depths+(depths==0))
            alt_freqs[alt_freqs<min_freq] = 0.0
            alt_freqs[alt_freqs>=(1-min_freq)] = 1.0
            
            passed_depths = (depths>0)[:,:,None]*(depths>0)[:,None,:]
    
            passed_depths = (depths>0)[:,0]*(depths>0)[:,1]
            
            delta_freqs = numpy.fabs(alt_freqs[:,1]-alt_freqs[:,0])
            delta_freqs[passed_depths==0] = 0
            delta_freqs[delta_freqs<min_change] = 0
    
            changed_sites = numpy.nonzero(delta_freqs)[0]
            
            if len(changed_sites)>0:
                # some fixations!
                
                for idx in changed_sites:
                    snp_changes.append((gene_name, allele_counts_map[gene_name][variant_type]['locations'][idx], (allele_counts[idx,0], depths[idx,0]), (allele_counts[idx,1],depths[idx,1]) ))
                        
    return snp_changes

