import numpy
from scipy.linalg import eigh
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster.hierarchy import cophenet
from scipy.cluster.hierarchy import fcluster
from numpy.random import shuffle, normal
import scipy.stats
from scipy.stats import binom
import config
from scipy.special import betainc
import sys
import parse_midas_data
import sample_utils
import stats_utils
import os.path
import sfs_utils

# Calls consensus genotypes from matrix of allele counts
#
# Returns: genotype matrix, passed_sitse matrix for polymorphic sites
#
def calculate_consensus_genotypes(allele_counts_matrix,lower_threshold=0.2,upper_threshold=0.8):
    
    num_sites, num_samples, num_alleles = allele_counts_matrix.shape
    
    depths = allele_counts_matrix.sum(axis=2)
    freqs = allele_counts_matrix[:,:,0]*1.0/(depths+(depths==0))
    passed_sites_matrix = (depths>0)*numpy.logical_or(freqs<=lower_threshold,freqs>=upper_threshold)
    # consensus approximation
    genotype_matrix = numpy.around(freqs)*passed_sites_matrix
    
    
    return genotype_matrix, passed_sites_matrix
    
    
def calculate_consensus_polymorphic_genotypes(allele_counts_matrix,lower_threshold=0.2,upper_threshold=0.8):
    
    genotype_matrix, passed_sites_matrix =  calculate_consensus_genotypes(allele_counts_matrix,lower_threshold,upper_threshold) 
    
    prevalences = (genotype_matrix*passed_sites_matrix).sum(axis=1)
    min_prevalences = 0.5
    max_prevalences = (passed_sites_matrix).sum(axis=1)-0.5
    
    polymorphic_sites = (prevalences>min_prevalences)*(prevalences<max_prevalences)
    
    return genotype_matrix[polymorphic_sites,:], passed_sites_matrix[polymorphic_sites,:]


# Calculates first two PCA coordinates for samples in allele_counts
# using the normalization scheme outlined in McVean (PLoS Genet, 2009).
#
# Returns: (vector of pca1 coords, vector of pca2 coords), (percent variance 1, percent variance 2)
#
def calculate_pca_coordinates(genotype_matrix, passed_sites_matrix):

    Zl = (genotype_matrix*passed_sites_matrix).sum(axis=1)/(passed_sites_matrix).sum(axis=1)

    Zli = (genotype_matrix-Zl[:,None])*passed_sites_matrix
    
    Mij = numpy.einsum('li,lj',Zli,Zli)/numpy.einsum('li,lj',passed_sites_matrix, passed_sites_matrix)

    # calculate eigenvectors & eigenvalues of the covariance matrix
    # use 'eigh' rather than 'eig' since R is symmetric, 
    # the performance gain is substantial
    evals, evecs = eigh(Mij)

    # sort eigenvalue in decreasing order
    idx = numpy.argsort(evals)[::-1]
    evals = evals[idx]
    evecs = evecs[:,idx]
    
    variances = evals/evals.sum()
    
    pca1_coords = evals[0]**0.5*evecs[:,0]
    pca2_coords = evals[1]**0.5*evecs[:,1]
    
    return (pca1_coords, pca2_coords), (variances[0],variances[1])
    
 

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


#####################################################################
def calculate_sigmasquared(allele_counts_1, allele_counts_2):
    # A standard measure of linkage disequilibrium:
    #
    # sigma_squared= E[X]/E[Y], where X=(p_ab-pa*pb)^2 and Y=(pa*(1-pa)*pb*(1-pb))
    # rsquared=E[X/Y]
    # see McVean 2002 for more notes on the difference. 

    # allele counts = 1 x samples x alleles vector
    
    freqs_1, passed_sites_1 = calculate_consensus_genotypes(allele_counts_1)
    freqs_2, passed_sites_2 = calculate_consensus_genotypes(allele_counts_2)
    
    # this asks which pairs of sites have depths >0 at BOTH sites
    # None here takes the product of the elements in the two vectors and returns a matrix. 
    joint_passed_sites = (passed_sites_1)[None,:,:]*(passed_sites_2)[:,None,:]
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


#####################################################################
def calculate_unbiased_sigmasquared(allele_counts_1, allele_counts_2):
    # An alternate version of a standard measure of linkage disequilibrium:
    #
    # sigma_squared= E[X]/E[Y], where X=(p_ab-pa*pb)^2 and Y=(pa*(1-pa)*pb*(1-pb))
    # rsquared=E[X/Y]
    # see McVean 2002 for more notes on the difference. 
    #
    # where we have corrected for finite sample effects

    genotypes_1, passed_sites_1 = calculate_consensus_genotypes(allele_counts_1)
    genotypes_2, passed_sites_2 = calculate_consensus_genotypes(allele_counts_2)
    
    
    # this asks which pairs of sites have depths >0 at BOTH sites
    # None here takes the product of the elements in the two vectors and returns a matrix. 
    joint_passed_sites = (passed_sites_1)[None,:,:]*(passed_sites_2)[:,None,:]
    # sites x sites x samples matrix
    
    # allele counts
    ns = joint_passed_sites.sum(axis=2)
    
    n11s = ((genotypes_1[None,:,:])*(genotypes_2[:,None,:])*joint_passed_sites).sum(axis=2)
    n10s = (genotypes_1[None,:,:]*(1-genotypes_2[:,None,:])*joint_passed_sites).sum(axis=2)
    n01s = ((1-genotypes_1[None,:,:])*(genotypes_2[:,None,:])*joint_passed_sites).sum(axis=2)
    n00s = ((1-genotypes_1[None,:,:])*(1-genotypes_2[:,None,:])*joint_passed_sites).sum(axis=2)
    
    #print "Gene:" 
    #print n11s
    #print n10s
    #print n01s
    #print n00s
    #print "--"

    
    # First calculate numerator
    rsquared_numerators = n11s*(n11s-1)*n00s*(n00s-1)
    rsquared_numerators -= 2*n10s*n01s*n11s*n00s
    rsquared_numerators += n10s*(n10s-1)*n01s*(n01s-1)
    
    #print "Before divide:"
    #print rsquared_numerators
    
    rsquared_numerators = rsquared_numerators*(ns>3.5)*1.0/(ns*(ns-1)*(ns-2)*(ns-3)+10*(ns<3.5))
    
    #print "After divide:"
    #print rsquared_numerators
    #print "---"
    # Now calculate denominator 
    # (more annoying... there are 16 terms rather than 4, so we will write them separately)
    
    #1
    rsquared_denominators = n10s*(n10s-1)*n01s*(n01s-1)    
    #2
    rsquared_denominators += n10s*n01s*(n01s-1)*n00s        
    #3
    rsquared_denominators += n10s*(n10s-1)*n01s*n11s         
    #4
    rsquared_denominators += n10s*n01s*n11s*n00s           
    #5
    rsquared_denominators += n10s*(n10s-1)*n01s*n00s       
    #6
    rsquared_denominators += n10s*n01s*n00s*(n00s-1)
    #7
    rsquared_denominators += n10s*(n10s-1)*n11s*n00s
    #8
    rsquared_denominators += n10s*n11s*n00s*(n00s-1)
    #9
    rsquared_denominators += n10s*n01s*(n01s-1)*n11s
    #10
    rsquared_denominators += n01s*(n01s-1)*n11s*n00s
    #11
    rsquared_denominators += n10s*n01s*n11s*(n11s-1)
    #12
    rsquared_denominators += n01s*n11s*(n11s-1)*n00s
    #13
    rsquared_denominators += n10s*n01s*n11s*n00s
    #14
    rsquared_denominators += n01s*n11s*n00s*(n00s-1)
    #15
    rsquared_denominators += n10s*n11s*(n11s-1)*n00s
    #16
    rsquared_denominators += n11s*(n11s-1)*n00s*(n00s-1)
    
    # divide by sample size
    rsquared_denominators = rsquared_denominators*(ns>3.5)*1.0/(ns*(ns-1)*(ns-2)*(ns-3)+10*(ns<3.5))
    
    return rsquared_numerators, rsquared_denominators




##################################
def generate_haplotype(allele_counts_4D, allele_counts_1D, location_dictionary, species_name):

    freqs={}
    depths={}

    depths['4D'] = allele_counts_4D.sum(axis=2)
    freqs['4D'] = allele_counts_4D[:,:,0]*1.0/(depths['4D']+(depths['4D']==0))

    depths['1D'] = allele_counts_1D.sum(axis=2)
    freqs['1D'] = allele_counts_1D[:,:,0]*1.0/(depths['1D']+(depths['1D']==0))
    
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
    outFile_consensus=open(os.path.expanduser('~/tmp_intermediate_files/tmp_consensus_%s.txt') % species_name ,'w')
    outFile_anno=open(os.path.expanduser('~/tmp_intermediate_files/tmp_anno_%s.txt') % species_name ,'w')

    for loc in range(0, len(locations)):
        location=str(int(locations[loc])) 
        index=location_dictionary[locations[loc]][0]
        variant_type=location_dictionary[locations[loc]][1]
        alleles=consensus[variant_type][index].tolist()
        annotation=freqs[variant_type][index].tolist()
        coverage=depths[variant_type][index].tolist() # if coverage ==0, then set to 'N' in both the consensus and annotation files. 

        for person in range(0, len(alleles)):
            alleles[person]=str(int(alleles[person]))
            if coverage[person] == 0.0:
                alleles[person]='N'
                annotation[person]='N'
            else:
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

####################################

def calculate_sample_freqs(allele_counts_map, passed_sites_map, variant_type='4D', allowed_genes=None, fold=True):

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
        if fold == True:
            freqs = numpy.fmin(freqs,1-freqs) #fold
        for sample_idx in xrange(0,freqs.shape[1]):
            gene_freqs = freqs[:,sample_idx]
            sample_freqs[sample_idx].extend( gene_freqs[gene_freqs>0])
            
        passed_sites += numpy.diagonal(passed_sites_map[gene_name][variant_type]['sites'])
        
    
    return sample_freqs, passed_sites




####################################


def calculate_temporal_sample_freqs(allele_counts_map, passed_sites_map, initial_sample_idx, final_sample_idx, allowed_variant_types=set(['1D','2D','3D','4D']), allowed_genes=None):

    desired_samples = numpy.array([initial_sample_idx, final_sample_idx])
    
    initial_freqs = []
    final_freqs = []
    gene_names = []
    chromosomes = []
    positions = []
    marginal_initial_depths = []
    marginal_final_depths = []
    
    
    if allowed_genes == None:
        allowed_genes = set(passed_sites_map.keys())
 

    for gene_name in allowed_genes:
        for variant_type in allele_counts_map[gene_name].keys():
            
            if variant_type not in allowed_variant_types:
                continue
                
            allele_counts = allele_counts_map[gene_name][variant_type]['alleles']

            if len(allele_counts)==0:
                continue

            chunk_chromosomes = numpy.array([chromosome for chromosome, position in allele_counts_map[gene_name][variant_type]['locations']])
            
            chunk_positions = numpy.array([position for chromosome, position in allele_counts_map[gene_name][variant_type]['locations']])
            
            allele_counts = allele_counts[:,desired_samples,:]            
            depths = allele_counts.sum(axis=2)
            freqs = allele_counts[:,:,0]*1.0/(depths+(depths==0))
            
            initial_depths = depths[:,0]
            final_depths = depths[:,1]
            
            marginal_passed_sites = numpy.logical_or((initial_depths>0), (final_depths>0))
            joint_passed_sites = numpy.logical_and((initial_depths>0), (final_depths>0))
            
            initial_depths = initial_depths[marginal_passed_sites]
            final_depths = final_depths[marginal_passed_sites]
            
            unpolarized_initial_freqs = freqs[marginal_passed_sites,0]
            unpolarized_final_freqs = freqs[marginal_passed_sites,1]
        
            unpolarized_dfs = unpolarized_final_freqs-unpolarized_initial_freqs+normal(0,1e-06,unpolarized_initial_freqs.shape)
        
            polarized_initial_freqs = unpolarized_initial_freqs + (1-2*unpolarized_initial_freqs)*(unpolarized_dfs<=0)
            polarized_final_freqs = unpolarized_final_freqs + (1-2*unpolarized_final_freqs)*(unpolarized_dfs<=0)
        
        
            #polarized_initial_freqs = unpolarized_initial_freqs + (1-2*unpolarized_initial_freqs)*(unpolarized_initial_freqs>0.5)
            #polarized_final_freqs = unpolarized_final_freqs + (1-2*unpolarized_final_freqs)*(unpolarized_initial_freqs>0.5)
        
            initial_freqs.extend(polarized_initial_freqs)
            final_freqs.extend(polarized_final_freqs)
            
            marginal_initial_depths.extend(initial_depths)
            marginal_final_depths.extend(final_depths)
            
            gene_names.extend([gene_name]*len(polarized_final_freqs))
            chromosomes.extend(chunk_chromosomes[marginal_passed_sites])
            positions.extend(chunk_positions[marginal_passed_sites])
        
    print len(gene_names), len(chromosomes), len(initial_freqs), len(initial_depths)
            
    return numpy.array(gene_names), numpy.array(chromosomes), numpy.array(positions), numpy.array(initial_freqs), numpy.array(final_freqs), numpy.array(marginal_initial_depths), numpy.array(marginal_final_depths)

def calculate_triplet_sample_freqs(allele_counts_map, passed_sites_map, i, j, k, allowed_variant_types=set(['1D','2D','3D','4D']), allowed_genes=None):

    desired_samples = numpy.array([i, j, k])
    
    initial_freqs = []
    middle_freqs = []
    final_freqs = []
    gene_names = []
    chromosomes = []
    positions = []
    
    if allowed_genes == None:
        allowed_genes = set(passed_sites_map.keys())
 

    for gene_name in allowed_genes:
        for variant_type in allele_counts_map[gene_name].keys():
            
            if variant_type not in allowed_variant_types:
                continue
                
            allele_counts = allele_counts_map[gene_name][variant_type]['alleles']

            if len(allele_counts)==0:
                continue

            chunk_chromosomes = numpy.array([chromosome for chromosome, position in allele_counts_map[gene_name][variant_type]['locations']])
            
            chunk_positions = numpy.array([position for chromosome, position in allele_counts_map[gene_name][variant_type]['locations']])
            
            allele_counts = allele_counts[:,desired_samples,:]            
            depths = allele_counts.sum(axis=2)
            freqs = allele_counts[:,:,0]*1.0/(depths+(depths==0))
        
            initial_depths = depths[:,0]
            middle_depths = depths[:,1]
            final_depths = depths[:,2]
            
            joint_passed_sites = (initial_depths>0)*(middle_depths>0)*(final_depths>0)
            
            unpolarized_initial_freqs = freqs[joint_passed_sites,0]
            unpolarized_middle_freqs = freqs[joint_passed_sites,1]
            unpolarized_final_freqs = freqs[joint_passed_sites,2]
            
            # polarize sites
            flipped_sites = (unpolarized_initial_freqs>0.5)
            
            polarized_initial_freqs = unpolarized_initial_freqs + (1-2*unpolarized_initial_freqs)*(flipped_sites)
            polarized_middle_freqs = unpolarized_middle_freqs + (1-2*unpolarized_middle_freqs)*(flipped_sites)
            polarized_final_freqs = unpolarized_final_freqs + (1-2*unpolarized_final_freqs)*(flipped_sites)
        
            # changed sites
            changed_sites = (polarized_initial_freqs<=0.2)*numpy.logical_or(polarized_final_freqs>=0.8, polarized_middle_freqs>=0.8)         
            
            if changed_sites.sum() > 0:
            
                initial_freqs.extend(polarized_initial_freqs[changed_sites])
                middle_freqs.extend(polarized_middle_freqs[changed_sites])
                final_freqs.extend(polarized_final_freqs[changed_sites])
                gene_names.extend( long(changed_sites.sum()) * [gene_name] )
            
    
    freqs = []
    for f0,f1,f2 in zip(initial_freqs, middle_freqs, final_freqs):
        freqs.append((f0,f1,f2))
        
    return freqs
    

####################

def calculate_sample_freqs_2D(allele_counts_map, passed_sites_map, desired_samples, variant_type='4D', allowed_genes=None, fold=True):

    
    if allowed_genes == None:
        allowed_genes = set(passed_sites_map.keys())
     
    num_samples=sum(desired_samples)
    sample_freqs = [[] for i in xrange(0, num_samples)]
    joint_passed_sites= [[] for i in xrange(0, num_samples)]
    passed_sites = numpy.zeros((num_samples, num_samples))*1.0
    

    for gene_name in allowed_genes:

        allele_counts = allele_counts_map[gene_name][variant_type]['alleles']

        if len(allele_counts)==0:
            continue

        allele_counts = allele_counts[:,desired_samples,:]            
        depths = allele_counts.sum(axis=2)
        freqs = allele_counts[:,:,0]*1.0/(depths+(depths==0))
        joint_passed_sites_tmp=(depths>0)[:,None,:]*(depths>0)[:,:,None]

        if fold== True:
            freqs = numpy.fmin(freqs,1-freqs) 
        
        for sample_idx in xrange(0,freqs.shape[1]):
            gene_freqs = freqs[:,sample_idx]
            sample_freqs[sample_idx].extend(gene_freqs)
            joint_passed_sites[sample_idx].extend(joint_passed_sites_tmp[:,0,sample_idx])
            idx=numpy.where(desired_samples==True)
        passed_sites += passed_sites_map[gene_name][variant_type]['sites'][:,idx[0]][idx[0],:]
    
    return sample_freqs, passed_sites, joint_passed_sites

####################

        
def calculate_pooled_freqs(allele_counts_map, passed_sites_map, allowed_sample_idxs=[], allowed_variant_types = set(['1D','2D','3D','4D']), allowed_genes=set([]), lower_threshold=0.2,upper_threshold=0.8):

    if len(allowed_sample_idxs)==0:
        # all samples are allowed
        allowed_sample_idxs = numpy.array([True for i in xrange(0,allele_counts_map.values()[0].values()[0]['alleles'].shape[1])])

    if len(allowed_genes)==0:
        allowed_genes = set(passed_sites_map.keys())
    allowed_genes = allowed_genes & set(passed_sites_map.keys())
     
    pooled_freqs = []
    
    for gene_name in allowed_genes:
        
        for variant_type in allele_counts_map[gene_name].keys():
            
            if variant_type not in allowed_variant_types:
                continue
            
                
            allele_counts = allele_counts_map[gene_name][variant_type]['alleles']
        
            if len(allele_counts)==0:
                continue
            
            #print allele_counts_map[gene_name][variant_type]['alleles'].shape, allowed_sample_idxs.shape
                
            allele_counts = allele_counts[:,allowed_sample_idxs,:]
            
            genotype_matrix, passed_sites_matrix = calculate_consensus_genotypes(allele_counts,lower_threshold,upper_threshold)
            prevalences = (genotype_matrix*passed_sites_matrix).sum(axis=1)
            min_prevalences = 0.5
            max_prevalences = (passed_sites_matrix).sum(axis=1)-0.5
    
            polymorphic_sites = (prevalences>min_prevalences)*(prevalences<max_prevalences)
    
            gene_pooled_freqs = prevalences*1.0/(passed_sites_matrix).sum(axis=1)
            gene_pooled_freqs = gene_pooled_freqs[polymorphic_sites]
            gene_pooled_freqs = numpy.fmin(gene_pooled_freqs,1-gene_pooled_freqs)
            pooled_freqs.extend(gene_pooled_freqs)

    pooled_freqs = numpy.array(pooled_freqs)
    return pooled_freqs

def calculate_pooled_counts(allele_counts_map, passed_sites_map, allowed_sample_idxs=[], allowed_variant_types = set(['1D','2D','3D','4D']), allowed_genes=set([]),pi_min_k=1,lower_threshold=0.2,upper_threshold=0.8):

    if len(allowed_sample_idxs)==0:
        # all samples are allowed
        allowed_sample_idxs = numpy.array([True for i in xrange(0,allele_counts_map.values()[0].values()[0]['alleles'].shape[1])])

    if len(allowed_genes)==0:
        allowed_genes = set(passed_sites_map.keys())
    allowed_genes = allowed_genes & set(passed_sites_map.keys())
     
    pi_weighted_number = 0
    pooled_counts = []
    
    for gene_name in allowed_genes:
        
        for variant_type in allele_counts_map[gene_name].keys():
            
            if variant_type not in allowed_variant_types:
                continue
                
            allele_counts = allele_counts_map[gene_name][variant_type]['alleles']
        
            if len(allele_counts)==0:
                continue
            
            #print allele_counts_map[gene_name][variant_type]['alleles'].shape, allowed_sample_idxs.shape
                
            allele_counts = allele_counts[:,allowed_sample_idxs,:]
            
            genotype_matrix, passed_sites_matrix = calculate_consensus_genotypes(allele_counts,lower_threshold,upper_threshold)
            prevalences = (genotype_matrix*passed_sites_matrix).sum(axis=1)
            min_prevalences = 0.5
            max_prevalences = (passed_sites_matrix).sum(axis=1)-0.5
    
            polymorphic_sites = (prevalences>min_prevalences)*(prevalences<max_prevalences)
    
            ks = prevalences[polymorphic_sites]
            ns = passed_sites_matrix.sum(axis=1)[polymorphic_sites]
            minor_ks = numpy.fmin(ks,ns-ks)
            pooled_counts.extend( minor_ks )
            
            pi_weighted_number += (ks*(ns-ks)*2.0/(ns*(ns-1))*(minor_ks>=pi_min_k)).sum()
            
    pooled_counts = numpy.array(pooled_counts)
    return pooled_counts, pi_weighted_number

def calculate_private_snvs(samples, allele_counts_map, passed_sites_map, allowed_variant_types=set([]), allowed_genes=set([]), lower_threshold=config.consensus_lower_threshold, 
upper_threshold=config.consensus_upper_threshold):

    # First 
    sample_host_matrix, hosts = sample_utils.calculate_sample_subject_matrix(samples)
    
    total_genes = set(passed_sites_map.keys())

    if len(allowed_genes)==0:
        allowed_genes = set(passed_sites_map.keys())
    
    allowed_genes = (allowed_genes & total_genes)     
    
    if len(allowed_variant_types)==0:
        allowed_variant_types = set(['1D','2D','3D','4D'])    
         
    private_snvs = []
     
    for gene_name in allowed_genes:
        
        for variant_type in passed_sites_map[gene_name].keys():
             
            if variant_type not in allowed_variant_types:
                continue
            
            passed_sites = passed_sites_map[gene_name][variant_type]['sites']
   
            allele_counts = allele_counts_map[gene_name][variant_type]['alleles']                        
            if len(allele_counts)==0:
                continue
                
            locations = allele_counts_map[gene_name][variant_type]['locations']
            
            depths = allele_counts.sum(axis=2)
            freqs = allele_counts[:,:,0]*1.0/(depths+(depths==0))
            
            derived_sites = (freqs>=upper_threshold)
            ancestral_sites = (freqs<=lower_threshold)
            
            # Sites where the major allele is at sufficiently high frequency
            high_freq_sites = numpy.logical_or(ancestral_sites, derived_sites)
            
            passed_depths = (depths>0)
            
            confident_sites = numpy.logical_and(high_freq_sites, passed_depths)
            
            #print confident_sites.shape
            
            # goes from L x n to L x h (sites across all hosts)
            host_confident_sites = (numpy.einsum('ij,jk', confident_sites, sample_host_matrix)>0.5)
            
            host_derived_sites = (numpy.einsum('ij,jk',derived_sites, sample_host_matrix)>0.5)
            
            #print host_confident_sites.shape
            
            host_sample_sizes = host_confident_sites.sum(axis=1)
            host_derived_counts = host_derived_sites.sum(axis=1)
            
            #print host_sample_sizes.shape
            #print host_derived_counts.shape
            
            private_idxs = numpy.nonzero((host_sample_sizes>3.5)*(host_derived_counts==1))[0]
            
            #print private_idxs.shape
            
            if len(private_idxs)>0:
                for snp_idx in private_idxs: 
                    host = hosts[numpy.nonzero(host_derived_sites[snp_idx])[0][0]]
                    
                    contig, location = allele_counts_map[gene_name][variant_type]['locations'][snp_idx]
                
                    private_snvs.append((contig, location, gene_name, variant_type, host))
            
    return private_snvs
     

def calculate_singleton_matrix(allele_counts_map, passed_sites_map, allowed_variant_types=set([]), allowed_genes=set([]), lower_threshold=config.consensus_lower_threshold, 
upper_threshold=config.consensus_upper_threshold):

    total_genes = set(passed_sites_map.keys())

    if len(allowed_genes)==0:
        allowed_genes = set(passed_sites_map.keys())
    
    allowed_genes = (allowed_genes & total_genes)     
    
    if len(allowed_variant_types)==0:
        allowed_variant_types = set(['1D','2D','3D','4D'])    
         
    doubleton_matrix = numpy.zeros_like( passed_sites_map.values()[0].values()[0]['sites'] )*1.0   
    singleton_matrix = numpy.zeros_like(doubleton_matrix)
    
    difference_matrix = numpy.zeros_like(doubleton_matrix)
    
    opportunity_matrix = numpy.zeros_like(doubleton_matrix)
    
     
    for gene_name in allowed_genes:
        
        for variant_type in passed_sites_map[gene_name].keys():
             
            if variant_type not in allowed_variant_types:
                continue
            
            passed_sites = passed_sites_map[gene_name][variant_type]['sites']
   
            allele_counts = allele_counts_map[gene_name][variant_type]['alleles']                        
            if len(allele_counts)==0:
                continue
            
            depths = allele_counts.sum(axis=2)
            freqs = allele_counts[:,:,0]*1.0/(depths+(depths==0))
            
            derived_sites = (freqs>=upper_threshold)
            ancestral_sites = (freqs<=lower_threshold)
            
            # Sites where the major allele is at sufficiently high frequency
            high_freq_sites = numpy.logical_or(ancestral_sites, derived_sites)
            # Those where it is not
            intermediate_freq_sites = numpy.logical_not(high_freq_sites)
            
            # site*sample*sample matrix of sites with sufficient coverage in both samples
            passed_depths = (depths>0)[:,:,None]*(depths>0)[:,None,:]

            # site*sample*sample matrix of sites where we can look for differences            
            confident_sites = numpy.logical_and(high_freq_sites[:,:,None], high_freq_sites[:,None,:])*passed_depths
            
            site_difference_matrix = numpy.logical_or(derived_sites[:,:,None]*ancestral_sites[:,None,:],  ancestral_sites[:,:,None]*derived_sites[:,None,:] )
            
            # this is really the only place you have to switch 
            # to within-between hosts
            
            site_sample_size_matrix = confident_sites.sum(axis=2)
            # Want at least a sample size of 4
            site_sufficient_sample_size_matrix = (site_sample_size_matrix>3.5)
            # 
            site_total_difference_matrix = (site_difference_matrix*confident_sites).sum(axis=2)
            
            potential_singletons = ((site_sample_size_matrix-site_total_difference_matrix)==1)
            
            potential_doubletons = ((site_sample_size_matrix-site_total_difference_matrix)==2)
            
            
            # Confident sites is now asymmetric
            confident_sites *= site_sufficient_sample_size_matrix[:,:,None]
            
            # Sites that we can't count (but we would have counted in passed_sites)
            non_confident_sites = numpy.logical_not(confident_sites)*passed_depths
            
            
            # total number of differences between i and j 
            # (regardless of singleton/doubleton status)
            differences = (site_difference_matrix*confident_sites).sum(axis=0)
            
            singletons = (confident_sites*site_difference_matrix*potential_singletons[:,:,None]).sum(axis=0)
            
            doubletons = (confident_sites * numpy.logical_not(site_difference_matrix) * potential_doubletons[:,:,None]).sum(axis=0)
            
            opportunities = passed_sites - (numpy.logical_not(confident_sites)*passed_depths).sum(axis=0) 
            
            singleton_matrix += singletons
            doubleton_matrix += doubletons
            difference_matrix += differences
            opportunity_matrix += opportunities
            
    
    return doubleton_matrix, singleton_matrix, difference_matrix, opportunity_matrix
    
def calculate_mutation_reversion_matrix(allele_counts_map, passed_sites_map, allowed_variant_types=set([]), allowed_genes=set([]), lower_threshold=config.consensus_lower_threshold, 
upper_threshold=config.consensus_upper_threshold, min_change=config.fixation_min_change):

    total_genes = set(passed_sites_map.keys())

    if len(allowed_genes)==0:
        allowed_genes = set(passed_sites_map.keys())
    
    allowed_genes = (allowed_genes & total_genes)     
    
    if len(allowed_variant_types)==0:
        allowed_variant_types = set(['1D','2D','3D','4D'])    
         
    mut_fixation_matrix = numpy.zeros_like( passed_sites_map.values()[0].values()[0]['sites'] )*1.0   
    rev_fixation_matrix = numpy.zeros_like(mut_fixation_matrix)
    
    mut_opportunity_matrix = numpy.zeros_like(mut_fixation_matrix)
    rev_opportunity_matrix = numpy.zeros_like(rev_fixation_matrix)
              
    for gene_name in allowed_genes:
        
        for variant_type in passed_sites_map[gene_name].keys():
             
            if variant_type not in allowed_variant_types:
                continue
        
        
            passed_sites = passed_sites_map[gene_name][variant_type]['sites']
   
            allele_counts = allele_counts_map[gene_name][variant_type]['alleles']                        
            if len(allele_counts)==0:
                continue
            
            depths = allele_counts.sum(axis=2)
            freqs = allele_counts[:,:,0]*1.0/(depths+(depths==0))
            
            derived_sites = (freqs>=upper_threshold)
            ancestral_sites = (freqs<=lower_threshold)
            
            # Sites where the major allele is at sufficiently high frequency
            high_freq_sites = numpy.logical_or(ancestral_sites, derived_sites)
            # Those where it is not
            intermediate_freq_sites = numpy.logical_not(high_freq_sites)
            
            # site*sample*sample matrix of sites with sufficient coverage in both samples
            passed_depths = (depths>0)[:,:,None]*(depths>0)[:,None,:]

            # site*sample*sample matrix of sites where we can look for differences            
            confident_sites = numpy.logical_and(high_freq_sites[:,:,None], high_freq_sites[:,None,:])*passed_depths
            

            # site*sample*sample matrix of sites that are missing data 
            # based on allele freqs, but which had sufficient coverage
            # (we need to remove these from opportunities below)
            missing_data_sites = numpy.logical_or(intermediate_freq_sites[:,:,None],intermediate_freq_sites[:,None,:])*passed_depths
              
            # Calculate mutations and reversions
            mutations = (ancestral_sites[:,:,None])*(derived_sites[:,None,:])*confident_sites
            
            reversions = (derived_sites[:,:,None])*(ancestral_sites[:,None,:])*confident_sites
            
            # sites were you could have had a reversion
            reversion_opportunities = derived_sites[:,:,None]*confident_sites
            
            mut_fixation_matrix += (mutations).sum(axis=0)
            rev_fixation_matrix += (reversions).sum(axis=0)
            
            rev_opportunity_matrix += (reversion_opportunities).sum(axis=0)
            mut_opportunity_matrix += (passed_sites - missing_data_sites.sum(axis=0) - reversion_opportunities.sum(axis=0) ) 
            
    return mut_fixation_matrix, rev_fixation_matrix, mut_opportunity_matrix, rev_opportunity_matrix    
    

def calculate_fixation_matrix(allele_counts_map, passed_sites_map, allowed_variant_types=set([]), allowed_genes=set([]), lower_threshold=config.consensus_lower_threshold, 
upper_threshold=config.consensus_upper_threshold, min_change=config.fixation_min_change):
    
    mut_fixation_matrix, rev_fixation_matrix, mut_opportunity_matrix, rev_opportunity_matrix = calculate_mutation_reversion_matrix(allele_counts_map, passed_sites_map, allowed_variant_types, allowed_genes, lower_threshold, 
upper_threshold, min_change)

    fixation_matrix = mut_fixation_matrix + rev_fixation_matrix
    opportunity_matrix = mut_opportunity_matrix + rev_opportunity_matrix
    
    return fixation_matrix, opportunity_matrix
    
    
# same as above, but returns two matrices with counts of 
# mutations (i->j away from consensus allele) and
# reversion (i->j toward consensus allele)
def calculate_fixation_matrix_mutation_reversion(allele_counts_map, passed_sites_map, allowed_variant_types=set([]), allowed_genes=set([]), lower_threshold=config.consensus_lower_threshold, 
upper_threshold=config.consensus_upper_threshold, min_change=config.fixation_min_change):

    total_genes = set(passed_sites_map.keys())

    if len(allowed_genes)==0:
        allowed_genes = set(passed_sites_map.keys())
    
    allowed_genes = (allowed_genes & total_genes)     
    
    if len(allowed_variant_types)==0:
        allowed_variant_types = set(['1D','2D','3D','4D'])    
                    
    fixation_matrix_mutation = numpy.zeros_like(passed_sites_map.values()[0].values()[0]['sites'])*1.0 
    fixation_matrix_reversion = numpy.zeros_like(fixation_matrix_mutation)*1.0
     
    passed_sites = numpy.zeros_like(fixation_matrix_mutation)*1.0
    
    for gene_name in allowed_genes:
        
        for variant_type in passed_sites_map[gene_name].keys():
             
            if variant_type not in allowed_variant_types:
                continue
        
            passed_sites += passed_sites_map[gene_name][variant_type]['sites']
   
            allele_counts = allele_counts_map[gene_name][variant_type]['alleles']                        
            if len(allele_counts)==0:
                continue
            
            depths = allele_counts.sum(axis=2)
            freqs = allele_counts[:,:,0]*1.0/(depths+(depths==0))
            
            intermediate_freq_sites = (freqs>lower_threshold)*(freqs<upper_threshold)
   
            passed_depths = (depths>0)[:,:,None]*(depths>0)[:,None,:]
            
            bad_sites = numpy.logical_or(intermediate_freq_sites[:,:,None],intermediate_freq_sites[:,None,:])*passed_depths
            
            delta_freqs = (freqs[:,:,None]-freqs[:,None,:])*passed_depths
            
            mutations = (delta_freqs>=min_change)
            reversions = (delta_freqs<=(-1*min_change))
            
            fixation_matrix_mutation += mutations.sum(axis=0) # sum over sites
            fixation_matrix_reversion += reversions.sum(axis=0) # sum over sites
            
            passed_sites -= bad_sites.sum(axis=0)
            
    return fixation_matrix_mutation, fixation_matrix_reversion, passed_sites  

# same as above, but returns two matrices with counts of 
# mutations (i->j away from consensus allele) and
# reversion (i->j toward consensus allele)
def calculate_preexisting_snps(allele_counts_map, passed_sites_map, allowed_variant_types=set([]), allowed_genes=set([]),min_freq=0.1):

    total_genes = set(passed_sites_map.keys())

    if len(allowed_genes)==0:
        allowed_genes = set(passed_sites_map.keys())
    
    allowed_genes = (allowed_genes & total_genes)     
    
    if len(allowed_variant_types)==0:
        allowed_variant_types = set(['1D','2D','3D','4D'])    
     
    snp_locations = [] 
                    
    for gene_name in allowed_genes:
        
        for variant_type in passed_sites_map[gene_name].keys():
             
            if variant_type not in allowed_variant_types:
                continue
            
            allele_counts = allele_counts_map[gene_name][variant_type]['alleles']                        
            if len(allele_counts)==0:
                continue
            
            depths = allele_counts.sum(axis=2)
            freqs = allele_counts[:,:,0]*1.0/(depths+(depths==0))
            
            site_raw_prevalence = (depths>0).sum(axis=1)
            snp_raw_prevalence = (freqs>min_freq).sum(axis=1)
            
            snp_prevalence = snp_raw_prevalence*1.0/(site_raw_prevalence+(site_raw_prevalence==0))            
            
            polymorphic_sites = (snp_raw_prevalence>0)
            if polymorphic_sites.sum()==0:
                continue
                
            # get locations
            polymorphic_site_idxs = numpy.nonzero(polymorphic_sites)[0]
            for idx in polymorphic_site_idxs:
                snp_locations.append( (allele_counts_map[gene_name][variant_type]['locations'][idx][0], allele_counts_map[gene_name][variant_type]['locations'][idx][1], snp_prevalence[idx] ) )
            
    return snp_locations

####
#
# Calculates the number of within-patient polymorphism differences between
# two samples. (e.g. something that is fixed in one timepoint and polymorphic
# in another. 
#
####
def calculate_new_snp_matrix(allele_counts_map, passed_sites_map, allowed_variant_types=set([]), allowed_genes=set([]), min_freq=0.05, max_freq=0.2):

    total_genes = set(passed_sites_map.keys())

    if len(allowed_genes)==0:
        allowed_genes = set(passed_sites_map.keys())
    
    allowed_genes = (allowed_genes & total_genes)     
    
    if len(allowed_variant_types)==0:
        allowed_variant_types = set(['1D','2D','3D','4D'])    
                    
    new_snp_matrix = numpy.zeros_like(passed_sites_map.values()[0].values()[0]['sites'])*1.0  
    passed_sites = numpy.zeros_like(new_snp_matrix)*1.0
    
    for gene_name in allowed_genes:
        
        for variant_type in passed_sites_map[gene_name].keys():
             
            if variant_type not in allowed_variant_types:
                continue
        
            passed_sites += passed_sites_map[gene_name][variant_type]['sites']
   
            allele_counts = allele_counts_map[gene_name][variant_type]['alleles']                        
            if len(allele_counts)==0:
                continue
            

            depths = allele_counts.sum(axis=2)
            freqs = allele_counts[:,:,0]/(depths+(depths==0))
            # turn into minor allele frequencies
            mafs = numpy.fmin(freqs,1-freqs)
            
            # Turn
            
            new_snps_1 = (mafs[:,:,None]<min_freq)*(mafs[:,None,:]>max_freq)
            new_snps_2 = (mafs[:,:,None]>max_freq)*(mafs[:,None,:]<min_freq)
            total_new_snps = new_snps_1+new_snps_2
             
            passed_depths = (depths>0)[:,:,None]*(depths>0)[:,None,:]
    
            total_new_snps[passed_depths==0] = 0
            
            new_snp_matrix += total_new_snps.sum(axis=0)
        
    return new_snp_matrix, passed_sites  


   
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
     
    # We used to normalize here    
    #pi_matrix = pi_matrix /(passed_sites+(passed_sites==0))
    #avg_pi_matrix = avg_pi_matrix/(passed_sites+(passed_sites==0))
    # Now we return passed sites
    
    return pi_matrix, avg_pi_matrix, passed_sites



    
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
# Returns list of differences, number of comparisons. Each difference is a tuple of form 
#
# (gene_name, (contig, location), (alt_i, depth_i), (alt_j, depth_j))
#
def calculate_snp_differences_between(i,j,allele_counts_map, passed_sites_map, avg_depth_i, avg_depth_j, allowed_variant_types=set([]), allowed_genes=set([]), lower_threshold=config.consensus_lower_threshold, 
upper_threshold=config.consensus_upper_threshold, log10_depth_ratio_threshold=config.fixation_log10_depth_ratio_threshold):

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
    
            safe_depths = depths+(depths==0)
            
            log10_depth_ratios = numpy.fabs(numpy.log10((safe_depths[:,0]/avg_depth_i)/(safe_depths[:,1]/avg_depth_j)))
                
            passed_depths = (depths>0)[:,0]*(depths>0)[:,1]*(log10_depth_ratios<log10_depth_ratio_threshold)

            mutations = (alt_freqs[:,0]<=lower_threshold)*(alt_freqs[:,1]>=upper_threshold)*passed_depths
            reversions = (alt_freqs[:,0]>=upper_threshold)*(alt_freqs[:,1]<=lower_threshold)*passed_depths
            
            changed_sites = numpy.nonzero( numpy.logical_or(mutations, reversions) )[0]
            
            if len(changed_sites)>0:
                # some fixations!
                
                for idx in changed_sites:
                    snp_changes.append((gene_name, allele_counts_map[gene_name][variant_type]['locations'][idx], variant_type, (allele_counts[idx,0,0], depths[idx,0]), (allele_counts[idx,1,0],depths[idx,1]) ))
                        
    return snp_changes

# Calculate polarized SNP changes from i to j that exceed threshold 
# Returns list of differences, number of comparisons. Each difference is a tuple of form 
#
# (gene_name, (contig, location), (alt_i, depth_i), (alt_j, depth_j))
#
def calculate_tracked_private_snvs(i,j,allele_counts_map, passed_sites_map, avg_depth_i, avg_depth_j, private_snv_map, allowed_variant_types=set([]), allowed_genes=set([]), lower_threshold=config.consensus_lower_threshold, 
upper_threshold=config.consensus_upper_threshold, log10_depth_ratio_threshold=config.fixation_log10_depth_ratio_threshold):

    if len(allowed_genes)==0:
        allowed_genes = set(passed_sites_map.keys())
        
    if len(allowed_variant_types)==0:
        allowed_variant_types = set(['1D','2D','3D','4D'])    
    
    tracked_private_snps = []    
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
    
            safe_depths = depths+(depths==0)
            
            log10_depth_ratios = numpy.fabs(numpy.log10((safe_depths[:,0]/avg_depth_i)/(safe_depths[:,1]/avg_depth_j)))
                
            passed_depths = (depths>0)[:,0]*(depths>0)[:,1]*(log10_depth_ratios<log10_depth_ratio_threshold)

            initial_high_freqs = alt_freqs[:,0]>=upper_threshold
            final_high_freqs = alt_freqs[:,1]>=upper_threshold
            final_low_freqs = alt_freqs[:,1]<=lower_threshold
            
            potential_private_snps = numpy.nonzero( initial_high_freqs*numpy.logical_or(final_high_freqs, final_low_freqs) )[0]
            
            if len(potential_private_snps)>0:
                # some candidates for private SNVs
                for idx in potential_private_snps:
                    # check to see if it is indeed a private SNV
                    location_tuple = allele_counts_map[gene_name][variant_type]['locations'][idx]
                    
                    if location_tuple in private_snv_map:
                        # it is indeed private! 
                        tracked_private_snps.append((gene_name, location_tuple, variant_type, (allele_counts[idx,0,0], depths[idx,0]), (allele_counts[idx,1,0],depths[idx,1]) ))
            
    return tracked_private_snps
    
    
    
#########################################

def calculate_mean_pi_matrix_per_pathway(pi_per_gene, avg_pi_per_gene, passed_sites_per_gene,num_people_with_data, kegg_ids,min_passed_sites_per_person=100):
    
    pi_per_pathway={}
    avg_pi_per_pathway={}
    passed_sites_per_pathway={}
    num_genes_per_pathway={}
    num_people_with_data_pathway={}
    gene_name=avg_pi_per_gene.keys()[0]
    
    pi_per_pathway['Annotated pathways'] = numpy.zeros_like(pi_per_gene[gene_name])
    avg_pi_per_pathway['Annotated pathways']=numpy.zeros_like(avg_pi_per_gene[gene_name])
    passed_sites_per_pathway['Annotated pathways']=numpy.zeros_like(passed_sites_per_gene[gene_name])
    num_genes_per_pathway['Annotated pathways']=0
    num_people_with_data_pathway['Annotated pathways']=0
    
    for gene_name in avg_pi_per_gene.keys():
        pathway=kegg_ids[gene_name][0][1]
        if pathway not in avg_pi_per_pathway.keys():
            pi_per_pathway[pathway]=pi_per_gene[gene_name]
            avg_pi_per_pathway[pathway]=avg_pi_per_gene[gene_name]
            passed_sites_per_pathway[pathway]=passed_sites_per_gene[gene_name]
            num_genes_per_pathway[pathway]=1
            num_people_with_data_pathway[pathway]=num_people_with_data[gene_name]
        else:
            pi_per_pathway[pathway]+=pi_per_gene[gene_name]
            avg_pi_per_pathway[pathway]+=avg_pi_per_gene[gene_name]
            passed_sites_per_pathway[pathway]+=passed_sites_per_gene[gene_name]  
            num_genes_per_pathway[pathway]+=1
            num_people_with_data_pathway[pathway]+=num_people_with_data[gene_name]       
        if pathway !='':
            pi_per_pathway['Annotated pathways'] += pi_per_gene[gene_name]
            avg_pi_per_pathway['Annotated pathways'] +=avg_pi_per_gene[gene_name]
            passed_sites_per_pathway['Annotated pathways'] +=passed_sites_per_gene[gene_name]
            num_genes_per_pathway['Annotated pathways']+=1
            num_people_with_data_pathway['Annotated pathways']+=num_people_with_data[gene_name]
            
    for pathway_name in avg_pi_per_pathway.keys():
        # we want to identify people that have few passed sites even after aggregating the data accross genes. Then set the values in these cells to zero because these data points are too noisy
        low_passed_sites_idxs=passed_sites_per_pathway[pathway_name]<min_passed_sites_per_person
        passed_sites_per_pathway[pathway_name][low_passed_sites_idxs]=0
        avg_pi_per_pathway[pathway_name][low_passed_sites_idxs]=0
        pi_per_pathway[pathway_name][low_passed_sites_idxs]=0
        # now compute pi/pathway.  
        avg_pi_per_pathway[pathway_name] = avg_pi_per_pathway[pathway_name]/(passed_sites_per_pathway[pathway_name]+(passed_sites_per_pathway[pathway_name]==0))     
        pi_per_pathway[pathway_name] = pi_per_pathway[pathway_name]/(passed_sites_per_pathway[pathway_name]+(passed_sites_per_pathway[pathway_name]==0))     
        #num_people_with_data_pathway[pathway_name]=sum(numpy.diagonal(passed_sites_per_pathway[pathway_name])>=min_passed_sites_per_person)
        num_people_with_data_pathway[pathway_name]= num_people_with_data_pathway[pathway_name]/num_genes_per_pathway[pathway_name]
    return pi_per_pathway,avg_pi_per_pathway,passed_sites_per_pathway,num_people_with_data_pathway, num_genes_per_pathway

#################################


def calculate_mean_fixation_matrix_per_pathway(fixation_per_gene, passed_sites_per_gene,num_people_with_data, kegg_ids, min_passed_sites_per_person=100):
    
    fixation_per_pathway={}
    passed_sites_per_pathway={}
    num_genes_per_pathway={}
    num_people_with_data_pathway={}

    gene_name=fixation_per_gene.keys()[0]
    fixation_per_pathway['Annotated pathways'] = numpy.zeros_like(fixation_per_gene[gene_name])
    passed_sites_per_pathway['Annotated pathways']=numpy.zeros_like(passed_sites_per_gene[gene_name])
    num_genes_per_pathway['Annotated pathways']=0
    num_people_with_data_pathway['Annotated pathways']=0


    for gene_name in fixation_per_gene.keys():
        pathway=kegg_ids[gene_name][0][1]
        if pathway not in fixation_per_pathway.keys():
            fixation_per_pathway[pathway]=fixation_per_gene[gene_name]
            passed_sites_per_pathway[pathway]=passed_sites_per_gene[gene_name]
            num_genes_per_pathway[pathway]=1
            num_people_with_data_pathway[pathway]=num_people_with_data[gene_name]
        else:
            fixation_per_pathway[pathway]+=fixation_per_gene[gene_name]
            passed_sites_per_pathway[pathway]+=passed_sites_per_gene[gene_name]  
            num_genes_per_pathway[pathway]+=1
            num_people_with_data_pathway[pathway]+=num_people_with_data[gene_name]
        if pathway !='':
            fixation_per_pathway['Annotated pathways'] += fixation_per_gene[gene_name]
            passed_sites_per_pathway['Annotated pathways'] +=passed_sites_per_gene[gene_name]
            num_genes_per_pathway['Annotated pathways']+=1
            num_people_with_data_pathway['Annotated pathways']+=num_people_with_data[gene_name]
            
            
    for pathway_name in fixation_per_pathway.keys():
       # we want to identify people that have few passed sites even after aggregating the data accross genes. Then set the values in these cells to zero because these data points are too noisy
        low_passed_sites_idxs=passed_sites_per_pathway[pathway_name]<min_passed_sites_per_person
        passed_sites_per_pathway[pathway_name][low_passed_sites_idxs]=0
        fixation_per_pathway[pathway_name][low_passed_sites_idxs]=0
        #now compute fixation/pathway
        fixation_per_pathway[pathway_name] = fixation_per_pathway[pathway_name]/(passed_sites_per_pathway[pathway_name]+(passed_sites_per_pathway[pathway_name]==0))     
        num_people_with_data_pathway[pathway_name]=num_people_with_data_pathway[pathway_name]/float(num_genes_per_pathway[pathway_name])
    return fixation_per_pathway, passed_sites_per_pathway, num_people_with_data_pathway, num_genes_per_pathway

#######################
#
# Calculate pi from SFS map
#
#######################
def calculate_pi_from_sfs_map(sfs_map):
    
    alts = []
    refs = []
    depths = []
    counts = []
    for key in sfs_map.keys():
        D,A = key
        n = sfs_map[key][0]
        
        alts.append(A)
        refs.append(D-A)
        depths.append(D)
        counts.append(n)
    
    alts = numpy.array(alts)
    refs = numpy.array(refs)
    depths = numpy.array(depths)
    counts = numpy.array(counts)
    
    alt_lower_threshold = numpy.ceil(depths*0.05)+0.5 #at least one read above 5%.
    alts[alts<alt_lower_threshold] = 0
    alt_upper_threshold = numpy.floor(depths*0.95)-0.5 #at least one read below 95%
    alts[alts>alt_upper_threshold] = depths[alts>alt_upper_threshold]
        
    total_pi = ((2*alts*(depths-alts)*1.0/(depths*(depths-1)+(depths<1.1)))*(counts)).sum()
    num_opportunities = counts.sum()
    
    return total_pi/num_opportunities
    
def calculate_polymorphism_rates_from_sfs_map(sfs_map,lower_threshold=0.2,upper_threshold=0.8):
    
    total_sites = 0
    within_sites = 0
    between_sites = 0
    for key in sfs_map.keys():
        D,A = key
        n = sfs_map[key][0]
        reverse_n = sfs_map[key][1]
        
        f = A*1.0/D
        
        total_sites += n
        
        if ((f>lower_threshold) and (f<upper_threshold)):
            # an intermediate frequency site
            within_sites += n
        else:    
            if f>0.5:
                between_sites += (n-reverse_n)
            else:
                between_sites += reverse_n
        
        
    between_polymorphism_rate = between_sites*1.0/total_sites
    within_polymorphism_rate = within_sites*1.0/total_sites
    
    return within_polymorphism_rate, between_polymorphism_rate
    
#######################
#
# Estimate smoothed within-person SFS with EM algorithm
#
#######################
def calculate_smoothed_sfs(sfs_map, num_iterations=100, perr=0.01, lower_threshold=config.consensus_lower_threshold, upper_threshold=config.consensus_upper_threshold):
    
    alts = []
    refs = []
    depths = []
    counts = []
    for key in sfs_map.keys():
        D,A = key
        n = sfs_map[key][0]
        
        alts.append(A)
        refs.append(D-A)
        depths.append(D)
        counts.append(n)
    
    alts = numpy.array(alts)
    refs = numpy.array(refs)
    depths = numpy.array(depths)
    counts = numpy.array(counts)
    weights = counts*1.0/counts.sum()
    
    # calculate median depth (or rough approximation)
    sorted_depths, sorted_counts = (numpy.array(x) for x in zip(*sorted(zip(depths, counts))))
    CDF = numpy.cumsum(sorted_counts)*1.0/sorted_counts.sum()
    Dbar = sorted_depths[CDF>0.5][0]
    #Dbar = min([Dbar,100])
    
    Abars = numpy.arange(0,Dbar+1)
    Rbars = Dbar-Abars
    fs = Abars*1.0/Dbar
    df = fs[1]-fs[0]
    flowers=  fs-df/2
    flowers[0] = 0-1e-10
    fuppers = fs+df/2
    fuppers[-1] = 1+1e-10
    
    pfs = numpy.zeros_like(fs)
    
    
    # first infer rate of polymorphisms (p_poly) using EM
    
    # Initial guess
    p_poly = 1e-04
    
    # calculate probability of data, conditioned on it not being polymorphic
    # (i.e., alt reads are sequencing errors)
    # (this doesn't depend on p_poly)
    pdata_errs = (betainc(alts+1,refs+1,perr)+betainc(refs+1,alts+1,perr))/(2*perr)
    pdata_intermediates = 1-(betainc(alts+1,refs+1, lower_threshold)+betainc(refs+1,alts+1,1-upper_threshold)) 
    # EM loop
    for iteration in xrange(0,num_iterations):
        posterior_polys = 1.0/(1.0+(1-p_poly)/(p_poly)*pdata_errs)
        p_poly = (posterior_polys*weights).sum()
    
    
    
    # Calculate avg posterior probability of freq being between lower and upper threshold
    p_intermediate = (posterior_polys*pdata_intermediates*weights).sum()
    
    # Now Calculate smoothed SFS estimate
    
    # Posterior method
    #posterior_frequencies = (betainc(alts[:,None]+1,refs[:,None]+1, fuppers[None,:])-betainc(alts[:,None]+1,refs[:,None]+1,flowers[None,:]))
    # The reason why we don't use this one is that it assumes a higher variance than our internal model. In reality, we believe that there are a few fixed frequencies, not that every one is independent. (Really we'd want to do some sort of EM, but it's slowly converging)
    
    
    
    # Bin overlap method
    
    freqs = alts*1.0/depths
    freqs_plushalf = numpy.clip((alts+0.5)*1.0/depths,0,1)
    freqs_minushalf = numpy.clip((alts-0.5)*1.0/depths,0,1)
    
    a = numpy.fmax(flowers[None,:],freqs_minushalf[:,None])
    b = numpy.fmin(fuppers[None,:],freqs_plushalf[:,None])
    
    posterior_frequencies = (b-a)*(b>a)/(freqs_plushalf-freqs_minushalf)[:,None]
    
    # Delta function method
    #posterior_frequencies = (freqs[:,None]>flowers[None,:])*(freqs[:,None]<=fuppers[None,:]) 
    # the reason why we don't use this one is that it suffers from binning artefacts 
    # though not *so* bad
    
    #pfs = ((posterior_frequencies)*((posterior_polys*weights)[:,None])).sum(axis=0)
    pfs = ((posterior_frequencies)*((weights)[:,None])).sum(axis=0)
    
    pfs /= pfs.sum()
    
    
    # Re-sampling method (too smooth)
    #prefactors = numpy.exp( loggamma(Abars[None,:]+alts[:,None]+1)+loggamma(Rbars[None,:]+refs[:,None]+1)+loggamma(Dbar+1)+loggamma(depths+1)[:,None]-loggamma(Dbar+depths+2)[:,None]-loggamma(Abars+1)[None,:]-loggamma(Rbars+1)[None,:]-loggamma(alts+1)[:,None]-loggamma(refs+1)[:,None])
    #pfs = ((prefactors*(p_poly+(1-p_poly)*(betainc(Abars[None,:]+alts[:,None]+1, Rbars[None,:]+refs[:,None]+1, perr)+betainc(Rbars[None,:]+refs[:,None]+1, Abars[None,:]+alts[:,None]+1, perr))/(2*perr)))*weights[:,None]).sum(axis=0)
    
    print p_poly, p_intermediate, Dbar
    
    return fs, pfs, p_intermediate, p_poly
    
    
    
#######################
#
# Estimate smoothed within-person SFS with EM algorithm
#
#######################
def calculate_smoothed_sfs_continuous_EM(sfs_map,fs=[],num_iterations=100):
    
    alts = []
    refs = []
    depths = []
    counts = []
    for key in sfs_map.keys():
        D,A = key
        n = sfs_map[key][0]
        
        alts.append(A)
        refs.append(D-A)
        depths.append(D)
        counts.append(n)
    
    alts = numpy.array(alts)
    refs = numpy.array(refs)
    depths = numpy.array(depths)
    counts = numpy.array(counts)
    
    
    weights = counts*1.0/counts.sum()
    
    if len(fs)==0:
        fs = numpy.linspace(0,1,101)[1:-1]
        
    dfs = fs[1]-fs[0]
    
    logfs = numpy.log(fs)
    log1minusfs = numpy.log(1-fs)
    
    # initial guess for pfs    
    pfs = numpy.zeros_like(fs)
    pfs[fs>=0.99] = 1e-02/(fs>=0.99).sum()
    pfs[(fs<0.99)*(fs>0.01)] = 1e-04/((fs<0.99)*(fs>0.01)).sum()
    pfs[fs<=0.01] = (1-1e-02-1e-04)/(fs<=0.01).sum()
    #print pfs.sum()
    pfs /= pfs.sum()
    
    # EM loop
    for iteration in xrange(0,num_iterations):
        log_pfs = numpy.log(pfs)
        
        log_posteriors = alts[:,None]*logfs[None,:]+refs[:,None]*log1minusfs[None,:]+numpy.log(pfs)[None,:]
        
        log_posteriors -= log_posteriors.max(axis=1)[:,None]
        
        posteriors = numpy.exp(log_posteriors)
        posteriors /= posteriors.sum(axis=1)[:,None]
        
        pfs = (posteriors*weights[:,None]).sum(axis=0)
        pfs = numpy.clip(pfs, 1e-100, 1e100)
    
        #print pfs.sum()
        
        # normalize
        pfs /= pfs.sum()
        
    return fs, pfs

def get_truong_pvalue(A,D):
    A = min([A,D-A])
    perr = 1e-02
    
    return scipy.stats.binom.sf(A,D,perr)+scipy.stats.binom.pmf(A,D,perr)
    
 
# definition of a polymorphic site according to Truong et al    
def is_polymorphic_truong(A,D):
    
    alpha = get_truong_pvalue(A,D)
    
    return alpha<0.05

def calculate_highcoverage_samples(species_name, min_coverage=config.min_median_coverage):
    
    # Load genomic coverage distributions
    sample_coverage_histograms, samples = parse_midas_data.parse_coverage_distribution(species_name)
    median_coverages = numpy.array([stats_utils.calculate_nonzero_median_from_histogram(sample_coverage_histogram) for sample_coverage_histogram in sample_coverage_histograms])
    sample_coverage_map = {samples[i]: median_coverages[i] for i in xrange(0,len(samples))}
    samples = numpy.array(samples)

    median_coverages = numpy.array([sample_coverage_map[samples[i]] for i in xrange(0,len(samples))])

    # Only plot samples above a certain depth threshold
    desired_samples = samples[(median_coverages>=min_coverage)]
    
    return desired_samples
    


def calculate_haploid_samples(species_name, min_coverage=config.min_median_coverage, threshold_pi=config.threshold_pi, threshold_within_between_fraction=config.threshold_within_between_fraction,debug=False):
    
    desired_samples = calculate_highcoverage_samples(species_name, min_coverage)
    
    if len(desired_samples)==0:
        return numpy.array([])
    
    # Old way, calculate pi_s
    # Load pi information for species_name
    # Load core gene set
    #sys.stderr.write("Loading core genes...\n")
    #core_genes = parse_midas_data.load_core_genes(species_name)
    #sys.stderr.write("Done! Core genome consists of %d genes\n" % len(core_genes))

    
    #sys.stderr.write("Loading within-sample diversity for %s...\n" % species_name)
    #samples, total_pis, total_pi_opportunities =     parse_midas_data.parse_within_sample_pi(species_name, allowed_genes=core_genes, debug=debug)
    #sys.stderr.write("Done!\n")
    #pis = total_pis/total_pi_opportunities

    #median_coverages = numpy.array([sample_coverage_map[samples[i]] for i in xrange(0,len(samples))])

    # Only plot samples above a certain depth threshold that are "haploids"
    #haploid_samples = samples[(median_coverages>=min_coverage)*(pis<=threshold_pi)]

    #return haploid_samples
    
    # New way with pre-computed SFS
    # Load SFS information for species_name
    import sfs_utils
    samples, sfs_map = parse_midas_data.parse_within_sample_sfs(species_name,     allowed_variant_types=set(['4D'])) 
    
    haploid_samples = []
    for sample in desired_samples:
        within_sites, between_sites, total_sites = sfs_utils.calculate_polymorphism_rates_from_sfs_map(sfs_map[sample])
    
        if within_sites <= threshold_within_between_fraction*between_sites:
            haploid_samples.append(sample)    
            
    return numpy.array(haploid_samples)

# Returns all high coverage samples that are involved in a temporal pair
def calculate_temporal_samples(species_name, min_coverage=config.min_median_coverage):

    highcoverage_samples = calculate_highcoverage_samples(species_name, min_coverage)
    
    sample_order_map = sample_utils.parse_sample_order_map()
    
    # Calculate which pairs of idxs belong to the same sample, which to the same subject
    # and which to different subjects
    same_sample_idxs, same_subject_idxs, diff_subject_idxs =     sample_utils.calculate_ordered_subject_pairs(sample_order_map, highcoverage_samples)


    temporal_samples = set()
    for sample_pair_idx in xrange(0,len(same_subject_idxs[0])):
   
        i = same_subject_idxs[0][sample_pair_idx]
        j = same_subject_idxs[1][sample_pair_idx]
        temporal_samples.add(highcoverage_samples[i])
        temporal_samples.add(highcoverage_samples[j])
        
    desired_samples = []
    for sample in highcoverage_samples:
        if sample in temporal_samples:
            desired_samples.append(sample)
            
    return numpy.array(desired_samples)

# Returns all high coverage samples that are involved in a temporal triplet
def calculate_triple_temporal_samples(species_name, min_coverage=config.min_median_coverage):

    highcoverage_samples = calculate_highcoverage_samples(species_name, min_coverage)
    
    sample_order_map = sample_utils.parse_sample_order_map()
    # Calculate which triplets of idxs belong to the same subject
    same_subject_idxs = parse_midas_data.calculate_ordered_subject_triplets(sample_order_map, highcoverage_samples)
    
    temporal_samples = set()
    for sample_triplet_idx in xrange(0,len(same_subject_idxs)):
        i,j,k = same_subject_idxs[sample_triplet_idx]
        
        temporal_samples.add(highcoverage_samples[i])
        temporal_samples.add(highcoverage_samples[j])
        temporal_samples.add(highcoverage_samples[k])
        
    desired_samples = []
    for sample in highcoverage_samples:
        if sample in temporal_samples:
            desired_samples.append(sample)
            
    return numpy.array(desired_samples)



def calculate_fixation_error_rate(sfs_map, sample_i, sample_j,dfs=[0.6], frequency_bins = numpy.linspace(0,1,21)):
    
    
    dfs = numpy.array(dfs)

    dummy_fs, pfs_i = sfs_utils.calculate_binned_sfs_from_sfs_map(sfs_map[sample_i],bins=frequency_bins)
    dummy_fs, pfs_j = sfs_utils.calculate_binned_sfs_from_sfs_map(sfs_map[sample_j],bins=frequency_bins)
    
    fs = frequency_bins[1:]-(frequency_bins[1]-frequency_bins[0])/2.0

    pfs = (pfs_i+pfs_j)/2.0
    # fold
    pfs = (pfs+pfs[::-1])/2
    
    # Calculate depth distributions
    dummy, D1s, pD1s = sfs_utils.calculate_binned_depth_distribution_from_sfs_map(sfs_map[sample_i])
    dummy, D2s, pD2s = sfs_utils.calculate_binned_depth_distribution_from_sfs_map(sfs_map[sample_j])
    
    fs = fs[pfs>0]
    pfs = pfs[pfs>0]
    
    D1s = D1s[pD1s>0]
    pD1s = pD1s[pD1s>0]
    D2s = D2s[pD2s>0]
    pD2s = pD2s[pD2s>0]
    
    
    perrs = {df:0 for df in dfs}
    for D1,pD1 in zip(D1s,pD1s):
        for D2,pD2 in zip(D2s,pD2s):
            for f,pf in zip(fs,pfs):
                for df in dfs:
                    perrs[df] += 2*binom.cdf(D1*(1-df)/2, D1, f)*binom.cdf(D2*(1-df)/2, D2, 1-f)*pD1*pD2*pf
    
    perrs = numpy.array([perrs[df] for df in dfs])
    return perrs




def find_snps_in_gene_pair(gene1_fasta, gene2_fasta):
    alignment={}
    # key=bp 
    # value=[B. vul, B. dorei]
    
    if len(gene1_fasta) == len(gene2_fasta):
        for bp in range(0, len(gene1_fasta)):
            if gene1_fasta[bp] != gene2_fasta[bp]:
                alignment[bp]=[gene1_fasta[bp],gene2_fasta[bp]]

    return alignment

    
