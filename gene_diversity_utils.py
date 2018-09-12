import numpy
import sys
from scipy.stats import poisson
from math import fabs
import config

# For each gene in gene_depth_matrix, calculates # of samples 
# in which it is "present". Returns vector of prevalences
def calculate_gene_prevalences(gene_depth_matrix, marker_coverages, min_copynum=0.3):

    gene_copynum_matrix = gene_depth_matrix * 1.0 / numpy.clip(marker_coverages,1,1e09)
    
    return (gene_copynum_matrix>=min_copynum).sum(axis=1)


def calculate_fractional_gene_prevalences(gene_depth_matrix, marker_coverages, min_copynum=0.3):    
    
    return calculate_gene_prevalences(gene_depth_matrix, marker_coverages, min_copynum)*1.0/len(marker_coverages)


# For each sample in gene_depth_matrix, calculates # of genes
# that are "present". Returns vector of gene numbers
def calculate_gene_numbers(gene_depth_matrix, marker_coverages, min_copynum=0.5):

    gene_copynum_matrix = gene_depth_matrix * 1.0 / numpy.clip(marker_coverages,1,1e09)
    
    return (gene_copynum_matrix>=min_copynum).sum(axis=0)

# Calculates the number of gene differences between pairs of samples
# Returns: matrix of # of gene differences
#          matrix of # of opportunities 
#          --> we should chat about what the denominator should be!
def calculate_coverage_based_gene_hamming_matrix(gene_reads_matrix, gene_depth_matrix, marker_coverages, absent_threshold=config.gainloss_max_absent_copynum, present_lower_threshold=config.gainloss_min_normal_copynum, present_upper_threshold= config.gainloss_max_normal_copynum):

    gene_hamming_matrix_gain, gene_hamming_matrix_loss, num_opportunities = calculate_coverage_based_gene_hamming_matrix_gain_loss(gene_reads_matrix, gene_depth_matrix, marker_coverages, absent_threshold=absent_threshold, present_lower_threshold=present_lower_threshold, present_upper_threshold=present_upper_threshold)
    
    gene_hamming_matrix = gene_hamming_matrix_gain + gene_hamming_matrix_loss 
    return gene_hamming_matrix, num_opportunities


def calculate_coverage_based_gene_hamming_matrix_gain_loss(gene_reads_matrix, gene_depth_matrix, marker_coverages, absent_threshold=config.gainloss_max_absent_copynum, present_lower_threshold=config.gainloss_min_normal_copynum, present_upper_threshold= config.gainloss_max_normal_copynum):
    # in this definition, we keep track of whether there was a gene 'gain' or 'loss' by removing numpy.fabs. This info will be used to plot the gain/loss events between two successive time pints. 

    gene_copynum_matrix = gene_depth_matrix*1.0/marker_coverages[None,:]
    #gene_copynum_matrix = numpy.clip(gene_copynum_matrix, 1e-09,1e09)

    # copynum is between 0.5 and 2
    is_present_copynum = (gene_copynum_matrix>present_lower_threshold)*(gene_copynum_matrix<present_upper_threshold)
    is_absent_copynum = (gene_copynum_matrix<=absent_threshold)
    is_low_copynum = (gene_copynum_matrix<present_upper_threshold)
  
  
  
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
        sub_is_present_copynum = is_present_copynum[lower_gene_idx:upper_gene_idx,:] 
        sub_is_absent_copynum = is_absent_copynum[lower_gene_idx:upper_gene_idx,:]
        sub_is_low_copynum = is_low_copynum[lower_gene_idx:upper_gene_idx,:]
    
        present_present = numpy.logical_and(sub_is_present_copynum[:,:,None], sub_is_present_copynum[:,None,:])
        
        present_absent = numpy.logical_and(sub_is_present_copynum[:,:,None], sub_is_absent_copynum[:,None,:])
    
        absent_present = numpy.logical_and(sub_is_absent_copynum[:,:,None], sub_is_present_copynum[:,None,:])
    
        gene_hamming_matrix_gain += (absent_present).sum(axis=0)
        gene_hamming_matrix_loss += (present_absent).sum(axis=0)
        num_opportunities += (present_present+absent_present+present_absent).sum(axis=0)
        
    return gene_hamming_matrix_gain, gene_hamming_matrix_loss, num_opportunities


# Calculate polarized gene copynum changes from i to j that exceed threshold 
# Returns list of differences. Each difference is a tuple of form 
#
# (gene_name, (cov_i, marker_cov_i), (cov_j, marker_cov_j))
#
def calculate_gene_differences_between(i, j, gene_reads_matrix, gene_depth_matrix, marker_coverages, absent_threshold=config.gainloss_max_absent_copynum, present_lower_threshold=config.gainloss_min_normal_copynum, present_upper_threshold= config.gainloss_max_normal_copynum):

    changed_genes = calculate_gene_differences_between_idxs(i, j, gene_reads_matrix, gene_depth_matrix, marker_coverages, absent_threshold=absent_threshold, present_lower_threshold=present_lower_threshold, present_upper_threshold=present_upper_threshold)

    gene_differences = []
    if len(changed_genes>0):
    
        # Look at these two samples
        gene_depth_matrix = gene_depth_matrix[:,[i,j]]
        marker_coverages = marker_coverages[[i,j]]
    
        for gene_idx in changed_genes:
            
            gene_differences.append( (gene_idx, (gene_depth_matrix[gene_idx,0], marker_coverages[0]), (gene_depth_matrix[gene_idx,1], marker_coverages[1])) )
            
    return gene_differences

def calculate_gene_differences_between_idxs(i, j, gene_reads_matrix, gene_depth_matrix, marker_coverages, absent_threshold=config.gainloss_max_absent_copynum, present_lower_threshold=config.gainloss_min_normal_copynum, present_upper_threshold= config.gainloss_max_normal_copynum):

    # Look at these two samples
    gene_depth_matrix = gene_depth_matrix[:,[i,j]]
    marker_coverages = marker_coverages[[i,j]]
    
    #marker_coverages = numpy.clip(marker_coverages,1,1e09)
    #gene_copynum_matrix = numpy.clip(gene_depth_matrix,1,1e09)*1.0/marker_coverages[None,:]

    gene_copynum_matrix = gene_depth_matrix*1.0/marker_coverages[None,:]
    gene_differences = []

    # copynum is between 0.5 and 2
    is_present_copynum = (gene_copynum_matrix>present_lower_threshold)*(gene_copynum_matrix<present_upper_threshold)
    is_absent_copynum = (gene_copynum_matrix<=absent_threshold)
    is_low_copynum = (gene_copynum_matrix<present_upper_threshold)
  
    present_present = numpy.logical_and(is_present_copynum[:,0], is_present_copynum[:,1])
        
    present_absent = numpy.logical_and(is_present_copynum[:,0], is_absent_copynum[:,1])
    
    absent_present = numpy.logical_and(is_absent_copynum[:,0], is_present_copynum[:,1])
    
    changed_genes = numpy.nonzero(numpy.logical_or(present_absent, absent_present))[0]

    return changed_genes
    


def calculate_triplet_gene_copynums(gene_depth_matrix, marker_coverages, i, j, k, absent_threshold=config.gainloss_max_absent_copynum, present_lower_threshold=config.gainloss_min_normal_copynum, present_upper_threshold= config.gainloss_max_normal_copynum):

    desired_samples = numpy.array([i, j, k])
    
    # Look at these three samples
    gene_depth_matrix = gene_depth_matrix[:,desired_samples]
    marker_coverages = marker_coverages[desired_samples]
    
    gene_copynum_matrix = gene_depth_matrix*1.0/marker_coverages[None,:]
    
    # copynum is between 0.5 and 2
    is_present_copynum = (gene_copynum_matrix>present_lower_threshold)*(gene_copynum_matrix<present_upper_threshold)
    is_absent_copynum = (gene_copynum_matrix<=absent_threshold)
    is_low_copynum = (gene_copynum_matrix<present_upper_threshold)
  
    changed_idxs = (is_low_copynum.all(axis=1))*(is_present_copynum.any(axis=1))*(is_absent_copynum.any(axis=1))*numpy.logical_or(is_present_copynum[:,0], is_absent_copynum[:,0])
  
    copynum_trajectories = []
    if changed_idxs.sum() > 0:
  
        gene_copynum_matrix = gene_copynum_matrix[changed_idxs]
  
        for gene_idx in xrange(0,gene_copynum_matrix.shape[0]):
            copynum_trajectories.append(gene_copynum_matrix[gene_idx,:])
        
    return copynum_trajectories


################################################################################
#
#  return a complete list of prevalences including all genes in the pangenome
#
################################################################################
def gene_prevalences_whole_pangenome(gene_names, gene_names_subset, prevalences):
    
    # first make a dictionary of gene names and prevalences:
    prevalence_dict={}
    for i in range(0, len(gene_names_subset)):
        gene=gene_names_subset[i]
        prevalence_dict[gene]=prevalences[i]

    gene_prevalences=[]
    for gene in gene_names:
        if gene in prevalence_dict.keys():
            gene_prevalences.append(prevalence_dict[gene])
        else:
            gene_prevalences.append(0)
            
    return numpy.asarray(gene_prevalences)


############################################
#                                          #
#  return a histogram  of kegg pathway IDs #
#                                          #  
############################################  

def kegg_pathways_histogram(kegg_ids, gene_names, gene_samples,gene_prevalences=[], spgenes=False):

    import pandas
  
    # if no gene_prevalences are provided, assume that every gene is in every person. 
    if len(gene_prevalences)==0:
        gene_prevalences=numpy.repeat(len(gene_samples),len(gene_names))/float(len(gene_samples))
    else:
        gene_prevalences=gene_prevalences/float(len(gene_samples))
        
  
    pathway_histogram={}
    pathway_description={}
    gene_idx=0
    for gene in gene_names:
        if gene in kegg_ids.keys():
            prevalence=gene_prevalences[gene_idx]
            pathways=kegg_ids[gene]
            for i in range(0, len(pathways)):
                pathway=pathways[i][0]
                description=pathways[i][1]
                if pathway not in pathway_histogram.keys():
                    pathway_histogram[pathway]=[prevalence]
                else:
                    pathway_histogram[pathway].append(prevalence)
                if spgenes==True:
                    pathway_description[pathway]=pathway
                else:
                    pathway_description[pathway]=description
            gene_idx +=1

    # create different prevalence bins for stacked histograms [100%, <0x<100, 0%]
    bins = numpy.asarray([0,0.1,0.5,0.9,1.0]) 
    pathway_counts_list={}
    for val in range(0, len(bins)): # the last value will be the total
        pathway_counts_list[val]=[]
    pathway_description_list=[]

    for pathway in pathway_histogram.keys():
        counts,dummy=numpy.histogram(pathway_histogram[pathway],bins)
        if pathway != '':
            pathway_description_list.append(pathway_description[pathway])
            for val in range(0, len(bins)-1):
                pathway_counts_list[val].append(counts[val])
            pathway_counts_list[len(bins)-1].append(sum(counts))

    # convert to dataframe:
    kegg_df={'total':pathway_counts_list[len(bins)-1],'names':pathway_description_list}
    for val in range(0, len(bins)-1):
        kegg_df[bins[val]]=pathway_counts_list[val]

    kegg_df=pandas.DataFrame(kegg_df)
    
#    sorted_kegg_df=pandas.DataFrame.sort(kegg_df, columns='total')
#    pathway_counts_list=sorted_kegg_df['counts'].tolist()
#    pathway_description_list=sorted_kegg_df['names'].tolist()   

#    return pathway_counts_list, pathway_description_list

    return kegg_df, pathway_description_list
    
    
def calculate_gene_error_rate(i, j, gene_reads_matrix, gene_depth_matrix, marker_coverages, absent_thresholds=[config.gainloss_max_absent_copynum], present_lower_threshold=config.gainloss_min_normal_copynum, present_upper_threshold= config.gainloss_max_normal_copynum):
   
    # Get reads, depths, and marker coverages in sample 1
    N1s = gene_reads_matrix[:,i]
    D1s = gene_depth_matrix[:,i]
    Dm1 = marker_coverages[i]
    
    # Get reads, depths, and marker coverages in sample 2
    N2s = gene_reads_matrix[:,j]
    D2s = gene_depth_matrix[:,j]
    Dm2 = marker_coverages[j]
    
    # Calculate gene copy numbers in both samples
    C1s = D1s*1.0/Dm1
    C2s = D2s*1.0/Dm2
    
    # Get list of genes that are in "normal" range in both samples
    good_idxs_1 = (C1s>=present_lower_threshold)*(C1s<=present_upper_threshold)
    good_idxs_2 = (C2s>=present_lower_threshold)*(C2s<=present_upper_threshold)
  
  
    # Calculate empirical prior distribution p(l) from SI 3.5.
    # For convenience, length_factor = 1/l
    length_factors = []
    if good_idxs_1.sum() > 0:
        length_factors.extend( D1s[good_idxs_1]*1.0/N1s[good_idxs_1] ) 
    if good_idxs_2.sum() > 0:
        length_factors.extend( D2s[good_idxs_2]*1.0/N2s[good_idxs_2] ) 
    length_factors = numpy.array(length_factors)
    # Will use whole array as prior distribution
     
    # Calculate empirical prior distribution p(c) from SI 3.5
    copynum_bins = numpy.linspace(0,2,21)
    Cs = copynum_bins[1:]-(copynum_bins[1]-copynum_bins[0])/2
    copynum_bins[0] = -1 # Just to make sure we include things with zero copynum
    copynum_bins[-1] = 1e09 # Assume things with c>2 have copynum 2 (conservative for detecting changes to lower values)
    
    # Prior distribution p(c) 
    pCs = numpy.histogram(numpy.hstack([C1s,C2s]), bins=copynum_bins)[0]
    pCs = pCs*1.0/pCs.sum()
    
    # Vectors for matrix calculation
    C_lowers = numpy.ones_like(Cs)*present_lower_threshold
    C_uppers = numpy.ones_like(Cs)*present_upper_threshold
    
    # The parameter of the poisson distribution of reads in sample 1
    Navg1s = Cs[:,None]*Dm1/length_factors[None,:]
    # Lower limit of number of reads for normal range
    Nlower1s = C_lowers[:,None]*Dm1/length_factors[None,:]
    # Upper limit of number of reads for normal range
    Nupper1s = C_uppers[:,None]*Dm1/length_factors[None,:]
    
    # Same for sample 2
    Navg2s = Cs[:,None]*Dm2/length_factors[None,:]
    Nlower2s = C_lowers[:,None]*Dm2/length_factors[None,:]
    Nupper2s = C_uppers[:,None]*Dm2/length_factors[None,:]
    
    perrs = []
    for absent_threshold in absent_thresholds:
        
        C_absents = numpy.ones_like(Cs)*absent_threshold
    
        # Upper limit of number of reads for "absent" range
        Nabsent1s = C_absents[:,None]*Dm1/length_factors[None,:]
        Nabsent2s = C_absents[:,None]*Dm2/length_factors[None,:]
        
        perr = 0
        perr += (poisson.cdf(Nabsent1s, Navg1s)*(poisson.cdf(Nupper2s,Navg2s)-poisson.cdf(Nlower2s,Navg2s))*pCs[:,None]).sum()
        perr += (poisson.cdf(Nabsent2s, Navg2s)*(poisson.cdf(Nupper1s,Navg1s)-poisson.cdf(Nlower1s,Navg1s))*pCs[:,None]).sum()
        # To emulate the integral over p(l)
        perr = perr/len(length_factors)
        
        # BG: 5/23/18. This is now done outside this function. 
        # To add up all the genes
        #perr = perr*len(C1s)
    
        perrs.append(perr)
        
    return numpy.array(perrs)


# Fuzzy matching of nearby genes
def is_nearby(gene_change_1, gene_change_2):
    
    gene_name_1 = gene_change_1[0]
    gene_name_2 = gene_change_2[0]
       
    gene_items_1 = gene_name_1.split(".")
    gene_items_2 = gene_name_2.split(".")
    
    genome_1 = ".".join([gene_items_1[0],gene_items_1[1]])
    gene_number_1 = long(gene_items_1[-1])
    genome_2 = ".".join([gene_items_2[0],gene_items_2[1]])
    gene_number_2 = long(gene_items_2[-1])
    
    if genome_1==genome_2:
        if fabs(gene_number_1-gene_number_2) < 6:
            return True        
    else:
        return False

def get_nearby_gene_idxs(gene_names, gene_idx, spacing=1, skip_target_gene=True):
    
    gene_name = gene_names[gene_idx]
    
    gene_items = gene_name.split(".")
    gene_id = long(gene_items[-1])
    
    
    idxs = []
    
    for i in xrange(-spacing,spacing+1):
    
        if skip_target_gene==True and i==0:
            continue
    
        gene_name = ".".join(gene_items[:-1]+[str(gene_id+i)])
        
        if gene_name in gene_names:
            idxs.append(gene_names.index(gene_name))
        
    return idxs
        
# Tries to merge nearby gene differences into blocks  
def merge_nearby_gene_differences(gene_differences):

    blocks = []
    
    for new_difference in gene_differences:
        print new_difference[0]
        matched=False
        for block_idx in xrange(0,len(blocks)):
            for old_difference in blocks[block_idx]:
                if is_nearby(new_difference, old_difference):
                    matched=True
                    break
            if matched:
                blocks[block_idx].append(new_difference)
                break
        
        if not matched:
            blocks.append([new_difference])
      
    print len(gene_differences), len(blocks)        
    return blocks
            
  
