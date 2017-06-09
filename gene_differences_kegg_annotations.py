import matplotlib  
matplotlib.use('Agg') 
import parse_midas_data
import pylab
import sys
import numpy

import diversity_utils
import gene_diversity_utils

import stats_utils
import matplotlib.colors as colors
import matplotlib.cm as cmx
from math import log10,ceil
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy.random import randint
import parse_patric
import pandas

################################################################################
#
# Standard header to read in argument information
#
################################################################################
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("species_name", help="name of species to process")
parser.add_argument("--debug", help="Loads only a subset of SNPs for speed", action="store_true")
parser.add_argument("--chunk-size", type=int, help="max number of records to load", default=1000000000)
args = parser.parse_args()

species_name = args.species_name
debug = args.debug
chunk_size = args.chunk_size
################################################################################

min_change = 0.8
min_coverage = 20
alpha = 0.5 # Confidence interval range for rate estimates
low_pi_threshold = 1e-03
low_divergence_threshold = 1e-03

# Load subject and sample metadata
sys.stderr.write("Loading HMP metadata...\n")
subject_sample_map = parse_midas_data.parse_subject_sample_map()
sys.stderr.write("Done!\n")
    
# Load core gene set
sys.stderr.write("Loading core genes...\n")
core_genes = parse_midas_data.load_core_genes(species_name)
sys.stderr.write("Done! Core genome consists of %d genes\n" % len(core_genes))
    
# Load genomic coverage distributions
sample_coverage_histograms, samples = parse_midas_data.parse_coverage_distribution(species_name)
median_coverages = numpy.array([stats_utils.calculate_nonzero_median_from_histogram(sample_coverage_histogram) for sample_coverage_histogram in sample_coverage_histograms])
sample_coverage_map = {samples[i]: median_coverages[i] for i in xrange(0,len(samples))}

# Load pi information for species_name
sys.stderr.write("Loading within-sample diversity for %s...\n" % species_name)
samples, total_pis, total_pi_opportunities = parse_midas_data.parse_within_sample_pi(species_name, allowed_genes=core_genes, debug=debug)
sys.stderr.write("Done!\n")
pis = total_pis/total_pi_opportunities

median_coverages = numpy.array([sample_coverage_map[samples[i]] for i in xrange(0,len(samples))])

# Only plot samples above a certain depth threshold that are "haploids"
snp_samples = samples[(median_coverages>=min_coverage)*(pis<=1e-03)]

# Analyze SNPs, looping over chunk sizes. 
# Clunky, but necessary to limit memory usage on cluster

# Load SNP information for species_name
sys.stderr.write("Loading SNPs for %s...\n" % species_name)
sys.stderr.write("(not just core genes...)\n")
pi_matrix_syn = numpy.array([])
avg_pi_matrix_syn = numpy.array([])
snp_difference_matrix = numpy.array([])
snp_opportunity_matrix = numpy.array([])
final_line_number = 0

while final_line_number >= 0:    
    sys.stderr.write("Loading chunk starting @ %d...\n" % final_line_number)
    dummy_samples, allele_counts_map, passed_sites_map, final_line_number = parse_midas_data.parse_snps(species_name, debug=debug, allowed_samples=snp_samples,chunk_size=chunk_size,initial_line_number=final_line_number)
    sys.stderr.write("Done! Loaded %d genes\n" % len(allele_counts_map.keys()))
#
    # Calculate fixation matrix
    sys.stderr.write("Calculating matrix of snp differences...\n")
    chunk_snp_difference_matrix, chunk_snp_opportunity_matrix =     diversity_utils.calculate_fixation_matrix(allele_counts_map, passed_sites_map, min_change=min_change)   # 
    sys.stderr.write("Done!\n")
#
    if snp_difference_matrix.shape[0]==0:
        snp_difference_matrix = numpy.zeros_like(chunk_snp_difference_matrix)*1.0
        snp_opportunity_matrix = numpy.zeros_like(snp_difference_matrix)*1.0
#    
    snp_difference_matrix += chunk_snp_difference_matrix
    snp_opportunity_matrix += chunk_snp_opportunity_matrix


snp_substitution_rate = snp_difference_matrix*1.0/(snp_opportunity_matrix+(snp_opportunity_matrix==0))
sys.stderr.write("Done!\n")   

# Load gene coverage information for species_name
sys.stderr.write("Loading pangenome data for %s...\n" % species_name)
gene_samples, gene_names, gene_presence_matrix, gene_depth_matrix, marker_coverages, gene_reads_matrix = parse_midas_data.parse_pangenome_data(species_name,allowed_samples=snp_samples,convert_centroid_names=False)
sys.stderr.write("Done!\n")

prevalence_idxs = (parse_midas_data.calculate_unique_samples(subject_sample_map, gene_samples))*(marker_coverages>=min_coverage)
    
prevalences = gene_diversity_utils.calculate_fractional_gene_prevalences(gene_depth_matrix[:,prevalence_idxs], marker_coverages[prevalence_idxs])

pangenome_prevalences = numpy.array(prevalences,copy=True)
pangenome_prevalences.sort()
        
# Calculate matrix of number of genes that differ
sys.stderr.write("Calculating matrix of gene differences...\n")
gene_difference_matrix, gene_opportunity_matrix = gene_diversity_utils.calculate_coverage_based_gene_hamming_matrix(gene_depth_matrix, marker_coverages, min_log2_fold_change=4)

##############################################################
# Now need to make the gene samples and snp samples match up #
# IDXs for same subject vs diff subject
##############################################################
desired_samples = gene_samples[marker_coverages>min_coverage]      

num_haploids = len(desired_samples)
     
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

typical_same_subject_snp_opportunities = numpy.median(snp_opportunity_matrix[same_subject_snp_idxs])
typical_diff_subject_snp_opportunities = numpy.median(snp_opportunity_matrix[diff_subject_snp_idxs])

typical_same_subject_gene_opportunities = numpy.median(gene_opportunity_matrix[same_subject_gene_idxs])
typical_diff_subject_gene_opportunities = numpy.median(gene_opportunity_matrix[diff_subject_gene_idxs])


###############################################
# Load kegg information 
##############################################
# load the kegg information for this species
kegg_ids=parse_patric.load_kegg_annotations(gene_names)


#######################################################################
# compute the kegg pathways in which genes that are changing are in   #
# first look at within-host changes                                   #
#######################################################################

pathways_within_host_changes={}
gene_count_within_host_changes={}
max_CNV_value_within=[] # I will store the higher CNV value for the pair of samples being compared
for sample_pair_idx in xrange(0,len(same_subject_snp_idxs[0])):
    i = same_subject_gene_idxs[0][sample_pair_idx]
    j = same_subject_gene_idxs[1][sample_pair_idx]
    gene_differences = gene_diversity_utils.calculate_gene_differences_between(i, j, gene_depth_matrix, marker_coverages)
    # iterate through and pull out the gene IDx and then map to kegg
    for gene_change in gene_differences:
        # get the gene index
        gene_idx = gene_change[0]
        # compute the CNV values
        gene_depth_i=gene_change[1][0]
        marker_coverage_i=gene_change[1][1]
        CNV_i=gene_depth_i/marker_coverage_i
        gene_depth_j=gene_change[2][0]
        marker_coverage_j=gene_change[2][1]
        CNV_j=gene_depth_j/marker_coverage_j
        if CNV_i > CNV_j: # store the higher CNV value for plotting later on. 
            max_CNV_value_within.append(CNV_i)
        else:
            max_CNV_value_within.append(CNV_j) 
        # get the gene name
        gene_name=gene_names[gene_idx]
        if gene_name not in gene_count_within_host_changes:
            gene_count_within_host_changes[gene_name]=1
        else:
            gene_count_within_host_changes[gene_name]+=1  
        kegg_id=kegg_ids[gene_name]
        for pathway_id in range(0, len(kegg_id)):
            pathway=kegg_id[pathway_id][1]
            if pathway=='':
                pathway='Unannotated pathways'
            if pathway not in pathways_within_host_changes:
                pathways_within_host_changes[pathway]=1
            else: 
                pathways_within_host_changes[pathway]+=1


#########################################
# Kegg pathways of between-host changes #
#########################################

pathways_between_host_changes={}
gene_count_between_host_changes={}
max_CNV_value_between=[] # I will store the higher CNV value for the pair of samples being compared
for sample_pair_idx in xrange(0,len(diff_subject_snp_idxs[0])): 
    i = diff_subject_gene_idxs[0][sample_pair_idx]
    j = diff_subject_gene_idxs[1][sample_pair_idx]
    gene_differences = gene_diversity_utils.calculate_gene_differences_between(i, j, gene_depth_matrix, marker_coverages)
    # iterate through and pull out the gene IDx and then map to kegg
    for gene_change in gene_differences:
        gene_idx = gene_change[0]
        # compute the CNV values
        gene_depth_i=gene_change[1][0]
        marker_coverage_i=gene_change[1][1]
        CNV_i=gene_depth_i/marker_coverage_i
        gene_depth_j=gene_change[2][0]
        marker_coverage_j=gene_change[2][1]
        CNV_j=gene_depth_j/marker_coverage_j
        if CNV_i > CNV_j: # store the higher CNV value for plotting later on. 
            max_CNV_value_between.append(CNV_i)
        else:
            max_CNV_value_between.append(CNV_j) 
        # get the gene name
        gene_name=gene_names[gene_idx]
        if gene_name not in gene_count_between_host_changes:
            gene_count_between_host_changes[gene_name]=1
        else:
            gene_count_between_host_changes[gene_name]+=1  
        kegg_id=kegg_ids[gene_name]
        for pathway_id in range(0, len(kegg_id)):
            pathway=kegg_id[pathway_id][1]
            if pathway!='': # exclude unannotated pathways because this dominates the plot
                if pathway not in pathways_between_host_changes:
                    pathways_between_host_changes[pathway]=1
                else: 
                    pathways_between_host_changes[pathway]+=1


###########################################################################
# Plot: histogram of pathways showing number of within-host gene changes  #
###########################################################################  

# try with dataframe
label_size = 10
matplotlib.rcParams['ytick.labelsize'] = label_size
pos = numpy.arange(len(pathways_within_host_changes))
width = 1.0 

df=pandas.DataFrame(pathways_within_host_changes, index=[0])
df=pandas.melt(df)
df=df.sort(['value'], ascending=[0])
ax =df.plot(kind='barh', stacked=False, figsize=(6,2), title=species_name, width=width)
ax.set_xlabel("Number of genes x number of hosts")
ax.set_ylabel("Kegg pathway")
ax.set_yticks(pos + (width / 2))
ax.set_yticklabels(df['variable'].tolist())

pylab.savefig('%s/%s_kegg_gene_changes_within.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight') 


# repeat for between
label_size = 8
matplotlib.rcParams['ytick.labelsize'] = label_size
pos = numpy.arange(len(pathways_between_host_changes))
width = 1.0 

df=pandas.DataFrame(pathways_between_host_changes, index=[0])
df=pandas.melt(df)
df=df.sort(['value'], ascending=[0])
ax =df.plot(kind='barh', stacked=False, figsize=(6,12), title=species_name, width=width)
ax.set_xlabel("Number of genes x number of host pairs")
ax.set_ylabel("Kegg pathway")
ax.set_yticks(pos + (width / 2))
ax.set_yticklabels(df['variable'].tolist())

pylab.savefig('%s/%s_kegg_gene_changes_between.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight') 



#########################################################################################################
# plot a histogram of the number of times a unique gene is found to be changed between different hosts  #
#########################################################################################################  

# within
pylab.figure(figsize=(6,6))
pylab.subplot(2, 1, 1)     
pylab.title(species_name +', within host')

plt.hist(gene_count_within_host_changes.values(), normed=False, bins=100)
plt.ylabel('Number of genes');

# between
pylab.subplot(2, 1, 2)
pylab.xlabel('Number of times a gene is changed')
pylab.title('Between host')

plt.hist(gene_count_between_host_changes.values(), normed=False, bins=100)
plt.ylabel('Number of genes');

pylab.savefig('%s/%s_frequency_of_gene_change.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight',dpi=300)


###########################################
# what is the distribution of CNV values? #
###########################################   
marker_coverages = numpy.clip(marker_coverages,1,1e09)
gene_copynum_matrix = numpy.clip(gene_depth_matrix,1,1e09)*1.0/marker_coverages[None,:]
copynum_distribution=list(gene_copynum_matrix.flat)    # if you want a list

# plot a histogram:
pylab.figure(figsize=(6,10))
pylab.subplot(3, 1, 1)
pylab.title(species_name + ', CNV distribution')
pylab.ylim(0,100)
pylab.xlim(0,210)
plt.hist(copynum_distribution, normed=False, bins=100)

# plot the distribution of CNV values for within host changes

pylab.subplot(3, 1, 2)
pylab.title('Within-host changes')
pylab.xlim(0,210)  
plt.hist(max_CNV_value_within, normed=False, bins=10)
plt.ylabel('Number of genes x hosts');

# plot the distribution of CNV values for between host changes
pylab.subplot(3, 1, 3)   
pylab.xlabel('CNV')
pylab.title('Between-host changes')
pylab.ylim(0,100) 
pylab.xlim(0,210) 
plt.hist(max_CNV_value_between, normed=False, bins=100)

pylab.savefig('%s/%s_CNV_distribution.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight',dpi=300)

#######################################################################################################################
# print the list of within-host gene changes. Outside of python, try blasting these genes to other species' genomes   #
#######################################################################################################################

outFN=('%s/%s_within_host_gene_changes.txt' % (parse_midas_data.analysis_directory,species_name))
outFile=open(outFN,'w')

for gene in gene_count_within_host_changes:
    print gene
    outFile.write(gene +'\n')
