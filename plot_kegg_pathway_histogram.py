import matplotlib  
matplotlib.use('Agg') 
import parse_midas_data
import pylab
import sys
import numpy
from numpy.random import normal
import diversity_utils
import gene_diversity_utils
import stats_utils
import os
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

#############################################################################

# Minimum median coverage of sample to look at
min_coverage = 20

#################
# Load metadata #
#################

# Load subject and sample metadata
sys.stderr.write("Loading HMP metadata...\n")
subject_sample_map = parse_midas_data.parse_subject_sample_map()
sys.stderr.write("Done!\n")
    
######################
# Load coverage data #
######################

# Load genomic coverage distributions
sample_coverage_histograms, samples = parse_midas_data.parse_coverage_distribution(species_name)
median_coverages = numpy.array([stats_utils.calculate_median_from_histogram(sample_coverage_histogram) for sample_coverage_histogram in sample_coverage_histograms])
sample_coverage_map = {samples[i]: median_coverages[i] for i in xrange(0,len(samples))}

###############################################################
# Compute Pi within patients to figure out which are haploid  #
###############################################################

# Load pi information for species_name
sys.stderr.write("Loading within-sample diversity for %s...\n" % species_name)
samples, total_pis, total_pi_opportunities = parse_midas_data.parse_within_sample_pi(species_name, debug)
sys.stderr.write("Done!\n")
pis = total_pis/total_pi_opportunities

######################
# compute median cov #
######################

median_coverages = numpy.array([sample_coverage_map[samples[i]] for i in xrange(0,len(samples))])

###############################################################
# Indexes for SNP samples that have high coverage #
###############################################################

# Only plot samples above a certain depth threshold that are "haploids"
low_pi_snp_samples = samples[(median_coverages>=min_coverage)*(pis<=1e-03)]
high_pi_snp_samples = samples[(median_coverages>=min_coverage)*(pis>1e-03)]

# Calculate which pairs of idxs belong to the same sample, which to the same subject
# and which to different subjects
#snp_same_sample_idxs, snp_same_subject_idxs, snp_diff_subject_idxs = parse_midas_data.calculate_subject_pairs(subject_sample_map, snp_samples)

####################################################
# Load gene coverage information for species_name
####################################################
sys.stderr.write("Loading pangenome data for %s...\n" % species_name)
low_pi_gene_samples, low_pi_gene_names, low_pi_gene_presence_matrix, low_pi_gene_depth_matrix, low_pi_marker_coverages, low_pi_gene_reads_matrix = parse_midas_data.parse_pangenome_data(species_name,allowed_samples=low_pi_snp_samples)

high_pi_gene_samples, high_pi_gene_names, high_pi_gene_presence_matrix, high_pi_gene_depth_matrix, high_pi_marker_coverages, high_pi_gene_reads_matrix = parse_midas_data.parse_pangenome_data(species_name,allowed_samples=high_pi_snp_samples)
sys.stderr.write("Done!\n")


###############################################
# Load kegg information 
##############################################

# redo loading of gene_names -- currently requires that gene is present in at least one ind. 

# load the kegg information for this species
kegg_ids=parse_midas_data.load_kegg_annotations(low_pi_gene_names)

# iterate through the kegg_ids and make a histogram of what is there. 
#kegg_df, pathway_description_list=gene_diversity_utils.kegg_pathways_histogram(kegg_ids, gene_names, gene_samples)


############################################################## 
# Distribution of presence/absence of genes accross samples  #
############################################################## 

# iterate through is_good_copynum to identify the number of samples in which a gene shows up
low_pi_gene_prevalences=gene_diversity_utils.calculate_gene_prevalences(low_pi_gene_depth_matrix, low_pi_marker_coverages, min_copynum=0.5)

high_pi_gene_prevalences=gene_diversity_utils.calculate_gene_prevalences(high_pi_gene_depth_matrix, high_pi_marker_coverages, min_copynum=0.5)

# using the gene_names list, figure out what pathways are present at which frequency across samples. 
low_pi_kegg_df, low_pi_pathway_description_list=gene_diversity_utils.kegg_pathways_histogram(kegg_ids, low_pi_gene_names, low_pi_gene_samples, low_pi_gene_prevalences)

high_pi_kegg_df, high_pi_pathway_description_list=gene_diversity_utils.kegg_pathways_histogram(kegg_ids, high_pi_gene_names, high_pi_gene_samples,high_pi_gene_prevalences)




##########################################
# plot stacked histogram with pandas df  #
########################################## 
# low pi

label_size = 8
matplotlib.rcParams['ytick.labelsize'] = label_size 
pos = numpy.arange(len(low_pi_pathway_description_list))
width = 1.0 

kegg_df=pandas.DataFrame.sort(low_pi_kegg_df, columns='total')
ax =kegg_df[[0.9,0.5,0.1,0.0]].plot(kind='barh', stacked=True, figsize=(10,18), title=species_name, color=['r','b','b','b'])
ax.set_xlabel("Number of genes")
ax.set_ylabel("Kegg pathway")
ax.set_yticks(pos + (width / 2))
ax.set_yticklabels(kegg_df['names'].tolist())

pylab.savefig('%s/%s_low_pi_kegg_histogram.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight') 



# plot for all prevalences:
label_size = 8
matplotlib.rcParams['ytick.labelsize'] = label_size 
pos = numpy.arange(len(low_pi_pathway_description_list))
width = 1.0 

kegg_df=pandas.DataFrame.sort(low_pi_kegg_df, columns='total')
ax =kegg_df[[0.9,0.5,0.1,0.0]].plot(kind='barh', stacked=True, figsize=(10,18), title=species_name, color=['r','y','g','b'])
ax.set_xlabel("Number of genes")
ax.set_ylabel("Kegg pathway")
ax.set_yticks(pos + (width / 2))
ax.set_yticklabels(kegg_df['names'].tolist())

pylab.savefig('%s/%s_low_pi_kegg_histogram_all_prevalences.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight') 


# high pi

label_size = 8
matplotlib.rcParams['ytick.labelsize'] = label_size 
pos = numpy.arange(len(high_pi_pathway_description_list))
width = 1.0 

kegg_df=pandas.DataFrame.sort(high_pi_kegg_df, columns='total')
ax =kegg_df[[0.9,0.5,0.1,0.0]].plot(kind='barh', stacked=True, figsize=(10,18), title=species_name, color=['r','b','b','b'])
ax.set_xlabel("Number of genes")
ax.set_ylabel("Kegg pathway")
ax.set_yticks(pos + (width / 2))
ax.set_yticklabels(kegg_df['names'].tolist())

pylab.savefig('%s/%s_high_pi_kegg_histogram.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight') 



# plot for all prevalences:
label_size = 8
matplotlib.rcParams['ytick.labelsize'] = label_size 
pos = numpy.arange(len(high_pi_pathway_description_list))
width = 1.0 

kegg_df=pandas.DataFrame.sort(high_pi_kegg_df, columns='total')
ax =kegg_df[[0.9,0.5,0.1,0.0]].plot(kind='barh', stacked=True, figsize=(10,18), title=species_name, color=['r','y','g','b'])
ax.set_xlabel("Number of genes")
ax.set_ylabel("Kegg pathway")
ax.set_yticks(pos + (width / 2))
ax.set_yticklabels(kegg_df['names'].tolist())

pylab.savefig('%s/%s_high_pi_kegg_histogram_all_prevalences.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight') 





# also plot later: what is the distribution of copy number values? Why is the max all the way up to 112? What properties does this gene have? I.e. could it be some weird pathway? 


