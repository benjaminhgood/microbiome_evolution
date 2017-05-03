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
import parse_patric
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

sys.stderr.write("Loading core genes...\n")
core_genes = parse_midas_data.load_core_genes(species_name)
sys.stderr.write("Done! %d core genes\n" % len(core_genes))

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
samples, total_pis, total_pi_opportunities = parse_midas_data.parse_within_sample_pi(species_name, allowed_variant_types=set(['4D']), allowed_genes=core_genes, debug=debug)
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

# Calculate which pairs of idxs belong to unique samples. Remove any samples that are duplicates (i.e. multiple time pts)
unique_idxs=parse_midas_data.calculate_unique_samples(subject_sample_map, low_pi_snp_samples)
low_pi_snp_samples=low_pi_snp_samples[unique_idxs]

unique_idxs=parse_midas_data.calculate_unique_samples(subject_sample_map, high_pi_snp_samples)
high_pi_snp_samples=low_pi_snp_samples[unique_idxs]

####################################################
# Load gene coverage information for species_name
####################################################
sys.stderr.write("Loading pangenome data for %s...\n" % species_name)
low_pi_gene_samples, low_pi_gene_names, low_pi_gene_presence_matrix, low_pi_gene_depth_matrix, low_pi_marker_coverages, low_pi_gene_reads_matrix = parse_midas_data.parse_pangenome_data(species_name,allowed_samples=low_pi_snp_samples)

high_pi_gene_samples, high_pi_gene_names, high_pi_gene_presence_matrix, high_pi_gene_depth_matrix, high_pi_marker_coverages, high_pi_gene_reads_matrix = parse_midas_data.parse_pangenome_data(species_name,allowed_samples=high_pi_snp_samples)
sys.stderr.write("Done!\n")

# this represents all gene names regardless of prevalences
gene_names, new_species_names=list(parse_midas_data.load_pangenome_genes(species_name))

# convert format of gene names from set to list:
gene_names=list(gene_names)


###############################################
# Load spgenes
###############################################  
sys.stderr.write("Loading spgenes for %s...\n" % species_name)
spgenes_ids, spgenes_set=parse_patric.load_spgenes_annotations(gene_names)
sys.stderr.write("Done!\n")

###############################################
# Load kegg information 
##############################################
# load the kegg information for this species
kegg_ids=parse_midas_data.load_kegg_annotations(gene_names)


############################################################## 
# Distribution of presence/absence of genes accross samples  #
############################################################## 

# iterate through is_good_copynum to identify the number of samples in which a gene shows up
low_pi_gene_prevalences=gene_diversity_utils.calculate_gene_prevalences(low_pi_gene_depth_matrix, low_pi_marker_coverages)

high_pi_gene_prevalences=gene_diversity_utils.calculate_gene_prevalences(high_pi_gene_depth_matrix, high_pi_marker_coverages)

# the prevalence given above is only for the subset where prev is at least 1. 
# want the full set, so merge low_pi_gene_names with gene_names and assign prev =0 to those gene names where there is no overlap.

low_pi_gene_prevalences_pangenome=gene_diversity_utils.gene_prevalences_whole_pangenome(gene_names, low_pi_gene_names, low_pi_gene_prevalences)
high_pi_gene_prevalences_pangenome=gene_diversity_utils.gene_prevalences_whole_pangenome(gene_names, high_pi_gene_names, high_pi_gene_prevalences)

# using the gene_names list, figure out what pathways are present at which frequency across samples. 

low_pi_kegg_df, low_pi_pathway_description_list=gene_diversity_utils.kegg_pathways_histogram(spgenes_ids, gene_names, low_pi_gene_samples, low_pi_gene_prevalences_pangenome)

high_pi_kegg_df, high_pi_pathway_description_list=gene_diversity_utils.kegg_pathways_histogram(kegg_ids, gene_names, high_pi_gene_samples,high_pi_gene_prevalences_pangenome)

# using the gene_names list, figure out what pathways are present at which frequency across samples in the spgenes 

low_pi_spgenes_df, low_pi_pathway_description_list=gene_diversity_utils.kegg_pathways_histogram(spgenes_ids, gene_names, low_pi_gene_samples, low_pi_gene_prevalences_pangenome, spgenes=True)

high_pi_spgenes_df, high_pi_pathway_description_list=gene_diversity_utils.kegg_pathways_histogram(spgenes_ids, gene_names, high_pi_gene_samples,high_pi_gene_prevalences_pangenome, spgenes=True)



#############################################################################
# plot stacked histogram with pandas df                                     #
# Each stack represents the numebr of genes in different prevalence classes #
########################################## ##################################
# low pi

# plot for all prevalences:
label_size = 10
matplotlib.rcParams['ytick.labelsize'] = label_size 
pos = numpy.arange(len(low_pi_pathway_description_list))
width = 1.0 

spgenes_df=pandas.DataFrame.sort(low_pi_spgenes_df, columns='total')
ax =spgenes_df[[0.9,0.5,0.1,0.0]].plot(kind='barh', stacked=True, figsize=(4,8), title=species_name, color=['r','y','g','b'])
ax.set_xlabel("Number of genes")
ax.set_ylabel("Type of special gene")
ax.set_yticks(pos + (width / 2))
ax.set_yticklabels(spgenes_df['names'].tolist())

pylab.savefig('%s/%s_low_pi_spgenes_histogram_all_prevalences.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight') 


# high pi

# plot for all prevalences:
label_size = 10
matplotlib.rcParams['ytick.labelsize'] = label_size 
pos = numpy.arange(len(high_pi_pathway_description_list))
width = 1.0 

spgenes_df=pandas.DataFrame.sort(high_pi_spgenes_df, columns='total')
ax =spgenes_df[[0.9,0.5,0.1,0.0]].plot(kind='barh', stacked=True, figsize=(4,8), title=species_name, color=['r','y','g','b'])
ax.set_xlabel("Number of genes")
ax.set_ylabel("Type of special gene")
ax.set_yticks(pos + (width / 2))
ax.set_yticklabels(spgenes_df['names'].tolist())

pylab.savefig('%s/%s_high_pi_spgenes_histogram_all_prevalences.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight') 


'''
##########################################################################
# plot side by side plot of number of genes in pathways for high vs low pis
#########################################################################

# create a new df with the low and high pi values for prevalence between 0.9 and 1
low_high_pi_kegg_df={'low_pis':low_pi_kegg_df[0.9], 'high_pis':high_pi_kegg_df[0.9], 'names':low_pi_kegg_df['names'], 'total':low_pi_kegg_df['total']}
low_high_pi_kegg_df = pandas.DataFrame(low_high_pi_kegg_df)


label_size = 8
matplotlib.rcParams['ytick.labelsize'] = label_size 
pos = numpy.arange(len(high_pi_pathway_description_list))
width = 1.0

kegg_df=pandas.DataFrame.sort(low_high_pi_kegg_df, columns='total')
ax =kegg_df[['low_pis','high_pis','total']].plot(kind='barh', stacked=False, figsize=(10,18), title=species_name, color=['r','y','b'], width=width)
ax.set_xlabel("Number of genes")
ax.set_ylabel("Kegg pathway")
ax.set_yticks(pos + (width / 2))
ax.set_yticklabels(kegg_df['names'].tolist())

pylab.savefig('%s/%s_spgenes_histogram_high_vs_low_pis_prevalence_0.9.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight') 


#######################################
# plot again with less prevalent genes
#######################################
low_high_pi_kegg_df={'low_pis':low_pi_kegg_df[0.5], 'high_pis':high_pi_kegg_df[0.5], 'names':low_pi_kegg_df['names'], 'total':low_pi_kegg_df['total']}
low_high_pi_kegg_df = pandas.DataFrame(low_high_pi_kegg_df)


label_size = 8
matplotlib.rcParams['ytick.labelsize'] = label_size 
pos = numpy.arange(len(high_pi_pathway_description_list))
width = 1.0

kegg_df=pandas.DataFrame.sort(low_high_pi_kegg_df, columns='total')
ax =kegg_df[['low_pis','high_pis','total']].plot(kind='barh', stacked=False, figsize=(10,18), title=species_name, color=['r','y','b'], width=width)
ax.set_xlabel("Number of genes")
ax.set_ylabel("Kegg pathway")
ax.set_yticks(pos + (width / 2))
ax.set_yticklabels(kegg_df['names'].tolist())

pylab.savefig('%s/%s_pi_kegg_histogram_high_vs_low_pis_prevalence_0.5.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight') 






low_high_pi_kegg_df={'low_pis':low_pi_kegg_df[0.1], 'high_pis':high_pi_kegg_df[0.1], 'names':low_pi_kegg_df['names'], 'total':low_pi_kegg_df['total']}
low_high_pi_kegg_df = pandas.DataFrame(low_high_pi_kegg_df)


label_size = 8
matplotlib.rcParams['ytick.labelsize'] = label_size 
pos = numpy.arange(len(high_pi_pathway_description_list))
width = 1.0

kegg_df=pandas.DataFrame.sort(low_high_pi_kegg_df, columns='total')
ax =kegg_df[['low_pis','high_pis','total']].plot(kind='barh', stacked=False, figsize=(10,18), title=species_name, color=['r','y','b'], width=width)
ax.set_xlabel("Number of genes")
ax.set_ylabel("Kegg pathway")
ax.set_yticks(pos + (width / 2))
ax.set_yticklabels(kegg_df['names'].tolist())

pylab.savefig('%s/%s_pi_kegg_histogram_high_vs_low_pis_prevalence_0.1.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight') 



####################################################
# how many species share a single pathway?         #
####################################################

'''
