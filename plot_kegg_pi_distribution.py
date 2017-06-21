import matplotlib  
matplotlib.use('Agg') 
import parse_midas_data
import pylab
from pylab import *
import sys
import numpy
from numpy.random import normal
import diversity_utils
import gene_diversity_utils
import stats_utils
import os
import pandas
import parse_patric
import pickle
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
min_passed_sites_per_gene=10
min_passed_sites_per_person=100

#core genes
sys.stderr.write("Loading core genes...\n")
core_genes = parse_midas_data.load_core_genes(species_name)
sys.stderr.write("Done! %d core genes\n" % len(core_genes))

# variable genes -- this is a list of anything that isn't a core gene but in the reference genome. 


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
# Indexes for SNP samples that have high coverage             #
###############################################################

# Only plot samples above a certain depth threshold that are "haploids"
low_pi_snp_samples = samples[(median_coverages>=min_coverage)*(pis<=1e-03)]
high_pi_snp_samples = samples[(median_coverages>=min_coverage)*(pis>1e-03)]

####################################################
# Load gene coverage information for species_name
####################################################

# Load all the genes from the reference genome regardless of prevalence. 
gene_names=parse_midas_data.load_reference_genes(species_name)
# convert format of gene names from set to list:
gene_names=list(gene_names)


# Variable genes: find the difference between gene_names and core_genes
gene_names_tmp=numpy.asarray(gene_names)
core_genes_tmp=numpy.asarray(list(core_genes))
variable_genes=set(numpy.asarray(list(numpy.setdiff1d(gene_names_tmp,core_genes_tmp))))

###############################################
# Load kegg information 
##############################################
# load the kegg information for this species
kegg_ids=parse_patric.load_kegg_annotations(gene_names)

##########################################################
# load SNP info
##########################################################
sys.stderr.write("Loading %s...\n" % species_name)
samples, allele_counts_map, passed_sites_map, final_line_number = parse_midas_data.parse_snps(species_name, debug, allowed_samples=low_pi_snp_samples)
sys.stderr.write("Done!\n")

#########################################################  
# figure out which samples belong to the same subject   #
######################################################### 
same_sample_idxs, same_subject_idxs, diff_subject_idxs = parse_midas_data.calculate_subject_pairs(subject_sample_map, samples)

###########################################################
# compute total pi genome-wide core genes
###########################################################
pi_matrix, avg_pi_matrix, passed_sites=diversity_utils.calculate_pi_matrix(allele_counts_map, passed_sites_map, variant_type='4D', allowed_genes=core_genes)
pi_matrix_core = pi_matrix /(passed_sites+(passed_sites==0))
avg_pi_matrix_core = avg_pi_matrix/(passed_sites+(passed_sites==0)) 

###########################################################
# compute total pi variable genes
###########################################################
pi_matrix, avg_pi_matrix, passed_sites=diversity_utils.calculate_pi_matrix(allele_counts_map, passed_sites_map, variant_type='4D', allowed_genes=variable_genes)
pi_matrix_variable = pi_matrix /(passed_sites+(passed_sites==0))
avg_pi_matrix_variable = avg_pi_matrix/(passed_sites+(passed_sites==0)) 

############################################################
# compute pi/pathway -- core genes
############################################################
avg_pi_per_gene={}
pi_per_gene={}
passed_sites_per_gene={}
num_people_with_data={}
for gene_name in core_genes:
    gene_pi_matrix, gene_avg_pi_matrix, gene_passed_sites= diversity_utils.calculate_pi_matrix(allele_counts_map, passed_sites_map, variant_type='4D', allowed_genes=[gene_name])
    # check if the values of gene_passed_sites are less than 5. If so, then zero out these idxs for gene_passed_sites, gene_pi_matrix, and gene_avg_pi_matrix. Basically all of these people or pairs of people have too few sites to compute realiable statistics. 
    low_passed_sites_idxs=(gene_passed_sites)<min_passed_sites_per_gene
    gene_passed_sites[low_passed_sites_idxs] =0
    gene_pi_matrix[low_passed_sites_idxs]=0
    gene_avg_pi_matrix[low_passed_sites_idxs]=0
    # put all the data into the dictionary for later aggregation. 
    # don't add any genes that have zero information.
    if sum(numpy.diagonal(gene_passed_sites)>min_passed_sites_per_gene)>0:
        passed_sites_per_gene[gene_name]=gene_passed_sites
        pi_per_gene[gene_name] = gene_pi_matrix
        avg_pi_per_gene[gene_name] = gene_avg_pi_matrix
        num_people_with_data[gene_name]=sum(numpy.diagonal(passed_sites_per_gene[gene_name])>min_passed_sites_per_gene)
#    # even though some genes are categorized as variable, they have data for all people. The core/var assignment came from the gene file, but the coverage/bp comes from the SNP file. 
#    if sum(numpy.diagonal(passed_sites_per_gene[gene_name])>5)>100:
#        print gene_name


pi_per_pathway_core, avg_pi_per_pathway_core,passed_sites_per_pathway_core, num_people_with_data_pathway_core, num_genes_per_pathway_core = diversity_utils.calculate_mean_pi_matrix_per_pathway(pi_per_gene, avg_pi_per_gene, passed_sites_per_gene,num_people_with_data,kegg_ids)


############################################################
# compute pi/pathway -- variable genes
############################################################
avg_pi_per_gene={}
pi_per_gene={}
passed_sites_per_gene={}
num_people_with_data={}
for gene_name in variable_genes:
    gene_pi_matrix, gene_avg_pi_matrix, gene_passed_sites= diversity_utils.calculate_pi_matrix(allele_counts_map, passed_sites_map, variant_type='4D', allowed_genes=[gene_name])
    # check if the values of gene_passed_sites are less than 5. If so, then zero out these idxs for gene_passed_sites, gene_pi_matrix, and gene_avg_pi_matrix. Basically all of these people or pairs of people have too few sites to compute realiable statistics. 
    low_passed_sites_idxs=(gene_passed_sites)<min_passed_sites_per_gene
    gene_passed_sites[low_passed_sites_idxs] =0
    gene_pi_matrix[low_passed_sites_idxs]=0
    gene_avg_pi_matrix[low_passed_sites_idxs]=0
    # put all the data into the dictionary for later aggregation. 
    # don't add any genes that have zero information.
    if sum(numpy.diagonal(gene_passed_sites)>min_passed_sites_per_gene)>0: 
        passed_sites_per_gene[gene_name]=gene_passed_sites
        pi_per_gene[gene_name] = gene_pi_matrix
        avg_pi_per_gene[gene_name] = gene_avg_pi_matrix
        num_people_with_data[gene_name]=sum(numpy.diagonal(passed_sites_per_gene[gene_name])>min_passed_sites_per_gene)
#    # even though some genes are categorized as variable, they have data for all people. The core/var assignment came from the gene file, but the coverage/bp comes from the SNP file. 
#    if sum(numpy.diagonal(passed_sites_per_gene[gene_name])>5)>100:
#        print gene_name


pi_per_pathway_variable, avg_pi_per_pathway_variable, passed_sites_per_pathway_variable,num_people_with_data_pathway_variable, num_genes_per_pathway_variable = diversity_utils.calculate_mean_pi_matrix_per_pathway(pi_per_gene, avg_pi_per_gene, passed_sites_per_gene,num_people_with_data,kegg_ids)



############################################################
# compute pi/pathway -- variable AND core genes
############################################################
avg_pi_per_gene={}
pi_per_gene={}
passed_sites_per_gene={}
num_people_with_data={}
for gene_name in gene_names: #gene_names has both core and variable genes. 
    gene_pi_matrix, gene_avg_pi_matrix, gene_passed_sites= diversity_utils.calculate_pi_matrix(allele_counts_map, passed_sites_map, variant_type='4D', allowed_genes=[gene_name])
    # check if the values of gene_passed_sites are less than 5. If so, then zero out these idxs for gene_passed_sites, gene_pi_matrix, and gene_avg_pi_matrix. Basically all of these people or pairs of people have too few sites to compute realiable statistics. 
    low_passed_sites_idxs=(gene_passed_sites)<min_passed_sites_per_gene
    gene_passed_sites[low_passed_sites_idxs] =0
    gene_pi_matrix[low_passed_sites_idxs]=0
    gene_avg_pi_matrix[low_passed_sites_idxs]=0
    # put all the data into the dictionary for later aggregation. 
    # don't add any genes that have zero information.
    if sum(numpy.diagonal(gene_passed_sites)>min_passed_sites_per_gene)>0: 
        passed_sites_per_gene[gene_name]=gene_passed_sites
        pi_per_gene[gene_name] = gene_pi_matrix
        avg_pi_per_gene[gene_name] = gene_avg_pi_matrix
        num_people_with_data[gene_name]=sum(numpy.diagonal(passed_sites_per_gene[gene_name])>min_passed_sites_per_gene)
#    # even though some genes are categorized as variable, they have data for all people. The core/var assignment came from the gene file, but the coverage/bp comes from the SNP file. 
#    if sum(numpy.diagonal(passed_sites_per_gene[gene_name])>5)>100:
#        print gene_name


pi_per_pathway_core_variable, avg_pi_per_pathway_core_variable, passed_sites_per_pathway_core_variable,num_people_with_data_pathway_core_variable, num_genes_per_pathway_core_variable = diversity_utils.calculate_mean_pi_matrix_per_pathway(pi_per_gene, avg_pi_per_gene, passed_sites_per_gene,num_people_with_data,kegg_ids)


##########################################################
# compute total fixations, genome-wide core genes
##########################################################
# Calculate fixation matrices
min_change=0.8
sys.stderr.write("Calculating 4D fixation matrix...\n")
fixation_matrix_syn, fixation_opportunities_syn = diversity_utils.calculate_fixation_matrix(allele_counts_map, passed_sites_map, allowed_variant_types=set(['4D']),allowed_genes=core_genes, min_change=min_change)
sys.stderr.write("Calculating 1D fixation matrix...\n")
fixation_matrix_non, fixation_opportunities_non = diversity_utils.calculate_fixation_matrix(allele_counts_map, passed_sites_map, allowed_variant_types=set(['1D']),allowed_genes=core_genes, min_change=min_change)
sys.stderr.write("Calculating total fixation matrix...\n")
fixation_matrix_all, fixation_opportunities_all = diversity_utils.calculate_fixation_matrix(allele_counts_map, passed_sites_map,allowed_genes=core_genes, min_change=min_change)

sys.stderr.write("Done!\n")

# Calculate fraction nonsynonymous  
dN = fixation_matrix_non/fixation_opportunities_non
dS = fixation_matrix_syn/fixation_opportunities_syn
dNplusdS = (dN+dS)
fraction_nonsynonymous_core = dN/(dNplusdS+(dNplusdS==0))

# Calculate total divergence
dtot_core = fixation_matrix_all/fixation_opportunities_all



##########################################################
# compute total fixations, genome-wide variable genes
##########################################################

# Calculate fixation matrices
min_change=0.8
sys.stderr.write("Calculating 4D fixation matrix...\n")
fixation_matrix_syn, fixation_opportunities_syn = diversity_utils.calculate_fixation_matrix(allele_counts_map, passed_sites_map, allowed_variant_types=set(['4D']),allowed_genes=variable_genes, min_change=min_change)
sys.stderr.write("Calculating 1D fixation matrix...\n")
fixation_matrix_non, fixation_opportunities_non = diversity_utils.calculate_fixation_matrix(allele_counts_map, passed_sites_map, allowed_variant_types=set(['1D']),allowed_genes=variable_genes, min_change=min_change)
sys.stderr.write("Calculating total fixation matrix...\n")
fixation_matrix_all, fixation_opportunities_all = diversity_utils.calculate_fixation_matrix(allele_counts_map, passed_sites_map,allowed_genes=variable_genes, min_change=min_change)
sys.stderr.write("Done!\n")

# Calculate fraction nonsynonymous  
dN = fixation_matrix_non/fixation_opportunities_non
dS = fixation_matrix_syn/fixation_opportunities_syn
dNplusdS = (dN+dS)
fraction_nonsynonymous_variable = dN/(dNplusdS+(dNplusdS==0))

# Calculate total divergence
dtot_variable = fixation_matrix_all/fixation_opportunities_all




###########################################################
# Compute fixations/pathway core genes
########################################################### 

min_change=0.8
fixation_matrix_per_gene_syn={}
fixation_opportunities_per_gene_syn={}
num_people_with_data_syn={}
fixation_matrix_per_gene_non={}
fixation_opportunities_per_gene_non={}
num_people_with_data_non={}
fixation_matrix_per_gene_all={}
fixation_opportunities_per_gene_all={}
num_people_with_data_all={}

for gene_name in core_genes:
    sys.stderr.write("Calculating 4D fixation matrix...\n")
    fixation_matrix, fixation_opportunities= diversity_utils.calculate_fixation_matrix(allele_counts_map, passed_sites_map, allowed_variant_types=set(['4D']),allowed_genes=set([gene_name]), min_change=min_change)
    # check if the values of gene_passed_sites are less than 5. If so, then zero out these idxs for gene_passed_sites, gene_pi_matrix, and gene_avg_pi_matrix. Basically all of these people or pairs of people have too few sites to compute realiable statistics. 
    low_passed_sites_idxs=(fixation_opportunities)<min_passed_sites_per_gene
    fixation_opportunities[low_passed_sites_idxs] =0
    fixation_matrix[low_passed_sites_idxs]=0
    # don't add any genes that have zero information. 
    if sum(numpy.diagonal(fixation_opportunities)>min_passed_sites_per_gene)>0:
        fixation_matrix_per_gene_syn[gene_name]=fixation_matrix
        fixation_opportunities_per_gene_syn[gene_name]=fixation_opportunities
        num_people_with_data_syn[gene_name]=sum(numpy.diagonal(fixation_opportunities)>min_passed_sites_per_gene)
    sys.stderr.write("Calculating 1D fixation matrix...\n")
    fixation_matrix, fixation_opportunities = diversity_utils.calculate_fixation_matrix(allele_counts_map, passed_sites_map, allowed_variant_types=set(['1D']),allowed_genes=set([gene_name]), min_change=min_change)
    low_passed_sites_idxs=(fixation_opportunities)<min_passed_sites_per_gene
    fixation_opportunities[low_passed_sites_idxs] =0
    fixation_matrix[low_passed_sites_idxs]=0
    # don't add any genes that have zero information. 
    if sum(numpy.diagonal(fixation_opportunities)>min_passed_sites_per_gene)>0:
        fixation_matrix_per_gene_non[gene_name]=fixation_matrix
        fixation_opportunities_per_gene_non[gene_name]=fixation_opportunities
        num_people_with_data_non[gene_name]=sum(numpy.diagonal(fixation_opportunities)>min_passed_sites_per_gene)
    sys.stderr.write("Calculating total fixation matrix...\n")
    fixation_matrix, fixation_opportunities = diversity_utils.calculate_fixation_matrix(allele_counts_map, passed_sites_map,allowed_genes=set([gene_name]), min_change=min_change)
    low_passed_sites_idxs=(fixation_opportunities)<min_passed_sites_per_gene
    fixation_opportunities[low_passed_sites_idxs] =0
    fixation_matrix[low_passed_sites_idxs]=0
    # don't add any genes that have zero information. 
    if sum(numpy.diagonal(fixation_opportunities)>min_passed_sites_per_gene)>0:
        fixation_matrix_per_gene_all[gene_name]=fixation_matrix
        fixation_opportunities_per_gene_all[gene_name]=fixation_opportunities
        num_people_with_data_all[gene_name]=sum(numpy.diagonal(fixation_opportunities)>min_passed_sites_per_gene)

sys.stderr.write("Done!\n")

# Aggregate by pathway
dS_per_pathway,fixation_opportunities_per_pathway_syn, num_people_with_data_per_pathway_syn, num_genes_per_pathway_syn=diversity_utils.calculate_mean_fixation_matrix_per_pathway(fixation_matrix_per_gene_syn,fixation_opportunities_per_gene_syn, num_people_with_data_syn,kegg_ids)
dN_per_pathway,fixation_opportunities_per_pathway_non, num_people_with_data_per_pathway_non, num_genes_per_pathway_non=diversity_utils.calculate_mean_fixation_matrix_per_pathway(fixation_matrix_per_gene_non,fixation_opportunities_per_gene_non,num_people_with_data_non, kegg_ids)
dtot_per_pathway_core,fixation_opportunities_per_pathway_all_core, num_people_with_data_per_pathway_tot_core, num_genes_per_pathway_tot_core =diversity_utils.calculate_mean_fixation_matrix_per_pathway(fixation_matrix_per_gene_all,fixation_opportunities_per_gene_all,num_people_with_data_all, kegg_ids)

#calculate fraction nonsynonymous
fraction_nonsynonymous_per_pathway_core={}
num_genes_per_pathway_syn_non_core={}
num_people_with_data_per_pathway_fixations_core={}
fixation_opportunities_per_pathway_syn_non_core={}
for pathway in dS_per_pathway.keys():
    if pathway in dN_per_pathway.keys():
        #check if both dN and dS have enough sites before adding in the data
        syn_non_idx=(fixation_opportunities_per_pathway_syn[pathway]==0)*(fixation_opportunities_per_pathway_non[pathway]==0)
        dS_per_pathway[pathway][syn_non_idx]=0
        dN_per_pathway[pathway][syn_non_idx]=0
        fixation_opportunities_per_pathway_syn[pathway][syn_non_idx]=0
        fixation_opportunities_per_pathway_non[pathway][syn_non_idx]=0
        dNplusdS=dS_per_pathway[pathway] + dN_per_pathway[pathway]
        fraction_nonsynonymous_per_pathway_core[pathway]=dN_per_pathway[pathway]/(dNplusdS+(dNplusdS==0))
        fixation_opportunities_per_pathway_syn_non_core[pathway]=fixation_opportunities_per_pathway_syn[pathway] + fixation_opportunities_per_pathway_non[pathway]
        num_people_with_data_per_pathway_fixations_core[pathway]=(num_people_with_data_per_pathway_syn[pathway]+ num_people_with_data_per_pathway_non[pathway])/2.0
        num_genes_per_pathway_syn_non_core[pathway]=(num_genes_per_pathway_syn[pathway] + num_genes_per_pathway_non[pathway])/2.0



######################################################
# Calculate fixation matrices/pathway variable genes #
######################################################
min_change=0.8
fixation_matrix_per_gene_syn={}
fixation_opportunities_per_gene_syn={}
num_people_with_data_syn={}
fixation_matrix_per_gene_non={}
fixation_opportunities_per_gene_non={}
num_people_with_data_non={}
fixation_matrix_per_gene_all={}
fixation_opportunities_per_gene_all={}
num_people_with_data_all={}

for gene_name in variable_genes:
    sys.stderr.write("Calculating 4D fixation matrix...\n")
    fixation_matrix, fixation_opportunities= diversity_utils.calculate_fixation_matrix(allele_counts_map, passed_sites_map, allowed_variant_types=set(['4D']),allowed_genes=set([gene_name]), min_change=min_change)
    # check if the values of gene_passed_sites are less than 5. If so, then zero out these idxs for gene_passed_sites, gene_pi_matrix, and gene_avg_pi_matrix. Basically all of these people or pairs of people have too few sites to compute realiable statistics. 
    low_passed_sites_idxs=(fixation_opportunities)<min_passed_sites_per_gene
    fixation_opportunities[low_passed_sites_idxs] =0
    fixation_matrix[low_passed_sites_idxs]=0
    # don't add any genes that have zero information. 
    if sum(numpy.diagonal(fixation_opportunities)>min_passed_sites_per_gene)>0:
        fixation_matrix_per_gene_syn[gene_name]=fixation_matrix
        fixation_opportunities_per_gene_syn[gene_name]=fixation_opportunities
        num_people_with_data_syn[gene_name]=sum(numpy.diagonal(fixation_opportunities)>min_passed_sites_per_gene)
    sys.stderr.write("Calculating 1D fixation matrix...\n")
    fixation_matrix, fixation_opportunities = diversity_utils.calculate_fixation_matrix(allele_counts_map, passed_sites_map, allowed_variant_types=set(['1D']),allowed_genes=set([gene_name]), min_change=min_change)
    low_passed_sites_idxs=(fixation_opportunities)<min_passed_sites_per_gene
    fixation_opportunities[low_passed_sites_idxs] =0
    fixation_matrix[low_passed_sites_idxs]=0
    # don't add any genes that have zero information. 
    if sum(numpy.diagonal(fixation_opportunities)>min_passed_sites_per_gene)>0:
        fixation_matrix_per_gene_non[gene_name]=fixation_matrix
        fixation_opportunities_per_gene_non[gene_name]=fixation_opportunities
        num_people_with_data_non[gene_name]=sum(numpy.diagonal(fixation_opportunities)>min_passed_sites_per_gene)
    sys.stderr.write("Calculating total fixation matrix...\n")
    fixation_matrix, fixation_opportunities = diversity_utils.calculate_fixation_matrix(allele_counts_map, passed_sites_map,allowed_genes=set([gene_name]), min_change=min_change)
    low_passed_sites_idxs=(fixation_opportunities)<min_passed_sites_per_gene
    fixation_opportunities[low_passed_sites_idxs] =0
    fixation_matrix[low_passed_sites_idxs]=0
    # don't add any genes that have zero information. 
    if sum(numpy.diagonal(fixation_opportunities)>min_passed_sites_per_gene)>0:
        fixation_matrix_per_gene_all[gene_name]=fixation_matrix
        fixation_opportunities_per_gene_all[gene_name]=fixation_opportunities
        num_people_with_data_all[gene_name]=sum(numpy.diagonal(fixation_opportunities)>min_passed_sites_per_gene)

sys.stderr.write("Done!\n")

# Aggregate by pathway
dS_per_pathway, fixation_opportunities_per_pathway_syn, num_people_with_data_per_pathway_syn, num_genes_per_pathway_syn=diversity_utils.calculate_mean_fixation_matrix_per_pathway(fixation_matrix_per_gene_syn,fixation_opportunities_per_gene_syn, num_people_with_data_syn,kegg_ids)
dN_per_pathway, fixation_opportunities_per_pathway_non, num_people_with_data_per_pathway_non, num_genes_per_pathway_non=diversity_utils.calculate_mean_fixation_matrix_per_pathway(fixation_matrix_per_gene_non,fixation_opportunities_per_gene_non,num_people_with_data_non, kegg_ids)
dtot_per_pathway_variable, fixation_opportunities_per_pathway_all_variable, num_people_with_data_per_pathway_tot_variable, num_genes_per_pathway_tot_variable =diversity_utils.calculate_mean_fixation_matrix_per_pathway(fixation_matrix_per_gene_all,fixation_opportunities_per_gene_all,num_people_with_data_all, kegg_ids)

#calculate fraction nonsynonymous
fraction_nonsynonymous_per_pathway_variable={}
fixation_opportunities_per_pathway_syn_non_variable={}
num_genes_per_pathway_syn_non_variable={}
num_people_with_data_per_pathway_fixations_variable={}
for pathway in dS_per_pathway.keys():
    if pathway in dN_per_pathway.keys():
        #check if both dN and dS have enough sites before adding in the data
        syn_non_idx=(fixation_opportunities_per_pathway_syn[pathway]==0)*(fixation_opportunities_per_pathway_non[pathway]==0)
        dS_per_pathway[pathway][syn_non_idx]=0
        dN_per_pathway[pathway][syn_non_idx]=0
        fixation_opportunities_per_pathway_syn[pathway][syn_non_idx]=0
        fixation_opportunities_per_pathway_non[pathway][syn_non_idx]=0
        dNplusdS=dS_per_pathway[pathway] + dN_per_pathway[pathway]
        fraction_nonsynonymous_per_pathway_variable[pathway]=dN_per_pathway[pathway]/(dNplusdS+(dNplusdS==0))
        fixation_opportunities_per_pathway_syn_non_variable[pathway]=fixation_opportunities_per_pathway_syn[pathway] + fixation_opportunities_per_pathway_non[pathway]
        num_people_with_data_per_pathway_fixations_variable[pathway]=(num_people_with_data_per_pathway_syn[pathway]+ num_people_with_data_per_pathway_non[pathway])/2.0
        num_genes_per_pathway_syn_non_variable[pathway]=(num_genes_per_pathway_syn[pathway] + num_genes_per_pathway_non[pathway])/2.0




###############################################################
# Calculate fixation matrices/pathway core AND variable genes #
###############################################################
min_change=0.8
fixation_matrix_per_gene_syn={}
fixation_opportunities_per_gene_syn={}
num_people_with_data_syn={}
fixation_matrix_per_gene_non={}
fixation_opportunities_per_gene_non={}
num_people_with_data_non={}
fixation_matrix_per_gene_all={}
fixation_opportunities_per_gene_all={}
num_people_with_data_all={}

for gene_name in gene_names:
    sys.stderr.write("Calculating 4D fixation matrix...\n")
    fixation_matrix, fixation_opportunities= diversity_utils.calculate_fixation_matrix(allele_counts_map, passed_sites_map, allowed_variant_types=set(['4D']),allowed_genes=set([gene_name]), min_change=min_change)
    # check if the values of gene_passed_sites are less than 5. If so, then zero out these idxs for gene_passed_sites, gene_pi_matrix, and gene_avg_pi_matrix. Basically all of these people or pairs of people have too few sites to compute realiable statistics. 
    low_passed_sites_idxs=(fixation_opportunities)<min_passed_sites_per_gene
    fixation_opportunities[low_passed_sites_idxs] =0
    fixation_matrix[low_passed_sites_idxs]=0
    # don't add any genes that have zero information. 
    if sum(numpy.diagonal(fixation_opportunities)>min_passed_sites_per_gene)>0:
        fixation_matrix_per_gene_syn[gene_name]=fixation_matrix
        fixation_opportunities_per_gene_syn[gene_name]=fixation_opportunities
        num_people_with_data_syn[gene_name]=sum(numpy.diagonal(fixation_opportunities)>min_passed_sites_per_gene)
    sys.stderr.write("Calculating 1D fixation matrix...\n")
    fixation_matrix, fixation_opportunities = diversity_utils.calculate_fixation_matrix(allele_counts_map, passed_sites_map, allowed_variant_types=set(['1D']),allowed_genes=set([gene_name]), min_change=min_change)
    low_passed_sites_idxs=(fixation_opportunities)<min_passed_sites_per_gene
    fixation_opportunities[low_passed_sites_idxs] =0
    fixation_matrix[low_passed_sites_idxs]=0
    # don't add any genes that have zero information. 
    if sum(numpy.diagonal(fixation_opportunities)>min_passed_sites_per_gene)>0:
        fixation_matrix_per_gene_non[gene_name]=fixation_matrix
        fixation_opportunities_per_gene_non[gene_name]=fixation_opportunities
        num_people_with_data_non[gene_name]=sum(numpy.diagonal(fixation_opportunities)>min_passed_sites_per_gene)
    sys.stderr.write("Calculating total fixation matrix...\n")
    fixation_matrix, fixation_opportunities = diversity_utils.calculate_fixation_matrix(allele_counts_map, passed_sites_map,allowed_genes=set([gene_name]), min_change=min_change)
    low_passed_sites_idxs=(fixation_opportunities)<min_passed_sites_per_gene
    fixation_opportunities[low_passed_sites_idxs] =0
    fixation_matrix[low_passed_sites_idxs]=0
    # don't add any genes that have zero information. 
    if sum(numpy.diagonal(fixation_opportunities)>min_passed_sites_per_gene)>0:
        fixation_matrix_per_gene_all[gene_name]=fixation_matrix
        fixation_opportunities_per_gene_all[gene_name]=fixation_opportunities
        num_people_with_data_all[gene_name]=sum(numpy.diagonal(fixation_opportunities)>min_passed_sites_per_gene)

sys.stderr.write("Done!\n")

# Aggregate by pathway
dS_per_pathway, fixation_opportunities_per_pathway_syn, num_people_with_data_per_pathway_syn, num_genes_per_pathway_syn=diversity_utils.calculate_mean_fixation_matrix_per_pathway(fixation_matrix_per_gene_syn,fixation_opportunities_per_gene_syn, num_people_with_data_syn,kegg_ids)
dN_per_pathway, fixation_opportunities_per_pathway_non, num_people_with_data_per_pathway_non, num_genes_per_pathway_non=diversity_utils.calculate_mean_fixation_matrix_per_pathway(fixation_matrix_per_gene_non,fixation_opportunities_per_gene_non,num_people_with_data_non, kegg_ids)
dtot_per_pathway_core_variable, fixation_opportunities_per_pathway_all_core_variable, num_people_with_data_per_pathway_tot_core_variable, num_genes_per_pathway_tot_core_variable =diversity_utils.calculate_mean_fixation_matrix_per_pathway(fixation_matrix_per_gene_all,fixation_opportunities_per_gene_all,num_people_with_data_all, kegg_ids)

#calculate fraction nonsynonymous
fraction_nonsynonymous_per_pathway_core_variable={}
fixation_opportunities_per_pathway_syn_non_core_variable={}
num_genes_per_pathway_syn_non_core_variable={}
num_people_with_data_per_pathway_fixations_core_variable={}
for pathway in dS_per_pathway.keys():
    if pathway in dN_per_pathway.keys():
        #check if both dN and dS have enough sites before adding in the data
        syn_non_idx=(fixation_opportunities_per_pathway_syn[pathway]==0)*(fixation_opportunities_per_pathway_non[pathway]==0)
        dS_per_pathway[pathway][syn_non_idx]=0
        dN_per_pathway[pathway][syn_non_idx]=0
        fixation_opportunities_per_pathway_syn[pathway][syn_non_idx]=0
        fixation_opportunities_per_pathway_non[pathway][syn_non_idx]=0
        dNplusdS=dS_per_pathway[pathway] + dN_per_pathway[pathway]
        fraction_nonsynonymous_per_pathway_core_variable[pathway]=dN_per_pathway[pathway]/(dNplusdS+(dNplusdS==0))
        fixation_opportunities_per_pathway_syn_non_core_variable[pathway]=fixation_opportunities_per_pathway_syn[pathway] + fixation_opportunities_per_pathway_non[pathway]
        num_people_with_data_per_pathway_fixations_core_variable[pathway]=(num_people_with_data_per_pathway_syn[pathway]+ num_people_with_data_per_pathway_non[pathway])/2.0
        num_genes_per_pathway_syn_non_core_variable[pathway]=(num_genes_per_pathway_syn[pathway] + num_genes_per_pathway_non[pathway])/2.0



###############################################################
# save variables so that I can plot cross-species comparisons #
###############################################################


saved_data=dict(
    avg_pi_matrix_core=avg_pi_matrix_core, 
    avg_pi_matrix_variable=avg_pi_matrix_variable, 
    avg_pi_per_pathway_core=avg_pi_per_pathway_core,
    avg_pi_per_pathway_variable=avg_pi_per_pathway_variable, 
    avg_pi_per_pathway_core_variable=avg_pi_per_pathway_core_variable, 
    passed_sites_per_pathway_core=passed_sites_per_pathway_core,
    passed_sites_per_pathway_variable=passed_sites_per_pathway_variable, 
    passed_sites_per_pathway_core_variable=passed_sites_per_pathway_core_variable, 
    num_genes_per_pathway_core=num_genes_per_pathway_core,
    num_genes_per_pathway_variable=num_genes_per_pathway_variable,
    num_genes_per_pathway_core_variable=num_genes_per_pathway_core_variable,
    num_people_with_data_pathway_core=num_people_with_data_pathway_core,
    num_people_with_data_pathway_variable=num_people_with_data_pathway_variable, 
    num_people_with_data_pathway_core_variable=num_people_with_data_pathway_core_variable, 
    fraction_nonsynonymous_core=fraction_nonsynonymous_core, 
    fraction_nonsynonymous_variable=fraction_nonsynonymous_variable, 
    fraction_nonsynonymous_per_pathway_core=fraction_nonsynonymous_per_pathway_core,
    fraction_nonsynonymous_per_pathway_variable=fraction_nonsynonymous_per_pathway_variable,
    fraction_nonsynonymous_per_pathway_core_variable=fraction_nonsynonymous_per_pathway_core_variable,
    fixation_opportunities_per_pathway_syn_non_core=fixation_opportunities_per_pathway_syn_non_core,
    fixation_opportunities_per_pathway_syn_non_variable=fixation_opportunities_per_pathway_syn_non_variable,
    fixation_opportunities_per_pathway_syn_non_core_variable=fixation_opportunities_per_pathway_syn_non_core_variable,
    num_genes_per_pathway_syn_non_core=num_genes_per_pathway_syn_non_core,
    num_genes_per_pathway_syn_non_variable=num_genes_per_pathway_syn_non_variable,
    num_genes_per_pathway_syn_non_core_variable=num_genes_per_pathway_syn_non_core_variable,
    num_people_with_data_per_pathway_fixations_core=num_people_with_data_per_pathway_fixations_core,
    num_people_with_data_per_pathway_fixations_variable=num_people_with_data_per_pathway_fixations_variable,
    num_people_with_data_per_pathway_fixations_core_variable=num_people_with_data_per_pathway_fixations_core_variable,
    dtot_core=dtot_core, 
    dtot_variable=dtot_variable,
    dtot_per_pathway_core=dtot_per_pathway_core,
    dtot_per_pathway_variable=dtot_per_pathway_variable,
    dtot_per_pathway_core_variable=dtot_per_pathway_core_variable,
    fixation_opportunities_per_pathway_all_core=fixation_opportunities_per_pathway_all_core,
    fixation_opportunities_per_pathway_all_variable=fixation_opportunities_per_pathway_all_variable,
    fixation_opportunities_per_pathway_all_corevariable=fixation_opportunities_per_pathway_all_core_variable,
    num_genes_per_pathway_tot_core=num_genes_per_pathway_tot_core,
    num_genes_per_pathway_tot_variable=num_genes_per_pathway_tot_variable,
    num_genes_per_pathway_tot_core_variable=num_genes_per_pathway_tot_core_variable,
    num_people_with_data_per_pathway_tot_core=num_people_with_data_per_pathway_tot_core,
    num_people_with_data_per_pathway_tot_variable=num_people_with_data_per_pathway_tot_variable,
    num_people_with_data_per_pathway_tot_core_variable=num_people_with_data_per_pathway_tot_core_variable,
    diff_subject_idxs=diff_subject_idxs, 
    same_sample_idxs=same_sample_idxs)

saved_data_file=os.path.expanduser('~/tmp_intermediate_files/kegg_pi_%s.dat' % species_name)
with open(saved_data_file, 'wb') as outfile:
    pickle.dump(saved_data, outfile, protocol=pickle.HIGHEST_PROTOCOL)

############################################################
# Plot
############################################################
###########################################
#### compare all core vs variable genes ###
###########################################

####################################################
# plot total divergence for core vs variable genes #
####################################################

pylab.figure(figsize=(6,2))
pylab.xlabel('Fixations/bp')
pylab.title(species_name+', ' + str(len(same_sample_idxs[0])) + ' subjects')
pylab.xlim(1e-7,1)

data=[]
labels=[]
# first add all data from entire genome
data.append(dtot_core[diff_subject_idxs])
labels.append('All core genes')
data.append(dtot_variable[diff_subject_idxs]) 
labels.append('All variable genes')
 
pylab.boxplot(data,0,'.',0, widths=0.75)
pylab.xscale('log')
locs, dummy_labels = pylab.yticks()
pylab.yticks(locs, labels, fontsize=9)

pylab.savefig('%s/%s_core_vs_variable_genes_fixations.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight',dpi=300)

###################################################################### 
# plot fraction nonsynonymous differences for core vs variable genes
###################################################################### 

pylab.figure(figsize=(6,2))
pylab.xlabel('Fraction nonsynonymous fixations')
pylab.title(species_name+', ' + str(len(same_sample_idxs[0])) + ' subjects')
pylab.xlim(0,1)

data=[]
labels=[]
# first add all data from entire genome
data.append(fraction_nonsynonymous_core[diff_subject_idxs])
labels.append('All core genes')
data.append(fraction_nonsynonymous_variable[diff_subject_idxs]) 
labels.append('All variable genes')
 
pylab.boxplot(data,0,'.',0, widths=0.75)
locs, dummy_labels = pylab.yticks()
pylab.yticks(locs, labels, fontsize=9)

pylab.savefig('%s/%s_core_vs_variable_genes_fraction_nonsynonymous_fixations.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight',dpi=300)

###########################################################
# plot distribution of pi within patients for core vs var #
###########################################################
colors=['#a1d99b','#c994c7']
pylab.figure(figsize=(6,2))
pylab.xlabel('Pi/bp')
pylab.title(species_name+', ' + str(len(same_sample_idxs[0])) + ' subjects')
pylab.xlim(1e-7,1)

data=[]
labels=[]
# first add all data from entire genome
data.append(avg_pi_matrix_core[same_sample_idxs])
labels.append('All core genes')
data.append(avg_pi_matrix_variable[same_sample_idxs])
labels.append('All variable genes')

bp=pylab.boxplot(data,0,'.',0, widths=0.75,patch_artist=True)
k=0
for patch in bp['boxes']:
    patch.set_facecolor(colors[k%2])
    k+=1
#for box in bp['boxes']:
#    box.set(facecolor='green')

pylab.xscale('log')
locs, dummy_labels = pylab.yticks()
pylab.yticks(locs, labels, fontsize=9)
pylab.axvline(x=0.001, ymin=0, ymax=1, hold=None)

pylab.savefig('%s/%s_core_vs_variable_genes_pi_within.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight',dpi=300)


##########################################
### Plot statistics per pathway        ###
##########################################

############################################################   
# plot distribution of pi within patients per pathway
############################################################   
colors=['#a1d99b','#c994c7']
# sort the pathways by number of genes
sorted_pathways=[]
for pathway in sorted(num_genes_per_pathway_core, key=num_genes_per_pathway_core.get, reverse=True):
    sorted_pathways.append(pathway)

pylab.figure(figsize=(4,15))
pylab.xlabel('$\\pi_s$')
pylab.title(species_name +', ' + str(len(same_sample_idxs[0])) + ' subjects')
pylab.xlim(1e-7,1)

data=[]
labels=[]
# first add all data from entire genome
data.append(avg_pi_matrix_core[same_sample_idxs])
data.append(avg_pi_matrix_variable[same_sample_idxs])
labels.append('All core genes, n=' +str(sum(num_genes_per_pathway_core.values()) - num_genes_per_pathway_core['Annotated pathways']))
labels.append('All variable genes, n=' +str(sum(num_genes_per_pathway_variable.values()) - num_genes_per_pathway_variable['Annotated pathways']))
for pathway in sorted_pathways:
    #check which people have enough data. Otherwise don't add this. ####FIX THIS#####
    high_num_sites_idx=numpy.where(passed_sites_per_pathway_core[pathway][same_sample_idxs]>=min_passed_sites_per_person)
    if len(high_num_sites_idx)>0:
        data.append(numpy.clip(avg_pi_per_pathway_core[pathway][same_sample_idxs][high_num_sites_idx],5*1e-7,1))
        num_genes=num_genes_per_pathway_core[pathway]
        if pathway !='':
            labels.append(pathway + ', n=' + str(num_genes)+ ', m='+str(num_people_with_data_pathway_core[pathway]))
        else: 
            labels.append('Unannotated Pathways' + ', n=' + str(num_genes)+ ', m='+str(num_people_with_data_pathway_core[pathway]))

bp=pylab.boxplot(data,0,'.',0, patch_artist=True)
for patch in bp['boxes']:
    patch.set_facecolor(colors[0]) 

pylab.xscale('log')
locs, dummy_labels = pylab.yticks()
pylab.yticks(locs, labels, fontsize=9)
pylab.axvline(x=0.001, ymin=0, ymax=1, hold=None)

pylab.savefig('%s/%s_pi_within_per_pathway_core.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight',dpi=300)

####################################
# Do the same for variable genes
####################################
colors=['#a1d99b','#c994c7']
# sort the pathways by number of genes
sorted_pathways=[]
for pathway in sorted(num_genes_per_pathway_variable, key=num_genes_per_pathway_variable.get, reverse=True):
    sorted_pathways.append(pathway)

pylab.figure(figsize=(4,15))
pylab.xlabel('$\\pi_s$')
pylab.title(species_name +', ' + str(len(same_sample_idxs[0])) + ' subjects')
pylab.xlim(1e-7,1)

data=[]
labels=[]
# first add all data from entire genome
data.append(avg_pi_matrix_core[same_sample_idxs])
data.append(avg_pi_matrix_variable[same_sample_idxs])
labels.append('All core genes, n=' +str(sum(num_genes_per_pathway_core.values()) - num_genes_per_pathway_core['Annotated pathways']))
labels.append('All variable genes, n=' +str(sum(num_genes_per_pathway_variable.values())- num_genes_per_pathway_variable['Annotated pathways']))

for pathway in sorted_pathways:
    #check which people have enough data. Otherwise don't add this.
    high_num_sites_idx=numpy.where(passed_sites_per_pathway_variable[pathway][same_sample_idxs]>=min_passed_sites_per_person)
    if len(high_num_sites_idx)>0:
        data.append(numpy.clip(avg_pi_per_pathway_variable[pathway][same_sample_idxs][high_num_sites_idx],5*1e-7,1))
        num_genes=num_genes_per_pathway_variable[pathway]
        if pathway !='':
            labels.append(pathway + ', n=' + str(num_genes)+ ', m='+str(num_people_with_data_pathway_variable[pathway]))
        else: 
            labels.append('Unannotated Pathways' + ', n=' + str(num_genes)+ ', m='+str(num_people_with_data_pathway_variable[pathway]))

bp=pylab.boxplot(data,0,'.',0, patch_artist=True)
for patch in bp['boxes']:
    patch.set_facecolor(colors[1]) 

pylab.xscale('log')
locs, dummy_labels = pylab.yticks()
pylab.yticks(locs, labels, fontsize=9)
pylab.axvline(x=0.001, ymin=0, ymax=1, hold=None)

pylab.savefig('%s/%s_pi_within_per_pathway_variable.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight',dpi=300)




#################################################################
# plot distribution of fraction DN within patients per pathway  #
# Core
#################################################################
colors=['#a1d99b','#c994c7']
sorted_pathways=[] 
for pathway in sorted(num_genes_per_pathway_syn_non_core, key=num_genes_per_pathway_syn_non_core.get, reverse=True):
    sorted_pathways.append(pathway)


pylab.figure(figsize=(6,15))
pylab.xlabel('Fraction nonsynonymous')
pylab.title(species_name+', ' + str(len(same_sample_idxs[0])) + ' subjects')

data=[]
labels=[]
# first add all data from entire genome
data.append(fraction_nonsynonymous_core[diff_subject_idxs])
data.append(fraction_nonsynonymous_variable[diff_subject_idxs])

labels.append('All core genes, n=' +str(sum(num_genes_per_pathway_syn_non_core.values())-num_genes_per_pathway_syn_non_core['Annotated pathways']))

labels.append('All variable genes, n=' +str(sum(num_genes_per_pathway_syn_non_variable.values())-num_genes_per_pathway_syn_non_variable['Annotated pathways']))

for pathway in sorted_pathways:
    #check which people have enough data. Otherwise don't add this. 
    high_num_sites_idx=numpy.where(fixation_opportunities_per_pathway_syn_non_core[pathway][diff_subject_idxs]>=2*min_passed_sites_per_person) # I'm multipleying min_passed_sites_per_person by 2 because I've summed syn and nonsyn together. 
    if len(high_num_sites_idx)>0:
        data.append(fraction_nonsynonymous_per_pathway_core[pathway][diff_subject_idxs][high_num_sites_idx])
        num_genes=num_genes_per_pathway_syn_non_core[pathway]
    if pathway !='':
        labels.append(pathway + ', n=' + str(int(num_genes)) + ', m='+str(int(num_people_with_data_per_pathway_fixations_core[pathway])))
    else: 
        labels.append('Unannotated Pathways' + ', n=' + str(int(num_genes)) + ', m='+str(int(num_people_with_data_per_pathway_fixations_core[pathway])))

bp=pylab.boxplot(data,0,'.',0, patch_artist=True)
for patch in bp['boxes']:
    patch.set_facecolor(colors[0]) 

locs, dummy_labels = pylab.yticks()
pylab.yticks(locs, labels, fontsize=9)
pylab.savefig('%s/%s_fraction_nonsyn_per_pathway_core.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight',dpi=300)






#################################################################
# plot distribution of fraction DN within patients per pathway  #
# Variable
#################################################################
colors=['#a1d99b','#c994c7']
sorted_pathways=[] 
for pathway in sorted(num_genes_per_pathway_syn_non_variable, key=num_genes_per_pathway_syn_non_variable.get, reverse=True):
    sorted_pathways.append(pathway)


pylab.figure(figsize=(6,15))
pylab.xlabel('Fraction nonsynonymous')
pylab.title(species_name+', ' + str(len(same_sample_idxs[0])) + ' subjects')

data=[]
labels=[]
# first add all data from entire genome
data.append(fraction_nonsynonymous_core[diff_subject_idxs])
data.append(fraction_nonsynonymous_variable[diff_subject_idxs])

labels.append('All core genes, n=' +str(sum(num_genes_per_pathway_syn_non_core.values())-num_genes_per_pathway_syn_non_core['Annotated pathways']))

labels.append('All variable genes, n=' +str(sum(num_genes_per_pathway_syn_non_variable.values())-num_genes_per_pathway_syn_non_variable['Annotated pathways']))

for pathway in sorted_pathways:
    #check which people have enough data. Otherwise don't add this. 
    high_num_sites_idx=numpy.where(fixation_opportunities_per_pathway_syn_non_variable[pathway][diff_subject_idxs]>=2*min_passed_sites_per_person) # I'm multipleying min_passed_sites_per_person by 2 because I've summed syn and nonsyn together. 
    if len(high_num_sites_idx)>0:
        data.append(fraction_nonsynonymous_per_pathway_variable[pathway][diff_subject_idxs][high_num_sites_idx])
        num_genes=num_genes_per_pathway_syn_non_variable[pathway]
    if pathway !='':
        labels.append(pathway + ', n=' + str(int(num_genes)) + ', m='+str(int(num_people_with_data_per_pathway_fixations_variable[pathway])))
    else: 
        labels.append('Unannotated Pathways' + ', n=' + str(int(num_genes)) + ', m='+str(int(num_people_with_data_per_pathway_fixations_variable[pathway])))

bp=pylab.boxplot(data,0,'.',0, patch_artist=True)
for patch in bp['boxes']:
    patch.set_facecolor(colors[1]) 

locs, dummy_labels = pylab.yticks()
pylab.yticks(locs, labels, fontsize=9)
pylab.savefig('%s/%s_fraction_nonsyn_per_pathway_variable.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight',dpi=300)



######################################
# plot total divergence per pathway  #
# core genes                         #
###################################### 
sorted_pathways=[] 
for pathway in sorted(num_genes_per_pathway_tot_core, key=num_genes_per_pathway_tot_core.get, reverse=True):
    sorted_pathways.append(pathway)

pylab.figure(figsize=(6,15))
pylab.xlabel('Total divergence')
pylab.title(species_name+', ' + str(len(same_sample_idxs[0])) + ' subjects')
pylab.xlim(1e-7,1)

data=[]
labels=[]
# first add all data from entire genome
data.append(dtot_core[diff_subject_idxs])
data.append(dtot_variable[diff_subject_idxs])
labels.append('All core genes, n=' +str(int(sum(num_genes_per_pathway_tot_core.values()) - num_genes_per_pathway_tot_core['Annotated pathways'])))
labels.append('All variable genes, n=' +str(int(sum(num_genes_per_pathway_tot_variable.values()) - num_genes_per_pathway_tot_variable['Annotated pathways'])))

for pathway in sorted_pathways:
    #check which people have enough data. Otherwise don't add this. 
    high_num_sites_idx=numpy.where(fixation_opportunities_per_pathway_all_core[pathway][diff_subject_idxs]>=min_passed_sites_per_person)  
    if len(high_num_sites_idx)>0: 
        data.append(numpy.clip(dtot_per_pathway_core[pathway][diff_subject_idxs][high_num_sites_idx],5*1e-7,1))
        num_genes=num_genes_per_pathway_tot_core[pathway]
    if pathway !='':
        labels.append(pathway + ', n=' + str(int(num_genes))+ ', m='+str(int(num_people_with_data_per_pathway_tot_core[pathway])))
    else: 
        labels.append('Unannotated Pathways' + ', n=' + str(int(num_genes))+ ', m='+str(int(num_people_with_data_per_pathway_tot_core[pathway])))

bp=pylab.boxplot(data,0,'.',0, patch_artist=True)

for patch in bp['boxes']:
    patch.set_facecolor(colors[0]) 

pylab.xscale('log')
locs, dummy_labels = pylab.yticks()
pylab.yticks(locs, labels, fontsize=9)

pylab.savefig('%s/%s_fixations_per_pathway_core.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight',dpi=300)



######################################
# plot total divergence per pathway  #
# variable genes                     #
###################################### 
sorted_pathways=[] 
for pathway in sorted(num_genes_per_pathway_tot_variable, key=num_genes_per_pathway_tot_variable.get, reverse=True):
    sorted_pathways.append(pathway)

pylab.figure(figsize=(6,15))
pylab.xlabel('Total divergence')
pylab.title(species_name+', ' + str(len(same_sample_idxs[0])) + ' subjects')
pylab.xlim(1e-7,1)

data=[]
labels=[]
# first add all data from entire genome
data.append(dtot_core[diff_subject_idxs])
data.append(dtot_variable[diff_subject_idxs])
labels.append('All core genes, n=' +str(int(sum(num_genes_per_pathway_tot_core.values()) - num_genes_per_pathway_tot_core['Annotated pathways'])))
labels.append('All variable genes, n=' +str(int(sum(num_genes_per_pathway_tot_variable.values()) - num_genes_per_pathway_tot_variable['Annotated pathways'])))

for pathway in sorted_pathways:
    #check which people have enough data. Otherwise don't add this. 
    high_num_sites_idx=numpy.where(fixation_opportunities_per_pathway_all_variable[pathway][diff_subject_idxs]>=min_passed_sites_per_person)  
    if len(high_num_sites_idx)>0: 
        data.append(numpy.clip(dtot_per_pathway_variable[pathway][diff_subject_idxs][high_num_sites_idx],5*1e-7,1))
        num_genes=num_genes_per_pathway_tot_variable[pathway]
    if pathway !='':
        labels.append(pathway + ', n=' + str(int(num_genes))+ ', m='+str(int(num_people_with_data_per_pathway_tot_variable[pathway])))
    else: 
        labels.append('Unannotated Pathways' + ', n=' + str(int(num_genes))+ ', m='+str(int(num_people_with_data_per_pathway_tot_variable[pathway])))

bp=pylab.boxplot(data,0,'.',0, patch_artist=True)

for patch in bp['boxes']:
    patch.set_facecolor(colors[0]) 

pylab.xscale('log')
locs, dummy_labels = pylab.yticks()
pylab.yticks(locs, labels, fontsize=9)

pylab.savefig('%s/%s_fixations_per_pathway_variable.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight',dpi=300)


