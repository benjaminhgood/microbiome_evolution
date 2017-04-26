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

same_sample_idxs, same_subject_idxs, diff_subject_idxs = parse_midas_data.calculate_subject_pairs(subject_sample_map, low_pi_snp_samples)

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

# Variable genes: find the difference between gene_names and core_genes
gene_names_tmp=numpy.asarray(gene_names)
core_genes_tmp=numpy.asarray(list(core_genes))
variable_genes=set(numpy.asarray(list(numpy.setdiff1d(gene_names_tmp,core_genes_tmp))))


###############################################
# Load kegg information 
##############################################
# load the kegg information for this species
kegg_ids=parse_midas_data.load_kegg_annotations(gene_names)


############################################################## 
# Distribution of presence/absence of genes accross samples  #
############################################################## 

# iterate through is_good_copynum to identify the number of samples in which a gene shows up
low_pi_gene_prevalences=gene_diversity_utils.calculate_gene_prevalences(low_pi_gene_depth_matrix, low_pi_marker_coverages, min_copynum=0.5)

high_pi_gene_prevalences=gene_diversity_utils.calculate_gene_prevalences(high_pi_gene_depth_matrix, high_pi_marker_coverages, min_copynum=0.5)

# the prevalence given above is only for the subset where prev is at least 1. 
# want the full set, so merge low_pi_gene_names with gene_names and assign prev =0 to those gene names where there is no overlap.

low_pi_gene_prevalences_pangenome=gene_diversity_utils.gene_prevalences_whole_pangenome(gene_names, low_pi_gene_names, low_pi_gene_prevalences)
high_pi_gene_prevalences_pangenome=gene_diversity_utils.gene_prevalences_whole_pangenome(gene_names, high_pi_gene_names, high_pi_gene_prevalences)

# using the gene_names list, figure out what pathways are present at which frequency across samples. 
low_pi_kegg_df, low_pi_pathway_description_list=gene_diversity_utils.kegg_pathways_histogram(kegg_ids, gene_names, low_pi_gene_samples, low_pi_gene_prevalences_pangenome)

high_pi_kegg_df, high_pi_pathway_description_list=gene_diversity_utils.kegg_pathways_histogram(kegg_ids, gene_names, high_pi_gene_samples,high_pi_gene_prevalences_pangenome)


##########################################################
# load SNP info
##########################################################

# note that this loads info for all samples. Later the desired samples are selected out. 
sys.stderr.write("Loading %s...\n" % species_name)
samples, allele_counts_map, passed_sites_map, final_line_number = parse_midas_data.parse_snps(species_name, debug, allowed_samples=low_pi_snp_samples)
sys.stderr.write("Done!\n")


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

pi_per_gene, avg_pi_per_gene, passed_sites_per_gene, num_people_with_data =  diversity_utils.calculate_pi_matrix_per_gene(allele_counts_map, passed_sites_map, variant_type='4D', allowed_genes=core_genes)

# aggregate pi/gene by pathway. 
pi_per_pathway_core, avg_pi_per_pathway_core, passed_sites_per_pathway_core, num_genes_per_pathway_core, num_people_with_data_per_pathway_core = diversity_utils.calculate_mean_pi_matrix_per_pathway(pi_per_gene, avg_pi_per_gene, passed_sites_per_gene,num_people_with_data,kegg_ids)

for pathway_name in avg_pi_per_pathway_core.keys():
    avg_pi_per_pathway_core[pathway_name] = avg_pi_per_pathway_core[pathway_name]/(passed_sites_per_pathway_core[pathway_name]+(passed_sites_per_pathway_core[pathway_name]==0)) 
    pi_per_pathway_core[pathway_name] = pi_per_pathway_core[pathway_name]/(passed_sites_per_pathway_core[pathway_name]+(passed_sites_per_pathway_core[pathway_name]==0))     

############################################################
# compute pi/pathway -- variable genes
############################################################

pi_per_gene, avg_pi_per_gene, passed_sites_per_gene, num_people_with_data =  diversity_utils.calculate_pi_matrix_per_gene(allele_counts_map, passed_sites_map, variant_type='4D', allowed_genes=variable_genes)

# aggregate pi/gene by pathway. 
pi_per_pathway_variable, avg_pi_per_pathway_variable, passed_sites_per_pathway_variable, num_genes_per_pathway_variable, num_people_with_data_per_pathway_variable = diversity_utils.calculate_mean_pi_matrix_per_pathway(pi_per_gene, avg_pi_per_gene, passed_sites_per_gene,num_people_with_data,kegg_ids)

for pathway_name in avg_pi_per_pathway_variable.keys():
    avg_pi_per_pathway_variable[pathway_name] = avg_pi_per_pathway_variable[pathway_name]/(passed_sites_per_pathway_variable[pathway_name]+(passed_sites_per_pathway_variable[pathway_name]==0))     
#debug:
#    pi_per_pathway_variable[pathway_name] = pi_per_pathway_variable[pathway_name]/(passed_sites_per_pathway_variable[pathway_name]+(passed_sites_per_pathway_variable[pathway_name]==0))     



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
fixation_matrix_syn, fixation_opportunities_syn= diversity_utils.calculate_fixation_matrix(allele_counts_map, passed_sites_map, allowed_variant_types=set(['4D']),allowed_genes=variable_genes, min_change=min_change)
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
# Compute fixations/pathway
########################################################### 
sys.stderr.write("Calculating 4D fixation matrix...\n")
fixation_matrix_syn, fixation_opportunities_syn, num_people_with_data_syn = diversity_utils.calculate_fixation_matrix_per_gene(allele_counts_map, passed_sites_map, allowed_variant_types=set(['4D']), allowed_genes=core_genes,min_change=min_change)
sys.stderr.write("Calculating 1D fixation matrix...\n")
fixation_matrix_non, fixation_opportunities_non, num_people_with_data_non = diversity_utils.calculate_fixation_matrix_per_gene(allele_counts_map, passed_sites_map, allowed_variant_types=set(['1D']), allowed_genes=core_genes,min_change=min_change)
sys.stderr.write("Calculating total fixation matrix...\n")
fixation_matrix_all, fixation_opportunities_all, num_people_with_data_all = diversity_utils.calculate_fixation_matrix_per_gene(allele_counts_map, passed_sites_map,allowed_genes=core_genes, min_change=min_change)

sys.stderr.write("Done!\n")



# Aggregate by pathway
dS_per_pathway, num_genes_per_pathway_syn, num_people_with_data_per_pathway_syn=diversity_utils.calculate_mean_fixation_matrix_per_pathway(fixation_matrix_syn,fixation_opportunities_syn, num_people_with_data_syn,kegg_ids)
dN_per_pathway, num_genes_per_pathway_non, num_people_with_data_per_pathway_non=diversity_utils.calculate_mean_fixation_matrix_per_pathway(fixation_matrix_non,fixation_opportunities_non,num_people_with_data_non, kegg_ids)
dtot_per_pathway, num_genes_per_pathway_tot, num_people_with_data_per_pathway_tot=diversity_utils.calculate_mean_fixation_matrix_per_pathway(fixation_matrix_all,fixation_opportunities_all,num_people_with_data_all, kegg_ids)

#calculate fraction nonsynonymous
fraction_nonsynonymous_per_pathway={}
num_genes_per_pathway_syn_non={}
for pathway in dS_per_pathway.keys():
    if pathway in dN_per_pathway.keys():
        dNplusdS=dS_per_pathway[pathway] + dN_per_pathway[pathway]
        fraction_nonsynonymous_per_pathway[pathway]=dN_per_pathway[pathway]/(dNplusdS+(dNplusdS==0))
        num_genes_per_pathway_syn_non[pathway]=(num_genes_per_pathway_syn[pathway] + num_genes_per_pathway_non[pathway])/2.0


############################################################
# Plot
############################################################
######################################
# compare all core vs variable genes
######################################

# plot total divergence for core vs variable genes

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


# plot fraction nonsynonymous differences for core vs variable genes

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


# plot distribution of pi within patients for core vs var
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

pylab.boxplot(data,0,'.',0, widths=0.75)
pylab.xscale('log')
locs, dummy_labels = pylab.yticks()
pylab.yticks(locs, labels, fontsize=9)
pylab.axvline(x=0.001, ymin=0, ymax=1, hold=None)

pylab.savefig('%s/%s_core_vs_variable_genes_pi_within.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight',dpi=300)


##########################################
# Return to this part of the code later  #
##########################################
'''

# plot distribution of pi within patients per pathway

# sort the pathways by number of genes
sorted_pathways=[]
for pathway in sorted(num_genes_per_pathway, key=num_genes_per_pathway.get, reverse=True):
    sorted_pathways.append(pathway)

pylab.figure(figsize=(4,15))
pylab.xlabel('$\\pi_s$')
pylab.title(species_name +', ' + str(len(same_sample_idxs[0])) + ' subjects')
pylab.xlim(1e-7,1)

data=[]
labels=[]
# first add all data from entire genome
data.append(avg_pi_matrix[same_sample_idxs])
#labels.append('All core genes, n=' +str(sum(num_genes_per_pathway.values())))
labels.append('All variable genes, n=' +str('x'))
for pathway in sorted_pathways:
    data.append(avg_pi_per_pathway[pathway][same_sample_idxs])
    #data.append(numpy.clip(avg_pi_per_pathway[pathway][same_sample_idxs],1e-10,1))
    num_genes=num_genes_per_pathway[pathway]
    if pathway !='':
        labels.append(pathway + ', n=' + str(num_genes))
    else: 
        labels.append('Unannotated Pathways' + ', n=' + str(num_genes))

pylab.boxplot(data,0,'.',0)
pylab.xscale('log')
locs, dummy_labels = pylab.yticks()
pylab.yticks(locs, labels, fontsize=9)
pylab.axvline(x=0.001, ymin=0, ymax=1, hold=None)

pylab.savefig('%s/%s_pi_within_per_pathway.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight',dpi=300)





# plot distribution of fraction DN within patients per pathway
sorted_pathways=[] 
for pathway in sorted(num_genes_per_pathway_syn_non, key=num_genes_per_pathway_syn_non.get, reverse=True):
    sorted_pathways.append(pathway)


pylab.figure(figsize=(6,15))
pylab.xlabel('Fraction nonsynonymous')
pylab.title(species_name+', ' + str(len(same_sample_idxs[0])) + ' subjects')

data=[]
labels=[]
# first add all data from entire genome
data.append(fraction_nonsynonymous[diff_subject_idxs])
#labels.append('All core genes, n=' +str(int(num_genes_per_pathway_syn_non['']+num_genes_per_pathway_syn_non['Annotated_pathways'])))
labels.append('All variable genes, n=' +str('x'))

for pathway in sorted_pathways:
    data.append(fraction_nonsynonymous_per_pathway[pathway][diff_subject_idxs])
    num_genes=num_genes_per_pathway_syn_non[pathway]
    if pathway !='':
        labels.append(pathway + ', n=' + str(int(num_genes)))
    else: 
        labels.append('Unannotated Pathways' + ', n=' + str(int(num_genes)))

pylab.boxplot(data,0,'.',0)
locs, dummy_labels = pylab.yticks()
pylab.yticks(locs, labels, fontsize=9)
pylab.savefig('%s/%s_fraction_nonsyn_per_pathway.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight',dpi=300)




# plot total divergence
sorted_pathways=[] 
for pathway in sorted(num_genes_per_pathway_tot, key=num_genes_per_pathway_tot.get, reverse=True):
    sorted_pathways.append(pathway)

pylab.figure(figsize=(6,15))
pylab.xlabel('Total divergence')
pylab.title(species_name+', ' + str(len(same_sample_idxs[0])) + ' subjects')
pylab.xlim(1e-7,1)

data=[]
labels=[]
# first add all data from entire genome
data.append(dtot[diff_subject_idxs])
#labels.append('All core genes, n=' +str(int(sum(num_genes_per_pathway_tot.values()))))
labels.append('All variable genes, n=' +str('x'))
for pathway in sorted_pathways:
    data.append(dtot_per_pathway[pathway][diff_subject_idxs])
    num_genes=num_genes_per_pathway_tot[pathway]
    if pathway !='':
        labels.append(pathway + ', n=' + str(int(num_genes)))
    else: 
        labels.append('Unannotated Pathways' + ', n=' + str(int(num_genes)))

pylab.boxplot(data,0,'.',0)
pylab.xscale('log')
locs, dummy_labels = pylab.yticks()
pylab.yticks(locs, labels, fontsize=9)

pylab.savefig('%s/%s_fixations_per_pathway.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight',dpi=300)

'''
