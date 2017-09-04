# Within-snp gene changes

import matplotlib  
matplotlib.use('Agg') 
import config
import parse_midas_data
import os
import parse_HMP_data


import pylab
import sys
import numpy
import random

import diversity_utils
import gene_diversity_utils
import sfs_utils
import calculate_substitution_rates
import calculate_temporal_changes
import parse_patric

import stats_utils
import matplotlib.colors as colors
import matplotlib.cm as cmx
from math import log10,ceil
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy.random import randint, choice

mpl.rcParams['font.size'] = 5
mpl.rcParams['lines.linewidth'] = 0.5
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'

#species_name = "Bacteroides_vulgatus_57955"

################################################################################
#
# Standard header to read in argument information
#
################################################################################
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--debug", help="Loads only a subset of SNPs for speed", action="store_true")
parser.add_argument("--chunk-size", type=int, help="max number of records to load", default=1000000000)
parser.add_argument('--other-species', type=str, help='Run the script for a different species')

args = parser.parse_args()

debug = args.debug
chunk_size = args.chunk_size
other_species = args.other_species

if other_species:
    species_name = other_species
    other_species_str = "_%s" % species_name
else:
    other_species_str = ""
    
################################################################################


if other_species=="":
    outFile=open('/pollard/home/ngarud/ben_nandita_hmp_analysis/gene_changes_shuffle_all_species.txt', 'w' )
else:
    outFile=open('/pollard/home/ngarud/ben_nandita_hmp_analysis/gene_changes_shuffle_%s.txt' % species_name , 'w' )   

modification_divergence_threshold = 1e-03 #the threshold for deciding when something is a modification vs a replacement. Like most other things, it is an arbitrary choice. 
min_coverage=20
clade_divergence_threshold = 1e-02

if other_species_str == "":
    good_species_list = parse_midas_data.parse_good_species_list()
else:
    good_species_list=[species_name]

# store all the species' data in a dictionary:
all_data={}

#key=species
#value={}, key=gene, valuee=num times gene shows up

for species_name in good_species_list: 
    all_species_gene_changes={}
    all_species_gene_changes_category={}
    all_species_null={}

    all_data[species_name]={}

    ####################
    # Analyze the data #
    ####################
    # find CPS samples:
    snp_samples = diversity_utils.calculate_haploid_samples(species_name, debug=debug)

    # load pre-computed data:
    sys.stderr.write("Loading pre-computed substitution rates for %s...\n" % species_name)
    substitution_rate_map = calculate_substitution_rates.load_substitution_rate_map(species_name)
    sys.stderr.write("Calculating matrix...\n")
    dummy_samples, snp_difference_matrix, snp_opportunity_matrix = calculate_substitution_rates.calculate_matrices_from_substitution_rate_map(substitution_rate_map, 'all', allowed_samples=snp_samples)
    snp_samples = dummy_samples
    sys.stderr.write("Done!\n")

    sys.stderr.write("Loading pre-computed temporal changes for %s...\n" % species_name)
    temporal_change_map = calculate_temporal_changes.load_temporal_change_map(species_name)
    sys.stderr.write("Done!\n")

    snp_substitution_rate = snp_difference_matrix*1.0/(snp_opportunity_matrix+(snp_opportunity_matrix==0))
    sys.stderr.write("Done!\n")   

    # get all genome ids for this species' pan genome:
    genome_ids=parse_midas_data.get_ref_genome_ids(species_name)

    # load the gene ids for all genomes coresponding to this speceis:
    gene_descriptions=parse_patric.load_patric_gene_descriptions(genome_ids)

    # create gene categories (poor proxy for GO terms):
    gene_categories, gene_category_map = parse_patric.cluster_patric_gene_descriptions(gene_descriptions)

    # load the kegg ids for all genomes corresponding to this species:
    kegg_ids=parse_patric.load_kegg_annotations(genome_ids)    
    
    ##################
    # pangenome null #
    ##################

    # load all pangenome genes for the species after clustering at 95% identity
    pangenome_gene_names, pangenome_new_species_names=parse_midas_data.load_pangenome_genes(species_name)
    

    ###########################################
    # load data for between host changes null #
    ###########################################

    # Load gene coverage information for species_name
    sys.stderr.write("Loading pangenome data for %s...\n" % species_name)
    gene_samples, gene_names, gene_presence_matrix, gene_depth_matrix, marker_coverages, gene_reads_matrix = parse_midas_data.parse_pangenome_data(species_name,allowed_samples=snp_samples)
    sys.stderr.write("Done!\n")
    
    # convert gene names to numpy array:
    gene_names=numpy.array(gene_names)

    # indexes for different subject pairs
    desired_samples = gene_samples

    sample_order_map = parse_HMP_data.parse_sample_order_map()
    desired_same_sample_idxs, desired_same_subject_idxs, desired_diff_subject_idxs = parse_midas_data.calculate_ordered_subject_pairs(sample_order_map, desired_samples)

    snp_sample_idx_map = parse_midas_data.calculate_sample_idx_map(desired_samples, snp_samples)
    gene_sample_idx_map = parse_midas_data.calculate_sample_idx_map(desired_samples, gene_samples)

    diff_subject_snp_idxs = parse_midas_data.apply_sample_index_map_to_indices(snp_sample_idx_map, desired_diff_subject_idxs)  
    diff_subject_gene_idxs = parse_midas_data.apply_sample_index_map_to_indices(gene_sample_idx_map, desired_diff_subject_idxs)  

    between_host_gene_idxs = [] # store idxs of genes that change between hosts
    for sample_pair_idx in xrange(0,len(diff_subject_snp_idxs[0])):
        snp_i = diff_subject_snp_idxs[0][sample_pair_idx]
        snp_j = diff_subject_snp_idxs[1][sample_pair_idx]

        i = diff_subject_gene_idxs[0][sample_pair_idx]
        j = diff_subject_gene_idxs[1][sample_pair_idx]
        if (marker_coverages[i]>min_coverage) and (marker_coverages[j]>min_coverage):
            if snp_substitution_rate[snp_i, snp_j] < clade_divergence_threshold:
                gene_idxs = gene_diversity_utils.calculate_gene_differences_between_idxs(i,j, gene_reads_matrix, gene_depth_matrix, marker_coverages)
                between_host_gene_idxs.extend(gene_idxs) # collect all gene changes occurring between hosts. Use this for the null.
    
    # get a unique set of genes that change between hosts:
    between_host_gene_idxs_tmp=between_host_gene_idxs
    between_host_gene_idxs=list(set(between_host_gene_idxs))
    

    #######################
    # within host changes #
    #######################

    # store the actual data in this:
    within_host_changes_gene_ids=[]
    
    # store null data in this to see how the actual data above compares. Draw from between_host_gene_idxs to populate this. 
    between_host_changes_gene_ids_null={} #dictionary which stores different trials (trial=key)
    pangenome_null={} # same except null is pangenome

    within_host_change_distribution=[] # use this for constructing the null
    for sample_pair in temporal_change_map.keys():
        sample_1=sample_pair[0]
        sample_2=sample_pair[1]
        # check if this pair underwent a replacement event. If so, ignore, otherwise stats get messed up. 
        # find the index of sample_1 and sample_2 in snp_samples to see if this is a CPS sample pair:
        include_sample_pair=False
        if sample_1 in snp_samples and sample_2 in snp_samples:
            sample_1_idx=snp_samples.index(sample_1)
            sample_2_idx=snp_samples.index(sample_2)
            if snp_substitution_rate[sample_1_idx, sample_2_idx] < modification_divergence_threshold:
                include_sample_pair=True
        if include_sample_pair == True:
            gene_perr, gains, losses=calculate_temporal_changes.calculate_gains_losses_from_temporal_change_map(temporal_change_map, sample_1, sample_2, lower_threshold=0.05)
            all_changes=gains+losses
            within_host_change_distribution.append(len(all_changes))
            #iterate through all_changes to store the gene_ids.
            for i in range(0, len(all_changes)):
                within_host_changes_gene_ids.append(all_changes[i][0])



    ###################################################################################
    # grab 100x the number of within-host gene changes from between host gene changes #
    # this is our null against we compare annotations to
    # make sure that the sampling is without replacment
    ##################################################################################
    
    num_trials=100
    for trial in range(0,num_trials):
        between_host_changes_gene_ids_null[trial]=[]
        for within_num in within_host_change_distribution:
            between_gene_idxs = choice(between_host_gene_idxs, within_num)
            # between_gene_idxs = choice(between_host_gene_idxs, len(within_host_changes_gene_ids))
            # get the gene_ids for these idxs:
            between_host_changes_gene_ids_null[trial].extend(gene_names[between_gene_idxs].tolist())

    #######################################################################################
    # pangenome null:
    # repeat similar procudure, using pangenome genes instead of between-host gene chagnes #
    #######################################################################################
    for trial in range(0,num_trials): 
        pangenome_null[trial]=[]
        for within_num in within_host_change_distribution:
            #pangenome_null[trial]=random.sample(pangenome_gene_names, len(within_host_changes_gene_ids))
            pangenome_null[trial].extend(random.sample(pangenome_gene_names, within_num))

    #########################################################
    # what is the gene change annotation
    # kegg pathway annotation
    # gene_cateogry annotation
    #########################################################

    gene_descriptions_gene_changes=[] # store the descriptions in this vector
    kegg_pathways_gene_changes=[] # store the pathways in this vector
    gene_categories_gene_changes=[] #stor categories in this vector

    for gene_id in within_host_changes_gene_ids:
        gene_descriptions_gene_changes.append(gene_descriptions[gene_id])
        gene_categories_gene_changes.append(gene_category_map[gene_id])
        for i in range(0, len(kegg_ids[gene_id])):
            kegg_pathways_gene_changes.append(kegg_ids[gene_id][i][1])

    gene_descriptions_gene_changes_null_between_host={}
    kegg_pathways_gene_changes_null_between_host={}
    gene_categories_gene_changes_null_between_host={}

    gene_descriptions_gene_changes_null_pangenome={}
    kegg_pathways_gene_changes_null_pangenome={}
    gene_categories_gene_changes_null_pangenome={}

    for trial in range(0, num_trials):
        gene_descriptions_gene_changes_null_between_host[trial]=[]
        kegg_pathways_gene_changes_null_between_host[trial]=[]
        gene_categories_gene_changes_null_between_host[trial]=[]
        for gene_id in between_host_changes_gene_ids_null[trial]:
            if gene_id in gene_descriptions:
                gene_descriptions_gene_changes_null_between_host[trial].append(gene_descriptions[gene_id])
                gene_categories_gene_changes_null_between_host[trial].append(gene_category_map[gene_id])
            for i in range(0, len(kegg_ids[gene_id])):
                if gene_id in kegg_ids:
                    kegg_pathways_gene_changes_null_between_host[trial].append(kegg_ids[gene_id][i][1])

        gene_descriptions_gene_changes_null_pangenome[trial]=[]
        kegg_pathways_gene_changes_null_pangenome[trial]=[]
        gene_categories_gene_changes_null_pangenome[trial]=[]
        for gene_id in pangenome_null[trial]:
            if gene_id in gene_descriptions:
                gene_descriptions_gene_changes_null_pangenome[trial].append(gene_descriptions[gene_id])
                gene_categories_gene_changes_null_pangenome[trial].append(gene_category_map[gene_id]) 
            for i in range(0, len(kegg_ids[gene_id])):
                if gene_id in kegg_ids:
                    kegg_pathways_gene_changes_null_pangenome[trial].append(kegg_ids[gene_id][i][1])
                
    




    #############################
    # compute empirical p-value #
    #############################

    # how likely is it to see a given gene change occur across and within people as many times as it does given all the nulls across people?

    # count the number of times a gene shows up in within-host changes
    all_gene_changes={}
    for gene in gene_descriptions_gene_changes:
        if gene not in all_gene_changes:
            all_gene_changes[gene]=1
        else:
            all_gene_changes[gene] +=1

    # count the number of times a gene shows up in between-host changes in the num_trails
    all_gene_changes_null_between_host={} # key == gene
    all_gene_changes_null_pangenome={}
    
    for trial in range(0,num_trials):
        for gene in gene_descriptions_gene_changes_null_between_host[trial]:
            all_gene_changes_null_between_host[gene]=[]
        for gene in gene_descriptions_gene_changes_null_pangenome[trial]: 
            all_gene_changes_null_pangenome[gene]=[]

    for gene in all_gene_changes_null_between_host.keys():
        for trial in range(0,num_trials): # iterate through the 100 trials
            num_occurrences=gene_descriptions_gene_changes_null_between_host[trial].count(gene)
            all_gene_changes_null_between_host[gene].append(num_occurrences)

    for gene in all_gene_changes_null_pangenome.keys():
        for trial in range(0,num_trials): # iterate through the 100 trials
            num_occurrences=gene_descriptions_gene_changes_null_pangenome[trial].count(gene)
            all_gene_changes_null_pangenome[gene].append(num_occurrences)



    # iterate through all gene changes and quantify observed vs expected:
    
    outFile.write(species_name + '\ngene_name\tprob_between_hosts\tprob_pangenome\n')
   
    for gene in all_gene_changes.keys():
        if gene in all_gene_changes_null_between_host:
            tmp_null = numpy.array(all_gene_changes_null_between_host[gene])
            prob_between_hosts=sum(all_gene_changes[gene] <= tmp_null)/100.0
        else:
            prob_between_hosts=0
        if gene in all_gene_changes_null_pangenome:
            tmp_null = numpy.array(all_gene_changes_null_pangenome[gene])
            prob_pangenome=sum(all_gene_changes[gene] <= tmp_null)/100.0
        else:
            prob_pangenome=0
        #store data in a dictionary for cross-species plotting after running on all species
        all_species_gene_changes[gene]=[all_gene_changes[gene],prob_between_hosts,prob_pangenome]
        
        outFile.write(gene + '\t' + str(all_gene_changes[gene]) +'\t' + str(prob_between_hosts) + '\t' + str(prob_pangenome) + '\n')
        
    outFile.write('\n')
    ###################################################################
    # check whether there is an enrichment of genes in a gene_category
    ##########################################################


    # count the number of times a gene shows up in within-host changes
    all_gene_changes_category={}
    for gene in gene_categories_gene_changes:
        if gene not in all_gene_changes_category:
            all_gene_changes_category[gene]=1
        else:
            all_gene_changes_category[gene] +=1

    # count the number of times a gene shows up in between-host changes in the num_trails
    all_gene_changes_category_null_between_host={} # key == gene
    all_gene_changes_category_null_pangenome={} # key == gene
    
    
    for trial in range(0,num_trials):
        for gene in gene_categories_gene_changes_null_between_host[trial]:
            all_gene_changes_category_null_between_host[gene]=[]
        for gene in gene_categories_gene_changes_null_pangenome[trial]:
            all_gene_changes_category_null_pangenome[gene]=[]

    for gene in all_gene_changes_category_null_between_host.keys():
        for trial in range(0,num_trials): # iterate through the 100 trials
            num_occurrences=gene_categories_gene_changes_null_between_host[trial].count(gene)
            all_gene_changes_category_null_between_host[gene].append(num_occurrences)
            
    for gene in all_gene_changes_category_null_pangenome.keys():
        for trial in range(0,num_trials): # iterate through the 100 trials
            num_occurrences=gene_categories_gene_changes_null_pangenome[trial].count(gene)
            all_gene_changes_category_null_pangenome[gene].append(num_occurrences)


    # iterate through all gene changes and quantify observed vs expected:
    outFile.write(species_name + '\ngene_name\tprob_between_hosts\tprob_pangenome\n')

    for gene in all_gene_changes_category.keys():
        if gene in all_gene_changes_category_null_between_host:
            tmp_null = numpy.array(all_gene_changes_category_null_between_host[gene])
            prob_between_hosts=sum(all_gene_changes_category[gene] <= tmp_null)/100.0
        else:
            prob_between_hosts=0
        if gene in all_gene_changes_category_null_pangenome:
            tmp_null = numpy.array(all_gene_changes_category_null_pangenome[gene])
            prob_pangenome=sum(all_gene_changes_category[gene] <= tmp_null)/100.0
        else:
            prob_pangenome=0
        all_species_gene_changes_category[gene]=[all_gene_changes_category[gene],prob_between_hosts,prob_pangenome]
        outFile.write(gene + '\t' + str(all_gene_changes_category[gene]) +'\t' + str(prob_between_hosts) + '\t' + str(prob_pangenome) + '\n')

    outFile.write('\n')

    # store the nulls:
    all_species_null={'between_host_genes':all_gene_changes_null_between_host,'pangenome_genes': all_gene_changes_null_pangenome, 'between_host_category': all_gene_changes_category_null_between_host, 'pangenome_category':all_gene_changes_category_null_pangenome}

    #store all the data 
    all_data[species_name]={'gene_changes':all_species_gene_changes, 'gene_changes_category':all_species_gene_changes_category, 'null':all_species_null}

import pickle
if other_species_str=="":
    pickle.dump( all_data, open( "/pollard/home/ngarud/tmp_intermediate_files/all_species_gene_changes.p", "wb" ) )
else:
    pickle.dump( all_data, open( "/pollard/home/ngarud/tmp_intermediate_files/%s_gene_changes.p" % species_name, "wb" ) )  


#######################################
# code for loading cross-species data #
#######################################

import pickle
good_species_list = parse_midas_data.parse_good_species_list() 


all_data={}
for species_name in good_species_list:
    print species_name
    if os.path.isfile("/pollard/home/ngarud/tmp_intermediate_files/%s_gene_changes.p" %species_name):
        all_data_species = pickle.load( open( "/pollard/home/ngarud/tmp_intermediate_files/%s_gene_changes.p" %species_name, "rb" ))
        all_data[species_name]=all_data_species[species_name]


# want to know what the sum of all gene changes is across species
all_data['all_species']={}
all_data['all_species']['gene_changes']={}
for species_name in all_data.keys():
    if species_name != 'all_species':
        for gene in all_data[species_name]['gene_changes'].keys():
            if gene not in all_data['all_species']['gene_changes'].keys():
                all_data['all_species']['gene_changes'][gene] = [all_data[species_name]['gene_changes'][gene][0],0,0]
            else:
                all_data['all_species']['gene_changes'][gene][0] +=all_data[species_name]['gene_changes'][gene][0]


# want to know what the null is combining all trials across speceis:
all_data['all_species']['null']={}
all_data['all_species']['null']['between_host_genes']={}

for species_name in all_data.keys(): 
    if species_name != 'all_species':
        for gene in all_data[species_name]['null']['between_host_genes'].keys():
            if gene not in  all_data['all_species']['null']['between_host_genes'].keys():
                all_data['all_species']['null']['between_host_genes'][gene]=all_data[species_name]['null']['between_host_genes'][gene]
            else:
                for i in range(0, num_trials):
                    all_data['all_species']['null']['between_host_genes'][gene][i] += all_data[species_name]['null']['between_host_genes'][gene][i]
                

# compute empirical p-value for all_species


# put everything in a 2D matrix so that I can plot a heatmap:
num_species=len(all_data.keys())
num_genes=len(all_data['all_species']['gene_changes'].keys())-1 # remove hypothetical protein
gene_occurrence_matrix=numpy.zeros((num_genes, num_species))

# sort genes based on how abundant they are accross species:
gene_order_sorted = sorted(all_data['all_species']['gene_changes'].items(), key=operator.itemgetter(1))
gene_order=[]
for i in range(0, len(gene_order_sorted)):
    if gene_order_sorted[i][0] !='hypothetical protein':
        gene_order.append(gene_order_sorted[i][0])



# fill in the matrix
species_idx=0
for species_name in all_data.keys():
    gene_idx=0
    for gene in gene_order:
        if gene != 'hypothetical protein':
            if gene in all_data[species_name]['gene_changes']:
                gene_occurrence_matrix[gene_idx, species_idx] = all_data[species_name]['gene_changes'][gene][0]
            gene_idx +=1
    species_idx +=1



species_idx=0
for species_name in all_data.keys():
    gene_idx=0
    for gene in gene_order:
        if gene != 'hypothetical protein':
            if gene in all_data[species_name]['gene_changes']:
                if all_data[species_name]['gene_changes'][gene][1] <=0.05:
                    gene_occurrence_matrix[gene_idx, species_idx] = 2
                else:
                    gene_occurrence_matrix[gene_idx, species_idx] = 1
            gene_idx +=1
    species_idx +=1



import seaborn as sns; sns.set()
pylab.figure(figsize=(10,15))   
sns.set(font_scale=0.4)
ax = sns.heatmap(gene_occurrence_matrix, yticklabels=gene_order, xticklabels=all_data.keys()) 
pylab.savefig('%s/gene_change_cross_species_distribution_betwn_host_prob.png' % (parse_midas_data.analysis_directory),bbox_inches='tight',dpi=300)


'''
    # put the results in a 2D matrix so that I can plot a heatmap.
    # size of the heatmap is num_samples x num_unique_pathways
    if len(kegg_pathways_gene_changes):
        first=len(kegg_pathways_gene_changes)
        second=len(kegg_pathways_gene_changes_dict.keys())
        kegg_occurrence_matrix=numpy.zeros((first,second))
        #fill in the matrix:
        sample_pair_idx=0
        for sample_pair in kegg_pathways_gene_changes_dict.keys():
            for pathway in kegg_pathways_gene_changes_dict[sample_pair]:
                pathway_idx=kegg_pathways_gene_changes.index(pathway)
                kegg_occurrence_matrix[pathway_idx,sample_pair_idx]+=1
            sample_pair_idx+=1

        # plot this as a heatmap:
        import seaborn as sns; sns.set()

        pylab.figure(figsize=(6,10)) 
        sns.set(font_scale=0.4) # set the font
        ax = sns.heatmap(kegg_occurrence_matrix,yticklabels=kegg_pathways_gene_changes)
        ax.set_title(species_name)
        pylab.savefig('%s/%s_gene_change_kegg_pathway_distribution.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight',dpi=300)


    # repeat for gene changes:
    if len(gene_descriptions_gene_changes) >0:
        first=len(gene_descriptions_gene_changes)
        second=len(gene_descriptions_gene_changes_dict.keys())
        gene_occurrence_matrix=numpy.zeros((first,second))
        #fill in the matrix:
        sample_pair_idx=0
        for sample_pair in gene_descriptions_gene_changes_dict.keys():
            for pathway in gene_descriptions_gene_changes_dict[sample_pair]:
                pathway_idx=gene_descriptions_gene_changes.index(pathway)
                gene_occurrence_matrix[pathway_idx,sample_pair_idx]+=1
            sample_pair_idx+=1

        # plot this as a heatmap:
        import seaborn as sns; sns.set()

        pylab.figure(figsize=(6,10)) 
        sns.set(font_scale=0.4) # set the font
        ax = sns.heatmap(gene_occurrence_matrix,yticklabels=gene_descriptions_gene_changes) # plot heatmap
        ax.set_title(species_name)
        pylab.savefig('%s/%s_gene_change_description_distribution.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight',dpi=300)
'''
