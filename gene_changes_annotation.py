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



modification_divergence_threshold = 1e-03 #the threshold for deciding when something is a modification vs a replacement. Like most other things, it is an arbitrary choice. 


good_species_list = parse_midas_data.parse_good_species_list()

#good_species_list=['Eubacterium_rectale_56927']
for species_name in good_species_list: 
    

    ####################
    # Analyze the data #
    ####################
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



    gene_descriptions_gene_changes=[] # store the descriptions in this vector
    kegg_pathways_gene_changes=[] # store the pathways in this vector
    gene_descriptions_gene_changes_dict={}
    kegg_pathways_gene_changes_dict={}

    for sample_pair in temporal_change_map.keys():
        sample_1=sample_pair[0]
        sample_2=sample_pair[1]
        # check if this pair underwent a replacement event. If so, ignore, otherwise stats get messed up. 
        # find the index of sample_1 and sample_2 in snp_samples if the sample was included in the SNP output:
        include_sample_pair=False
        if sample_1 in snp_samples and sample_2 in snp_samples:
            sample_1_idx=snp_samples.index(sample_1)
            sample_2_idx=snp_samples.index(sample_2)
            if snp_substitution_rate[sample_1_idx, sample_2_idx] < modification_divergence_threshold:
                include_sample_pair=True
        else:
            include_sample_pair=True # not sure what to do if the snp matrix doesn't include this pair
        if include_sample_pair == True:
            gene_perr, gains, losses=calculate_temporal_changes.calculate_gains_losses_from_temporal_change_map(temporal_change_map, sample_1, sample_2, lower_threshold=0.05)
            all_changes=gains+losses
            gene_ids=[] # store the gene ids from all changes in this.
            #iterate through all_changes to pull out the gene_ids.
            for i in range(0, len(all_changes)):
                gene_ids.append(all_changes[i][0])
                
            # what kegg pathways are the genes in? Load all kegg_ids for all genomes in which the gene changes are in. 
            kegg_ids=parse_patric.load_kegg_annotations(gene_ids)    
            # what are the patric gene_descriptions for all genes that are changing? Load all gene_descritpions for the genomes in which the gene changes are in
            gene_descriptions=parse_patric.load_patric_gene_descriptions(gene_ids) 
            # now iterate again and pull out the pathway and gene_descriptions for hte genes in the gene_ids vector:
            kegg_pathways_gene_changes_dict[(sample_1,sample_2)]=[]
            gene_descriptions_gene_changes_dict[(sample_1,sample_2)]=[]
            for gene_id in gene_ids:
                #if kegg_ids[gene_id][0][0] != '':
                #print sample_1 +'\t'+sample_2 +'\t' + gene_id +'\t' + str(kegg_ids[gene_id])
                for i in range(0, len(kegg_ids[gene_id])):
                    if kegg_ids[gene_id][i][1] !='':
                        if kegg_ids[gene_id][i][1] not in kegg_pathways_gene_changes:
                            kegg_pathways_gene_changes.append(kegg_ids[gene_id][i][1])
                        kegg_pathways_gene_changes_dict[(sample_1,sample_2)].append(kegg_ids[gene_id][i][1])
                if 'hypothetical protein' not in gene_descriptions[gene_id]:
                    if gene_descriptions[gene_id] not in gene_descriptions_gene_changes:
                        gene_descriptions_gene_changes.append(gene_descriptions[gene_id]) 
                    gene_descriptions_gene_changes_dict[(sample_1,sample_2)].append(gene_descriptions[gene_id])


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
