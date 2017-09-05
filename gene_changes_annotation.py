# Within-host gene change annotation

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
import species_phylogeny_utils

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
    outFile=open('%sgene_changes_shuffle_all_species.txt' % parse_midas_data.analysis_directory, 'w' )
else:
    outFile=open('%sgene_changes_shuffle_%s.txt' %(parse_midas_data.analysis_directory, species_name) , 'w' )   

modification_divergence_threshold = 1e-03 #the threshold for deciding when something is a modification vs a replacement. Like most other things, it is an arbitrary choice. 
min_coverage=20
clade_divergence_threshold = 1e-02
num_trials=100

if other_species_str == "":
    good_species_list = parse_midas_data.parse_good_species_list()
else:
    good_species_list=[species_name]

# store all the species' data in a dictionary:
all_data={} 
#key=species
#value={}, key=gene, valuee=num times gene shows up

for species_name in good_species_list: 
    
    # data structures for storing information for pickling later on
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
    

    # store null data in this to see how the actual data compares. 
    between_host_changes_gene_ids_null={} #dictionary which stores different trials (trial=key)
    present_gene_null={}
    pangenome_null={} # same except null is pangenome
    
    for trial in range(0,num_trials):
        between_host_changes_gene_ids_null[trial]=[]
        present_gene_null[trial]=[]        
        pangenome_null[trial]=[]
        

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
    
    # compute gene cnv for constructing a null based on which genes are present later on.
    gene_copynum_matrix = gene_depth_matrix*1.0/(marker_coverages+(marker_coverages==0))
    
    # convert gene_samples to list:
    gene_samples=gene_samples.tolist()

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
    
    

    #######################
    # within host changes #
    # present gene null -- construct a null consisting of any gene present at either time pt
    #######################

    # store the actual data in this:
    within_host_changes_gene_ids=[]
    

    within_host_change_distribution=[] # use this for constructing the null
    for sample_pair in temporal_change_map.keys():
        sample_1=sample_pair[0]
        sample_2=sample_pair[1]
        if sample_1 in gene_samples and sample_2 in gene_samples:
            sample_1_gene_idx=gene_samples.index(sample_1)
            sample_2_gene_idx=gene_samples.index(sample_2)
            #check if the marker coverages are in range:
            if (marker_coverages[sample_2_gene_idx]>min_coverage) and (marker_coverages[sample_2_gene_idx]>min_coverage): 
                # check if this pair underwent a replacement event. If so, ignore, otherwise stats get messed up. 
                # find the index of sample_1 and sample_2 in snp_samples to see if this is a CPS sample pair:
                include_sample_pair=False
                if sample_1 in snp_samples and sample_2 in snp_samples:
                    sample_1_idx=snp_samples.index(sample_1)
                    sample_2_idx=snp_samples.index(sample_2)
                    if snp_substitution_rate[sample_1_idx, sample_2_idx] < modification_divergence_threshold:
                        gene_perr, gains, losses=calculate_temporal_changes.calculate_gains_losses_from_temporal_change_map(temporal_change_map, sample_1, sample_2, lower_threshold=0.05)
                        all_changes=gains+losses
                    #within_host_change_distribution.append(len(all_changes))
                    #iterate through all_changes to store the gene_ids.
                    for i in range(0, len(all_changes)):
                        within_host_changes_gene_ids.append(all_changes[i][0])
                    #
                    # construct a null comprising of all genes present at either time point:
                    present_gene_idxs = []
                    present_gene_idxs.extend( numpy.nonzero( (gene_copynum_matrix[:,sample_1_gene_idx]>0.5)*(gene_copynum_matrix[:,sample_1_gene_idx]<2))[0] )
                    present_gene_idxs.extend( numpy.nonzero( (gene_copynum_matrix[:,sample_1_gene_idx]>0.5)*(gene_copynum_matrix[:,sample_1_gene_idx]<2))[0] )
                    #
                    # sample 100x the number of within-host gene changes from the three different nulls:
                    for trial in range(0,num_trials):
                        between_gene_null_idxs = choice(between_host_gene_idxs, len(all_changes))
                        present_gene_null_idxs = choice(present_gene_idxs, len(all_changes))
                        # get the gene_ids for these idxs:
                        between_host_changes_gene_ids_null[trial].extend(gene_names[between_gene_null_idxs].tolist())
                        present_gene_null[trial].extend(gene_names[present_gene_null_idxs].tolist())                
                        pangenome_null[trial].extend(random.sample(pangenome_gene_names, len(all_changes)))
                


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

    gene_descriptions_gene_changes_null_present={}
    kegg_pathways_gene_changes_null_present={}
    gene_categories_gene_changes_null_present={}

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

        gene_descriptions_gene_changes_null_present[trial]=[]
        kegg_pathways_gene_changes_null_present[trial]=[]
        gene_categories_gene_changes_null_present[trial]=[]
        for gene_id in present_gene_null[trial]:
            if gene_id in gene_descriptions:
                gene_descriptions_gene_changes_null_present[trial].append(gene_descriptions[gene_id])
                gene_categories_gene_changes_null_present[trial].append(gene_category_map[gene_id])
            for i in range(0, len(kegg_ids[gene_id])):
                if gene_id in kegg_ids:
                    kegg_pathways_gene_changes_null_present[trial].append(kegg_ids[gene_id][i][1])

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

    # how likely is it to see a given gene change as many times as it does given all the nulls?

    # count the number of times a gene shows up in within-host changes
    all_gene_changes={}
    for gene in gene_descriptions_gene_changes:
        if gene not in all_gene_changes:
            all_gene_changes[gene]=1
        else:
            all_gene_changes[gene] +=1

    # count the number of times a gene shows up in between-host changes in the num_trails
    all_gene_changes_null_between_host={} # key == gene
    all_gene_changes_null_present={}
    all_gene_changes_null_pangenome={}
    
    for trial in range(0,num_trials):
        for gene in gene_descriptions_gene_changes_null_between_host[trial]:
            all_gene_changes_null_between_host[gene]=[]
        for gene in gene_descriptions_gene_changes_null_present[trial]: 
            all_gene_changes_null_present[gene]=[]
        for gene in gene_descriptions_gene_changes_null_pangenome[trial]: 
            all_gene_changes_null_pangenome[gene]=[]

    for gene in all_gene_changes_null_between_host.keys():
        for trial in range(0,num_trials): # iterate through the 100 trials
            num_occurrences=gene_descriptions_gene_changes_null_between_host[trial].count(gene)
            all_gene_changes_null_between_host[gene].append(num_occurrences)

    for gene in all_gene_changes_null_present.keys():
        for trial in range(0,num_trials): # iterate through the 100 trials
            num_occurrences=gene_descriptions_gene_changes_null_present[trial].count(gene)
            all_gene_changes_null_present[gene].append(num_occurrences)

    for gene in all_gene_changes_null_pangenome.keys():
        for trial in range(0,num_trials): # iterate through the 100 trials
            num_occurrences=gene_descriptions_gene_changes_null_pangenome[trial].count(gene)
            all_gene_changes_null_pangenome[gene].append(num_occurrences)



    # iterate through all gene changes and quantify observed vs expected:
    
    outFile.write(species_name + '\ngene_name\tprob_between_hosts\tprob_present\tprob_pangenome\n')
   
    for gene in all_gene_changes.keys():
        if gene in all_gene_changes_null_between_host:
            tmp_null = numpy.array(all_gene_changes_null_between_host[gene])
            prob_between_hosts=sum(all_gene_changes[gene] <= tmp_null)/float(num_trials)
        else:
            prob_between_hosts=0
        if gene in all_gene_changes_null_present:
            tmp_null = numpy.array(all_gene_changes_null_present[gene])
            prob_present=sum(all_gene_changes[gene] <= tmp_null)/float(num_trials)
        else:
            prob_present=0
        if gene in all_gene_changes_null_pangenome:
            tmp_null = numpy.array(all_gene_changes_null_pangenome[gene])
            prob_pangenome=sum(all_gene_changes[gene] <= tmp_null)/float(num_trials)
        else:
            prob_pangenome=0
        #store data in a dictionary for cross-species plotting after running on all species
        all_species_gene_changes[gene]=[all_gene_changes[gene],prob_between_hosts,prob_present,prob_pangenome]
        
        outFile.write(gene + '\t' + str(all_gene_changes[gene]) +'\t' + str(prob_between_hosts) + '\t'+ str(prob_present)+ '\t' + str(prob_pangenome) + '\n')
        
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
    all_gene_changes_category_null_present={}
    all_gene_changes_category_null_pangenome={} # key == gene
    
    
    for trial in range(0,num_trials):
        for gene in gene_categories_gene_changes_null_between_host[trial]:
            all_gene_changes_category_null_between_host[gene]=[]
        for gene in gene_categories_gene_changes_null_present[trial]:
            all_gene_changes_category_null_present[gene]=[] 
        for gene in gene_categories_gene_changes_null_pangenome[trial]:
            all_gene_changes_category_null_pangenome[gene]=[]

    for gene in all_gene_changes_category_null_between_host.keys():
        for trial in range(0,num_trials): # iterate through the 100 trials
            num_occurrences=gene_categories_gene_changes_null_between_host[trial].count(gene)
            all_gene_changes_category_null_between_host[gene].append(num_occurrences)
            
    for gene in all_gene_changes_category_null_present.keys():
        for trial in range(0,num_trials): # iterate through the 100 trials
            num_occurrences=gene_categories_gene_changes_null_present[trial].count(gene)
            all_gene_changes_category_null_present[gene].append(num_occurrences)

    for gene in all_gene_changes_category_null_pangenome.keys():
        for trial in range(0,num_trials): # iterate through the 100 trials
            num_occurrences=gene_categories_gene_changes_null_pangenome[trial].count(gene)
            all_gene_changes_category_null_pangenome[gene].append(num_occurrences)


    # iterate through all gene changes and quantify observed vs expected:
    outFile.write(species_name + '\ngene_name\tprob_between_hosts\tprob_present\tprob_pangenome\n')

    for gene in all_gene_changes_category.keys():
        if gene in all_gene_changes_category_null_between_host:
            tmp_null = numpy.array(all_gene_changes_category_null_between_host[gene])
            prob_between_hosts=sum(all_gene_changes_category[gene] <= tmp_null)/100.0
        else:
            prob_between_hosts=0
        if gene in all_gene_changes_category_null_present:
            tmp_null = numpy.array(all_gene_changes_category_null_present[gene])
            prob_present=sum(all_gene_changes_category[gene] <= tmp_null)/100.0
        else:
            prob_present=0
        if gene in all_gene_changes_category_null_pangenome:
            tmp_null = numpy.array(all_gene_changes_category_null_pangenome[gene])
            prob_pangenome=sum(all_gene_changes_category[gene] <= tmp_null)/100.0
        else:
            prob_pangenome=0
        all_species_gene_changes_category[gene]=[all_gene_changes_category[gene],prob_between_hosts,prob_present,prob_pangenome]
        outFile.write(gene + '\t' + str(all_gene_changes_category[gene]) +'\t' + str(prob_between_hosts) + '\t' + str(prob_present)+ '\t' +  str(prob_pangenome) + '\n')

    outFile.write('\n')

    # store the nulls:
    all_species_null={'between_host_genes':all_gene_changes_null_between_host,'present_genes': all_gene_changes_null_present, 'pangenome_genes': all_gene_changes_null_pangenome, 'between_host_category': all_gene_changes_category_null_between_host, 'present_category':all_gene_changes_category_null_present, 'pangenome_category':all_gene_changes_category_null_pangenome}

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


#  sum  all gene changes  across species
all_data['all_species']={}
for gene_type in ['gene_changes','gene_changes_category']:
    all_data['all_species']['gene_changes']={}
    all_data['all_species']['gene_changes_category']={}

for species_name in all_data.keys():
    if species_name != 'all_species':
        for gene_type in ['gene_changes','gene_changes_category']:
            for gene in all_data[species_name][gene_type].keys():
                if gene not in all_data['all_species'][gene_type].keys():
                    all_data['all_species'][gene_type][gene] = [all_data[species_name][gene_type][gene][0],0,0,0]
                else:
                    all_data['all_species'][gene_type][gene][0] +=all_data[species_name][gene_type][gene][0]
                                                                 

# compute a null for all species by combining all trials across speceis:
all_data['all_species']['null']={}
for null_type in ['between_host_genes', 'present_genes', 'pangenome_genes', 'between_host_category','present_category','pangenome_category']:
    all_data['all_species']['null'][null_type]={}

for species_name in all_data.keys(): 
    print species_name
    if species_name != 'all_species' and species_name !='Escherichia_coli_58110':
        for null_type in ['between_host_genes', 'present_genes', 'pangenome_genes', 'between_host_category','present_category','pangenome_category']:
            for gene in all_data[species_name]['null'][null_type].keys():
                if gene not in  all_data['all_species']['null'][null_type].keys():
                    all_data['all_species']['null'][null_type][gene]=all_data[species_name]['null'][null_type][gene]
                else:
                    for i in range(0, num_trials):
                        all_data['all_species']['null'][null_type][gene][i] += all_data[species_name]['null'][null_type][gene][i]

                

# compute empirical p-value for all_species 
for gene in all_data['all_species']['gene_changes'].keys():
    for null_type in ['between_host_genes', 'present_genes', 'pangenome_genes']:
        if gene in all_data['all_species']['null'][null_type].keys():
            tmp_null = numpy.array(all_data['all_species']['null'][null_type][gene])
            prob=sum(all_data['all_species']['gene_changes'][gene][0] <= tmp_null)/float(num_trials)
        else:
            prob=0
        if null_type=='between_host_genes':
            all_data['all_species']['gene_changes'][gene][1]=prob
        elif null_type=='present_genes':
            all_data['all_species']['gene_changes'][gene][2]=prob
        else:
            all_data['all_species']['gene_changes'][gene][3]=prob
        


# repeat for gene_category

for gene in all_data['all_species']['gene_changes_category'].keys():
    for null_type in ['between_host_category','present_category','pangenome_category']:
        if gene in all_data['all_species']['null'][null_type].keys():
            tmp_null = numpy.array(all_data['all_species']['null'][null_type][gene])
            prob=sum(all_data['all_species']['gene_changes_category'][gene][0] <= tmp_null)/float(num_trials)
        else:
            prob=0
        if null_type=='between_host_category':
            all_data['all_species']['gene_changes_category'][gene][1]=prob
        elif null_type=='present_category':
            all_data['all_species']['gene_changes_category'][gene][2]=prob
        else:
            all_data['all_species']['gene_changes_category'][gene][3]=prob
        

# print the p-values of the genes in the all species gene changes:
outFile_gene_change=open('%sgene_changes_accross_species.txt' %  parse_midas_data.analysis_directory,'w')

# store data for plotting cdf:
p_val_between_array=[]
p_val_present_array=[]
p_val_pangenome_array=[]

gene_order_sorted = sorted(all_data['all_species']['gene_changes'].items(), key=operator.itemgetter(1))
for i in range(0, len(gene_order_sorted)):
    gene=gene_order_sorted[i][0]
    counts=all_data['all_species']['gene_changes'][gene][0]
    p_val_between=all_data['all_species']['gene_changes'][gene][1]
    p_val_present=all_data['all_species']['gene_changes'][gene][2]
    p_val_pangenome=all_data['all_species']['gene_changes'][gene][3]
    p_val_between_array.append(p_val_between)
    p_val_present_array.append(p_val_present)
    p_val_pangenome_array.append(p_val_pangenome)
    print gene + '\t' + str(counts) + '\t' +str(p_val_between) +'\t' + str(p_val_present) +'\t' + str(p_val_pangenome) 
    outFile_gene_change.write(gene + '\t' +str(counts) +'\t'+ str(p_val_between) +'\t' + str(p_val_present) +'\t' + str(p_val_pangenome) + '\n' )


outFile_gene_change.close()

p_val_between_array=numpy.asarray(sorted(p_val_between_array))
p_val_present_array=numpy.asarray(sorted(p_val_present_array))
p_val_pangenome_array=numpy.asarray(sorted(p_val_pangenome_array))

# print the p-values of the genes in the all species category:
outFile_category=open('%sgene_changes_accross_species_category.txt' %  parse_midas_data.analysis_directory,'w')

gene_order_sorted = sorted(all_data['all_species']['gene_changes_category'].items(), key=operator.itemgetter(1))
for i in range(0, len(gene_order_sorted)):
    gene=gene_order_sorted[i][0]
    counts=all_data['all_species']['gene_changes_category'][gene][0]
    p_val_between=all_data['all_species']['gene_changes_category'][gene][1]
    p_val_present=all_data['all_species']['gene_changes_category'][gene][2]
    p_val_pangenome=all_data['all_species']['gene_changes_category'][gene][3]
    print gene + '\t' + str(counts) + '\t' +str(p_val_between) +'\t' + str(p_val_present) +'\t' + str(p_val_pangenome) 
    outFile_category.write(gene + '\t' +str(counts) +'\t'+ str(p_val_between) +'\t' + str(p_val_present) +'\t' + str(p_val_pangenome) + '\n' )

outFile_category.close()


##################################################
# list of common genes that show up              #
##################################################

keywords={}
keywords['hypothetical'] = 'ypothetical'
keywords['conjugative']='onjugati'
keywords['transposase']='anspos'
keywords['tRNA']='tRNA'
keywords['recombinase']='ecombinase'
keywords['integrase']='ntegrase'
keywords['ABC transporter']='ABC'
keywords['mob']='mob'
keywords['ATP']='ATP'
keywords['excisionase']='xcisionase'
keywords['mobile']='obile'
keywords['transmembrane']='embrane'
keywords['replication']='eplication'
keywords['phase']='hage'
keywords['mobilization']='obilization'
keywords['regulator']='egulator'
keywords['transcription']='anscription'
keywords['endonuclease']='ndonuclease'
keywords['DNA']='DNA'
keywords['RNA']='RNA'
keywords['toxin']='toxin'
keywords['ferredoxin']='erredoxin'
keywords['restriction']='estriction'
keywords['replication']='eplication'
keywords['transferase']='ansferase'
keywords['reductase']='eductase'
keywords['phosphatase']='phosphatase'
keywords['plasmid']='plasmid'
keywords['helicase']='elicase'
keywords['kinase']='kinase'
keywords['dehydrogenase']='dehydrogenase'
keywords['drug']='drug'
keywords['cell wall']='ell wall'
keywords['other']='other'

common_genes={}

# key=keyword
# value={}
# num={}
# genes={}
for keyword in keywords:
    common_genes[keyword]={'num':0, 'genes':[]}

others=[]
for gene in all_data['all_species']['gene_changes']:
    keyword_found=False
    for keyword in keywords.keys():
        if keywords[keyword] in gene:
            common_genes[keyword]['num']+=all_data['all_species']['gene_changes'][gene][0]
            common_genes[keyword]['genes'].append(gene)
            keyword_found=True
    if keyword_found==False:
        common_genes['other']['num']+=all_data['all_species']['gene_changes'][gene][0]
        common_genes['other']['genes'].append(gene)

outFile_keywords=open('%sgene_changes_accross_species_keywords.txt' %  parse_midas_data.analysis_directory,'w')

for gene in common_genes.keys():
    outFile_keywords.write(gene + '\t' + str(common_genes[gene]['num']) +'\t' + ','.join(common_genes[gene]['genes']) +'\n')

outFile_keywords.close()
    

###################################################
# plot CDF of p-values for the different nulls 
###################################################

pylab.figure(figsize=(6,6))
prevalence_axis = pylab.subplot(111)

prevalence_axis.set_ylabel('Fraction genes $\leq p$',labelpad=2)
prevalence_axis.set_xlabel('P-val, $p$',labelpad=2)
prevalence_axis.set_xlim([0,1.05])
prevalence_axis.set_ylim([0,1.1])

prevalence_axis.spines['top'].set_visible(False)
prevalence_axis.spines['right'].set_visible(False)
prevalence_axis.get_xaxis().tick_bottom()
prevalence_axis.get_yaxis().tick_left()

xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(p_val_present_array)
prevalence_axis.step(xs,1-ns*1.0/ns[0],'b-',label='Within-host present genes',zorder=2)

xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(p_val_between_array)
prevalence_axis.step(xs,1-ns*1.0/ns[0],'r-',label='Between-host gene changes',zorder=1)

xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(p_val_pangenome_array)
prevalence_axis.step(xs,1-ns*1.0/ns[0],'k-',label='Pangenome',zorder=0)

#prevalence_axis.set_ylim([0,1.1])
#prevalence_axis.set_xlim([0,1.05])

prevalence_axis.legend(loc='upper right',frameon=False,fontsize=4)

pylab.savefig('%s/p_vals_gene_change_annotation.png' % (parse_midas_data.analysis_directory),bbox_inches='tight',dpi=300)







##############################################################
# 2D matrix of the p-values (genes x species) for a heatmap:
##############################################################
num_species=len(all_data.keys())
num_genes=len(all_data['all_species']['gene_changes'].keys())-1 # remove hypothetical protein
gene_occurrence_matrix=numpy.zeros((num_genes, num_species))

# sort genes based on how abundant they are accross species:
gene_order_sorted = sorted(all_data['all_species']['gene_changes'].items(), key=operator.itemgetter(1))
gene_order=[]
for i in range(0, len(gene_order_sorted)):
    if gene_order_sorted[i][0] !='hypothetical protein':
        gene_order.append(gene_order_sorted[i][0])

species_order=species_phylogeny_utils.sort_phylogenetically(all_data.keys())
species_order.append('all_species')



# fill matrix with probability under the present gene null 
species_idx=0
for species_name in species_order:
    gene_idx=0
    for gene in gene_order:
        if gene != 'hypothetical protein':
            if gene in all_data[species_name]['gene_changes']:
                if all_data[species_name]['gene_changes'][gene][2] <=0.05:
                    gene_occurrence_matrix[gene_idx, species_idx]=2
                else:
                    gene_occurrence_matrix[gene_idx, species_idx]=1
            gene_idx +=1
    species_idx +=1



import seaborn as sns; sns.set()
pylab.figure(figsize=(10,15))   
sns.set(font_scale=0.4)
ax = sns.heatmap(gene_occurrence_matrix, yticklabels=gene_order, xticklabels=species_order) 
pylab.savefig('%s/gene_change_cross_species_distribution_present_null.png' % (parse_midas_data.analysis_directory),bbox_inches='tight',dpi=300)



gene_occurrence_matrix=numpy.zeros((num_genes, num_species))

# fill matrix with probability under the between-host gene changes null 
species_idx=0
for species_name in species_order:
    gene_idx=0
    for gene in gene_order:
        if gene != 'hypothetical protein':
            if gene in all_data[species_name]['gene_changes']:
                if all_data[species_name]['gene_changes'][gene][1] <=0.05:
                    gene_occurrence_matrix[gene_idx, species_idx]=2
                else:
                    gene_occurrence_matrix[gene_idx, species_idx]=1
            gene_idx +=1
    species_idx +=1



import seaborn as sns; sns.set()
pylab.figure(figsize=(10,15))   
sns.set(font_scale=0.4)
ax = sns.heatmap(gene_occurrence_matrix, yticklabels=gene_order, xticklabels=species_order) 
pylab.savefig('%s/gene_change_cross_species_distribution_between_host_null.png' % (parse_midas_data.analysis_directory),bbox_inches='tight',dpi=300)


#####
# repeat for gene_changes_category

# put everything in a 2D matrix so that I can plot a heatmap:
num_species=len(all_data.keys())
num_genes=len(all_data['all_species']['gene_changes_category'].keys())-1 # remove hypothetical protein
gene_occurrence_matrix=numpy.zeros((num_genes, num_species))

# sort genes based on how abundant they are accross species:
gene_order_sorted = sorted(all_data['all_species']['gene_changes_category'].items(), key=operator.itemgetter(1))
gene_order=[]
for i in range(0, len(gene_order_sorted)):
    if gene_order_sorted[i][0] !='hypothetical protein':
        gene_order.append(gene_order_sorted[i][0])

species_order=species_phylogeny_utils.sort_phylogenetically(all_data.keys())
species_order.append('all_species')



# fill based on the probability under the present gene model
species_idx=0
for species_name in species_order:
    gene_idx=0
    for gene in gene_order:
        if gene != 'hypothetical protein':
            if gene in all_data[species_name]['gene_changes_category']:
                if all_data[species_name]['gene_changes_category'][gene][2] <=0.05:
                    gene_occurrence_matrix[gene_idx, species_idx] = 2
                else:
                    gene_occurrence_matrix[gene_idx, species_idx] = 1
            gene_idx +=1
    species_idx +=1



import seaborn as sns; sns.set()
pylab.figure(figsize=(10,15))   
sns.set(font_scale=0.4)
ax = sns.heatmap(gene_occurrence_matrix, yticklabels=gene_order, xticklabels=species_order) 
pylab.savefig('%s/gene_change_cross_species_distribution_category_present_null.png' % (parse_midas_data.analysis_directory),bbox_inches='tight',dpi=300)


# fill based on the probability under the between-host model
species_idx=0
for species_name in species_order:
    gene_idx=0
    for gene in gene_order:
        if gene != 'hypothetical protein':
            if gene in all_data[species_name]['gene_changes_category']:
                if all_data[species_name]['gene_changes_category'][gene][1] <=0.05:
                    gene_occurrence_matrix[gene_idx, species_idx] = 2
                else:
                    gene_occurrence_matrix[gene_idx, species_idx] = 1
            gene_idx +=1
    species_idx +=1



import seaborn as sns; sns.set()
pylab.figure(figsize=(10,15))   
sns.set(font_scale=0.4)
ax = sns.heatmap(gene_occurrence_matrix, yticklabels=gene_order, xticklabels=species_order) 
pylab.savefig('%s/gene_change_cross_species_distribution_category_between_host_null.png' % (parse_midas_data.analysis_directory),bbox_inches='tight',dpi=300)





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
