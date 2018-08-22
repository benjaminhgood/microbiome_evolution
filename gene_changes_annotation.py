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
import core_gene_utils

import stats_utils
import matplotlib.colors as colors
import matplotlib.cm as cmx
from math import log10,ceil
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy.random import randint, choice
import pickle


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


if other_species_str=="":
    outFile=open('%sgene_changes_shuffle_all_species.txt' % parse_midas_data.analysis_directory, 'w' )
else:
    outFile=open('%sgene_changes_shuffle_%s.txt' % (parse_midas_data.analysis_directory, species_name) , 'w' )   

modification_difference_threshold = config.modification_difference_threshold
min_coverage = config.min_median_coverage
clade_divergence_threshold = 1e-02 # TODO: change to top level clade definition later
num_trials=100
min_sample_size = 5

within_host_classes = ['gains','losses','all','snps']

# Load subject and sample metadata
sys.stderr.write("Loading sample metadata...\n")
subject_sample_map = parse_HMP_data.parse_subject_sample_map()
sample_order_map = parse_HMP_data.parse_sample_order_map()
sys.stderr.write("Done!\n")


if other_species_str == "":
    good_species_list = parse_midas_data.parse_good_species_list()
else:
    good_species_list=[species_name]

# store all the species' data in a dictionary:
all_data={} 
#key=species
#value={}, key=gene, valuee=num times gene shows up


for species_name in good_species_list: 
    dummy_samples, sfs_map = parse_midas_data.parse_within_sample_sfs(species_name, allowed_variant_types=set(['1D','2D','3D','4D'])) 
    #
    # data structures for storing information for pickling later on
    all_species_gene_changes={}
    #all_species_gene_changes_category={}
    all_species_null={}
    all_data[species_name]={}
    #
    ####################
    # Analyze the data #
    ####################
    #
    # Only plot samples above a certain depth threshold that are "haploids"
    haploid_samples = diversity_utils.calculate_haploid_samples(species_name, debug=debug)
    #
    if len(haploid_samples) < min_sample_size:
        continue
    #
    same_sample_idxs, same_subject_idxs, diff_subject_idxs = parse_midas_data.calculate_ordered_subject_pairs(sample_order_map, haploid_samples)
    #
    snp_samples = set()
    sample_size = 0        
    for sample_pair_idx in xrange(0,len(same_subject_idxs[0])):
        #
        i = same_subject_idxs[0][sample_pair_idx]
        j = same_subject_idxs[1][sample_pair_idx]
        #
        snp_samples.add(haploid_samples[i])
        snp_samples.add(haploid_samples[j])
        #    
        sample_size += 1
        #    
    snp_samples = list(snp_samples)
    allowed_sample_set = set(snp_samples)
    #
    if sample_size < min_sample_size:
        continue
    #
    # load pre-computed data:
    sys.stderr.write("Loading pre-computed substitution rates for %s...\n" % species_name)
    substitution_rate_map = calculate_substitution_rates.load_substitution_rate_map(species_name)
    sys.stderr.write("Calculating matrix...\n")
    dummy_samples, snp_difference_matrix, snp_opportunity_matrix = calculate_substitution_rates.calculate_matrices_from_substitution_rate_map(substitution_rate_map, 'all', allowed_samples=snp_samples)
    snp_samples = dummy_samples
    sys.stderr.write("Done!\n")
    #
    sys.stderr.write("Loading pre-computed temporal changes for %s...\n" % species_name)
    temporal_change_map = calculate_temporal_changes.load_temporal_change_map(species_name)
    sys.stderr.write("Done!\n")
    #
    snp_substitution_rate = snp_difference_matrix*1.0/(snp_opportunity_matrix+(snp_opportunity_matrix==0))
    sys.stderr.write("Done!\n")   
    #
    # get all genome ids for this species' pan genome:
    genome_ids=parse_midas_data.get_ref_genome_ids(species_name)
    #
    # Load the non-shared genes (whitelisted genes):
    non_shared_genes = core_gene_utils.parse_non_shared_reference_genes(species_name)
    
    # load the gene descriptions for all genomes coresponding to this speceis:
    gene_descriptions=parse_patric.load_patric_gene_descriptions(genome_ids, non_shared_genes)
   #
    # create gene categories (poor proxy for GO terms):
    gene_categories, gene_category_map = parse_patric.cluster_patric_gene_descriptions(gene_descriptions)
    #
    # load the kegg ids for all genomes corresponding to this species:
    kegg_ids=parse_patric.load_kegg_annotations(genome_ids)    
    #
    # store null data in this to see how the actual data compares. 
    between_host_changes_gene_ids_null={} #dictionary which stores different trials (trial=key)
    present_gene_null={}
    pangenome_null={}
    for change_type in within_host_classes:
        between_host_changes_gene_ids_null[change_type]={}
        present_gene_null[change_type]={}
        pangenome_null[change_type]={}
        for trial in range(0,num_trials):
            between_host_changes_gene_ids_null[change_type][trial]=[]
            present_gene_null[change_type][trial]=[]        
            pangenome_null[change_type][trial]=[]
     #   
     #
    ##################
    # pangenome null #
    ##################
    #
    # load all pangenome genes for the species after clustering at 95% identity
    pangenome_gene_names, pangenome_new_species_names=parse_midas_data.load_pangenome_genes(species_name, non_shared_genes)
    #exclude any genes that are in the whitelisted set from pangenome_gene_names (the pangenome_new_species_names is not used):
    #
    #
    ###########################################
    # load data for between host changes null #
    ###########################################
    #
    # Load gene coverage information for species_name
    sys.stderr.write("Loading pangenome data for %s...\n" % species_name)
    gene_samples, gene_names, gene_presence_matrix, gene_depth_matrix, marker_coverages, gene_reads_matrix = parse_midas_data.parse_pangenome_data(species_name,allowed_samples=snp_samples)
    sys.stderr.write("Done!\n")
    #
    # compute gene cnv for constructing a null based on which genes are present later on.
    gene_copynum_matrix = gene_depth_matrix*1.0/(marker_coverages+(marker_coverages==0))
    #
    # convert gene_samples to list:
    gene_samples=gene_samples.tolist()
    #
    # convert gene names to numpy array:
    gene_names=numpy.array(gene_names)
    #
    # indexes for different subject pairs
    desired_samples = gene_samples
    #
    desired_same_sample_idxs, desired_same_subject_idxs, desired_diff_subject_idxs = parse_midas_data.calculate_ordered_subject_pairs(sample_order_map, desired_samples)
    #
    snp_sample_idx_map = parse_midas_data.calculate_sample_idx_map(desired_samples, snp_samples)
    gene_sample_idx_map = parse_midas_data.calculate_sample_idx_map(desired_samples, gene_samples)
    #
    same_subject_snp_idxs = parse_midas_data.apply_sample_index_map_to_indices(snp_sample_idx_map, desired_same_subject_idxs)  
    same_subject_gene_idxs = parse_midas_data.apply_sample_index_map_to_indices(gene_sample_idx_map, desired_same_subject_idxs)  
    #
    diff_subject_snp_idxs = parse_midas_data.apply_sample_index_map_to_indices(snp_sample_idx_map, desired_diff_subject_idxs)  
    diff_subject_gene_idxs = parse_midas_data.apply_sample_index_map_to_indices(gene_sample_idx_map, desired_diff_subject_idxs)  
    #
    between_host_gene_idxs = [] # store idxs of genes that change between hosts
    for sample_pair_idx in xrange(0,len(diff_subject_snp_idxs[0])):
        snp_i = diff_subject_snp_idxs[0][sample_pair_idx]
        snp_j = diff_subject_snp_idxs[1][sample_pair_idx]
        #
        i = diff_subject_gene_idxs[0][sample_pair_idx]
        j = diff_subject_gene_idxs[1][sample_pair_idx]
        if (marker_coverages[i]>min_coverage) and (marker_coverages[j]>min_coverage):
            if snp_substitution_rate[snp_i, snp_j] < clade_divergence_threshold:
                gene_idxs = gene_diversity_utils.calculate_gene_differences_between_idxs(i,j, gene_reads_matrix, gene_depth_matrix, marker_coverages)
                between_host_gene_idxs.extend(gene_idxs) # collect all gene changes occurring between hosts. Use this for the null.
    #
    #
    #
    #######################
    # within host changes #
    # present gene null -- construct a null consisting of any gene present at either time pt
    #######################
    #
    # store the actual data in this:
    within_host_changes_gene_ids={type:[] for type in within_host_classes}
    #
    # BG: Can't do it this way! Will pick up lots of diploids!
    #for sample_pair in temporal_change_map.keys():
    #    sample_1=sample_pair[0]
    #    sample_2=sample_pair[1]
    # 
    for sample_pair_idx in xrange(0,len(same_subject_snp_idxs[0])):
        #    
        i = same_subject_snp_idxs[0][sample_pair_idx]
        j = same_subject_snp_idxs[1][sample_pair_idx]
        #
        sample_i = snp_samples[i] 
        sample_j = snp_samples[j]
        #
        if not ((sample_i in allowed_sample_set) and (sample_j in allowed_sample_set)):
            continue
        #
        # Load SNP and gene changes!
        #
        # First SNP changes
        snp_opportunities, perr, mutations, reversions = calculate_temporal_changes.calculate_mutations_reversions_from_temporal_change_map(temporal_change_map, sample_i, sample_j)
        #
        # Look at higher threshold if error rate is too high
        if perr>=0.5:
            #
            # Calculate a more fine grained value!
            #
            dfs = numpy.array([0.6,0.7,0.8,0.9])
            perrs = diversity_utils.calculate_fixation_error_rate(sfs_map, sample_i, sample_j,dfs=dfs) * snp_opportunity_matrix[i, j]
            #
            if (perrs<0.5).any():
                # take most permissive one!
                perr_idx = numpy.nonzero(perrs<0.5)[0][0]
                df = dfs[perr_idx]
                perr = perrs[perr_idx]
                #
                # recalculate stuff!    
                perr, mutations, reversions = calculate_temporal_changes.calculate_mutations_reversions_from_temporal_change_map(temporal_change_map, sample_i, sample_j,lower_threshold=(1-df)/2.0, upper_threshold=(1+df)/2.0)
            #    
            else:
                df = 2
                perr = 1
                mutations = None
                reversions = None
                #
        if mutations==None or perr>=0.5:
            num_mutations = 0
            num_reversions = 0
            num_snp_changes = -1
        else:
            num_mutations = len(mutations)
            num_reversions = len(reversions)
            num_snp_changes = num_mutations+num_reversions
            #
        # Now do gene changes
        gene_opportunities, gene_perr, gains, losses = calculate_temporal_changes.calculate_gains_losses_from_temporal_change_map(temporal_change_map, sample_i, sample_j)
        all_changes=gains+losses
        #
        if (gains==None) or (gene_perr<-0.5) or (gene_perr>0.5):
            num_gains = 0
            num_losses = 0
            num_gene_changes = -1
        else:
            num_gains = len(gains)
            num_losses = len(losses)
            num_gene_changes = num_gains+num_losses
            #
        # Don't want to look at modifications or things with high error rates!
        if num_snp_changes<0 or num_snp_changes>=modification_difference_threshold:
            continue
        if (num_snp_changes<=0) and (num_gene_changes<=0):
            continue
        if num_snp_changes <0.5 and num_gene_changes <0.5:
            continue
        gene_change_dictionary={'gains':gains, 'losses':losses, 'all':all_changes}
        #                
        #iterate through all_changes to store the gene_ids.
        for change_type in ['gains','losses','all']:
            for i in range(0, len(gene_change_dictionary[change_type])):
                within_host_changes_gene_ids[change_type].append( gene_change_dictionary[change_type][i][0])
        #        
        # do same thing for SNP changes
        all_snp_changes = mutations+reversions
        snp_genes = set() # don't double count genes w/ 2 snps. probably same transfer event
        for snp_change in all_snp_changes:
            snp_genes.add(snp_change[0])
            #
        #print len(snp_genes)
        within_host_changes_gene_ids['snps'].extend(list(snp_genes))
        gene_change_dictionary['snps']=list(snp_genes) # added @12:08pm
        #
        #        
        # construct a null comprising of all genes present at either time point:
        sample_1_gene_idx = same_subject_gene_idxs[0][sample_pair_idx]
        sample_2_gene_idx = same_subject_gene_idxs[1][sample_pair_idx]
        #
        present_gene_idxs = []
        present_gene_idxs.extend( numpy.nonzero( (gene_copynum_matrix[:,sample_1_gene_idx]>0.5)*(gene_copynum_matrix[:,sample_1_gene_idx]<2))[0] )
        present_gene_idxs.extend( numpy.nonzero( (gene_copynum_matrix[:,sample_1_gene_idx]>0.5)*(gene_copynum_matrix[:,sample_1_gene_idx]<2))[0] )
        #
        #
        # sample 100x the number of within-host gene changes from the three different nulls:
        for change_type in within_host_classes:
            for trial in range(0,num_trials):
                between_gene_null_idxs = choice(between_host_gene_idxs, len(gene_change_dictionary[change_type]))
                present_gene_null_idxs = choice(present_gene_idxs, len(gene_change_dictionary[change_type]))
                # get the gene_ids for these idxs:
                between_host_changes_gene_ids_null[change_type][trial].extend( gene_names[between_gene_null_idxs].tolist() )
                present_gene_null[change_type][trial].extend( gene_names[present_gene_null_idxs].tolist())                
                pangenome_null[change_type][trial].extend( random.sample(pangenome_gene_names, len(all_changes)))
                


    #########################################################
    # what is the gene change annotation
    # kegg pathway annotation
    # gene_cateogry annotation
    #########################################################
    
    # annotate the real data:
    gene_descriptions_gene_changes={} # store the descriptions in this vector
    kegg_pathways_gene_changes={} # store the pathways in this vector
    gene_categories_gene_changes={} #stor categories in this vector
    
    for change_type in within_host_classes:
        gene_descriptions_gene_changes[change_type]=[]
        kegg_pathways_gene_changes[change_type]=[]
        gene_categories_gene_changes[change_type]=[]
        #
        #
        for gene_id in within_host_changes_gene_ids[change_type]:
            gene_descriptions_gene_changes[change_type].append(gene_descriptions[gene_id])
            gene_categories_gene_changes[change_type].append(gene_category_map[gene_id])
            for i in range(0, len(kegg_ids[gene_id])):
                kegg_pathways_gene_changes[change_type].append(kegg_ids[gene_id][i][1])
    
    # annotate the random draws:

    gene_descriptions_gene_changes_null_between_host={}
    kegg_pathways_gene_changes_null_between_host={}
    gene_categories_gene_changes_null_between_host={}

    gene_descriptions_gene_changes_null_present={}
    kegg_pathways_gene_changes_null_present={}
    gene_categories_gene_changes_null_present={}

    gene_descriptions_gene_changes_null_pangenome={}
    kegg_pathways_gene_changes_null_pangenome={}
    gene_categories_gene_changes_null_pangenome={}

    for change_type in within_host_classes:
        gene_descriptions_gene_changes_null_between_host[change_type]={}
        kegg_pathways_gene_changes_null_between_host[change_type]={}
        gene_categories_gene_changes_null_between_host[change_type]={}
        
        gene_descriptions_gene_changes_null_present[change_type]={}
        kegg_pathways_gene_changes_null_present[change_type]={}
        gene_categories_gene_changes_null_present[change_type]={}
        
        gene_descriptions_gene_changes_null_pangenome[change_type]={}
        kegg_pathways_gene_changes_null_pangenome[change_type]={}
        gene_categories_gene_changes_null_pangenome[change_type]={}
        
        for trial in range(0, num_trials):
            gene_descriptions_gene_changes_null_between_host[change_type][trial]=[]
            kegg_pathways_gene_changes_null_between_host[change_type][trial]=[]
            gene_categories_gene_changes_null_between_host[change_type][trial]=[]
            
            gene_descriptions_gene_changes_null_present[change_type][trial]=[]
            kegg_pathways_gene_changes_null_present[change_type][trial]=[]
            gene_categories_gene_changes_null_present[change_type][trial]=[]
            
            gene_descriptions_gene_changes_null_pangenome[change_type][trial]=[]
            kegg_pathways_gene_changes_null_pangenome[change_type][trial]=[]
            gene_categories_gene_changes_null_pangenome[change_type][trial]=[]
            
            for gene_id in between_host_changes_gene_ids_null[change_type][trial]:
                if gene_id in gene_descriptions:
                    gene_descriptions_gene_changes_null_between_host[change_type][trial].append(gene_descriptions[gene_id])
                    gene_categories_gene_changes_null_between_host[change_type][trial].append(gene_category_map[gene_id])
                for i in range(0, len(kegg_ids[gene_id])):
                    if gene_id in kegg_ids:
                        kegg_pathways_gene_changes_null_between_host[change_type][trial].append(kegg_ids[gene_id][i][1])
            
            for gene_id in present_gene_null[change_type][trial]:
                if gene_id in gene_descriptions:
                    gene_descriptions_gene_changes_null_present[change_type][trial].append(gene_descriptions[gene_id])
                    gene_categories_gene_changes_null_present[change_type][trial].append(gene_category_map[gene_id])
                for i in range(0, len(kegg_ids[gene_id])):
                    if gene_id in kegg_ids:
                        kegg_pathways_gene_changes_null_present[change_type][trial].append(kegg_ids[gene_id][i][1])
            
            for gene_id in pangenome_null[change_type][trial]:
                if gene_id in gene_descriptions:
                    gene_descriptions_gene_changes_null_pangenome[change_type][trial].append(gene_descriptions[gene_id])
                    gene_categories_gene_changes_null_pangenome[change_type][trial].append(gene_category_map[gene_id]) 
                for i in range(0, len(kegg_ids[gene_id])):
                    if gene_id in kegg_ids:
                        kegg_pathways_gene_changes_null_pangenome[change_type][trial].append(kegg_ids[gene_id][i][1])
                



    #############################
    # compute empirical p-value #
    #############################

    # how likely is it to see a given gene change as many times as it does given all the nulls?

    # count the number of times a gene shows up in within-host changes
    all_gene_changes={}
    for change_type in within_host_classes:
        if change_type=='snps':
            print "Getting ready for snps!"
            
        for gene in gene_descriptions_gene_changes[change_type]:
            print gene, "in snps"
            if gene not in all_gene_changes:
                #all_gene_changes[gene]={0 for c in within_host_classes}
                all_gene_changes[gene]={}
                for c in within_host_classes:
                    all_gene_changes[gene][c]=0
            all_gene_changes[gene][change_type] +=1

    # count the number of times a gene shows up in between-host changes in the num_trials
    all_gene_changes_null_between_host={} # key == gene
    all_gene_changes_null_present={}
    all_gene_changes_null_pangenome={}
    

    for trial in range(0,num_trials):
        for gene in gene_descriptions_gene_changes_null_between_host['all'][trial]:
            all_gene_changes_null_between_host[gene] = {c : [] for c in within_host_classes}
        for gene in gene_descriptions_gene_changes_null_present['all'][trial]: 
            all_gene_changes_null_present[gene] = {c : [] for c in within_host_classes}
        for gene in gene_descriptions_gene_changes_null_pangenome['all'][trial]: 
            all_gene_changes_null_pangenome[gene] = {c : [] for c in within_host_classes}

    for change_type in within_host_classes:
        for trial in range(0,num_trials): # iterate through the 100 trials
            for gene in all_gene_changes_null_between_host.keys():
                num_occurrences=gene_descriptions_gene_changes_null_between_host[change_type][trial].count(gene)
                all_gene_changes_null_between_host[gene][change_type].append(num_occurrences)
                #
            for gene in all_gene_changes_null_present.keys():        
                num_occurrences=gene_descriptions_gene_changes_null_present[change_type][trial].count(gene)
                all_gene_changes_null_present[gene][change_type].append(num_occurrences)
                #
            for gene in all_gene_changes_null_pangenome.keys():
                num_occurrences=gene_descriptions_gene_changes_null_pangenome[change_type][trial].count(gene)
                all_gene_changes_null_pangenome[gene][change_type].append(num_occurrences)



    # iterate through all gene changes and quantify observed vs expected:
    
    #outFile.write(species_name + '\ngene_name\tprob_between_hosts\tprob_present\tprob_pangenome\n')
   
    #for gene in all_gene_changes.keys():
        #if gene in all_gene_changes_null_between_host:
        #    tmp_null = numpy.array(all_gene_changes_null_between_host[gene])
        #    prob_between_hosts=sum(all_gene_changes[gene] <= tmp_null)/float(num_trials)
        #else:
        #    prob_between_hosts=0
        #if gene in all_gene_changes_null_present:
        #    tmp_null = numpy.array(all_gene_changes_null_present[gene])
        #    prob_present=sum(all_gene_changes[gene] <= tmp_null)/float(num_trials)
        #else:
        #    prob_present=0
        #if gene in all_gene_changes_null_pangenome:
        #    tmp_null = numpy.array(all_gene_changes_null_pangenome[gene])
        #    prob_pangenome=sum(all_gene_changes[gene] <= tmp_null)/float(num_trials)
        #else:
        #    prob_pangenome=0
        #store data in a dictionary for cross-species plotting after running on all species
        #all_species_gene_changes[gene]=[all_gene_changes[gene],prob_between_hosts,prob_present,prob_pangenome]
        #all_species_gene_changes[gene]=all_gene_changes[gene]
        #outFile.write(gene + '\t' + str(all_gene_changes[gene]) +'\t' + str(prob_between_hosts) + '\t'+ str(prob_present)+ '\t' + str(prob_pangenome) + '\n')
        
    #outFile.write('\n')
    
    ###################################################################
    # check whether there is an enrichment of genes in a gene_category
    ##########################################################

    '''
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
    '''
    # store the nulls:
    #all_species_null={'between_host_genes':all_gene_changes_null_between_host,'present_genes': all_gene_changes_null_present, 'pangenome_genes': all_gene_changes_null_pangenome, 'between_host_category': all_gene_changes_category_null_between_host, 'present_category':all_gene_changes_category_null_present, 'pangenome_category':all_gene_changes_category_null_pangenome}
    all_species_null={'between_host_genes':all_gene_changes_null_between_host,'present_genes': all_gene_changes_null_present, 'pangenome_genes': all_gene_changes_null_pangenome}

    #store all the data 
    #all_data[species_name]={'gene_changes':all_gene_changes, 'gene_changes_category':all_species_gene_changes_category, 'null':all_species_null}
    all_data[species_name]={'gene_changes':all_gene_changes, 'null':all_species_null}

if other_species_str=="":
    pickle.dump( all_data, open( "/pollard/home/ngarud/tmp_intermediate_files/all_species_gene_changes.p", "wb" ) )
else:
    pickle.dump( all_data, open( "/pollard/home/ngarud/tmp_intermediate_files/%s_gene_changes.p" % species_name, "wb" ) )  





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
