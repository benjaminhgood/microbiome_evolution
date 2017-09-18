import config
import parse_midas_data
import parse_HMP_data
import os.path
import pylab
import sys
import numpy

import diversity_utils
import gene_diversity_utils

import stats_utils
from math import log10,ceil
from numpy.random import randint

intermediate_filename = '%ssubstitution_rates.txt' % (parse_midas_data.data_directory)
    

min_coverage = config.min_median_coverage
alpha = 0.5 # Confidence interval range for rate estimates
low_pi_threshold = 1e-03
low_divergence_threshold = 1e-03
min_change = 0.8
min_sample_size = 10
allowed_variant_types = set(['1D','2D','3D','4D'])


def load_substitution_rate_map(species_name):
# This definition is called whenever another script downstream uses the output of this data.

    substitution_rate_map = {}

    file = open(intermediate_filename,"r")
    file.readline() # header
    for line in file:
        items = line.split(",")
        if items[0].strip()!=species_name:
            continue
            
        sample_1 = items[1].strip()
        sample_2 = items[2].strip()
        type = items[3].strip()
        num_changes = float(items[4])
        num_opportunities = float(items[5])
        
        sample_pair = (sample_1, sample_2)
        if sample_pair not in substitution_rate_map:
            substitution_rate_map[sample_pair] = {}
            
        substitution_rate_map[(sample_1,sample_2)][type] = (num_changes, num_opportunities)
    
    return substitution_rate_map
    
def calculate_matrices_from_substitution_rate_map(substitution_rate_map, type, allowed_samples=[]):
# once the map is loaded, then we can compute rate matrices in this definition (so, it relies on the previous def)    

    sample_set = set([])
    for sample_1, sample_2 in substitution_rate_map.keys():
        sample_set.add(sample_1)
        sample_set.add(sample_2)
    
    if len(allowed_samples)>0:
        allowed_sample_set = set(allowed_samples)    
    else:
        allowed_sample_set = sample_set
        
    sample_set = sample_set & allowed_sample_set
    samples = list(sample_set)
    
    sample_idx_map = {samples[i]: i for i in xrange(0,len(samples))}
    
    difference_matrix = numpy.zeros((len(samples), len(samples)))*1.0
    opportunity_matrix = numpy.zeros_like(difference_matrix)
    
    for sample_pair in substitution_rate_map.keys():
        
        sample_i = sample_pair[0]
        sample_j = sample_pair[1]
        
        if not (sample_i in sample_set and sample_j in sample_set):
            continue
        
        i = sample_idx_map[sample_i]
        j = sample_idx_map[sample_j]
        
        num_changes, num_opportunities = substitution_rate_map[sample_pair][type]
        
        difference_matrix[i,j] = num_changes
        difference_matrix[j,i] = difference_matrix[i,j]
        
        opportunity_matrix[i,j] = num_opportunities
        opportunity_matrix[j,i] = opportunity_matrix[i,j]
    
    return samples, difference_matrix, opportunity_matrix




if __name__=='__main__':


    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--debug", help="Loads only a subset of SNPs for speed", action="store_true")
    parser.add_argument("--chunk-size", type=int, help="max number of records to load", default=1000000000)
    parser.add_argument("--species", help="Name of specific species to run code on", default="all")
    args = parser.parse_args()

    debug = args.debug
    chunk_size = args.chunk_size
    species=args.species

    # Load subject and sample metadata
    sys.stderr.write("Loading sample metadata...\n")
    subject_sample_map = parse_HMP_data.parse_subject_sample_map()
    sys.stderr.write("Done!\n")
    
    # get a list of specis to run this script on. 
    good_species_list = parse_midas_data.parse_good_species_list()
    if debug:
        good_species_list = good_species_list[:3]
    elif species !='all':
        good_species_list = [species]

    # header for the output file.
    record_strs = [", ".join(['Species', 'Sample1', 'Sample2', 'Type', 'Num_changes', 'Num_opportunities'])]
    #good_species_list=['Bacteroides_vulgatus_57955']
    for species_name in good_species_list:

        sys.stderr.write("Loading haploid samples...\n")

        # Only plot samples above a certain depth threshold that are confidently phaseable.
        snp_samples = diversity_utils.calculate_haploid_samples(species_name, debug=debug)
        
        if len(snp_samples) < min_sample_size:
            sys.stderr.write("Not enough haploid samples!\n")
            continue
        
        sys.stderr.write("Proceeding with %d haploid samples!\n" % len(snp_samples))

        sys.stderr.write("Loading core genes...\n")
        core_genes = parse_midas_data.load_core_genes(species_name)
        sys.stderr.write("Done! Core genome consists of %d genes\n" % len(core_genes))

        # Analyze SNPs, looping over chunk sizes. 
        # Clunky, but necessary to limit memory usage on cluster

        # Load SNP information for species_name
        sys.stderr.write("Loading SNPs for %s...\n" % species_name)
        sys.stderr.write("(core genes only...)\n")
        pi_matrix_syn = numpy.array([])
        avg_pi_matrix_syn = numpy.array([])
    
        syn_difference_matrix = numpy.array([]) # 4d sites in core genes
        syn_opportunity_matrix = numpy.array([])
    
        non_difference_matrix = numpy.array([]) # 1d sites in core genes
        non_opportunity_matrix = numpy.array([])
    
        core_difference_matrix = numpy.array([]) # all sites in core genes
        core_opportunity_matrix = numpy.array([])
    
        snp_difference_matrix = numpy.array([]) # all sites in all genes
        snp_opportunity_matrix = numpy.array([])
    
        final_line_number = 0
        while final_line_number >= 0:
    
            sys.stderr.write("Loading chunk starting @ %d...\n" % final_line_number)
            dummy_samples, allele_counts_map, passed_sites_map, final_line_number = parse_midas_data.parse_snps(species_name, debug=debug, allowed_samples=snp_samples, chunk_size=chunk_size,initial_line_number=final_line_number)
            sys.stderr.write("Done! Loaded %d genes\n" % len(allele_counts_map.keys()))
            
            # Calculate fixation matrix
            sys.stderr.write("Calculating matrix of snp differences...\n")
            # Synonymous (4D)
            chunk_syn_difference_matrix, chunk_syn_opportunity_matrix = diversity_utils.calculate_fixation_matrix(allele_counts_map, passed_sites_map, min_change=min_change, allowed_genes=core_genes, allowed_variant_types=set(['4D']))
            # Nonsynonymous (1D) 
            chunk_non_difference_matrix, chunk_non_opportunity_matrix = diversity_utils.calculate_fixation_matrix(allele_counts_map, passed_sites_map, min_change=min_change, allowed_genes=core_genes, allowed_variant_types=set(['1D']))
            # Core (all)
            chunk_core_difference_matrix, chunk_core_opportunity_matrix = diversity_utils.calculate_fixation_matrix(allele_counts_map, passed_sites_map, min_change=min_change, allowed_genes=core_genes)
            # All
            chunk_snp_difference_matrix, chunk_snp_opportunity_matrix = diversity_utils.calculate_fixation_matrix(allele_counts_map, passed_sites_map, min_change=min_change)   # 
            sys.stderr.write("Done!\n")
    
            if snp_difference_matrix.shape[0]==0:
                snp_difference_matrix = numpy.zeros_like(chunk_snp_difference_matrix)*1.0
                snp_opportunity_matrix = numpy.zeros_like(snp_difference_matrix)*1.0
                syn_difference_matrix = numpy.zeros_like(snp_difference_matrix)*1.0
                syn_opportunity_matrix = numpy.zeros_like(snp_difference_matrix)*1.0
                non_difference_matrix = numpy.zeros_like(snp_difference_matrix)*1.0
                non_opportunity_matrix = numpy.zeros_like(snp_difference_matrix)*1.0
                core_difference_matrix = numpy.zeros_like(snp_difference_matrix)*1.0
                core_opportunity_matrix = numpy.zeros_like(snp_difference_matrix)*1.0
    
            # Add syn
            syn_difference_matrix += chunk_syn_difference_matrix
            syn_opportunity_matrix += chunk_syn_opportunity_matrix

            # Add non
            non_difference_matrix += chunk_non_difference_matrix
            non_opportunity_matrix += chunk_non_opportunity_matrix

            # Add core
            core_difference_matrix += chunk_core_difference_matrix
            core_opportunity_matrix += chunk_core_opportunity_matrix

            # Add all
            snp_difference_matrix += chunk_snp_difference_matrix
            snp_opportunity_matrix += chunk_snp_opportunity_matrix

            snp_samples = dummy_samples
    
    
        # Now calculate gene differences
        # Load gene coverage information for species_name
        sys.stderr.write("Loading pangenome data for %s...\n" % species_name)
        gene_samples, gene_names, gene_presence_matrix, gene_depth_matrix, marker_coverages,     gene_reads_matrix = parse_midas_data.parse_pangenome_data(species_name,allowed_samples=snp_samples)
        sys.stderr.write("Done!\n")
    
        gene_sample_list = list(gene_samples)
        gene_sample_set = set(gene_samples)
        
        # Calculate matrix of number of genes that differ
        sys.stderr.write("Calculating matrix of gene differences...\n")
        gene_difference_matrix, gene_opportunity_matrix =     gene_diversity_utils.calculate_coverage_based_gene_hamming_matrix(gene_reads_matrix, gene_depth_matrix, marker_coverages)

        good_marker_coverages = (marker_coverages>=min_coverage)

        gene_difference_matrix = gene_difference_matrix*good_marker_coverages[:,None]*good_marker_coverages[None,:]
        gene_opportunity_matrix = gene_opportunity_matrix*good_marker_coverages[:,None]*good_marker_coverages[None,:]
    
        # Add records to output
        # Calculate which pairs of idxs belong to the same sample, which to the same subject
        # and which to different subjects
        same_sample_idxs, same_subject_idxs, diff_subject_idxs = parse_midas_data.calculate_subject_pairs(subject_sample_map, snp_samples)

        for idxs in [same_subject_idxs, diff_subject_idxs]:
            for sample_pair_idx in xrange(0,len(idxs[0])):
   
                i = idxs[0][sample_pair_idx]
                j = idxs[1][sample_pair_idx]
        
                sample_i = snp_samples[i]
                sample_j = snp_samples[j]
        
                record_str_items = [species_name, sample_i, sample_j, '4D', str(syn_difference_matrix[i,j]), str(syn_opportunity_matrix[i,j])]
                record_strs.append( ", ".join(record_str_items) )
        
                record_str_items = [species_name, sample_i, sample_j, '1D', str(non_difference_matrix[i,j]), str(non_opportunity_matrix[i,j])]
                record_strs.append( ", ".join(record_str_items) )
        
                record_str_items = [species_name, sample_i, sample_j, 'core', str(core_difference_matrix[i,j]), str(core_opportunity_matrix[i,j])]
                record_strs.append( ", ".join(record_str_items) )
        
                record_str_items = [species_name, sample_i, sample_j, 'all', str(snp_difference_matrix[i,j]), str(snp_opportunity_matrix[i,j])]
                record_strs.append( ", ".join(record_str_items) )
        
                if (sample_i in gene_sample_set) and (sample_j in gene_sample_set):
                    gene_i = gene_sample_list.index(sample_i)
                    gene_j = gene_sample_list.index(sample_j)
            
                    record_str_items = [species_name, sample_i, sample_j, 'genes', str(gene_difference_matrix[gene_i, gene_j]), str(gene_opportunity_matrix[gene_i, gene_j])]
                else:
                    record_str_items = [species_name, sample_i, sample_j, 'genes', '0', '0']
                
                record_strs.append( ", ".join(record_str_items) )
        
        sys.stderr.write("Done with %s!\n" % species_name) 
    
    sys.stderr.write("Done looping over species!\n")
    
    sys.stderr.write("Writing intermediate file...\n")
    file = open(intermediate_filename,"w")
    record_str = "\n".join(record_strs)
    file.write(record_str)
    file.close()
    sys.stderr.write("Done!\n")

 
