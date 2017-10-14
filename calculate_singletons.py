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

intermediate_filename = '%ssingleton_rates.txt' % (parse_midas_data.data_directory)
    

min_coverage = config.min_median_coverage
alpha = 0.5 # Confidence interval range for rate estimates
low_pi_threshold = 1e-03
low_divergence_threshold = 1e-03
min_change = 0.8
min_sample_size = 10
allowed_variant_types = set(['1D','2D','3D','4D'])


def load_singleton_rate_map(species_name):
# This definition is called whenever another script downstream uses the output of this data.

    singleton_rate_map = {}

    file = open(intermediate_filename,"r")
    file.readline() # header
    for line in file:
        items = line.split(",")
        if items[0].strip()!=species_name:
            continue
            
        sample = items[1].strip()
        type = items[2].strip()
        num_singletons = float(items[3])
        num_opportunities = float(items[4])
        
        if type not in singleton_rate_map:
            singleton_rate_map[type] = {}
          
        singleton_rate_map[type][sample] = (num_singletons, num_opportunities)
        
    return singleton_rate_map

def calculate_matrices_from_singleton_rate_map(singleton_rate_map, type, allowed_samples=[]):
# once the map is loaded, then we can compute rate matrices in this definition (so, it relies on the previous def)    

    sample_set = set([])
    for sample in singleton_rate_map[type].keys():
        sample_set.add(sample)
    
    if len(allowed_samples)>0:
        allowed_sample_set = set(allowed_samples)    
    else:
        allowed_sample_set = sample_set
        
    sample_set = sample_set & allowed_sample_set
    samples = list(sorted(sample_set))
    
    sample_idx_map = {samples[i]: i for i in xrange(0,len(samples))}
    
    singleton_matrix = numpy.zeros(len(samples))*1.0
    opportunity_matrix = numpy.zeros_like(singleton_matrix)
    
    for sample in singleton_rate_map[type].keys():
        
        if not (sample in sample_set):
            continue
        
        i = sample_idx_map[sample]
        
        num_singletons, num_opportunities = singleton_rate_map[type][sample]
        
        singleton_matrix[i] = num_singletons
        opportunity_matrix[i] = num_opportunities
     
    #print singleton_matrix, opportunity_matrix    
    return samples, singleton_matrix, opportunity_matrix

    
def calculate_matrices_from_substitution_rate_map(substitution_rate_map, type, allowed_samples=[]):
# once the map is loaded, then we can compute rate matrices in this definition (so, it relies on the previous def)    

    samples, mut_difference_matrix, rev_difference_matrix, mut_opportunity_matrix, rev_opportunity_matrix = calculate_mutrev_matrices_from_substitution_rate_map( substitution_rate_map, type, allowed_samples)

    difference_matrix = mut_difference_matrix+rev_difference_matrix
    opportunity_matrix = mut_opportunity_matrix+rev_opportunity_matrix
    
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
    record_strs = [", ".join(['Species', 'Sample1', 'Sample2', 'Type', 'Num_muts', 'Num_revs', 'Num_mut_opportunities', 'Num_rev_opportunities'])]
    #good_species_list=['Bacteroides_vulgatus_57955']
    for species_name in good_species_list:

        sys.stderr.write("Loading haploid samples...\n")

        # Only plot samples above a certain depth threshold that are confidently phaseable.
        snp_samples = diversity_utils.calculate_haploid_samples(species_name, debug=debug)
        
        if len(snp_samples) < min_sample_size:
            sys.stderr.write("Not enough haploid samples!\n")
            continue
        
        # Only consider one sample per person
        snp_samples = snp_samples[parse_midas_data.calculate_unique_samples(subject_sample_map, sample_list=snp_samples)]
        
        sys.stderr.write("Proceeding with %d haploid samples!\n" % len(snp_samples))

        sys.stderr.write("Loading core genes...\n")
        core_genes = parse_midas_data.load_core_genes(species_name)
        sys.stderr.write("Done! Core genome consists of %d genes\n" % len(core_genes))

        # Analyze SNPs, looping over chunk sizes. 
        # Clunky, but necessary to limit memory usage on cluster

        # Load SNP information for species_name
        sys.stderr.write("Loading SNPs for %s...\n" % species_name)
        sys.stderr.write("(core genes only...)\n")
        
        snp_singleton_count_matrix = numpy.array([]) 
        snp_singleton_opportunity_matrix = numpy.array([])
                
        syn_singleton_count_matrix = numpy.array([])
        syn_singleton_opportunity_matrix = numpy.array([])
                
        non_singleton_count_matrix = numpy.array([])
        non_singleton_opportunity_matrix = numpy.array([])
                
        core_singleton_count_matrix = numpy.array([])
        core_singleton_opportunity_matrix = numpy.array([])
                
        final_line_number = 0
        while final_line_number >= 0:
    
            sys.stderr.write("Loading chunk starting @ %d...\n" % final_line_number)
            dummy_samples, allele_counts_map, passed_sites_map, final_line_number = parse_midas_data.parse_snps(species_name, debug=debug, allowed_samples=snp_samples, chunk_size=chunk_size,initial_line_number=final_line_number)
            sys.stderr.write("Done! Loaded %d genes\n" % len(allele_counts_map.keys()))
            
            # Calculate fixation matrix
            sys.stderr.write("Calculating matrix of singletons...\n")
            # Synonymous (4D)
            
            chunk_syn_singleton_count_matrix, chunk_syn_singleton_opportunity_matrix = diversity_utils.calculate_singleton_matrix(allele_counts_map, passed_sites_map, allowed_genes=core_genes, allowed_variant_types=set(['4D']))
            
            # Nonsynonymous (1D) 
            chunk_non_singleton_count_matrix, chunk_non_singleton_opportunity_matrix = diversity_utils.calculate_singleton_matrix(allele_counts_map, passed_sites_map, allowed_genes=core_genes, allowed_variant_types=set(['1D']))
            
            # Core (all)
            chunk_core_singleton_count_matrix, chunk_core_singleton_opportunity_matrix = diversity_utils.calculate_singleton_matrix(allele_counts_map, passed_sites_map, allowed_genes=core_genes)
            
            # All
            chunk_snp_singleton_count_matrix, chunk_snp_singleton_opportunity_matrix = diversity_utils.calculate_singleton_matrix(allele_counts_map, passed_sites_map)
            
            sys.stderr.write("Done!\n")
    
            if len(snp_singleton_count_matrix)==0:
                snp_singleton_count_matrix = numpy.zeros_like(chunk_snp_singleton_count_matrix)*1.0
                snp_singleton_opportunity_matrix = numpy.zeros_like(chunk_snp_singleton_count_matrix)*1.0
                
                syn_singleton_count_matrix = numpy.zeros_like(chunk_snp_singleton_count_matrix)*1.0
                syn_singleton_opportunity_matrix = numpy.zeros_like(chunk_snp_singleton_count_matrix)*1.0
                
                non_singleton_count_matrix = numpy.zeros_like(chunk_snp_singleton_count_matrix)*1.0
                non_singleton_opportunity_matrix = numpy.zeros_like(chunk_snp_singleton_count_matrix)*1.0
                
                core_singleton_count_matrix = numpy.zeros_like(chunk_snp_singleton_count_matrix)*1.0
                core_singleton_opportunity_matrix = numpy.zeros_like(chunk_snp_singleton_count_matrix)*1.0
                
            # Add syn
            syn_singleton_count_matrix += chunk_syn_singleton_count_matrix
            syn_singleton_opportunity_matrix += chunk_syn_singleton_opportunity_matrix
            
            # Add non
            non_singleton_count_matrix += chunk_non_singleton_count_matrix
            non_singleton_opportunity_matrix += chunk_non_singleton_opportunity_matrix
            
            # Add core
            core_singleton_count_matrix += chunk_core_singleton_count_matrix
            core_singleton_opportunity_matrix += chunk_core_singleton_opportunity_matrix
            
            # Add all
            snp_singleton_count_matrix += chunk_snp_singleton_count_matrix
            snp_singleton_opportunity_matrix += chunk_snp_singleton_opportunity_matrix
            
            snp_samples = dummy_samples
    
        # Add records to output
        for i in xrange(0,len(snp_samples)):
            
            sample_i = snp_samples[i]
            
            record_str_items = [species_name, sample_i, '4D', str(syn_singleton_count_matrix[i]), str(syn_singleton_opportunity_matrix[i])]
            record_strs.append( ", ".join(record_str_items) )
        
            record_str_items = [species_name, sample_i, '1D', str(non_singleton_count_matrix[i]), str(non_singleton_opportunity_matrix[i])]
            record_strs.append( ", ".join(record_str_items) )
        
            record_str_items = [species_name, sample_i, 'core', str(core_singleton_count_matrix[i]), str(core_singleton_opportunity_matrix[i])]
            record_strs.append( ", ".join(record_str_items) )
        
            record_str_items = [species_name, sample_i, 'all', str(snp_singleton_count_matrix[i]), str(snp_singleton_opportunity_matrix[i])]
            record_strs.append( ", ".join(record_str_items) )
        

        sys.stderr.write("Done with %s!\n" % species_name) 
    
    sys.stderr.write("Done looping over species!\n")
    
    sys.stderr.write("Writing intermediate file...\n")
    file = open(intermediate_filename,"w")
    record_str = "\n".join(record_strs)
    file.write(record_str)
    file.close()
    sys.stderr.write("Done!\n")

 
