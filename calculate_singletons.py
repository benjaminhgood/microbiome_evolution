import sample_utils
import config
import parse_midas_data
import os.path
import pylab
import sys
import numpy

import diversity_utils
import gene_diversity_utils

import stats_utils
from math import log10,ceil
from numpy.random import randint

import core_gene_utils
import gzip

singleton_directory = '%ssingleton_rates/' % (parse_midas_data.data_directory)
intermediate_filename_template = '%s%s.txt.gz'  
   

min_coverage = config.min_median_coverage
min_sample_size = 10

def load_singleton_rate_map(species_name):
# This definition is called whenever another script downstream uses the output of this data.

    intermediate_filename = intermediate_filename_template % (singleton_directory, species_name)

    singleton_rate_map = {}

    if not os.path.isfile(intermediate_filename):
        return singleton_rate_map
    
    file = gzip.open(intermediate_filename,"r")
    file.readline() # header
    for line in file:
        items = line.split(",")
        if items[0].strip()!=species_name:
            continue
            
        sample_i = items[1].strip()
        sample_j = items[2].strip()
        type = items[3].strip()
        num_singletons = float(items[4])
        num_doubletons = float(items[5])
        num_differences = float(items[6])
        num_opportunities = float(items[7])
        
        if type not in singleton_rate_map:
            singleton_rate_map[type] = {}
        
        if sample_i==sample_j:
            num_singletons = 0
            num_doubletons = 0
            num_differences = 0
          
        singleton_rate_map[type][sample_i, sample_j] = (num_singletons, num_doubletons, num_differences, num_opportunities)
        
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
    
    sample_set = set()
    for sample_i, sample_j in singleton_rate_map[type]:
        sample_set.add(sample_i)
        sample_set.add(sample_j)    
         
    if len(allowed_samples)==0:
        allowed_samples = list(sorted(allowed_sample_set))
    
    
    samples = []
    # preserve same order as allowed samples
    for sample in allowed_samples:
        if sample in sample_set:
            samples.append(sample)
               
    singleton_matrix = numpy.zeros((len(samples),len(samples)))*1.0
    doubleton_matrix = numpy.zeros_like(singleton_matrix)
    difference_matrix = numpy.zeros_like(singleton_matrix)
    opportunity_matrix = numpy.zeros_like(singleton_matrix)
    
    for i in xrange(0,len(samples)):
        for j in xrange(0,len(samples)):
            
            num_singletons, num_doubletons, num_differences, num_opportunities = singleton_rate_map[type][(samples[i], samples[j])]
        
            if i==j:
                num_doubletons = 0
        
            singleton_matrix[i,j] = num_singletons
            doubleton_matrix[i,j] = num_doubletons
            difference_matrix[i,j] = num_differences
            opportunity_matrix[i,j] = num_opportunities
     
    #print singleton_matrix, opportunity_matrix    
    return samples, singleton_matrix, doubleton_matrix, difference_matrix, opportunity_matrix

    
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
    parser.add_argument("species", help="Name of specific species to run code on")
    args = parser.parse_args()

    debug = args.debug
    chunk_size = args.chunk_size
    species_name=args.species
    good_species_list = [species_name]
    
    os.system('mkdir -p %s' % singleton_directory)

    # Load subject and sample metadata
    sys.stderr.write("Loading sample metadata...\n")
    subject_sample_map = sample_utils.parse_subject_sample_map()
    sys.stderr.write("Done!\n")
    
    # header for the output file.
    record_strs = [", ".join(['Species', 'Sample1', 'Sample2', 'Type', 'Num_muts', 'Num_revs', 'Num_mut_opportunities', 'Num_rev_opportunities'])]
    for species_name in good_species_list:

        sys.stderr.write("Loading haploid samples...\n")

        # Only plot samples above a certain depth threshold that are confidently phaseable.
        snp_samples = diversity_utils.calculate_haploid_samples(species_name, debug=debug)
        
        if len(snp_samples) < min_sample_size:
            sys.stderr.write("Not enough haploid samples!\n")
            continue
        
        # Only consider one sample per person
        snp_samples = snp_samples[sample_utils.calculate_unique_samples(subject_sample_map, sample_list=snp_samples)]
        
        sys.stderr.write("Proceeding with %d haploid samples!\n" % len(snp_samples))

        sys.stderr.write("Loading whitelisted genes...\n")
        core_genes = core_gene_utils.parse_core_genes(species_name)
        non_shared_genes = core_gene_utils.parse_non_shared_reference_genes(species_name)
        shared_pangenome_genes = core_gene_utils.parse_shared_genes(species_name)
        sys.stderr.write("Done! %d core genes and %d shared genes and %d non-shared genes\n" % (len(core_genes), len(shared_pangenome_genes), len(non_shared_genes)))


        # Analyze SNPs, looping over chunk sizes. 
        # Clunky, but necessary to limit memory usage on cluster

        # Load SNP information for species_name
        sys.stderr.write("Loading SNPs for %s...\n" % species_name)
        sys.stderr.write("(core genes only...)\n")
        
        snp_doubleton_count_matrix = numpy.array([])
        snp_difference_count_matrix = numpy.array([])         
        snp_singleton_count_matrix = numpy.array([]) 
        snp_singleton_opportunity_matrix = numpy.array([])
        
        syn_doubleton_count_matrix = numpy.array([])
        syn_difference_count_matrix = numpy.array([])                 
        syn_singleton_count_matrix = numpy.array([])
        syn_singleton_opportunity_matrix = numpy.array([])
        
        non_doubleton_count_matrix = numpy.array([])
        non_difference_count_matrix = numpy.array([])                 
        non_singleton_count_matrix = numpy.array([])
        non_singleton_opportunity_matrix = numpy.array([])
        
        core_doubleton_count_matrix = numpy.array([])
        core_difference_count_matrix = numpy.array([])             
        core_singleton_count_matrix = numpy.array([])
        core_singleton_opportunity_matrix = numpy.array([])
                
        final_line_number = 0
        while final_line_number >= 0:
    
            sys.stderr.write("Loading chunk starting @ %d...\n" % final_line_number)
            dummy_samples, allele_counts_map, passed_sites_map, final_line_number = parse_midas_data.parse_snps(species_name, debug=debug, allowed_samples=snp_samples, chunk_size=chunk_size,initial_line_number=final_line_number, allowed_genes=non_shared_genes)
            sys.stderr.write("Done! Loaded %d genes\n" % len(allele_counts_map.keys()))
            
            # Calculate fixation matrix
            sys.stderr.write("Calculating matrix of singletons...\n")
            # Synonymous (4D)
            
            chunk_syn_doubleton_count_matrix, chunk_syn_singleton_count_matrix, chunk_syn_difference_count_matrix, chunk_syn_singleton_opportunity_matrix = diversity_utils.calculate_singleton_matrix(allele_counts_map, passed_sites_map, allowed_genes=core_genes, allowed_variant_types=set(['4D']))
            
            # Nonsynonymous (1D) 
            chunk_non_doubleton_count_matrix, chunk_non_singleton_count_matrix, chunk_non_difference_count_matrix, chunk_non_singleton_opportunity_matrix = diversity_utils.calculate_singleton_matrix(allele_counts_map, passed_sites_map, allowed_genes=core_genes, allowed_variant_types=set(['1D']))
            
            # Core (all)
            chunk_core_doubleton_count_matrix, chunk_core_singleton_count_matrix, chunk_core_difference_count_matrix, chunk_core_singleton_opportunity_matrix = diversity_utils.calculate_singleton_matrix(allele_counts_map, passed_sites_map, allowed_genes=core_genes)
            
            # All
            chunk_snp_doubleton_count_matrix, chunk_snp_singleton_count_matrix, chunk_snp_difference_count_matrix, chunk_snp_singleton_opportunity_matrix = diversity_utils.calculate_singleton_matrix(allele_counts_map, passed_sites_map)
            
            sys.stderr.write("Done!\n")
    
            if len(snp_singleton_count_matrix)==0:
                snp_doubleton_count_matrix = numpy.zeros_like(chunk_snp_singleton_count_matrix)*1.0
                snp_difference_count_matrix = numpy.zeros_like(chunk_snp_singleton_count_matrix)*1.0
                
                snp_singleton_count_matrix = numpy.zeros_like(chunk_snp_singleton_count_matrix)*1.0
                snp_singleton_opportunity_matrix = numpy.zeros_like(chunk_snp_singleton_count_matrix)*1.0
                
                syn_doubleton_count_matrix = numpy.zeros_like(chunk_snp_singleton_count_matrix)*1.0
                syn_difference_count_matrix = numpy.zeros_like(chunk_snp_singleton_count_matrix)*1.0
                
                syn_singleton_count_matrix = numpy.zeros_like(chunk_snp_singleton_count_matrix)*1.0
                syn_singleton_opportunity_matrix = numpy.zeros_like(chunk_snp_singleton_count_matrix)*1.0
                
                non_doubleton_count_matrix = numpy.zeros_like(chunk_snp_singleton_count_matrix)*1.0
                non_difference_count_matrix = numpy.zeros_like(chunk_snp_singleton_count_matrix)*1.0
                
                non_singleton_count_matrix = numpy.zeros_like(chunk_snp_singleton_count_matrix)*1.0
                non_singleton_opportunity_matrix = numpy.zeros_like(chunk_snp_singleton_count_matrix)*1.0
                
                core_doubleton_count_matrix = numpy.zeros_like(chunk_snp_singleton_count_matrix)*1.0
                core_difference_count_matrix = numpy.zeros_like(chunk_snp_singleton_count_matrix)*1.0
                
                core_singleton_count_matrix = numpy.zeros_like(chunk_snp_singleton_count_matrix)*1.0
                core_singleton_opportunity_matrix = numpy.zeros_like(chunk_snp_singleton_count_matrix)*1.0
                
            # Add syn
            syn_singleton_count_matrix += chunk_syn_singleton_count_matrix
            syn_singleton_opportunity_matrix += chunk_syn_singleton_opportunity_matrix
            syn_doubleton_count_matrix += chunk_syn_doubleton_count_matrix
            syn_difference_count_matrix += chunk_syn_difference_count_matrix
            
            # Add non
            non_singleton_count_matrix += chunk_non_singleton_count_matrix
            non_singleton_opportunity_matrix += chunk_non_singleton_opportunity_matrix
            non_doubleton_count_matrix += chunk_non_doubleton_count_matrix
            non_difference_count_matrix += chunk_non_difference_count_matrix
            
            # Add core
            core_singleton_count_matrix += chunk_core_singleton_count_matrix
            core_singleton_opportunity_matrix += chunk_core_singleton_opportunity_matrix
            core_doubleton_count_matrix += chunk_core_doubleton_count_matrix
            core_difference_count_matrix += chunk_core_difference_count_matrix
            
            # Add all
            snp_singleton_count_matrix += chunk_snp_singleton_count_matrix
            snp_singleton_opportunity_matrix += chunk_snp_singleton_opportunity_matrix
            snp_doubleton_count_matrix += chunk_snp_doubleton_count_matrix
            snp_difference_count_matrix += chunk_snp_difference_count_matrix
            
            snp_samples = dummy_samples
    
        # Add records to output
        for i in xrange(0,len(snp_samples)):
            for j in xrange(0,len(snp_samples)):
            
                sample_i = snp_samples[i]
                sample_j = snp_samples[j]
            
                record_str_items = [species_name, sample_i, sample_j, '4D', str(syn_singleton_count_matrix[i,j]), str(syn_doubleton_count_matrix[i,j]), str(syn_difference_count_matrix[i,j]), str(syn_singleton_opportunity_matrix[i,j])]
                record_strs.append( ", ".join(record_str_items) )
        
                record_str_items = [species_name, sample_i, sample_j, '1D', str(non_singleton_count_matrix[i,j]), str(non_doubleton_count_matrix[i,j]), str(non_difference_count_matrix[i,j]), str(non_singleton_opportunity_matrix[i,j])]
                record_strs.append( ", ".join(record_str_items) )
        
                record_str_items = [species_name, sample_i, sample_j, 'core', str(core_singleton_count_matrix[i,j]), str(core_doubleton_count_matrix[i,j]), str(core_difference_count_matrix[i,j]), str(core_singleton_opportunity_matrix[i,j])]
                record_strs.append( ", ".join(record_str_items) )
        
                record_str_items = [species_name, sample_i, sample_j, 'all', str(snp_singleton_count_matrix[i,j]), str(snp_doubleton_count_matrix[i,j]), str(snp_difference_count_matrix[i,j]), str(snp_singleton_opportunity_matrix[i,j])]
                record_strs.append( ", ".join(record_str_items) )
        

        sys.stderr.write("Done with %s!\n" % species_name) 
    
        sys.stderr.write("Writing intermediate file...\n")
        intermediate_filename = intermediate_filename_template % (singleton_directory, species_name)
        file = gzip.open(intermediate_filename,"w")
        record_str = "\n".join(record_strs)
        file.write(record_str)
        file.close()
        sys.stderr.write("Done!\n")
    
        sys.stderr.write("Testing loading...\n")
        singleton_rate_map = load_singleton_rate_map(species_name)
        sys.stderr.write("Done!\n")

 
