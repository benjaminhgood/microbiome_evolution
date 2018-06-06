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

import core_gene_utils
import gzip
import calculate_substitution_rates
import clade_utils

private_snv_directory = '%ssnv_distances/' % (parse_midas_data.data_directory)
intermediate_filename_template = '%s%s.txt.gz'  

min_coverage = config.min_median_coverage
alpha = 0.5 # Confidence interval range for rate estimates
low_pi_threshold = 1e-03
low_divergence_threshold = 1e-03
min_change = 0.8
min_sample_size = 10
allowed_variant_types = set(['1D','2D','3D','4D'])

def load_snv_distance_map(species_name):
# This definition is called whenever another script downstream uses the output of this data.

    intermediate_filename = intermediate_filename_template % (private_snv_directory, species_name)

    snv_distance_map = {}

    file = gzip.open(intermediate_filename,"r")
    file.readline() # header
    for line in file:
     
        items = line.split(",")
        
        contig = items[0].strip()
        location = long(items[1])
        variant_type = items[2].strip()
        derived_allele_count = long(items[3])
        ancestral_allele_count = long(items[4])
        min_between_d = float(items[5])
        max_within_d1 = float(items[6])
        max_within_d2 = float(items[7])
        
        snv_distance_map[(contig, location)] = (variant_type, derived_allele_count, ancestral_allele_count, min_between_d, max_within_d1, max_within_d2)
        
    return snv_distance_map




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

    os.system('mkdir -p %s' % private_snv_directory)

    for species_name in good_species_list:

        

        # Only plot samples above a certain depth threshold that are "haploids"
        snp_samples = diversity_utils.calculate_haploid_samples(species_name, debug=debug)
        # Only consider one sample per person
        snp_samples =     snp_samples[parse_midas_data.calculate_unique_samples(subject_sample_map, sample_list=snp_samples)]
        sys.stderr.write("Proceeding with %d haploid samples!\n" % len(snp_samples))
        
        if len(snp_samples) < min_sample_size:
            sys.stderr.write("Not enough haploid samples!\n")
            continue
                
        sys.stderr.write("Proceeding with %d haploid samples!\n" % len(snp_samples))

        sys.stderr.write("Loading core genes...\n")
        core_genes = core_gene_utils.parse_core_genes(species_name)
        non_shared_genes = core_gene_utils.parse_non_shared_reference_genes(species_name)
        shared_pangenome_genes = core_gene_utils.parse_shared_genes(species_name)
        sys.stderr.write("Done! Core genome consists of %d genes\n" % len(core_genes))
        sys.stderr.write("%d shared genes and %d non-shared genes\n" % (len(shared_pangenome_genes), len(non_shared_genes)))

        sys.stderr.write("Loading pre-computed substitution rates for %s...\n" % species_name)
        substitution_rate_map = calculate_substitution_rates.load_substitution_rate_map(species_name)
        sys.stderr.write("Calculating matrix...\n")
        dummy_samples, snp_difference_matrix, snp_opportunity_matrix = calculate_substitution_rates.calculate_matrices_from_substitution_rate_map(substitution_rate_map, 'core', allowed_samples=snp_samples)
        snp_samples = dummy_samples
        sys.stderr.write("Done!\n")

        snp_substitution_rate = snp_difference_matrix*1.0/(snp_opportunity_matrix+(snp_opportunity_matrix==0))

        # Analyze SNPs, looping over chunk sizes. 
        # Clunky, but necessary to limit memory usage on cluster

        # Load SNP information for species_name
        sys.stderr.write("Loading SNPs for %s...\n" % species_name)
        sys.stderr.write("(core genes only...)\n")
        snp_data = []
    
        final_line_number = 0
        while final_line_number >= 0:
    
            sys.stderr.write("Loading chunk starting @ %d...\n" % final_line_number)
            dummy_samples, allele_counts_map, passed_sites_map, final_line_number = parse_midas_data.parse_snps(species_name, debug=debug, allowed_samples=snp_samples, chunk_size=chunk_size,initial_line_number=final_line_number, allowed_genes=core_genes)
            sys.stderr.write("Done! Loaded %d genes\n" % len(allele_counts_map.keys()))
        
            if not (dummy_samples==snp_samples).all():
                sys.stderr.write("Ordering problem!\n")
            
            # Calculate fixation matrix
            sys.stderr.write("Calculating snp distances...\n")
            chunk_snp_data = clade_utils.calculate_snp_distances(allele_counts_map, passed_sites_map, snp_substitution_rate)
            
            sys.stderr.write("Done!\n")
    
            snp_data.extend(chunk_snp_data)
        
        sys.stderr.write("Done!\n")
     
        if len(snp_data)>0:
            
            intermediate_filename = intermediate_filename_template % (private_snv_directory, species_name)
            
            # Now add records
            output_file = gzip.open(intermediate_filename,"w")
            # Header
            output_file.write("contig, location, var_type, derived_allele_count, ancestral_allele_count, min_between_d, max_within_derived_d, max_within_ancestral_d\n")
            for location_tuple, variant_type, derived_allele_count, ancestral_allele_count, min_between_d, max_within_d1, avg_within_d1, max_within_d2, avg_within_d2 in snp_data:
            
                contig, location = location_tuple
                
                record_str_items = [contig, str(location), variant_type, str(derived_allele_count), str(ancestral_allele_count), str(min_between_d), str(max_within_d1), str(max_within_d2)]           
                record_str = ", ".join(record_str_items)
                output_file.write(record_str)
                output_file.write("\n")
            
            output_file.close()

        sys.stderr.write("Done with %s!\n" % species_name) 
    
    sys.stderr.write("Done looping over species!\n")
    
    sys.stderr.write("Testing loading...\n")
    snv_distance_map = load_snv_distance_map(good_species_list[0])
    sys.stderr.write("Done!\n")

 
