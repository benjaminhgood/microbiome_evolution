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

private_snv_directory = '%sprivate_snvs/' % (parse_midas_data.data_directory)
intermediate_filename_template = '%s%s.txt.gz'  

min_coverage = config.min_median_coverage
min_sample_size = 5
allowed_variant_types = set(['1D','2D','3D','4D'])


def load_private_snv_map(species_name):
# This definition is called whenever another script downstream uses the output of this data.

    intermediate_filename = intermediate_filename_template % (private_snv_directory, species_name)

    private_snv_map = {}

    if not os.path.isfile(intermediate_filename):
        return private_snv_map
    
    file = gzip.open(intermediate_filename,"r")
    file.readline() # header
    for line in file:
     
        items = line.split(",")
        
        contig = items[0].strip()
        location = long(items[1])
        gene_name = items[2].strip()
        variant_type = items[3].strip()
        host = items[4].strip()
        
        private_snv_map[(contig, location)] = (gene_name, variant_type, host)
        
    return private_snv_map


if __name__=='__main__':


    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--debug", help="Loads only a subset of SNPs for speed", action="store_true")
    parser.add_argument("--chunk-size", type=int, help="max number of records to load", default=1000000000)
    parser.add_argument("species", help="Name of specific species to run code on", default="all")
    args = parser.parse_args()

    debug = args.debug
    chunk_size = args.chunk_size
    species_name=args.species
    good_species_list = [species_name]
    
    # Load subject and sample metadata
    sys.stderr.write("Loading sample metadata...\n")
    subject_sample_map = sample_utils.parse_subject_sample_map()
    sys.stderr.write("Done!\n")
    
    os.system('mkdir -p %s' % private_snv_directory)

    for species_name in good_species_list:

        sys.stderr.write("Loading haploid samples...\n")

        # Only plot samples above a certain depth threshold that are confidently phaseable.
        snp_samples = diversity_utils.calculate_haploid_samples(species_name, debug=debug)
        
        if len(snp_samples) < min_sample_size:
            sys.stderr.write("Not enough haploid samples!\n")
            continue
                
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
        
        private_snvs = []
        
        final_line_number = 0
        while final_line_number >= 0:
    
            sys.stderr.write("Loading chunk starting @ %d...\n" % final_line_number)
            dummy_samples, allele_counts_map, passed_sites_map, final_line_number = parse_midas_data.parse_snps(species_name, debug=debug, allowed_samples=snp_samples, chunk_size=chunk_size,initial_line_number=final_line_number, allowed_genes=core_genes)
            snp_samples = dummy_samples
            sys.stderr.write("Done! Loaded %d genes\n" % len(allele_counts_map.keys()))
            
            chunk_private_snvs = diversity_utils.calculate_private_snvs(snp_samples, allele_counts_map, passed_sites_map)
            
            private_snvs.extend(chunk_private_snvs)
        
        if len(private_snvs)>0:
            
            intermediate_filename = intermediate_filename_template % (private_snv_directory, species_name)
            
            # Now add records
            output_file = gzip.open(intermediate_filename,"w")
            # Header
            output_file.write("contig, location, gene_name, var_type, host_id\n")
            for contig, location, gene_name, variant_type, host in private_snvs:
            
                record_str_items = [contig, str(location), gene_name, variant_type, host]           
                record_str = ", ".join(record_str_items)
                output_file.write(record_str)
                output_file.write("\n")
            
            output_file.close()

        sys.stderr.write("Done with %s!\n" % species_name) 
    
    sys.stderr.write("Done looping over species!\n")
    
    sys.stderr.write("Testing loading...\n")
    private_snv_map = load_private_snv_map(good_species_list[0])
    sys.stderr.write("Done!\n")

 
