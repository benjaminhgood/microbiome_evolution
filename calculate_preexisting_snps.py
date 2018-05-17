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
import gzip

intermediate_filename = '%spreexisting_snps.txt.gz' % (parse_midas_data.data_directory)
    

min_coverage = 5
allowed_variant_types = set(['1D','2D','3D','4D'])

def parse_preexisting_snps(species_name):
    
    preexisting_snps = {}  
    file = gzip.GzipFile(intermediate_filename,"r")
    for line in file:
        if line.startswith(species_name):
            contig_items = line.split(";")[1:]
            for contig_item in contig_items:
                contig_item = contig_item.strip()
                if contig_item=="":
                    continue
                contig_subitems = contig_item.split(":")
                contig = contig_subitems[0].strip()
                snp_items = contig_subitems[1].split()
                for snp_item in snp_items:
                    snp_subitems = snp_item.split(",")
                    location = long(snp_subitems[0])
                    prevalence = float(snp_subitems[1])
                    
                    preexisting_snps[(contig,location)] = prevalence
                    
    file.close()
    
    return preexisting_snps
    
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

    file = gzip.GzipFile(intermediate_filename,"w")

    # Load subject and sample metadata
    sys.stderr.write("Loading sample metadata...\n")
    subject_sample_map = parse_HMP_data.parse_subject_sample_map()
    sys.stderr.write("Done!\n")
    
    # get a list of specis to run this script on. 
    good_species_list = parse_midas_data.parse_good_species_list()
    if species!='all':
        good_species_list = [species]
    else:    
        if debug:
            good_species_list = good_species_list[:3]
    
    # header for the output file.
    record_strs = []
    
    for species_name in good_species_list:

        sys.stderr.write("Loading samples...\n")

        # Only plot samples above a certain depth threshold that are confidently phaseable.
        snp_samples = diversity_utils.calculate_highcoverage_samples(species_name, min_coverage=min_coverage)
        
        if len(snp_samples)<2:
            continue
            
        sys.stderr.write("found %d samples\n" % len(snp_samples))
        
        # Analyze SNPs, looping over chunk sizes. 
        # Clunky, but necessary to limit memory usage on cluster

        # Load SNP information for species_name
        sys.stderr.write("Loading SNPs for %s...\n" % species_name)
        
        snps = []
        snp_map = {} # contig: list of locations map
        
        final_line_number = 0
        while final_line_number >= 0:
    
            sys.stderr.write("Loading chunk starting @ %d...\n" % final_line_number)
            dummy_samples, allele_counts_map, passed_sites_map, final_line_number = parse_midas_data.parse_snps(species_name, debug=debug, chunk_size=chunk_size,initial_line_number=final_line_number,allowed_samples=snp_samples)
            sys.stderr.write("Done! Loaded %d genes\n" % len(allele_counts_map.keys()))
            snp_samples = dummy_samples
            
            # Calculate fixation matrix
            sys.stderr.write("Calculating matrix of snp differences...\n")
            
            chunk_snps = diversity_utils.calculate_preexisting_snps(allele_counts_map, passed_sites_map, allowed_variant_types=set(['1D','2D','3D','4D']))
            
            sys.stderr.write("Done!\n")
            for contig,location,prevalence in chunk_snps:
                if contig not in snp_map:
                    snp_map[contig] = {}
                snp_map[contig][location] = prevalence 
                
            #snps.extend(chunk_snps)
        
        contig_strs = [species_name]
        for contig in snp_map.keys():
            contig_str = contig + ": "+(" ".join(["%d,%0.2f" % (location, snp_map[contig][location]) for location in sorted(snp_map[contig])]))
            contig_strs.append(contig_str)
        
        record_str = "; ".join(contig_strs)
            
        #record_str_items = [species_name]+["%s,%d" % (contig, location) for contig, location in snps]   
        #record_str = " ".join(record_str_items)
        file.write(record_str)
        file.write("\n")
            
        sys.stderr.write("Done with %s!\n" % species_name) 
    
    sys.stderr.write("Done looping over species!\n")
    
    file.close()
    sys.stderr.write("Done!\n")