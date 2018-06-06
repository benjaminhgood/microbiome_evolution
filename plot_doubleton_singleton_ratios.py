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
from math import log10,ceil,exp
from numpy.random import randint, normal

import core_gene_utils

intermediate_filename = '%ssubstitution_rates.txt' % (parse_midas_data.data_directory)
    

min_coverage = config.min_median_coverage
alpha = 0.5 # Confidence interval range for rate estimates
low_pi_threshold = 1e-03
low_divergence_threshold = 1e-03
min_change = 0.8
min_sample_size = 10
allowed_variant_types = set(['1D','2D','3D','4D'])


if __name__=='__main__':


    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--debug", help="Loads only a subset of SNPs for speed", action="store_true")
    parser.add_argument("--chunk-size", type=int, help="max number of records to load", default=1000000000)
    parser.add_argument("--species", help="Name of specific species to run code on", default="Bacteroides_vulgatus_57955")
    args = parser.parse_args()

    debug = args.debug
    chunk_size = args.chunk_size
    species_name=args.species

    # Load subject and sample metadata
    sys.stderr.write("Loading sample metadata...\n")
    subject_sample_map = parse_HMP_data.parse_subject_sample_map()
    sys.stderr.write("Done!\n")
    
    # Only plot samples above a certain depth threshold that are "haploids"
    snp_samples = diversity_utils.calculate_haploid_samples(species_name, debug=debug)
    # Only consider one sample per person
    snp_samples =     snp_samples[parse_midas_data.calculate_unique_samples(subject_sample_map, sample_list=snp_samples)]
    sys.stderr.write("Proceeding with %d haploid samples!\n" % len(snp_samples))
        
    if len(snp_samples) < min_sample_size:
        sys.stderr.write("Not enough haploid samples!\n")
        sys.exit(1)
                
    sys.stderr.write("Proceeding with %d haploid samples!\n" % len(snp_samples))

    sys.stderr.write("Loading core genes...\n")
    core_genes = core_gene_utils.parse_core_genes(species_name)
    non_shared_genes = core_gene_utils.parse_non_shared_reference_genes(species_name)
    shared_pangenome_genes = core_gene_utils.parse_shared_genes(species_name)
    sys.stderr.write("Done! Core genome consists of %d genes\n" % len(core_genes))
    sys.stderr.write("%d shared genes and %d non-shared genes\n" % (len(shared_pangenome_genes), len(non_shared_genes)))

        
    # Analyze SNPs, looping over chunk sizes. 
    # Clunky, but necessary to limit memory usage on cluster

    # Load SNP information for species_name
    sys.stderr.write("Loading SNPs for %s...\n" % species_name)
    sys.stderr.write("(core genes only...)\n")
    pi_matrix_syn = numpy.array([])
    avg_pi_matrix_syn = numpy.array([])

    syn_mut_difference_matrix = numpy.array([]) # 4d sites in core genes
    syn_mut_opportunity_matrix = numpy.array([])
    syn_rev_difference_matrix = numpy.array([]) # 4d sites in core genes
    syn_rev_opportunity_matrix = numpy.array([])
    
    non_mut_difference_matrix = numpy.array([]) # 1d sites in core genes
    non_mut_opportunity_matrix = numpy.array([])
    non_rev_difference_matrix = numpy.array([]) # 1d sites in core genes
    non_rev_opportunity_matrix = numpy.array([])
    
    core_mut_difference_matrix = numpy.array([]) # all sites in core genes
    core_mut_opportunity_matrix = numpy.array([])
    core_rev_difference_matrix = numpy.array([]) # all sites in core genes
    core_rev_opportunity_matrix = numpy.array([])
    
    core_doubleton_matrix = numpy.array([]) # all sites in core genes
    core_doubleton_opportunity_matrix = numpy.array([])
    
    final_line_number = 0
    while final_line_number >= 0:
    
        sys.stderr.write("Loading chunk starting @ %d...\n" % final_line_number)
        dummy_samples, allele_counts_map, passed_sites_map, final_line_number = parse_midas_data.parse_snps(species_name, debug=debug, allowed_samples=snp_samples, chunk_size=chunk_size,initial_line_number=final_line_number, allowed_genes=core_genes)
        sys.stderr.write("Done! Loaded %d genes\n" % len(allele_counts_map.keys()))
            
        # Calculate fixation matrix
        sys.stderr.write("Calculating matrix of snp differences...\n")
        # Synonymous (4D)
            
        chunk_syn_mut_difference_matrix, chunk_syn_rev_difference_matrix, chunk_syn_mut_opportunity_matrix, chunk_syn_rev_opportunity_matrix = diversity_utils.calculate_mutation_reversion_matrix(allele_counts_map, passed_sites_map, allowed_genes=core_genes, allowed_variant_types=set(['4D']))
            
        # Nonsynonymous (1D) 
        chunk_non_mut_difference_matrix, chunk_non_rev_difference_matrix, chunk_non_mut_opportunity_matrix, chunk_non_rev_opportunity_matrix = diversity_utils.calculate_mutation_reversion_matrix(allele_counts_map, passed_sites_map, allowed_genes=core_genes, allowed_variant_types=set(['1D']))
            
        # Core (all)
        chunk_core_mut_difference_matrix, chunk_core_rev_difference_matrix, chunk_core_mut_opportunity_matrix, chunk_core_rev_opportunity_matrix = diversity_utils.calculate_mutation_reversion_matrix(allele_counts_map, passed_sites_map, allowed_genes=core_genes)
        
        chunk_doubleton_matrix, chunk_doubleton_opportunity_matrix = diversity_utils.calculate_doubleton_matrix(allele_counts_map, passed_sites_map, allowed_genes=core_genes)   
            
        sys.stderr.write("Done!\n")
    
        if syn_mut_difference_matrix.shape[0]==0:
            syn_mut_difference_matrix = numpy.zeros_like(chunk_syn_mut_difference_matrix)*1.0
            syn_mut_opportunity_matrix = numpy.zeros_like(syn_mut_difference_matrix)*1.0
            syn_rev_difference_matrix = numpy.zeros_like(syn_mut_difference_matrix)*1.0
            syn_rev_opportunity_matrix = numpy.zeros_like(syn_mut_difference_matrix)*1.0
                
                
            non_mut_difference_matrix = numpy.zeros_like(syn_mut_difference_matrix)*1.0
            non_mut_opportunity_matrix = numpy.zeros_like(syn_mut_difference_matrix)*1.0
            non_rev_difference_matrix = numpy.zeros_like(syn_mut_difference_matrix)*1.0
            non_rev_opportunity_matrix = numpy.zeros_like(syn_mut_difference_matrix)*1.0
                
            core_mut_difference_matrix = numpy.zeros_like(syn_mut_difference_matrix)*1.0
            core_mut_opportunity_matrix = numpy.zeros_like(syn_mut_difference_matrix)*1.0
            core_rev_difference_matrix = numpy.zeros_like(syn_mut_difference_matrix)*1.0
            core_rev_opportunity_matrix = numpy.zeros_like(syn_mut_difference_matrix)*1.0
            
            core_doubleton_matrix = numpy.zeros_like(chunk_doubleton_matrix)
            core_doubleton_opportunity_matrix = numpy.zeros_like(core_doubleton_matrix)
                
        # Add syn
        syn_mut_difference_matrix += chunk_syn_mut_difference_matrix
        syn_mut_opportunity_matrix += chunk_syn_mut_opportunity_matrix
        syn_rev_difference_matrix += chunk_syn_rev_difference_matrix
        syn_rev_opportunity_matrix += chunk_syn_rev_opportunity_matrix

        # Add non
        non_mut_difference_matrix += chunk_non_mut_difference_matrix
        non_mut_opportunity_matrix += chunk_non_mut_opportunity_matrix
        non_rev_difference_matrix += chunk_non_rev_difference_matrix
        non_rev_opportunity_matrix += chunk_non_rev_opportunity_matrix

        # Add core
        core_mut_difference_matrix += chunk_core_mut_difference_matrix
        core_mut_opportunity_matrix += chunk_core_mut_opportunity_matrix
        core_rev_difference_matrix += chunk_core_rev_difference_matrix
        core_rev_opportunity_matrix += chunk_core_rev_opportunity_matrix

        # Add doubletons
        core_doubleton_matrix += chunk_doubleton_matrix
        core_doubleton_opportunity_matrix += chunk_doubleton_opportunity_matrix
            
        snp_samples = dummy_samples
    
    sys.stderr.write("Done!\n")
        
    doubleton_fractions = core_doubleton_matrix*1.0/(core_doubleton_opportunity_matrix+(core_doubleton_opportunity_matrix==0))
        
    difference_matrix = core_mut_difference_matrix+core_rev_difference_matrix
    opportunity_matrix = core_mut_opportunity_matrix + core_rev_opportunity_matrix
        
    substitution_rates = difference_matrix*1.0/(opportunity_matrix+(opportunity_matrix==0))

    pylab.figure()
    pylab.xlabel('Synonymous divergence, $d_S$')
    pylab.ylabel('Private SNV sharing')
    pylab.semilogx([1e-07,1e-01],[0,0],'k:')
    pylab.xlim([3e-07,1e-01])
    pylab.ylim([-0.05,1.05])
        
    # Add records to output
    # Calculate which pairs of idxs belong to the same sample, which to the same subject
    # and which to different subjects
    same_sample_idxs, same_subject_idxs, diff_subject_idxs = parse_midas_data.calculate_subject_pairs(subject_sample_map, snp_samples)

    for idxs in [diff_subject_idxs]:
        for sample_pair_idx in xrange(0,len(idxs[0])):
                
            # do both order
            for i,j in [(idxs[0][sample_pair_idx], idxs[1][sample_pair_idx]), (idxs[1][sample_pair_idx], idxs[0][sample_pair_idx])]:
                
                if opportunity_matrix[i,j]==0:
                    continue
                        
                if core_doubleton_opportunity_matrix[i,j]<0.5:
                    continue
                
                ds = max([substitution_rates[i,j],1e-06])
                df = doubleton_fractions[i,j]
                
                #print ds, df    
                
                ds = ds*exp(normal(0,0.2))
                    
                pylab.plot([ds],[df],'b.',alpha=0.5,markersize=3)


    pylab.savefig('doubleton_fraction.pdf',bbox_inches='tight')