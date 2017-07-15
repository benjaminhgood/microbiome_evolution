import matplotlib  
matplotlib.use('Agg') 
import config
import parse_midas_data
import parse_HMP_data
import os.path
import pylab
import sys
import numpy
import sfs_utils
        

import diversity_utils
import gene_diversity_utils

import stats_utils
from math import log10,ceil
from numpy.random import randint

min_coverage = config.min_median_coverage
alpha = 0.5 # Confidence interval range for rate estimates
low_pi_threshold = 1e-03
low_divergence_threshold = 1e-03
min_change = 0.8
min_sample_size = 5

intermediate_filename = '%stemporal_changes.txt' % (parse_midas_data.analysis_directory)

def load_temporal_change_map(species_name):
    
    temporal_change_map = {}

    file = open(intermediate_filename,"r")
    file.readline() # header
    for line in file:
        items = line.split(",")
        if items[0].strip()!=species_name:
            continue
            
        sample_1 = items[1].strip()
        sample_2 = items[2].strip()
        type = items[3].strip()
        perr = float(items[4])
        sample_pair = (sample_1, sample_2)
        if sample_pair not in temporal_change_map:
            temporal_change_map[sample_pair] = {}
        
        changes = []
        if len(items)<6:
            pass
        else:
            change_strs = items[5:]
            for change_str in change_strs:
            
                subitems = change_str.split(";")
                
                # switch on type of change
                if type=='snps':    
                    gene_name = subitems[0].strip()
                    contig = subitems[1].strip()
                    position = long(subitems[2])
                    variant_type = subitems[3].strip()
                    A1 = float(subitems[4])
                    D1 = float(subitems[5])
                    A2 = float(subitems[6])
                    D2 = float(subitems[7])
                    changes.append( (gene_name, contig, position, variant_type, A1, D1, A2, D2) )
                            
                elif type=='genes':
                    gene_name = subitems[0].strip()
                    D1 = float(subitems[1])
                    Dm1 = float(subitems[2])
                    D2 = float(subitems[3])
                    Dm2 = float(subitems[4])
                    changes.append( (gene_name, D1, Dm1, D2, Dm2) )
                
        temporal_change_map[sample_pair][type] = perr, changes
    
    return temporal_change_map

def calculate_mutations_reversions_from_temporal_change_map(temporal_change_map, sample_1, sample_2, lower_threshold=config.consensus_lower_threshold, 
upper_threshold=config.consensus_upper_threshold):

    sample_pair = sample_1, sample_2
    if sample_pair not in temporal_change_map:
        return -1, None, None
        
    if 'snps' not in temporal_change_map[sample_pair]:
        return -1, None, None
        
    # otherwise, some hope! 
    snp_perr, snp_changes = temporal_change_map[sample_pair]['snps']
    
    mutations = []
    reversions = []
    for snp_change in snp_changes:
    
        a,b,c,d,A1,D1,A2,D2 = snp_change
        
        f1 = A1*1.0/D1
        f2 = A2*1.0/D2
        
        if (f1<=lower_threshold) and (f2>=upper_threshold):
            mutations.append(snp_change)
        elif (f1>=upper_threshold) and (f2<=lower_threshold):
            reversions.append(snp_change)
            
    
    return snp_perr, mutations, reversions


def calculate_gains_losses_from_temporal_change_map(temporal_change_map, sample_1, sample_2, lower_threshold=0.05):

    sample_pair = sample_1, sample_2
    if sample_pair not in temporal_change_map:
        return -1, None, None
        
    if 'genes' not in temporal_change_map[sample_pair]:
        return -1, None, None
        
    # otherwise, some hope! 
    gene_perr, gene_changes = temporal_change_map[sample_pair]['genes']
    
    gains = []
    losses = []
    for gene_change in gene_changes:
    
        gene_name, D1, Dm1, D2, Dm2 = gene_change
        
        copynum_1 = D1/Dm1
        copynum_2 = D2/Dm2
        
        if (copynum_1<=lower_threshold):
            gains.append(gene_change)
        elif (copynum_2<=lower_threshold):
            losses.append(gene_change)
            
    
    return gene_perr, gains, losses


if __name__=='__main__':

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--debug", help="Loads only a subset of SNPs for speed", action="store_true")
    parser.add_argument("--chunk-size", type=int, help="max number of records to load", default=1000000000)

    args = parser.parse_args()

    debug = args.debug
    chunk_size = args.chunk_size
         
    # Load subject and sample metadata
    sys.stderr.write("Loading sample metadata...\n")
    subject_sample_map = parse_HMP_data.parse_subject_sample_map()
    sample_order_map = parse_HMP_data.parse_sample_order_map()
    sys.stderr.write("Done!\n")
    
    good_species_list = parse_midas_data.parse_good_species_list()
    if debug:
        good_species_list = good_species_list[:3]

    record_strs = [", ".join(['Species', 'Sample1', 'Sample2', 'Type', 'FalsePositives', 'Change1', '...'])]

    for species_name in good_species_list:

        sys.stderr.write("Loading SFSs for %s...\t" % species_name)
        samples, sfs_map = parse_midas_data.parse_within_sample_sfs(species_name, allowed_variant_types=set(['1D','2D','3D','4D'])) 
        sys.stderr.write("Done!\n")


        sys.stderr.write("Loading temporal samples...\n")
        # Only plot samples above a certain depth threshold that are involved in timecourse
        snp_samples = diversity_utils.calculate_temporal_samples(species_name)
    
        same_sample_idxs, same_subject_idxs, diff_subject_idxs = parse_midas_data.calculate_ordered_subject_pairs(sample_order_map, snp_samples)
        sample_size = len(same_subject_idxs[0])
    
        if sample_size < min_sample_size:
            sys.stderr.write("Not enough temporal samples!\n")
            continue
        
        sys.stderr.write("Proceeding with %d temporal samples!\n" % sample_size)

        # Analyze SNPs, looping over chunk sizes. 
        # Clunky, but necessary to limit memory usage on cluster

    
        # Load SNP information for species_name
        sys.stderr.write("Loading SNPs for %s...\n" % species_name)    
        snp_changes = {}
        gene_changes = {}
        
        snp_perrs = {}
        gene_perrs = {}
        
        snp_difference_matrix = numpy.array([]) # all sites in all genes
        snp_opportunity_matrix = numpy.array([])
    
        final_line_number = 0

        while final_line_number >= 0:
    
            sys.stderr.write("Loading chunk starting @ %d...\n" % final_line_number)
            dummy_samples, allele_counts_map, passed_sites_map, final_line_number = parse_midas_data.parse_snps(species_name, debug=debug, allowed_samples=snp_samples, chunk_size=chunk_size,initial_line_number=final_line_number)
            sys.stderr.write("Done! Loaded %d genes\n" % len(allele_counts_map.keys()))
            snp_samples = dummy_samples
        
            # All
            chunk_snp_difference_matrix, chunk_snp_opportunity_matrix = diversity_utils.calculate_fixation_matrix(allele_counts_map, passed_sites_map, min_change=min_change)   # 
    
            if snp_difference_matrix.shape[0]==0:
                snp_difference_matrix = numpy.zeros_like(chunk_snp_difference_matrix)*1.0
                snp_opportunity_matrix = numpy.zeros_like(snp_difference_matrix)*1.0
                
            # Add all
            snp_difference_matrix += chunk_snp_difference_matrix
            snp_opportunity_matrix += chunk_snp_opportunity_matrix

        
            same_sample_idxs, same_subject_idxs, diff_subject_idxs = parse_midas_data.calculate_ordered_subject_pairs(sample_order_map, snp_samples)
        
            for sample_pair_idx in xrange(0,len(same_subject_idxs[0])):
    
                i = same_subject_idxs[0][sample_pair_idx]
                j = same_subject_idxs[1][sample_pair_idx]

                sample_i = snp_samples[i]
                sample_j = snp_samples[j]
                chunk_snp_changes = diversity_utils.calculate_snp_differences_between(i, j, allele_counts_map, passed_sites_map)
        
                sample_pair = (sample_i, sample_j)
        
                if sample_pair not in snp_changes:
                    snp_changes[sample_pair] = []
                    gene_changes[sample_pair] = []
                    snp_perrs[sample_pair] = -1
                    gene_perrs[sample_pair] = -1
            
                snp_changes[sample_pair].extend(chunk_snp_changes)
    
        # Calculate SNP error rate
        same_sample_idxs, same_subject_idxs, diff_subject_idxs = parse_midas_data.calculate_ordered_subject_pairs(sample_order_map, snp_samples)   
        for sample_pair_idx in xrange(0,len(same_subject_idxs[0])):
    
            i = same_subject_idxs[0][sample_pair_idx]
            j = same_subject_idxs[1][sample_pair_idx]

            sample_i = snp_samples[i]
            sample_j = snp_samples[j]
            sample_pair = (sample_i, sample_j)
        
            perr = diversity_utils.calculate_fixation_error_rate(sfs_map, sample_i, sample_j)[0] * snp_opportunity_matrix[i, j]
            
            snp_perrs[sample_pair] = perr
    
        # Now calculate gene differences
        # Load gene coverage information for species_name
        sys.stderr.write("Loading pangenome data for %s...\n" % species_name)
        gene_samples, gene_names, gene_presence_matrix, gene_depth_matrix, marker_coverages, gene_reads_matrix = parse_midas_data.parse_pangenome_data(species_name,allowed_samples=snp_samples)
        sys.stderr.write("Done!\n")
    
        same_sample_idxs, same_subject_idxs, diff_subject_idxs = parse_midas_data.calculate_ordered_subject_pairs(sample_order_map, gene_samples)
    
        sys.stderr.write("Calculating gene changes...\n")
        for sample_pair_idx in xrange(0,len(same_subject_idxs[0])):
    
            i = same_subject_idxs[0][sample_pair_idx]
            j = same_subject_idxs[1][sample_pair_idx]
        
            if not (marker_coverages[i]>=min_coverage)*(marker_coverages[j]>=min_coverage):
                continue
            
            sample_i = gene_samples[i]
            sample_j = gene_samples[j]
        
            sample_pair = (sample_i, sample_j)
        
            if sample_pair not in gene_changes:
                snp_changes[sample_pair] = []
                gene_changes[sample_pair] = []
    
            gene_changes[sample_pair].extend( gene_diversity_utils.calculate_gene_differences_between(i, j, gene_reads_matrix, gene_depth_matrix, marker_coverages) )
            
            gene_perr = gene_diversity_utils.calculate_gene_error_rate(i, j, gene_reads_matrix, gene_depth_matrix, marker_coverages)[0]
            
            gene_perrs[sample_pair] = gene_perr
    
        sys.stderr.write("Done!\n") 
    
        for sample_i, sample_j in snp_changes.keys():
        
            snp_strs = []
        
            for snp_change in snp_changes[(sample_i, sample_j)]:
        
            
                gene_name, location, variant_type, allele_counts_1, allele_counts_2 = snp_change
                contig = location[0]
                position = location[1]
            
                A1,D1 = allele_counts_1
                A2,D2 = allele_counts_2
            
                snp_str = ('%s;%s;%d;%s;%d;%d;%d;%d' % (gene_name, contig, position, variant_type, A1, D1, A2, D2))
            
                snp_strs.append(snp_str)
            
            gene_strs = []
            for gene_change in gene_changes[(sample_i, sample_j)]:
                gene_idx, coverages_1, coverages_2 = gene_change
                gene_name = gene_names[gene_idx]
                D1,Dm1 = coverages_1
                D2,Dm2 = coverages_2
            
                gene_str = ('%s;%0.2f;%0.2f;%0.2f;%0.2f' % (gene_name, D1, Dm1, D2, Dm2))
                gene_strs.append(gene_str)
            
        
        
            record_str_items = [species_name, sample_i, sample_j, 'snps', "%g" % snp_perrs[(sample_i, sample_j)]] + snp_strs
            record_strs.append(", ".join(record_str_items))
        
            record_str_items = [species_name, sample_i, sample_j, 'genes', "%g" % gene_perrs[(sample_i, sample_j)]] + gene_strs
            record_strs.append(", ".join(record_str_items))
    
        sys.stderr.write("Done with %s!\n" % species_name) 
    
    sys.stderr.write("Done looping over species!\n")
    
    sys.stderr.write("Writing intermediate file...\n")
    file = open(intermediate_filename,"w")
    record_str = "\n".join(record_strs)
    file.write(record_str)
    file.close()
    sys.stderr.write("Done!\n")

 