import matplotlib  
matplotlib.use('Agg') 
import sample_utils
import config
import parse_midas_data
import os.path
import pylab
import sys
import numpy
import sfs_utils
        

import diversity_utils
import gene_diversity_utils
import core_gene_utils
import gzip
import os

temporal_change_directory = '%stemporal_changes/' % (parse_midas_data.data_directory)
intermediate_filename_template = '%s%s.txt.gz'  

min_coverage = config.min_median_coverage
min_sample_size = 2

import stats_utils
from math import log10,ceil
from numpy.random import randint

def load_temporal_change_map(species_name):
    
    intermediate_filename = intermediate_filename_template % (temporal_change_directory, species_name)

    temporal_change_map = {}


    if not os.path.isfile(intermediate_filename):
        return temporal_change_map
    
    file = gzip.open(intermediate_filename,"r")
    file.readline() # header
    for line in file:
        items = line.split(",")
        if items[0].strip()!=species_name:
            continue
            
        sample_1 = items[1].strip()
        sample_2 = items[2].strip()
        type = items[3].strip()
        num_opportunities = float(items[4])
        perr = float(items[5])
        sample_pair = (sample_1, sample_2)
        if sample_pair not in temporal_change_map:
            temporal_change_map[sample_pair] = {}
        
        changes = []
        if len(items)<7:
            pass
        else:
            change_strs = items[6:]
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
                    
                elif type=='private_snps':
                    
                    gene_name = subitems[0].strip()
                    contig = subitems[1].strip()
                    position = long(subitems[2])
                    variant_type = subitems[3].strip()
                    A1 = float(subitems[4])
                    D1 = float(subitems[5])
                    A2 = float(subitems[6])
                    D2 = float(subitems[7])
                    changes.append( (gene_name, contig, position, variant_type, A1, D1, A2, D2) )
                    
        temporal_change_map[sample_pair][type] = num_opportunities, perr, changes
    
    return temporal_change_map

def calculate_private_reversions_from_temporal_change_map(temporal_change_map, sample_1, sample_2, lower_threshold=config.consensus_lower_threshold, 
upper_threshold=config.consensus_upper_threshold):
    
    sample_pair = sample_1, sample_2
    if sample_pair not in temporal_change_map:
        return -1, None, None
        
    if 'private_snps' not in temporal_change_map[sample_pair]:
        return -1, None, None
        
    # otherwise, some hope! 
    
    private_snp_opportunities, private_snp_perr, private_snps = temporal_change_map[sample_pair]['private_snps']
    
    mutations = []
    private_snp_reversions = []
    for snp_change in private_snps:
    
        a,b,c,d,A1,D1,A2,D2 = snp_change
        
        if D1==0 or D2==0:
            private_snp_opportunities-=1
            continue
        
        f1 = A1*1.0/D1
        f2 = A2*1.0/D2
        
        if f1>=upper_threshold and f2<=lower_threshold:
            private_snp_reversions.append(snp_change)
        if f1<=upper_threshold and f2>=upper_threshold:
            mutations.append(snp_change)        
    
    return private_snp_opportunities, private_snp_perr, private_snp_reversions


def calculate_mutations_reversions_from_temporal_change_map(temporal_change_map, sample_1, sample_2, lower_threshold=config.consensus_lower_threshold, 
upper_threshold=config.consensus_upper_threshold):

    sample_pair = sample_1, sample_2
    if sample_pair not in temporal_change_map:
        return -1, -1, [], []
        
    if 'snps' not in temporal_change_map[sample_pair]:
        return -1, -1, [], []
        
    # otherwise, some hope! 
    snp_opportunities, snp_perr, snp_changes = temporal_change_map[sample_pair]['snps']
    
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
            
    
    return snp_opportunities, snp_perr, mutations, reversions


def calculate_gains_losses_from_temporal_change_map(temporal_change_map, sample_1, sample_2, max_absent_copynum=config.gainloss_max_absent_copynum, min_normal_copynum=config.gainloss_min_normal_copynum, max_normal_copynum=config.gainloss_max_normal_copynum):


    sample_pair = sample_1, sample_2
    if sample_pair not in temporal_change_map:
        return -1, -1, [], []
        
    if 'genes' not in temporal_change_map[sample_pair]:
        return -1, -1, [], []
        
    # otherwise, some hope! 
    gene_opportunities, gene_perr, gene_changes = temporal_change_map[sample_pair]['genes']
    
    gains = []
    losses = []
    for gene_change in gene_changes:
    
        gene_name, D1, Dm1, D2, Dm2 = gene_change
        
        copynum_1 = D1/Dm1
        copynum_2 = D2/Dm2
        
        if (copynum_1<=max_absent_copynum) and (copynum_2>=min_normal_copynum) and (copynum_2<=max_normal_copynum):
            gains.append(gene_change)
        elif (copynum_2<=max_absent_copynum) and (copynum_1>=min_normal_copynum) and (copynum_1<=max_normal_copynum):
            losses.append(gene_change)
            
    
    return gene_opportunities, gene_perr, gains, losses


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

    os.system('mkdir -p %s' % temporal_change_directory)
 
    # Load subject and sample metadata
    sys.stderr.write("Loading sample metadata...\n")
    subject_sample_map = sample_utils.parse_subject_sample_map()
    sample_order_map = sample_utils.parse_sample_order_map()
    sys.stderr.write("Done!\n")

    intermediate_filename = intermediate_filename_template % (temporal_change_directory, species_name)

    output_file = gzip.open(intermediate_filename,"w")
    # header!
    output_file.write(", ".join(['Species', 'Sample1', 'Sample2', 'Type', 'L','Perr', 'Change1', '...']))
    output_file.write("\n")
    
    for species_name in good_species_list:
        
        sample_coverage_map = parse_midas_data.parse_sample_coverage_map(species_name)

        sys.stderr.write("Loading SFSs for %s...\t" % species_name)
        samples, sfs_map = parse_midas_data.parse_within_sample_sfs(species_name, allowed_variant_types=set(['1D','2D','3D','4D'])) 
        sys.stderr.write("Done!\n")


        sys.stderr.write("Loading temporal samples...\n")
        # Only plot samples above a certain depth threshold that are involved in timecourse
        snp_samples = diversity_utils.calculate_temporal_samples(species_name)
    
        # On purpose looking at non-consecutive pairs too
        # (restriction to consecutive pairs is later)
        same_sample_idxs, same_subject_idxs, diff_subject_idxs = sample_utils.calculate_nonconsecutive_ordered_subject_pairs(sample_order_map, snp_samples)
        
        if len(same_subject_idxs[0]) < min_sample_size:
            sys.stderr.write("Not enough temporal samples!\n")
            continue
          
        
        sys.stderr.write("Proceeding with %d comparisons of %d temporal samples!\n" % (len(same_subject_idxs[0]), len(snp_samples)))

        # Analyze SNPs, looping over chunk sizes. 
        # Clunky, but necessary to limit memory usage on cluster

        sys.stderr.write("Loading whitelisted genes...\n")
        non_shared_genes = core_gene_utils.parse_non_shared_reference_genes(species_name)
        shared_pangenome_genes = core_gene_utils.parse_shared_genes(species_name)
        sys.stderr.write("Done! %d shared genes and %d non-shared genes\n" % (len(shared_pangenome_genes), len(non_shared_genes)))

        # Now calculate gene differences
        # Load gene coverage information for species_name
        sys.stderr.write("Loading pangenome data for %s...\n" % species_name)
        gene_samples, gene_names, gene_presence_matrix, gene_depth_matrix, marker_coverages, gene_reads_matrix = parse_midas_data.parse_pangenome_data(species_name,allowed_samples=snp_samples, disallowed_genes=shared_pangenome_genes)
        sys.stderr.write("Done!\n")
    
        snp_samples = gene_samples
    
        same_sample_idxs, same_subject_idxs, diff_subject_idxs = sample_utils.calculate_nonconsecutive_ordered_subject_pairs(sample_order_map, snp_samples)
    
        if len(same_subject_idxs[0]) < min_sample_size:
            sys.stderr.write("Not enough temporal samples!\n")
            continue
        
        # Use if HMP    
        import calculate_private_snvs
        private_snv_map = calculate_private_snvs.load_private_snv_map(species_name)
    
        # If other dataset, use this (and uncomment private snv line below)
        #import calculate_snp_prevalences
        #snv_freq_map = calculate_snp_prevalences.parse_population_freqs(species_name,polarize_by_consensus=True)
        
        
        # Load SNP information for species_name
        sys.stderr.write("Loading SNPs for %s...\n" % species_name)    
        snp_changes = {}
        gene_changes = {}
        tracked_private_snps = {}
        snp_opportunities = {}
        gene_opportunities = {}
        tracked_private_snp_opportunities = {}
        
        snp_perrs = {}
        gene_perrs = {}
        tracked_private_snp_perrs = {}
        
        snp_difference_matrix = numpy.array([]) # all sites in all genes
        snp_opportunity_matrix = numpy.array([])
    
        final_line_number = 0

        while final_line_number >= 0:
    
            sys.stderr.write("Loading chunk starting @ %d...\n" % final_line_number)
            dummy_samples, allele_counts_map, passed_sites_map, final_line_number = parse_midas_data.parse_snps(species_name, debug=debug, allowed_samples=snp_samples, chunk_size=chunk_size,initial_line_number=final_line_number,allowed_genes=non_shared_genes)
            sys.stderr.write("Done! Loaded %d genes\n" % len(allele_counts_map.keys()))
            
            # Calculate private snvs
            # (use this if applying to a different dataset)
            # private_snv_map = calculate_private_snvs.calculate_private_snv_map_from_external_snv_prevalence_map(allele_counts_map, snv_freq_map)
        
            # All
            chunk_snp_difference_matrix, chunk_snp_opportunity_matrix = diversity_utils.calculate_fixation_matrix(allele_counts_map, passed_sites_map)   # 
    
            if snp_difference_matrix.shape[0]==0:
                snp_difference_matrix = numpy.zeros_like(chunk_snp_difference_matrix)*1.0
                snp_opportunity_matrix = numpy.zeros_like(snp_difference_matrix)*1.0
                
            # Add all
            snp_difference_matrix += chunk_snp_difference_matrix
            snp_opportunity_matrix += chunk_snp_opportunity_matrix

        
            #same_sample_idxs, same_subject_idxs, diff_subject_idxs = sample_utils.calculate_ordered_subject_pairs(sample_order_map, snp_samples)
        
            for sample_pair_idx in xrange(0,len(same_subject_idxs[0])):
    
                i = same_subject_idxs[0][sample_pair_idx]
                j = same_subject_idxs[1][sample_pair_idx]

                sample_i = snp_samples[i]
                sample_j = snp_samples[j]
                
                avg_depth_i = sample_coverage_map[sample_i]
                avg_depth_j = sample_coverage_map[sample_j]
                
                chunk_tracked_private_snps = diversity_utils.calculate_tracked_private_snvs(i, j, allele_counts_map, passed_sites_map, avg_depth_i, avg_depth_j, private_snv_map)
                
                chunk_snp_changes = diversity_utils.calculate_snp_differences_between(i, j, allele_counts_map, passed_sites_map, avg_depth_i, avg_depth_j)
        
                sample_pair = (sample_i, sample_j)
        
                if sample_pair not in snp_changes:
                    snp_changes[sample_pair] = []
                    gene_changes[sample_pair] = []
                    snp_opportunities[sample_pair] = 0
                    gene_opportunities[sample_pair] = 0
                    snp_perrs[sample_pair] = -1
                    gene_perrs[sample_pair] = -1
                    
                    tracked_private_snps[sample_pair] = []
                    tracked_private_snp_opportunities[sample_pair] = 0
                    tracked_private_snp_perrs[sample_pair] = -1
                    
                snp_changes[sample_pair].extend(chunk_snp_changes)
                snp_opportunities[sample_pair] += chunk_snp_opportunity_matrix[i,j]
                
                tracked_private_snps[sample_pair].extend( chunk_tracked_private_snps)
                
                tracked_private_snp_opportunities[sample_pair] += len(chunk_tracked_private_snps)
                               
        # Calculate SNP error rate
        for sample_pair_idx in xrange(0,len(same_subject_idxs[0])):
    
            i = same_subject_idxs[0][sample_pair_idx]
            j = same_subject_idxs[1][sample_pair_idx]

            sample_i = snp_samples[i]
            sample_j = snp_samples[j]
            sample_pair = (sample_i, sample_j)
        
            perr = diversity_utils.calculate_fixation_error_rate(sfs_map, sample_i, sample_j)[0]
            
            snp_perrs[sample_pair] = perr
            tracked_private_snp_perrs[sample_pair] = perr
    
            gene_changes[sample_pair].extend( gene_diversity_utils.calculate_gene_differences_between(i, j, gene_reads_matrix, gene_depth_matrix, marker_coverages) )
            
            gene_perr = gene_diversity_utils.calculate_gene_error_rate(i, j, gene_reads_matrix, gene_depth_matrix, marker_coverages)[0]
            
            gene_opportunities[sample_pair] = gene_depth_matrix.shape[0]
            
            gene_perrs[sample_pair] = gene_perr
    
        sys.stderr.write("Done!\n") 
    
        for sample_i, sample_j in snp_changes.keys():
        
            # First output SNPs
            snp_strs = []
            for snp_change in snp_changes[(sample_i, sample_j)]:
        
            
                gene_name, location, variant_type, allele_counts_1, allele_counts_2 = snp_change
                contig = location[0]
                position = location[1]
            
                A1,D1 = allele_counts_1
                A2,D2 = allele_counts_2
            
                snp_str = ('%s;%s;%d;%s;%d;%d;%d;%d' % (gene_name, contig, position, variant_type, A1, D1, A2, D2))
            
                snp_strs.append(snp_str)
            
            record_str_items = [species_name, sample_i, sample_j, 'snps', "%g" % snp_opportunities[(sample_i, sample_j)], "%g" % snp_perrs[(sample_i, sample_j)]] + snp_strs
            record_str = ", ".join(record_str_items)
            output_file.write(record_str)
            output_file.write("\n")
            
            # Now output genes
            gene_strs = []
            for gene_change in gene_changes[(sample_i, sample_j)]:
                gene_idx, coverages_1, coverages_2 = gene_change
                gene_name = gene_names[gene_idx]
                D1,Dm1 = coverages_1
                D2,Dm2 = coverages_2
            
                gene_str = ('%s;%0.2f;%0.2f;%0.2f;%0.2f' % (gene_name, D1, Dm1, D2, Dm2))
                gene_strs.append(gene_str)
        
            record_str_items = [species_name, sample_i, sample_j, 'genes', "%g" % gene_opportunities[(sample_i, sample_j)], "%g" % gene_perrs[(sample_i, sample_j)]] + gene_strs
            record_str = ", ".join(record_str_items)
            output_file.write(record_str)
            output_file.write("\n")
    
            # Now output private SNPS
            private_snp_strs = []
            for snp_change in tracked_private_snps[(sample_i, sample_j)]:
        
            
                gene_name, location, variant_type, allele_counts_1, allele_counts_2 = snp_change
                contig = location[0]
                position = location[1]
            
                A1,D1 = allele_counts_1
                A2,D2 = allele_counts_2
            
                snp_str = ('%s;%s;%d;%s;%d;%d;%d;%d' % (gene_name, contig, position, variant_type, A1, D1, A2, D2))
            
                private_snp_strs.append(snp_str)
            
            record_str_items = [species_name, sample_i, sample_j, 'private_snps', "%g" % tracked_private_snp_opportunities[(sample_i, sample_j)], "%g" % tracked_private_snp_perrs[(sample_i, sample_j)]] + private_snp_strs
            record_str = ", ".join(record_str_items)
            output_file.write(record_str)
            output_file.write("\n")
            
    
        sys.stderr.write("Done with %s!\n" % species_name) 
    
    sys.stderr.write("Done looping over species!\n")
    output_file.close()
    sys.stderr.write("Done!\n")
    
    # testing loading of intermediate file
    temporal_change_map = load_temporal_change_map(good_species_list[0])
 
