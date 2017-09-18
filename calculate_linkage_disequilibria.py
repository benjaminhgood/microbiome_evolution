import config
import parse_midas_data
import parse_HMP_data
import os.path
import pylab
import sys
import numpy

import diversity_utils
import gene_diversity_utils
import calculate_substitution_rates
import clade_utils
import stats_utils
from math import log10,ceil,fabs
from numpy.random import randint, choice

intermediate_filename = '%slinkage_disequilibria.txt' % (parse_midas_data.data_directory)
old_intermediate_filename = '%slinkage_disequilibria.txt.old' % (parse_midas_data.data_directory)
        

min_coverage = config.min_median_coverage
alpha = 0.5 # Confidence interval range for rate estimates
low_pi_threshold = 1e-03
low_divergence_threshold = 1e-03 #NRG: how was this picked?
min_change = 0.8
min_sample_size = 10
allowed_variant_types = set(['1D','4D'])

def load_ld_map(species_name):

    ld_map = {}

    file = open(intermediate_filename,"r")
    header_line = file.readline() # header
    header_items = header_line.split(",")
    
    distance_items = header_items[4:-1]
    distances = numpy.array([float(item.split(":")[-1]) for item in distance_items])
    
    for line in file:
        items = line.split(",")
        if items[0].strip()!=species_name:
            continue
        
        clade_type = items[1].strip()
        variant_type = items[2].strip()
        pi = float(items[3])
        
        rsquared_numerators = []
        rsquared_denominators = []
        lds = []
        counts = []
        for item in items[4:]:
            subitems = item.split(":")
            rsquared_numerators.append(float(subitems[0]))
            rsquared_denominators.append(float(subitems[1]))
            counts.append(float(subitems[2]))
        
        rsquared_numerators = numpy.array(rsquared_numerators)
        rsquared_denominators = numpy.array(rsquared_denominators)
        counts = numpy.array(counts)
        
        
        lds = rsquared_numerators/rsquared_denominators
        
        control_numerator = rsquared_numerators[-1]
        control_denominator = rsquared_denominators[-1]
        control_count = counts[-1]
        
        control_ld = control_numerator/control_denominator
        
        rsquared_numerators = rsquared_numerators[:-1]
        rsquared_denominators = rsquared_denominators[:-1]
        counts = counts[:-1]
        
        ld_map[(clade_type, variant_type)] = (distances, rsquared_numerators, rsquared_denominators, counts, control_numerator, control_denominator, control_count, pi)
        
    return ld_map


def load_ld_map_old(species_name):

    ld_map = {}

    file = open(old_intermediate_filename,"r")
    header_line = file.readline() # header
    header_items = header_line.split(",")
    
    distance_items = header_items[3:-1]
    distances = numpy.array([float(item.split(":")[1]) for item in distance_items])
    
    for line in file:
        items = line.split(",")
        if items[0].strip()!=species_name:
            continue
        
        variant_type = items[1].strip()
        pi = float(items[2])
        
        lds = []
        counts = []
        for item in items[3:]:
            subitems = item.split(":")
            lds.append(float(subitems[0]))
            counts.append(float(subitems[1]))
            
        lds = numpy.array(lds)
        counts = numpy.array(counts)
        
        ld_map[variant_type] = (distances, lds[:-1], counts[:-1], lds[-1], counts[-1], pi)
        
    return ld_map




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
    
    good_species_list = parse_midas_data.parse_good_species_list()
    if debug:
        good_species_list = good_species_list[:3]
    elif species !='all':
        good_species_list = [species]

    #good_species_list=['Bacteroides_vulgatus_57955'] 
    # better binning scheme (multiple of 3)
    distance_bin_locations = numpy.arange(0,1002)*3.0
    distance_bins = numpy.arange(-1,1002)*3+1.5
    distance_bins[0] = 0 # no such thing as negative distance
    distance_bins[1] = 2.5 # want at least one codon separation
    distance_bins[-1] = 1e09 # catch everything
    
    distance_strs = ["LD_N:LD_D:%g" % d for d in distance_bin_locations[1:-1]] # N=numerator and D=denominator
    distance_strs = distance_strs+["LD_N:LD_D:intergene"]
    

    # header of the output file.
    record_strs = [", ".join(['Species', 'CladeType', 'VariantType', 'Pi']+distance_strs)]
    
    for species_name in good_species_list:

        sys.stderr.write("Loading haploid samples...\n")
        # Only plot samples above a certain depth threshold that are "haploids"
        snp_samples = diversity_utils.calculate_haploid_samples(species_name, debug=debug)
    
        if len(snp_samples) < min_sample_size:
            sys.stderr.write("Not enough haploid samples!\n")
            continue
        
        sys.stderr.write("Calculating unique samples...\n")
        # Only consider one sample per person
        snp_samples = snp_samples[parse_midas_data.calculate_unique_samples(subject_sample_map, sample_list=snp_samples)]


        if len(snp_samples) < min_sample_size:
            sys.stderr.write("Not enough unique samples!\n")
            continue

        # Load divergence matrices 
        sys.stderr.write("Loading pre-computed substitution rates for %s...\n" % species_name)
        substitution_rate_map = calculate_substitution_rates.load_substitution_rate_map(species_name)
        sys.stderr.write("Calculating matrix...\n")
        dummy_samples, snp_difference_matrix, snp_opportunity_matrix = calculate_substitution_rates.calculate_matrices_from_substitution_rate_map(substitution_rate_map, 'core', allowed_samples=snp_samples)
        snp_samples = numpy.array(dummy_samples)
        substitution_rate = snp_difference_matrix*1.0/(snp_opportunity_matrix+(snp_opportunity_matrix==0))
        sys.stderr.write("Done!\n")

        sys.stderr.write("Clustering samples with low divergence...\n")
        coarse_grained_idxs, coarse_grained_cluster_list = clade_utils.cluster_samples(substitution_rate, min_d=low_divergence_threshold) # NRG: what is this returning?

        coarse_grained_samples = snp_samples[coarse_grained_idxs]
        clade_sets = clade_utils.load_manual_clades(species_name)

        if len(clade_sets)==0:
            continue

        clade_idxss = clade_utils.calculate_clade_idxs_from_clade_sets(coarse_grained_samples, clade_sets)

        clade_sizes = numpy.array([clade_idxs.sum() for clade_idxs in clade_idxss])

        largest_clade_samples = coarse_grained_samples[ clade_idxss[clade_sizes.argmax()] ]

        largest_clade_set = set(largest_clade_samples)

        sys.stderr.write("Top level: %d clades, %s\n" % (len(clade_sets), str(clade_sizes)))
        sys.stderr.write("Max: %d\n" % len(largest_clade_samples))

        snp_samples = coarse_grained_samples

        if len(largest_clade_samples) < min_sample_size:
            sys.stderr.write("Not enough haploid samples!\n")
            continue
                
        
        # Analyze SNPs, looping over chunk sizes. 
        # Clunky, but necessary to limit memory usage on cluster

        sys.stderr.write("Loading core genes...\n")
        core_genes = parse_midas_data.load_core_genes(species_name)
        sys.stderr.write("Done! Core genome consists of %d genes\n" % len(core_genes))

        # Load SNP information for species_name
        sys.stderr.write("Loading SNPs for %s...\n" % species_name)
        sys.stderr.write("(core genes only...)\n")
        
        clade_types = ['all','largest_clade']
        variant_types = ['4D','1D']
        
        binned_rsquared_numerators = {}
        binned_rsquared_denominators = {}
        binned_counts = {}

        # NRG: total_control=between genes? 
        total_control_rsquared_numerators = {}
        total_control_rsquared_denominators = {}
        total_control_counts = {}
        
        for clade_type in clade_types:
            for variant_type in variant_types:
        
                binned_rsquared_numerators[(clade_type,variant_type)] = numpy.zeros_like(distance_bin_locations)
                binned_rsquared_denominators[(clade_type,variant_type)] = numpy.zeros_like(distance_bin_locations)  
                binned_counts[(clade_type,variant_type)] = numpy.zeros_like(distance_bin_locations)  
         
                total_control_rsquared_numerators[(clade_type,variant_type)] = 0
                total_control_rsquared_denominators[(clade_type,variant_type)] = 0
                total_control_counts[(clade_type,variant_type)] = 0
            
        final_line_number = 0
        while final_line_number >= 0:
    
            sys.stderr.write("Loading chunk starting @ %d...\n" % final_line_number)
            snp_samples, allele_counts_map, passed_sites_map, final_line_number = parse_midas_data.parse_snps(species_name, debug=debug,         allowed_variant_types=allowed_variant_types, allowed_samples=snp_samples,allowed_genes=core_genes, chunk_size=chunk_size,initial_line_number=final_line_number)
            sys.stderr.write("Done! Loaded %d genes\n" % len(allele_counts_map.keys()))
    
            largest_clade_idxs = numpy.array([sample in largest_clade_set for sample in snp_samples])
    
            sys.stderr.write("Calculating intra-gene synonymous LD...\n")
            for clade_type in clade_types:
    
                for variant_type in variant_types:
            
                    for gene_name in allele_counts_map.keys():
            
                        if gene_name not in core_genes:
                            continue
        
                        locations = numpy.array([location for chromosome, location in allele_counts_map[gene_name][variant_type]['locations']])*1.0
                        allele_counts = allele_counts_map[gene_name][variant_type]['alleles']
        
        
                        
        
                        if len(allele_counts)==0:
                            # no diversity to look at!
                            continue
        
                        if clade_type=='largest_clade':        
                            # Now restrict to largest clade
                            allele_counts = allele_counts[:,largest_clade_idxs,:]
            

                        #compute the distances between all pairs of sites 
                        # None in the two index positions results in a transpose of the vector relative to each other
                        # Subtraction between the two vectors results in pairwise subtraction of each element in each vector.
                        distances = numpy.fabs(locations[:,None]-locations[None,:])
    
                        rsquared_numerators, rsquared_denominators = diversity_utils.calculate_unbiased_sigmasquared(allele_counts, allele_counts)
                    
                    
                        # pick a random gene somewhere else as a control
                        # 10 to 1 control to regular
                        control_rsquared_numerators = []
                        control_rsquared_denominators = []
                        gene_peg_number = long(gene_name.split(".")[-1])
                        for control_idx in xrange(0,10):
                    
                            control_gene_name = gene_name
                            control_allele_counts = []
                    
                            # get the right gene name
                            while True:
                        
                                control_gene_name = choice(allele_counts_map.keys())
                            
                                control_gene_peg_number = long(control_gene_name.split(".")[-1])
                            
                                control_allele_counts = allele_counts_map[control_gene_name][variant_type]['alleles']
                            
                                if len(control_allele_counts)==0:
                                    continue
                                    
                                if (fabs(control_gene_peg_number - gene_peg_number) < 5.5):
                                    continue
                            
                                if clade_type=='largest_clade':        
                                # Now restrict to largest clade
                                    control_allele_counts = control_allele_counts[:,largest_clade_idxs,:]
                                
                                break
            
                            
                                
                             
        
                            control_gene_rsquared_numerators, control_gene_rsquared_denominators = diversity_utils.calculate_unbiased_sigmasquared(allele_counts, control_allele_counts)
                        
                            control_rsquared_numerators.extend( control_gene_rsquared_numerators.flatten() )
                            control_rsquared_denominators.extend( control_gene_rsquared_denominators.flatten() )
                    
                        control_rsquared_numerators = numpy.array( control_rsquared_numerators )
                        control_rsquared_denominators = numpy.array( control_rsquared_denominators )  
        
                        # get the indices of the upper diagonal of the distance matrix
                        # numpy triu_indices returns upper diagnonal including diagonal
                        # the 1 inside the function excludes diagonal. Diagnonal has distance of zero.
                        desired_idxs = numpy.triu_indices(distances.shape[0],1)
        
                        #print distances.shape, rsquared_numerators.shape
        
                        # fetch the distances and rsquared vals corresponding to the upper diagonal. 
                        distances = distances[desired_idxs]
                        rsquared_numerators = rsquared_numerators[desired_idxs]
                        rsquared_denominators = rsquared_denominators[desired_idxs]
        
                        # fetch entries where denominator != 0 (remember, denominator=pa*(1-pa)*pb*(1-pb). If zero, then at least one site is invariant)
                        distances = distances[rsquared_denominators>1e-09]
                        rsquared_numerators = rsquared_numerators[rsquared_denominators>1e-09] 
                        rsquared_denominators = rsquared_denominators[rsquared_denominators>1e-09]
        
                        if len(distances) == 0:
                            continue

                        # numpy.digitize: For each distance value, return the bin index it belongs to in distances_bins. 
                        bin_idxs = numpy.digitize(distances,bins=distance_bins)-1
            
                        for i in xrange(0,len(bin_idxs)):
                    
                            binned_counts[(clade_type,variant_type)][bin_idxs[i]] += 1
                            binned_rsquared_numerators[(clade_type,variant_type)][bin_idxs[i]] += rsquared_numerators[i]
                            binned_rsquared_denominators[(clade_type,variant_type)][bin_idxs[i]] += rsquared_denominators[i]
        
                        total_control_counts[(clade_type,variant_type)] += (control_rsquared_denominators>1e-09).sum()
                        total_control_rsquared_numerators[(clade_type,variant_type)] += control_rsquared_numerators[control_rsquared_denominators>1e-09].sum()
                        total_control_rsquared_denominators[(clade_type,variant_type)] += control_rsquared_denominators[control_rsquared_denominators>1e-09].sum()
    
        for clade_type in clade_types:
            for variant_type in variant_types:
            
                desired_samples = snp_samples
                if clade_type=='largest_clade':
                    desired_samples = largest_clade_samples
            
                # Calculate pi!
                dummy_samples, snp_difference_matrix, snp_opportunity_matrix = calculate_substitution_rates.calculate_matrices_from_substitution_rate_map(substitution_rate_map, variant_type, allowed_samples=desired_samples)
                
                substitution_rate = snp_difference_matrix*1.0/(snp_opportunity_matrix+(snp_opportunity_matrix==0))
        
                iu = numpy.triu_indices(substitution_rate.shape[0], 1)
        
                pi = numpy.median(substitution_rate[iu])
        
                binned_rsquareds = binned_rsquared_numerators[(clade_type,variant_type)]*1.0/(binned_rsquared_denominators[(clade_type,variant_type)] + (binned_rsquared_denominators[(clade_type,variant_type)] == 0))
            
                control_rsquareds = total_control_rsquared_numerators[(clade_type,variant_type)]*1.0/(total_control_rsquared_denominators[(clade_type,variant_type)]+(total_control_rsquared_denominators[(clade_type,variant_type)]==0))   
            
                rsquared_strs = ["%g:%g:%d" % (rsquared_numerator, rsquared_denominator, count) for rsquared_numerator, rsquared_denominator, count in zip(binned_rsquared_numerators[(clade_type,variant_type)], binned_rsquared_denominators[(clade_type,variant_type)], binned_counts[(clade_type,variant_type)])[1:-1]]
            
                control_rsquared_str = "%g:%g:%d" % (total_control_rsquared_numerators[(clade_type,variant_type)], total_control_rsquared_denominators[(clade_type,variant_type)], total_control_counts[(clade_type,variant_type)])
            
                pi_str = str(pi)
            
                record_str = ", ".join([species_name, clade_type, variant_type, pi_str]+rsquared_strs+[control_rsquared_str])
            
                record_strs.append(record_str)
            
        sys.stderr.write("Done with %s!\n" % species_name) 
    
    sys.stderr.write("Done looping over species!\n")
    
    sys.stderr.write("Writing intermediate file...\n")
    file = open(intermediate_filename,"w")
    record_str = "\n".join(record_strs)
    file.write(record_str)
    file.close()
    sys.stderr.write("Done!\n")

 
