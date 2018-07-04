import sample_utils
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
        

low_divergence_threshold = 5e-04 # this was picked by looking at inflection point of dN/dS vs dS plot 
min_sample_size = 10
allowed_variant_types = set(['1D','4D'])

def load_ld_map(species_name):

    ld_map = {}

    file = open(intermediate_filename,"r")
    header_line = file.readline() # header
    header_items = header_line.split(",")
    
    distance_strs = [item.split(":")[-1] for item in header_items[4:]]
    
    distances = []
    intragene_idxs = []
    
    intergene_distances = []
    intergene_idxs = []
    
    control_idx = -1
    
    for i in xrange(0,len(distance_strs)-1):
        
        if distance_strs[i].startswith('g'):
            # an intergene distance
            intergene_idxs.append(i)
            intergene_distances.append(long(distance_strs[i][1:]))
        else:
            # an intragene distance
            intragene_idxs.append(i)
            distances.append(float(distance_strs[i]))
            
    distances = numpy.array(distances)
    intragene_idxs = numpy.array(intragene_idxs)
    
    intergene_distances = numpy.array(intergene_distances)
    intergene_idxs = numpy.array(intergene_idxs)
    
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
        
        control_numerator = rsquared_numerators[control_idx]
        control_denominator = rsquared_denominators[control_idx]
        control_count = counts[control_idx]
        
        control_ld = control_numerator/control_denominator
        
        intragene_rsquared_numerators = rsquared_numerators[intragene_idxs]
        intragene_rsquared_denominators = rsquared_denominators[intragene_idxs]
        intragene_counts = counts[intragene_idxs]
        
        intergene_rsquared_numerators = rsquared_numerators[intergene_idxs]
        intergene_rsquared_denominators = rsquared_denominators[intergene_idxs]
        intergene_counts = counts[intergene_idxs]
        
        ld_map[(clade_type, variant_type)] = (distances, intragene_rsquared_numerators, intragene_rsquared_denominators, intragene_counts, intergene_distances, intergene_rsquared_numerators, intergene_rsquared_denominators, intergene_counts, control_numerator, control_denominator, control_count, pi)
        
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
    subject_sample_map = sample_utils.parse_subject_sample_map()
    sys.stderr.write("Done!\n")

    # load the identity of isolates and mixtures so that I can filter them
    isolates, mixtures=parse_HMP_data.list_of_isolates_and_mixtures()
    

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
    
    neighbor_distances = numpy.array([1,2,3,4,5])
                        
    
    distance_strs = ["LD_N:LD_D:%g" % d for d in distance_bin_locations[1:-1]] # N=numerator and D=denominator
    distance_strs = distance_strs+["LD_N:LD_D:g%d" % nd for nd in neighbor_distances]+["LD_N:LD_D:intergene"]
    
    # header of the output file.
    record_strs = [", ".join(['Species', 'CladeType', 'VariantType', 'Pi']+distance_strs)]
    
    for species_name in good_species_list:

        sys.stderr.write("Loading haploid samples...\n")
        # Only plot samples above a certain depth threshold that are "haploids"
        snp_samples = diversity_utils.calculate_haploid_samples(species_name, debug=debug)
    
        # Only consider samples from isolates
        snp_samples_isolates=[]
        for sample in snp_samples:
            if sample in isolates:
                snp_samples_isolates.append(sample)

        snp_samples=numpy.asarray(snp_samples_isolates)

        if len(snp_samples) < min_sample_size:
            sys.stderr.write("Not enough haploid samples!\n")
            continue
        else:
            sys.stderr.write("Found %d haploid samples!\n" % len(snp_samples))
        
        sys.stderr.write("Calculating unique hosts...\n")
        # Only consider one sample per person
        snp_samples = snp_samples[sample_utils.calculate_unique_samples(subject_sample_map, sample_list=snp_samples)]


        if len(snp_samples) < min_sample_size:
            sys.stderr.write("Not enough hosts!\n")
            continue
        else:
            sys.stderr.write("Found %d unique hosts!\n" % len(snp_samples)) 

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
        # CHANGE THIS FOR SIMULATIONS:
        #clade_sets = clade_utils.load_manual_clades(species_name)
        clade_sets=[set(coarse_grained_samples)]

        sys.stderr.write("%d samples remaining after clustering!\n" % len(coarse_grained_samples))

        if len(clade_sets)==0:
            continue

        clade_idxss = clade_utils.calculate_clade_idxs_from_clade_sets(coarse_grained_samples, clade_sets)

        clade_sizes = numpy.array([clade_idxs.sum() for clade_idxs in clade_idxss])

        largest_clade_samples = coarse_grained_samples[ clade_idxss[clade_sizes.argmax()] ]

        largest_clade_set = set(largest_clade_samples)

        sys.stderr.write("Top level clades: %d clades, sizes: %s\n" % (len(clade_sets), str(clade_sizes)))
        sys.stderr.write("Max clade size: %d\n" % len(largest_clade_samples))

        snp_samples = coarse_grained_samples

        if len(largest_clade_samples) < min_sample_size:
            sys.stderr.write("Not enough haploid samples!\n")
            continue
        else:
            sys.stderr.write("Proceeding with %d coarse-grained samples in largest clade!\n" % len(largest_clade_samples))
                
        
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

        neighboring_gene_rsquared_numerators = {}
        neighboring_gene_rsquared_denominators = {}
        neighboring_gene_counts = {}

        # total_control=between genes.  
        total_control_rsquared_numerators = {}
        total_control_rsquared_denominators = {}
        total_control_counts = {}
        
        for clade_type in clade_types:
            for variant_type in variant_types:
        
                binned_rsquared_numerators[(clade_type,variant_type)] = numpy.zeros_like(distance_bin_locations)
                binned_rsquared_denominators[(clade_type,variant_type)] = numpy.zeros_like(distance_bin_locations)  
                binned_counts[(clade_type,variant_type)] = numpy.zeros_like(distance_bin_locations)  
         
                neighboring_gene_rsquared_numerators[(clade_type,variant_type)] = numpy.zeros_like(neighbor_distances)*1.0
                neighboring_gene_rsquared_denominators[ (clade_type,variant_type)] = numpy.zeros_like(neighbor_distances)*1.0
                neighboring_gene_counts[(clade_type,variant_type)] = numpy.zeros_like(neighbor_distances)*1.0  
         
                total_control_rsquared_numerators[(clade_type,variant_type)] = 0
                total_control_rsquared_denominators[(clade_type,variant_type)] = 0
                total_control_counts[(clade_type,variant_type)] = 0
                
                
            
        final_line_number = 0
        while final_line_number >= 0:
    
            sys.stderr.write("Loading chunk starting @ %d...\n" % final_line_number)
            snp_samples, allele_counts_map, passed_sites_map, final_line_number = parse_midas_data.parse_snps(species_name, debug=debug,         allowed_variant_types=allowed_variant_types, allowed_samples=snp_samples,allowed_genes=core_genes, chunk_size=chunk_size,initial_line_number=final_line_number)
            sys.stderr.write("Done! Loaded %d genes\n" % len(allele_counts_map.keys()))
    
            largest_clade_idxs = numpy.array([sample in largest_clade_set for sample in snp_samples])
    
            sys.stderr.write("Calculating LD...\n")
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
        
                        target_chromosome = allele_counts_map[gene_name][variant_type]['locations'][0][0]
        
                        if clade_type=='largest_clade':        
                            # Now restrict to largest clade
                            allele_counts = allele_counts[:,largest_clade_idxs,:]
            

                        #compute the distances between all pairs of sites 
                        # None in the two index positions results in a transpose of the vector relative to each other
                        # Subtraction between the two vectors results in pairwise subtraction of each element in each vector.
                        distances = numpy.fabs(locations[:,None]-locations[None,:])
    
                        rsquared_numerators, rsquared_denominators = diversity_utils.calculate_unbiased_sigmasquared(allele_counts, allele_counts)
                    
                        neighbor_rsquared_numeratorss = [[] for d in neighbor_distances]
                        neighbor_rsquared_denominatorss = [[] for d in neighbor_distances]
                        
                        for neighbor_distance_idx in xrange(0,len(neighbor_distances)):
                            
                            neighbor_distance = neighbor_distances[neighbor_distance_idx]     
                        
                            gene_name_items = gene_name.split(".")
                            gene_peg_number = long(gene_name.split(".")[-1])
                            nearest_gene_peg_numbers = [gene_peg_number-neighbor_distance,gene_peg_number+neighbor_distance]
                            neighboring_genes = [".".join(gene_name_items[:-1]+[str(n)]) for n in nearest_gene_peg_numbers]
                        
                            for neighboring_gene_name in neighboring_genes:
                    
                                # first make sure it's a real gene
                                if neighboring_gene_name not in allele_counts_map:
                                    continue
                                
                                if neighboring_gene_name not in core_genes:
                                    continue
        
                                neighboring_allele_counts = allele_counts_map[neighboring_gene_name][variant_type]['alleles']
                            
                                # then make sure it has some variants
                                if len(neighboring_allele_counts)==0:
                                    continue
                                
                                neighboring_chromosome = allele_counts_map[neighboring_gene_name][variant_type]['locations'][0][0]
                                
                                if neighboring_chromosome!=target_chromosome:
                                    continue
                                
                                if clade_type=='largest_clade':        
                                # Now restrict to largest clade
                                    neighboring_allele_counts = neighboring_allele_counts[:,largest_clade_idxs,:]
                                
                            
                                chunk_rsquared_numerators, chunk_rsquared_denominators = diversity_utils.calculate_unbiased_sigmasquared(allele_counts, neighboring_allele_counts)
                        
                                neighbor_rsquared_numeratorss[ neighbor_distance_idx].extend( chunk_rsquared_numerators.flatten() )
                                neighbor_rsquared_denominatorss[ neighbor_distance_idx].extend( chunk_rsquared_denominators.flatten() )
                            
                            neighbor_rsquared_numeratorss[ neighbor_distance_idx] = numpy.array( neighbor_rsquared_numeratorss[neighbor_distance_idx] )
                                
                            neighbor_rsquared_denominatorss[ neighbor_distance_idx] = numpy.array( neighbor_rsquared_denominatorss[neighbor_distance_idx] )
                               
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
                                
                                if control_gene_name not in core_genes:
                                    continue
        
                            
                                control_gene_peg_number = long(control_gene_name.split(".")[-1])
                            
                                control_allele_counts = allele_counts_map[control_gene_name][variant_type]['alleles']
                            
                                if len(control_allele_counts)==0:
                                    continue
                                
                                # make sure you don't have one too close by!     
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
        
                        for i in xrange(0,len(neighbor_distances)):
                            good_idxs = (neighbor_rsquared_denominatorss[i]>1e-09)
                            neighboring_gene_counts[(clade_type,variant_type)][i] += good_idxs.sum() 
                            
                            neighboring_gene_rsquared_numerators[ (clade_type,variant_type)][i] += neighbor_rsquared_numeratorss[i][good_idxs].sum()  
                            
                            neighboring_gene_rsquared_denominators[ (clade_type,variant_type)][i] += neighbor_rsquared_denominatorss[i][good_idxs].sum() 
                            
                                    
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
            
                gene_rsquared_strs = ["%g:%g:%d" % (rsquared_numerator, rsquared_denominator, count) for rsquared_numerator, rsquared_denominator, count in zip(neighboring_gene_rsquared_numerators[(clade_type,variant_type)], neighboring_gene_rsquared_denominators[(clade_type,variant_type)], neighboring_gene_counts[(clade_type,variant_type)])]
            
                control_rsquared_str = "%g:%g:%d" % (total_control_rsquared_numerators[(clade_type,variant_type)], total_control_rsquared_denominators[(clade_type,variant_type)], total_control_counts[(clade_type,variant_type)])
            
                pi_str = str(pi)
            
                record_str = ", ".join([species_name, clade_type, variant_type, pi_str]+rsquared_strs+gene_rsquared_strs+[control_rsquared_str])
            
                record_strs.append(record_str)
            
        sys.stderr.write("Done with %s!\n" % species_name) 
    
    sys.stderr.write("Done looping over species!\n")
    
    sys.stderr.write("Writing intermediate file...\n")
    file = open(intermediate_filename,"w")
    record_str = "\n".join(record_strs)
    file.write(record_str)
    file.close()
    sys.stderr.write("Done!\n")

 
