import matplotlib  
matplotlib.use('Agg') 
import sample_utils
import config
import parse_midas_data
import os.path
import pylab
import sys
import numpy
import diversity_utils
import gene_diversity_utils
import calculate_temporal_changes
import calculate_substitution_rates
import stats_utils
import sfs_utils
    
    
import matplotlib.colors as colors
import matplotlib.cm as cmx
from math import log10,ceil,log,exp
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy.random import randint, random, choice
import matplotlib.colors as mcolors
import matplotlib.patheffects as pe

from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster.hierarchy import cophenet
from scipy.cluster.hierarchy import fcluster


mpl.rcParams['font.size'] = 7
mpl.rcParams['lines.linewidth'] = 0.5
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'

species_name = "Bacteroides_vulgatus_57955"

################################################################################
#
# Standard header to read in argument information
#
################################################################################
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--debug", help="Loads only a subset of SNPs for speed", action="store_true")
parser.add_argument("--memoize", help="Loads stuff from disk", action="store_true")
parser.add_argument("--chunk-size", type=int, help="max number of records to load", default=1000000000)
parser.add_argument("--modification-threshold", type=int, help="max number of SNV differences before calling a modification", default=config.modification_difference_threshold)


args = parser.parse_args()

debug = args.debug
chunk_size = args.chunk_size
memoize = args.memoize
modification_difference_threshold = args.modification_threshold
replacement_difference_threshold = config.replacement_difference_threshold
################################################################################

min_coverage = config.min_median_coverage
min_sample_size = 3
min_haploid_sample_size = 10

variant_types = ['1D','4D']


# Must compete divergence matrix on the fly! 
            
# Load subject and sample metadata
sys.stderr.write("Loading sample metadata...\n")
sample_utils.parse_subject_sample_map
subject_sample_map = sample_utils.parse_subject_sample_map()
sample_order_map = sample_utils.parse_sample_order_map()
sample_country_map = sample_utils.parse_sample_country_map()
sys.stderr.write("Done!\n")
    
good_species_list = parse_midas_data.parse_good_species_list()
if debug:
    good_species_list = ["Bacteroides_vulgatus_57955", "Bacteroides_uniformis_57318"]

num_passed_species = 0

num_temporal_change_map = {}

total_snp_modification_map = {}
total_null_snp_modification_map = {}

total_gene_modification_map = {}
total_null_gene_modification_map = {}

# Variant type distribution of SNPs
total_snps = {var_type:0 for var_type in variant_types} # observed variant_type distribution
total_random_null_snps = {var_type:0 for var_type in variant_types} # expectation from randomly drawing sites on genome (mutation)
total_between_null_snps = {var_type:0 for var_type in variant_types} # expectation from randomly drawing sites that vary between samples (recombination) 

total_reversion_snps = {var_type: 0 for var_type in variant_types}
total_mutation_snps = {var_type: 0 for var_type in variant_types}
total_private_snps = {var_type: 0 for var_type in variant_types}

#derived_freq_bins = numpy.array([-1,0,1.0/32,1.0/16,1.0/8,1.0/4,1.0/2,3.0/4,7.0/8,15.0/16, 31.0/32,1,2])
#derived_virtual_xticks = list(derived_virtual_freqs[:-1]+0.5)
#derived_virtual_xticklabels = ['0','1/32','1/16','1/8','1/4','1/2','3/4','7/8','15/16','31/32','1']

derived_freq_bins = numpy.array([-1,0,0.01,0.1,0.5,0.9,0.99,1,2])
derived_virtual_freqs = numpy.arange(0,len(derived_freq_bins)-1)

derived_virtual_xticks = list(derived_virtual_freqs[:-1]+0.5)
derived_virtual_xticklabels = ['0','.01','.1','.5','.9','.99','1']


total_freq_snps = {var_type: numpy.zeros_like(derived_virtual_freqs) for var_type in variant_types}
total_freq_all_snps = numpy.zeros_like(derived_virtual_freqs)
total_null_freq_all_snps = numpy.zeros_like(derived_virtual_freqs)*1.0

total_snp_mutrevs = {'muts': 0, 'revs':0}
total_random_null_snp_mutrevs = {'muts':0, 'revs':0}
total_between_null_snp_mutrevs = {'muts': 0, 'revs':0}

total_gene_gainlosses = {'gains':0, 'losses':0}
total_between_null_gene_gainlosses = {'gains':0, 'losses':0}

gene_freq_bins = numpy.array([-1,0.1,0.9,2])
gene_freq_xticks      = [-3,  -2,   -1,   0,   1,    2,   3]
gene_freq_xticklabels = ['0','0.1','0.9','1','0.9','0.1','0']

total_freq_gains = numpy.zeros(len(gene_freq_bins)-1)*1.0
total_freq_losses = numpy.zeros_like(total_freq_gains)
total_null_freq_losses = numpy.zeros_like(total_freq_gains)

gene_gain_virtual_freqs = numpy.array([2.5,1.5,0.5])
gene_loss_virtual_freqs = numpy.array([-2.5,-1.5,-0.5])

species_snp_change_distribution = {}
species_gene_change_distribution = {}

# observed within host value for twins
pooled_young_twin_snp_change_distribution = []
pooled_young_twin_gene_change_distribution = []

# observed within host value for twins
pooled_twin_snp_change_distribution = []
pooled_twin_gene_change_distribution = []

# observed within host value
pooled_snp_change_distribution = []
pooled_gene_change_distribution = []

pooled_replacement_gene_change_distribution = []

# typical value, median other sample
pooled_between_snp_change_distribution = []
pooled_between_gene_change_distribution = []

# closest other sample
pooled_min_between_snp_change_distribution = []
pooled_min_between_gene_change_distribution = []

replacement_map = {}

countries = ["United States", "United Kingdom", "Western Europe"]

young_twin_output_strs = []
twin_output_strs = []

for species_name in good_species_list:

    sys.stderr.write("\nProcessing %s...\n" % species_name)

    # all samples
    all_samples = sample_order_map.keys()

    # list of samples that meet coverage criteria for this species
    highcoverage_samples = set(diversity_utils.calculate_highcoverage_samples(species_name))
    
    # list of samples that meet QP criteria for this species
    haploid_samples = set(diversity_utils.calculate_haploid_samples(species_name))
    
    #print len(all_samples), len(highcoverage_samples), len(haploid_samples)
       
    if len(haploid_samples) < min_haploid_sample_size:
        continue

    same_sample_idxs, same_subject_idxs, diff_subject_idxs = sample_utils.calculate_ordered_subject_pairs(sample_order_map, all_samples)

    hmp_sample_size = 0        
       
    hmp_samples = set()
    twin_samples = set()
    young_twin_samples = set()
    
    qp_counts = {country:[0,0,0,0] for country in countries}
    
    #print len(all_samples), len(highcoverage_samples), len(haploid_samples)
    #print len(highcoverage_samples - haploid_samples)
    #print len(haploid_samples - highcoverage_samples)

    for sample_pair_idx in xrange(0,len(same_subject_idxs[0])):
   
        i = same_subject_idxs[0][sample_pair_idx]
        j = same_subject_idxs[1][sample_pair_idx]

        sample_i = all_samples[i]
        sample_j = all_samples[j]
        
        country = sample_country_map[sample_i]
        
        #print sample_i, sample_j, country
            
        if country not in countries:
            continue
                
        # Figure out what kind of pair it is

        if not ((sample_i in highcoverage_samples) and (sample_j in highcoverage_samples)):
            
            # In one but not the other
            if ((sample_i in highcoverage_samples) or (sample_j in highcoverage_samples)):
            
                qp_counts[country][0] += 1
            
        else:
            
            # both are highcoverage samples
            
            if (sample_i in haploid_samples) and (sample_j in haploid_samples):
                
                qp_counts[country][1] += 1
                    
                if country==countries[0]:
                        
                    # An HMP pair
                
                    hmp_samples.add(sample_i)
                    hmp_samples.add(sample_j)
                
                    hmp_sample_size += 1
                        
                elif country==countries[1]:
                        
                    # A UK Twin pair
                    twin_samples.add(sample_i)
                    twin_samples.add(sample_j)
                        
                        
                elif country==countries[2]:
                    young_twin_samples.add(sample_i)
                    young_twin_samples.add(sample_j)
                else:
                    pass
                        
            
            elif (sample_i not in haploid_samples) and (sample_j not in haploid_samples):
                
                qp_counts[country][2] += 1
            
            else:
                
                qp_counts[country][3] += 1
        
        
    print qp_counts
        
    sample_size = hmp_sample_size
    if sample_size < min_sample_size:
        continue
            
    species_snp_change_distribution[species_name] = []
    species_gene_change_distribution[species_name] = []
                     
    snp_samples = list(hmp_samples)
    allowed_sample_set = set(hmp_samples)

    twin_samples = list(twin_samples)
    allowed_twin_sample_set = set(twin_samples)
    
    young_twin_samples = list(young_twin_samples)
    allowed_young_twin_sample_set = set(young_twin_samples)
    
        
    combined_samples = (snp_samples+twin_samples+young_twin_samples)
    allowed_combined_sample_set = set(combined_samples)
    
    sys.stderr.write("Proceeding with %d longitudinal comparisons with %d samples!\n" % (sample_size, len(snp_samples)))
    
    import calculate_private_snvs
    private_snv_map = calculate_private_snvs.load_private_snv_map(species_name)
    
    import calculate_snp_prevalences
    snv_freq_map = calculate_snp_prevalences.parse_population_freqs(species_name,polarize_by_consensus=True)
    snv_freq_values = snv_freq_map.values()
    
    import core_gene_utils
    gene_freq_map = core_gene_utils.parse_gene_freqs(species_name)
    gene_freq_values = numpy.array(gene_freq_map.values())
    gene_freq_weights = gene_freq_values*1.0/gene_freq_values.sum()
    
    sys.stderr.write("Loading SFSs for %s...\t" % species_name)
    dummy_samples, sfs_map = parse_midas_data.parse_within_sample_sfs(species_name, allowed_variant_types=set(['1D','2D','3D','4D'])) 
    sys.stderr.write("Done!\n")

    sys.stderr.write("Loading pre-computed substitution rates for %s...\n" % species_name)
    substitution_rate_map = calculate_substitution_rates.load_substitution_rate_map(species_name)
    sys.stderr.write("Calculating SNV matrix...\n")
    dummy_samples, snp_mut_difference_matrix, snp_rev_difference_matrix, snp_mut_opportunity_matrix, snp_rev_opportunity_matrix = calculate_substitution_rates.calculate_mutrev_matrices_from_substitution_rate_map(substitution_rate_map, 'all', allowed_samples=snp_samples)
    snp_samples = dummy_samples
    
    gene_samples, gene_loss_difference_matrix, gene_gain_difference_matrix, gene_loss_opportunity_matrix, gene_gain_opportunity_matrix = calculate_substitution_rates.calculate_mutrev_matrices_from_substitution_rate_map(substitution_rate_map, 'genes', allowed_samples=snp_samples)
    
    gene_difference_matrices = {'gains': gene_gain_difference_matrix, 'losses': gene_loss_difference_matrix}
    gene_opportunity_matrix = gene_loss_opportunity_matrix
    
    opportunity_matrices = {}
    difference_matrices = {}

    
    for var_type in variant_types:
        
        dummy_samples, difference_matrix, opportunity_matrix =    calculate_substitution_rates.calculate_matrices_from_substitution_rate_map(substitution_rate_map, var_type, allowed_samples=snp_samples)
    
        difference_matrices[var_type] = difference_matrix
        opportunity_matrices[var_type] = opportunity_matrix

    difference_matrices['muts'] = snp_mut_difference_matrix
    difference_matrices['revs'] = snp_rev_difference_matrix
    opportunity_matrices['muts'] = snp_mut_opportunity_matrix
    opportunity_matrices['revs'] = snp_rev_opportunity_matrix
    
    snp_difference_matrix = snp_mut_difference_matrix+snp_rev_difference_matrix
    snp_opportunity_matrix = snp_mut_opportunity_matrix+snp_rev_opportunity_matrix
    
    gene_difference_matrix = gene_gain_difference_matrix + gene_loss_difference_matrix
        
    snp_substitution_rate =     snp_difference_matrix*1.0/(snp_opportunity_matrix+(snp_opportunity_matrix==0))
    sys.stderr.write("Done!\n")

    sys.stderr.write("Loading pre-computed temporal changes for %s...\n" % species_name)
    temporal_change_map = calculate_temporal_changes.load_temporal_change_map(species_name)
    sys.stderr.write("Done!\n")

    
    total_snp_modification_map[species_name] = 0
    total_null_snp_modification_map[species_name] = 0
    total_gene_modification_map[species_name] = 0
    total_null_gene_modification_map[species_name] = 0

    temporal_changes = []
    
    ############################
    #
    # Calculate between-twin changes
    #
    ############################
    same_sample_idxs, same_subject_idxs, diff_subject_idxs = sample_utils.calculate_ordered_subject_pairs(sample_order_map, twin_samples)
    
    for sample_pair_idx in xrange(0,len(same_subject_idxs[0])):
    #    
        i = same_subject_idxs[0][sample_pair_idx]
        j = same_subject_idxs[1][sample_pair_idx]
    
        sample_i = twin_samples[i] 
        sample_j = twin_samples[j]
        
        if not ((sample_i in allowed_twin_sample_set) and (sample_j in allowed_twin_sample_set)):
            continue
        
        L, perr, mutations, reversions = calculate_temporal_changes.calculate_mutations_reversions_from_temporal_change_map(temporal_change_map, sample_i, sample_j)
        
        nerr = L*perr
        
        num_mutations = len(mutations)
        num_reversions = len(reversions)
        num_snp_changes = num_mutations+num_reversions
        
        gene_L, gene_perr, gains, losses = calculate_temporal_changes.calculate_gains_losses_from_temporal_change_map(temporal_change_map, sample_i, sample_j)
        
        gene_nerr = gene_L*gene_perr
        num_gains = len(gains)
        num_losses = len(losses)
        num_gene_changes = num_gains+num_losses
        
        if (perr<-0.5) or (gene_perr < -0.5):
            continue
        
        if (nerr > max([0.5, 0.1*num_snp_changes])) or (gene_nerr > max([0.5, 0.1*num_gene_changes])):     
            continue 
            
        pooled_twin_snp_change_distribution.append(num_snp_changes)
        pooled_twin_gene_change_distribution.append(num_gene_changes)    
    
        if num_snp_changes<replacement_difference_threshold:
            twin_output_strs.append("%s, %s: n_snv=%d, n_gene=%d" % (species_name, sample_i, num_snp_changes, num_gene_changes))
        
    
    ############################
    #
    # Calculate between-twin changes for young twins from Korpela et al
    #
    ############################
    same_sample_idxs, same_subject_idxs, diff_subject_idxs = sample_utils.calculate_ordered_subject_pairs(sample_order_map, young_twin_samples)
    
    for sample_pair_idx in xrange(0,len(same_subject_idxs[0])):
    #    
        i = same_subject_idxs[0][sample_pair_idx]
        j = same_subject_idxs[1][sample_pair_idx]
    
        sample_i = young_twin_samples[i] 
        sample_j = young_twin_samples[j]
        
        if not ((sample_i in allowed_young_twin_sample_set) and (sample_j in allowed_young_twin_sample_set)):
            print "Not in right sample set!", sample_i, sample_j
            continue
        
        L, perr, mutations, reversions = calculate_temporal_changes.calculate_mutations_reversions_from_temporal_change_map(temporal_change_map, sample_i, sample_j)
        
        nerr = L*perr
        
        num_mutations = len(mutations)
        num_reversions = len(reversions)
        num_snp_changes = num_mutations+num_reversions
        
        gene_L, gene_perr, gains, losses = calculate_temporal_changes.calculate_gains_losses_from_temporal_change_map(temporal_change_map, sample_i, sample_j)
        
        gene_nerr = gene_L*gene_perr
        num_gains = len(gains)
        num_losses = len(losses)
        num_gene_changes = num_gains+num_losses
        
        if (perr<-0.5): # or (gene_perr < -0.5):
            print "Invalid error rates!", sample_i, sample_j
            continue
        
        if (nerr > max([0.5, 0.1*num_snp_changes])): # or (gene_nerr > max([0.5, 0.1*num_gene_changes])):    
            print "Bad error rates!", sample_i, sample_j, nerr, gene_nerr
            continue 
            
        young_twin_output_strs.append("%s, %s: n_snv=%d, n_gene=%d" % (species_name, sample_i, num_snp_changes, num_gene_changes))
        
        pooled_young_twin_snp_change_distribution.append(num_snp_changes)
        pooled_young_twin_gene_change_distribution.append(num_gene_changes)            
 
 
    ###
    #
    # Now do it for the HMP samples
    #
    ###
    same_sample_idxs, same_subject_idxs, diff_subject_idxs = sample_utils.calculate_ordered_subject_pairs(sample_order_map, snp_samples)
    
    for sample_pair_idx in xrange(0,len(same_subject_idxs[0])):
#    
        i = same_subject_idxs[0][sample_pair_idx]
        j = same_subject_idxs[1][sample_pair_idx]
    
        sample_i = snp_samples[i] 
        sample_j = snp_samples[j]
        
        if not ((sample_i in allowed_sample_set) and (sample_j in allowed_sample_set)):
            continue
        
        L, perr, mutations, reversions = calculate_temporal_changes.calculate_mutations_reversions_from_temporal_change_map(temporal_change_map, sample_i, sample_j)
        
        nerr = L*perr
        
        num_mutations = len(mutations)
        num_reversions = len(reversions)
        num_snp_changes = num_mutations+num_reversions
        
        gene_L, gene_perr, gains, losses = calculate_temporal_changes.calculate_gains_losses_from_temporal_change_map(temporal_change_map, sample_i, sample_j)
        
        gene_nerr = gene_L*gene_perr
        num_gains = len(gains)
        num_losses = len(losses)
        num_gene_changes = num_gains+num_losses
        
        if (perr<-0.5) or (gene_perr < -0.5):
            continue
        
        if (nerr > max([0.5, 0.1*num_snp_changes])) or (gene_nerr > max([0.5, 0.1*num_gene_changes])):
            continue # Only take things with low-ish FPR
    
        
        if not ((sample_i in allowed_sample_set) and (sample_j in allowed_sample_set)):
            continue
        
        pooled_snp_change_distribution.append(num_snp_changes)
                
        species_snp_change_distribution[species_name].append( num_snp_changes)
            
        good_idxs = sample_utils.calculate_samples_in_different_subjects( subject_sample_map, snp_samples, sample_i)

        # typical
        pooled_between_snp_change_distribution.append( choice( snp_difference_matrix[i, good_idxs] ) )
        
        #numpy.median(snp_difference_matrix[i, good_idxs]) )
            
        # minimum
        pooled_min_between_snp_change_distribution.append( snp_difference_matrix[i, good_idxs].min() )
            
            
        if (num_snp_changes>=replacement_difference_threshold):
            sample_pair = (sample_i, sample_j)
            if sample_pair not in replacement_map:
                replacement_map[sample_pair] = []
            replacement_map[sample_pair].append(species_name)
            pooled_replacement_gene_change_distribution.append( num_gene_changes )   
            
        if (num_snp_changes<=modification_difference_threshold):
        
            total_snp_modification_map[species_name] += num_snp_changes
            total_null_snp_modification_map[species_name] += perr
            
            # Count up variant types of observed mutations
            variant_type_counts = {var_type: 0 for var_type in variant_types}
            mutation_variant_type_counts = {var_type: 0 for var_type in variant_types}
            private_variant_type_counts = {var_type: 0 for var_type in variant_types}
            reversion_variant_type_counts = {var_type: 0 for var_type in variant_types}
            
            
            
            for snp_change in mutations:
            
                if snp_change[3] in variant_types:
                    variant_type_counts[snp_change[3]] += 1
                    mutation_variant_type_counts[snp_change[3]] += 1
                else:
                    pass
                    
            for snp_change in reversions:
                if snp_change[3] in variant_types:
                    variant_type_counts[snp_change[3]] += 1
                    reversion_variant_type_counts[snp_change[3]] += 1
                else:
                    pass
            
            for snp_change in (mutations+reversions):        
                gene_name, contig, position, variant_type, A1, D1, A2, D2 = snp_change
                
                f1 = A1*1.0/D2
                f2 = A2*1.0/D2
                
                is_reversion = (f1>f2)
                
                location_tuple = (contig, position)
                
                is_private_snv = (location_tuple in private_snv_map)
                
                if is_private_snv and (variant_type in variant_types):
                    private_variant_type_counts[variant_type] += 1   
                    
                # Now calculate frequency-stratified version
                
                if location_tuple in snv_freq_map:
                    f = snv_freq_map[location_tuple]
                else:
                    sys.stderr.write("SNP not in map. Shouldn't happen!\n")
                    f = -0.5
                
                # Let's impose that private snvs have zero freq        
                if is_private_snv:
                    f = -0.5
                
                # Change f so that it represents
                # frequency of allele at second timepoint
                if is_reversion:
                    f = 1-f
                
                f_idx = ((f>derived_freq_bins[:-1])*(f<=derived_freq_bins[1:])).argmax()    
                
                if variant_type in variant_types:
                    total_freq_snps[variant_type][f_idx] += 1
                total_freq_all_snps[f_idx] += 1
                
                # Now draw a null from genome
                L = snp_opportunity_matrix[i,j]
                L_snv = len(snv_freq_map) # A slight overestimate
                snv_fraction = L_snv*1.0/L
                num_bootstraps = 10
                for bootstrap_idx in xrange(0,num_bootstraps):
                    
                    if random()<snv_fraction:
                        f = snv_freq_values[randint(0,len(snv_freq_values))]
                        rev_f = 1-f
                    
                        f_idx = ((f>derived_freq_bins[:-1])*(f<=derived_freq_bins[1:])).argmax()    
                    
                        rev_f_idx = ((rev_f>derived_freq_bins[:-1])*(rev_f<=derived_freq_bins[1:])).argmax()  
                    
                        total_null_freq_all_snps[f_idx] += (1-f)*1.0/num_bootstraps
                        total_null_freq_all_snps[rev_f_idx] += f*1.0/num_bootstraps
                        
                    else:
                    
                        total_null_freq_all_snps[0] += 1.0/num_bootstraps
                           
            # Add them to running total
            # & form running total for this sample only
            observed_sample_size = 0
            for var_type in variant_types:
                total_snps[var_type] += variant_type_counts[var_type]
                total_mutation_snps[var_type] += mutation_variant_type_counts[var_type]
                total_private_snps[var_type] += private_variant_type_counts[var_type]
                total_reversion_snps[var_type] += reversion_variant_type_counts[var_type]
                observed_sample_size += variant_type_counts[var_type]
            
            # Opportunities for these two samples
            total_opportunities = sum([opportunity_matrices[var_type][i,j] for var_type in variant_types])
            
            for var_type in variant_types:
                total_random_null_snps[var_type] += observed_sample_size*opportunity_matrices[var_type][i,j]*1.0/total_opportunities  
            
            # Now get a null from between-host changes
            good_comparison_idxs = (snp_opportunity_matrix[i,:]>0.5)
            good_comparison_idxs *= sample_utils.calculate_samples_in_different_subjects( subject_sample_map, snp_samples, sample_i)

            
            total_between_host_changes = sum([numpy.median(difference_matrices[var_type][i,good_comparison_idxs]) for var_type in variant_types])
            
            for var_type in variant_types:
                total_between_null_snps[var_type] += observed_sample_size*(numpy.median(difference_matrices[var_type][i, good_comparison_idxs]))*1.0/total_between_host_changes  
            
            # Now do the same thing, except for SNV mutations/reversions
            
            # Tally mutations & reversions    
            total_snp_mutrevs['muts'] += num_mutations
            total_snp_mutrevs['revs'] += num_reversions
            
            observed_sample_size = num_mutations+num_reversions
            
            # Now get a null from randomly drawing from genome
            total_opportunities = opportunity_matrices['muts'][i,j]+opportunity_matrices['revs'][i,j]
            for type in ['muts','revs']:
                total_random_null_snp_mutrevs[type] += observed_sample_size*opportunity_matrices[type][i,j]*1.0/total_opportunities
                
            # Now get a null from between-host changes
            total_between_host_changes = sum([numpy.median(difference_matrices[type][i,good_comparison_idxs]) for type in ['muts','revs']])
            for type in ['muts','revs']:
                total_between_null_snp_mutrevs[type] += observed_sample_size*(numpy.median(difference_matrices[type][i,good_comparison_idxs]))*1.0/total_between_host_changes  
            
            # If there are gene changes to look at: 
            if num_gene_changes > -0.5:
            
                for gene_change in gains:
                    gene_name = gene_change[0]
                    
                    if gene_name in gene_freq_map:
                        f = gene_freq_map[gene_name]
                    else:
                        f = 0
                    
                    f_idx = ((f>gene_freq_bins[:-1])*(f<=gene_freq_bins[1:])).argmax()    
                    total_freq_gains[f_idx] += 1
                    
                        
                    
                for gene_change in losses:
                    gene_name = gene_change[0]
                    
                    if gene_name in gene_freq_map:
                        f = gene_freq_map[gene_name]
                    else:
                        f = 0
                    
                    f_idx = ((f>gene_freq_bins[:-1])*(f<=gene_freq_bins[1:])).argmax()    
                    total_freq_losses[f_idx] += 1
                    
                    num_bootstraps = 10
                    fs = choice(gene_freq_values, size=num_bootstraps, p=gene_freq_weights)
                    for f in fs:
                        f_idx = ((f>gene_freq_bins[:-1])*(f<=gene_freq_bins[1:])).argmax()    
                        total_null_freq_losses[f_idx] += 1.0/num_bootstraps
                
            
            
                gene_i = i
                #print sample_i, sample_j
                good_comparison_idxs = (gene_opportunity_matrix[i,:]>0.5)
                
                #print gene_opportunity_matrix[i,:].shape
                #print good_comparison_idxs.sum() 
                 
                good_comparison_idxs *= sample_utils.calculate_samples_in_different_subjects( subject_sample_map, gene_samples, sample_i)
                
                #print good_comparison_idxs.sum() 
                            
                pooled_gene_change_distribution.append(num_gene_changes) 
                species_gene_change_distribution[species_name].append(num_gene_changes) 
                
                # Typical value
                pooled_between_gene_change_distribution.append( numpy.median(gene_difference_matrix[gene_i, good_comparison_idxs]) )      
                # Minimum value
                pooled_min_between_gene_change_distribution.append( gene_difference_matrix[gene_i, good_comparison_idxs].min() )
                
                total_gene_modification_map[species_name] += num_gene_changes
                total_null_gene_modification_map[species_name] += gene_perr   
                
                total_gene_gainlosses['gains'] += num_gains
                total_gene_gainlosses['losses'] += num_losses
                
                if num_gene_changes>0.5: # you actually have some genes to draw a null ffrom...
                    
                    if good_comparison_idxs.sum() < 0.5:
                        print sample_i, "no gene comparisons!" 
                    else:
                        observed_sample_size = num_gains+num_losses
                
                        # Now get a null from between-host changes
                        total_between_host_changes = sum([numpy.median(gene_difference_matrices[type][gene_i,good_comparison_idxs]) for type in ['gains','losses']])
                
                        if total_between_host_changes < 0.5:
                            print sample_i, num_gene_changes, gene_difference_matrices['gains'][gene_i,:], gene_difference_matrices['losses'][gene_i,:], gene_difference_matrices['losses'][gene_i,:], gene_opportunity_matrix[gene_i,:]
                
                        for type in ['gains','losses']:
                            total_between_null_gene_gainlosses[type] += observed_sample_size*(numpy.median(gene_difference_matrices[type][gene_i,:]))*1.0/total_between_host_changes  
    
        temporal_changes.append((sample_i, sample_j, num_snp_changes, num_gene_changes))        

    num_passed_species+=1
    
    if len(temporal_changes) > 0:
        num_temporal_change_map[species_name] = temporal_changes
        
    sys.stderr.write("Done with %s!\n" % species_name) 
    
sys.stderr.write("Done looping over species!\n")    

print "Young twins:"
print "\n".join(young_twin_output_strs)

print "All twins:"
print "\n".join(twin_output_strs)


species_names = []
sample_sizes = []


for species_name in num_temporal_change_map.keys():
    species_names.append(species_name)
    sample_sizes.append( len(num_temporal_change_map[species_name]) )
    
# sort in descending order of sample size
# Sort by num haploids    
sample_sizes, species_names = zip(*sorted(zip(sample_sizes, species_names),reverse=True))
     
sys.stderr.write("Postprocessing %d species!\n" % len(species_names))


cmap_str = 'YlGnBu'
vmin = -2
vmax = 3
cmap = pylab.get_cmap(cmap_str) 

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=-1):
    if n == -1:
        n = cmap.N
    new_cmap = mcolors.LinearSegmentedColormap.from_list(
         'trunc({name},{a:.2f},{b:.2f})'.format(name=cmap.name, a=minval, b=maxval),
         cmap(numpy.linspace(minval, maxval, n)))
    return new_cmap

cmap = truncate_colormap(cmap, 0.25, 1.0)
cNorm  = colors.Normalize(vmin=0, vmax=vmax)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cmap)



####################################################
#
# Set up Figure (1 panels, arranged in 1x1 grid)
#
####################################################

####################################################
#
# Set up Figure (1 panels, arranged in 1x1 grid)
#
####################################################

pylab.figure(1,figsize=(5,3.5))
fig = pylab.gcf()
# make three panels panels
outer_grid  = gridspec.GridSpec(1,2,width_ratios=[50,1],wspace=0.05)

change_axis = plt.Subplot(fig, outer_grid[0])
fig.add_subplot(change_axis)

change_axis.set_ylabel('Number of samples')
change_axis.set_ylim([-75,75])
change_axis.set_xlim([-1,len(species_names)])
change_axis.plot([-1,len(species_names)],[0,0],'k-')

change_axis.set_yticks([-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70])
change_axis.set_yticklabels(['70','60','50','40','30','20','10','0','10','20','30','40','50','60','70'])

xticks = numpy.arange(0,len(species_names))
#xticklabels = ["%s (%d)" % (species_names[i],sample_sizes[i]) for i in xrange(0,len(sample_sizes))]
xticklabels = ["%s" % (species_names[i]) for i in xrange(0,len(sample_sizes))]

change_axis.set_xticks(xticks)
change_axis.set_xticklabels(xticklabels, rotation='vertical',fontsize=4)

#change_axis.spines['top'].set_visible(False)
#change_axis.spines['right'].set_visible(False)
change_axis.get_xaxis().tick_bottom()
change_axis.get_yaxis().tick_left()

cax = plt.Subplot(fig, outer_grid[1])
fig.add_subplot(cax)


##############
#
# Real figure
#
###############
pylab.figure(2,figsize=(7,4))
fig2 = pylab.gcf()
# make three panels panels
outer_outer_grid_2 = gridspec.GridSpec(1,1) #2, width_ratios=[1,0.2],wspace=0.1) 
outer_grid_2 = gridspec.GridSpecFromSubplotSpec(2,1, height_ratios=[1,0.7],hspace=0.6, subplot_spec=outer_outer_grid_2[0])

pylab.figure(3,figsize=(5,1.5))
fig3 = pylab.gcf()
outer_grid_3 = gridspec.GridSpec(1,1)

prevalence_grid = gridspec.GridSpecFromSubplotSpec(1,3, width_ratios=[1,1,0.3],wspace=0.5,subplot_spec=outer_grid_2[1])

upper_grid = gridspec.GridSpecFromSubplotSpec(1,4, width_ratios=[1,1,1,0.6],wspace=0.45,subplot_spec=outer_grid_3[0])

dnds_axis = plt.Subplot(fig3, upper_grid[0])
fig3.add_subplot(dnds_axis)
dnds_axis.set_ylabel('# changes')

dnds_axis.spines['top'].set_visible(False)
dnds_axis.spines['right'].set_visible(False)
dnds_axis.get_xaxis().tick_bottom()
dnds_axis.get_yaxis().tick_left()

dnds_axis.set_xlim([0.3,2.7])
dnds_axis.set_xticks([1,2])
dnds_axis.set_xticklabels(['non','syn'])
dnds_axis.set_ylim([0,300])
dnds_axis.set_yticks([0,100,200,300])

# TODO: significance of DNDS < 1!

# Mutation / reversion
mutrev_axis = plt.Subplot(fig3, upper_grid[1])
fig3.add_subplot(mutrev_axis)

mutrev_axis.spines['top'].set_visible(False)
mutrev_axis.spines['right'].set_visible(False)
mutrev_axis.get_xaxis().tick_bottom()
mutrev_axis.get_yaxis().tick_left()

mutrev_axis.set_xlim([0.3,2.7])
#mutrev_axis.set_ylim([0,645])
mutrev_axis.set_yticks([0,200,400,600])
mutrev_axis.set_xticks([1,2])
#mutrev_axis.set_xticklabels(['mut','rev'])
mutrev_axis.set_xticklabels(['away\nfrom\nref','toward\nref'],fontsize=6)
#mutrev_axis.set_xticklabels(['alt','ref'],fontsize=6)

#mutrev_axis.set_yticklabels([])



# Gain / loss
gainloss_axis = plt.Subplot(fig3, upper_grid[2])
fig3.add_subplot(gainloss_axis)
gainloss_axis.spines['top'].set_visible(False)
gainloss_axis.spines['right'].set_visible(False)
gainloss_axis.get_xaxis().tick_bottom()
gainloss_axis.get_yaxis().tick_left()
gainloss_axis.set_xlim([0.3,2.7])
#gainloss_axis.set_ylim([0,2100])
gainloss_axis.set_yticks([0,1000,2000,3000])
gainloss_axis.set_xticks([1,2])
gainloss_axis.set_xticklabels(['loss','gain'])
#gainloss_axis.set_yticklabels([])


legend_axis = plt.Subplot(fig3, upper_grid[3])
fig3.add_subplot(legend_axis)

legend_axis.set_ylim([0,1])
legend_axis.set_xlim([0,1])

legend_axis.spines['top'].set_visible(False)
legend_axis.spines['right'].set_visible(False)
legend_axis.spines['left'].set_visible(False)
legend_axis.spines['bottom'].set_visible(False)

legend_axis.set_xticks([])
legend_axis.set_yticks([])


pooled_grid = gridspec.GridSpecFromSubplotSpec(1,3,width_ratios=[1,1,0.2],wspace=0.15,subplot_spec=outer_grid_2[0])

pooled_snp_axis = plt.Subplot(fig2, pooled_grid[0])
fig2.add_subplot(pooled_snp_axis)
pooled_snp_axis.set_ylabel('Fraction comparisons $\geq n$')
pooled_snp_axis.set_xlabel('# SNV changes')
#pooled_axis.set_ylim([-35,35])
#pooled_snp_axis.set_xlim([2e-01,1e05])
pooled_snp_axis.set_xlim([1,1e05])

pooled_snp_axis.set_xticklabels([])

pooled_snp_axis.spines['top'].set_visible(False)
pooled_snp_axis.spines['right'].set_visible(False)
pooled_snp_axis.get_xaxis().tick_bottom()
pooled_snp_axis.get_yaxis().tick_left()
 
pooled_gene_axis = plt.Subplot(fig2, pooled_grid[1])
fig2.add_subplot(pooled_gene_axis)
#pooled_gene_axis.set_ylabel('Number of samples')
pooled_gene_axis.set_xlabel('# gene changes')
#pooled_axis.set_ylim([-35,35])
#pooled_gene_axis.set_xlim([2e-01,1e04])
pooled_gene_axis.set_xlim([1,1e04])

pooled_gene_axis.spines['top'].set_visible(False)
pooled_gene_axis.spines['right'].set_visible(False)
pooled_gene_axis.get_xaxis().tick_bottom()
pooled_gene_axis.get_yaxis().tick_left()
 
legend2_axis = plt.Subplot(fig2, pooled_grid[2])
fig2.add_subplot(legend2_axis)

legend2_axis.set_ylim([0,1])
legend2_axis.set_xlim([0,1])

legend2_axis.spines['top'].set_visible(False)
legend2_axis.spines['right'].set_visible(False)
legend2_axis.spines['left'].set_visible(False)
legend2_axis.spines['bottom'].set_visible(False)

legend2_axis.set_xticks([])
legend2_axis.set_yticks([])

legend2_axis.plot([-2,-1],[-1,-1],'-',linewidth=1, color='#08519c',label='Within-host')
legend2_axis.plot([-2,-1],[-1,-1],'-',color='#08519c',linewidth=1, label='modification',zorder=2,path_effects=[pe.Stroke(linewidth=5, foreground='#9ecae1'), pe.Normal()])
legend2_axis.plot([-2,-1],[-1,-1],'-',color='#08519c', label='replacement',linewidth=1,path_effects=[pe.Stroke(linewidth=5, foreground='#fee0d2'), pe.Normal()])
legend2_axis.plot([-2,-1],[-1,-1], '-',linewidth=1,color='w', alpha=0.5, label=' ')
legend2_axis.plot([-2,-1],[-1,-1], '-',linewidth=1,color='r', alpha=0.5, label='Between-host\n(unrelated)')
legend2_axis.plot([-2,-1],[-1,-1],'-',linewidth=1,color='#8856a7', label='Between-host\n(adult twins)')

legend2_axis.legend(loc='upper center',frameon=False,fontsize=5,numpoints=1,ncol=1,handlelength=1)   


frequency_axis = plt.Subplot(fig2, prevalence_grid[0])
fig2.add_subplot(frequency_axis)

frequency_axis.spines['top'].set_visible(False)
frequency_axis.spines['right'].set_visible(False)
frequency_axis.get_xaxis().tick_bottom()
frequency_axis.get_yaxis().tick_left()
 
frequency_axis.set_xlabel('Derived allele prevalence\nacross hosts')
frequency_axis.set_ylabel('# SNV changes')

frequency_axis.set_xticks(derived_virtual_xticks)
frequency_axis.set_xticklabels(derived_virtual_xticklabels) #,rotation='vertical')

frequency_axis.set_ylim([0,200])

gene_frequency_axis = plt.Subplot(fig2, prevalence_grid[1])
fig2.add_subplot(gene_frequency_axis)

gene_frequency_axis.spines['top'].set_visible(False)
gene_frequency_axis.spines['right'].set_visible(False)
gene_frequency_axis.get_xaxis().tick_bottom()
gene_frequency_axis.get_yaxis().tick_left()
 
gene_frequency_axis.set_xlabel('Gene prevalence across hosts')
gene_frequency_axis.set_ylabel('# gene changes')

gene_frequency_axis.set_xlim([gene_freq_xticks[0],gene_freq_xticks[-1]])
gene_frequency_axis.set_xticks(gene_freq_xticks)
gene_frequency_axis.set_xticklabels(gene_freq_xticklabels) #,rotation='vertical')

gene_frequency_axis.plot([0,0],[100,100],'k-')
gene_frequency_axis.set_ylim([0,100])

gene_legend_axis = plt.Subplot(fig2, prevalence_grid[2])
fig2.add_subplot(gene_legend_axis)

gene_legend_axis.set_ylim([0,1])
gene_legend_axis.set_xlim([0,1])

gene_legend_axis.spines['top'].set_visible(False)
gene_legend_axis.spines['right'].set_visible(False)
gene_legend_axis.spines['left'].set_visible(False)
gene_legend_axis.spines['bottom'].set_visible(False)

gene_legend_axis.set_xticks([])
gene_legend_axis.set_yticks([])

gene_legend_axis.bar([-2],[-1],width=0.2, linewidth=0,facecolor='#b3de69',label='gain')
gene_legend_axis.bar([-2],[-1],width=0.2, linewidth=0,facecolor='#ff7f00',label='loss')
gene_legend_axis.bar([-2],[-1],width=0.2, linewidth=0,facecolor='0.7',label='de novo\nexpectation')

gene_legend_axis.legend(loc='center left',frameon=False,fontsize=5,numpoints=1,ncol=1,handlelength=1)   


pylab.figure(4,figsize=(2,2))
fig4 = pylab.gcf()
outer_grid_4 = gridspec.GridSpec(1,1)

avg_axis = plt.Subplot(fig4, outer_grid_4[0])
fig4.add_subplot(avg_axis)

#avg_axis.set_ylabel('Microbiome-wide avg')
avg_axis.set_xlim([0.5,2.5])
avg_axis.set_xticks([])

avg_axis.spines['top'].set_visible(False)
avg_axis.spines['right'].set_visible(False)
avg_axis.get_xaxis().tick_bottom()
avg_axis.get_yaxis().tick_left()


##############################################################################
#
# Plot results
#
##############################################################################

 
for species_idx in xrange(0,len(species_names)):

    species_name = species_names[species_idx]
    temporal_changes = num_temporal_change_map[species_name]
    
    
    if total_snp_modification_map[species_name] > 0:
        snp_is_significant = ((total_null_snp_modification_map[species_name]*1.0/total_snp_modification_map[species_name]) < 0.1)
    else:
        snp_is_significant = 0
    
    if total_gene_modification_map[species_name] > 0:
        gene_is_significant = ((total_null_gene_modification_map[species_name]*1.0/total_gene_modification_map[species_name]) < 0.1)
    else:
        gene_is_significant = 0
    
    print species_name, total_snp_modification_map[species_name], total_null_snp_modification_map[species_name], snp_is_significant, total_gene_modification_map[species_name], total_null_gene_modification_map[species_name], gene_is_significant
    
     
    snp_changes = []
    gene_changes = []
    for sample_1, sample_2, num_snps, num_genes in temporal_changes:
         
         snp_changes.append(num_snps)
         
         if num_genes>=0:
             gene_changes.append(num_genes)
             
    snp_changes = numpy.array(snp_changes)
    snp_changes.sort()
    gene_changes = numpy.array(gene_changes)
    gene_changes.sort()
    
    for idx in xrange(0,len(snp_changes)):
        
        if snp_changes[idx]<0.5:
            colorVal='0.7'
        else:
        
            colorVal = scalarMap.to_rgba(log10(snp_changes[idx]))
        
        change_axis.fill_between([species_idx-0.3,species_idx+0.3], [idx,idx],[idx+1.05,idx+1.05],color=colorVal,linewidth=0)
    
        
        if snp_is_significant:
            
            change_axis.text(species_idx, len(snp_changes),'*',fontsize=4)
    
    for idx in xrange(0,len(gene_changes)):
        
        if gene_changes[idx]<0.5:
            colorVal='0.7'
        else:
            colorVal = scalarMap.to_rgba(log10(gene_changes[idx]))
        
        change_axis.fill_between([species_idx-0.3,species_idx+0.3], [-idx-1.05,-idx-1.05],[-idx,-idx],color=colorVal,linewidth=0)
        
        if gene_is_significant:
            
            change_axis.text(species_idx, -len(gene_changes)-3,'*',fontsize=4)
            

m = change_axis.scatter([200],[1],c=[0], vmin=0, vmax=vmax, cmap=cmap, marker='^')


cbar = fig.colorbar(m,cax=cax,orientation='vertical', ticks=[0,1,2,3])
cbar.set_ticklabels(['$1$','$10$','$10^{2}$','$10^{3}$'])
cbar.set_label('Number of changes',rotation=270,labelpad=10)
cl = pylab.getp(cbar.ax, 'ymajorticklabels')
pylab.setp(cl, fontsize=9) 
#fig.text(0.945,0.05,'$\\pi/\\pi_0$',fontsize=12)

cbar.ax.tick_params(labelsize=5)
change_axis.text(20,25,'SNVs',fontsize=5)
change_axis.text(20,-20,'Genes',fontsize=5)

######
#
# Distribution of nucleotide changes
#
######

pooled_snp_change_distribution = numpy.array(pooled_snp_change_distribution)
pooled_twin_snp_change_distribution = numpy.array(pooled_twin_snp_change_distribution)
pooled_between_snp_change_distribution = numpy.array(pooled_between_snp_change_distribution)
pooled_min_between_snp_change_distribution = numpy.array(pooled_min_between_snp_change_distribution)

print pooled_young_twin_snp_change_distribution

print "Mean within host snps =", pooled_snp_change_distribution.mean()
print "Median withon host snps =", numpy.median(pooled_snp_change_distribution)

num_replacements = (pooled_snp_change_distribution>replacement_difference_threshold).sum()
num_twin_replacements = (pooled_twin_snp_change_distribution>1000).sum()


print num_replacements, "replacements", "(%g)" % (num_replacements*1.0/len(pooled_snp_change_distribution))

print num_twin_replacements, "old twin replacements", "(%g)" % (num_twin_replacements*1.0/len(pooled_twin_snp_change_distribution))




pooled_snp_change_distribution = numpy.clip(pooled_snp_change_distribution, 1e-06,1e08)
pooled_twin_snp_change_distribution = numpy.clip(pooled_twin_snp_change_distribution, 1e-06,1e08)
pooled_between_snp_change_distribution = numpy.clip(pooled_between_snp_change_distribution, 1e-06,1e08)
pooled_min_between_snp_change_distribution = numpy.clip(pooled_min_between_snp_change_distribution, 1e-06,1e08)

# Plot SNP change averages

print pooled_snp_change_distribution.mean()
print pooled_between_snp_change_distribution.mean()

avg_axis.bar([0.75],[pooled_snp_change_distribution.mean()],width=0.5,facecolor='#08519c',linewidth=0,log=True)
avg_axis.bar([1.75],[pooled_between_snp_change_distribution.mean()],width=0.5, facecolor='r',alpha=0.5,linewidth=0,log=True)
#avg_axis.semilogy([-1],[0],'k.')

avg_axis.set_ylim([1,1e05])

fig4.savefig('%s/supplemental_within_between_avg.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight',transparent=True)


# Now start plotting the SNP change distribution

# First within hosts

xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(pooled_snp_change_distribution, min_x=1e-02, max_x=1e09)

pooled_snp_axis.step(xs,ns/ns[0],'-',color='#08519c',linewidth=1, label=('Within-host (n=%d)' % ns[0]), where='mid',zorder=4)

# This lets you set the y range

ymin = 1.0/ns[0]
ymax = 1.3

# Put on log scale
pooled_snp_axis.loglog([0.1],[1],'k.')
pooled_gene_axis.loglog([0.1],[1],'k.')
pooled_gene_axis.set_yticklabels([])

pooled_snp_axis.set_ylim([ymin,ymax])

# Save intermediate version (for Keynote animations)
fig2.savefig('%s/figure_6.1.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight',transparent=True)

# Now the typical between host

# Random between host
xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(pooled_between_snp_change_distribution, min_x=1e-02, max_x=1e09)

# Min between host
#xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(pooled_min_between_snp_change_distribution, min_x=1e-02, max_x=1e09)

pooled_snp_axis.step(xs,ns/ns[0],'-',color='r',linewidth=0.5, alpha=0.5, label='Between-host', where='mid',zorder=2)

# Save intermediate version (for Keynote animations)
fig2.savefig('%s/figure_6.2.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight',transparent=True)

# Now do twins

xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(pooled_twin_snp_change_distribution, min_x=1e-02, max_x=1e09)

pooled_snp_axis.step(xs,ns/ns[0],'-',color='#8856a7',linewidth=1, label=('Twins (n=%d)' % ns[0]), where='mid',zorder=2)

# Save intermediate version (for Keynote animations)
fig2.savefig('%s/figure_6.3.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight',transparent=True)

# Now fill in the graphics

pooled_snp_axis.fill_between([1e-01,modification_difference_threshold],[ymin,ymin],[ymax,ymax],color='#deebf7',zorder=1)
pooled_snp_axis.fill_between([replacement_difference_threshold,1e05],[ymin,ymin],[ymax,ymax],color='#fee0d2',zorder=1)

pooled_snp_axis.text(exp((log(1e05)+log(replacement_difference_threshold))/2), ymax*1.2, 'putative\nreplacement',fontsize=6,fontstyle='italic',ha='center',color='#fc9272',zorder=1)
pooled_snp_axis.text(exp((log(1)+log(modification_difference_threshold))/2), ymax*1.2, 'putative\nmodification',fontsize=6,fontstyle='italic',ha='center',color='#9ecae1',zorder=1)
#pooled_snp_axis.text(exp((log(modification_difference_threshold)+log(replacement_difference_threshold))/2), ymax*1.2, 'unclassified',fontsize=6,fontstyle='italic',ha='center')

# Now do same thing for genes

pooled_twin_gene_change_distribution = numpy.array(pooled_twin_gene_change_distribution)
pooled_gene_change_distribution = numpy.array(pooled_gene_change_distribution)
pooled_replacement_gene_change_distribution = numpy.array(pooled_replacement_gene_change_distribution)
pooled_between_gene_change_distribution = numpy.array(pooled_between_gene_change_distribution)
pooled_min_between_gene_change_distribution = numpy.array(pooled_min_between_gene_change_distribution)

pooled_gene_change_distribution = numpy.clip(pooled_gene_change_distribution, 1e-01,1e08)
pooled_replacement_gene_change_distribution = numpy.clip(pooled_replacement_gene_change_distribution, 1e-01,1e08)

pooled_twin_gene_change_distribution = numpy.clip(pooled_twin_gene_change_distribution, 1e-01,1e08)
pooled_between_gene_change_distribution = numpy.clip(pooled_between_gene_change_distribution, 1e-01,1e08)
pooled_min_between_gene_change_distribution = numpy.clip(pooled_min_between_gene_change_distribution, 1e-01,1e08)

xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(pooled_gene_change_distribution, min_x=1e-02, max_x=1e09)

print "regular", xs, ns

pooled_gene_axis.step(xs,ns/ns[0],'-',color='#08519c',linewidth=1, label='Within-host',zorder=2,where='mid',path_effects=[pe.Stroke(linewidth=5, foreground='#9ecae1'), pe.Normal()])

pooled_gene_axis.loglog([1e-01,1e05],[1.0/ns[0],1.0/ns[0]],'k:')

pooled_gene_axis.set_ylim([1.0/ns[0],1.3])
pooled_gene_axis.set_yticklabels([])
#pooled_snp_axis.set_yticklabels(['0.01','0.1','1'])


xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(pooled_replacement_gene_change_distribution, min_x=1e-02, max_x=1e09)

print "replacement", xs, ns

pooled_gene_axis.step(xs,ns/ns[0],'-',color='#08519c',linewidth=1, label='Within-host',zorder=1,where='mid',path_effects=[pe.Stroke(linewidth=5, foreground='#fee0d2'), pe.Normal()])



xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(pooled_twin_gene_change_distribution, min_x=1e-02, max_x=1e09)

pooled_gene_axis.step(xs,ns/ns[0],'-',color='#8856a7',linewidth=1, label='Twin',zorder=1,where='mid')

xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(pooled_min_between_gene_change_distribution, min_x=1e-02, max_x=1e09)

pooled_gene_axis.step(xs,ns/ns[0],'-',color='r',linewidth=0.5, label='Between-host',zorder=0,alpha=0.5, where='mid')


#pooled_gene_axis.legend(loc='lower left', frameon=False, fontsize=5, numpoints=1, handlelength=1)

# Plot dNdS and expected version

observed_nonsynonymous = total_snps['1D']
expected_nonsynonymous = total_snps['4D']/total_random_null_snps['4D']*total_random_null_snps['1D']

within_dnds = observed_nonsynonymous/expected_nonsynonymous
print "Within-host dNdS =", within_dnds


observed_nonsynonymous = total_between_null_snps['1D']
expected_nonsynonymous = total_between_null_snps['4D']/total_random_null_snps['4D']*total_random_null_snps['1D']

between_dnds = observed_nonsynonymous/expected_nonsynonymous
print "Between-host dNdS =", between_dnds


observed_totals = numpy.array([total_snps[var_type] for var_type in variant_types])*1.0
random_totals = numpy.array([total_random_null_snps[var_type] for var_type in variant_types])*1.0 
between_totals = numpy.array([total_between_null_snps[var_type] for var_type in variant_types])*1.0 

print observed_totals
print random_totals
print between_totals

legend_axis.bar([-2],[-1],width=0.2, linewidth=0,color='#08519c',label='Within-host\n(modification)')
legend_axis.bar([-2],[-1],width=0.2, linewidth=0,color='r', alpha=0.5, label='Between-host\n(unrelated)')
legend_axis.bar([-2],[-1],width=0.2, linewidth=0,color='#8856a7', label='Between-host\n(adult twins)')
legend_axis.bar([-2],[-1],width=0.2, linewidth=0,color='k', alpha=0.5, label='De novo\nexpectation')

legend_axis.legend(loc='upper center',frameon=False,fontsize=5,numpoints=1,ncol=1,handlelength=1)   



dnds_axis.bar(numpy.arange(1,3)-0.1, observed_totals, width=0.2, linewidth=0, color='#08519c',label='Obs')

dnds_axis.bar(numpy.arange(1,3)-0.3, random_totals, width=0.2, linewidth=0, color='k',alpha=0.5,label='Null 1')

dnds_axis.bar(numpy.arange(1,3)+0.1, between_totals, width=0.2, linewidth=0, color='r',alpha=0.5,label='Null 2')

#dnds_axis.legend(loc='upper right',frameon=False, handlelength=1)

# Plot mutations, reversions

mutrev_axis.bar([1-0.3, 2-0.3], [total_random_null_snp_mutrevs['muts'], total_random_null_snp_mutrevs['revs']], width=0.2, linewidth=0, color='k',alpha=0.5,label='Null 1')
 
mutrev_axis.bar([1-0.1, 2-0.1], [total_snp_mutrevs['muts'], total_snp_mutrevs['revs']], width=0.2, linewidth=0,color='#08519c')

mutrev_axis.bar([1+0.1, 2+0.1], [total_between_null_snp_mutrevs['muts'], total_between_null_snp_mutrevs['revs']], width=0.2, linewidth=0, color='r',alpha=0.5,label='Null 2')


# Plot gains and losses

gainloss_axis.bar([1-0.3, 2-0.3], [total_gene_gainlosses['losses']+total_gene_gainlosses['gains'],0], width=0.2, linewidth=0,color='k',alpha=0.5)

gainloss_axis.bar([1-0.1, 2-0.1], [total_gene_gainlosses['losses'], total_gene_gainlosses['gains']], width=0.2, linewidth=0,color='#08519c')

gainloss_axis.bar([1+0.1, 2+0.1], [total_between_null_gene_gainlosses['losses'],total_between_null_gene_gainlosses['gains']], width=0.2, linewidth=0, color='r',alpha=0.5)

print "Printing replacement map!"
for sample_pair in replacement_map.keys():
    print sample_pair, len(replacement_map[sample_pair]), replacement_map[sample_pair]


print 'Reversions', total_reversion_snps, total_reversion_snps['1D']*1.0/(total_reversion_snps['4D']+(total_reversion_snps['4D']==0))/(total_random_null_snps['1D']*1.0/total_random_null_snps['4D']) 

print 'Mutations', total_mutation_snps, total_mutation_snps['1D']*1.0/(total_mutation_snps['4D']+(total_mutation_snps['4D']==0))/(total_random_null_snps['1D']*1.0/total_random_null_snps['4D']) 

print 'Private', total_private_snps, total_private_snps['1D']*1.0/(total_private_snps['4D']+(total_private_snps['4D']==0))/(total_random_null_snps['1D']*1.0/total_random_null_snps['4D']) 

print "All", total_snps, total_snps['1D']*1.0/total_snps['4D']/(total_random_null_snps['1D']*1.0/total_random_null_snps['4D']) 

print "Random", total_random_null_snps, 1.0

print "Between", total_between_null_snps, total_between_null_snps['1D']*1.0/total_between_null_snps['4D'] /(total_random_null_snps['1D']*1.0/total_random_null_snps['4D']) 

print '1D', total_freq_snps['1D']
print '4D', total_freq_snps['4D']
print 'all', total_freq_all_snps
print 'null', total_null_freq_all_snps

frequency_axis.bar(derived_virtual_freqs, total_freq_snps['4D'],width=0.3,linewidth=0,facecolor='#b3de69',label='syn (4D)',zorder=3)

frequency_axis.bar(derived_virtual_freqs, total_freq_snps['1D']+total_freq_snps['4D'],width=0.3,linewidth=0,facecolor='#ff7f00',label='non (1D)',zorder=2)

frequency_axis.bar(derived_virtual_freqs, total_freq_all_snps,width=0.3,linewidth=0,facecolor='#b15928',label='(2D & 3D)',zorder=1)


frequency_axis.bar(derived_virtual_freqs-0.3, total_null_freq_all_snps,width=0.3,linewidth=0,facecolor='0.7',label='de novo\nexpectation',zorder=0)

frequency_axis.legend(loc='upper center',frameon=False,fontsize=5,numpoints=1,ncol=2,handlelength=1)   

gene_frequency_axis.bar(gene_gain_virtual_freqs, total_freq_gains,width=0.3,linewidth=0,facecolor='#b3de69',label='gain')

gene_frequency_axis.bar(gene_loss_virtual_freqs, total_freq_losses,width=0.3,linewidth=0, facecolor='#ff7f00',label='loss')

gene_frequency_axis.bar(gene_loss_virtual_freqs-0.3, total_null_freq_losses,width=0.3,linewidth=0, facecolor='0.7',label='de novo\nexpectation')

#gene_frequency_axis.legend(loc='upper center',frameon=False,fontsize=5,numpoints=1,ncol=3,handlelength=1)   

sys.stderr.write("Saving figure...\t")
fig2.savefig('%s/figure_6.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')
fig2.savefig('%s/figure_6.svg' % (parse_midas_data.analysis_directory),bbox_inches='tight')
fig.savefig('%s/supplemental_within_across_species.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')
sys.stderr.write("Done!\n")


fig3.savefig('%s/within_allele_freqs.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')

from scipy.stats import fisher_exact

for i in xrange(1,len(total_freq_all_snps)/2+1):
    
    non_counts = total_freq_snps['1D'][:i].sum()+total_freq_snps['1D'][-i:].sum()
    syn_counts = total_freq_snps['4D'][:i].sum()+total_freq_snps['4D'][-i:].sum()
    
    
    print derived_freq_bins[i], non_counts, syn_counts, non_counts*1.0/syn_counts/(total_random_null_snps['1D']*1.0/total_random_null_snps['4D']),  fisher_exact([[non_counts, total_random_null_snps['1D']], [syn_counts, total_random_null_snps['4D']]])


       