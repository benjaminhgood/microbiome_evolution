import matplotlib  
matplotlib.use('Agg') 
import parse_midas_data
import parse_HMP_data
import pylab
import sys
import numpy

import diversity_utils
import gene_diversity_utils

import stats_utils
import matplotlib.colors as colors
import matplotlib.cm as cmx
from math import log10,ceil
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy.random import randint, choice

from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster.hierarchy import cophenet
from scipy.cluster.hierarchy import fcluster

mpl.rcParams['font.size'] = 7
mpl.rcParams['lines.linewidth'] = 0.5
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'


################################################################################
#
# Standard header to read in argument information
#
################################################################################
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("species_name", help="name of species to process")
parser.add_argument("--debug", help="Loads only a subset of SNPs for speed", action="store_true")
parser.add_argument("--include-china", help="Includes Chinese subjects from Qin et al (Nature, 2012)", action="store_true")
parser.add_argument("--chunk-size", type=int, help="max number of records to load", default=1000000000)
args = parser.parse_args()

species_name = args.species_name
debug = args.debug
chunk_size = args.chunk_size
include_china = args.include_china
################################################################################

min_change = 0.8
min_coverage = 20
alpha = 0.5 # Confidence interval range for rate estimates
allowed_variant_types = set(['1D','2D','3D','4D'])
max_clade_d = 1e-02

# Load subject and sample metadata
sys.stderr.write("Loading HMP metadata...\n")
subject_sample_map = parse_HMP_data.parse_subject_sample_map()
sample_country_map = parse_HMP_data.parse_sample_country_map()
sys.stderr.write("Done!\n")

# Load core gene set
sys.stderr.write("Loading core genes...\n")
core_genes = parse_midas_data.load_core_genes(species_name)
sys.stderr.write("Done! Core genome consists of %d genes\n" % len(core_genes))
 
 # Only plot samples above a certain depth threshold that are "haploids"
snp_samples = diversity_utils.calculate_haploid_samples(species_name, debug=debug)
snp_samples = snp_samples[ diversity_utils.parse_midas_data.calculate_unique_samples(subject_sample_map, snp_samples)]

# Analyze SNPs, looping over chunk sizes. 
# Clunky, but necessary to limit memory usage on cluster

# Load SNP information for species_name
sys.stderr.write("Loading %d samples for %s...\n" % (len(snp_samples), species_name))

snp_difference_matrix = numpy.array([])
snp_opportunity_matrix = numpy.array([])

final_line_number = 0
while final_line_number >= 0:
    
    sys.stderr.write("Loading chunk starting @ %d...\n" % final_line_number)
    dummy_samples, allele_counts_map, passed_sites_map, final_line_number = parse_midas_data.parse_snps(species_name, debug=debug, allowed_variant_types=allowed_variant_types, allowed_genes=core_genes, allowed_samples=snp_samples,chunk_size=chunk_size,initial_line_number=final_line_number)
    sys.stderr.write("Done! Loaded %d genes\n" % len(allele_counts_map.keys()))
    
    # Calculate fixation matrix
    sys.stderr.write("Calculating matrix of snp differences...\n")
    chunk_snp_difference_matrix, chunk_snp_opportunity_matrix = diversity_utils.calculate_fixation_matrix(allele_counts_map, passed_sites_map, allowed_variant_types=allowed_variant_types, allowed_genes=core_genes, min_change=min_change)    
    sys.stderr.write("Done!\n")
    
    if snp_difference_matrix.shape[0]==0:
        snp_difference_matrix = numpy.zeros_like(chunk_snp_difference_matrix)*1.0
        snp_opportunity_matrix = numpy.zeros_like(snp_difference_matrix)*1.0
    
    snp_difference_matrix += chunk_snp_difference_matrix
    snp_opportunity_matrix += chunk_snp_opportunity_matrix
    
     
substitution_rate = snp_difference_matrix*1.0/snp_opportunity_matrix 

coarse_grained_idxs, coarse_grained_cluster_list = diversity_utils.cluster_samples(substitution_rate, min_d=5e-04, max_ds=[max_clade_d])
coarse_grained_cluster_idxss = coarse_grained_cluster_list[0] 
coarse_grained_cluster_sizes = numpy.array([cluster_idxs.sum() for cluster_idxs in coarse_grained_cluster_idxss])

print "Top level:", len(coarse_grained_cluster_idxss), coarse_grained_cluster_sizes



# only focus on the members of the largest clade
remapped_cluster_idxss = [cluster_idxs[coarse_grained_idxs] for cluster_idxs in coarse_grained_cluster_idxss]
largest_clade_idxs = remapped_cluster_idxss[0]
largest_clade_size = largest_clade_idxs.sum()

coarse_grained_samples = snp_samples[coarse_grained_idxs]
largest_clade_samples = set(coarse_grained_samples[largest_clade_idxs])

# Load gene coverage information for species_name
sys.stderr.write("Loading pangenome data for %s...\n" % species_name)
gene_samples, gene_names, gene_presence_matrix, gene_depth_matrix, marker_coverages, gene_reads_matrix = parse_midas_data.parse_pangenome_data(species_name,allowed_samples=snp_samples)
sys.stderr.write("Done!\n")

desired_gene_sample_idxs = numpy.array([sample in largest_clade_samples for sample in gene_samples])
gene_prevalences = gene_diversity_utils.calculate_fractional_gene_prevalences(gene_depth_matrix[:,desired_gene_sample_idxs], marker_coverages[desired_gene_sample_idxs])

gene_prevalence_map = {gene_name: prevalence for gene_name, prevalence in zip(gene_names, gene_prevalences)} 

if len(coarse_grained_samples)>2:

    sys.stderr.write("Continuing with %d samples...\n" % len(coarse_grained_samples))

    # Load SNP information for species_name
    sys.stderr.write("Re-loading %s...\n" % species_name)

    snp_difference_matrix = numpy.array([])
    snp_opportunity_matrix = numpy.array([])
    
    polymorphic_freqs = [] 
    inconsistent_freqs = [] 
    null_inconsistent_freqs = []

    # initialize prevalence bins 
    prevalence_bins = numpy.linspace(0,1,11)
    prevalence_bins[-1] = 1.01
    
    prevalence_locations = numpy.arange(0,10)*0.1+0.05
    
    prevalence_synonymous_differences = { i: 0.0 for i in xrange(0,len(prevalence_locations)) }
    prevalence_synonymous_opportunities = { i: 0.0 for i in xrange(0,len(prevalence_locations)) }
    prevalence_nonsynonymous_differences = { i: 0.0 for i in xrange(0,len(prevalence_locations)) }
    prevalence_nonsynonymous_opportunities = { i: 0.0 for i in xrange(0,len(prevalence_locations)) }


    # initialize distance bins for LD computations
    distance_bins = numpy.logspace(0,4,20) # bins start from 1 to 10^4 and there are 20 evenly spaced bins log(1)=0, log(10^4)-4
    distance_bin_locations = numpy.array(distance_bins[:-1],copy=True) # shifted one to avoid edge effects for plotting.
    distance_bins[0] = 0.5 # made smallest bin 0.5 to avoid edge effects
    distance_bins[-1] = 1e09 # made largest bin very large to catch anything >10^4. 
    
    binned_rsquared_numerators = numpy.zeros_like(distance_bin_locations)
    binned_rsquared_denominators = numpy.zeros_like(distance_bin_locations)   
    total_control_rsquared_numerators = 0
    total_control_rsquared_denominators = 0
    
    nonsynonymous_binned_rsquared_numerators = numpy.zeros_like(distance_bin_locations)
    nonsynonymous_binned_rsquared_denominators = numpy.zeros_like(distance_bin_locations)   
    nonsynonymous_total_control_rsquared_numerators = 0
    nonsynonymous_total_control_rsquared_denominators = 0
    
    all_binned_rsquared_numerators = numpy.zeros_like(distance_bin_locations)
    all_binned_rsquared_denominators = numpy.zeros_like(distance_bin_locations)   
    all_total_control_rsquared_numerators = 0
    all_total_control_rsquared_denominators = 0
    
    maf_bins = numpy.arange(1,largest_clade_size+1)*1.0/largest_clade_size
    maf_bins -= (maf_bins[1]-maf_bins[0])/2
    maf_bins[0]=-0.1
    maf_bins[-1] = 1.1
    mafs = numpy.arange(1,largest_clade_size)*1.0/largest_clade_size
    
    synonymous_sfs = numpy.zeros_like(mafs)
    nonsynonymous_sfs = numpy.zeros_like(mafs)
    
    final_line_number = 0
    while final_line_number >= 0:
    
        sys.stderr.write("Loading chunk starting @ %d...\n" % final_line_number)
        dummy_samples, allele_counts_map, passed_sites_map, final_line_number = parse_midas_data.parse_snps(species_name, debug=debug, allowed_variant_types=allowed_variant_types, allowed_samples=coarse_grained_samples,chunk_size=chunk_size,initial_line_number=final_line_number)
        sys.stderr.write("Done! Loaded %d genes\n" % len(allele_counts_map.keys()))
    
        # Calculate fixation matrix
        sys.stderr.write("Calculating matrix of snp differences...\n")
        chunk_snp_difference_matrix, chunk_snp_opportunity_matrix = diversity_utils.calculate_fixation_matrix(allele_counts_map, passed_sites_map, allowed_genes=core_genes, min_change=min_change)    
        sys.stderr.write("Done!\n")
    
        if snp_difference_matrix.shape[0]==0:
            snp_difference_matrix = numpy.zeros_like(chunk_snp_difference_matrix)*1.0
            snp_opportunity_matrix = numpy.zeros_like(snp_difference_matrix)*1.0
    
        snp_difference_matrix += chunk_snp_difference_matrix
        snp_opportunity_matrix += chunk_snp_opportunity_matrix
    
        sys.stderr.write("Calculating phylogenetic consistency...\n")
        chunk_polymorphic_freqs, chunk_inconsistent_freqs, chunk_null_inconsistent_freqs = diversity_utils.calculate_phylogenetic_consistency(allele_counts_map, passed_sites_map, [largest_clade_idxs], allowed_genes=core_genes)    
        polymorphic_freqs.extend(chunk_polymorphic_freqs) 
        inconsistent_freqs.extend(chunk_inconsistent_freqs)
        null_inconsistent_freqs.extend(chunk_null_inconsistent_freqs)
        sys.stderr.write("Done!\n")
    
        sys.stderr.write("Calculating the SFS...\n")
        chunk_synonymous_freqs = diversity_utils.calculate_pooled_freqs(allele_counts_map, passed_sites_map, allowed_sample_idxs=largest_clade_idxs, allowed_variant_types = set(['4D']), allowed_genes=core_genes)
        chunk_nonsynonymous_freqs = diversity_utils.calculate_pooled_freqs(allele_counts_map, passed_sites_map, allowed_sample_idxs=largest_clade_idxs, allowed_variant_types = set(['1D']), allowed_genes=core_genes)
        
        chunk_synonymous_sfs, dummy = numpy.histogram(chunk_synonymous_freqs, bins=maf_bins) 
        synonymous_sfs += chunk_synonymous_sfs
        
        chunk_nonsynonymous_sfs, dummy = numpy.histogram(chunk_nonsynonymous_freqs, bins=maf_bins) 
        nonsynonymous_sfs += chunk_nonsynonymous_sfs
        
        sys.stderr.write("Calculating intra-gene synonymous LD...\n")
        for gene_name in allele_counts_map.keys():
            
            if gene_name not in core_genes:
                continue
        
            locations = numpy.array([location for chromosome, location in allele_counts_map[gene_name]['4D']['locations']])*1.0
            allele_counts = allele_counts_map[gene_name]['4D']['alleles']
        
            if len(allele_counts)==0:
                # no diversity to look at!
                continue
        
            # pick a random gene somewhere else as a control
            control_gene_name = gene_name
            control_allele_counts = []
            while gene_name==control_gene_name or len(control_allele_counts)==0:
                control_gene_name = choice(allele_counts_map.keys())
                control_allele_counts = allele_counts_map[control_gene_name]['4D']['alleles']
        

            #compute the distances between all pairs of sites 
            # None in the two index positions results in a transpose of the vector relative to each other
            # Subtraction between the two vectors results in pairwise subtraction of each element in each vector.
            distances = numpy.fabs(locations[:,None]-locations[None,:])
    
            rsquared_numerators, rsquared_denominators = diversity_utils.calculate_unbiased_sigmasquared(allele_counts, allele_counts)
            control_rsquared_numerators, control_rsquared_denominators = diversity_utils.calculate_unbiased_sigmasquared(allele_counts, control_allele_counts)
        
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
            distances = distances[rsquared_denominators>0]
            rsquared_numerators = rsquared_numerators[rsquared_denominators>0] 
            rsquared_denominators = rsquared_denominators[rsquared_denominators>0]
        
            if len(distances) == 0:
                continue

            # numpy.digitize: For each distance value, return the bin index it belongs to in distances_bins. 
            bin_idxs = numpy.digitize(distances,bins=distance_bins)-1
            
            for i in xrange(0,len(bin_idxs)):
        
                all_binned_rsquared_numerators[bin_idxs[i]] += rsquared_numerators[i]
            
                all_binned_rsquared_denominators[bin_idxs[i]] += rsquared_denominators[i]
        
            control_rsquared_numerators = control_rsquared_numerators[control_rsquared_denominators>0]
            control_rsquared_denominators = control_rsquared_denominators[control_rsquared_denominators>0]
        
            all_total_control_rsquared_numerators += (control_rsquared_numerators).sum()
            all_total_control_rsquared_denominators += (control_rsquared_denominators).sum()
           
        
            # Now restrict to largest clade
            allele_counts = allele_counts[:,largest_clade_idxs,:]
            control_allele_counts = control_allele_counts[:,largest_clade_idxs,:]
    
            #compute the distances between all pairs of sites 
            # None in the two index positions results in a transpose of the vector relative to each other
            # Subtraction between the two vectors results in pairwise subtraction of each element in each vector.
            distances = numpy.fabs(locations[:,None]-locations[None,:])
    
            rsquared_numerators, rsquared_denominators = diversity_utils.calculate_unbiased_sigmasquared(allele_counts, allele_counts)
            control_rsquared_numerators, control_rsquared_denominators = diversity_utils.calculate_unbiased_sigmasquared(allele_counts, control_allele_counts)
        
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
            distances = distances[rsquared_denominators>0]
            rsquared_numerators = rsquared_numerators[rsquared_denominators>0] 
            rsquared_denominators = rsquared_denominators[rsquared_denominators>0]
        
            if len(distances) == 0:
                continue

            # numpy.digitize: For each distance value, return the bin index it belongs to in distances_bins. 
            bin_idxs = numpy.digitize(distances,bins=distance_bins)-1
            
            for i in xrange(0,len(bin_idxs)):
        
                binned_rsquared_numerators[bin_idxs[i]] += rsquared_numerators[i]
            
                binned_rsquared_denominators[bin_idxs[i]] += rsquared_denominators[i]
        
            control_rsquared_numerators = control_rsquared_numerators[control_rsquared_denominators>0]
            control_rsquared_denominators = control_rsquared_denominators[control_rsquared_denominators>0]
        
            total_control_rsquared_numerators += (control_rsquared_numerators).sum()
            total_control_rsquared_denominators += (control_rsquared_denominators).sum()
            
        # Now repeat for nonsynonymous ones
        for gene_name in allele_counts_map.keys():
            
            if gene_name not in core_genes:
                continue
            
            locations = numpy.array([location for chromosome, location in allele_counts_map[gene_name]['1D']['locations']])*1.0
            allele_counts = allele_counts_map[gene_name]['1D']['alleles']
        
            if len(allele_counts)==0:
                # no diversity to look at!
                continue
            
            # pick a random gene somewhere else as a control
            control_gene_name = gene_name
            control_allele_counts = []
            while gene_name==control_gene_name or len(control_allele_counts)==0:
                control_gene_name = choice(allele_counts_map.keys())
                control_allele_counts = allele_counts_map[control_gene_name]['1D']['alleles']
         
        
            # Now restrict to largest clade
            allele_counts = allele_counts[:,largest_clade_idxs,:]
            control_allele_counts = control_allele_counts[:,largest_clade_idxs,:]
    
            #compute the distances between all pairs of sites 
            # None in the two index positions results in a transpose of the vector relative to each other
            # Subtraction between the two vectors results in pairwise subtraction of each element in each vector.
            distances = numpy.fabs(locations[:,None]-locations[None,:])
    
            rsquared_numerators, rsquared_denominators = diversity_utils.calculate_unbiased_sigmasquared(allele_counts, allele_counts)
            control_rsquared_numerators, control_rsquared_denominators = diversity_utils.calculate_unbiased_sigmasquared(allele_counts, control_allele_counts)
        
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
            distances = distances[rsquared_denominators>0]
            rsquared_numerators = rsquared_numerators[rsquared_denominators>0] 
            rsquared_denominators = rsquared_denominators[rsquared_denominators>0]
        
            if len(distances) == 0:
                continue

            # numpy.digitize: For each distance value, return the bin index it belongs to in distances_bins. 
            bin_idxs = numpy.digitize(distances,bins=distance_bins)-1
            
            for i in xrange(0,len(bin_idxs)):
        
                nonsynonymous_binned_rsquared_numerators[bin_idxs[i]] += rsquared_numerators[i]
            
                nonsynonymous_binned_rsquared_denominators[bin_idxs[i]] += rsquared_denominators[i]
        
            control_rsquared_numerators = control_rsquared_numerators[control_rsquared_denominators>0]
            control_rsquared_denominators = control_rsquared_denominators[control_rsquared_denominators>0]
        
            nonsynonymous_total_control_rsquared_numerators += (control_rsquared_numerators).sum()
            nonsynonymous_total_control_rsquared_denominators += (control_rsquared_denominators).sum()    
        #sys.stderr.write("Calculating something else?...\n")
        
        # Construct indices of pairs in largest subclade
        largest_clade_idx_numbers = numpy.nonzero(largest_clade_idxs)[0]
        pair_idxs_1 = []
        pair_idxs_2 = []
        for i in xrange(0,len(largest_clade_idx_numbers)):
            for j in xrange(i+1,len(largest_clade_idx_numbers)):
                pair_idxs_1.append(largest_clade_idx_numbers[i])
                pair_idxs_2.append(largest_clade_idx_numbers[j])
        
        pair_idxs = (numpy.array(pair_idxs_1),numpy.array(pair_idxs_2))
        
        # Now repeat for prevalence-specific dNdSs
        for gene_name in allele_counts_map.keys():
        
            if gene_name in gene_prevalence_map:
                prevalence = gene_prevalence_map[gene_name]
                prevalence_idx = numpy.digitize([prevalence], prevalence_bins)[0]-1
            else:
                sys.stderr.write("No prevalence found: %s!\n" % gene_name)
                prevalence_idx = 0
                
            #print prevalence, prevalence_idx
            #print gene_snp_difference_matrix[pair_idxs].sum()
            # Calculate fixation matrix
            gene_snp_difference_matrix, gene_snp_opportunity_matrix = diversity_utils.calculate_fixation_matrix(allele_counts_map, passed_sites_map, allowed_variant_types=set(['4D']), allowed_genes=set([gene_name]), min_change=min_change)     
            prevalence_synonymous_differences[prevalence_idx] += gene_snp_difference_matrix[pair_idxs].sum()
            prevalence_synonymous_opportunities[prevalence_idx] += gene_snp_opportunity_matrix[pair_idxs].sum()
            
            # Calculate fixation matrix
            gene_snp_difference_matrix, gene_snp_opportunity_matrix = diversity_utils.calculate_fixation_matrix(allele_counts_map, passed_sites_map, allowed_variant_types=set(['1D']), allowed_genes=set([gene_name]), min_change=min_change)     
            prevalence_nonsynonymous_differences[prevalence_idx] += gene_snp_difference_matrix[pair_idxs].sum()
            prevalence_nonsynonymous_opportunities[prevalence_idx] += gene_snp_opportunity_matrix[pair_idxs].sum()
            
         
        
    
    substitution_rate = snp_difference_matrix*1.0/snp_opportunity_matrix 
    
    polymorphic_freqs = numpy.array(polymorphic_freqs)
    inconsistent_freqs = numpy.array(inconsistent_freqs)
    null_inconsistent_freqs = numpy.array(null_inconsistent_freqs)
    
    all_binned_rsquareds = all_binned_rsquared_numerators/(all_binned_rsquared_denominators+(all_binned_rsquared_denominators==0))    
    all_control_rsquareds = all_total_control_rsquared_numerators/(all_total_control_rsquared_denominators+(all_total_control_rsquared_denominators==0))
    
    binned_rsquareds = binned_rsquared_numerators/(binned_rsquared_denominators+(binned_rsquared_denominators==0))    
    control_rsquareds = total_control_rsquared_numerators/(total_control_rsquared_denominators+(total_control_rsquared_denominators==0))

    nonsynonymous_binned_rsquareds = nonsynonymous_binned_rsquared_numerators/(nonsynonymous_binned_rsquared_denominators+(nonsynonymous_binned_rsquared_denominators==0))    
    nonsynonymous_control_rsquareds = nonsynonymous_total_control_rsquared_numerators/(nonsynonymous_total_control_rsquared_denominators+(nonsynonymous_total_control_rsquared_denominators==0))

    pylab.figure(1,figsize=(3.42,2))
    pylab.suptitle(species_name)
    pylab.xlabel('Within-clade MAF, $f$')
    pylab.ylabel('Scaled SFS, $f(1-f) P(f)$')
    
    pylab.plot(mafs, synonymous_sfs*mafs*(1-mafs)/(synonymous_sfs*mafs*(1-mafs)).sum(), 'b.-',label='4D')
    pylab.plot(mafs, nonsynonymous_sfs*mafs*(1-mafs)/(nonsynonymous_sfs*mafs*(1-mafs)).sum(),'r.-',label='1D')
    
    pylab.xlim([0,0.5])
    pylab.legend(loc='upper right',frameon=False,fontsize=6)
    pylab.savefig('%s/%s_pooled_sfs.pdf' % (parse_midas_data.analysis_directory, species_name), bbox_inches='tight')
    pylab.savefig('%s/%s_pooled_sfs.png' % (parse_midas_data.analysis_directory, species_name), bbox_inches='tight', dpi=300)
 
    pylab.figure(2,figsize=(3.42,2))
    pylab.suptitle(species_name)

    #xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(polymorphic_freqs)
    #pylab.step(xs,ns*1.0/ns[0],'b-',label='All polymorphisms')
    
    if len(null_inconsistent_freqs)>0:
        xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(null_inconsistent_freqs)
        pylab.step(xs,ns*1.0/ns[0],'-',color='0.7',linewidth=0.5, label=('Unlinked expectation'))
     
    
    if len(inconsistent_freqs)>0:
        xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(inconsistent_freqs)
        pylab.step(xs,ns*1.0/ns[0],'r-',label=('Inconsistent ($d=%g$)' % max_clade_d))
       
    pylab.xlim([0,0.5])
    pylab.ylim([0,1])
    pylab.xlabel('Within-clade MAF, $f$')
    pylab.ylabel('SNPs $\geq f$')
    pylab.legend(loc='upper right', frameon=False,fontsize=6)
    
    pylab.savefig('%s/%s_phylogenetically_inconsistent_sfs.pdf' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight')
    pylab.savefig('%s/%s_phylogenetically_inconsistent_sfs.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight',dpi=300)

    
    pylab.figure(3,figsize=(3.42,2))
    pylab.suptitle(species_name)
    pylab.xlabel('Distance between SNPs')
    pylab.ylabel("Linkage disequilibrium, $\\sigma^2_d$")

    pylab.gca().spines['top'].set_visible(False)
    pylab.gca().spines['right'].set_visible(False)
    pylab.gca().get_xaxis().tick_bottom()
    pylab.gca().get_yaxis().tick_left()


    control_x = 2e04

    pylab.xlim([1,3e04])
    pylab.ylim([1e-02,1])
    pylab.loglog(distance_bin_locations[all_binned_rsquared_denominators>0], all_binned_rsquareds[all_binned_rsquared_denominators>0],'b.-',alpha=0.5,label='All (4D)')
    pylab.loglog([distance_bin_locations[all_binned_rsquared_denominators>0][-1], control_x], [all_binned_rsquareds[all_binned_rsquared_denominators>0][-1], all_control_rsquareds], 'b:',alpha=0.5)
    pylab.loglog([control_x], [all_control_rsquareds], 'b.',alpha=0.5)
    
    
    pylab.loglog(distance_bin_locations[binned_rsquared_denominators>0], binned_rsquareds[binned_rsquared_denominators>0],'b.-',label='Largest clade (4D)')
    pylab.loglog([distance_bin_locations[binned_rsquared_denominators>0][-1], control_x], [binned_rsquareds[binned_rsquared_denominators>0][-1], control_rsquareds], 'b:')
    pylab.loglog([control_x], [control_rsquareds], 'b.')
    
    pylab.loglog(distance_bin_locations[nonsynonymous_binned_rsquared_denominators>0], nonsynonymous_binned_rsquareds[nonsynonymous_binned_rsquared_denominators>0],'r.-',label='Nonsynonymous (1D)')
    pylab.loglog([distance_bin_locations[nonsynonymous_binned_rsquared_denominators>0][-1], control_x], [nonsynonymous_binned_rsquareds[nonsynonymous_binned_rsquared_denominators>0][-1], nonsynonymous_control_rsquareds], 'r:')
    pylab.loglog([control_x], [nonsynonymous_control_rsquareds], 'r.')
    
    
    pylab.legend(loc='lower left',frameon=False,fontsize=6)
    pylab.savefig('%s/%s_intragene_ld.pdf' % (parse_midas_data.analysis_directory, species_name), bbox_inches='tight')
    pylab.savefig('%s/%s_intragene_ld.png' % (parse_midas_data.analysis_directory, species_name), bbox_inches='tight', dpi=300)

    pylab.figure(4,figsize=(3.42,4))
    pylab.suptitle(species_name)
    
    outer_grid = gridspec.GridSpec(2, 1, height_ratios=[1,1], hspace=0.05)
    
    fig = pylab.gcf()

    prevalence_dS_axis = plt.Subplot(fig, outer_grid[0])
    fig.add_subplot(prevalence_dS_axis)

    prevalence_dNdS_axis = plt.Subplot(fig, outer_grid[1])
    fig.add_subplot(prevalence_dNdS_axis)

    prevalence_dNdS_axis.spines['top'].set_visible(False)
    prevalence_dNdS_axis.spines['right'].set_visible(False)
    prevalence_dNdS_axis.get_xaxis().tick_bottom()
    prevalence_dNdS_axis.get_yaxis().tick_left()

    prevalence_dS_axis.spines['top'].set_visible(False)
    prevalence_dS_axis.spines['right'].set_visible(False)
    prevalence_dS_axis.get_xaxis().tick_bottom()
    prevalence_dS_axis.get_yaxis().tick_left()


    prevalence_dNdS_axis.set_xlim([0,1])
    prevalence_dNdS_axis.set_ylim([0,1.1])
    
    prevalence_dNdS_axis.plot([0,1],[1,1],'k:')
    
    prevalence_dS_axis.set_xlim([0,1])
     
    prevalence_dNdS_axis.set_xlabel('Gene prevalence (%)')
    prevalence_dNdS_axis.set_ylabel("dN/dS")
    prevalence_dS_axis.set_ylabel("dS")
    prevalence_dS_axis.set_xticklabels([])
    
    prevalence_dNdSs = []
    prevalence_dSs = []
    for prevalence_idx in xrange(0,len(prevalence_locations)):
        
        
        dNdS = ( (prevalence_nonsynonymous_differences[prevalence_idx]+1.0)/(prevalence_nonsynonymous_opportunities[prevalence_idx]+1.0) ) / ( (prevalence_synonymous_differences[prevalence_idx]+1.0)/(prevalence_synonymous_opportunities[prevalence_idx]+1.0) )
        prevalence_dNdSs.append(dNdS)
    
        prevalence_dSs.append( (prevalence_synonymous_differences[prevalence_idx])*1.0/(prevalence_synonymous_opportunities[prevalence_idx]+(prevalence_synonymous_opportunities[prevalence_idx]<1)) )
    
    prevalence_dS_axis.semilogy(prevalence_locations, prevalence_dSs, 'b.-')
    prevalence_dNdS_axis.plot(prevalence_locations, prevalence_dNdSs, 'b.-')
    
    pylab.savefig('%s/%s_diversity_vs_prevalence.pdf' % (parse_midas_data.analysis_directory, species_name), bbox_inches='tight')
    pylab.savefig('%s/%s_diversity_vs_prevalence.png' % (parse_midas_data.analysis_directory, species_name), bbox_inches='tight', dpi=300)

else:
    
    sys.stderr.write("No clades!\n")
    
