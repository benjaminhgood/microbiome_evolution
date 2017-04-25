import matplotlib  
matplotlib.use('Agg') 
import parse_midas_data
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
parser.add_argument("--chunk-size", type=int, help="max number of records to load", default=1000000000)
args = parser.parse_args()

species_name = args.species_name
debug = args.debug
chunk_size = args.chunk_size
################################################################################

min_change = 0.8
min_coverage = 20
alpha = 0.5 # Confidence interval range for rate estimates
allowed_variant_types = set(['1D','2D','3D','4D'])
max_clade_d = 1e-02

# Load subject and sample metadata
sys.stderr.write("Loading HMP metadata...\n")
subject_sample_map = parse_midas_data.parse_subject_sample_map()
sys.stderr.write("Done!\n")

# Load core gene set
sys.stderr.write("Loading core genes...\n")
core_genes = parse_midas_data.load_core_genes(species_name)
sys.stderr.write("Done! Core genome consists of %d genes\n" % len(core_genes))
 
    
# Load genomic coverage distributions
sample_coverage_histograms, samples = parse_midas_data.parse_coverage_distribution(species_name)
median_coverages = numpy.array([stats_utils.calculate_nonzero_median_from_histogram(sample_coverage_histogram) for sample_coverage_histogram in sample_coverage_histograms])
sample_coverage_map = {samples[i]: median_coverages[i] for i in xrange(0,len(samples))}

# Load pi information for species_name
sys.stderr.write("Loading within-sample diversity for %s...\n" % species_name)
samples, total_pis, total_pi_opportunities = parse_midas_data.parse_within_sample_pi(species_name, allowed_genes=core_genes, debug=debug)
sys.stderr.write("Done!\n")
pis = total_pis/total_pi_opportunities

median_coverages = numpy.array([sample_coverage_map[samples[i]] for i in xrange(0,len(samples))])

# Only plot samples above a certain depth threshold that are "haploids"
snp_samples = samples[(median_coverages>=min_coverage)*(pis<=1e-03)]
# Restrict to single timepoint single timepoints per person
unique_subject_idxs = parse_midas_data.calculate_unique_samples(subject_sample_map, snp_samples)
snp_samples = snp_samples[unique_subject_idxs]



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


# calculate compressed distance matrix suitable for agglomerative clustering
#Y = []
#for i in xrange(0,substitution_rate.shape[0]):
#    for j in xrange(i+1,substitution_rate.shape[1]):
#        Y.append(substitution_rate[i,j]) 
#Y = numpy.array(Y) 
    
#Z = linkage(Y, method='average')        
    
#c, coph_dists = cophenet(Z, Y)
#sys.stderr.write("cophenetic correlation: %g\n" % c)


# Make a zoomed dendrogram
#pylab.figure(2, figsize=(15, 5))
#pylab.title('Zoomed hierarchical clustering dendrogram for %s' % species_name)
#pylab.xlabel('sample index')
#pylab.ylabel('distance')
#dendrogram(
#    Z,
#    leaf_rotation=90.,  # rotates the x axis labels
#    leaf_font_size=8.,  # font size for the x axis labels
#)
#pylab.ylim([0,5e-04])
#pylab.savefig('%s/%s_zoomed_dendrogram.pdf' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight')
#pylab.savefig('%s/%s_zoomed_dendrogram.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight',dpi=300)

# only focus on the members of the largest clade
remapped_cluster_idxss = [cluster_idxs[coarse_grained_idxs] for cluster_idxs in coarse_grained_cluster_idxss]
largest_clade_idxs = remapped_cluster_idxss[0]
largest_clade_size = largest_clade_idxs.sum()
coarse_grained_samples = snp_samples[coarse_grained_idxs]

if len(coarse_grained_samples)>2:

    sys.stderr.write("Continuing with %d samples...\n" % len(coarse_grained_samples))

    # Load SNP information for species_name
    sys.stderr.write("Re-loading %s...\n" % species_name)

    snp_difference_matrix = numpy.array([])
    snp_opportunity_matrix = numpy.array([])
    
    polymorphic_freqs = [] 
    inconsistent_freqs = [] 
    null_inconsistent_freqs = []

    # initialize distance bins for LD computations
    distance_bins = numpy.logspace(0,4,20) # bins start from 1 to 10^4 and there are 20 evenly spaced bins log(1)=0, log(10^4)-4
    distance_bin_locations = numpy.array(distance_bins[:-1],copy=True) # shifted one to avoid edge effects for plotting.
    distance_bins[0] = 0.5 # made smallest bin 0.5 to avoid edge effects
    distance_bins[-1] = 1e09 # made largest bin very large to catch anything >10^4. 
    
    binned_rsquared_numerators = numpy.zeros_like(distance_bin_locations)
    binned_rsquared_denominators = numpy.zeros_like(distance_bin_locations)   
    total_control_rsquared_numerators = 0
    total_control_rsquared_denominators = 0
    
    
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
        dummy_samples, allele_counts_map, passed_sites_map, final_line_number = parse_midas_data.parse_snps(species_name, debug=debug, allowed_variant_types=allowed_variant_types, allowed_samples=coarse_grained_samples,chunk_size=chunk_size,initial_line_number=final_line_number, allowed_genes=core_genes)
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
            
    
        #sys.stderr.write("Calculating something else?...\n")
        
        
    
    substitution_rate = snp_difference_matrix*1.0/snp_opportunity_matrix 
    
    polymorphic_freqs = numpy.array(polymorphic_freqs)
    inconsistent_freqs = numpy.array(inconsistent_freqs)
    null_inconsistent_freqs = numpy.array(null_inconsistent_freqs)
    
    all_binned_rsquareds = all_binned_rsquared_numerators/(all_binned_rsquared_denominators+(all_binned_rsquared_denominators==0))    
    all_control_rsquareds = all_total_control_rsquared_numerators/(all_total_control_rsquared_denominators+(all_total_control_rsquared_denominators==0))
    
    print binned_rsquared_numerators
    print binned_rsquared_denominators
    
    binned_rsquareds = binned_rsquared_numerators/(binned_rsquared_denominators+(binned_rsquared_denominators==0))    
    control_rsquareds = total_control_rsquared_numerators/(total_control_rsquared_denominators+(total_control_rsquared_denominators==0))

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
    pylab.ylabel("Ohta and Kimura's $\\sigma^2_d$")

    pylab.xlim([1,1e04])
    pylab.ylim([1e-02,1])
    pylab.loglog(distance_bin_locations[all_binned_rsquared_denominators>0], all_binned_rsquareds[all_binned_rsquared_denominators>0],'.-',color='0.7',label='All')
    pylab.loglog(distance_bin_locations[binned_rsquared_denominators>0], binned_rsquareds[binned_rsquared_denominators>0],'b.-',label='Largest clade')
    pylab.loglog(distance_bin_locations, numpy.ones_like(distance_bin_locations)*control_rsquareds,'b:')
    pylab.legend(loc='upper right',frameon=False,fontsize=6)
    pylab.savefig('%s/%s_intragene_ld.pdf' % (parse_midas_data.analysis_directory, species_name), bbox_inches='tight')
    pylab.savefig('%s/%s_intragene_ld.png' % (parse_midas_data.analysis_directory, species_name), bbox_inches='tight', dpi=300)


else:
    
    sys.stderr.write("No clades!\n")
    
