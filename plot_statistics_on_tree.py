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
from numpy.random import randint

from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster.hierarchy import cophenet
from scipy.cluster.hierarchy import fcluster

mpl.rcParams['font.size'] = 8
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


max_ds = numpy.logspace(-4,-1,19)
fine_grained_idxs, fine_grained_cluster_list = diversity_utils.cluster_samples(substitution_rate, min_d=1e-07, max_ds=max_ds)


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

# Plot phylogenetic consistency vs time
remapped_cluster_list = []
for cluster_idxss in fine_grained_cluster_list:

    remapped_cluster_idxss = [cluster_idxs[fine_grained_idxs] for cluster_idxs in cluster_idxss]
    remapped_cluster_list.append(remapped_cluster_idxss)

fine_grained_samples = snp_samples[fine_grained_idxs]

if len(fine_grained_samples)>2:

    sys.stderr.write("Continuing with %d samples...\n" % len(fine_grained_samples))

    # Load SNP information for species_name
    sys.stderr.write("Re-loading %s...\n" % species_name)

    snp_difference_matrix = numpy.array([])
    snp_opportunity_matrix = numpy.array([])

    total_polymorphic_sites = [0 for max_d in max_ds]
    total_inconsistent_sites = [0 for max_d in max_ds]
    total_null_inconsistent_sites = [0 for max_d in max_ds]
 
    final_line_number = 0
    while final_line_number >= 0:
    
        sys.stderr.write("Loading chunk starting @ %d...\n" % final_line_number)
        dummy_samples, allele_counts_map, passed_sites_map, final_line_number = parse_midas_data.parse_snps(species_name, debug=debug, allowed_variant_types=allowed_variant_types, allowed_samples=fine_grained_samples,chunk_size=chunk_size,initial_line_number=final_line_number, allowed_genes=core_genes)
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
        for i in xrange(0,len(max_ds)):
        
            cluster_idxss = remapped_cluster_list[i]
        
            chunk_polymorphic_freqs, chunk_inconsistent_freqs, chunk_null_inconsistent_freqs = diversity_utils.calculate_phylogenetic_consistency(allele_counts_map, passed_sites_map, cluster_idxss, allowed_genes=core_genes)    
        
            total_polymorphic_sites[i] += len(chunk_polymorphic_freqs) 
            total_inconsistent_sites[i] += len(chunk_inconsistent_freqs)
            total_null_inconsistent_sites[i] += len(chunk_null_inconsistent_freqs)
    
        sys.stderr.write("Done!\n")
    
    substitution_rate = snp_difference_matrix*1.0/snp_opportunity_matrix 
    
        
    total_polymorphic_sites = numpy.array(total_polymorphic_sites)
    total_inconsistent_sites = numpy.array(total_inconsistent_sites)
    total_null_inconsistent_sites = numpy.array(total_null_inconsistent_sites)
  
    fraction_inconsistent = total_inconsistent_sites*1.0/total_polymorphic_sites
    null_fraction_inconsistent = total_null_inconsistent_sites*1.0/total_polymorphic_sites
    
    print total_polymorphic_sites
    print total_inconsistent_sites
    print total_null_inconsistent_sites
    
    # Set up figure
    pylab.figure(3,figsize=(3.42,4))
    fig = pylab.gcf()
    # Set up grids to hold figure panels
    outer_grid = gridspec.GridSpec(2, 1, height_ratios=[1,1], hspace=0.05)


    absolute_axis = plt.Subplot(fig, outer_grid[0])
    fig.add_subplot(absolute_axis)

    relative_axis = plt.Subplot(fig, outer_grid[1])
    fig.add_subplot(relative_axis)

    absolute_axis.set_title(species_name,fontsize=7)
    absolute_axis.set_ylabel('% inconsistent SNPs')
    
    absolute_axis.set_xlim([max_ds[0]/1.1,max_ds[-1]*1.1])
    absolute_axis.set_ylim([0,1.1])
    
    relative_axis.set_ylabel('Relative to unlinked')
    relative_axis.set_xlabel('Divergence threshold, $d$')
    relative_axis.set_xlim([max_ds[0]/1.1,max_ds[-1]*1.1])
    relative_axis.set_ylim([0,1.1])
    
    absolute_axis.semilogx(max_ds, fraction_inconsistent, 'b.-',label='Observed')
    absolute_axis.plot(max_ds, null_fraction_inconsistent, 'k.-',linewidth=0.25,color='0.7',label='Unlinked')

    relative_axis.semilogx(max_ds, fraction_inconsistent/null_fraction_inconsistent, 'b.-')
    relative_axis.plot(max_ds, numpy.ones_like(max_ds), 'k-',linewidth=0.25,color='0.7')
    
    absolute_axis.set_xticklabels([])
    
    absolute_axis.legend(loc='upper left', frameon=False,fontsize=6)
    
    pylab.savefig('%s/%s_phylogenetic_inconsistency.pdf' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight')
    pylab.savefig('%s/%s_phylogenetic_inconsistency.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight',dpi=300)

    sys.exit(0)
    
    polymorphic_freqs = numpy.array(polymorphic_freqs)
    inconsistent_freqs = numpy.array(inconsistent_freqs)
    
    total_polymorphic_sites = len(polymorphic_freqs)
    total_inconsistent_sites = len(inconsistent_freqs)
    null_inconsistent_sites = len(null_inconsistent_freqs)
    
    print total_polymorphic_sites, total_inconsistent_sites
    
    sys.stderr.write('%d SNPs, %d inconsistent (%g)\n' % (total_polymorphic_sites, total_inconsistent_sites, total_inconsistent_sites*1.0/total_polymorphic_sites))
    sys.stderr.write('Null = %d (%g)\n' % (null_inconsistent_sites, null_inconsistent_sites*1.0/total_polymorphic_sites))
    pylab.figure(3,figsize=(3.42,2))
    pylab.suptitle(species_name)

    xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(polymorphic_freqs)
    pylab.step(xs,ns*1.0/ns[0],'b-',label='All polymorphisms')
    if len(inconsistent_freqs)>0:
        xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(inconsistent_freqs)
        pylab.step(xs,ns*1.0/ns[0],'r-',label=('Inconsistent ($d=%g$)' % max_clade_d))
    pylab.xlim([0,0.5])
    pylab.ylim([0,1])
    pylab.xlabel('Within-clade MAF, $f$')
    pylab.ylabel('SNPs $\geq f$')
    pylab.legend(loc='upper right', frameon=False,fontsize=6)
    
    pylab.savefig('%s/%s_phylogenetic_inconsistency.pdf' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight')
    pylab.savefig('%s/%s_phylogenetic_inconsistency.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight',dpi=300)



else:
    
    sys.stderr.write("No clades!\n")

    
