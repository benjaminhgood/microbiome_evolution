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
    
    snp_samples = dummy_samples
     
substitution_rate = snp_difference_matrix*1.0/snp_opportunity_matrix 

max_ds = numpy.array([1e-05, 1e-04, 1e-03, 1e-02])

#max_ds = numpy.logspace(-4,-1,19)
fine_grained_idxs, fine_grained_cluster_list = diversity_utils.cluster_samples(substitution_rate, min_d=1e-07, max_ds=max_ds)


# Plot phylogenetic consistency vs time
remapped_cluster_list = []
for cluster_idxss in fine_grained_cluster_list:

    remapped_cluster_idxss = [cluster_idxs[fine_grained_idxs] for cluster_idxs in cluster_idxss]
    remapped_cluster_list.append(remapped_cluster_idxss)

fine_grained_samples = snp_samples[fine_grained_idxs]

#maf_bins = numpy.arange(1,len(fine_grained_samples)+1)*1.0/len(fine_grained_samples)
#maf_bins -= (maf_bins[1]-maf_bins[0])/2
#maf_bins[0]=-0.1
#maf_bins[-1] = 1.1
#mafs = numpy.arange(1,len(fine_grained_samples))*1.0/len(fine_grained_samples) 


if len(fine_grained_samples)>2:

    sys.stderr.write("Continuing with %d samples...\n" % len(fine_grained_samples))

    # Load SNP information for species_name
    sys.stderr.write("Re-loading %s...\n" % species_name)

    snp_difference_matrix = numpy.array([])
    snp_opportunity_matrix = numpy.array([])

    total_singleton_sites = [0 for max_d in max_ds]
    total_polymorphic_sites = [0 for max_d in max_ds]
    total_inconsistent_sites = [0 for max_d in max_ds]
    total_null_inconsistent_sites = [0 for max_d in max_ds]
 
    polymorphic_variant_types = [{'1D':0, '2D':0, '3D':0, '4D':0} for max_d in max_ds]
    inconsistent_variant_types = [{'1D':0, '2D':0, '3D':0, '4D':0} for max_d in max_ds]
 
    polymorphic_freqs = [[] for max_d in max_ds]
    inconsistent_freqs = [[] for max_d in max_ds]
 
    #polymorphic_sfs = [numpy.zeros_like(mafs) for max_d in max_ds]
    #inconsistent_sfs = [numpy.zeros_like(mafs) for max_d in max_ds]
 
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
        
            chunk_singleton_freqs, chunk_polymorphic_freqs, chunk_inconsistent_freqs, chunk_null_inconsistent_freqs, chunk_singleton_variant_types, chunk_polymorphic_variant_types, chunk_inconsistent_variant_types, chunk_null_variant_types = diversity_utils.calculate_phylogenetic_consistency(allele_counts_map, passed_sites_map, cluster_idxss, allowed_genes=core_genes)    
        
            total_singleton_sites[i] += len(chunk_singleton_freqs)
            total_polymorphic_sites[i] += len(chunk_polymorphic_freqs) 
            total_inconsistent_sites[i] += len(chunk_inconsistent_freqs)
            total_null_inconsistent_sites[i] += len(chunk_null_inconsistent_freqs)
    
            polymorphic_freqs[i].extend(chunk_polymorphic_freqs)
            inconsistent_freqs[i].extend(chunk_inconsistent_freqs)
            
            print "Singleton:", chunk_singleton_variant_types
            print "Polymorphic:", chunk_polymorphic_variant_types
            print "Inconsistent:", chunk_inconsistent_variant_types
                
            for variant_type in polymorphic_variant_types[i].keys():
                singleton_variant_types[i][variant_type] += chunk_singleton_variant_types[variant_type]
                polymorphic_variant_types[i][variant_type] += chunk_polymorphic_variant_types[variant_type]
                inconsistent_variant_types[i][variant_type] += chunk_inconsistent_variant_types[variant_type]
                
            
            
        sys.stderr.write("Done!\n")
    
    substitution_rate = snp_difference_matrix*1.0/snp_opportunity_matrix 
    
        
    total_polymorphic_sites = numpy.array(total_polymorphic_sites)
    total_inconsistent_sites = numpy.array(total_inconsistent_sites)
    total_null_inconsistent_sites = numpy.array(total_null_inconsistent_sites)
  
    fraction_inconsistent = total_inconsistent_sites*1.0/total_polymorphic_sites
    null_fraction_inconsistent = total_null_inconsistent_sites*1.0/total_polymorphic_sites
    
    print total_singleton_sites
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

    for i in xrange(0,len(polymorphic_freqs)):
        polymorphic_freqs[i] = numpy.array(polymorphic_freqs[i])
        inconsistent_freqs[i] = numpy.array(inconsistent_freqs[i])
        print "d=", max_ds[i]
        print "Site", "Polymorphic", "Inconsistent"
        for variant_type in sorted(polymorphic_variant_types[i].keys()):
            variant_type, polymorphic_variant_types[i][variant_type], inconsistent_variant_types[i][variant_type] 
        print ""
    
    pylab.figure(4,figsize=(3.42,2))
    pylab.suptitle(species_name)

    for i in xrange(0,len(polymorphic_freqs)):
    
        if len(polymorphic_freqs[i])==0:
            continue
    
        xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(polymorphic_freqs[i])
        pylab.step(xs,ns*1.0/ns[0],'-',label='Polymorphic ($d=%g$)' % max_ds[i])
        print 'Polymorphic (d=%g), n=%g' % (max_ds[i], ns[0])
    
    if len(inconsistent_freqs[1])>0:
        xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(inconsistent_freqs[1])
        pylab.step(xs,ns*1.0/ns[0],'r-',linewidth=2,label=('Inconsistent ($d=%g$)' % max_ds[1]))
        
        print "Inconsistent (d=%g), n=%g" %  (max_ds[1], ns[0])
     
        
    pylab.xlim([0,0.5])
    pylab.ylim([0,1])
    pylab.xlabel('Minor allele frequency, $f$')
    pylab.ylabel('SNPs $\geq f$')
    pylab.legend(loc='upper right', frameon=False,fontsize=6)
    
    pylab.savefig('%s/%s_phylogenetic_inconsistency_sfs.pdf' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight')
    


else:
    
    sys.stderr.write("No clades!\n")

    
