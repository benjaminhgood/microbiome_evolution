import matplotlib  
matplotlib.use('Agg') 
import config
import parse_midas_data
import sample_utils
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
import figure_utils
    
import matplotlib.colors as colors
import matplotlib.cm as cmx
from math import log10,ceil
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy.random import randint
import matplotlib.colors as mcolors

from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster.hierarchy import cophenet
from scipy.cluster.hierarchy import fcluster


mpl.rcParams['font.size'] = 6
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

args = parser.parse_args()

debug = args.debug
chunk_size = args.chunk_size
memoize = args.memoize

################################################################################

min_coverage = config.min_median_coverage
min_sample_size = 5

# Must compete divergence matrix on the fly! 
            
# Load subject and sample metadata
sys.stderr.write("Loading sample metadata...\n")
subject_sample_map = sample_utils.parse_subject_sample_map()
sample_order_map = sample_utils.parse_sample_order_map()
sys.stderr.write("Done!\n")

examples = [['Bacteroides_vulgatus_57955',  '700021876'], 
            ['Bacteroides_uniformis_57318', '700016456'],
            ['Bacteroides_caccae_53434',    '700024998']]

for example_idx in xrange(0,len(examples)):
    species_name = examples[example_idx][0]
    initial_sample = examples[example_idx][1]

    # Only plot samples above a certain depth threshold that are "haploids"
    haploid_samples = diversity_utils.calculate_haploid_samples(species_name, debug=debug)

    sample_order_map = sample_utils.parse_sample_order_map()
    # Calculate which triplets of idxs belong to the same subject
    same_subject_idxs = parse_midas_data.calculate_ordered_subject_triplets(sample_order_map, haploid_samples)
    
    temporal_samples = set()
    for sample_triplet_idx in xrange(0,len(same_subject_idxs)):
        i,j,k = same_subject_idxs[sample_triplet_idx]
        
        if haploid_samples[i] == initial_sample:
            examples[example_idx].append( (haploid_samples[i],haploid_samples[j],haploid_samples[k]) )
        
# Now can plot them the same way as before    
total_panels = len(examples)

# Make plot. title species           
 
####################################################
#
# Set up Figure (3 main panels, arranged in 1x3 grid)
#
####################################################

pylab.figure(1,figsize=(7,2.75))
fig = pylab.gcf()
# make three panels panels

outer_grid  = gridspec.GridSpec(2, total_panels, height_ratios=[1,1], hspace=0.05, width_ratios=[1 for i in xrange(0,total_panels)], wspace=0.05)


snp_axes = []
gene_axes = []

for example_idx in xrange(0,total_panels):

    snp_axis = plt.Subplot(fig, outer_grid[0, example_idx])
    fig.add_subplot(snp_axis)
    line, = snp_axis.plot([0.5,3.5],[0.2,0.2],':',color='0.7')
    line.set_dashes((1,1))
    line, = snp_axis.plot([0.5,3.5],[0.8,0.8],':',color='0.7')
    line.set_dashes((1,1))
    snp_axis.set_xlim([0.6,3.4])
    snp_axis.set_xticks([1,2,3])
    snp_axis.set_ylim([-0.02,1.02])
    snp_axis.set_xticklabels([])
    snp_axes.append(snp_axis)
    
    
    gene_axis = plt.Subplot(fig, outer_grid[1, example_idx])
    fig.add_subplot(gene_axis)
    line, = gene_axis.semilogy([0.5,3.5],[0.5,0.5],':',color='0.7')
    line.set_dashes((1,1))
    line, = gene_axis.semilogy([0.5,3.5],[0.05,0.05],':',color='0.7')
    line.set_dashes((1,1))
    gene_axis.set_xlim([0.6,3.4])
    gene_axis.set_xticks([1,2,3])
    gene_axis.set_ylim([0.02,2])
    gene_axes.append(gene_axis)
    
    snp_axis.set_title(figure_utils.get_pretty_species_name(examples[example_idx][0]),fontsize=6)
    gene_axis.set_xlabel('Visit number')
    if example_idx==0:
        snp_axis.set_ylabel('SNV frequency')
        gene_axis.set_ylabel('Gene copynum')
    else:
        snp_axis.set_yticklabels([])
        gene_axis.set_yticklabels([])
    

# Now actually plot
for example_idx in xrange(0,len(examples)):

    species_name = examples[example_idx][0]
    
    import core_gene_utils
    sys.stderr.write("Loading whitelisted genes...\n")
    non_shared_genes = core_gene_utils.parse_non_shared_reference_genes(species_name)
    shared_pangenome_genes = core_gene_utils.parse_shared_genes(species_name)
    sys.stderr.write("Done! %d shared genes and %d non-shared genes\n" % (len(shared_pangenome_genes), len(non_shared_genes)))
    
    initial_sample, middle_sample, final_sample = examples[example_idx][2]
    snp_samples = [initial_sample, middle_sample, final_sample]
    
    triplet_freqs = []
    # Analyze SNPs, looping over chunk sizes. 
    # Clunky, but necessary to limit memory usage on cluster
    # Load SNP information for species_name
    sys.stderr.write("Loading SNPs for %s...\n" % species_name)
    sys.stderr.write("(not just core genes...)\n") 
    final_line_number = 0
    while final_line_number >= 0:
    
        sys.stderr.write("Loading chunk starting @ %d...\n" % final_line_number)
        dummy_samples, allele_counts_map, passed_sites_map, final_line_number = parse_midas_data.parse_snps(species_name, debug=debug, allowed_samples=snp_samples,chunk_size=chunk_size,initial_line_number=final_line_number, allowed_genes=non_shared_genes)
        sys.stderr.write("Done! Loaded %d genes\n" % len(allele_counts_map.keys()))
        print len(dummy_samples), "dummy samples!"
    
    
        initial_sample_idx = numpy.nonzero(dummy_samples==initial_sample)[0][0]
        middle_sample_idx = numpy.nonzero(dummy_samples==middle_sample)[0][0]
        final_sample_idx = numpy.nonzero(dummy_samples==final_sample)[0][0]
    
        chunk_freqs = diversity_utils.calculate_triplet_sample_freqs(allele_counts_map, passed_sites_map, initial_sample_idx, middle_sample_idx, final_sample_idx)
        triplet_freqs.extend( chunk_freqs ) 
        
    sys.stderr.write("Done!\n")   

        
    visnos = numpy.array([1,2,3])
    
    for i in xrange(0,len(triplet_freqs)):
        snp_axes[example_idx].plot(visnos, triplet_freqs[i],'.-',alpha=0.5,markersize=2,markeredgewidth=0,rasterized=True)
    
    
    # Now do genes
    # Load gene coverage information for species_name
    sys.stderr.write("Loading pangenome data for %s...\n" % species_name)
    dummy_samples, gene_names, gene_presence_matrix, gene_depth_matrix, marker_coverages, gene_reads_matrix =     parse_midas_data.parse_pangenome_data(species_name,allowed_samples=snp_samples, disallowed_genes=shared_pangenome_genes)
    sys.stderr.write("Done!\n")
    
    
    initial_idx = numpy.nonzero(dummy_samples==initial_sample)[0][0]
    middle_idx = numpy.nonzero(dummy_samples==middle_sample)[0][0]
    final_idx = numpy.nonzero(dummy_samples==final_sample)[0][0]
    
    copynum_trajectories = gene_diversity_utils.calculate_triplet_gene_copynums(gene_depth_matrix, marker_coverages, initial_idx, middle_idx, final_idx)
    
    for i in xrange(0,len(copynum_trajectories)):
            
        cs = numpy.clip(copynum_trajectories[i],1e-03,2)
                
        gene_axes[example_idx].plot(visnos, cs,'.-',alpha=0.5,markersize=2,markeredgewidth=0,rasterized=True)
    
    
fig.savefig('%s/supplemental_triplet_trajectories.png' % (parse_midas_data.analysis_directory),bbox_inches='tight',dpi=300)
sys.stderr.write("Done!\n")
