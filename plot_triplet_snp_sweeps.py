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

variant_types = ['1D','4D']

modification_difference_threshold = 1000

# Must compete divergence matrix on the fly! 
            
# Load subject and sample metadata
sys.stderr.write("Loading sample metadata...\n")
subject_sample_map = sample_utils.parse_subject_sample_map()
sample_order_map = sample_utils.parse_sample_order_map()
sys.stderr.write("Done!\n")
    
good_species_list = parse_midas_data.parse_good_species_list()
if debug:
    good_species_list = good_species_list[0:2]

triplet_map = {}
for species_name in good_species_list:

    haploid_triplets = []

    # Only plot samples above a certain depth threshold that are "haploids"
    haploid_samples = diversity_utils.calculate_haploid_samples(species_name, debug=debug)

    sample_order_map = sample_utils.parse_sample_order_map()
    # Calculate which triplets of idxs belong to the same subject
    same_subject_idxs = sample_utils.calculate_ordered_subject_triplets(sample_order_map, haploid_samples)
    
    temporal_samples = set()
    for sample_triplet_idx in xrange(0,len(same_subject_idxs)):
        i,j,k = same_subject_idxs[sample_triplet_idx]
        
        temporal_samples.add(haploid_samples[i])
        temporal_samples.add(haploid_samples[j])
        temporal_samples.add(haploid_samples[k])
        
        haploid_triplets.append( (species_name, haploid_samples[i], haploid_samples[j], haploid_samples[k]) )
     
    if len(haploid_triplets) > 4:
        print species_name, len(haploid_triplets)
        triplet_map[species_name] = haploid_triplets
        
# Now can plot them the same way as before    
total_panels = 0
for species_name in triplet_map.keys():
    total_panels+=len(triplet_map[species_name])
   
print total_panels

# Make plot. title species           
 
####################################################
#
# Set up Figure (3 main panels, arranged in 1x3 grid)
#
####################################################

pylab.figure(1,figsize=(3.42,2*total_panels))
fig = pylab.gcf()
# make three panels panels

outer_grid  = gridspec.GridSpec(total_panels, 1, height_ratios=[1 for i in xrange(0,total_panels)],hspace=0.2)


sfs_axes = []


for panel_idx in xrange(0,total_panels):

    sfs_axis = plt.Subplot(fig, outer_grid[panel_idx])
    fig.add_subplot(sfs_axis)
 
    sfs_axis.set_xlim([0.5,3.5])
    sfs_axis.set_ylim([-0.05,1.05])
    sfs_axes.append(sfs_axis)
    
# Now actually plot
panel_idx=-1
for species_name in good_species_list:

    if species_name not in triplet_map:
        continue
    
    import core_gene_utils
    sys.stderr.write("Loading whitelisted genes...\n")
    non_shared_genes = core_gene_utils.parse_non_shared_reference_genes(species_name)
    shared_pangenome_genes = core_gene_utils.parse_shared_genes(species_name)
    sys.stderr.write("Done! %d shared genes and %d non-shared genes\n" % (len(shared_pangenome_genes), len(non_shared_genes)))
        
    
    triplet_freqs = []    
    snp_samples = set()
    for dummy1, sample_i, sample_j, sample_k in triplet_map[species_name]:
        snp_samples.add(sample_i)
        snp_samples.add(sample_j)
        snp_samples.add(sample_k)
        triplet_freqs.append([])

    snp_samples = list(snp_samples)
    
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
    
    
        desired_triplet_idxs = []
        for dummy1, initial_sample, middle_sample, final_sample in triplet_map[species_name]:
            initial_idx = numpy.nonzero(dummy_samples==initial_sample)[0][0]
            middle_idx = numpy.nonzero(dummy_samples==middle_sample)[0][0]
            final_idx = numpy.nonzero(dummy_samples==final_sample)[0][0]
    
            desired_triplet_idxs.append((initial_idx, middle_idx, final_idx))
        
    
        # Calculate fixation matrix
        sys.stderr.write("Calculating joint freqs...\n")
        for pair_idx in xrange(0,len(triplet_map[species_name])):
        
            initial_sample_idx, middle_sample_idx, final_sample_idx = desired_triplet_idxs[pair_idx]
        
            chunk_freqs = diversity_utils.calculate_triplet_sample_freqs(allele_counts_map, passed_sites_map, initial_sample_idx, middle_sample_idx, final_sample_idx)
            triplet_freqs[pair_idx].extend( chunk_freqs ) 
            
        sys.stderr.write("Done!\n")
    
    sys.stderr.write("Done!\n")   

    # Plot joint SFS!
    for pair_idx in xrange(0,len(triplet_map[species_name])):
    
        
        dummy, sample_i, sample_j, sample_k = triplet_map[species_name][pair_idx]
    
        panel_idx+=1 

        freqs = triplet_freqs[pair_idx]
        
        visnos = numpy.array([1,2,3])
    
        for i in xrange(0,len(freqs)):
            sfs_axes[panel_idx].plot(visnos, freqs[i],'.-',alpha=0.5,markersize=2,markeredgewidth=0)
    
    
        sfs_axes[panel_idx].set_ylabel('%s\n%s' % (figure_utils.get_abbreviated_species_name(species_name), sample_i))
    
fig.savefig('%s/haploid_triplet_sweeps.png' % (parse_midas_data.analysis_directory),bbox_inches='tight',dpi=300)
sys.stderr.write("Done!\n")
