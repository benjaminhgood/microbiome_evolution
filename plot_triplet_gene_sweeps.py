import matplotlib  
matplotlib.use('Agg') 
import config
import parse_midas_data
import parse_HMP_data
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
subject_sample_map = parse_HMP_data.parse_subject_sample_map()
sample_order_map = parse_HMP_data.parse_sample_order_map()
sys.stderr.write("Done!\n")
    
good_species_list = parse_midas_data.parse_good_species_list()
if debug:
    good_species_list = good_species_list[0:2]

triplet_map = {}
for species_name in good_species_list:

    

    haploid_triplets = []

    # Only plot samples above a certain depth threshold that are "haploids"
    haploid_samples = diversity_utils.calculate_haploid_samples(species_name, debug=debug)

    sample_order_map = parse_HMP_data.parse_sample_order_map()
    # Calculate which triplets of idxs belong to the same subject
    same_subject_idxs = parse_midas_data.calculate_ordered_subject_triplets(sample_order_map, haploid_samples)
    
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
    line, = sfs_axis.semilogy([0.5,3.5],[0.5,0.5],':',color='0.7')
    line.set_dashes((1,1))
    line, = sfs_axis.semilogy([0.5,3.5],[0.05,0.05],':',color='0.7')
    line.set_dashes((1,1))
    
    sfs_axis.set_xlim([0.5,3.5])
    sfs_axis.set_ylim([0.02,2])
    sfs_axes.append(sfs_axis)
    
# Now actually plot
panel_idx=-1
for species_name in good_species_list:

    if species_name not in triplet_map:
        continue
        
    
    triplet_freqs = []    
    snp_samples = set()
    for dummy1, sample_i, sample_j, sample_k in triplet_map[species_name]:
        snp_samples.add(sample_i)
        snp_samples.add(sample_j)
        snp_samples.add(sample_k)
        triplet_freqs.append([])

    snp_samples = list(snp_samples)
    
    import core_gene_utils
    sys.stderr.write("Loading whitelisted genes...\n")
    non_shared_genes = core_gene_utils.parse_non_shared_reference_genes(species_name)
    shared_pangenome_genes = core_gene_utils.parse_shared_genes(species_name)
    sys.stderr.write("Done! %d shared genes and %d non-shared genes\n" % (len(shared_pangenome_genes), len(non_shared_genes)))
    
    
    # Analyze SNPs, looping over chunk sizes. 
    # Clunky, but necessary to limit memory usage on cluster


    # Load gene coverage information for species_name
    sys.stderr.write("Loading pangenome data for %s...\n" % species_name)
    dummy_samples, gene_names, gene_presence_matrix, gene_depth_matrix, marker_coverages, gene_reads_matrix =     parse_midas_data.parse_pangenome_data(species_name,allowed_samples=snp_samples, disallowed_genes=shared_pangenome_genes)
    sys.stderr.write("Done!\n")
    
    desired_triplet_idxs = []
    for dummy1, initial_sample, middle_sample, final_sample in triplet_map[species_name]:
        panel_idx+=1 

        if (initial_sample in dummy_samples) and (middle_sample in dummy_samples) and (final_sample in dummy_samples):
            
            initial_idx = numpy.nonzero(dummy_samples==initial_sample)[0][0]
            middle_idx = numpy.nonzero(dummy_samples==middle_sample)[0][0]
            final_idx = numpy.nonzero(dummy_samples==final_sample)[0][0]
    
            copynum_trajectories = gene_diversity_utils.calculate_triplet_gene_copynums(gene_depth_matrix, marker_coverages, initial_idx, middle_idx, final_idx)
    
            visnos = numpy.array([1,2,3])
    
            for i in xrange(0,len(copynum_trajectories)):
            
                cs = numpy.clip(copynum_trajectories[i],5e-02,2)
                
                sfs_axes[panel_idx].plot(visnos, cs,'.-',alpha=0.5,markersize=2,markeredgewidth=0)
    
    
            sfs_axes[panel_idx].set_ylabel('%s\n%s' % (figure_utils.get_abbreviated_species_name(species_name), initial_sample))
    
fig.savefig('%s/haploid_triplet_gene_sweeps.png' % (parse_midas_data.analysis_directory),bbox_inches='tight',dpi=300)
sys.stderr.write("Done!\n")
