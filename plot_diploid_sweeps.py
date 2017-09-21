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

modification_difference_threshold = 100

# Must compete divergence matrix on the fly! 
            
# Load subject and sample metadata
sys.stderr.write("Loading sample metadata...\n")
subject_sample_map = parse_HMP_data.parse_subject_sample_map()
sample_order_map = parse_HMP_data.parse_sample_order_map()
sys.stderr.write("Done!\n")
    
good_species_list = parse_midas_data.parse_good_species_list()
if debug:
    good_species_list = good_species_list[0:2]

diploid_pair_map = {}
for species_name in good_species_list:

    diploid_pairs = []

    # Only plot samples above a certain depth threshold that are "haploids"
    haploid_samples = diversity_utils.calculate_haploid_samples(species_name, debug=debug)

    if len(haploid_samples) < min_sample_size:
        continue
        
    # all samples
    temporal_samples = diversity_utils.calculate_temporal_samples(species_name)
    
    # temporal changes    
    sys.stderr.write("Loading pre-computed temporal changes for %s...\n" % species_name)
    temporal_change_map = calculate_temporal_changes.load_temporal_change_map(species_name)
    sys.stderr.write("Done!\n")

    same_sample_idxs, same_subject_idxs, diff_subject_idxs = parse_midas_data.calculate_ordered_subject_pairs(sample_order_map, temporal_samples)
    
    snp_samples = set()
    sample_size = 0        
    for sample_pair_idx in xrange(0,len(same_subject_idxs[0])):
   
        i = same_subject_idxs[0][sample_pair_idx]
        j = same_subject_idxs[1][sample_pair_idx]
        
        sample_i = temporal_samples[i]
        sample_j = temporal_samples[j]

        # Don't look at haploids here!
        if sample_i in haploid_samples:
            continue

        if sample_j in haploid_samples:
            final_state = "haploid"
        else:
            final_state = "diploid"

        perr, mutations, reversions = calculate_temporal_changes.calculate_mutations_reversions_from_temporal_change_map(temporal_change_map, sample_i, sample_j)
    
        if perr<-0.5:
            continue
    
        highconfidence_snps = []
        highlighted_gene_set = set()
        gene_fixation_positions = {}
        for snp_change in (mutations+reversions):
            gene_name = snp_change[0]
            position = snp_change[2]
        
            if gene_name not in gene_fixation_positions:
                gene_fixation_positions[gene_name] = []
            gene_fixation_positions[gene_name].append(position) 
        
        for gene_name in gene_fixation_positions:
            if len(gene_fixation_positions[gene_name]) >= 2:
                # Calculate max position difference between SNPs.
            
                positions = numpy.array(gene_fixation_positions[gene_name])
            
                max_distance = numpy.fabs(positions[:,None]-positions[None,:]).max()
            
                if max_distance>100:
                    highlighted_gene_set.add(gene_name)
    
        for snp_change in (mutations+reversions):
            if snp_change[0] in highlighted_gene_set:
                highconfidence_snps.append(snp_change)
        
        if len(highconfidence_snps)>0:        
            
            # If more than 10 genes, probably a replacement event
            if len(highlighted_gene_set) < 10:
                diploid_pairs.append((species_name, sample_i, sample_j, final_state, highlighted_gene_set, highconfidence_snps))
                print species_name, sample_i, sample_j, final_state, len(highlighted_gene_set)
                for snp in highconfidence_snps:
                    print snp
            else:
                pass
    if len(diploid_pairs)>0:
        diploid_pair_map[species_name] = diploid_pairs
        
# Now can plot them the same way as before    
total_panels = 0
for species_name in diploid_pair_map.keys():
    total_panels+=len(diploid_pair_map[species_name])
   
print len(diploid_pair_map.keys()), total_panels 
    
# Make plot. title species, sample_i, sample_j           
 
####################################################
#
# Set up Figure (3 main panels, arranged in 1x3 grid)
#
####################################################

pylab.figure(1,figsize=(3,2*total_panels))
fig = pylab.gcf()
# make three panels panels

outer_grid  = gridspec.GridSpec(total_panels, 1, height_ratios=[1 for i in xrange(0,total_panels)],hspace=0.2)


sfs_axes = []


for panel_idx in xrange(0,total_panels):

    sfs_axis = plt.Subplot(fig, outer_grid[panel_idx])
    fig.add_subplot(sfs_axis)
 
    sfs_axis.set_xlim([-0.05,1.05])
    sfs_axis.set_ylim([-0.05,1.05])
    #sfs_axis.plot([0,1],[1,1],'k-')
    #sfs_axis_1.plot([0,0],[0,1],'k-')
    #sfs_axis_1.plot([0,0.2],[0.8,0.8],'k-')
    #sfs_axis_1.plot([0.2,0.2],[0.8,1.0],'k-')

    sfs_axis.plot([0,1],[0,1.0],'k-')
    
    sfs_axes.append(sfs_axis)
    
# Now actually plot
panel_idx=-1
for species_name in good_species_list:

    if species_name not in diploid_pair_map:
        continue
        
    gene_names = [[] for i in xrange(0,len(diploid_pair_map[species_name]))]
    initial_freqs = [[] for i in xrange(0,len(diploid_pair_map[species_name]))]
    final_freqs = [[] for i in xrange(0,len(diploid_pair_map[species_name]))]

    snp_samples = set()
    for dummy1, sample_i, sample_j, dummy2, dummy3, dummy4 in diploid_pair_map[species_name]:
        snp_samples.add(sample_i)
        snp_samples.add(sample_j)

    snp_samples = list(snp_samples)
    
    # Analyze SNPs, looping over chunk sizes. 
    # Clunky, but necessary to limit memory usage on cluster

    # Load SNP information for species_name
    sys.stderr.write("Loading SNPs for %s...\n" % species_name)
    sys.stderr.write("(not just core genes...)\n")
    
    final_line_number = 0
    while final_line_number >= 0:
    
        sys.stderr.write("Loading chunk starting @ %d...\n" % final_line_number)
        dummy_samples, allele_counts_map, passed_sites_map, final_line_number = parse_midas_data.parse_snps(species_name, debug=debug, allowed_samples=snp_samples,chunk_size=chunk_size,initial_line_number=final_line_number)
        sys.stderr.write("Done! Loaded %d genes\n" % len(allele_counts_map.keys()))
        print len(dummy_samples), "dummy samples!"
    
    
        desired_sample_pair_idxs = []
        for dummy1, initial_sample, final_sample, dummy2, dummy3, dummy4 in diploid_pair_map[species_name]:
            initial_idx = numpy.nonzero(dummy_samples==initial_sample)[0][0]
            final_idx = numpy.nonzero(dummy_samples==final_sample)[0][0]
    
            desired_sample_pair_idxs.append((initial_idx, final_idx))
        
    
        # Calculate fixation matrix
        sys.stderr.write("Calculating joint freqs...\n")
        for pair_idx in xrange(0,len(diploid_pair_map[species_name])):
        
            initial_sample_idx, final_sample_idx = desired_sample_pair_idxs[pair_idx]
        
            chunk_gene_names, chunk_chromosomes, chunk_positions, chunk_initial_freqs, chunk_final_freqs, chunk_initial_depths, chunk_final_depths = diversity_utils.calculate_temporal_sample_freqs(allele_counts_map, passed_sites_map, initial_sample_idx, final_sample_idx)  # 
    
            joint_passed_sites = (chunk_initial_depths>0)*(chunk_final_depths>0)
    
            gene_names[pair_idx].extend(chunk_gene_names[joint_passed_sites])
            initial_freqs[pair_idx].extend( chunk_initial_freqs[joint_passed_sites])
            final_freqs[pair_idx].extend( chunk_final_freqs[joint_passed_sites])
         
        sys.stderr.write("Done!\n")
    
    for pair_idx in xrange(0,len(diploid_pair_map[species_name])):
        gene_names[pair_idx] = numpy.array(gene_names[pair_idx])
        initial_freqs[pair_idx] = numpy.array(initial_freqs[pair_idx])
        final_freqs[pair_idx] = numpy.array(final_freqs[pair_idx])
    
    sys.stderr.write("Done!\n")   

    # Plot joint SFS!
    for pair_idx in xrange(0,len(diploid_pair_map[species_name])):
    
        dummy, sample_i, sample_j, final_state, highlighted_gene_names, highconfidence_snps = diploid_pair_map[species_name][pair_idx]
    
        panel_idx+=1 


        f0s = initial_freqs[pair_idx]
        f1s = final_freqs[pair_idx]
    
        dfs = (f1s-f0s)
    
        abs_dfs = numpy.fabs(dfs)
    
        major_freqs = numpy.fmax(f1s,1-f1s)
    
        # Don't plot a bunch of crap near 0
        good_sites = numpy.logical_or(f0s>0.05, f1s>0.05)
    
    
        sfs_axes[panel_idx].plot(f0s[good_sites], f1s[good_sites],'.',alpha=0.1,markersize=2,markeredgewidth=0,color='0.7')
    
        # Now plot ones in fixed gene set with colors
        for gene_name in sorted(highlighted_gene_names):
        
            in_gene = (gene_names[pair_idx]==gene_name)*good_sites
         
            sfs_axes[panel_idx].plot(f0s[in_gene], f1s[in_gene],'.',alpha=0.75,markersize=4,markeredgewidth=0)
    
        
        sfs_axes[panel_idx].set_ylabel('%s\n%s,%s' % (figure_utils.get_abbreviated_species_name(species_name), sample_i, sample_j))
    
fig.savefig('%s/all_diploid_sweeps.png' % (parse_midas_data.analysis_directory),bbox_inches='tight',dpi=300)
sys.stderr.write("Done!\n")
