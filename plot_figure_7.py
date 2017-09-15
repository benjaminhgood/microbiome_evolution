# Within-snp gene changes

import matplotlib  
matplotlib.use('Agg') 
import config
import parse_midas_data
###
#
# For today while the new data processes
#
import os
#parse_midas_data.data_directory = os.path.expanduser("~/ben_nandita_hmp_data_062517/")
#########################################
import parse_HMP_data


import pylab
import sys
import numpy

import diversity_utils
import gene_diversity_utils
import calculate_temporal_changes

import stats_utils
import matplotlib.colors as colors
import matplotlib.cm as cmx
from math import log10,ceil
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy.random import randint, choice

mpl.rcParams['font.size'] = 5
mpl.rcParams['lines.linewidth'] = 0.5
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'

species_name = "Bacteroides_vulgatus_57955"
#species_name = 'Bacteroides_uniformis_57318'
#species_name = 'Prevotella_copri_61740'
################################################################################
#
# Standard header to read in argument information
#
################################################################################
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--debug", help="Loads only a subset of SNPs for speed", action="store_true")
parser.add_argument("--chunk-size", type=int, help="max number of records to load", default=1000000000)
parser.add_argument('--other-species', type=str, help='Run the script for a different species')

args = parser.parse_args()

debug = args.debug
chunk_size = args.chunk_size
other_species = args.other_species

if other_species:
    species_name = other_species
    other_species_str = "_%s" % species_name
else:
    other_species_str = ""
    
################################################################################

min_coverage = config.min_median_coverage
alpha = 0.5 # Confidence interval range for rate estimates
low_pi_threshold = 1e-03
clade_divergence_threshold = 1e-02
modification_divergence_threshold = 1e-03
min_change = 0.8
include_high_copynum = False
#include_high_copynum = True

snv_lower_threshold = 0.2
snv_upper_threshold = 0.8

# half-width of window to zoom in on for gene panels
window_size = 5000

# Load subject and sample metadata
sys.stderr.write("Loading sample metadata...\n")
subject_sample_map = parse_HMP_data.parse_subject_sample_map()
sample_country_map = parse_HMP_data.parse_sample_country_map()
sample_order_map = parse_HMP_data.parse_sample_order_map()
sys.stderr.write("Done!\n")
  
  
desired_sample_pairs = {}

desired_sample_pairs['Bacteroides_uniformis_57318'] = []
desired_sample_pairs['Bacteroides_uniformis_57318'].append(('700023113', '700023720'))
desired_sample_pairs['Bacteroides_uniformis_57318'].append(('700015245',  '700099803'))
desired_sample_pairs['Bacteroides_uniformis_57318'].append(('700033502', '700102356'))

desired_sample_pairs['Prevotella_copri_61740'] = []
desired_sample_pairs['Prevotella_copri_61740'].append(('700024437',  '700106663'))
desired_sample_pairs['Prevotella_copri_61740'].append(('700024437',  '700106663'))
desired_sample_pairs['Prevotella_copri_61740'].append(('700024437',  '700106663'))
  

desired_sample_pairs['Bacteroides_vulgatus_57955'] = []
desired_sample_pairs['Bacteroides_vulgatus_57955'].append(('700024318', '700105372')) # diploid->diploid example # this first one is ignored
desired_sample_pairs['Bacteroides_vulgatus_57955'].append(('700016960', '700101581')) # diploid->haploid example?
desired_sample_pairs['Bacteroides_vulgatus_57955'].append(('700014837', '700098561')) # diploid->diploid example
#desired_sample_pairs['Bacteroides_vulgatus_57955'].append(('700024318', '700105372')) # diploid->diploid example

#435590.9.peg.1205       Putative outer membrane protein, probably involved in nutrient binding
#435590.9.peg.2541       Thiol-activated cytolysin


# Locations of SNP changes to zoom in on
desired_diploid_genes = ['435590.9.peg.1205',  '435590.9.peg.2541'] 
desired_locations = [('NC_009614', 1567636L), ('NC_009614', 3163600)]

#desired_locations = [('NZ_DS362239', 189868L), ('NZ_DS362220', 116578L)]

desired_sample_pairs = desired_sample_pairs[species_name]
snp_samples = set()
for initial_sample, final_sample in desired_sample_pairs:
    snp_samples.add(initial_sample)
    snp_samples.add(final_sample)
    
snp_samples = list(snp_samples)

sys.stderr.write("Loading pre-computed temporal changes for %s...\n" % species_name)
temporal_change_map = calculate_temporal_changes.load_temporal_change_map(species_name)
sys.stderr.write("Done!\n") 
       
####################################################
#
# Set up Figure (3 main panels, arranged in 1x3 grid)
#
####################################################

pylab.figure(1,figsize=(7,1.5))
fig = pylab.gcf()
# make three panels panels

outer_grid  = gridspec.GridSpec(1,2,width_ratios=[1.3,1],wspace=0.1)

transition_grid = gridspec.GridSpecFromSubplotSpec(1,2, width_ratios=[1,1], wspace=0.05, subplot_spec=outer_grid[0])

gene_grid = gridspec.GridSpecFromSubplotSpec(4,1, height_ratios=[1,1,1,1], hspace=0.3, subplot_spec=outer_grid[1])

sfs_axes = []

sfs_axis_1 = plt.Subplot(fig, transition_grid[0])
fig.add_subplot(sfs_axis_1)
sfs_axis_1.set_ylabel('Final frequency')
sfs_axis_1.set_xlabel('Initial frequency')
sfs_axis_1.set_xlim([-0.05,0.51])
sfs_axis_1.set_ylim([-0.05,1.05])
sfs_axis_1.plot([0,1],[1,1],'k-')
sfs_axis_1.plot([0,0],[0,1],'k-')
#sfs_axis_1.plot([0,0.2],[0.8,0.8],'k-')
#sfs_axis_1.plot([0.2,0.2],[0.8,1.0],'k-')

sfs_axis_1.plot([0,1],[0,1.0],'k-')

sfs_axis_1.spines['top'].set_visible(False)
sfs_axis_1.spines['right'].set_visible(False)
sfs_axis_1.get_xaxis().tick_bottom()
sfs_axis_1.get_yaxis().tick_left()
sfs_axes.append(sfs_axis_1)

sfs_axis_2 = plt.Subplot(fig, transition_grid[1])
fig.add_subplot(sfs_axis_2)
sfs_axis_2.set_xlabel('Initial frequency')
sfs_axis_2.set_xlim([-0.05,0.51])
sfs_axis_2.set_ylim([-0.05,1.05])
sfs_axis_2.plot([0,1],[1,1],'k-')
sfs_axis_2.plot([0,0],[0,1],'k-')
#sfs_axis_2.plot([0,0.2],[0.8,1.0],'k-')
sfs_axis_2.plot([0,1],[0,1.0],'k-')

sfs_axis_2.set_yticklabels([])
sfs_axis_2.spines['top'].set_visible(False)
sfs_axis_2.spines['right'].set_visible(False)
sfs_axis_2.get_xaxis().tick_bottom()
sfs_axis_2.get_yaxis().tick_left()
sfs_axes.append(sfs_axis_2)

depth_axes = []
gene_axes = []
for idx in xrange(0,len(desired_locations)):
    
    
    depth_axis = plt.Subplot(fig, gene_grid[2*idx])
    fig.add_subplot(depth_axis)
    depth_axis.set_ylabel('Depth')
    depth_axis.spines['top'].set_visible(False)
    depth_axis.spines['right'].set_visible(False)
    depth_axis.get_xaxis().tick_bottom()
    depth_axis.get_yaxis().tick_left()
    
    #depth_axis.set_ylim([-0.05,1.05])
    #depth_axis.set_yticks([0,0.5,1])
    depth_axis.set_xlim([desired_locations[idx][1]-window_size, desired_locations[idx][1]+window_size])
    depth_axis.set_xticklabels([])
    depth_axes.append(depth_axis)
    
    gene_axis = plt.Subplot(fig, gene_grid[2*idx+1])
    fig.add_subplot(gene_axis)
    gene_axis.set_ylabel('Freq')
    gene_axis.spines['top'].set_visible(False)
    gene_axis.spines['right'].set_visible(False)
    gene_axis.get_xaxis().tick_bottom()
    gene_axis.get_yaxis().tick_left()
    
    gene_axis.set_ylim([-0.05,1.05])
    gene_axis.set_yticks([0,0.5,1])
    gene_axis.set_xlim([desired_locations[idx][1]-window_size, desired_locations[idx][1]+window_size])
    gene_axis.set_xticklabels([])
    gene_axes.append(gene_axis)
    

gene_axes[1].set_xlabel('Position')

highlighted_gene_names = [[] for i in xrange(0,len(desired_sample_pairs))]
gene_names = [[] for i in xrange(0,len(desired_sample_pairs))]
initial_freqs = [[] for i in xrange(0,len(desired_sample_pairs))]
final_freqs = [[] for i in xrange(0,len(desired_sample_pairs))]

snp_changes = [[] for i in xrange(0,len(desired_sample_pairs))]


if debug==True: # or debug==False:
    spatial_gene_names = [['test'] for loc in desired_locations]
    spatial_chromosomes = [['test'] for loc in desired_locations]
    spatial_positions = [[loc[1]] for loc in desired_locations]
    spatial_initial_freqs = [[0.1] for loc in desired_locations]
    spatial_final_freqs = [[1] for loc in desired_locations]
    spatial_initial_depths = [[10] for loc in desired_locations]
    spatial_final_depths = [[10] for loc in desired_locations]
else:
    spatial_gene_names = [[] for loc in desired_locations]
    spatial_chromosomes = [[] for loc in desired_locations]
    spatial_positions = [[] for loc in desired_locations]
    spatial_initial_freqs = [[] for loc in desired_locations]
    spatial_final_freqs = [[] for loc in desired_locations]
    spatial_initial_depths = [[] for loc in desired_locations]
    spatial_final_depths = [[] for loc in desired_locations]


# Analyze SNPs, looping over chunk sizes. 
# Clunky, but necessary to limit memory usage on cluster

# Load SNP information for species_name
sys.stderr.write("Loading SNPs for %s...\n" % species_name)
sys.stderr.write("(not just core genes...)\n")
initial_freqs = [[] for i in xrange(0,len(desired_sample_pairs))]
final_freqs = [[] for i in xrange(0,len(desired_sample_pairs))]

# First get list of genes with >= 2 fixations
for pair_idx in xrange(0,len(desired_sample_pairs)):
    
    perr, mutations, reversions = calculate_temporal_changes.calculate_mutations_reversions_from_temporal_change_map(temporal_change_map, desired_sample_pairs[pair_idx][0], desired_sample_pairs[pair_idx][1])
    
    print pair_idx, "SNP perr =", perr
    if perr<-0.5:
        mutations = []
        reversions = []    
    
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
            
                print gene_name
            
                highlighted_gene_set.add(gene_name)
    
    for snp_change in (mutations+reversions):
        if snp_change[0] in highlighted_gene_set:
            print snp_change
                
    highlighted_gene_names[pair_idx] = highlighted_gene_set
    
final_line_number = 0
while final_line_number >= 0:
    
    sys.stderr.write("Loading chunk starting @ %d...\n" % final_line_number)
    dummy_samples, allele_counts_map, passed_sites_map, final_line_number = parse_midas_data.parse_snps(species_name, debug=debug, allowed_samples=snp_samples,chunk_size=chunk_size,initial_line_number=final_line_number)
    sys.stderr.write("Done! Loaded %d genes\n" % len(allele_counts_map.keys()))
    print len(dummy_samples), "dummy samples!"
    
    
    desired_sample_pair_idxs = []
    for initial_sample, final_sample in desired_sample_pairs:
        initial_idx = numpy.nonzero(dummy_samples==initial_sample)[0][0]
        final_idx = numpy.nonzero(dummy_samples==final_sample)[0][0]
    
        desired_sample_pair_idxs.append((initial_idx, final_idx))
        
    
    # Calculate fixation matrix
    sys.stderr.write("Calculating joint freqs...\n")
    for pair_idx in xrange(0,len(desired_sample_pairs)):
        
        initial_sample_idx, final_sample_idx = desired_sample_pair_idxs[pair_idx]
        
        chunk_gene_names, chunk_chromosomes, chunk_positions, chunk_initial_freqs, chunk_final_freqs, chunk_initial_depths, chunk_final_depths = diversity_utils.calculate_temporal_sample_freqs(allele_counts_map, passed_sites_map, initial_sample_idx, final_sample_idx)  # 
    
        joint_passed_sites = (chunk_initial_depths>0)*(chunk_final_depths>0)
    
        gene_names[pair_idx].extend(chunk_gene_names[joint_passed_sites])
        initial_freqs[pair_idx].extend(chunk_initial_freqs[joint_passed_sites])
        final_freqs[pair_idx].extend(chunk_final_freqs[joint_passed_sites])
        
        if pair_idx==2:
            for spatial_idx in xrange(0,len(desired_locations)):
                desired_chromosome, desired_position = desired_locations[spatial_idx]
            
                desired_sites = (chunk_chromosomes==desired_chromosome)*(numpy.fabs(chunk_positions-desired_position)<window_size)
            
                if desired_sites.sum() > 0:
                    
                    spatial_gene_names[spatial_idx].extend( chunk_gene_names[desired_sites] )
                    spatial_positions[spatial_idx].extend( chunk_positions[desired_sites])
                    
                    spatial_initial_freqs[spatial_idx].extend( chunk_initial_freqs[desired_sites] )
                    
                    spatial_final_freqs[spatial_idx].extend( chunk_final_freqs[desired_sites] )
                    
                    spatial_initial_depths[spatial_idx].extend( chunk_initial_depths[desired_sites])
                    
                    spatial_final_depths[spatial_idx].extend( chunk_final_depths[desired_sites])
                    
        
        
        snp_changes[pair_idx].extend( diversity_utils.calculate_snp_differences_between(initial_sample_idx, final_sample_idx, allele_counts_map, passed_sites_map,lower_threshold=snv_lower_threshold, 
upper_threshold=snv_upper_threshold) )
         
    sys.stderr.write("Done!\n")
    
for pair_idx in xrange(0,len(desired_sample_pairs)):
    gene_names[pair_idx] = numpy.array(gene_names[pair_idx])
    initial_freqs[pair_idx] = numpy.array(initial_freqs[pair_idx])
    final_freqs[pair_idx] = numpy.array(final_freqs[pair_idx])
    
    
sys.stderr.write("Done!\n")   

# Plot joint SFS!

for pair_idx in xrange(1,len(desired_sample_pairs)):

    initial_sample, final_sample = desired_sample_pairs[pair_idx]

    
    f0s = initial_freqs[pair_idx]
    f1s = final_freqs[pair_idx]
    
    dfs = (f1s-f0s)
    
    abs_dfs = numpy.fabs(dfs)
    
    major_freqs = numpy.fmax(f1s,1-f1s)
    
    # Don't plot a bunch of crap near 0
    good_sites = numpy.logical_or(f0s>0.05, f1s>0.05)
    
    #sfs_axes[pair_idx-1].plot(f0s[good_sites], f1s[good_sites],'.',alpha=0.1,markersize=2,markeredgewidth=0,color='0.7')
    
    sfs_axes[pair_idx-1].plot(f0s[good_sites], f1s[good_sites],'.',alpha=0.1,markersize=2,markeredgewidth=0,color='0.7')
    
    # Now plot ones in fixed gene set with colors
    for gene_name in sorted(highlighted_gene_names[pair_idx]):
        
        in_gene = (gene_names[pair_idx]==gene_name)*good_sites
         
        sfs_axes[pair_idx-1].plot(f0s[in_gene], f1s[in_gene],'.',alpha=0.75,markersize=4,markeredgewidth=0)
    
    #print initial_sample, final_sample
    #for snp_change in snp_changes[pair_idx]:
    #    print snp_change
    
    
    sfs_axes[pair_idx-1].set_xlim([0,1.05])
    sfs_axes[pair_idx-1].set_ylim([0,1.05])
    

# Plot spatial snp profiles near fixations
for spatial_idx in xrange(0,len(desired_locations)):
    
    # Sort snps by position on chromosome
    
    if len(spatial_positions[spatial_idx]) > 0:
        positions, gene_names, initial_freqs, final_freqs, initial_depths, final_depths = (numpy.array(x) for x in zip(*sorted(zip(spatial_positions[spatial_idx], spatial_gene_names[spatial_idx], spatial_initial_freqs[spatial_idx], spatial_final_freqs[spatial_idx],spatial_initial_depths[spatial_idx], spatial_final_depths[spatial_idx]))))
    else:
        positions, gene_names, initial_freqs, final_freqs, initial_depths, final_depths = (numpy.array([0]) for x in xrange(0,6))

    # Only looks at sites above a minimum freq threshold
    #good_sites = numpy.logical_or(initial_freqs>0.05, final_freqs>0.05)
    # Does nothing:
    #good_sites = numpy.logical_or(initial_freqs>-1, final_freqs>-1)
    good_sites = numpy.logical_and(initial_depths>0, final_depths>0)
    
    initial_hs = 4*initial_freqs*(1-initial_freqs)
    final_hs = 4*final_freqs*(1-final_freqs)
    
 
    for position, f0, ff, h0, hf in zip(positions[good_sites], initial_freqs[good_sites], final_freqs[good_sites], initial_hs[good_sites], final_hs[good_sites]):
    
        # Polarize such that everything is given by the major allele at end. 
        if ff<0.5:
            f0 = 1-f0
            ff = 1-ff 
        
        ff += 0.01
        
        if (ff>0.8 and f0<0.2) or (ff<0.2 and f0>0.8):
            color='r'
        else:
            color='0.7'
        
        gene_axes[spatial_idx].plot([position, position], [f0,ff],'-',color=color,linewidth=0.25)
        gene_axes[spatial_idx].plot([position],[ff],'.',color=color,markeredgewidth=0, markersize=1.5)
        
    # Plot depths
    for position, D0, Df in zip(positions, initial_depths, final_depths):
        
        depth_axes[spatial_idx].plot([position, position], [D0,Df],'b-',linewidth=0.25)
        depth_axes[spatial_idx].plot([position],[Df],'b.',markeredgewidth=0, markersize=1.5)
    
    

sys.stderr.write("Saving figure...\t")
fig.savefig('%s/figure_7%s.png' % (parse_midas_data.analysis_directory, other_species_str),bbox_inches='tight',dpi=600)
sys.stderr.write("Done!\n")

    
