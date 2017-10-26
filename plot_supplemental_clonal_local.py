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
from numpy.random import randint, choice,normal

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

snv_lower_threshold = 0.2
snv_upper_threshold = 0.8

# Load subject and sample metadata
sys.stderr.write("Loading sample metadata...\n")
subject_sample_map = parse_HMP_data.parse_subject_sample_map()
sample_country_map = parse_HMP_data.parse_sample_country_map()
sample_order_map = parse_HMP_data.parse_sample_order_map()
sys.stderr.write("Done!\n")

species_name = 'Bacteroides_vulgatus_57955'
desired_sample_pairs = {}
desired_sample_pairs['Bacteroides_vulgatus_57955'] = []
# Clonal sweep example
desired_sample_pairs['Bacteroides_vulgatus_57955'].append(('700110812', '700123827')) 
# Local sweep example
desired_sample_pairs['Bacteroides_vulgatus_57955'].append(('700014837', '700098561')) 
# Other locals weep example
desired_sample_pairs['Bacteroides_vulgatus_57955'].append(('700106170', '700163772')) 

num_examples = len(desired_sample_pairs[species_name])

#435590.9.peg.1205       Putative outer membrane protein, probably involved in nutrient binding
#435590.9.peg.2541       Thiol-activated cytolysin

# Locations of SNP changes to zoom in on
desired_diploid_genes = ['435590.9.peg.1205',  '435590.9.peg.2541'] 
desired_locations = [('NC_009614', 1567636L), ('NC_009614', 3163600)]

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

pylab.figure(1,figsize=(7,1.25))
fig = pylab.gcf()
# make three panels panels

outer_grid  = gridspec.GridSpec(1,2,width_ratios=[num_examples,0.6],wspace=0.2)

transition_grid = gridspec.GridSpecFromSubplotSpec(1,num_examples, width_ratios=[1 for i in xrange(0,num_examples)], wspace=0.05, subplot_spec=outer_grid[0])

gene_grid = gridspec.GridSpecFromSubplotSpec(num_examples-1,1, height_ratios=[1 for i in xrange(0,num_examples-1)], hspace=0.05, subplot_spec=outer_grid[1])

sfs_axes = []
for example_idx in xrange(0,num_examples):


    sfs_axis = plt.Subplot(fig, transition_grid[example_idx])
    fig.add_subplot(sfs_axis)
    sfs_axis.set_xlim([-0.05,1.05])
    sfs_axis.set_ylim([-0.05,1.05])
    line, = sfs_axis.plot([0,1],[0.8, 0.8],'-',color='0.7')
    line.set_dashes((1,1))
    line, = sfs_axis.plot([0.2,0.2],[0,1],'-',color='0.7')
    line.set_dashes((1,1))
    
    if example_idx==0:
        sfs_axis.set_ylabel('Final frequency')
    else:
        sfs_axis.set_yticklabels([])
    sfs_axis.set_xlabel('Initial frequency')
    
    
    #sfs_axis.plot([0,1],[0,1.0],'k-')

    sfs_axis.spines['top'].set_visible(False)
    sfs_axis.spines['right'].set_visible(False)
    sfs_axis.get_xaxis().tick_bottom()
    sfs_axis.get_yaxis().tick_left()
    sfs_axes.append(sfs_axis)

gene_axes = []
for example_idx in xrange(1, num_examples):
    
    gene_axis = plt.Subplot(fig, gene_grid[example_idx-1])
    fig.add_subplot(gene_axis)
    
    
    gene_axis.semilogx([1e-03],[1],'k.')
    gene_axis.set_xlim([5e-02,2])
    
    if example_idx==1:
        gene_axis.set_xticklabels([])    
    else:
        gene_axis.set_xlabel('Relative coverage')
        gene_axis.set_ylabel('                       Fraction of sites')
        
    gene_axis.spines['top'].set_visible(False)
    gene_axis.spines['right'].set_visible(False)
    gene_axis.get_xaxis().tick_bottom()
    gene_axis.get_yaxis().tick_left()
    
    gene_axis.set_xticklabels([])
    gene_axes.append(gene_axis)
    
highlighted_gene_names = [[] for i in xrange(0,len(desired_sample_pairs))]
gene_names = [[] for i in xrange(0,len(desired_sample_pairs))]
initial_freqs = [[] for i in xrange(0,len(desired_sample_pairs))]
final_freqs = [[] for i in xrange(0,len(desired_sample_pairs))]

snp_changes = [[] for i in xrange(0,len(desired_sample_pairs))]

# Analyze SNPs, looping over chunk sizes. 
# Clunky, but necessary to limit memory usage on cluster

# Load SNP information for species_name
sys.stderr.write("Loading SNPs for %s...\n" % species_name)
sys.stderr.write("(not just core genes...)\n")
initial_freqs = [[] for i in xrange(0,len(desired_sample_pairs))]
final_freqs = [[] for i in xrange(0,len(desired_sample_pairs))]
initial_depths = [[] for i in xrange(0,len(desired_sample_pairs))]

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
        initial_depths[pair_idx].extend(chunk_initial_depths[joint_passed_sites])
        final_freqs[pair_idx].extend(chunk_final_freqs[joint_passed_sites])
        
        snp_changes[pair_idx].extend( diversity_utils.calculate_snp_differences_between(initial_sample_idx, final_sample_idx, allele_counts_map, passed_sites_map,lower_threshold=snv_lower_threshold, 
upper_threshold=snv_upper_threshold) )
         
    sys.stderr.write("Done!\n")
    
for pair_idx in xrange(0,len(desired_sample_pairs)):
    gene_names[pair_idx] = numpy.array(gene_names[pair_idx])
    initial_freqs[pair_idx] = numpy.array(initial_freqs[pair_idx])
    final_freqs[pair_idx] = numpy.array(final_freqs[pair_idx])
    
    
sys.stderr.write("Done!\n")   

highlighted_gene_colors = {}
highlighted_gene_depths = {}

# Plot joint SFS!

for pair_idx in xrange(0,num_examples):

    initial_sample, final_sample = desired_sample_pairs[pair_idx]

    D0s = numpy.array(initial_depths[pair_idx])
    f0s = initial_freqs[pair_idx]
    f1s = final_freqs[pair_idx]
    
    dfs = (f1s-f0s)
    
    abs_dfs = numpy.fabs(dfs)
    
    major_freqs = numpy.fmax(f1s,1-f1s)
    
    # Don't plot a bunch of crap near 0
    good_sites = numpy.logical_or(f0s>0.05, f1s>0.05)
    
    sfs_axes[pair_idx].plot(f0s[good_sites], f1s[good_sites],'.',alpha=0.1,markersize=2,markeredgewidth=0,color='0.7')
    
    # Now plot ones in fixed gene set with colors
    highlighted_gene_colors[pair_idx] = {}
    highlighted_gene_depths[pair_idx] = {}
    
    for gene_name in sorted(highlighted_gene_names[pair_idx]):
        
        in_gene = (gene_names[pair_idx]==gene_name)*good_sites
         
        line, = sfs_axes[pair_idx].plot(f0s[in_gene], f1s[in_gene],'.',alpha=0.75,markersize=4,markeredgewidth=0)
        
        highlighted_gene_colors[pair_idx][gene_name] = pylab.getp(line,'color')
        if in_gene.sum() > 0:
            highlighted_gene_depths[pair_idx][gene_name] = numpy.median( D0s[in_gene] )
        else:
            # No such sites. Plot value off screen
            highlighted_gene_depths[pair_idx][gene_name] = 1000
        
    #print initial_sample, final_sample
    #for snp_change in snp_changes[pair_idx]:
    #    print snp_change
    
    
    sfs_axes[pair_idx].set_xlim([0,1.05])
    sfs_axes[pair_idx].set_ylim([0,1.05])
    

# Calculate depth distributions
import sfs_utils
samples, sfs_map = parse_midas_data.parse_within_sample_sfs(species_name,     allowed_variant_types=set(['1D','2D','3D','4D'])) 

for example_idx in xrange(1,num_examples):
    
    axis = gene_axes[example_idx-1]
    
    initial_sample = desired_sample_pairs[example_idx][0]


 
    bins, Ds, pDs =     sfs_utils.calculate_binned_depth_distribution_from_sfs_map(sfs_map[initial_sample])
    
    
    Dbar = (Ds*pDs).sum()/pDs.sum()
    
    bins = bins/Dbar
    
    ymax = pDs.max()*1.1
    
    haploid_color = '#08519c'
    
    #heights,bins,patches = axis.hist(gene_copynums,bins=bins,zorder=1,color=haploid_color,edgecolor=haploid_color)
    axis.bar(bins[0:-1], pDs,width=(bins[1:]-bins[0:-1]), zorder=1,color=haploid_color,edgecolor=haploid_color)
    
    axis.fill_between([0.01, 0.05],[0,0],[ymax,ymax],color='0.8',zorder=0)
    axis.fill_between([0.5, 2],[0,0],[ymax,ymax],color='0.8',zorder=0)
    axis.semilogx([2],[-1],'k.')
    axis.set_ylim([0,ymax])
    #axis.set_xlim([1e-02,4e00])   
    axis.set_xlim([5e-02,5])

    for gene_name in sorted(highlighted_gene_names[example_idx]):
        copynum = highlighted_gene_depths[example_idx][gene_name]/Dbar
        color = highlighted_gene_colors[example_idx][gene_name]
        height = ymax/2+(ymax/2)*0.1*normal(0,1)
        axis.plot([copynum, copynum],[0,height],'-',color=color)
        axis.plot([copynum],[height],'.',color=color,markersize=3)     

sys.stderr.write("Saving figure...\t")
fig.savefig('%s/supplemental_clonal_local_sweeps%s.png' % (parse_midas_data.analysis_directory, other_species_str),bbox_inches='tight',dpi=600)
sys.stderr.write("Done!\n")

    
