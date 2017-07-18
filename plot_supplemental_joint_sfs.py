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

# Load subject and sample metadata
sys.stderr.write("Loading sample metadata...\n")
subject_sample_map = parse_HMP_data.parse_subject_sample_map()
sample_country_map = parse_HMP_data.parse_sample_country_map()
sample_order_map = parse_HMP_data.parse_sample_order_map()
sys.stderr.write("Done!\n")
  
  
desired_sample_pairs = {}
desired_sample_pairs['Bacteroides_vulgatus_57955'] = []
#desired_sample_pairs['Bacteroides_vulgatus_57955'].append(('700016765', '700100312')) # haploid->haploid example
#desired_sample_pairs['Bacteroides_vulgatus_57955'].append(('700016960', '700101581')) # haploid->diploid example
#desired_sample_pairs['Bacteroides_vulgatus_57955'].append(('700014837', '700098561')) # diploid->diploid example
desired_sample_pairs['Bacteroides_vulgatus_57955'].append(('700024318', '700105372')) # diploid->diploid example
desired_sample_pairs['Bacteroides_vulgatus_57955'].append(('700016960', '700101581')) # diploid->diploid example
desired_sample_pairs['Bacteroides_vulgatus_57955'].append(('700014837', '700098561')) # diploid->diploid example


desired_sample_pairs['Eubacterium_rectale_56927'] = []
desired_sample_pairs['Eubacterium_rectale_56927'].append(('700098429', '700110812')) # haploid->haploid example
desired_sample_pairs['Eubacterium_rectale_56927'].append(('700098429', '700110812')) # haploid->haploid example 2
desired_sample_pairs['Eubacterium_rectale_56927'].append(('700021306', '700037632')) # diploid->haploid example

desired_sample_pairs['Bacteroides_uniformis_57318'] = []
desired_sample_pairs['Bacteroides_uniformis_57318'].append(('700021902', '700105210')) # haploid->haploid example
desired_sample_pairs['Bacteroides_uniformis_57318'].append(('700023113', '700023720')) # diploid->haploid example
desired_sample_pairs['Bacteroides_uniformis_57318'].append(('700015415', '700101243')) # diploid->diploid example

desired_sample_pairs['Bacteroides_ovatus_58035'] = []
desired_sample_pairs['Bacteroides_ovatus_58035'].append(('700037042', '700103621')) # haploid->haploid example
desired_sample_pairs['Bacteroides_ovatus_58035'].append(('700023578', '700103446')) # diploid->diploid example1
desired_sample_pairs['Bacteroides_ovatus_58035'].append(('700024615', '700105306')) # diploid->diploid example2


desired_sample_pairs = desired_sample_pairs[species_name]
snp_samples = set()
for initial_sample, final_sample in desired_sample_pairs:
    snp_samples.add(initial_sample)
    snp_samples.add(final_sample)
    
snp_samples = list(snp_samples)
       
####################################################
#
# Set up Figure (3 panels, arranged in 1x3 grid)
#
####################################################

pylab.figure(1,figsize=(3.42,1.5))
fig = pylab.gcf()
# make three panels panels
outer_grid  = gridspec.GridSpec(1,2, width_ratios=[1,1], wspace=0.05)

sfs_axes = []

sfs_axis_1 = plt.Subplot(fig, outer_grid[0])
fig.add_subplot(sfs_axis_1)
sfs_axis_1.set_ylabel('Final frequency')
sfs_axis_1.set_xlabel('Initial frequency')
sfs_axis_1.set_xlim([-0.05,0.51])
sfs_axis_1.set_ylim([-0.05,1.05])
sfs_axis_1.plot([0,1],[1,1],'k-')
sfs_axis_1.plot([0,0],[0,1],'k-')
sfs_axis_1.plot([0,0.2],[0.8,1.0],'k-')
sfs_axis_1.plot([0,1],[0,1.0],'k-')

sfs_axis_1.spines['top'].set_visible(False)
sfs_axis_1.spines['right'].set_visible(False)
sfs_axis_1.get_xaxis().tick_bottom()
sfs_axis_1.get_yaxis().tick_left()
sfs_axes.append(sfs_axis_1)

sfs_axis_2 = plt.Subplot(fig, outer_grid[1])
fig.add_subplot(sfs_axis_2)
sfs_axis_2.set_xlabel('Initial frequency')
sfs_axis_2.set_xlim([-0.05,0.51])
sfs_axis_2.set_ylim([-0.05,1.05])
sfs_axis_2.plot([0,1],[1,1],'k-')
sfs_axis_2.plot([0,0],[0,1],'k-')
sfs_axis_2.plot([0,0.2],[0.8,1.0],'k-')
sfs_axis_2.plot([0,1],[0,1.0],'k-')

sfs_axis_2.set_yticklabels([])
sfs_axis_2.spines['top'].set_visible(False)
sfs_axis_2.spines['right'].set_visible(False)
sfs_axis_2.get_xaxis().tick_bottom()
sfs_axis_2.get_yaxis().tick_left()
sfs_axes.append(sfs_axis_2)

#sfs_axis_3 = plt.Subplot(fig, outer_grid[2])
#fig.add_subplot(sfs_axis_3)
#sfs_axis_3.set_xlabel('Initial frequency')
#sfs_axis_3.set_xlim([-0.05,0.51])
#sfs_axis_3.set_ylim([-0.05,1.05])
#sfs_axis_3.plot([0,1],[1,1],'k-')
#sfs_axis_3.plot([0,0],[0,1],'k-')
#sfs_axis_3.plot([0,0.2],[0.8,1.0],'k-')
#sfs_axis_3.plot([0,1],[0,1.0],'k-')

#sfs_axis_3.set_yticklabels([])
#sfs_axis_3.spines['top'].set_visible(False)
#sfs_axis_3.spines['right'].set_visible(False)
#sfs_axis_3.get_xaxis().tick_bottom()
#sfs_axis_3.get_yaxis().tick_left()
#sfs_axes.append(sfs_axis_3)

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
        
        chunk_initial_freqs, chunk_final_freqs = diversity_utils.calculate_temporal_sample_freqs(allele_counts_map, passed_sites_map, initial_sample_idx, final_sample_idx)  # 
    
        initial_freqs[pair_idx].extend(chunk_initial_freqs)
        final_freqs[pair_idx].extend(chunk_final_freqs)
        
        
        snp_changes[pair_idx].extend( diversity_utils.calculate_snp_differences_between(initial_sample_idx, final_sample_idx, allele_counts_map, passed_sites_map) )
         
    sys.stderr.write("Done!\n")
    
for pair_idx in xrange(0,len(desired_sample_pairs)):
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
    
    good_sites = numpy.logical_or(f0s>0.05, f1s>0.05)
    fixed_sites = good_sites*(dfs>0.8)
    
    print initial_sample, final_sample
    for snp_change in snp_changes[pair_idx]:
        print snp_change
    #sfs_axes[pair_idx].plot(f0s[good_sites],f1s[good_sites],'b.',alpha=0.1,markersize=2,markeredgewidth=0)
    #sfs_axes[pair_idx].plot(f0s[fixed_sites],f1s[fixed_sites],'b.',alpha=1,markersize=3,markeredgewidth=0)
    
    sfs_axes[pair_idx-1].plot(f0s[good_sites],f1s[good_sites],'b.',alpha=0.1,markersize=2,markeredgewidth=0)
    sfs_axes[pair_idx-1].plot(f0s[fixed_sites],f1s[fixed_sites],'b.',alpha=1,markersize=3,markeredgewidth=0)
    sfs_axes[pair_idx-1].set_xlim([0,1])
    sfs_axes[pair_idx-1].set_ylim([0,1])
    

sys.stderr.write("Saving figure...\t")
fig.savefig('%s/supplemental_joint_sfs%s.png' % (parse_midas_data.analysis_directory, other_species_str),bbox_inches='tight',dpi=600)
sys.stderr.write("Done!\n")

    
