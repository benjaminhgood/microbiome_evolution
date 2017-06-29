# within-host changes across species

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

################################################################################
#
# Standard header to read in argument information
#
################################################################################
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--debug", help="Loads only a subset of SNPs for speed", action="store_true")
parser.add_argument("--chunk-size", type=int, help="max number of records to load", default=1000000000)

args = parser.parse_args()

debug = args.debug
chunk_size = args.chunk_size
    
################################################################################

min_coverage = config.min_median_coverage
alpha = 0.5 # Confidence interval range for rate estimates
low_pi_threshold = 1e-03
clade_divergence_threshold = 1e-02
modification_divergence_threshold = 1e-03
min_change = 0.8


####################################################
#
# Set up Figure (4 panels, arranged in 2x2 grid)
#
####################################################

pylab.figure(1,figsize=(5,2))
fig = pylab.gcf()
# make three panels panels
outer_grid  = gridspec.GridSpec(2,2, height_ratios=[1,1], hspace=0.05, width_ratios=[10,1],wspace=0.15)

snp_axis = plt.Subplot(fig, outer_grid[0,0])
fig.add_subplot(snp_axis)

snp_axis.set_ylabel('Number of samples')
snp_axis.set_ylim([0,50])

snp_axis.spines['top'].set_visible(False)
snp_axis.spines['right'].set_visible(False)
snp_axis.get_xaxis().tick_bottom()
snp_axis.get_yaxis().tick_left()

gene_axis = plt.Subplot(fig, outer_grid[1,0])
fig.add_subplot(gene_axis)

gene_axis.set_ylabel('Number of samples')
gene_axis.set_ylim([0,50])

gene_axis.spines['top'].set_visible(False)
gene_axis.spines['right'].set_visible(False)
gene_axis.get_xaxis().tick_bottom()
gene_axis.get_yaxis().tick_left()

# Load subject and sample metadata
sys.stderr.write("Loading sample metadata...\n")
subject_sample_map = parse_HMP_data.parse_subject_sample_map()
sample_country_map = parse_HMP_data.parse_sample_country_map()
sample_order_map = parse_HMP_data.parse_sample_order_map()
sys.stderr.write("Done!\n")
       

sys.stderr.write("Saving figure...\t")
fig.savefig('%s/figure_7.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')
sys.stderr.write("Done!\n")
