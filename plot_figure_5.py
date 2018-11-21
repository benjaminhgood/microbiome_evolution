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
    
    
import matplotlib.colors as colors
import matplotlib.cm as cmx
from math import log10,ceil,log,exp
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy.random import randint, random, choice, multinomial, shuffle
import matplotlib.colors as mcolors
import matplotlib.patheffects as pe

from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster.hierarchy import cophenet
from scipy.cluster.hierarchy import fcluster


mpl.rcParams['font.size'] = 7
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
parser.add_argument("--modification-threshold", type=int, help="max number of SNV differences before calling a modification", default=config.modification_difference_threshold)


args = parser.parse_args()

debug = args.debug
chunk_size = args.chunk_size
memoize = args.memoize
modification_difference_threshold = args.modification_threshold
replacement_difference_threshold = config.replacement_difference_threshold
twin_modification_difference_threshold = config.twin_modification_difference_threshold
twin_replacement_difference_threshold = config.twin_replacement_difference_threshold
default_num_bootstraps = 10000
################################################################################

#####################
#
# Settings for calculation:
#
#####################

output_filename = ('%s/within_host_output.txt' % (parse_midas_data.analysis_directory))
output_strs = []

min_coverage = config.min_median_coverage
min_sample_size = 3
min_haploid_sample_size = 10

variant_types = ['1D','4D']

within_host_type = 'consecutive' # consecutive timepoints
#within_host_type = 'longest' # longest interval available

# For partitioning SNVs according to prevalence
derived_freq_bins = numpy.array([-1,0,0.01,0.1,0.5,0.9,0.99,1,2])
derived_virtual_freqs = numpy.arange(0,len(derived_freq_bins)-1)
derived_virtual_xticks = list(derived_virtual_freqs[:-1]+0.5)
derived_virtual_xticklabels = ['0','.01','.1','.5','.9','.99','1']

# For partitioning genes into different prevalence classes
gene_freq_bins = numpy.array([-1,0.1,0.5,0.9,2])
gene_freq_xticks      = [-4, -3,  -2,   -1,   0,   1,    2,   3, 4]
gene_freq_xticklabels = ['0','0.1','0.5', '0.9','1','0.9','0.5', '0.1','0']
gene_gain_virtual_freqs = numpy.array([3.5,2.5,1.5,0.5])
gene_loss_virtual_freqs = numpy.array([-3.5,-2.5,-1.5,-0.5])

#gene_freq_bins = numpy.array([-1,0.1,0.9,2])
#gene_freq_xticks      = [-3,  -2,   -1,   0,   1,    2,   3]
#gene_freq_xticklabels = ['0','0.1', '0.9','1','0.9', '0.1','0']
#gene_gain_virtual_freqs = numpy.array([2.5,1.5,0.5])
#gene_loss_virtual_freqs = numpy.array([-2.5,-1.5,-0.5])


#####
#
# Settings for different cohorts we are looking at 
#
#####
cohorts = ["hmp", "twins", "young_twins"]
countries = ["United States", "United Kingdom", "Western Europe"]
country_cohort_map = {country: cohort for country,cohort in zip(countries,cohorts)}

modification_difference_thresholds = {"hmp": modification_difference_threshold, "twins": 1e06, "young_twins": twin_modification_difference_threshold}

replacement_difference_thresholds = {"hmp": replacement_difference_threshold, "twins": twin_replacement_difference_threshold, "young_twins": twin_replacement_difference_threshold}

################################
#
# Set up figures
#
################################

####################################################
#
# Distribution of changes across individual species
#
####################################################

pylab.figure(1,figsize=(7,5))
fig = pylab.gcf()
# make three panels panels
outer_grid = gridspec.GridSpec(2,1,height_ratios=[1,1],hspace=0.25)

upper_grid  = gridspec.GridSpecFromSubplotSpec(1,2,width_ratios=[1,1],wspace=0.1, subplot_spec=outer_grid[0])

species_snp_axis = plt.Subplot(fig, upper_grid[0])
fig.add_subplot(species_snp_axis)

species_snp_axis.spines['top'].set_visible(False)
species_snp_axis.spines['right'].set_visible(False)
species_snp_axis.get_xaxis().tick_bottom()
species_snp_axis.get_yaxis().tick_left()

species_snp_axis.set_ylabel('Fraction comparisons $\geq n$')
species_snp_axis.set_xlabel('# SNV changes')


species_gene_axis = plt.Subplot(fig, upper_grid[1])
fig.add_subplot(species_gene_axis)

species_gene_axis.spines['top'].set_visible(False)
species_gene_axis.spines['right'].set_visible(False)
species_gene_axis.get_xaxis().tick_bottom()
species_gene_axis.get_yaxis().tick_left()

species_snp_axis.loglog([0.01],[1],'k.')
species_gene_axis.loglog([0.01],[1],'k.')

species_snp_axis.set_xlim([3e-01,1e05])
species_gene_axis.set_xlim([3e-01,1e04])

species_gene_axis.set_xlabel('# gene changes')

species_legend_axis = plt.Subplot(fig, outer_grid[1])
fig.add_subplot(species_legend_axis)

species_legend_axis.set_ylim([0,1])
species_legend_axis.set_xlim([0,1])

species_legend_axis.spines['top'].set_visible(False)
species_legend_axis.spines['right'].set_visible(False)
species_legend_axis.spines['left'].set_visible(False)
species_legend_axis.spines['bottom'].set_visible(False)

species_legend_axis.set_xticks([])
species_legend_axis.set_yticks([])

##############
#
# Main text figure
#
###############
pylab.figure(2,figsize=(7,4))
fig2 = pylab.gcf()
# make three panels panels
outer_outer_grid_2 = gridspec.GridSpec(1,1) #2, width_ratios=[1,0.2],wspace=0.1) 
outer_grid_2 = gridspec.GridSpecFromSubplotSpec(2,1, height_ratios=[1,0.7],hspace=0.6, subplot_spec=outer_outer_grid_2[0])

prevalence_outer_grid = gridspec.GridSpecFromSubplotSpec(1,2, width_ratios=[1,0.1],wspace=0.2,subplot_spec=outer_grid_2[1])

prevalence_grid = gridspec.GridSpecFromSubplotSpec(1,2, width_ratios=[1,1],wspace=0.5,subplot_spec=prevalence_outer_grid[0])

pooled_grid = gridspec.GridSpecFromSubplotSpec(1,3,width_ratios=[1,1,0.2],wspace=0.15,subplot_spec=outer_grid_2[0])

pooled_snp_axis = plt.Subplot(fig2, pooled_grid[0])
fig2.add_subplot(pooled_snp_axis)
pooled_snp_axis.set_ylabel('Fraction comparisons $\geq n$')
pooled_snp_axis.set_xlabel('# SNV changes')
#pooled_axis.set_ylim([-35,35])
#pooled_snp_axis.set_xlim([2e-01,1e05])
pooled_snp_axis.set_xlim([0.6,1e05])

pooled_snp_axis.set_xticklabels([])

pooled_snp_axis.spines['top'].set_visible(False)
pooled_snp_axis.spines['right'].set_visible(False)
pooled_snp_axis.get_xaxis().tick_bottom()
pooled_snp_axis.get_yaxis().tick_left()
 
pooled_gene_axis = plt.Subplot(fig2, pooled_grid[1])
fig2.add_subplot(pooled_gene_axis)
#pooled_gene_axis.set_ylabel('Number of samples')
pooled_gene_axis.set_xlabel('# gene changes')
#pooled_axis.set_ylim([-35,35])
#pooled_gene_axis.set_xlim([2e-01,1e04])
pooled_gene_axis.set_xlim([0.6,1e04])

pooled_gene_axis.spines['top'].set_visible(False)
pooled_gene_axis.spines['right'].set_visible(False)
pooled_gene_axis.get_xaxis().tick_bottom()
pooled_gene_axis.get_yaxis().tick_left()

pooled_snp_axis.loglog([0.1],[1],'k.')
pooled_gene_axis.loglog([0.1],[1],'k.')
 
legend2_axis = plt.Subplot(fig2, pooled_grid[2])
fig2.add_subplot(legend2_axis)

legend2_axis.set_ylim([0,1])
legend2_axis.set_xlim([0,1])

legend2_axis.spines['top'].set_visible(False)
legend2_axis.spines['right'].set_visible(False)
legend2_axis.spines['left'].set_visible(False)
legend2_axis.spines['bottom'].set_visible(False)

legend2_axis.set_xticks([])
legend2_axis.set_yticks([])

legend2_axis.plot([-2,-1],[-1,-1],'-',linewidth=1, color='#08519c',label='Within-host')
legend2_axis.plot([-2,-1],[-1,-1],'-',color='#08519c',linewidth=1, label='modification',zorder=2,path_effects=[pe.Stroke(linewidth=5, foreground='#9ecae1'), pe.Normal()])
legend2_axis.plot([-2,-1],[-1,-1],'-',color='#08519c', label='replacement',linewidth=1,path_effects=[pe.Stroke(linewidth=5, foreground='#fee0d2'), pe.Normal()])
legend2_axis.plot([-2,-1],[-1,-1],'-',color='#08519c', label='no SNVs',linewidth=1,path_effects=[pe.Stroke(linewidth=5, foreground='0.8'), pe.Normal()])
legend2_axis.plot([-2,-1],[-1,-1], '-',linewidth=1,color='w', alpha=0.5, label=' ')
legend2_axis.plot([-2,-1],[-1,-1], '-',linewidth=1,color='r', alpha=0.5, label='Between-host\n(unrelated)')
legend2_axis.plot([-2,-1],[-1,-1],'-',linewidth=1,color='#8856a7', label='Between-host\n(adult twins)')

legend2_axis.legend(loc='upper center',frameon=False,fontsize=5,numpoints=1,ncol=1,handlelength=1)   


hmp_frequency_axis = plt.Subplot(fig2, prevalence_grid[0])
fig2.add_subplot(hmp_frequency_axis)

hmp_frequency_axis.spines['top'].set_visible(False)
hmp_frequency_axis.spines['right'].set_visible(False)
hmp_frequency_axis.get_xaxis().tick_bottom()
hmp_frequency_axis.get_yaxis().tick_left()
 
hmp_frequency_axis.set_xlabel('Derived allele prevalence\nacross hosts')
hmp_frequency_axis.set_ylabel('# SNV changes')

hmp_frequency_axis.set_xticks(derived_virtual_xticks)
hmp_frequency_axis.set_xticklabels(derived_virtual_xticklabels) #,rotation='vertical')

hmp_frequency_axis.set_ylim([0,200])

hmp_gene_frequency_axis = plt.Subplot(fig2, prevalence_grid[1])
fig2.add_subplot(hmp_gene_frequency_axis)

hmp_gene_frequency_axis.spines['top'].set_visible(False)
hmp_gene_frequency_axis.spines['right'].set_visible(False)
hmp_gene_frequency_axis.get_xaxis().tick_bottom()
hmp_gene_frequency_axis.get_yaxis().tick_left()
 
hmp_gene_frequency_axis.set_xlabel('Gene prevalence across hosts')
hmp_gene_frequency_axis.set_ylabel('# gene changes')

hmp_gene_frequency_axis.set_xlim([gene_freq_xticks[0],gene_freq_xticks[-1]])
hmp_gene_frequency_axis.set_xticks(gene_freq_xticks)
hmp_gene_frequency_axis.set_xticklabels(gene_freq_xticklabels) #,rotation='vertical')

hmp_gene_frequency_axis.plot([0,0],[100,100],'k-')
hmp_gene_frequency_axis.set_ylim([0,60])

hmp_gene_legend_axis = plt.Subplot(fig2, prevalence_outer_grid[1])
fig2.add_subplot(hmp_gene_legend_axis)

hmp_gene_legend_axis.set_ylim([0,1])
hmp_gene_legend_axis.set_xlim([0,1])

hmp_gene_legend_axis.spines['top'].set_visible(False)
hmp_gene_legend_axis.spines['right'].set_visible(False)
hmp_gene_legend_axis.spines['left'].set_visible(False)
hmp_gene_legend_axis.spines['bottom'].set_visible(False)

hmp_gene_legend_axis.set_xticks([])
hmp_gene_legend_axis.set_yticks([])

hmp_gene_legend_axis.bar([-2],[-1],width=0.2, linewidth=0,facecolor='#b3de69',label='gain')
hmp_gene_legend_axis.bar([-2],[-1],width=0.2, linewidth=0,facecolor='#ff7f00',label='loss')
hmp_gene_legend_axis.bar([-2],[-1],width=0.2, linewidth=0,facecolor='0.7',label='de novo\nexpectation')

hmp_gene_legend_axis.legend(loc='center left',frameon=False,fontsize=5,numpoints=1,ncol=1,handlelength=1)   

##############
#
# Young twin version of SNV and gene changes
#
###############

pylab.figure(3,figsize=(7,2))
fig3 = pylab.gcf()
outer_grid_3 = gridspec.GridSpec(1,2,width_ratios=[1,1],wspace=0.15)

young_snp_axis = plt.Subplot(fig3, outer_grid_3[0])
fig3.add_subplot(young_snp_axis)

young_snp_axis.set_ylabel('Fraction comparisons $\geq n$')
young_snp_axis.set_xlabel('# SNV changes')
#young_snp_axis.set_xlim([1,1e05])
young_snp_axis.set_xlim([0.6,1e05])
young_snp_axis.spines['top'].set_visible(False)
young_snp_axis.spines['right'].set_visible(False)
young_snp_axis.get_xaxis().tick_bottom()
young_snp_axis.get_yaxis().tick_left()

young_gene_axis = plt.Subplot(fig3, outer_grid_3[1])
fig3.add_subplot(young_gene_axis)

young_gene_axis.set_xlim([0.6,1e04])
young_gene_axis.set_xlabel('# gene changes')

young_gene_axis.spines['top'].set_visible(False)
young_gene_axis.spines['right'].set_visible(False)
young_gene_axis.get_xaxis().tick_bottom()
young_gene_axis.get_yaxis().tick_left()

young_snp_axis.loglog([0.1],[1],'k.')
young_gene_axis.loglog([0.1],[1],'k.')


##############
#
# Avg # of differences within hosts
#
###############

pylab.figure(4,figsize=(2,2))
fig4 = pylab.gcf()
outer_grid_4 = gridspec.GridSpec(1,1)

avg_axis = plt.Subplot(fig4, outer_grid_4[0])
fig4.add_subplot(avg_axis)

#avg_axis.set_ylabel('Microbiome-wide avg')
avg_axis.set_xlim([0.5,2.5])
avg_axis.set_xticks([])

avg_axis.spines['top'].set_visible(False)
avg_axis.spines['right'].set_visible(False)
avg_axis.get_xaxis().tick_bottom()
avg_axis.get_yaxis().tick_left()

##############
#
# Twin version of SNV & gene prevalence distribution
#
###############
pylab.figure(5,figsize=(7,2))
fig5 = pylab.gcf()
# make three panels panels
twin_prevalence_grid = gridspec.GridSpec(1,3, width_ratios=[1,1,0.3],wspace=0.5)

twin_frequency_axis = plt.Subplot(fig5, twin_prevalence_grid[0])
fig5.add_subplot(twin_frequency_axis)

twin_frequency_axis.spines['top'].set_visible(False)
twin_frequency_axis.spines['right'].set_visible(False)
twin_frequency_axis.get_xaxis().tick_bottom()
twin_frequency_axis.get_yaxis().tick_left()
 
twin_frequency_axis.set_xlabel('Derived allele prevalence\nacross hosts')
twin_frequency_axis.set_ylabel('# SNV changes')

twin_frequency_axis.set_xticks(derived_virtual_xticks)
twin_frequency_axis.set_xticklabels(derived_virtual_xticklabels) #,rotation='vertical')

#twin_frequency_axis.set_ylim([0,200])

twin_gene_frequency_axis = plt.Subplot(fig5, twin_prevalence_grid[1])
fig5.add_subplot(twin_gene_frequency_axis)

twin_gene_frequency_axis.spines['top'].set_visible(False)
twin_gene_frequency_axis.spines['right'].set_visible(False)
twin_gene_frequency_axis.get_xaxis().tick_bottom()
twin_gene_frequency_axis.get_yaxis().tick_left()
 
twin_gene_frequency_axis.set_xlabel('Gene prevalence across hosts')
twin_gene_frequency_axis.set_ylabel('# gene changes')

twin_gene_frequency_axis.set_xlim([gene_freq_xticks[0],gene_freq_xticks[-1]])
twin_gene_frequency_axis.set_xticks(gene_freq_xticks)
twin_gene_frequency_axis.set_xticklabels(gene_freq_xticklabels) #,rotation='vertical')

twin_gene_frequency_axis.plot([0,0],[100,100],'k-')
#twin_gene_frequency_axis.set_ylim([0,100])

twin_gene_legend_axis = plt.Subplot(fig5, twin_prevalence_grid[2])
fig5.add_subplot(twin_gene_legend_axis)

twin_gene_legend_axis.set_ylim([0,1])
twin_gene_legend_axis.set_xlim([0,1])

twin_gene_legend_axis.spines['top'].set_visible(False)
twin_gene_legend_axis.spines['right'].set_visible(False)
twin_gene_legend_axis.spines['left'].set_visible(False)
twin_gene_legend_axis.spines['bottom'].set_visible(False)

twin_gene_legend_axis.set_xticks([])
twin_gene_legend_axis.set_yticks([])

twin_gene_legend_axis.bar([-2],[-1],width=0.2, linewidth=0,facecolor='#b3de69',label='gain')
twin_gene_legend_axis.bar([-2],[-1],width=0.2, linewidth=0,facecolor='#ff7f00',label='loss')
twin_gene_legend_axis.bar([-2],[-1],width=0.2, linewidth=0,facecolor='0.7',label='de novo\nexpectation')

twin_gene_legend_axis.legend(loc='center left',frameon=False,fontsize=5,numpoints=1,ncol=1,handlelength=1)   

##############
#
# Young twin version of SNV & gene prevalence distribution
#
###############
pylab.figure(7,figsize=(7,2))
fig7 = pylab.gcf()
# make three panels panels
young_twin_prevalence_grid = gridspec.GridSpec(1,3, width_ratios=[1,1,0.3],wspace=0.5)

young_twin_frequency_axis = plt.Subplot(fig7, young_twin_prevalence_grid[0])
fig7.add_subplot(young_twin_frequency_axis)

young_twin_gene_frequency_axis = plt.Subplot(fig7, young_twin_prevalence_grid[1])
fig7.add_subplot(young_twin_gene_frequency_axis)

young_twin_gene_legend_axis = plt.Subplot(fig7, young_twin_prevalence_grid[2])
fig7.add_subplot(young_twin_gene_legend_axis)

####################################################
#
# Set up Suppplemental Fig (temporal haploid classification)
#
####################################################
# This figure spreads them all out

pylab.figure(6,figsize=(5,6))
fig6 = pylab.gcf()
# make three panels panels
outer_grid6  = gridspec.GridSpec(1,2,width_ratios=[1,1],wspace=0.2)

hmp_haploid_axis = plt.Subplot(fig6, outer_grid6[0])
fig6.add_subplot(hmp_haploid_axis)
hmp_haploid_axis.set_xlabel('HMP timepoint pairs')

twin_haploid_axis = plt.Subplot(fig6, outer_grid6[1])
fig6.add_subplot(twin_haploid_axis)
twin_haploid_axis.set_xlabel('Twin pairs')

####################################################
#
# Supplemental figure w/ colorbars
#
####################################################

cmap_str = 'YlGnBu'
vmin = -2
vmax = 3
cmap = pylab.get_cmap(cmap_str) 

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=-1):
    if n == -1:
        n = cmap.N
    new_cmap = mcolors.LinearSegmentedColormap.from_list(
         'trunc({name},{a:.2f},{b:.2f})'.format(name=cmap.name, a=minval, b=maxval),
         cmap(numpy.linspace(minval, maxval, n)))
    return new_cmap

cmap = truncate_colormap(cmap, 0.25, 1.0)
cNorm  = colors.Normalize(vmin=0, vmax=vmax)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cmap)


pylab.figure(8,figsize=(5,4))
fig8 = pylab.gcf()
# make three panels panels
outer_grid8  = gridspec.GridSpec(1,2,width_ratios=[50,1],wspace=0.05)

change_axis = plt.Subplot(fig8, outer_grid8[0])
fig8.add_subplot(change_axis)

change_axis.set_ylabel('HMP timepoint pairs')

#change_axis.spines['top'].set_visible(False)
#change_axis.spines['right'].set_visible(False)
change_axis.get_xaxis().tick_bottom()
change_axis.get_yaxis().tick_left()

cax = plt.Subplot(fig8, outer_grid8[1])
fig8.add_subplot(cax)


##############
#
# Full distribution of SNV and gene prevalences
#
###############
pylab.figure(9,figsize=(7,2))
fig9 = pylab.gcf()
# make three panels panels
ks_prevalence_grid = gridspec.GridSpec(1,3, width_ratios=[1,1,1],wspace=0.1)

snv_ks_axis = plt.Subplot(fig9, ks_prevalence_grid[0])
fig9.add_subplot(snv_ks_axis)
 
snv_ks_axis.set_xlabel('Derived allele prevalence\nacross hosts, $p$')
snv_ks_axis.set_ylabel('Fraction changes $\geq p$')
snv_ks_axis.set_xlim([-0.05,1.05])
snv_ks_axis.set_ylim([0,1])

dnds_ks_axis = plt.Subplot(fig9, ks_prevalence_grid[1])
fig9.add_subplot(dnds_ks_axis)
dnds_ks_axis.set_xlim([-0.05,1.05])
dnds_ks_axis.set_ylim([0,1])
dnds_ks_axis.set_yticklabels([])

dnds_ks_axis.set_xlabel('Derived allele prevalence\nacross hosts, $p$')

gene_ks_axis = plt.Subplot(fig9, ks_prevalence_grid[2])
fig9.add_subplot(gene_ks_axis)
gene_ks_axis.set_xlabel('Gene prevalence\nacross hosts, $p$')
gene_ks_axis.set_xlim([-0.05,1.05])
gene_ks_axis.set_ylim([0,1])
gene_ks_axis.set_yticklabels([])

##############
#
# Number of sites and genes retained
#
###############
pylab.figure(10,figsize=(6,2))
fig10 = pylab.gcf()
# make three panels panels
genome_length_grid = gridspec.GridSpec(1,2, width_ratios=[1,1], wspace=0.1)

genome_length_axis = plt.Subplot(fig10, genome_length_grid[0])
fig10.add_subplot(genome_length_axis)

genome_length_axis.set_xlabel('Number of sites compared')
genome_length_axis.set_ylabel('Fraction timepoint pairs $\geq n$')

genome_length_axis.set_ylim([0,1.1])

genome_length_axis.spines['top'].set_visible(False)
genome_length_axis.spines['right'].set_visible(False)
genome_length_axis.get_xaxis().tick_bottom()
genome_length_axis.get_yaxis().tick_left()


pangenome_length_axis = plt.Subplot(fig10, genome_length_grid[1])
fig10.add_subplot(pangenome_length_axis)

pangenome_length_axis.set_xlabel('Number of genes compared')

pangenome_length_axis.set_ylim([0,1.1])
pangenome_length_axis.set_yticklabels([])

pangenome_length_axis.spines['top'].set_visible(False)
pangenome_length_axis.spines['right'].set_visible(False)
pangenome_length_axis.get_xaxis().tick_bottom()
pangenome_length_axis.get_yaxis().tick_left()


################################
#
# Now do calculation
#
################################

hmp_species_qp_counts = {}
twin_species_qp_counts = {}

species_snp_change_distribution = {cohort: {} for cohort in cohorts}
species_snp_nerrs = {cohort: {} for cohort in cohorts}
species_gene_change_distribution = {cohort: {} for cohort in cohorts}
species_gene_nerrs = {cohort: {} for cohort in cohorts}

# observed within host value
pooled_snp_change_distribution = {cohort : [] for cohort in cohorts}
pooled_gene_change_distribution = {cohort : [] for cohort in cohorts} # for modifications

pooled_snp_length_distribution = {cohort : [] for cohort in cohorts}
pooled_gene_length_distribution = {cohort : [] for cohort in cohorts}

# typical value, median other sample
pooled_between_snp_change_distribution = {cohort : [] for cohort in cohorts}
pooled_between_gene_change_distribution = {cohort : [] for cohort in cohorts}

# closest other sample
pooled_min_between_snp_change_distribution = {cohort : [] for cohort in cohorts}
pooled_min_between_gene_change_distribution = {cohort : [] for cohort in cohorts}

replacement_map = {cohort: {} for cohort in cohorts}

total_freq_snps = {cohort: {} for cohort in cohorts}
total_null_freq_snps = {cohort: {} for cohort in cohorts}
for cohort in cohorts:
    
    # This holds the # of SNVs in each prevalence class of each var type
    total_freq_snps[cohort] = {var_type: numpy.zeros_like(derived_virtual_freqs) for var_type in variant_types}
    
    # This holds the expectation of the # of SNVs in each prevalence class of 1D and 4D
    # (i.e., relative fraction of opportunities for different species), conditioned on prevalence
    total_null_freq_snps[cohort] = {var_type: numpy.zeros_like(derived_virtual_freqs) for var_type in variant_types}


# This sums up across the different var_type categories
total_freq_all_snps = {cohort: numpy.zeros_like(derived_virtual_freqs) for cohort in cohorts}

# This is the null distribution of prevalence (conditioned on total # of SNVs)
total_null_freq_all_snps = {cohort: numpy.zeros_like(derived_virtual_freqs)*1.0 for cohort in cohorts}

total_freq_gains = {cohort: numpy.zeros(len(gene_freq_bins)-1)*1.0 for cohort in cohorts}
total_freq_losses = {cohort: numpy.zeros_like(total_freq_gains[cohort]) for cohort in cohorts}
total_null_freq_losses = {cohort: numpy.zeros_like(total_freq_gains[cohort]) for cohort in cohorts}


# SNV and gene prevalences resolved by subject
# so that we can do bootstrapping
snv_prevalence_count_map = {cohort: {} for cohort in cohorts}
gene_gain_count_map = {cohort: {} for cohort in cohorts}
gene_loss_count_map = {cohort: {} for cohort in cohorts}

snv_prevalence_map = {cohort : {} for cohort in cohorts}
gene_gain_prevalence_map = {cohort: {} for cohort in cohorts}
gene_loss_prevalence_map = {cohort: {} for cohort in cohorts}

variant_type_prevalence_map = {cohort: {'1D':[], '4D':[]} for cohort in cohorts}
cohort_output_strs = {cohort : [] for cohort in cohorts}

######
#
# Helper function for calculating the prevalence of the sweeping allele in the larger cohort
#
######
def get_sweep_prevalence(snp_change, snv_freq_map, private_snv_map):
 
    gene_name, contig, position, variant_type, A1, D1, A2, D2 = snp_change
                
    f1 = A1*1.0/D2
    f2 = A2*1.0/D2
                
    is_reversion = (f1>f2)
                
    location_tuple = (contig, position)
                
    is_private_snv = (location_tuple in private_snv_map)
                
                    
    # Now calculate frequency-stratified version
                
    if location_tuple in snv_freq_map:
        f = snv_freq_map[location_tuple]
    else:
        sys.stderr.write("SNP not in map. Shouldn't happen!\n")
        f = -0.5
                
    # Let's impose that private snvs have zero freq (specifically, lim 0^-)        
    if is_private_snv:
        f = -0.5
                
    # Change f so that it represents
    # frequency of allele at second timepoint
    if is_reversion:
        f = 1-f
        
    return f
    
# Helper functions for stats
from scipy.stats import fisher_exact, ks_2samp, anderson_ksamp
from numpy.random import multinomial as sample_multinomial
from numpy.random import binomial as sample_binomial

# First we can look at prevalence distribution
from scipy.special import gammaln as loggamma

# Helper function to calculate loglikelihood from a multinomial distribution
def multinomial_loglikelihood(ns, ntot, ps):

    return loggamma(ntot+1)-loggamma(ns+1).sum() + (ns*numpy.log(ps+(ns==0))).sum()    

def binomial_loglikelihoods(ns,ntots,ps):
    
    anti_ns = ntots-ns
    
    return loggamma(ntots+1).sum()-loggamma(ns+1).sum() - loggamma(anti_ns+1).sum() + (ns*numpy.log(ps+(ns==0))).sum() + (anti_ns*numpy.log(1-ps+(anti_ns==0))).sum()   

def reversal_loglikelihoods(ns):
    ntots = ns + ns[::-1]
    ps = numpy.ones_like(ns)*0.5    
    return binomial_loglikelihoods(ns,ntots,ps)

def ks_distance(n1s, n2s):
    if len(n1s)==0 or len(n2s)==0:
        return 1e06
        
    ks, dummy = ks_2samp(n1s, n2s) # Kolmogorov-Smirnov version
    #ks = anderson_ksamp([n1s,n2s]) # Anderson-Darling version
    return ks
    
def symmetrized_ks_distance(prevalences):
    prevalences = numpy.array(prevalences)
    symmetrized_prevalences = numpy.hstack([prevalences, (1-prevalences)])
    return ks_distance(prevalences, symmetrized_prevalences)
    


#####################
#
# Do calculation
#
#####################

            
# Load subject and sample metadata
sys.stderr.write("Loading sample metadata...\n")
sample_utils.parse_subject_sample_map
subject_sample_map = sample_utils.parse_subject_sample_map()
sample_order_map = sample_utils.parse_sample_order_map()
sample_country_map = sample_utils.parse_sample_country_map()
sys.stderr.write("Done!\n")
    
good_species_list = parse_midas_data.parse_good_species_list()
if debug:
    good_species_list = good_species_list[0:4]
    default_num_bootstraps = 100
    
num_passed_species = 0

for species_name in good_species_list:

    sys.stderr.write("\nProcessing %s...\n" % species_name)
    
    # First we have to enumerate QP pairs in each cohort
    sys.stderr.write("Enumerating QP pairs...\n")

    # all samples
    all_samples = sample_order_map.keys()

    # list of samples that meet coverage criteria for this species
    highcoverage_samples = set(diversity_utils.calculate_highcoverage_samples(species_name))
    
    # list of samples that meet QP criteria for this species
    haploid_samples = set(diversity_utils.calculate_haploid_samples(species_name))
    
    #print len(all_samples), len(highcoverage_samples), len(haploid_samples)
       
    if len(haploid_samples) < min_haploid_sample_size:
        continue

    all_samples = list(haploid_samples)
    same_sample_idxs, same_subject_idxs, diff_subject_idxs = sample_utils.calculate_ordered_subject_pairs(sample_order_map, all_samples, within_host_type=within_host_type)

    hmp_sample_size = 0        

    qp_sample_sets = {cohort: set() for cohort in cohorts}    
    qp_counts = {cohort:[0,0,0,0] for cohort in cohorts}
    
    for sample_pair_idx in xrange(0,len(same_subject_idxs[0])):
   
        i = same_subject_idxs[0][sample_pair_idx]
        j = same_subject_idxs[1][sample_pair_idx]

        sample_i = all_samples[i]
        sample_j = all_samples[j]
        
        country = sample_country_map[sample_i]
        
        if country not in countries:
            continue
        
        # Figure out cohort
        cohort = country_cohort_map[country]
                
        # Figure out QP status of pair

        
        if not ((sample_i in highcoverage_samples) and (sample_j in highcoverage_samples)):
            # Both are not highcoverage samples
            
            if ((sample_i in highcoverage_samples) or (sample_j in highcoverage_samples)):
                # One sample is high coverage
                qp_counts[cohort][0] += 1
            else:
                # Neither sample is high coverage, ignore
                pass
            
        else:
            
            # Both are highcoverage samples
            
            if (sample_i in haploid_samples) and (sample_j in haploid_samples):
                
                # Both are QP samples!
                
                qp_counts[cohort][1] += 1
                qp_sample_sets[cohort].add(sample_i)
                qp_sample_sets[cohort].add(sample_j) 
                #print sample_i, sample_j   
            
            elif (sample_i not in haploid_samples) and (sample_j not in haploid_samples):
                # pair that is non-QP at both timepoints
                qp_counts[cohort][2] += 1
            
            else:
                # pair that is QP at one timepoint and non-QP at another
                qp_counts[cohort][3] += 1
    
    sys.stderr.write("Done!\n")
    
    output_strs.append(species_name)
    for cohort in cohorts:    
        output_strs.append("%s: %d temporal QP pairs" % (cohort, qp_counts[cohort][1]))
        
    combined_sample_set = set()
    for cohort in cohorts:
        combined_sample_set.update(qp_sample_sets[cohort])
    combined_samples = list(sorted(combined_sample_set))
    combined_sample_idx_map = {combined_samples[i] : i for i in xrange(0,len(combined_samples))}    
    
    qp_sample_lists = {cohort: list(sorted(qp_sample_sets[cohort])) for cohort in cohorts}
    
    sample_size = len(qp_sample_sets['hmp'])
        
    if sample_size < min_sample_size:
        continue
    
    hmp_species_qp_counts[species_name] = qp_counts['hmp']
    twin_species_qp_counts[species_name] = qp_counts['twins']
    
    for cohort in cohorts:
        species_snp_change_distribution[cohort][species_name] = []
        species_snp_nerrs[cohort][species_name] = []
        species_gene_change_distribution[cohort][species_name] = []
        species_gene_nerrs[cohort][species_name] = []
    
    sys.stderr.write("Proceeding with %d HMP longitudinal comparisons!\n" % (sample_size))
    
    import calculate_private_snvs
    private_snv_map = calculate_private_snvs.load_private_snv_map(species_name)
    
    import calculate_snp_prevalences
    snv_freq_map = calculate_snp_prevalences.parse_population_freqs(species_name,polarize_by_consensus=True)
    snv_freq_keys = snv_freq_map.keys()
    snv_freq_values = snv_freq_map.values()
    
    import core_gene_utils
    gene_freq_map = core_gene_utils.parse_gene_freqs(species_name)
    gene_freq_values = numpy.array(gene_freq_map.values())
    gene_freq_weights = gene_freq_values*1.0/gene_freq_values.sum()
    
    sys.stderr.write("Loading pre-computed substitution rates for %s...\n" % species_name)
    substitution_rate_map = calculate_substitution_rates.load_substitution_rate_map(species_name)
    sys.stderr.write("Calculating SNV matrix...\n")
    dummy_samples, snp_mut_difference_matrix, snp_rev_difference_matrix, snp_mut_opportunity_matrix, snp_rev_opportunity_matrix = calculate_substitution_rates.calculate_mutrev_matrices_from_substitution_rate_map(substitution_rate_map, 'all', allowed_samples=combined_samples)
    
    snp_difference_matrix = snp_mut_difference_matrix+snp_rev_difference_matrix
    snp_opportunity_matrix = snp_mut_opportunity_matrix+snp_rev_opportunity_matrix    
    snp_substitution_rate =     snp_difference_matrix*1.0/(snp_opportunity_matrix+(snp_opportunity_matrix==0))
    sys.stderr.write("Done!\n")
    
    sys.stderr.write("Loading gene matrix...\n")
    gene_samples, gene_loss_difference_matrix, gene_gain_difference_matrix, gene_loss_opportunity_matrix, gene_gain_opportunity_matrix = calculate_substitution_rates.calculate_mutrev_matrices_from_substitution_rate_map(substitution_rate_map, 'genes', allowed_samples=combined_samples)
    gene_difference_matrix = gene_gain_difference_matrix + gene_loss_difference_matrix
    gene_opportunity_matrix = gene_loss_opportunity_matrix
    gene_difference_matrices = {'gains': gene_gain_difference_matrix, 'losses': gene_loss_difference_matrix}
    sys.stderr.write("Done!\n")

    
    sys.stderr.write("Loading 1D & 4D opportunity matrices...\n")
    
    opportunity_matrices = {}
    difference_matrices = {}

    twin_opportunity_matrices = {}
    twin_difference_matrices = {}

    for var_type in variant_types:
        
        # First do HMP
        dummy_samples, difference_matrix, opportunity_matrix =    calculate_substitution_rates.calculate_matrices_from_substitution_rate_map(substitution_rate_map, var_type, allowed_samples=combined_samples)
    
        difference_matrices[var_type] = difference_matrix
        opportunity_matrices[var_type] = opportunity_matrix
        
    sys.stderr.write("Done!\n")
    
    sys.stderr.write("Loading pre-computed temporal changes...\n")
    temporal_change_map = calculate_temporal_changes.load_temporal_change_map(species_name)
    sys.stderr.write("Done!\n")


    ### Now loop over different cohorts
    for cohort in cohorts:
        
        modification_difference_threshold = modification_difference_thresholds[cohort]
        replacement_difference_threshold = replacement_difference_thresholds[cohort]
        
        desired_samples = qp_sample_lists[cohort]
        
        same_sample_idxs, same_subject_idxs, diff_subject_idxs = sample_utils.calculate_ordered_subject_pairs(sample_order_map, desired_samples, within_host_type=within_host_type)
        
        #apply_sample_index_map_to_indices(sample_idx_map, idxs):
        #new_idxs = (numpy.array([sample_idx_map[i] for i in idxs[0]]), numpy.array([sample_idx_map[i] for i in idxs[1]]))
    
        for sample_pair_idx in xrange(0,len(same_subject_idxs[0])):
        
            
            sample_i = desired_samples[same_subject_idxs[0][sample_pair_idx]] 
            sample_j = desired_samples[same_subject_idxs[1][sample_pair_idx]]
    
            i = combined_sample_idx_map[sample_i]
            j = combined_sample_idx_map[sample_j]
            
            good_idxs = sample_utils.calculate_samples_in_different_subjects( subject_sample_map, combined_samples, sample_i)
            good_idxs *= ( (snp_opportunity_matrix[i,:]>0.5) * (gene_opportunity_matrix[i,:]>0.5) )
            
            if good_idxs.sum() < 1:
                continue
                
            L, perr, mutations, reversions = calculate_temporal_changes.calculate_mutations_reversions_from_temporal_change_map(temporal_change_map, sample_i, sample_j)
        
            if L<config.min_opportunities:
                continue
        
            nerr = L*perr
        
            num_mutations = len(mutations)
            num_reversions = len(reversions)
            num_snp_changes = num_mutations+num_reversions
        
            gene_L, gene_perr, gains, losses = calculate_temporal_changes.calculate_gains_losses_from_temporal_change_map(temporal_change_map, sample_i, sample_j) #, min_normal_copynum = 0.6, max_normal_copynum = 1.2)
        
            gene_nerr = gene_L*gene_perr
            num_gains = len(gains)
            num_losses = len(losses)
            num_gene_changes = num_gains+num_losses
        
            if (perr<-0.5) or (gene_perr < -0.5):
                continue
        
            if (nerr > max([0.5, 0.1*num_snp_changes])) or (gene_nerr > max([0.5, 0.1*num_gene_changes])):
                continue # Only take things with low-ish FPR
        
            # Species specific distributions
            species_snp_change_distribution[cohort][species_name].append( num_snp_changes)
            species_snp_nerrs[cohort][species_name].append( nerr)
            species_gene_change_distribution[cohort][species_name].append( num_gene_changes)
            
            species_gene_nerrs[cohort][species_name].append( gene_nerr)
            
            # Pooled distributions
            pooled_snp_change_distribution[cohort].append(num_snp_changes)
            pooled_gene_change_distribution[cohort].append(num_gene_changes)
            
            pooled_snp_length_distribution[cohort].append(L)
            pooled_gene_length_distribution[cohort].append(gene_L)
            
            # Matched between-host samples
            # typical
            pooled_between_snp_change_distribution[cohort].append( choice( snp_difference_matrix[i, good_idxs] ) ) 
            pooled_between_gene_change_distribution[cohort].append( choice(gene_difference_matrix[i, good_idxs]) )           
            # minimum
            pooled_min_between_snp_change_distribution[cohort].append( snp_difference_matrix[i, good_idxs].min() )
            pooled_min_between_gene_change_distribution[cohort].append( gene_difference_matrix[i, good_idxs].min() )
            
            if (cohort=='young_twins'):
            
                cohort_output_strs[cohort].append("%s, %s, %s: n_snv=%d, n_gene=%d" % (species_name, sample_i, sample_j, num_snp_changes, num_gene_changes))
            
            # Store sample names for replacement to see if many species are
            # replaced in the same individual
            if (num_snp_changes>=replacement_difference_threshold):
                sample_pair = (sample_i, sample_j)
                if sample_pair not in replacement_map[cohort]:
                    replacement_map[cohort][sample_pair] = []
                replacement_map[cohort][sample_pair].append(species_name)
        
        
            if cohort=='twins' and (num_snp_changes<replacement_difference_threshold):
                cohort_output_strs[cohort].append("%s, %s, %s: n_snv=%d, n_gene=%d" % (species_name, sample_i, sample_j, num_snp_changes, num_gene_changes))
            
        
            # If deemed a modification, investigate properties of SNVs and genes        
            if (num_snp_changes<=modification_difference_threshold):
        
                    
                for snp_change in (mutations+reversions):        
                
                    variant_type = snp_change[3]
                
                    f = get_sweep_prevalence(snp_change, snv_freq_map, private_snv_map)
                
                    f_idx = ((f>derived_freq_bins[:-1])*(f<=derived_freq_bins[1:])).argmax()    
                
                    if variant_type in variant_types:
                        total_freq_snps[cohort][variant_type][f_idx] += 1
                    
                        # Calculate null version based on # of opportunities
                        total_opportunities = 0.0
                        for other_variant_type in variant_types:
                            total_opportunities = opportunity_matrices[other_variant_type][i,j]
                        
                        for other_variant_type in variant_types:
                            total_null_freq_snps[cohort][other_variant_type][f_idx] += opportunity_matrices[other_variant_type][i,j]/total_opportunities
                    
                    total_freq_all_snps[cohort][f_idx] += 1
                
                    if (sample_i, sample_j, species_name) not in snv_prevalence_count_map[cohort]:
                        snv_prevalence_count_map[cohort][(sample_i, sample_j, species_name)] = numpy.zeros_like(total_freq_all_snps[cohort])
                        snv_prevalence_map[cohort][(sample_i, sample_j, species_name)] = []
                        
                    snv_prevalence_count_map[cohort][(sample_i, sample_j, species_name)][f_idx] += 1
                    
                    snv_prevalence_map[cohort][(sample_i, sample_j, species_name)].append(f) 
                    
                    if variant_type not in variant_type_prevalence_map[cohort]:
                        variant_type_prevalence_map[cohort][variant_type] = []
                    
                    variant_type_prevalence_map[cohort][variant_type].append(f) 
                    
                    # Now draw a null prevalence from the genome
                    L = snp_opportunity_matrix[i,j]
                    L_snv = len(snv_freq_map) # A slight overestimate
                    snv_fraction = L_snv*1.0/L
                    num_bootstraps = 10
                    for bootstrap_idx in xrange(0,num_bootstraps):
                    
                        if random()<snv_fraction:
                            # A polymorphic site
                        
                            random_snv_idx = randint(0,len(snv_freq_keys))
                            random_snv_location = snv_freq_keys[random_snv_idx]
                            f = snv_freq_values[random_snv_idx]
                            
                            rev_f = 1-f
                            
                            if random_snv_location in private_snv_map:
                                # A private SNV. Use private bins
                                # use private         
                                f_idx = 0
                                rev_f_idx = -1
                            else:
                            
                                f_idx = ((f>derived_freq_bins[:-1])*(f<=derived_freq_bins[1:])).argmax()    
                    
                                rev_f_idx = ((rev_f>derived_freq_bins[:-1])*(rev_f<=derived_freq_bins[1:])).argmax()  
                    
                    
                            # Now add in probability weight
                            total_null_freq_all_snps[cohort][f_idx] += (1-f)*1.0/num_bootstraps
                            total_null_freq_all_snps[cohort][rev_f_idx] += f*1.0/num_bootstraps
                        
                        else:
                            # A truly invariant site
                            total_null_freq_all_snps[cohort][0] += 1.0/num_bootstraps
                           
            
                for gene_change in gains:
                    gene_name = gene_change[0]
                    
                    if gene_name in gene_freq_map:
                        f = gene_freq_map[gene_name]
                    else:
                        f = 0
                    
                    f_idx = ((f>gene_freq_bins[:-1])*(f<=gene_freq_bins[1:])).argmax()    
                    total_freq_gains[cohort][f_idx] += 1
                    
                    if (sample_i, sample_j, species_name) not in gene_gain_count_map[cohort]:
                        
                        gene_gain_count_map[cohort][(sample_i, sample_j, species_name)] = numpy.zeros_like(total_freq_gains[cohort])
                        gene_loss_count_map[cohort][(sample_i, sample_j, species_name)] = numpy.zeros_like(total_freq_losses[cohort])
                        gene_gain_prevalence_map[cohort][(sample_i, sample_j, species_name)] = []
                        gene_loss_prevalence_map[cohort][(sample_i, sample_j, species_name)] = []     
                        
                        
                    gene_gain_count_map[cohort][(sample_i, sample_j, species_name)][f_idx] += 1
                    gene_gain_prevalence_map[cohort][(sample_i, sample_j, species_name)].append(f)
                    
                
                        
                    
                for gene_change in losses:
                    gene_name = gene_change[0]
                    
                    if gene_name in gene_freq_map:
                        f = gene_freq_map[gene_name]
                    else:
                        f = 0
                    
                    f_idx = ((f>gene_freq_bins[:-1])*(f<=gene_freq_bins[1:])).argmax()    
                    total_freq_losses[cohort][f_idx] += 1
                    
                    if (sample_i, sample_j, species_name) not in gene_gain_count_map[cohort]:
                        
                        gene_gain_count_map[cohort][(sample_i, sample_j, species_name)] = numpy.zeros_like(total_freq_gains[cohort])
                        gene_loss_count_map[cohort][(sample_i, sample_j, species_name)] = numpy.zeros_like(total_freq_losses[cohort])
                        gene_gain_prevalence_map[cohort][(sample_i, sample_j, species_name)] = []
                        gene_loss_prevalence_map[cohort][(sample_i, sample_j, species_name)] = []     
                        
                    gene_loss_count_map[cohort][(sample_i, sample_j, species_name)][f_idx] += 1
                    
                    gene_loss_prevalence_map[cohort][(sample_i, sample_j, species_name)].append(f)
                
                    
                    num_bootstraps = 10
                    fs = choice(gene_freq_values, size=num_bootstraps, p=gene_freq_weights)
                    for f in fs:
                        f_idx = ((f>gene_freq_bins[:-1])*(f<=gene_freq_bins[1:])).argmax()    
                        total_null_freq_losses[cohort][f_idx] += 1.0/num_bootstraps
                
        
sys.stderr.write("Done looping over species!\n")    

output_strs.append("Young twins data:")
output_strs.append("\n".join(cohort_output_strs['young_twins']))
output_strs.append("Old twins modifications:")
output_strs.append("\n".join(cohort_output_strs['twins']))

#print "All twins modifications:"
#print "\n".join(output_strs['twins'])

output_strs.append("Printing replacement map!")
for sample_pair in replacement_map['hmp'].keys():
    output_strs.append(" ".join([str(sample_pair), str(len(replacement_map['hmp'][sample_pair])), str(replacement_map['hmp'][sample_pair])]))


species_names = []
sample_sizes = []

for species_name in species_snp_change_distribution['hmp']:
    sample_size = len(species_snp_change_distribution['hmp'][species_name])
    if sample_size > 0:
        species_names.append(species_name)
        sample_sizes.append( sample_size )
    
# sort in descending order of sample size
# Sort by num haploids    
sample_sizes, species_names = zip(*sorted(zip(sample_sizes, species_names),reverse=True))

max_sample_size = 10**(ceil(log10(sample_sizes[0]*1.0)))

     
sys.stderr.write("Postprocessing %d species!\n" % len(species_names))

# Clip things so that we can plot on loglog axis
for cohort in cohorts:
    pooled_snp_change_distribution[cohort] = numpy.clip(numpy.array(pooled_snp_change_distribution[cohort]), 1e-06, 1e08)
    pooled_gene_change_distribution[cohort] = numpy.clip(numpy.array(pooled_gene_change_distribution[cohort]), 1e-06, 1e08)
    
    pooled_between_snp_change_distribution[cohort] = numpy.clip(numpy.array(pooled_between_snp_change_distribution[cohort]), 1e-06, 1e08)
    pooled_between_gene_change_distribution[cohort] = numpy.clip(numpy.array(pooled_between_gene_change_distribution[cohort]), 1e-06, 1e08)
    
    pooled_min_between_snp_change_distribution[cohort] = numpy.clip(numpy.array(pooled_min_between_snp_change_distribution[cohort]), 1e-06, 1e08)
    pooled_min_between_gene_change_distribution[cohort] = numpy.clip(numpy.array(pooled_min_between_gene_change_distribution[cohort]), 1e-06, 1e08)
    

##############################################################################
#
# Plot results
#
##############################################################################

ymin=0
ymax=0
for cohort in cohorts:
    
    if cohort=='hmp':
    
        # If HMP, plot species-specific snp and gene change distribution
        
        change_axis_labels = []
        
        for species_idx in xrange(0,len(species_names)):

            species_name = species_names[species_idx]
    
            # first figure out whether # of SNV changes is significant
            # at 10% species-wide FDR
            
            species_snp_change_distribution[cohort][species_name] = numpy.array(species_snp_change_distribution[cohort][species_name])
            species_snp_nerrs[cohort][species_name] = numpy.array(species_snp_nerrs[cohort][species_name])
            
            species_gene_change_distribution[cohort][species_name] = numpy.array(species_gene_change_distribution[cohort][species_name])
            species_gene_nerrs[cohort][species_name] = numpy.array(species_gene_nerrs[cohort][species_name])
            
            modification_idxs = species_snp_change_distribution[cohort][species_name]<=modification_difference_thresholds[cohort]
            
            total_snv_changes = species_snp_change_distribution[cohort][species_name][modification_idxs].sum()
            total_snv_nerr = species_snp_nerrs[cohort][species_name][modification_idxs].sum()
            snv_is_significant = (total_snv_nerr < 0.1*total_snv_changes)
            
            
            # then figure out whether # of gene changes is significant
            total_gene_changes = species_gene_change_distribution[cohort][species_name][modification_idxs].sum()
            total_gene_nerr = species_gene_nerrs[cohort][species_name][modification_idxs].sum()
            gene_is_significant = (total_gene_nerr < 0.1*total_gene_changes)
            
            label = species_name
            change_axis_labels.append(species_name)
            if snv_is_significant and gene_is_significant:
                label += " (S,G)"
                pass
            elif snv_is_significant and (not gene_is_significant):
                label += " (S)"
                pass
            elif gene_is_significant and (not snv_is_significant):
                label += " (G)"
                pass
            else:
                pass 
            
            # Plot SNVs
            xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(numpy.clip(species_snp_change_distribution[cohort][species_name],5e-01,1e09), min_x=1e-07, max_x=1e09)

            line, = species_snp_axis.step(xs,ns,'-',linewidth=1, where='pre',zorder=4)
            color = pylab.getp(line,'color')
            
            # Plot genes
            xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(numpy.clip(species_gene_change_distribution[cohort][species_name],5e-01,1e09), min_x=1e-07, max_x=1e09)

            line, = species_gene_axis.step(xs,ns,'-',color=color,linewidth=1, where='pre',zorder=4)
            
            species_legend_axis.plot([-1],[-1],'-',color=color, linewidth=1,label=label)
            
            snp_changes = list(sorted(species_snp_change_distribution[cohort][species_name]))
            
            gene_changes = list(sorted(species_gene_change_distribution[cohort][species_name]))
            
            # Plot SNP change colors
            for idx in xrange(0,len(snp_changes)):
        
                if snp_changes[idx]<0.5:
                    colorVal='0.7'
                else:
                    colorVal = scalarMap.to_rgba(log10(snp_changes[idx]))
        
                change_axis.fill_between([species_idx-0.3,species_idx+0.3], [idx,idx],[idx+1.05,idx+1.05],color=colorVal,linewidth=0)
    
        
                if snv_is_significant:
            
                    change_axis.text(species_idx, len(snp_changes),'*',fontsize=4)
    
            for idx in xrange(0,len(gene_changes)):
        
                if gene_changes[idx]<0.5:
                    colorVal='0.7'
                else:
                    colorVal = scalarMap.to_rgba(log10(gene_changes[idx]))
        
                change_axis.fill_between([species_idx-0.3,species_idx+0.3], [-idx-1.05,-idx-1.05],[-idx,-idx],color=colorVal,linewidth=0)
        
                if gene_is_significant:
            
                    change_axis.text(species_idx, -len(gene_changes)-3,'*',fontsize=4)
            

        change_axis.set_ylim([-90,90])
        change_axis.set_xlim([-1,len(change_axis_labels)])
        change_axis.plot([-1,len(change_axis_labels)],[0,0],'k-')

        change_axis.set_yticks([-90,-80,-70, -60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80,90])
        change_axis.set_yticklabels(['90','80','70', '60','50','40','30','20','10','0','10','20','30','40','50','60','70','80','90'])

        xticks = numpy.arange(0,len(change_axis_labels))
        xticklabels = change_axis_labels

        change_axis.set_xticks(xticks)
        change_axis.set_xticklabels(xticklabels, rotation='vertical',fontsize=4)

        m = change_axis.scatter([200],[1],c=[0], vmin=0, vmax=vmax, cmap=cmap, marker='^')

        cbar = fig.colorbar(m,cax=cax,orientation='vertical', ticks=[0,1,2,3])
        cbar.set_ticklabels(['$1$','$10$','$10^{2}$','$10^{3}$'])
        cbar.set_label('Number of changes',rotation=270,labelpad=10)
        cl = pylab.getp(cbar.ax, 'ymajorticklabels')
        pylab.setp(cl, fontsize=9) 
        #fig.text(0.945,0.05,'$\\pi/\\pi_0$',fontsize=12)

        cbar.ax.tick_params(labelsize=5)
        change_axis.text(20,25,'SNVs',fontsize=5)
        change_axis.text(20,-20,'Genes',fontsize=5)
        
            
        species_gene_axis.set_ylim([0.3,max_sample_size])
        species_snp_axis.set_ylim([0.3,max_sample_size])
        species_gene_axis.set_yticklabels([])        
        species_legend_axis.legend(loc='center', frameon=False,fontsize=5,numpoints=1,ncol=3,handlelength=1)   
    
    # Now plot pooled versions
        
    # First Break distributions into replacements and modifications
    modification_idxs = pooled_snp_change_distribution[cohort]<=modification_difference_thresholds[cohort]
    replacement_idxs = pooled_snp_change_distribution[cohort]>=replacement_difference_thresholds[cohort]
    non_replacement_idxs = numpy.logical_not(replacement_idxs)
    
    no_snv_change_idxs = numpy.logical_and(modification_idxs, pooled_snp_change_distribution[cohort]<0.5)
    true_snv_modification_idxs = numpy.logical_and(modification_idxs, pooled_snp_change_distribution[cohort]>0.5)
    
    true_gene_modification_idxs = numpy.logical_and(modification_idxs, pooled_gene_change_distribution[cohort]>0.5)
    
    true_gene_and_snv_modification_idxs = numpy.logical_and(true_snv_modification_idxs, true_gene_modification_idxs)
    true_gene_and_no_snv_modification_idxs = numpy.logical_and(no_snv_change_idxs, true_gene_modification_idxs)
    
      
    true_modification_idxs = numpy.logical_or(true_snv_modification_idxs, true_gene_modification_idxs)
    
            
    modification_snp_change_distribution = pooled_snp_change_distribution[cohort][modification_idxs]
    replacement_snp_change_distribution = pooled_snp_change_distribution[cohort][replacement_idxs]
       
    no_snv_gene_change_distribution = pooled_gene_change_distribution[cohort][no_snv_change_idxs] 
        
    true_modification_gene_change_distribution = pooled_gene_change_distribution[cohort][true_snv_modification_idxs]
    
    modification_gene_change_distribution = pooled_gene_change_distribution[cohort][modification_idxs]
    
    replacement_gene_change_distribution = pooled_gene_change_distribution[cohort][replacement_idxs]
        
    current_num_species = 0
    for species in species_snp_change_distribution[cohort]:
        if len(species_snp_change_distribution[cohort][species])>0:
            current_num_species += 1
    
    snv_modification_species = {}
    for sample_pair in snv_prevalence_map[cohort]:
        if len(snv_prevalence_map[cohort][sample_pair])>0:
            species_name = sample_pair[-1]
            if species_name not in snv_modification_species:
                snv_modification_species[species_name] = 0
                
            snv_modification_species[species_name] += len(snv_prevalence_map[cohort][sample_pair])

    
    # sort species by descending number of changes
    sorted_snv_modification_species = []
    for species_name in sorted(snv_modification_species.keys(), key=lambda x: snv_modification_species[x], reverse=True):
        sorted_snv_modification_species.append((species_name, snv_modification_species[species_name]))

    
    gene_modification_species = {}
    for sample_pair in gene_gain_prevalence_map[cohort]:
        num_gene_changes = len(gene_gain_prevalence_map[cohort][sample_pair])+len(gene_loss_prevalence_map[cohort][sample_pair])
        if num_gene_changes>0:
        
            species_name = sample_pair[-1]
            
            if species_name not in gene_modification_species:
                gene_modification_species[species_name]=0
                
            gene_modification_species[species_name] += num_gene_changes
            
    # sort species by descending number of changes
    sorted_gene_modification_species = []
    for species_name in sorted(gene_modification_species.keys(), key=lambda x: gene_modification_species[x], reverse=True):
        sorted_gene_modification_species.append((species_name, gene_modification_species[species_name]))
            
    output_strs.append("Cohort: %s" % cohort)
    output_strs.append("".join(['-']*80))
    output_strs.append("%d total comparisons across %d species" % (len(pooled_snp_change_distribution[cohort]), current_num_species))
    
    output_strs.append("%d replacements (%g, %d total SNVs, %d total genes)" % (replacement_idxs.sum(), replacement_idxs.sum()*1.0/len(replacement_idxs), replacement_snp_change_distribution.sum(), replacement_gene_change_distribution.sum()))
    
    if cohort=='twins':
        output_strs.append(" ".join(["Non-replacement SNVs", str(numpy.median(pooled_snp_change_distribution[cohort][non_replacement_idxs])), str(pooled_snp_change_distribution[cohort][non_replacement_idxs])]))
        output_strs.append(" ".join(["Non-replacement genes", str(numpy.median(pooled_gene_change_distribution[cohort][non_replacement_idxs])), str(pooled_gene_change_distribution[cohort][non_replacement_idxs])]))
        
            
    output_strs.append("%d total modifications (%g)" % (true_modification_idxs.sum(), true_modification_idxs.sum()*1.0/len(true_modification_idxs)))
    
    output_strs.append("%d SNV modifications (%g, %d total SNVs in %d pairs across %d species)" % (true_snv_modification_idxs.sum(), true_snv_modification_idxs.sum()*1.0/len(true_snv_modification_idxs), modification_snp_change_distribution.sum(), len(snv_prevalence_count_map[cohort]), len(snv_modification_species)))
    output_strs.append(str(sorted_snv_modification_species))
    
    output_strs.append( "%d gene modifications (%g, %d total genes in %d pairs across %d species)" % (true_gene_modification_idxs.sum(), true_gene_modification_idxs.sum()*1.0/len(true_gene_modification_idxs), modification_gene_change_distribution.sum(), len(gene_gain_count_map[cohort]), len(gene_modification_species)) )
    output_strs.append(str(sorted_gene_modification_species))
    
    output_strs.append( "p(gene change | snv change) = %g" % (true_gene_and_snv_modification_idxs.sum()*1.0/true_snv_modification_idxs.sum()))
    
    output_strs.append( "p(gene change | no snv change) = %g" % (true_gene_and_no_snv_modification_idxs.sum()*1.0/no_snv_change_idxs.sum()))
    
    n11 = true_gene_and_snv_modification_idxs.sum()
    n10 = true_gene_and_no_snv_modification_idxs.sum()
    n01 = true_snv_modification_idxs.sum()-n11
    n00 = no_snv_change_idxs.sum()-n10
    
    oddsratio, pvalue = fisher_exact([[n11, n10],[n01, n00]])
    output_strs.append( "Fisher exact pvalue = %g" % pvalue)
    # comment
    
    if cohort=='hmp':
        
        # Plot snp and gene length distribution
        # Plot within and between for snvs
        
        xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(pooled_snp_length_distribution[cohort], min_x=1e05, max_x=1e07)

        genome_length_axis.step(xs,ns/ns[0],'-',color='#08519c',linewidth=1, label=('Within-host (n=%d)' % ns[0]), where='pre',zorder=4)
        
        genome_length_axis.semilogx([1e02],[1],'k.')
        genome_length_axis.set_xlim([1e05,1e07])
        genome_length_axis.set_ylim([0,1.1])
        
        xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(pooled_gene_length_distribution[cohort], min_x=3e02, max_x=3e04)

        pangenome_length_axis.step(xs,ns/ns[0],'-',color='#08519c',linewidth=1, where='pre',zorder=4)
        
        pangenome_length_axis.semilogx([1e01],[1],'k.')
        pangenome_length_axis.set_xlim([3e02,3e04])
        pangenome_length_axis.set_ylim([0,1.1])
        
        # Plot average within and between
        # Plot SNP change averages

        avg_axis.bar([0.75], [pooled_snp_change_distribution[cohort].mean()],width=0.5,facecolor='#08519c',linewidth=0,log=True)
        avg_axis.bar([1.75], [pooled_between_snp_change_distribution[cohort].mean()],width=0.5, facecolor='r',alpha=0.5,linewidth=0,log=True)

        avg_axis.set_ylim([1,1e05])
        avg_axis.set_ylabel('Avg # SNV differences')
        avg_axis.set_xticks([1,2])
        avg_axis.set_xticklabels(['Within\nhost','Between\nhost'])
    
        # Plot within and between for snvs
        
        xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(pooled_snp_change_distribution[cohort], min_x=1e-07, max_x=1e09)

        good_idxs = (pooled_snp_change_distribution[cohort]>0.5)*(pooled_snp_change_distribution[cohort]<1.5)

        print good_idxs.sum(), "one change"

        pooled_snp_axis.step(xs,ns/ns[0],'-',color='#08519c',linewidth=1, label=('Within-host (n=%d)' % ns[0]), where='pre',zorder=4)
        #pooled_snp_axis.plot(xs,ns/ns[0],'.-',color='#08519c',linewidth=1, label=('Within-host (n=%d)' % ns[0]),zorder=4)

      

        # This lets you set the y range

        ymin = 1.0/ns[0]
        ymax = 1.3

        # Put on log scale
        pooled_snp_axis.set_ylim([ymin,ymax])
        pooled_gene_axis.set_ylim([ymin,ymax])
        
        pooled_gene_axis.set_yticklabels([])

        


        # Save intermediate version (for Keynote animations)
        fig2.savefig('%s/figure_6.1.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight',transparent=True)

        # Now the typical between host

        # Random between host
        xs, ns =       stats_utils.calculate_unnormalized_survival_from_vector(pooled_between_snp_change_distribution[cohort], min_x=1e-07, max_x=1e09)

        # Min between host
        #xs, ns =     stats_utils.calculate_unnormalized_survival_from_vector(pooled_min_between_snp_change_distribution, min_x=1e-02, max_x=1e09)

        pooled_snp_axis.step(xs,ns/ns[0],'-',color='r',linewidth=0.5, alpha=0.5, label='Between-host', where='pre',zorder=2)
        #pooled_snp_axis.plot(xs,ns/ns[0],'.-',color='r',linewidth=0.5, alpha=0.5, label='Between-host',zorder=2)

        # Save intermediate version (for Keynote animations)
        fig2.savefig('%s/figure_6.2.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight',transparent=True)

        # Genes: Plot modification and replacement within (separately)
        
        
        xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(pooled_between_gene_change_distribution[cohort], min_x=1e-07, max_x=1e09)

        pooled_gene_axis.step(xs,ns/ns[0],linewidth=1,color='r', alpha=0.5, label='Between-host',zorder=1)
            
        xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(true_modification_gene_change_distribution, min_x=1e-07, max_x=1e09)

        pooled_gene_axis.step(xs,ns/ns[0],'-',color='#08519c',linewidth=1, label='Within-host',zorder=3,where='pre',path_effects=[pe.Stroke(linewidth=5, foreground='#9ecae1'), pe.Normal()])

        #pooled_gene_axis.loglog([1e-01,1e05],[1.0/ns[0],1.0/ns[0]],'k:')

        #pooled_gene_axis.set_ylim([1.0/ns[0],1.3])
        pooled_gene_axis.set_yticklabels([])
        #pooled_snp_axis.set_yticklabels(['0.01','0.1','1'])


        xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(replacement_gene_change_distribution, min_x=1e-07, max_x=1e09)

        pooled_gene_axis.step(xs,ns/ns[0],'-',color='#08519c',linewidth=1, label='Within-host',zorder=2,where='pre',path_effects=[pe.Stroke(linewidth=5, foreground='#fee0d2'), pe.Normal()])
    
        xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(no_snv_gene_change_distribution, min_x=1e-07, max_x=1e09)

        pooled_gene_axis.step(xs,ns/ns[0],'-',color='#08519c',linewidth=1, label='Within-host',zorder=2,where='pre',path_effects=[pe.Stroke(linewidth=5, foreground='0.8'), pe.Normal()])
            
    elif cohort=='twins':
            
        # Plot all for snvs, all for genes 
        
        # Now do twins

        xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(pooled_snp_change_distribution[cohort], min_x=1e-07, max_x=1e09)

        pooled_snp_axis.step(xs,ns/ns[0],'-',color='#8856a7',linewidth=1, label=('Twins (n=%d)' % ns[0]), where='pre',zorder=4)

        young_snp_axis.step(xs,ns/ns[0],'-',color='#8856a7',linewidth=1, label=('Adult twins (n=%d)' % ns[0]), where='pre',zorder=4)

        #pooled_snp_axis.plot(xs,ns/ns[0],'.-',color='#8856a7',linewidth=1, label=('Twins (n=%d)' % ns[0]),zorder=2)



        # Save intermediate version (for Keynote animations)
        fig2.savefig('%s/figure_6.3.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight',transparent=True)

        # Now fill in the graphics

        pooled_snp_axis.fill_between([1e-01,1], [ymin,ymin],[ymax,ymax],color='0.8',zorder=1)
        pooled_snp_axis.fill_between([ 1e0,modification_difference_thresholds['hmp']],[ymin,ymin],[ymax,ymax],color='#deebf7',zorder=1)
        pooled_snp_axis.fill_between([ replacement_difference_thresholds['hmp'],1e05],[ymin,ymin],[ymax,ymax],color='#fee0d2',zorder=1)

        pooled_snp_axis.text( exp((log(1e05)+log(replacement_difference_thresholds['hmp']))/2), ymax*1.2, 'putative\nreplacement',fontsize=6,fontstyle='italic',ha='center',color='#fc9272',zorder=1)
        pooled_snp_axis.text( exp((log(1)+log(modification_difference_thresholds['hmp']))/2), ymax*1.2, 'putative\nmodification',fontsize=6,fontstyle='italic',ha='center',color='#9ecae1',zorder=1)
        
        
        xs, ns = stats_utils.calculate_unnormalized_survival_from_vector( pooled_gene_change_distribution[cohort], min_x=1e-07, max_x=1e09)

        pooled_gene_axis.step(xs,ns/ns[0],'-',color='#8856a7',linewidth=1, label='Twin',zorder=3,where='pre')

        young_gene_axis.step(xs,ns/ns[0],'-',color='#8856a7',linewidth=1,zorder=3,where='pre')

        
    elif cohort=='young_twins':
        
        # Plot all for snvs, all for genes
        xs, ns = stats_utils.calculate_unnormalized_survival_from_vector( pooled_snp_change_distribution[cohort], min_x=1e-07, max_x=1e09)

        young_ymin = 1.0/ns[0]
        young_ymax = 1.3

        print pooled_snp_change_distribution[cohort]
        print xs
        print ns
        print ns/ns[0]

        young_snp_axis.step(xs,ns/ns[0],'-',color='#8856a7',linewidth=1, label=('Younger twins (n=%d)' % ns[0]), where='pre',zorder=2,alpha=0.5)
    
        xs, ns = stats_utils.calculate_unnormalized_survival_from_vector( pooled_gene_change_distribution[cohort], min_x=1e-07, max_x=1e09)

        young_gene_axis.step(xs,ns/ns[0],'-',color='#8856a7',linewidth=1, where='pre',zorder=2,alpha=0.5)
        
        young_snp_axis.set_ylim([young_ymin, young_ymax])
        young_gene_axis.set_ylim([young_ymin, young_ymax])
        young_gene_axis.set_yticklabels([])
    
        young_snp_axis.legend(loc='lower left',frameon=False,fontsize=6,numpoints=1,handlelength=1)   


    # Now to SNV and gene prevalence
    if cohort=='hmp':
        frequency_axis = hmp_frequency_axis
        gene_frequency_axis = hmp_gene_frequency_axis
        gene_legend_axis = hmp_gene_legend_axis
    elif cohort=='twins':
        frequency_axis = twin_frequency_axis
        gene_frequency_axis = twin_gene_frequency_axis
        gene_legend_axis = twin_gene_legend_axis
    else:
        frequency_axis = young_twin_frequency_axis
        gene_frequency_axis = young_twin_gene_frequency_axis
        gene_legend_axis = young_twin_gene_legend_axis

    frequency_axis.bar(derived_virtual_freqs, total_freq_snps[cohort]['4D'],width=0.3,linewidth=0,facecolor='#b3de69',label='syn (4D)',zorder=3)

    frequency_axis.bar(derived_virtual_freqs, total_freq_snps[cohort]['1D']+total_freq_snps[cohort]['4D'],width=0.3,linewidth=0,facecolor='#ff7f00',label='non (1D)',zorder=2)

    frequency_axis.bar(derived_virtual_freqs, total_freq_all_snps[cohort],width=0.3,linewidth=0,facecolor='#b15928',label='(2D & 3D)',zorder=1)


    frequency_axis.bar(derived_virtual_freqs-0.3, total_null_freq_all_snps[cohort],width=0.3,linewidth=0,facecolor='0.7',label='de novo\nexpectation',zorder=0)

    frequency_axis.legend(loc='upper center',frameon=False,fontsize=5,numpoints=1,ncol=2,handlelength=1)   

    gene_frequency_axis.bar(gene_gain_virtual_freqs, total_freq_gains[cohort],width=0.3,linewidth=0,facecolor='#b3de69',label='gain')

    gene_frequency_axis.bar(gene_loss_virtual_freqs, total_freq_losses[cohort],width=0.3,linewidth=0, facecolor='#ff7f00',label='loss')

    gene_frequency_axis.bar(gene_loss_virtual_freqs-0.3, total_null_freq_losses[cohort],width=0.3,linewidth=0, facecolor='0.7',label='de novo\nexpectation')

    #gene_frequency_axis.legend(loc='upper center',frameon=False,fontsize=5,numpoints=1,ncol=3,handlelength=1) 
    
    # Now calculate statistical tests 
    
    # First check if prevalence is higher than null distribution
    num_bootstraps = default_num_bootstraps
    observed_counts = total_freq_all_snps[cohort] #[:-1]
    sample_size = observed_counts.sum()
    null_counts = total_null_freq_all_snps[cohort] #[:-1]
    prevalence_weights = null_counts*1.0/null_counts.sum()

    observed_loglikelihood = multinomial_loglikelihood(observed_counts, sample_size, prevalence_weights)

    bootstrapped_loglikelihoods = []
    for bootstrap_idx in xrange(0,num_bootstraps):
    
        bootstrapped_counts = sample_multinomial(sample_size, prevalence_weights)
    
        bootstrapped_loglikelihoods.append(multinomial_loglikelihood(bootstrapped_counts, sample_size, prevalence_weights))
    
    bootstrapped_loglikelihoods = numpy.array(bootstrapped_loglikelihoods)

    p_value = ((bootstrapped_loglikelihoods<=observed_loglikelihood).sum()+1.0)/(len(bootstrapped_loglikelihoods)+1.0)
    
    output_strs.append( "Prevalence P-value: %g" % p_value )
    #output_strs.append( "Observed loglikelihood: %g" % observed_loglikelihood )

    # Then check whether SNV prevalence distribution is asymmetric
    observed_counts = total_freq_all_snps[cohort]
    sample_size = observed_counts.sum()
    
    null_counts = (observed_counts+observed_counts[::-1])/2.0
    prevalence_weights = null_counts*1.0/null_counts.sum()

    #observed_loglikelihood = multinomial_loglikelihood(observed_counts, sample_size, prevalence_weights)  
    observed_loglikelihood = reversal_loglikelihoods(observed_counts)
    
    #print snv_prevalence_count_map[cohort]

    observed_prevalences = []
    for sample_pair in snv_prevalence_map[cohort]:
        snv_prevalence_map[cohort][sample_pair] = numpy.array(snv_prevalence_map[cohort][sample_pair])
        if len(snv_prevalence_map[cohort][sample_pair]>100):
            snv_prevalence_map[cohort][sample_pair] = choice(snv_prevalence_map[cohort][sample_pair],100)
              
        observed_prevalences.extend(snv_prevalence_map[cohort][sample_pair])
    observed_prevalences = numpy.clip(observed_prevalences,0,1)
    
    observed_ks = symmetrized_ks_distance(observed_prevalences)    

    bootstrapped_loglikelihoods = []
    bootstrapped_kss = []
    for bootstrap_idx in xrange(0,num_bootstraps):
    
        # Old method: time-reverses each SNV independently. (not conservative)
        # bootstrapped_counts = sample_multinomial(sample_size, prevalence_weights)
    
        # New method: time reverses sample pair
        bootstrapped_counts = numpy.zeros_like(observed_counts)
        bootstrapped_prevalences = []
        for sample_pair in snv_prevalence_count_map[cohort]:
            
            if random()<0.5:
                # normal order
                bootstrapped_counts += snv_prevalence_count_map[cohort][sample_pair]
            else:
                # reverse
                bootstrapped_counts += snv_prevalence_count_map[cohort][sample_pair][::-1]
                
            if random()<0.5:
                bootstrapped_prevalences.extend( snv_prevalence_map[cohort][sample_pair] )
            else:
                bootstrapped_prevalences.extend( 1-snv_prevalence_map[cohort][sample_pair] )       
                
        bootstrapped_sample_size = bootstrapped_counts.sum()
        
        #bootstrapped_loglikelihood = multinomial_loglikelihood(bootstrapped_counts, bootstrapped_sample_size, prevalence_weights))
        bootstrapped_loglikelihood = reversal_loglikelihoods(bootstrapped_counts)
        bootstrapped_loglikelihoods.append(bootstrapped_loglikelihood)
        bootstrapped_kss.append( symmetrized_ks_distance(bootstrapped_prevalences) ) 

            
    bootstrapped_loglikelihoods = numpy.array(bootstrapped_loglikelihoods)
    bootstrapped_kss = numpy.array(bootstrapped_kss)
    
    p_value = ((bootstrapped_loglikelihoods<=observed_loglikelihood).sum()+1.0)/(len(bootstrapped_loglikelihoods)+1.0)
    
    ks_p_value = ((bootstrapped_kss>=observed_ks).sum()+1.0)/(len(bootstrapped_kss)+1.0)
    
    output_strs.append( "SNV time-reversal P-value: %g" % p_value)
    output_strs.append( "SNV ks P-value: %g" % ks_p_value)
    # Then check whether SNV prevalence distribution is asymmetric
    
    observed_counts = numpy.hstack([total_freq_losses[cohort], total_freq_gains[cohort][::-1]])
    sample_size = observed_counts.sum()
    
    null_counts = (observed_counts+observed_counts[::-1])*1.0/2.0
    prevalence_weights = null_counts*1.0/null_counts.sum()

    #observed_loglikelihood = multinomial_loglikelihood(observed_counts, sample_size, prevalence_weights)
    observed_loglikelihood = reversal_loglikelihoods(observed_counts)
    
    observed_gain_prevalences = []
    observed_loss_prevalences = []
    for sample_pair in gene_gain_count_map[cohort]:
            
        sample_gain_prevalences = gene_gain_prevalence_map[cohort][sample_pair]
        sample_loss_prevalences = gene_loss_prevalence_map[cohort][sample_pair]

        observed_gain_prevalences.extend( sample_gain_prevalences )
        observed_loss_prevalences.extend( sample_loss_prevalences )
    
    if len(observed_gain_prevalences)==0:
        observed_gain_prevalences.append(0.5)
    if len(observed_loss_prevalences)==0:
        observed_loss_prevalences.append(0.5)
        
    observed_ks, dummy = ks_2samp(observed_gain_prevalences, observed_loss_prevalences)
        
            
    bootstrapped_loglikelihoods = []
    bootstrapped_kss = []
    other_bootstrapped_kss = [] # bootstrapped using other method, useful for checking donating but not invasion of strains
    for bootstrap_idx in xrange(0,num_bootstraps):
    
        # Old method: time-reverses each gene independently. (not conservative)
        #bootstrapped_counts = sample_multinomial(sample_size, prevalence_weights)
        
        # New method: time reverses sample pair
        bootstrapped_counts = numpy.zeros_like(observed_counts)
        bootstrapped_gain_prevalences = []
        bootstrapped_loss_prevalences = []
        other_bootstrapped_gain_prevalences = []
        other_bootstrapped_loss_prevalences = []
        
        for sample_pair in gene_gain_count_map[cohort]:
            
            bootstrapped_gains = gene_gain_count_map[cohort][sample_pair]
            bootstrapped_losses = gene_loss_count_map[cohort][sample_pair]
            
            sample_gain_prevalences = gene_gain_prevalence_map[cohort][sample_pair]
            sample_loss_prevalences = gene_loss_prevalence_map[cohort][sample_pair]
            
            if random()<0.5:
                bootstrapped_gains, bootstrapped_losses = bootstrapped_losses, bootstrapped_gains
                sample_gain_prevalences, sample_loss_prevalences = sample_loss_prevalences, sample_gain_prevalences
                    
            bootstrapped_counts += numpy.hstack([bootstrapped_losses, bootstrapped_gains[::-1]])   
            
            bootstrapped_gain_prevalences.extend( sample_gain_prevalences )
            bootstrapped_loss_prevalences.extend( sample_loss_prevalences )
            
            # In other bootstrap method, flip gains and losses independently of each other
            # (ok null for read donating, bad null for invasion of strains)
            if random()<0.5:
                other_bootstrapped_gain_prevalences.extend( sample_gain_prevalences )
            else:
                other_bootstrapped_loss_prevalences.extend( sample_gain_prevalences ) 
                
            if random()<0.5:
                other_bootstrapped_gain_prevalences.extend( sample_loss_prevalences )
            else:
                other_bootstrapped_loss_prevalences.extend( sample_loss_prevalences )   
                 
        bootstrapped_sample_size = bootstrapped_counts.sum()
        
        #print "O:", observed_counts
        #print "B:", bootstrapped_counts
        
        
        #bootstrapped_loglikelihood = multinomial_loglikelihood(bootstrapped_counts, bootstrapped_sample_size, prevalence_weights))
        bootstrapped_loglikelihood = reversal_loglikelihoods(bootstrapped_counts)
        
        bootstrapped_loglikelihoods.append(bootstrapped_loglikelihood)
        
        if (len(bootstrapped_gain_prevalences)==0) or (len(bootstrapped_loss_prevalences)==0):
            continue
        
        bootstrapped_ks, dummy = ks_2samp(bootstrapped_gain_prevalences, bootstrapped_loss_prevalences)
        
        bootstrapped_kss.append(bootstrapped_ks)
        
        other_bootstrapped_ks, dummy = ks_2samp(other_bootstrapped_gain_prevalences, other_bootstrapped_loss_prevalences)
        
        other_bootstrapped_kss.append(other_bootstrapped_ks)
        
        
    bootstrapped_loglikelihoods = numpy.array(bootstrapped_loglikelihoods)
    bootstrapped_kss = numpy.array(bootstrapped_kss)
    other_bootstrapped_kss = numpy.array(other_bootstrapped_kss)
    
    
    p_value = ((bootstrapped_loglikelihoods<=observed_loglikelihood).sum()+1.0)/(len(bootstrapped_loglikelihoods)+1.0)
    
    ks_p_value = ((bootstrapped_kss>=observed_ks).sum()+1.0)/(len(bootstrapped_kss)+1.0)
    other_ks_p_value = ((other_bootstrapped_kss>=observed_ks).sum()+1.0)/(len(other_bootstrapped_kss)+1.0)


    output_strs.append( "gene time-reversal P-value: %g" % p_value)
    output_strs.append( "gene ks P-value: %g" % ks_p_value)
    output_strs.append( "other gene ks P-value (gains and losses reversed independently): %g" % other_ks_p_value)
    

    # Now do dN/dSs

    non_ns = total_freq_snps[cohort]['1D']
    syn_ns = total_freq_snps[cohort]['4D']
    total_ns = non_ns + syn_ns
    ps = non_ns*1.0/total_ns

    non_opportunities = total_null_freq_snps[cohort]['1D']
    syn_opportunities = total_null_freq_snps[cohort]['4D']
    total_opportunities = non_opportunities + syn_opportunities

    observed_dNdSs = non_ns*1.0/syn_ns / (non_opportunities*1.0/syn_opportunities)

    bootstrapped_dNdSs = []

    for bootstrap_idx in xrange(0,num_bootstraps):

        bootstrapped_non_ns = sample_binomial(total_ns, ps)
        bootstrapped_syn_ns = total_ns - bootstrapped_non_ns
   
        bootstrapped_dNdSs.append( bootstrapped_non_ns*1.0/bootstrapped_syn_ns / (non_opportunities*1.0/syn_opportunities) )

    bootstrapped_dNdSs = numpy.sort( bootstrapped_dNdSs, axis=0 )  
   
    output_strs.append( "Observed dNdS " + str(observed_dNdSs) )
    output_strs.append( "Lower CI " + str(bootstrapped_dNdSs[long(0.025*num_bootstraps),:]))
    output_strs.append( "Upper CI " + str(bootstrapped_dNdSs[long(0.975*num_bootstraps),:]))

    # Now do pooled version
    pooled_p = non_ns.sum()*1.0/(total_ns.sum())
    pooled_ps = numpy.ones_like(total_ns)*pooled_p

    observed_loglikelihood = binomial_loglikelihoods(non_ns, total_ns, pooled_ps)

    # 
    for variant_type in variant_type_prevalence_map[cohort]:
        if len(variant_type_prevalence_map[cohort][variant_type])>1000:
            # downsample (For speed). Only used for twins
            variant_type_prevalence_map[cohort][variant_type] = choice(variant_type_prevalence_map[cohort][variant_type], 1000)
        
        variant_type_prevalence_map[cohort][variant_type] = numpy.clip(variant_type_prevalence_map[cohort][variant_type],0,1)
            
    
    observed_ks = ks_distance(variant_type_prevalence_map[cohort]['1D'], variant_type_prevalence_map[cohort]['4D']) 


    bootstrapped_loglikelihoods = []
    bootstrapped_kss = []
    for bootstrap_idx in xrange(0,num_bootstraps):
    
        bootstrapped_non_ns = sample_binomial(total_ns, pooled_ps)
    
        bootstrapped_loglikelihoods.append(binomial_loglikelihoods(bootstrapped_non_ns, total_ns, pooled_ps))
        
        nonsynonymous_idxs = numpy.hstack([numpy.ones_like(variant_type_prevalence_map[cohort]['1D']), numpy.zeros_like(variant_type_prevalence_map[cohort]['4D'])])
        
        prevalences = numpy.hstack([variant_type_prevalence_map[cohort]['1D'], variant_type_prevalence_map[cohort]['4D']])
        shuffle(prevalences)
        
        bootstrapped_nonsynonymous_prevalences = prevalences[nonsynonymous_idxs>0.5]
        bootstrapped_synonymous_prevalences = prevalences[nonsynonymous_idxs<0.5]
        bootstrapped_kss.append( ks_distance(bootstrapped_nonsynonymous_prevalences, bootstrapped_synonymous_prevalences) )
    
    bootstrapped_loglikelihoods = numpy.array(bootstrapped_loglikelihoods)
    bootstrapped_kss = numpy.array(bootstrapped_kss)
    
    
    p_value = ((bootstrapped_loglikelihoods<=observed_loglikelihood).sum()+1.0)/(len(bootstrapped_loglikelihoods)+1.0)
    ks_p_value = ((bootstrapped_kss>=observed_ks).sum()+1.0)/(len(bootstrapped_kss)+1.0)
    

    output_strs.append( "dNdS uniformity P-value: %g" % p_value)
    output_strs.append( "dNdS ks P-value: %g" % ks_p_value)
    
    for i in xrange(1,len(total_freq_all_snps[cohort])/2+1):
    
        non_counts = total_freq_snps[cohort]['1D'][:i].sum() #+total_freq_snps[cohort]['1D'][-i:].sum()
        syn_counts = total_freq_snps[cohort]['4D'][:i].sum() #+total_freq_snps[cohort]['4D'][-i:].sum()
        total_counts = non_counts + syn_counts
        p = non_counts*1.0/(total_counts+(total_counts==0))
    
    
        non_opportunities = total_null_freq_snps[cohort]['1D'][:i].sum() #+ total_null_freq_snps[cohort]['1D'][-i:].sum() 
        syn_opportunities = total_null_freq_snps[cohort]['4D'][:i].sum() #+ total_null_freq_snps[cohort]['4D'][-i:].sum()
        total_opportunities = non_opportunities+syn_opportunities
    
        observed_dNdS = (non_counts*1.0/syn_counts) / (non_opportunities*1.0/syn_opportunities)
    
        bootstrapped_non_ns = sample_binomial(total_counts,p,num_bootstraps)
        bootstrapped_syn_ns = total_counts - bootstrapped_non_ns
    
        bootstrapped_dNdSs = (bootstrapped_non_ns*1.0/bootstrapped_syn_ns) / (non_opportunities*1.0/syn_opportunities)
    
        bootstrapped_dNdSs = numpy.sort(bootstrapped_dNdSs)
    
        output_strs.append("%g %g %g %g %g %g" % (derived_freq_bins[i], non_counts, syn_counts, observed_dNdS, bootstrapped_dNdSs[long(0.025*num_bootstraps)], bootstrapped_dNdSs[long(0.975*num_bootstraps)]))
 

    # Now let's go the opposite direction   
    for i in xrange(0,len(total_freq_all_snps[cohort])/2+1):
    
        n = len(total_freq_all_snps[cohort])
    
        non_counts = total_freq_snps[cohort]['1D'][i:n-i].sum()
        syn_counts = total_freq_snps[cohort]['4D'][i:n-i].sum()
        total_counts = non_counts + syn_counts
        p = non_counts*1.0/(total_counts+(total_counts==0))
    
    
        non_opportunities = total_null_freq_snps[cohort]['1D'][i:n-i].sum()
        syn_opportunities = total_null_freq_snps[cohort]['4D'][i:n-i].sum()
        total_opportunities = non_opportunities+syn_opportunities
    
        observed_dNdS = (non_counts*1.0/syn_counts) / (non_opportunities*1.0/syn_opportunities)
    
        bootstrapped_non_ns = sample_binomial(total_counts,p,num_bootstraps)
        bootstrapped_syn_ns = total_counts - bootstrapped_non_ns
    
        bootstrapped_dNdSs = (bootstrapped_non_ns*1.0/bootstrapped_syn_ns) / (non_opportunities*1.0/syn_opportunities)
    
        bootstrapped_dNdSs = numpy.sort(bootstrapped_dNdSs)
    
        output_strs.append("%g %g %g %g %g %g" % (derived_freq_bins[i], non_counts, syn_counts, observed_dNdS, bootstrapped_dNdSs[long(0.025*num_bootstraps)], bootstrapped_dNdSs[long(0.975*num_bootstraps)]))
    
    output_strs.append("")
    output_strs.append("")
     
     
    if cohort=='hmp':
        # Plot ks figures
         
        # SNV one
        xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(observed_prevalences, min_x=-0.01, max_x=1.01, min_p=0)
        snv_ks_axis.step(xs,ns*1.0/ns[0],'-',linewidth=1, where='pre',zorder=4,color='#08519c',label='Observed')
        
          
        symmetrized_prevalences = numpy.hstack([observed_prevalences, 1-observed_prevalences])
        xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(symmetrized_prevalences, min_x=-0.01, max_x=1.01, min_p=0)
        snv_ks_axis.step(xs,ns*1.0/ns[0],'-',linewidth=0.5, color='0.7', where='pre',zorder=3,label='Time-reversal\nsymmetric')
        
        snv_ks_axis.legend(loc='upper right',frameon=False,fontsize=5,numpoints=1,ncol=1,handlelength=1)   
        
        # dN/dS one
        xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(variant_type_prevalence_map[cohort]['1D'], min_x=-0.01, max_x=1.01)
        dnds_ks_axis.step(xs,ns*1.0/ns[0],'-',linewidth=1, where='pre',zorder=4,color='#ff7f00',label='non (1D)')
        
        xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(variant_type_prevalence_map[cohort]['4D'], min_x=-0.01, max_x=1.01)
        dnds_ks_axis.step(xs,ns*1.0/ns[0],'-',linewidth=1, where='pre',zorder=3,color='#b3de69',label='syn (4D)')

        dnds_ks_axis.legend(loc='upper right',frameon=False,fontsize=5,numpoints=1,ncol=1,handlelength=1)   
       

        # gene one

        xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(observed_loss_prevalences, min_x=-0.01, max_x=1.01)
        gene_ks_axis.step(xs,ns*1.0/ns[0],'-',linewidth=1, where='pre',zorder=3,color='#ff7f00',label='loss')
 
        xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(observed_gain_prevalences, min_x=-0.01, max_x=1.01)
        
        gene_ks_axis.step(xs,ns*1.0/ns[0],'-',linewidth=1, where='pre',zorder=4,color='#b3de69',label='gain')
        
        gene_ks_axis.legend(loc='upper right',frameon=False,fontsize=5,numpoints=1,ncol=1,handlelength=1)   
       
        
        
            

        

### Now plot temporal qp figures
species_names = hmp_species_qp_counts.keys()

species_names = list(sorted(species_names, key=lambda s: sum(hmp_species_qp_counts[s])))

ys = numpy.arange(0,len(species_names))

yticklabels = []

for y,species_name in zip(ys,species_names):

    yticklabels.append(species_name)
    

    total_samples = sum(hmp_species_qp_counts[species_name])
    
    if total_samples>0:
    
        qp_samples = hmp_species_qp_counts[species_name][1]
        non_qp_samples = hmp_species_qp_counts[species_name][2]
        mixed_samples = hmp_species_qp_counts[species_name][3]
        dropout_samples = hmp_species_qp_counts[species_name][0]
        
        hmp_haploid_axis.barh([y],[qp_samples],linewidth=0, color='#08519c')
        
        hmp_haploid_axis.barh([y],[non_qp_samples], left=[qp_samples], linewidth=0, color='#de2d26')
        
        hmp_haploid_axis.barh([y],[mixed_samples], left=[qp_samples+non_qp_samples], linewidth=0, color='#8856a7')
        
        hmp_haploid_axis.barh([y],[dropout_samples], left=[qp_samples+non_qp_samples+mixed_samples], linewidth=0, color='0.7')
        
    total_samples = sum(twin_species_qp_counts[species_name])
    
    if total_samples>0:
    
        qp_samples = twin_species_qp_counts[species_name][1]
        non_qp_samples = twin_species_qp_counts[species_name][2]
        mixed_samples = twin_species_qp_counts[species_name][3]
        dropout_samples = twin_species_qp_counts[species_name][0]
        
        twin_haploid_axis.barh([y],[qp_samples],linewidth=0,label='QP->QP', color='#08519c')
        
        twin_haploid_axis.barh([y],[non_qp_samples], left=[qp_samples], linewidth=0, color='#de2d26')
        
        twin_haploid_axis.barh([y],[mixed_samples], left=[qp_samples+non_qp_samples], linewidth=0, color='#8856a7')
        
        twin_haploid_axis.barh([y],[dropout_samples], left=[qp_samples+non_qp_samples+mixed_samples], linewidth=0, color='0.7')
        
hmp_haploid_axis.yaxis.tick_left()
hmp_haploid_axis.xaxis.tick_bottom()
  
twin_haploid_axis.yaxis.tick_left()
twin_haploid_axis.xaxis.tick_bottom()

hmp_haploid_axis.set_yticks(ys+0.5)
twin_haploid_axis.set_yticks(ys+0.5)
hmp_haploid_axis.set_ylim([-1,len(ys)])
twin_haploid_axis.set_ylim([-1,len(ys)])

hmp_haploid_axis.set_yticklabels(yticklabels,fontsize=5)
twin_haploid_axis.set_yticklabels([])

hmp_haploid_axis.tick_params(axis='y', direction='out',length=3,pad=1)
twin_haploid_axis.tick_params(axis='y', direction='out',length=3,pad=1)

hmp_haploid_axis.set_xlim([0,200])
twin_haploid_axis.set_xlim([0,200])


### Do stuff for legend
hmp_haploid_axis.barh([-10],[1],linewidth=0,label='QP->QP', color='#08519c')
hmp_haploid_axis.barh([-10],[1],linewidth=0,label='non->non', color='#de2d26')
hmp_haploid_axis.barh([-10],[1],linewidth=0,label='mixed', color='#8856a7')
hmp_haploid_axis.barh([-10],[1],linewidth=0,label='dropout', color='0.7')
hmp_haploid_axis.legend(loc='lower right',frameon=False)

### Calculate KS test between twins and young twins

observed_twin_distribution = pooled_snp_change_distribution['twins']
observed_young_twin_distribution = pooled_snp_change_distribution['young_twins']

joint_distribution = numpy.hstack([observed_young_twin_distribution, observed_twin_distribution])

observed_ks, dummy = ks_2samp(observed_twin_distribution, observed_young_twin_distribution)

print observed_ks, dummy

bootstrapped_kss = []
num_bootstraps = default_num_bootstraps
for bootstrap_idx in xrange(0,num_bootstraps):

    shuffle(joint_distribution)

    bootstrapped_young_twin_distribution = joint_distribution[0:len(observed_young_twin_distribution)]
    
    bootstrapped_twin_distribution = joint_distribution[len(observed_young_twin_distribution):]
    
    bootstrapped_ks, dummy = ks_2samp(bootstrapped_twin_distribution, bootstrapped_young_twin_distribution)

    bootstrapped_kss.append(bootstrapped_ks)
    
bootstrapped_kss = numpy.array(bootstrapped_kss)

pvalue = ((bootstrapped_kss>=observed_ks).sum()+1.0)/(len(bootstrapped_kss)+1.0)

output_strs.append("Twin/Young twin KS test: Pvalue=%g" % pvalue)

output_file = open(output_filename,"w")
output_file.write("\n".join(output_strs))
output_file.write("\n")
output_file.close()

sys.stderr.write("Saving figures...\t")
fig2.savefig('%s/figure_5.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')
fig.savefig('%s/supplemental_within_across_species.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')
sys.stderr.write("Done!\n")
fig3.savefig('%s/supplemental_young_twin_within.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight',transparent=True)
fig4.savefig('%s/supplemental_within_between_avg.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight',transparent=True)
fig5.savefig('%s/supplemental_twin_modification_frequency.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight',transparent=True)
#fig6.savefig('%s/supplemental_temporal_qp_sample_size.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight',transparent=True)
fig8.savefig('%s/supplemental_within_across_species.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight',transparent=True)
fig9.savefig('%s/supplemental_within_ks.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight',transparent=True)
fig10.savefig('%s/supplemental_sites_retained.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight',transparent=True)



    
