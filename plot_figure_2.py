import matplotlib  
matplotlib.use('Agg') 
import sample_utils
import config
import parse_midas_data
import os.path
import pylab
import sys
import numpy
from numpy.random import choice

import species_phylogeny_utils
import diversity_utils
import gene_diversity_utils
import calculate_substitution_rates
import figure_utils

import stats_utils
import matplotlib.colors as colors
import matplotlib.cm as cmx
from math import log10,ceil,log,exp
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy.random import randint
from mpl_toolkits.axes_grid.inset_locator import inset_axes

from numpy.random import randint, binomial, choice
from scipy.stats import poisson as poisson_distribution
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster.hierarchy import cophenet
from scipy.cluster.hierarchy import fcluster

from scipy.stats import gaussian_kde

mpl.rcParams['font.size'] = 5
mpl.rcParams['axes.labelpad'] = 2
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

num_bootstraps = 10000

min_coverage = config.min_median_coverage
low_divergence_threshold = config.between_low_divergence_threshold
min_sample_size = config.between_host_min_sample_size # 46 gives at least 1000 pairs, 33 gives at least 500 (actually 528)
allowed_variant_types = set(['1D','2D','3D','4D'])

output_filename = ('%s/between_host_output.txt' % (parse_midas_data.analysis_directory))
output_strs = []

divergence_matrices = {}
low_divergence_pair_counts = {}
null_low_divergence_pair_counts = [{} for i in xrange(0,num_bootstraps)]

low_divergence_same_continent_counts = {True:0, False:0}
null_low_divergence_same_continent_counts = [{True:0, False:0} for i in xrange(0,num_bootstraps)]


low_divergence_snp_differences = []
low_divergence_gene_differences = []
low_divergence_clock_null_gene_differences = []
normal_divergence_gene_differences = []

good_species_list = parse_midas_data.parse_good_species_list()
if debug:
    good_species_list = ["Bacteroides_vulgatus_57955", "Bacteroides_uniformis_57318"]

sys.stderr.write("Loading sample metadata...\n")
subject_sample_map = sample_utils.parse_subject_sample_map()
sample_continent_map = sample_utils.parse_sample_continent_map()
sys.stderr.write("Done!\n")

for species_name in good_species_list:

    sys.stderr.write("Loading haploid samples...\n")
    # Only plot samples above a certain depth threshold that are "haploids"
    snp_samples = diversity_utils.calculate_haploid_samples(species_name, debug=debug)
    
    if len(snp_samples) < min_sample_size:
        sys.stderr.write("Not enough haploid samples!\n")
        continue
        
    sys.stderr.write("Calculating unique samples...\n")
    # Only consider one sample per person
    snp_samples = snp_samples[sample_utils.calculate_unique_samples(subject_sample_map, sample_list=snp_samples)]


    if len(snp_samples) < min_sample_size:
        sys.stderr.write("Not enough unique samples!\n")
        continue
        
    # Load divergence matrices 
    sys.stderr.write("Loading pre-computed substitution rates for %s...\n" % species_name)
    substitution_rate_map = calculate_substitution_rates.load_substitution_rate_map(species_name)
    sys.stderr.write("Calculating matrix...\n")
    dummy_samples, snp_difference_matrix, snp_opportunity_matrix = calculate_substitution_rates.calculate_matrices_from_substitution_rate_map(substitution_rate_map, 'core', allowed_samples=snp_samples)
    snp_samples = dummy_samples
    
    dummy_samples, gene_difference_matrix, gene_opportunity_matrix = calculate_substitution_rates.calculate_matrices_from_substitution_rate_map(substitution_rate_map, 'genes', allowed_samples=snp_samples)
    
    sys.stderr.write("Done!\n")

    snp_substitution_matrix = snp_difference_matrix*1.0/(snp_opportunity_matrix+(snp_opportunity_matrix==0))
    
    gene_differences = []
    for i in xrange(0, gene_difference_matrix.shape[0]):
        for j in xrange(i+1, gene_difference_matrix.shape[0]):
            if gene_opportunity_matrix[i,j]>0.5:
                gene_differences.append(gene_difference_matrix[i,j])
                
    gene_differences = numpy.array(gene_differences)
    
    
    snp_substitution_rates = []
    for i in xrange(0, snp_opportunity_matrix.shape[0]):
        for j in xrange(i+1, snp_opportunity_matrix.shape[0]):
            if snp_opportunity_matrix[i,j]>0.5:
                snp_substitution_rates.append(snp_substitution_matrix[i,j]) 
    snp_substitution_rates = numpy.array(snp_substitution_rates)
    
    median_gene_difference = numpy.median(gene_differences)
    
    scaled_gene_differences = gene_differences*1.0/median_gene_difference
    
    scaled_gene_difference_matrix = gene_difference_matrix*1.0/median_gene_difference
    
    # Find closely related samples
    for i in xrange(0, snp_substitution_matrix.shape[0]):
        for j in xrange(i+1, snp_substitution_matrix.shape[0]):
            
            if snp_opportunity_matrix[i,j] > 0.5:
                
                if snp_substitution_matrix[i,j] <= low_divergence_threshold:
                
                    low_divergence_snp_differences.append(snp_difference_matrix[i,j])
                
                    if gene_opportunity_matrix[i,j] > 0.5:
                        low_divergence_gene_differences.append(gene_difference_matrix[i,j])
                        
                        # get clock-like null
                        low_divergence_clock_null_gene_differences.append( numpy.median(gene_differences) * snp_substitution_matrix[i,j] / numpy.median(snp_substitution_rates) )
                        
                        
                        # Choosing it here ensures that the species contribute according to how
                        # many low divergence samples they have
                        normal_divergence_gene_differences.extend(choice(gene_differences,size=100))
                
                    sample_pair = frozenset([snp_samples[i],snp_samples[j]])
                    
                    if sample_pair not in low_divergence_pair_counts:
                        low_divergence_pair_counts[sample_pair] = 0
                        
                    low_divergence_pair_counts[sample_pair] += 1
                    
                    same_continent = (sample_continent_map[snp_samples[i]] == sample_continent_map[snp_samples[j]])
                    
                    low_divergence_same_continent_counts[same_continent] += 1
                    
                    
                    for bootstrap_idx in xrange(0,num_bootstraps):
                        
                        # now draw null pair from 
                        null_samples = choice(snp_samples,size=2,replace=False)
                        null_pair = frozenset([null_samples[0], null_samples[1]])
                    
                        if null_pair not in null_low_divergence_pair_counts[bootstrap_idx]:
                            null_low_divergence_pair_counts[bootstrap_idx][null_pair] = 0
                        
                        null_low_divergence_pair_counts[bootstrap_idx][null_pair] += 1
                        
                        same_continent = (sample_continent_map[null_samples[0]] == sample_continent_map[null_samples[1]])
                    
                        null_low_divergence_same_continent_counts[bootstrap_idx][same_continent] += 1
                    
                        
                 
                
    divergence_matrices[species_name] = snp_substitution_matrix
 
# # low divergence strains across species
# # samples
        
species_names = []
sample_sizes = []

for species_name in divergence_matrices.keys():
    species_names.append(species_name)
    
    if species_name=='Bacteroides_vulgatus_57955':
        sample_sizes.append(-1000 )
    else:
        sample_sizes.append( -divergence_matrices[species_name].shape[0] )
    
sorted_species_names = species_phylogeny_utils.sort_phylogenetically(divergence_matrices.keys(),first_entry='Bacteroides_vulgatus_57955',second_sorting_attribute=sample_sizes)

species_names = []
sample_sizes = []
for species_name in reversed(sorted_species_names):
    species_names.append(species_name)
    sample_sizes.append( divergence_matrices[species_name].shape[0] )
    


    
# sort in descending order of sample size
# Sort by num haploids    
#sample_sizes, species_names = zip(*sorted(zip(sample_sizes, species_names),reverse=True))
    
sys.stderr.write("Postprocessing %d species...\n" % len(species_names))
print "Analyzing %d species with %d or more QP samples" % (len(species_names), min_sample_size)
output_strs.append("Analyzing %d species with %d or more QP samples" % (len(species_names), min_sample_size))

####################################################
#
# Set up Figure (1 panels, arranged in 1x1 grid)
#
####################################################

haploid_color = '#08519c'

pylab.figure(1,figsize=(6, 3))
fig = pylab.gcf()
# make three panels panels
outer_grid  = gridspec.GridSpec(1,2, width_ratios=[1.2,1],wspace=0.2)

pylab.figure(2,figsize=(6, 2))
fig2 = pylab.gcf()
# make three panels panels
outer_grid_2  = gridspec.GridSpec(1,2,width_ratios=[1,1.2])


left_grid = gridspec.GridSpecFromSubplotSpec(1,2, width_ratios=[0.7,1],wspace=0,subplot_spec=outer_grid[1])

right_grid = gridspec.GridSpecFromSubplotSpec(2,1, height_ratios=[2,1],hspace=0.5,subplot_spec=outer_grid[0])

bottom_grid = gridspec.GridSpecFromSubplotSpec(1,2, width_ratios=[0.5, 0.7], wspace=0.6,subplot_spec=right_grid[1])

bottom_right_grid = gridspec.GridSpecFromSubplotSpec(1,2, width_ratios=[1,0.1],wspace=0.01,subplot_spec=outer_grid_2[1])

divergence_axis = plt.Subplot(fig, left_grid[0])
fig.add_subplot(divergence_axis)

divergence_axis.set_xlabel('Divergence, $d$')
divergence_axis.set_xlim([1e-06,1e-01])
divergence_axis.set_ylim([-1.5,len(species_names)+0.5])
#line, = divergence_axis.plot([low_divergence_threshold, low_divergence_threshold],[-1.5,len(species_names)+0.5],'-',color='k',linewidth=0.25)
#line.set_dashes((1,1))

divergence_axis.fill_between([1e-06,low_divergence_threshold],[-1.5,-1.5],[len(species_names)+0.5, len(species_names)+0.5],color='#fee0d2')
divergence_axis.text(exp((log(1e-06)+log(low_divergence_threshold))/2), len(species_names)+0.5 + (len(species_names)+2)*0.02, "'closely\nrelated'",fontsize=5,fontstyle='italic',ha='center',color='#fc9272')

# get better haploid species names
pretty_species_names = []
for species_name in species_names:

    base_name = figure_utils.get_pretty_species_name(species_name)
    pretty_name = base_name
    if pretty_name in pretty_species_names:
        idx = 1
        while pretty_name in pretty_species_names:
            idx += 1
            pretty_name = base_name + (" %d" % (idx))
    
    pretty_species_names.append(pretty_name)
    

yticks = numpy.arange(0,len(species_names))
yticklabels = ["%s, n=%d" % (pretty_species_names[i],sample_sizes[i]) for i in xrange(0,len(sample_sizes))]

divergence_axis.set_yticks(yticks)
divergence_axis.set_yticklabels(yticklabels,fontsize=4)
divergence_axis.tick_params(axis='y', direction='out',length=3,pad=1)


#divergence_axis.spines['top'].set_visible(False)
#divergence_axis.spines['right'].set_visible(False)
#divergence_axis.get_xaxis().tick_bottom()
divergence_axis.get_yaxis().tick_right()




# Plot percentiles of divergence distribution
for species_idx in xrange(0,len(species_names)):

    species_name = species_names[species_idx]
    
    sys.stderr.write("Postprocessing %s (%d samples)...\n" % (species_name, divergence_matrices[species_name].shape[0]))
    divergence_matrix = divergence_matrices[species_name]
    divergences = []
    for i in xrange(0, divergence_matrix.shape[0]):
        for j in xrange(i+1, divergence_matrix.shape[0]):
            
            if divergence_matrix[i,j] >= 0:
                
                divergences.append(divergence_matrix[i,j])
    
    divergences = numpy.array(divergences)
    divergences = numpy.clip(divergences,1e-06,1)
    divergences.sort() # ascending by default
    
    min_divergence = divergences[0]
    nextmin_divergence = divergences[1]
    
    log_divergences = numpy.log(divergences)
    
    kernel = gaussian_kde(log_divergences)
    
    n = len(divergences)
    
    
    
    #percentiles = numpy.array([0.001,0.01,0.1,0.25,0.5,0.75,0.9,0.99,0.999])
    percentiles = numpy.array([0.001,0.01,0.5])
    
    quantiles = numpy.array([divergences[long(n*p)] for p in percentiles])
    quantiles = numpy.clip(quantiles,1e-06,1)
    
    # Use second smallest value for robustness. 
    if quantiles[0]<nextmin_divergence:
        quantiles[0] = nextmin_divergence
    
    theory_log_divergences = numpy.linspace(log_divergences.min(), log_divergences.max()+1,100)
    theory_divergences = numpy.exp(theory_log_divergences)
    theory_pdf = kernel(theory_log_divergences)
    theory_pdf = theory_pdf / theory_pdf.max() * 0.45
    
    
    divergence_axis.fill_between(theory_divergences, species_idx-theory_pdf, species_idx+theory_pdf,linewidth=0,facecolor='0.3') #facecolor='#1f2c88') #,facecolor='0.7') #facecolor='#de2d26') # red color
    #divergence_axis.semilogy(numpy.ones_like(quantiles)*species_idx, quantiles,'k_',markersize=3)
    
    # Median
    divergence_axis.semilogx([quantiles[-1]], [species_idx], '|',markersize=3,color='#de2d26')
    # 1%-tile
    divergence_axis.semilogx([quantiles[1]], [species_idx], '.',markersize=2.5,color='#de2d26',markeredgewidth=0)
    # 0.1%-tile
    divergence_axis.semilogx([quantiles[0]], [species_idx], '.',markersize=4,color='#de2d26',markeredgewidth=0)
    # Line connecting them
    divergence_axis.semilogx([quantiles[0],quantiles[-1]],[species_idx,species_idx], '-',color='#de2d26')


histogram_axis = plt.Subplot(fig, bottom_grid[1])
fig.add_subplot(histogram_axis)

histogram_axis.set_ylabel('Host pairs')
histogram_axis.set_xlabel('# closely related\nstrains / pair')
histogram_axis.set_xlim([0,4])

histogram_axis.spines['top'].set_visible(False)
histogram_axis.spines['right'].set_visible(False)
histogram_axis.get_xaxis().tick_bottom()
histogram_axis.get_yaxis().tick_left()

ks = numpy.arange(1,4)

# Calculate observed histogram
low_divergence_counts = numpy.array(low_divergence_pair_counts.values())

observed_histogram = numpy.array([(low_divergence_counts==k).sum() for k in ks])*1.0

null_histogram = numpy.zeros_like(ks)*1.0

pvalue = 0

# Calculate null histogram
for bootstrap_idx in xrange(0,num_bootstraps):
    
    bootstrapped_low_divergence_counts = numpy.array(null_low_divergence_pair_counts[bootstrap_idx].values())
    
    bootstrapped_histogram = numpy.array([(bootstrapped_low_divergence_counts==k).sum() for k in ks])*1.0

    null_histogram += bootstrapped_histogram*1.0/num_bootstraps

    pvalue += (bootstrapped_histogram[1:].sum() >= observed_histogram[1:].sum())

pvalue = (pvalue+1)/(num_bootstraps+1.0)

output_strs.append("pvalue for closely related pair distribution = %g" % pvalue)

# Plot histograms

histogram_axis.bar(ks-0.3, observed_histogram, width=0.3, linewidth=0, color='r',label='Obs',bottom=1e-03)

histogram_axis.bar(ks, null_histogram, width=0.3, linewidth=0, color='0.7',label='Null',bottom=1e-03)

histogram_axis.semilogy([1e-03,1e-03],[0,4],'k-')
histogram_axis.set_ylim([3e-01,3e03])
histogram_axis.set_xticks([1,2,3])
histogram_axis.legend(loc='upper right',frameon=False,fontsize=4,numpoints=1, handlelength=1)

snp_difference_axis = plt.Subplot(fig2, outer_grid_2[0])
fig2.add_subplot(snp_difference_axis)

snp_difference_axis.set_xlabel('# SNV differences')
snp_difference_axis.set_ylabel('Fraction $\leq n$')

snp_difference_axis.spines['top'].set_visible(False)
snp_difference_axis.spines['right'].set_visible(False)
snp_difference_axis.get_xaxis().tick_bottom()
snp_difference_axis.get_yaxis().tick_left()

snp_difference_axis.semilogx([1,1])
snp_difference_axis.set_xlim([1,1e03])
snp_difference_axis.set_ylim([0,1.174])

gene_difference_axis = plt.Subplot(fig2, bottom_right_grid[0])
fig2.add_subplot(gene_difference_axis)

gene_difference_axis.set_xlabel('# gene differences')
gene_difference_axis.set_ylabel('Fraction $\leq n$')

gene_difference_axis.spines['top'].set_visible(False)
gene_difference_axis.spines['right'].set_visible(False)
gene_difference_axis.get_xaxis().tick_bottom()
gene_difference_axis.get_yaxis().tick_left()

gene_difference_axis.semilogx([1,1])
gene_difference_axis.set_xlim([1,1e04])
gene_difference_axis.set_ylim([0,1.174])

low_divergence_snp_differences = numpy.array(low_divergence_snp_differences)
low_divergence_gene_differences = numpy.array(low_divergence_gene_differences)
low_divergence_clock_null_gene_differences = numpy.array(low_divergence_clock_null_gene_differences)
normal_divergence_gene_differences = numpy.array(normal_divergence_gene_differences)

xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(low_divergence_gene_differences, min_x=0.1,max_x=1e04)
snp_difference_axis.step(xs,1-ns*1.0/ns[0],'r-',label='Closely\nrelated',zorder=1)


xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(low_divergence_gene_differences, min_x=0.1,max_x=1e04)
gene_difference_axis.step(xs,1-ns*1.0/ns[0],'r-',label='Closely\nrelated',zorder=1)

xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(normal_divergence_gene_differences,min_x=0.1,max_x=1e04)
gene_difference_axis.step(xs,1-ns*1.0/ns[0],'k-',label='All',zorder=2)

xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(low_divergence_clock_null_gene_differences, min_x=0.1,max_x=1e04)
gene_difference_axis.step(xs,1-ns*1.0/ns[0],'-', color='0.7', label='Scaled',zorder=0)





#gene_difference_axis.legend(loc=(0.01,0.92),frameon=False,fontsize=4, ncol=3, numpoints=1, handlelength=1)
gene_difference_axis.legend(loc=(0.9,0.15),frameon=False,fontsize=4, ncol=1, numpoints=1, handlelength=1)

output_strs.append("%d same continent closely related strains, %d different continent closely related strains" % (low_divergence_same_continent_counts[True], low_divergence_same_continent_counts[False]))

###
#
# Same / different continent distribution
#
###

# Calculate observed histogram
low_divergence_counts = numpy.array(low_divergence_pair_counts.values())

observed_same = low_divergence_same_continent_counts[True] 
observed_different = low_divergence_same_continent_counts[False]

null_same = 0
null_different = 0

pvalue = 0

# Calculate null histogram
for bootstrap_idx in xrange(0,num_bootstraps):
    
    null_same += null_low_divergence_same_continent_counts[bootstrap_idx][True]*1.0/num_bootstraps 
    
    null_different += null_low_divergence_same_continent_counts[bootstrap_idx][False]*1.0/num_bootstraps 
    
    
    pvalue += (null_low_divergence_same_continent_counts[bootstrap_idx][True] >= observed_same)

pvalue = (pvalue+1.0)/(num_bootstraps+1.0)

output_strs.append( "pvalue for closely related continent distribution = %g" %  pvalue)

output_strs.append("%g expected same continent closely related strains, %g different continent" % (null_same, null_different))

continent_axis = plt.Subplot(fig, bottom_grid[0])
fig.add_subplot(continent_axis)

continent_axis.set_ylabel('Closely related strains')
continent_axis.set_xlim([0,3])

continent_axis.spines['top'].set_visible(False)
continent_axis.spines['right'].set_visible(False)
continent_axis.get_xaxis().tick_bottom()
continent_axis.get_yaxis().tick_left()

continent_axis.bar([1-0.3, 2-0.3], [observed_same, observed_different], width=0.3, linewidth=0, color='r',label='Obs',bottom=1e-03)

continent_axis.bar([1, 2], [null_same, null_different], width=0.3, linewidth=0, color='0.7',label='Null',bottom=1e-03)

continent_axis.set_ylim([0,1.6e03])
continent_axis.set_xticks([1,2])
continent_axis.set_xticklabels(['Same\ncontinent','Diff\ncontinent'], rotation='vertical',fontsize=4)

continent_axis.legend(loc='upper right',frameon=False,fontsize=4,numpoints=1, handlelength=1)

output_file = open(output_filename,"w")
output_file.write("\n".join(output_strs))
output_file.write("\n")
output_file.close()

sys.stderr.write("Saving figure...\t")
fig.savefig('%s/figure_2.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')
fig2.savefig('%s/supplemental_closely_related_differences.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')


sys.stderr.write("Done!\n")

 