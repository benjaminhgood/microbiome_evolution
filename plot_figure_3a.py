import matplotlib  
matplotlib.use('Agg') 
import config
import parse_midas_data
import parse_HMP_data
import os.path
import pylab
import sys
import numpy
from numpy.random import choice

import species_phylogeny_utils
import diversity_utils
import gene_diversity_utils
import calculate_substitution_rates

import stats_utils
import matplotlib.colors as colors
import matplotlib.cm as cmx
from math import log10,ceil
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy.random import randint

from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster.hierarchy import cophenet
from scipy.cluster.hierarchy import fcluster

from scipy.stats import gaussian_kde

mpl.rcParams['font.size'] = 6
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

num_bootstraps = 1000

min_coverage = config.min_median_coverage
alpha = 0.5 # Confidence interval range for rate estimates
low_pi_threshold = 1e-03
low_divergence_threshold = 2e-04
min_change = 0.8
min_sample_size = 33 # 46 gives at least 1000 pairs
allowed_variant_types = set(['1D','2D','3D','4D'])

divergence_matrices = {}
low_divergence_pair_counts = {}
null_low_divergence_pair_counts = [{} for i in xrange(0,num_bootstraps)]

good_species_list = parse_midas_data.parse_good_species_list()

sys.stderr.write("Loading sample metadata...\n")
subject_sample_map = parse_HMP_data.parse_subject_sample_map()
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
    snp_samples = snp_samples[parse_midas_data.calculate_unique_samples(subject_sample_map, sample_list=snp_samples)]


    if len(snp_samples) < min_sample_size:
        sys.stderr.write("Not enough unique samples!\n")
        continue
        
    # Load divergence matrices 
    sys.stderr.write("Loading pre-computed substitution rates for %s...\n" % species_name)
    substitution_rate_map = calculate_substitution_rates.load_substitution_rate_map(species_name)
    sys.stderr.write("Calculating matrix...\n")
    dummy_samples, snp_difference_matrix, snp_opportunity_matrix = calculate_substitution_rates.calculate_matrices_from_substitution_rate_map(substitution_rate_map, 'core', allowed_samples=snp_samples)
    snp_samples = dummy_samples
    sys.stderr.write("Done!\n")

    snp_substitution_matrix = snp_difference_matrix*1.0/(snp_opportunity_matrix+(snp_opportunity_matrix==0))
    
    # Find closely related samples
    for i in xrange(0, snp_substitution_matrix.shape[0]):
        for j in xrange(i+1, snp_substitution_matrix.shape[0]):
            
            if snp_substitution_matrix[i,j] >= 0:
                
                if snp_substitution_matrix[i,j] <= low_divergence_threshold:
                
                    sample_pair = frozenset([snp_samples[i],snp_samples[j]])
                    
                    if sample_pair not in low_divergence_pair_counts:
                        low_divergence_pair_counts[sample_pair] = 0
                        
                    low_divergence_pair_counts[sample_pair] += 1
                    
                    for bootstrap_idx in xrange(0,num_bootstraps):
                        
                        # now draw null pair from 
                        null_samples = choice(snp_samples,size=2,replace=False)
                        null_pair = frozenset([null_samples[0], null_samples[1]])
                    
                        if null_pair not in null_low_divergence_pair_counts[bootstrap_idx]:
                            null_low_divergence_pair_counts[bootstrap_idx][null_pair] = 0
                        
                        null_low_divergence_pair_counts[bootstrap_idx][null_pair] += 1
                        
    
    divergence_matrices[species_name] = snp_substitution_matrix
 
# # low divergence strains across species
# # samples
        
species_names = []
sample_sizes = []

for species_name in species_phylogeny_utils.sort_phylogenetically(divergence_matrices.keys()):
    species_names.append(species_name)
    sample_sizes.append( divergence_matrices[species_name].shape[0] )
    
# sort in descending order of sample size
# Sort by num haploids    
#sample_sizes, species_names = zip(*sorted(zip(sample_sizes, species_names),reverse=True))
    
sys.stderr.write("Postprocessing %d species...\n" % len(species_names))
        

####################################################
#
# Set up Figure (1 panels, arranged in 1x1 grid)
#
####################################################

haploid_color = '#08519c'

pylab.figure(1,figsize=(4,1))
fig = pylab.gcf()
# make three panels panels
outer_grid  = gridspec.GridSpec(1,1)

divergence_axis = plt.Subplot(fig, outer_grid[0])
fig.add_subplot(divergence_axis)

divergence_axis.set_ylabel('Divergence, $d$')
divergence_axis.set_ylim([1e-06,1e-01])
divergence_axis.set_xlim([-1,len(species_names)])

xticks = numpy.arange(0,len(species_names))
xticklabels = ["%s (%d)" % (species_names[i],sample_sizes[i]) for i in xrange(0,len(sample_sizes))]
#xticklabels = ["%s" % (species_names[i]) for i in xrange(0,len(sample_sizes))]

divergence_axis.set_xticks(xticks)
divergence_axis.set_xticklabels(xticklabels, rotation='vertical',fontsize=4)

divergence_axis.spines['top'].set_visible(False)
divergence_axis.spines['right'].set_visible(False)
divergence_axis.get_xaxis().tick_bottom()
divergence_axis.get_yaxis().tick_left()


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
    
    log_divergences = numpy.log(divergences)
    
    kernel = gaussian_kde(log_divergences)
    
    n = len(divergences)
    
    
    
    #percentiles = numpy.array([0.001,0.01,0.1,0.25,0.5,0.75,0.9,0.99,0.999])
    percentiles = numpy.array([0.001,0.01,0.5])
    
    quantiles = numpy.array([divergences[long(n*p)] for p in percentiles])
    quantiles = numpy.clip(quantiles,1e-06,1)
    
    theory_log_divergences = numpy.linspace(log_divergences.min(), log_divergences.max()+1,100)
    theory_divergences = numpy.exp(theory_log_divergences)
    theory_pdf = kernel(theory_log_divergences)
    theory_pdf = theory_pdf / theory_pdf.max() * 0.45
    
    
    divergence_axis.fill_betweenx(theory_divergences, species_idx-theory_pdf, species_idx+theory_pdf,linewidth=0,facecolor='0.7') #facecolor='#de2d26') # red color
    #divergence_axis.semilogy(numpy.ones_like(quantiles)*species_idx, quantiles,'k_',markersize=3)
    
    # Median
    divergence_axis.semilogy([species_idx], [quantiles[-1]],'_',markersize=3,color='#de2d26')
    # 1%-tile
    divergence_axis.semilogy([species_idx], [quantiles[1]],'.',markersize=2.5,color='#de2d26',markeredgewidth=0)
    # 0.1%-tile
    divergence_axis.semilogy([species_idx], [quantiles[0]],'.',markersize=4,color='#de2d26',markeredgewidth=0)
    # Line connecting them
    divergence_axis.semilogy([species_idx,species_idx], [quantiles[0],quantiles[-1]],'-',color='#de2d26')

# Now plot distribution of closely related species per pair
pylab.figure(2,figsize=(1.5,1))
fig2 = pylab.gcf()
# make three panels panels
outer_grid  = gridspec.GridSpec(1,1)

histogram_axis = plt.Subplot(fig2, outer_grid[0])
fig2.add_subplot(histogram_axis)

histogram_axis.set_ylabel('Sample pairs')
histogram_axis.set_xlabel('Low divergence strains')
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

    null_histogram += bootstrapped_histogram/num_bootstraps

    pvalue += (bootstrapped_histogram[1:].sum() >= observed_histogram[1:].sum())

pvalue = (pvalue+1)/(num_bootstraps+1.0)

print "pvalue for closely related pair distribution =", pvalue

# Plot histograms

histogram_axis.bar(ks-0.3, observed_histogram, width=0.3, linewidth=0, color='r',label='Observed',bottom=1e-03)

histogram_axis.bar(ks, null_histogram, width=0.3, linewidth=0, color='0.7',label='Null',bottom=1e-03)

histogram_axis.semilogy([1e-03,1e-03],[0,4],'k-')
histogram_axis.set_ylim([1e-01,1e03])
histogram_axis.set_xticks([1,2,3])
histogram_axis.legend(loc='upper right',frameon=False)

sys.stderr.write("Saving figure...\t")
fig.savefig('%s/figure_3a.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')
fig2.savefig('%s/supplemental_low_divergence_pair_distribution.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')

sys.stderr.write("Done!\n")

 