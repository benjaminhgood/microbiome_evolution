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
import calculate_substitution_rates

import stats_utils
import matplotlib.colors as colors
import matplotlib.cm as cmx
from math import log10,ceil
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy.random import randint, binomial
from scipy.stats import poisson as poisson_distribution
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

min_coverage = config.min_median_coverage
alpha = 0.5 # Confidence interval range for rate estimates
low_pi_threshold = 1e-03
low_divergence_threshold = 1e-03
min_change = 0.8
min_sample_size = 33 # 46 gives at least 1000 pairs
allowed_variant_types = set(['1D','2D','3D','4D'])

syn_differences = {}
syn_opportunities = {}
syn_pseudocounts = {}
non_differences = {}
non_pseudocounts = {}
non_opportunities = {}

good_species_list = parse_midas_data.parse_good_species_list()
if debug:
    good_species_list = good_species_list[0:2]

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
    sys.stderr.write("Calculating matrices...\n")
    dummy_samples, syn_difference_matrix, syn_opportunity_matrix = calculate_substitution_rates.calculate_matrices_from_substitution_rate_map(substitution_rate_map, '4D', allowed_samples=snp_samples)
    dummy_samples, non_difference_matrix, non_opportunity_matrix = calculate_substitution_rates.calculate_matrices_from_substitution_rate_map(substitution_rate_map, '1D', allowed_samples=snp_samples)
    snp_samples = dummy_samples
    
    syn_differences[species_name] = []
    syn_pseudocounts[species_name] = []
    syn_opportunities[species_name] = []
    
    non_differences[species_name] = []
    non_pseudocounts[species_name] = []
    non_opportunities[species_name] = []
    
    
    for i in xrange(0, syn_difference_matrix.shape[0]):
        for j in xrange(i+1, syn_difference_matrix.shape[0]):
            
            if syn_opportunity_matrix[i,j]>0 and non_opportunity_matrix[i,j]>0:
            
                syn_differences[species_name].append(syn_difference_matrix[i,j]+1)
                syn_pseudocounts[species_name].append(1)
                syn_opportunities[species_name].append(syn_opportunity_matrix[i,j])
                
                non_differences[species_name].append( non_difference_matrix[i,j] )
                non_pseudocounts[species_name].append( non_opportunity_matrix[i,j]*1.0/syn_opportunity_matrix[i,j] )
                non_opportunities[species_name].append(non_opportunity_matrix[i,j])
        
                
            
    syn_differences[species_name] = numpy.array(syn_differences[species_name])
    syn_pseudocounts[species_name] = numpy.array(syn_pseudocounts[species_name])
    
    syn_opportunities[species_name] = numpy.array(syn_opportunities[species_name])
    
    non_differences[species_name] = numpy.array(non_differences[species_name])
    non_pseudocounts[species_name] = numpy.array(non_pseudocounts[species_name])
    
    non_opportunities[species_name] = numpy.array(non_opportunities[species_name])    

species_names = []

for species_name in syn_differences.keys():
    species_names.append(species_name)
    
        
sys.stderr.write("Postprocessing %d species...\n" % len(species_names))
        

####################################################
#
# Set up Figure (1 panels, arranged in 1x1 grid)
#
####################################################

haploid_color = '#08519c'

pylab.figure(1,figsize=(3.42,1.7))
fig = pylab.gcf()
# make three panels panels
outer_grid  = gridspec.GridSpec(1,1)

divergence_axis = plt.Subplot(fig, outer_grid[0])
fig.add_subplot(divergence_axis)

divergence_axis.set_ylabel('$d_N/d_S$')
divergence_axis.set_xlabel('Synonymous divergence, $d_S$')

divergence_axis.spines['top'].set_visible(False)
divergence_axis.spines['right'].set_visible(False)
divergence_axis.get_xaxis().tick_bottom()
divergence_axis.get_yaxis().tick_left()

pylab.figure(2,figsize=(3,2))
fig2 = pylab.gcf()
# make three panels panels
outer_grid2  = gridspec.GridSpec(1,1)

short_divergence_axis = plt.Subplot(fig2, outer_grid2[0])
fig2.add_subplot(short_divergence_axis)

short_divergence_axis.plot([0,1e-04],[0,1e-04],'k:')


all_syn_thinned_differences = []
all_syn_thinned_opportunities = []
all_non_differences = []
all_non_opportunities = []
all_pSs = []

median_pNs = []
median_pSs = []


# Plot percentiles of divergence distribution
for species_idx in xrange(0,len(species_names)):

    species_name = species_names[species_idx]

    thinned_syn_differences_1 = binomial(numpy.array(syn_differences[species_name],dtype=numpy.int32),0.5)
    thinned_syn_differences_2 = syn_differences[species_name]-thinned_syn_differences_1
    
    
    pS1s = thinned_syn_differences_1*1.0/syn_opportunities[species_name]*2
    pS2s = thinned_syn_differences_2*1.0/syn_opportunities[species_name]*2
    pSs = syn_differences[species_name]*1.0/syn_opportunities[species_name]
    pNs = non_differences[species_name]*1.0/non_opportunities[species_name]
    
    pseudo_pSs = syn_pseudocounts[species_name]*1.0/syn_opportunities[species_name]
    pseudo_pNs = non_pseudocounts[species_name]*1.0/non_opportunities[species_name]
    
    pS2s = numpy.clip(pS2s,1e-06,1)
    
    pNpSs = ((pseudo_pNs+pNs)/(pseudo_pSs+pS1s) )
    
    median_pSs.append( numpy.median(pSs) )
    median_pNs.append( numpy.median(pNs) )
    
    all_syn_thinned_differences.extend( thinned_syn_differences_1 )
    all_syn_thinned_opportunities.extend( syn_opportunities[species_name]/2.0 )
    all_non_differences.extend( non_differences[species_name] )
    all_non_opportunities.extend( non_opportunities[species_name] )
    all_pSs.extend(pS2s)
    
    if species_name.startswith('Bacteroides_vulgatus'):
        divergence_axis.loglog(pS2s, pNpSs, 'r.', markersize=2,markeredgewidth=0,zorder=1,label=("%s" % species_name))
    
    divergence_axis.loglog(pSs, pNpSs, '.', color='0.7', markersize=2,alpha=0.5,markeredgewidth=0,zorder=0)
    
    short_divergence_axis.plot(pNs, pSs, '.', color='b', markersize=3,alpha=0.5,markeredgewidth=0,zorder=0)
    
all_syn_thinned_differences = numpy.array(all_syn_thinned_differences)
all_syn_thinned_opportunities = numpy.array(all_syn_thinned_opportunities)
all_non_differences = numpy.array(all_non_differences)
all_non_opportunities = numpy.array(all_non_opportunities)
all_pSs = numpy.array(all_pSs)

    

# Calculate pvalues
pvalues = []
pS_thresholds = numpy.logspace(-6,-1,20)
for pSstar in pS_thresholds:

    good_idxs = (all_pSs<=pSstar)

    if not good_idxs.any():
        pvalues.append(1)
        continue

    nobs = all_non_differences[good_idxs].sum()
    lam = (all_non_opportunities[good_idxs].sum())*1.0/(all_syn_thinned_opportunities[good_idxs].sum())*(all_syn_thinned_differences[good_idxs].sum())
    
    pNpS = nobs/lam
    
    pvalue = 1-poisson_distribution.cdf(nobs-1, lam)
    
    print pNpS, nobs, lam, pvalue
    pvalues.append( pvalue )
    
      
median_pSs = numpy.array(median_pSs)
median_pNs = numpy.array(median_pNs)   

divergence_axis.plot([1e-09],[100], '.', color='0.7', markersize=2,alpha=0.5,markeredgewidth=0,zorder=0,label='All species')
       
divergence_axis.loglog(median_pSs, median_pNs*1.0/median_pSs, 'kx',markersize=2,label='Species median',alpha=0.5)

divergence_axis.legend(loc='upper right',frameon=False,numpoints=1)

divergence_axis.set_ylim([1e-02,10])    
divergence_axis.set_xlim([1e-06,1e-01])

theory_ds = numpy.logspace(-6,-1,100)

asymptotic_dNdS = 0.12
dStar = 3e-04
sbymu = 1/dStar/asymptotic_dNdS
print "s/u =", sbymu
print "s =", sbymu*1e-09
theory_dNdSs = asymptotic_dNdS+(1-asymptotic_dNdS)*(1-numpy.exp(-sbymu*theory_ds))/(theory_ds*sbymu)

divergence_axis.loglog(theory_ds, theory_dNdSs,'k-')

short_divergence_axis.set_xlim([-1e-06,3e-05])
short_divergence_axis.set_ylim([-1e-06,4e-05])

sys.stderr.write("Saving figure...\t")
fig.savefig('%s/figure_5.5.png' % (parse_midas_data.analysis_directory),bbox_inches='tight', dpi=600)
fig2.savefig('%s/figure_5.5.2.png' % (parse_midas_data.analysis_directory),bbox_inches='tight', dpi=600)

sys.stderr.write("Done!\n")

 