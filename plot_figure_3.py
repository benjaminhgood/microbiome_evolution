import matplotlib  
matplotlib.use('Agg') 
import sample_utils
import config
import parse_midas_data
import sample_utils
import os.path
import pylab
import sys
import numpy
from numpy.random import choice, binomial

import species_phylogeny_utils
import diversity_utils
import gene_diversity_utils
import calculate_substitution_rates
import calculate_singletons
import figure_utils

import stats_utils
import matplotlib.colors as colors
import matplotlib.cm as cmx
from math import log10,ceil
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy.random import randint
from mpl_toolkits.axes_grid.inset_locator import inset_axes

from numpy.random import randint, binomial, choice, poisson
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

num_bootstraps = 1000

min_coverage = config.min_median_coverage
alpha = 0.05 # Confidence interval range for rate estimates
low_divergence_threshold = 2e-04
min_change = 0.8
min_sample_size = 33 # 46 gives at least 1000 pairs
allowed_variant_types = set(['1D','2D','3D','4D'])

# Purifying selection model
# Fitted manually
asymptotic_dNdS = 0.12
dStar = 3e-04
sbymu = 1/dStar/asymptotic_dNdS
print "s/u =", sbymu
print "s =", sbymu*1e-09

def theory_dN(dS):
    return (asymptotic_dNdS+(1-asymptotic_dNdS)*(1-numpy.exp(-sbymu*dS))/(theory_ds*sbymu))*dS

####

divergence_matrices = {}
low_divergence_pair_counts = {}
null_low_divergence_pair_counts = [{} for i in xrange(0,num_bootstraps)]

low_divergence_gene_differences = []
low_divergence_clock_null_gene_differences = []
normal_divergence_gene_differences = []

good_species_list = parse_midas_data.parse_good_species_list()
if debug:
    good_species_list = ["Bacteroides_vulgatus_57955", "Bacteroides_uniformis_57318"]

sys.stderr.write("Loading sample metadata...\n")
subject_sample_map = sample_utils.parse_subject_sample_map()
sys.stderr.write("Done!\n")
       

####################################################
#
# Set up Figure (1 panels, arranged in 1x1 grid)
#
####################################################

haploid_color = '#08519c'

#pylab.figure(1,figsize=(7, 2))
#fig = pylab.gcf()
# make three panels panels
#outer_grid  = gridspec.GridSpec(1,2, width_ratios=[1,1],wspace=0.275)

# When it 1 then 2
#pylab.figure(1,figsize=(4.5, 2))
#fig = pylab.gcf()
#outer_grid  = gridspec.GridSpec(1,2, width_ratios=[3,1.3],wspace=0.4)
#right_grid = gridspec.GridSpecFromSubplotSpec(2,1, height_ratios=[1,1],hspace=0.5,subplot_spec=outer_grid[1])

pylab.figure(1,figsize=(5, 3))
fig = pylab.gcf()
outer_grid  = gridspec.GridSpec(1,1)

pylab.figure(2,figsize=(5, 3))
fig2 = pylab.gcf()
outer_grid2  = gridspec.GridSpec(1,1)


###########################################
#
# Do calculations and plotting for panel D (dN/dS vs dS)
#
###########################################


syn_differences = {}
syn_opportunities = {}
syn_pseudocounts = {}
non_differences = {}
non_pseudocounts = {}
non_opportunities = {}

core_differences = {}
core_opportunities = {}

data = {}

for species_name in good_species_list:

    #if not species_name.startswith('Bacteroides'):
    #    continue

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
    #sys.stderr.write("Loading pre-computed substitution rates for %s...\n" % species_name)
    #substitution_rate_map = calculate_substitution_rates.load_substitution_rate_map(species_name)
    #sys.stderr.write("Calculating matrices...\n")
    #dummy_samples, syn_difference_matrix, syn_opportunity_matrix = calculate_substitution_rates.calculate_matrices_from_substitution_rate_map(substitution_rate_map, '4D', allowed_samples=snp_samples)
    #dummy_samples, non_difference_matrix, non_opportunity_matrix = calculate_substitution_rates.calculate_matrices_from_substitution_rate_map(substitution_rate_map, '1D', allowed_samples=snp_samples)
    #snp_samples = dummy_samples
    
    # Load singleton matrices 
    sys.stderr.write("Loading pre-computed singleton rates for %s...\n" % species_name)
    singleton_rate_map = calculate_singletons.load_singleton_rate_map(species_name)
    sys.stderr.write("Calculating matrix...\n")
    
    
    snp_samples, syn_singleton_matrix, syn_doubleton_matrix, syn_difference_matrix, syn_opportunity_matrix  = calculate_singletons.calculate_matrices_from_singleton_rate_map(singleton_rate_map, '4D', allowed_samples=snp_samples)
    snp_samples, non_singleton_matrix, non_doubleton_matrix, non_difference_matrix, non_opportunity_matrix  = calculate_singletons.calculate_matrices_from_singleton_rate_map(singleton_rate_map, '1D', allowed_samples=snp_samples)
    snp_samples, core_singleton_matrix, core_doubleton_matrix, core_difference_matrix, core_opportunity_matrix  = calculate_singletons.calculate_matrices_from_singleton_rate_map(singleton_rate_map, '1D', allowed_samples=snp_samples)
    
    #snp_samples, singleton_matrix, doubleton_matrix, difference_matrix, opportunity_matrix = calculate_singletons.calculate_matrices_from_singleton_rate_map(singleton_rate_map, '4D', allowed_samples=snp_samples)
    
    singleton_matrix = syn_singleton_matrix + non_singleton_matrix
    doubleton_matrix = syn_doubleton_matrix + non_doubleton_matrix
    difference_matrix = syn_difference_matrix + non_difference_matrix
    opportunity_matrix = syn_opportunity_matrix + non_opportunity_matrix

    
    substitution_rate_matrix = difference_matrix*1.0/(opportunity_matrix+(opportunity_matrix==0))
    syn_substitution_rate_matrix = syn_difference_matrix*1.0/(syn_opportunity_matrix+(syn_opportunity_matrix==0))
    
    sys.stderr.write("Done!\n")
    
    
    syn_differences[species_name] = []
    syn_pseudocounts[species_name] = []
    syn_opportunities[species_name] = []
    
    non_differences[species_name] = []
    non_pseudocounts[species_name] = []
    non_opportunities[species_name] = []
    
    core_differences[species_name] = []
    core_opportunities[species_name] = []
    
    
    for i in xrange(0, syn_difference_matrix.shape[0]):
        for j in xrange(i+1, syn_difference_matrix.shape[0]):
            
            if syn_opportunity_matrix[i,j]>0 and non_opportunity_matrix[i,j]>0:
            
                syn_differences[species_name].append(syn_difference_matrix[i,j]+1)
                syn_pseudocounts[species_name].append(syn_opportunity_matrix[i,j]*1.0/(syn_opportunity_matrix[i,j]+non_opportunity_matrix[i,j]))
                syn_opportunities[species_name].append(syn_opportunity_matrix[i,j])
                
                non_differences[species_name].append( non_difference_matrix[i,j] )
                non_pseudocounts[species_name].append( non_opportunity_matrix[i,j]*1.0/(syn_opportunity_matrix[i,j]+non_opportunity_matrix[i,j]) )
                non_opportunities[species_name].append(non_opportunity_matrix[i,j])
                
                core_differences[species_name].append( core_difference_matrix[i,j] )
                core_opportunities[species_name].append( core_opportunity_matrix[i,j] )
        
    
    # Calculate singletons for closest sample
    
    closest_singleton_vector = []
    closest_difference_vector = []
    closest_opportunity_vector = []
    
    closest_syn_singleton_vector = []
    closest_syn_difference_vector = []
    closest_syn_opportunity_vector = []
    
    closest_non_singleton_vector = []
    closest_non_difference_vector = []
    closest_non_opportunity_vector = [] 
    
    for i in xrange(0, syn_substitution_rate_matrix.shape[0]):
    
        substitution_rates = syn_substitution_rate_matrix[i,:]
    
        # LEFT OFF HERE
    
        bad_idxs = (syn_opportunity_matrix[i,:]<0.5) # Make sure some things to compare to
        bad_idxs[i] = True # Don't compare it to itself
        #numpy.logical_or((difference_matrix[i,:]<0.5), (opportunity_matrix[i,:]<0.5))
        
        substitution_rates[bad_idxs] = 2 # so that it can't be too close
        
        closest_idx = substitution_rates.argmin()
        
        #print singleton_matrix[i,closest_idx], difference_matrix[i,closest_idx], opportunity_matrix[i,closest_idx]
        
        closest_singleton_vector.append( core_singleton_matrix[i,closest_idx] )
        closest_difference_vector.append( core_difference_matrix[i,closest_idx] )
        closest_opportunity_vector.append( core_opportunity_matrix[i,closest_idx]  )
        
        closest_syn_singleton_vector.append( syn_singleton_matrix[i,closest_idx] )
        closest_syn_difference_vector.append( syn_difference_matrix[i,closest_idx] )
        #closest_syn_difference_vector.append( syn_difference_matrix[i,closest_idx]-syn_singleton_matrix[i,closest_idx] + syn_singleton_matrix[closest_idx, i]  )
        closest_syn_opportunity_vector.append( syn_opportunity_matrix[i,closest_idx] )
        
        closest_non_singleton_vector.append( non_singleton_matrix[i,closest_idx] )
        closest_non_difference_vector.append( non_difference_matrix[i,closest_idx] )
        closest_non_opportunity_vector.append( non_opportunity_matrix[i,closest_idx] )
    
                
            
    syn_differences[species_name] = numpy.array(syn_differences[species_name])
    syn_pseudocounts[species_name] = numpy.array(syn_pseudocounts[species_name])
    
    syn_opportunities[species_name] = numpy.array(syn_opportunities[species_name])
    
    non_differences[species_name] = numpy.array(non_differences[species_name])
    non_pseudocounts[species_name] = numpy.array(non_pseudocounts[species_name])
    
    non_opportunities[species_name] = numpy.array(non_opportunities[species_name])    

    core_differences[species_name] = numpy.array(core_differences[species_name])
    core_opportunities[species_name] = numpy.array(core_opportunities[species_name])

    data[species_name] = numpy.array(closest_singleton_vector), numpy.array(closest_difference_vector), numpy.array(closest_opportunity_vector), numpy.array(closest_syn_singleton_vector), numpy.array(closest_syn_difference_vector), numpy.array(closest_syn_opportunity_vector), numpy.array(closest_non_singleton_vector), numpy.array(closest_non_difference_vector), numpy.array(closest_non_opportunity_vector)

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

divergence_axis = plt.Subplot(fig, outer_grid[0])
fig.add_subplot(divergence_axis)

divergence_axis.set_ylabel('Nonsynonymous ratio, $d_N/d_S$')
divergence_axis.set_xlabel('Synonymous divergence, $d_S$')

divergence_axis.spines['top'].set_visible(False)
divergence_axis.spines['right'].set_visible(False)
divergence_axis.get_xaxis().tick_bottom()
divergence_axis.get_yaxis().tick_left()

divergence_axis.set_ylabel('Nonsynonymous ratio, $d_N/d_S$')
divergence_axis.set_xlabel('Synonymous divergence, $d_S$')

divergence_axis.spines['top'].set_visible(False)
divergence_axis.spines['right'].set_visible(False)
divergence_axis.get_xaxis().tick_bottom()
divergence_axis.get_yaxis().tick_left()



cumulative_axis = inset_axes(divergence_axis, width="25%", height="25%", borderpad=0, bbox_to_anchor=(-0.01,0,1, 1), bbox_transform=divergence_axis.transAxes)

#cumulative_axis = plt.Subplot(fig, right_grid[0])
#fig.add_subplot(cumulative_axis)


cumulative_axis.spines['top'].set_visible(False)
cumulative_axis.spines['right'].set_visible(False)
cumulative_axis.get_xaxis().tick_bottom()
cumulative_axis.get_yaxis().tick_left()

cumulative_axis.set_ylabel('Cumulative $d_N/d_S$')
cumulative_axis.set_xlabel('Synonymous divergence, $d_S$')

line, = cumulative_axis.loglog([1e-05,1e-02],[1,1],'k:',linewidth=0.25,zorder=1)
line.set_dashes((1,1))

singleton_axis = plt.Subplot(fig2, outer_grid2[0])
fig2.add_subplot(singleton_axis)


singleton_axis.spines['top'].set_visible(False)
singleton_axis.spines['right'].set_visible(False)
singleton_axis.get_xaxis().tick_bottom()
singleton_axis.get_yaxis().tick_left()

singleton_axis.set_xlabel('Closest synonymous divergence, $d^*$')
singleton_axis.set_ylabel('Private nonsynonymous ratio, $d_N/d_S$')

cumulative_singleton_axis = inset_axes(singleton_axis, width="25%", height="25%", borderpad=0, bbox_to_anchor=(-0.01,0,1, 1), bbox_transform=singleton_axis.transAxes)
#-0.025

#cumulative_singleton_axis = plt.Subplot(fig, right_grid[1])
#fig.add_subplot(cumulative_singleton_axis)

line, = cumulative_singleton_axis.loglog([1e-05,1e-02],[1,1],'k:',linewidth=0.25,zorder=1)
line.set_dashes((1,1))

cumulative_singleton_axis.spines['top'].set_visible(False)
cumulative_singleton_axis.spines['right'].set_visible(False)
cumulative_singleton_axis.get_xaxis().tick_bottom()
cumulative_singleton_axis.get_yaxis().tick_left()

cumulative_singleton_axis.set_ylabel('Cumulative private $d_N/d_S$')
cumulative_singleton_axis.set_xlabel('Closest synonymous \n divergence, $d_S^*$')


all_syn_differences = []
all_syn_opportunities = []
all_non_differences = []
all_non_opportunities = []
all_core_differences = []
all_core_opportunities = []
median_pNs = []
median_pSs = []


# Plot percentiles of divergence distribution
for species_idx in xrange(0,len(species_names)):

    species_name = species_names[species_idx]

    # Use the Poisson thinning theorem to cut down on 
    # non-biological correlations between dS and dN/dS 
    # (i.e., fact that dS is in denominator of dN/dS.
    thinned_syn_differences_1 = binomial(numpy.array(syn_differences[species_name],dtype=numpy.int32),0.5)
    thinned_syn_differences_2 = syn_differences[species_name]-thinned_syn_differences_1
    
    
    pS1s = thinned_syn_differences_1*1.0/(syn_opportunities[species_name]/2.0)
    pS2s = thinned_syn_differences_2*1.0/(syn_opportunities[species_name]/2.0)
    pSs = syn_differences[species_name]*1.0/syn_opportunities[species_name]
    pNs = non_differences[species_name]*1.0/non_opportunities[species_name]
    ptots = (syn_differences[species_name]+non_differences[species_name])*1.0/(syn_opportunities[species_name] + non_opportunities[species_name])
    
    pseudo_pSs = 1.0/(syn_opportunities[species_name]/2.0+non_opportunities[species_name])
    pseudo_pNs = 1.0/(syn_opportunities[species_name]/2.0+non_opportunities[species_name])
    
    pS2s = numpy.clip(pS2s,1e-06,1)
    
    pNpSs = ((pseudo_pNs+pNs)/(pseudo_pSs+pS1s) )
    
    median_pSs.append( numpy.median(pSs) )
    median_pNs.append( numpy.median(pNs) )
    
    all_syn_differences.extend( syn_differences[species_name] )
    all_syn_opportunities.extend( syn_opportunities[species_name] )
        
    all_non_differences.extend( non_differences[species_name] )
    all_non_opportunities.extend( non_opportunities[species_name] )
    all_core_differences.extend( core_differences[species_name] )
    all_core_opportunities.extend( core_opportunities[species_name] )
    
    
    good_idxs = ((syn_differences[species_name]+non_differences[species_name])>=10)
    bad_idxs = numpy.logical_not(good_idxs)
    
    
    if False and species_name.startswith('Bacteroides'):
        divergence_axis.loglog(pS2s[good_idxs], pNpSs[good_idxs], 'r.', markersize=2,markeredgewidth=0,zorder=1,rasterized=True)
    else:
        divergence_axis.loglog(pS2s[good_idxs], pNpSs[good_idxs], '.', color='0.7', markersize=2,alpha=0.5,markeredgewidth=0,zorder=0,rasterized=True)
        #divergence_axis.loglog(ptots[good_idxs], pNpSs[good_idxs], '.', color='0.7', markersize=2,alpha=0.5,markeredgewidth=0,zorder=0,rasterized=True)
        
all_syn_differences = numpy.array(all_syn_differences,dtype=numpy.int32)
all_syn_opportunities = numpy.array(all_syn_opportunities,dtype=numpy.int32)
all_non_differences = numpy.array(all_non_differences,dtype=numpy.int32)
all_non_opportunities = numpy.array(all_non_opportunities,dtype=numpy.int32)
all_core_differences = numpy.array(all_non_differences,dtype=numpy.int32)
all_core_opportunities = numpy.array(all_non_opportunities,dtype=numpy.int32)

all_core_divergence = all_core_differences*1.0/all_core_opportunities
all_syn_divergence = all_syn_differences*1.0/all_syn_opportunities

all_NS_differences = all_syn_differences + all_non_differences
all_NS_opportunities = all_syn_opportunities + all_non_opportunities
all_fractions = all_non_differences*1.0/(all_NS_differences+(all_NS_differences==0))
    
ds = numpy.logspace(-5,-2,20)

cf_ratios = [] # cumulative estimates <= total d
sf_ratios = [] # cumulative estimates >= total d

sys.stderr.write("Bootstrapping dN/dS...\n")
num_bootstraps = 500
for bootstrap_idx in xrange(0,num_bootstraps):
    
    lower_pNpSs = []
    upper_pNpSs = []
    
    # bootstrap dataset using poisson resampling
    # Pseudocounts are added so that things w/ zero counts are not "stuck" there in resampling
    # Pseudocounts are chosen w/ dN/dS=1, so should be conservative?
    # (alternatively, we could choose dN/dS=0.1, but that seems a little unfair)
    pseudocount = 0 #1.0
    bootstrapped_non_differences = poisson(all_non_differences+pseudocount) 
    bootstrapped_syn_differences = poisson(all_syn_differences+all_syn_opportunities*pseudocount/all_non_opportunities)
    bootstrapped_NS_differences = bootstrapped_non_differences + bootstrapped_syn_differences
    bootstrapped_thinned_syn_differences_1 = binomial(bootstrapped_syn_differences,0.5)
    bootstrapped_thinned_syn_differences_2 = bootstrapped_syn_differences-bootstrapped_thinned_syn_differences_1
    
    bootstrapped_divergence = bootstrapped_thinned_syn_differences_1 / (all_syn_opportunities/2.0)
    
    for d in ds:
        
        lower_idxs = (bootstrapped_divergence <= d)*(all_NS_differences>0.5)*(bootstrapped_NS_differences>0.5)
        upper_idxs = (bootstrapped_divergence > d)*(all_NS_differences>0.5)*(bootstrapped_NS_differences>0.5)
        
        if lower_idxs.sum()<1.5:
            lower_pNpSs.append(-1)
        else:
            
            lower_cumulative_non_differences = (bootstrapped_non_differences)[lower_idxs].sum()
            lower_cumulative_expected_non_differences = (bootstrapped_thinned_syn_differences_2[lower_idxs]*2.0/all_syn_opportunities[lower_idxs]*all_non_opportunities[lower_idxs]).sum() 
            lower_pNpSs.append( (lower_cumulative_non_differences)/(lower_cumulative_expected_non_differences) )
        
        
        if upper_idxs.sum()<1.5:
            upper_pNpSs.append(-1)
        else:
            upper_cumulative_non_differences = (bootstrapped_non_differences[upper_idxs]).sum()
            upper_cumulative_expected_non_differences = (bootstrapped_thinned_syn_differences_2[upper_idxs]*2.0/all_syn_opportunities[upper_idxs]*all_non_opportunities[upper_idxs]).sum() 
            upper_pNpSs.append( (upper_cumulative_non_differences)/(upper_cumulative_expected_non_differences) )
        
        
    cf_ratios.append(lower_pNpSs)
    sf_ratios.append(upper_pNpSs)
    
    
cf_ratios = numpy.array(cf_ratios)
sf_ratios = numpy.array(sf_ratios)

avg_cf_ratios = []
std_cf_ratios = []
median_cf_ratios = []
lower_cf_ratios = []
upper_cf_ratios = []

avg_sf_ratios = [] 
std_sf_ratios = [] 

for i in xrange(0,len(ds)):
    
    ratios = numpy.sort(cf_ratios[:,i])
    good_idxs = (ratios>-0.5)
    if good_idxs.sum()<1.5:
        avg_cf_ratios.append(-1)
        std_cf_ratios.append(0)
    else:
    
        median_cf_ratios.append(numpy.median(ratios[good_idxs]))
        idx = long(0.025*good_idxs.sum())
        lower_cf_ratios.append( ratios[good_idxs][idx] )
        upper_cf_ratios.append(ratios[good_idxs][-idx-1])
    
        avg_cf_ratios.append( ratios[good_idxs].mean() )
        std_cf_ratios.append( ratios[good_idxs].std() )
    
    ratios = sf_ratios[:,i]
    good_idxs = (ratios>-0.5)
    if good_idxs.sum()<1.5:
        avg_sf_ratios.append(-1)
        std_sf_ratios.append(0)
    else:
        avg_sf_ratios.append( ratios[good_idxs].mean() )
        std_sf_ratios.append( ratios[good_idxs].std() )        
    
avg_cf_ratios = numpy.array(avg_cf_ratios)
std_cf_ratios = numpy.array(std_cf_ratios)
median_cf_ratios = numpy.array(median_cf_ratios)
upper_cf_ratios = numpy.array(upper_cf_ratios)
lower_cf_ratios = numpy.array(lower_cf_ratios)

avg_sf_ratios = numpy.array(avg_sf_ratios)
std_sf_ratios = numpy.array(std_sf_ratios)

#good_idxs = (avg_sf_ratios>-0.5)
#cumulative_axis.fill_between(ds[good_idxs], avg_sf_ratios[good_idxs]-2*std_sf_ratios[good_idxs], avg_sf_ratios[good_idxs]+2*std_sf_ratios[good_idxs],color='b',linewidth=0)
#cumulative_axis.loglog(ds[good_idxs], avg_sf_ratios[good_idxs],'b-')

good_idxs = (avg_cf_ratios>-0.5)
cumulative_axis.fill_between(ds[good_idxs], lower_cf_ratios[good_idxs], upper_cf_ratios[good_idxs],color='0.7',linewidth=0,zorder=0)
cumulative_axis.loglog(ds[good_idxs], avg_cf_ratios[good_idxs],'k-',zorder=2)
      
median_pSs = numpy.array(median_pSs)
median_pNs = numpy.array(median_pNs)   

divergence_axis.plot([1e-09],[100], 'o', color='0.7', markersize=2,markeredgewidth=0,zorder=0,label='(Species x host x host)')
       
divergence_axis.loglog(median_pSs, median_pNs*1.0/median_pSs, 'kx',markersize=2,label='Median of each species',alpha=0.5)


divergence_axis.set_ylim([1e-02,10])    
divergence_axis.set_xlim([1e-06,1e-01])

theory_ds = numpy.logspace(-6,-1,100)

#theory_dNdSs = asymptotic_dNdS+(1-asymptotic_dNdS)*(1-numpy.exp(-sbymu*theory_ds))/(theory_ds*sbymu)
theory_dNdSs = theory_dN(theory_ds)/theory_ds

line, = divergence_axis.loglog([1e-06,1e-03],[1,1],'k:',linewidth=0.25,label='Neutral model')
line.set_dashes((1,1))
divergence_axis.loglog(theory_ds, theory_dNdSs,'r-',label='Purifying selection model')

divergence_axis.legend(loc='lower left',frameon=False,numpoints=1)
#divergence_axis.legend(loc='upper right',frameon=False,numpoints=1)


cumulative_axis.set_xlim([1e-05,1e-02])
cumulative_axis.set_ylim([5e-02,2])
cumulative_singleton_axis.set_ylim([5e-02,2])


singleton_axis.set_xlim([1e-06,1e-02])
singleton_axis.set_ylim([1e-01,1e01])
singleton_axis.loglog([1e-06,1e-01],[1,1],'k:')
#cumulative_axis.set_xticklabels([])

all_closest_differences = []
all_closest_opportunities = []

all_syn_singletons = []
all_syn_differences = []
all_syn_opportunities = []
all_non_singletons = []
all_non_opportunities = []


for species_idx in xrange(0,len(species_names)):

    species_name = species_names[species_idx]

    closest_singleton_vector, closest_difference_vector, closest_opportunity_vector, closest_syn_singleton_vector, closest_syn_difference_vector, closest_syn_opportunity_vector, closest_non_singleton_vector, closest_non_difference_vector, closest_non_opportunity_vector = data[species_name]
    
    all_closest_differences.extend(closest_difference_vector)
    all_closest_opportunities.extend(closest_opportunity_vector)
    
    all_syn_singletons.extend(closest_syn_singleton_vector)
    all_syn_differences.extend(closest_syn_difference_vector)
    all_syn_opportunities.extend(closest_syn_opportunity_vector)
    
    if not (closest_syn_opportunity_vector>=0).all():
        print species_name, "bad closest syn opportunities"
    
    all_non_singletons.extend(closest_non_singleton_vector)
    all_non_opportunities.extend(closest_non_opportunity_vector)
    
    syn_pseudocount_vector = closest_syn_opportunity_vector*1.0/(closest_syn_opportunity_vector + closest_non_opportunity_vector)
    non_pseudocount_vector = closest_non_opportunity_vector*1.0/(closest_syn_opportunity_vector + closest_non_opportunity_vector)
    
     
    closest_divergence_vector = closest_difference_vector*1.0/closest_opportunity_vector
    singleton_pNpS_vector = (closest_non_singleton_vector+non_pseudocount_vector)/(closest_non_opportunity_vector)*1.0/((closest_syn_singleton_vector+syn_pseudocount_vector)*1.0/(closest_syn_opportunity_vector))
    
    singleton_axis.semilogx(closest_divergence_vector, singleton_pNpS_vector, '.', color='0.7', markersize=2,alpha=0.5,markeredgewidth=0,zorder=0,rasterized=True)


all_closest_differences = numpy.array(all_closest_differences,dtype=numpy.int32)
all_closest_opportunities = numpy.array(all_closest_opportunities,dtype=numpy.int32)
all_syn_singletons = numpy.array(all_syn_singletons,dtype=numpy.int32)
all_syn_differences = numpy.array(all_syn_differences,dtype=numpy.int32)
all_syn_opportunities = numpy.array(all_syn_opportunities,dtype=numpy.int32)
all_non_singletons = numpy.array(all_non_singletons,dtype=numpy.int32)
all_non_opportunities = numpy.array(all_non_opportunities,dtype=numpy.int32)

avg_singleton_dnds = all_non_singletons.sum()*1.0/all_non_opportunities.sum()/(all_syn_singletons.sum()*1.0/all_syn_opportunities.sum())

all_core_divergence = all_closest_differences*1.0/all_closest_opportunities

all_NS_singletons = all_syn_singletons + all_non_singletons
all_NS_opportunities = all_syn_opportunities + all_non_opportunities
all_fractions = all_non_singletons*1.0/(all_NS_singletons+(all_NS_singletons==0))

####
## Need to fix this!
##
all_syn_other = all_syn_differences-all_syn_singletons
#all_syn_other = all_syn_differences # now that we've removed singletons from here

    
print (all_syn_opportunities>=0).all()
    
ds = numpy.logspace(-5,-2,20)

cf_ratios = [] # cumulative estimates <= total d
sf_ratios = [] # cumulative estimates >= total d

sys.stderr.write("Bootstrapping singleton dN/dS...\n")
num_bootstraps = 10000
for bootstrap_idx in xrange(0,num_bootstraps):
    
    # bootstrap dataset using poisson resampling
    # Pseudocounts are added so that things w/ zero counts are not "stuck" there in resampling
    # Pseudocounts are chosen w/ dN/dS=1, so should be conservative?
    # (alternatively, we could choose dN/dS=0.1, but that seems a little unfair)
    pseudocount = 0.0
    bootstrapped_non_singletons = poisson(all_non_singletons+pseudocount) 
    bootstrapped_syn_singletons = poisson(all_syn_singletons+all_syn_opportunities*pseudocount/all_non_opportunities)
    bootstrapped_syn_other = poisson(all_syn_other)

    bootstrapped_NS_singletons = bootstrapped_non_singletons + bootstrapped_syn_singletons
    
    bootstrapped_thinned_syn_singletons_1 = binomial(bootstrapped_syn_singletons,0.5)
    bootstrapped_thinned_syn_singletons_2 = bootstrapped_syn_singletons-bootstrapped_thinned_syn_singletons_1
    bootstrapped_thinned_syn_other_1 = binomial(bootstrapped_syn_other,0.5)
    bootstrapped_divergence = (bootstrapped_thinned_syn_singletons_1 + bootstrapped_thinned_syn_other_1) / (all_syn_opportunities/2.0)
    
    # Now that we've removed singletons from divergence estimate
    #bootstrapped_thinned_syn_singletons_2 = bootstrapped_syn_singletons/2.0
    #bootstrapped_divergence = bootstrapped_syn_other * 1.0 / (all_syn_opportunities)
    
    #print (bootstrapped_non_singletons>=0).all()
    #print (bootstrapped_thinned_syn_singletons_2>=0).all()
    
    
    lower_pNpSs = []
    upper_pNpSs = []
    
    for d in ds:
        
        lower_idxs = (bootstrapped_divergence <= d)*(all_NS_singletons>0.5)*(bootstrapped_NS_singletons>0.5)
        upper_idxs = (bootstrapped_divergence > d)*(all_NS_singletons>0.5)*(bootstrapped_NS_singletons>0.5)
        
        if lower_idxs.sum()<1.5:
            lower_pNpSs.append(-1)
        else:
        
            numerators = (bootstrapped_non_singletons[lower_idxs])
            denominators = (bootstrapped_thinned_syn_singletons_2[lower_idxs]*2.0/all_syn_opportunities[lower_idxs]*all_non_opportunities[lower_idxs]) 
        
            lower_cumulative_non_singletons = numerators.sum() 
            lower_cumulative_expected_non_singletons = denominators.sum()
            lower_pNpSs.append( (lower_cumulative_non_singletons)/(lower_cumulative_expected_non_singletons) )
        
        
        
        if upper_idxs.sum()<1.5:
            upper_pNpSs.append(-1)
        else:
        
            
            upper_cumulative_non_singletons = (bootstrapped_non_singletons[upper_idxs]).sum()
            upper_cumulative_expected_non_singletons = (bootstrapped_thinned_syn_singletons_2[upper_idxs]*2.0/all_syn_opportunities[upper_idxs]*all_non_opportunities[upper_idxs]).sum() 
            upper_pNpSs.append( (upper_cumulative_non_singletons)/(upper_cumulative_expected_non_singletons) )
        

        
        
    cf_ratios.append(lower_pNpSs)
    sf_ratios.append(upper_pNpSs)
    
    
cf_ratios = numpy.array(cf_ratios)
sf_ratios = numpy.array(sf_ratios)

sys.stderr.write("Done!\n")

avg_cf_ratios = []
std_cf_ratios = []
median_cf_ratios = []
lower_cf_ratios = []
upper_cf_ratios = []

avg_sf_ratios = [] 
std_sf_ratios = [] 

for i in xrange(0,len(ds)):
    
    ratios = numpy.sort(cf_ratios[:,i])
    good_idxs = (ratios>-0.5)
    if good_idxs.sum()<1.5:
        avg_cf_ratios.append(-1)
        std_cf_ratios.append(0)
    else:
    
        median_cf_ratios.append(numpy.median(ratios[good_idxs]))
        idx = long(0.025*good_idxs.sum())
        lower_cf_ratios.append( ratios[good_idxs][idx] )
        upper_cf_ratios.append(ratios[good_idxs][-idx-1])
    
        avg_cf_ratios.append( ratios[good_idxs].mean() )
        std_cf_ratios.append( ratios[good_idxs].std() )
    
    ratios = sf_ratios[:,i]
    good_idxs = (ratios>-0.5)
    if good_idxs.sum()<1.5:
        avg_sf_ratios.append(-1)
        std_sf_ratios.append(0)
    else:
        avg_sf_ratios.append( ratios[good_idxs].mean() )
        std_sf_ratios.append( ratios[good_idxs].std() )        
    
avg_cf_ratios = numpy.array(avg_cf_ratios)
std_cf_ratios = numpy.array(std_cf_ratios)
median_cf_ratios = numpy.array(median_cf_ratios)
upper_cf_ratios = numpy.array(upper_cf_ratios)
lower_cf_ratios = numpy.array(lower_cf_ratios)

avg_sf_ratios = numpy.array(avg_sf_ratios)
std_sf_ratios = numpy.array(std_sf_ratios)

#good_idxs = (avg_sf_ratios>-0.5)
#cumulative_axis.fill_between(ds[good_idxs], avg_sf_ratios[good_idxs]-2*std_sf_ratios[good_idxs], avg_sf_ratios[good_idxs]+2*std_sf_ratios[good_idxs],color='b',linewidth=0)
#cumulative_axis.loglog(ds[good_idxs], avg_sf_ratios[good_idxs],'b-')

good_idxs = (avg_cf_ratios>-0.5)
cumulative_singleton_axis.fill_between(ds[good_idxs], lower_cf_ratios[good_idxs], upper_cf_ratios[good_idxs],color='0.7',linewidth=0,zorder=0)
cumulative_singleton_axis.loglog(ds[good_idxs], avg_cf_ratios[good_idxs],'k-',zorder=1)



sys.stderr.write("Saving figure...\t")
fig.savefig('%s/figure_3.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight',dpi=600)
fig2.savefig('%s/supplemental_singleton_dNdS.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight',dpi=600)
sys.stderr.write("Done!\n")

 