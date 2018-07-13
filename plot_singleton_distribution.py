import matplotlib  
matplotlib.use('Agg') 
import sample_utils
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
import calculate_singletons

import stats_utils
import matplotlib.colors as colors
import matplotlib.cm as cmx
from math import log10,ceil,log
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy.random import randint, shuffle, poisson

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

data = {}

good_species_list = parse_midas_data.parse_good_species_list()
if debug:
    good_species_list = ["Bacteroides_vulgatus_57955", "Bacteroides_uniformis_57318"]

sys.stderr.write("Loading sample metadata...\n")
subject_sample_map = sample_utils.parse_subject_sample_map()
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

    # Load singleton matrices 
    sys.stderr.write("Loading pre-computed singleton rates for %s...\n" % species_name)
    singleton_rate_map = calculate_singletons.load_singleton_rate_map(species_name)
    sys.stderr.write("Calculating matrix...\n")
    snp_samples, singleton_matrix, doubleton_matrix, difference_matrix, opportunity_matrix = calculate_singletons.calculate_matrices_from_singleton_rate_map(singleton_rate_map, '4D', allowed_samples=snp_samples)
    
    snp_samples, syn_singleton_matrix, syn_doubleton_matrix, syn_difference_matrix, syn_opportunity_matrix  = calculate_singletons.calculate_matrices_from_singleton_rate_map(singleton_rate_map, '4D', allowed_samples=snp_samples)
    snp_samples, non_singleton_matrix, non_doubleton_matrix, non_difference_matrix, non_opportunity_matrix  = calculate_singletons.calculate_matrices_from_singleton_rate_map(singleton_rate_map, '1D', allowed_samples=snp_samples)
    
    substitution_rate_matrix = difference_matrix*1.0/(opportunity_matrix+(opportunity_matrix==0))
    
    sys.stderr.write("Done!\n")
    
    # Calculate relative singleton matrices (# singletons / # differences to closest sample)
    
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
    
    for i in xrange(0, substitution_rate_matrix.shape[0]):
    
        substitution_rates = substitution_rate_matrix[i,:]
    
        bad_idxs = numpy.logical_or((difference_matrix[i,:]<0.5), (opportunity_matrix[i,:]<0.5))
        
        substitution_rates[bad_idxs] = 2 # so that it can't be too close
        
        closest_idx = substitution_rates.argmin()
        
        #print singleton_matrix[i,closest_idx], difference_matrix[i,closest_idx], opportunity_matrix[i,closest_idx]
        
        closest_singleton_vector.append( singleton_matrix[i,closest_idx] )
        closest_difference_vector.append( difference_matrix[i,closest_idx] )
        closest_opportunity_vector.append(opportunity_matrix[i,closest_idx]  )
        
        closest_syn_singleton_vector.append( syn_singleton_matrix[i,closest_idx] )
        closest_syn_difference_vector.append( syn_difference_matrix[i,closest_idx] )
        closest_syn_opportunity_vector.append( syn_opportunity_matrix[i,closest_idx] )
        
        closest_non_singleton_vector.append( non_singleton_matrix[i,closest_idx] )
        closest_non_difference_vector.append( non_difference_matrix[i,closest_idx] )
        closest_non_opportunity_vector.append( non_opportunity_matrix[i,closest_idx] )
        

    #print closest_divergence_vector, relative_singleton_rate_vector
    data[species_name] = numpy.array(closest_singleton_vector), numpy.array(closest_difference_vector), numpy.array(closest_opportunity_vector), numpy.array(closest_syn_singleton_vector), numpy.array(closest_syn_difference_vector), numpy.array(closest_syn_opportunity_vector), numpy.array(closest_non_singleton_vector), numpy.array(closest_non_difference_vector), numpy.array(closest_non_opportunity_vector)
    
    #print numpy.array(closest_opportunity_vector)
   
# # low divergence strains across species
# # samples
        
species_names = []
sample_sizes = []

for species_name in species_phylogeny_utils.sort_phylogenetically(data.keys()):
    species_names.append(species_name)
    sample_sizes.append( len(data[species_name][0]) )
    
# sort in descending order of sample size
# Sort by num haploids    
#sample_sizes, species_names = zip(*sorted(zip(sample_sizes, species_names),reverse=True))
    
sys.stderr.write("Postprocessing %d species...\n" % len(species_names))

divergences = numpy.logspace(-4,-1,10)
        

####################################################
#
# Set up Figure (1 panels, arranged in 1x1 grid)
#
####################################################

haploid_color = '#08519c'

pylab.figure(1,figsize=(3.42,2))
fig = pylab.gcf()
# make three panels panels
outer_grid = gridspec.GridSpec(1,1)
#outer_grid  = gridspec.GridSpec(2,1, height_ratios=[1,1],hspace=0.2)

singleton_axis = plt.Subplot(fig, outer_grid[0])
fig.add_subplot(singleton_axis)

singleton_axis.set_ylabel('Fraction singletons')
singleton_axis.set_xlabel('Minimum synonymous divergence, $d_S$') 
singleton_axis.set_xlim([1e-06,1e-01])
#singleton_axis.set_ylim([0,1])

singleton_axis.spines['top'].set_visible(False)
singleton_axis.spines['right'].set_visible(False)
singleton_axis.get_xaxis().tick_bottom()
singleton_axis.get_yaxis().tick_left()

pylab.figure(2,figsize=(3.42,2))
fig2 = pylab.gcf()
# make three panels panels
outer_grid = gridspec.GridSpec(1,1)

dnds_axis = plt.Subplot(fig2, outer_grid[0])
fig2.add_subplot(dnds_axis)

dnds_axis.set_xlabel('Minimum synonymous divergence, $d_S$')
dnds_axis.set_ylabel('Singleton $dN/dS$')
dnds_axis.set_xlim([1e-06,1e-01])
#dnds_axis.set_ylim([0,1])

dnds_axis.spines['top'].set_visible(False)
dnds_axis.spines['right'].set_visible(False)
dnds_axis.get_xaxis().tick_bottom()
dnds_axis.get_yaxis().tick_left()

total_closest_differences = []
total_closest_opportunities = []
total_closest_singletons = []
total_closest_singleton_opportunities = []
ds = numpy.logspace(-5,-2,50)

cumulative_singleton_differences = numpy.zeros_like(ds)
cumulative_singleton_opportunities = numpy.zeros_like(ds)

cumulative_syn_singleton_differences = numpy.zeros_like(ds)
cumulative_syn_singleton_opportunities = numpy.zeros_like(ds)
cumulative_syn_opportunities = numpy.zeros_like(ds)

cumulative_non_singleton_differences = numpy.zeros_like(ds)
cumulative_non_singleton_opportunities = numpy.zeros_like(ds)
cumulative_non_opportunities = numpy.zeros_like(ds)

cumulative_closest_differences = numpy.zeros_like(ds)
cumulative_closest_opportunities = numpy.zeros_like(ds)

for species_idx in xrange(0,len(species_names)):

    species_name = species_names[species_idx]

    closest_singleton_vector, closest_difference_vector, closest_opportunity_vector, closest_syn_singleton_vector, closest_syn_difference_vector, closest_syn_opportunity_vector, closest_non_singleton_vector, closest_non_difference_vector, closest_non_opportunity_vector = data[species_name]
    
    total_closest_differences.extend(closest_difference_vector)
    total_closest_opportunities.extend(closest_opportunity_vector)
    total_closest_singletons.extend(closest_singleton_vector)
    total_closest_singleton_opportunities.extend(closest_difference_vector)
     
    closest_divergence_vector = closest_difference_vector*1.0/closest_opportunity_vector
     
    # Add to cumulative distribution
    for d_idx in xrange(0,len(ds)):
        d = ds[d_idx]
        
        lower_idxs = (closest_divergence_vector<=d)
        
        cumulative_singleton_differences[d_idx] += closest_singleton_vector[lower_idxs].sum()
        cumulative_singleton_opportunities[d_idx] += closest_difference_vector[lower_idxs].sum()
        
        cumulative_syn_singleton_differences[d_idx] += closest_syn_singleton_vector[lower_idxs].sum()
        cumulative_syn_singleton_opportunities[d_idx] += closest_syn_difference_vector[lower_idxs].sum()
        cumulative_syn_opportunities[d_idx] += closest_syn_opportunity_vector[lower_idxs].sum()
        
        cumulative_non_singleton_differences[d_idx] += closest_non_singleton_vector[lower_idxs].sum()
        cumulative_non_singleton_opportunities[d_idx] += closest_non_difference_vector[lower_idxs].sum()
        cumulative_non_opportunities[d_idx] += closest_non_opportunity_vector[lower_idxs].sum()
        
        cumulative_closest_differences[d_idx] += closest_difference_vector[lower_idxs].sum()
        cumulative_closest_opportunities[d_idx] += closest_opportunity_vector[lower_idxs].sum() 
    

print ds
print cumulative_non_singleton_differences
    
cumulative_singleton_rates = cumulative_singleton_differences*1.0/cumulative_singleton_opportunities

cumulative_closest_rates = cumulative_closest_differences*1.0/cumulative_closest_opportunities 

cumulative_non_singleton_rates = cumulative_non_singleton_differences*1.0/cumulative_non_opportunities

cumulative_syn_singleton_rates = cumulative_syn_singleton_differences*1.0/cumulative_syn_opportunities

cumulative_relative_singleton_rates = cumulative_singleton_rates

cumulative_pNpS =  (cumulative_non_singleton_differences/cumulative_non_opportunities)/(cumulative_syn_singleton_differences/cumulative_syn_opportunities)

bootstrapped_relative_singleton_ratess = []
bootstrapped_pNpSs = []
num_bootstraps = 1000
for bootstrap_idx in xrange(0,num_bootstraps):

    # resampe everything at known rates
    bootstrapped_cumulative_singleton_rates = poisson(cumulative_singleton_rates*cumulative_singleton_opportunities)*1.0/cumulative_singleton_opportunities
    
    bootstrapped_cumulative_closest_rates = poisson(cumulative_closest_rates*cumulative_closest_opportunities)*1.0/ cumulative_closest_opportunities
    
    bootstrapped_relative_singleton_ratess.append( bootstrapped_cumulative_singleton_rates ) 
    
    bootstrapped_syn_singleton_rates = poisson(cumulative_syn_singleton_rates*cumulative_syn_opportunities)*1.0/cumulative_syn_opportunities

    bootstrapped_non_singleton_rates = poisson(cumulative_non_singleton_rates*cumulative_non_opportunities)*1.0/cumulative_non_opportunities

    bootstrapped_pNpSs.append( bootstrapped_non_singleton_rates/bootstrapped_syn_singleton_rates )

bootstrapped_relative_singleton_ratess = numpy.array(bootstrapped_relative_singleton_ratess)
    
bootstrapped_pNpSs = numpy.array(bootstrapped_pNpSs) 

avg_rates = bootstrapped_relative_singleton_ratess.mean(axis=0)
std_rates = bootstrapped_relative_singleton_ratess.std(axis=0)

avg_pNpSs = bootstrapped_pNpSs.mean(axis=0)
std_pNpSs = bootstrapped_pNpSs.std(axis=0)


singleton_axis.fill_between(ds, avg_rates-2*std_rates, avg_rates+2*std_rates,color='0.7',linewidth=0)

singleton_axis.loglog(ds, avg_rates, 'k-')

dnds_axis.fill_between(ds, avg_pNpSs-2*std_pNpSs, avg_pNpSs+2*std_pNpSs,color='0.7',linewidth=0)

dnds_axis.loglog(ds, avg_pNpSs, 'k-')

sys.stderr.write("Saving figure...\t")
fig.savefig('%s/supplemental_singleton_distribution.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight',dpi=600)
sys.stderr.write("Done!\n")
sys.stderr.write("Saving figure...\t")
fig2.savefig('%s/supplemental_singleton_pnps.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight',dpi=600)
sys.stderr.write("Done!\n")
sys.exit(0)

## Do bootstrapping and significance testing for singleton ratios
total_closest_differences = numpy.array(total_closest_differences)
total_closest_opportunities = numpy.array(total_closest_opportunities)

total_closest_singleton_differences = numpy.array(total_closest_singleton_differences)
total_closest_singleton_opportunities = numpy.array(total_closest_singleton_opportunities)

total_closest_rates = total_closest_differences*1.0/total_closest_opportunities

pavg = total_closest_singleton_differences.sum()*1.0/total_closest_singleton_opportunities.sum()

dstars = numpy.array([2e-04])

plowers = []
puppers = []
LRTs = []
for dstar in dstars:
    
    lower_idxs = (total_closest_rates<=dstar)
    upper_idxs = (total_closest_rates>dstar)

    print total_singleton_differences[lower_idxs].sum(), total_closest_differences[lower_idxs].sum()

    plower = total_closest_singleton_differences[lower_idxs].sum()*1.0/total_closest_singleton_opportunities[lower_idxs].sum()
    
    pupper = total_closest_singleton_differences[upper_idxs].sum()*1.0/total_closest_singleton_opportunities[upper_idxs].sum()

    LRT = total_closest_singleton_differences[upper_idxs].sum()*log(pupper/pavg) + (pavg-pupper)*(total_closest_differences[upper_idxs].sum()) + total_closest_singleton_differences[lower_idxs].sum()*log(plower/pavg) + (pavg-plower)*(total_closest_differences[lower_idxs].sum())

    plowers.append( plower )
    
    puppers.append( pupper )
    
    LRTs.append( LRT ) 

plowers = numpy.array(plowers)
puppers = numpy.array(puppers)
LRTs = numpy.array(LRTs)
   
print plowers
    

bootstrapped_plowerss = []
bootstrapped_pupperss = []
bootstrapped_pavgss = []
bootstrapped_LRTss = []

num_bootstraps = 1000

for bootstrap_idx in xrange(0,num_bootstraps):

    bootstrapped_closest_rates = numpy.array(total_closest_rates,copy=True)
    shuffle(bootstrapped_closest_rates)
    
    bootstrapped_closest_differences = poisson(bootstrapped_closest_rates*total_closest_opportunities)
    bootstrapped_closest_opportunities = total_closest_opportunities
    bootstrapped_singleton_differences = poisson(pavg*bootstrapped_closest_rates*total_singleton_opportunities)
    bootstrapped_singleton_opportunities = total_singleton_opportunities
    bootstrapped_pavg = bootstrapped_singleton_differences.sum()*1.0/bootstrapped_singleton_opportunities.sum()/(bootstrapped_closest_differences.sum()*1.0/bootstrapped_closest_opportunities.sum())
    
    bootstrapped_plowers = []
    bootstrapped_puppers = []
    bootstrapped_LRTs = []
    for dstar in dstars:
    
        lower_idxs = (bootstrapped_closest_rates<=dstar)
        upper_idxs = (bootstrapped_closest_rates>dstar)

        bootstrapped_plower = bootstrapped_singleton_differences[lower_idxs].sum()*1.0/bootstrapped_singleton_opportunities[lower_idxs].sum()/(bootstrapped_closest_differences[lower_idxs].sum()*1.0/bootstrapped_closest_opportunities[lower_idxs].sum())

        bootstrapped_pupper = bootstrapped_singleton_differences[upper_idxs].sum()*1.0/bootstrapped_singleton_opportunities[upper_idxs].sum()/(bootstrapped_closest_differences[upper_idxs].sum()*1.0/bootstrapped_closest_opportunities[upper_idxs].sum())

        bootstrapped_LRT = bootstrapped_singleton_differences[upper_idxs].sum()*log(bootstrapped_pupper/bootstrapped_pavg) + (bootstrapped_pavg-bootstrapped_pupper)*(bootstrapped_closest_differences[upper_idxs].sum()) + bootstrapped_singleton_differences[lower_idxs].sum()*log(bootstrapped_plower/bootstrapped_pavg) + (bootstrapped_pavg-bootstrapped_plower)*(bootstrapped_closest_differences[lower_idxs].sum())

        bootstrapped_plowers.append( bootstrapped_plower )
    
        bootstrapped_puppers.append( bootstrapped_pupper  )
        
        bootstrapped_LRTs.append(bootstrapped_LRT)
        
    bootstrapped_plowerss.append(bootstrapped_plowers)
    bootstrapped_pupperss.append(bootstrapped_puppers)
    bootstrapped_pavgss.append(bootstrapped_pavg*numpy.ones(len(bootstrapped_plowers)))
    bootstrapped_LRTss.append( bootstrapped_LRTs )

bootstrapped_plowerss = numpy.array(bootstrapped_plowerss)    
bootstrapped_pupperss = numpy.array(bootstrapped_pupperss)
bootstrapped_pavgss = numpy.array(bootstrapped_pavgss)   
bootstrapped_LRTss = numpy.array(bootstrapped_LRTs)

print bootstrapped_plowerss.mean(axis=0)

print bootstrapped_LRTss.mean()
print LRTs

pvalues = (bootstrapped_LRTss>=(LRTs[None,:])).sum(axis=0)*1.0/num_bootstraps   
    
print pvalues
    



 