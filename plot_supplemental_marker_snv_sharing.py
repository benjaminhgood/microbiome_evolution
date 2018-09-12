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
import calculate_singletons
import calculate_temporal_changes

import stats_utils
import matplotlib.colors as colors
import matplotlib.cm as cmx
from math import log10,ceil,log
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy.random import randint, shuffle, poisson, binomial, choice, hypergeometric

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
within_host_data = {}

good_species_list = parse_midas_data.parse_good_species_list()
if debug:
    good_species_list = good_species_list[0:2]

sys.stderr.write("Loading sample metadata...\n")
subject_sample_map = parse_HMP_data.parse_subject_sample_map()
sample_order_map = parse_HMP_data.parse_sample_order_map()
sys.stderr.write("Done!\n")

for species_name in good_species_list:

    sys.stderr.write("Loading haploid samples...\n")
    # Only plot samples above a certain depth threshold that are "haploids"
    haploid_samples = diversity_utils.calculate_haploid_samples(species_name, debug=debug)
    
    if len(haploid_samples) < min_sample_size:
        sys.stderr.write("Not enough haploid samples!\n")
        continue
        
    sys.stderr.write("Calculating unique samples...\n")
    # Only consider one sample per person
    snp_samples = haploid_samples[parse_midas_data.calculate_unique_samples(subject_sample_map, sample_list=haploid_samples)]

    if len(snp_samples) < min_sample_size:
        sys.stderr.write("Not enough unique samples!\n")
        continue

    # Load singleton matrices 
    sys.stderr.write("Loading pre-computed singleton rates for %s...\n" % species_name)
    singleton_rate_map = calculate_singletons.load_singleton_rate_map(species_name)
    sys.stderr.write("Calculating matrix...\n")
    snp_samples, singleton_matrix, doubleton_matrix, difference_matrix, opportunity_matrix = calculate_singletons.calculate_matrices_from_singleton_rate_map(singleton_rate_map, 'core', allowed_samples=snp_samples)
    snp_samples, syn_singleton_matrix, syn_doubleton_matrix, syn_difference_matrix, syn_opportunity_matrix  = calculate_singletons.calculate_matrices_from_singleton_rate_map(singleton_rate_map, '4D', allowed_samples=snp_samples)
    snp_samples, non_singleton_matrix, non_doubleton_matrix, non_difference_matrix, non_opportunity_matrix  = calculate_singletons.calculate_matrices_from_singleton_rate_map(singleton_rate_map, '1D', allowed_samples=snp_samples)
    
    doubleton_opportunity_matrix = singleton_matrix+doubleton_matrix
    substitution_rate_matrix = syn_difference_matrix*1.0/(syn_opportunity_matrix+(syn_opportunity_matrix==0))
    sys.stderr.write("Done!\n")
    
    good_idxs = (doubleton_opportunity_matrix>0.5)
    
    doubletons = doubleton_matrix[good_idxs]
    doubleton_opportunities = doubleton_opportunity_matrix[good_idxs]
    substitution_rates = substitution_rate_matrix[good_idxs]
    
    data[species_name] = doubletons, doubleton_opportunities, substitution_rates

    shared_snps = []
    shared_snp_opportunities = []
    replacement_shared_snps = []
    replacement_shared_snp_opportunities = []
    
    # get temporal samples
    same_sample_idxs, same_subject_idxs, diff_subject_idxs = parse_midas_data.calculate_ordered_subject_pairs(sample_order_map, haploid_samples)
    
    sys.stderr.write("Loading pre-computed temporal changes for %s...\n" % species_name)
    temporal_change_map = calculate_temporal_changes.load_temporal_change_map(species_name)
    sys.stderr.write("Done!\n")
    
    if len(temporal_change_map)==0:
        continue
    
    for sample_pair_idx in xrange(0,len(same_subject_idxs[0])):
   
        i = same_subject_idxs[0][sample_pair_idx]
        j = same_subject_idxs[1][sample_pair_idx]

        sample_i = haploid_samples[i]
        sample_j = haploid_samples[j]
        
        L, perr, mutations, reversions = calculate_temporal_changes.calculate_mutations_reversions_from_temporal_change_map(temporal_change_map, sample_i, sample_j)
        
        total_num_changes = len(mutations)+len(reversions)
        
        private_L, private_perr, private_reversions = calculate_temporal_changes.calculate_private_reversions_from_temporal_change_map(temporal_change_map, sample_i, sample_j)
        
        if L==0 or private_L==0:
            continue
            
        if total_num_changes==0:
            continue
            
        if private_L*private_perr > 0.5:
            continue
        
        print private_L, len(private_reversions) 
            
        if total_num_changes>config.modification_difference_threshold:
            #print "Replacement!"
            replacement_shared_snp_opportunities.append(private_L)
            replacement_shared_snps.append(private_L-len(private_reversions))
        else:   
            shared_snp_opportunities.append(private_L)
            shared_snps.append(private_L-len(private_reversions))
            
            
    within_host_data[species_name] = shared_snps, shared_snp_opportunities, replacement_shared_snps, replacement_shared_snp_opportunities
    
    
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

doubleton_axis = plt.Subplot(fig, outer_grid[0])
fig.add_subplot(doubleton_axis)

doubleton_axis.set_ylabel('Mean marker SNV sharing')
doubleton_axis.set_xlabel('Synonymous divergence, $d_S$') 
doubleton_axis.set_xlim([1e-06,1e-01])
doubleton_axis.set_ylim([0,1])

doubleton_axis.spines['top'].set_visible(False)
doubleton_axis.spines['right'].set_visible(False)
doubleton_axis.get_xaxis().tick_bottom()
doubleton_axis.get_yaxis().tick_left()

pylab.figure(2,figsize=(3.42,2.5))
fig2 = pylab.gcf()
# make three panels panels
outer_grid = gridspec.GridSpec(1,1)
#outer_grid  = gridspec.GridSpec(2,1, height_ratios=[1,1],hspace=0.2)

sharing_axis = plt.Subplot(fig2, outer_grid[0])
fig2.add_subplot(sharing_axis)

sharing_axis.set_ylabel('Fraction host pairs $\geq p$')
sharing_axis.set_xlabel('Marker SNV sharing, $p$') 
sharing_axis.set_xlim([-0.025,1.025])
sharing_axis.set_ylim([0,1])

sharing_axis.spines['top'].set_visible(False)
sharing_axis.spines['right'].set_visible(False)
sharing_axis.get_xaxis().tick_bottom()
sharing_axis.get_yaxis().tick_left()

ds = numpy.logspace(-5,-1,50)

cumulative_doubletons = numpy.zeros_like(ds)
cumulative_doubleton_opportunities = numpy.zeros_like(ds)

low_doubletons = []
low_doubleton_opportunities = []
all_doubletons = []
all_doubleton_opportunities = []

within_shared_snps = []
within_shared_snp_opportunities = []
replacement_shared_snps = []
replacement_shared_snp_opportunities = []

for species_idx in xrange(0,len(species_names)):

    species_name = species_names[species_idx]

    doubletons, doubleton_opportunities, substitution_rates = data[species_name]
    
    substitution_rates = substitution_rates*1e-02/(substitution_rates.mean())
     
    # Add to cumulative distribution
    for d_idx in xrange(0,len(ds)):
        d = ds[d_idx]
        
        lower_idxs = (substitution_rates<=d)
        
        cumulative_doubletons[d_idx] += doubletons[lower_idxs].sum()
        cumulative_doubleton_opportunities[d_idx] += doubleton_opportunities[lower_idxs].sum()
        
        
        
    cumulative_doubleton_rates = cumulative_doubletons*1.0/(cumulative_doubleton_opportunities+(cumulative_doubleton_opportunities==0))
    
    lower_idxs = (substitution_rates<=low_divergence_threshold)
    higher_idxs = numpy.logical_not(lower_idxs)
    
    low_doubletons.extend(doubletons[lower_idxs])
    low_doubleton_opportunities.extend(doubleton_opportunities[lower_idxs] )
    all_doubletons.extend(doubletons[higher_idxs])
    all_doubleton_opportunities.extend(doubleton_opportunities[higher_idxs] )
    
    if species_name in within_host_data:
        
        shared_snps, shared_snp_opportunities, replacement_snps, replacement_opportunities = within_host_data[species_name]    

        #print replacement_snps
        #print replacement_opportunities

        within_shared_snps.extend(shared_snps)
        within_shared_snp_opportunities.extend(shared_snp_opportunities)
        replacement_shared_snps.extend(replacement_snps)
        replacement_shared_snp_opportunities.extend(replacement_opportunities)
        

low_doubletons = numpy.array(low_doubletons)
low_doubleton_opportunities = numpy.array(low_doubleton_opportunities)
all_doubletons = numpy.array(all_doubletons)
all_doubleton_opportunities = numpy.array(all_doubleton_opportunities)

within_shared_snps = numpy.array(within_shared_snps)
within_shared_snp_opportunities = numpy.array(within_shared_snp_opportunities)

replacement_shared_snps = numpy.array(replacement_shared_snps)
replacement_shared_snp_opportunities = numpy.array(replacement_shared_snp_opportunities)

#print replacement_shared_snps
#print replacement_shared_snp_opportunities    
    
bootstrapped_cumulative_doubleton_ratess = []
num_bootstraps = 1000
for bootstrap_idx in xrange(0,num_bootstraps):

    # resampe everything at known rates
    
    bootstrapped_doubletons = poisson(cumulative_doubletons)*1.0
    bootstrapped_singletons = poisson(cumulative_doubleton_opportunities-cumulative_doubletons)*1.0
    bootstrapped_doubleton_opportunities = bootstrapped_doubletons+bootstrapped_singletons
    
    
    bootstrapped_cumulative_doubleton_ratess.append( bootstrapped_doubletons/(bootstrapped_doubleton_opportunities+(bootstrapped_doubleton_opportunities==0)) )
    
bootstrapped_cumulative_doubleton_ratess = numpy.array(bootstrapped_cumulative_doubleton_ratess)
avg_rates = bootstrapped_cumulative_doubleton_ratess.mean(axis=0)
std_rates = bootstrapped_cumulative_doubleton_ratess.std(axis=0)

doubleton_axis.fill_between(ds, avg_rates-2*std_rates, avg_rates+2*std_rates,color='0.7',linewidth=0)

doubleton_axis.semilogx(ds, avg_rates, 'k-')
doubleton_axis.semilogx(ds[cumulative_doubleton_opportunities>0], cumulative_doubleton_rates[cumulative_doubleton_opportunities>0], 'r-')

bootstrapped_low_ps = []
bootstrapped_all_ps = []
bootstrapped_fake_low_ps = []
bootstrapped_fake_all_ps = []
real_all_ps = all_doubletons*1.0/all_doubleton_opportunities 
real_low_ps = low_doubletons*1.0/low_doubleton_opportunities 
within_ps = within_shared_snps*1.0/within_shared_snp_opportunities
replacement_ps = replacement_shared_snps*1.0/replacement_shared_snp_opportunities

num_bootstraps = 10
for bootstrap_idx in xrange(0,num_bootstraps):

    # resampe everything at known rates
    
    idxs = choice(numpy.arange(0,len(all_doubletons)),size=len(low_doubletons))
    
    sample_sizes = numpy.fmin(low_doubleton_opportunities, all_doubleton_opportunities[idxs]).astype(numpy.int32)
    
    low_ngood = low_doubletons.astype(numpy.int32)
    low_nbad = (low_doubleton_opportunities-low_doubletons).astype(numpy.int32)
    
    low_p = low_doubletons.sum()*1.0/low_doubleton_opportunities.sum()
    
    all_ngood = all_doubletons[idxs].astype(numpy.int32)
    all_nbad = (all_doubleton_opportunities[idxs] - all_ngood).astype(numpy.int32)
    all_p = all_doubletons.sum()*1.0/all_doubleton_opportunities.sum()
    
    bootstrapped_low_ps.extend( hypergeometric(low_ngood, low_nbad, sample_sizes)*1.0/sample_sizes )
    bootstrapped_all_ps.extend( hypergeometric(all_ngood, all_nbad, sample_sizes)*1.0/sample_sizes )
    bootstrapped_fake_low_ps.extend( binomial(sample_sizes, low_p)*1.0/sample_sizes )
    bootstrapped_fake_all_ps.extend( binomial(sample_sizes, all_p)*1.0/sample_sizes )

xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(bootstrapped_low_ps, min_x=0,max_x=2)
sharing_axis.step(xs,ns*1.0/ns[0],'r-',label='Low $d_S$ (matched)',zorder=3)

xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(bootstrapped_all_ps, min_x=0,max_x=2)
sharing_axis.step(xs,ns*1.0/ns[0],'k-',label='All (matched)',zorder=2)


#xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(bootstrapped_fake_low_ps, min_x=0,max_x=1)
#sharing_axis.step(xs,ns*1.0/ns[0],'r-',label='Low $d_S$ (pooled)',zorder=1,alpha=0.5)

xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(real_low_ps, min_x=0,max_x=2)
sharing_axis.step(xs,ns*1.0/ns[0],'r-',label='Low $d_S$',zorder=1,alpha=0.5)

xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(real_all_ps, min_x=0,max_x=2)
sharing_axis.step(xs,ns*1.0/ns[0],'k-',label='All',zorder=1,alpha=0.5)

xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(replacement_ps, min_x=0,max_x=2)
sharing_axis.step(xs,ns*1.0/ns[0],'b-',label='Within-host (rep)',zorder=3,alpha=0.5)


xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(within_ps, min_x=0,max_x=2)
sharing_axis.step(xs,ns*1.0/ns[0],'b-',label='Within-host (mod)',zorder=3)



#xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(bootstrapped_fake_all_ps, min_x=0,max_x=1)
#sharing_axis.step(xs,ns*1.0/ns[0],'k-',label='All (pooled)',zorder=0,alpha=0.5)



sharing_axis.legend(loc='upper center',frameon=False)

sys.stderr.write("Saving figure...\t")
fig.savefig('%s/supplemental_marker_sharing.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')
fig2.savefig('%s/supplemental_low_divergence_marker_sharing.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')
sys.stderr.write("Done!\n")
sys.exit(0)


 