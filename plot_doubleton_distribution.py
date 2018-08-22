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
import matplotlib.colors as mcolors
import matplotlib.patheffects as pe

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

min_opportunities = 10

data = {}
within_host_data = {}

good_species_list = parse_midas_data.parse_good_species_list()
if debug:
    good_species_list = good_species_list[0:2]

# Load subject and sample metadata
sys.stderr.write("Loading sample metadata...\n")
sample_utils.parse_subject_sample_map
subject_sample_map = sample_utils.parse_subject_sample_map()
sample_order_map = sample_utils.parse_sample_order_map()
sample_country_map = sample_utils.parse_sample_country_map()
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
    snp_samples = haploid_samples[sample_utils.calculate_unique_samples(subject_sample_map, sample_list=haploid_samples)]

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
    substitution_rate_matrix = difference_matrix*1.0/(opportunity_matrix+(opportunity_matrix==0))
    sys.stderr.write("Done!\n")
    
    good_idxs = (doubleton_opportunity_matrix>0.5)
    
    doubletons = doubleton_matrix[good_idxs]
    doubleton_opportunities = doubleton_opportunity_matrix[good_idxs]
    substitution_rates = substitution_rate_matrix[good_idxs]
    
    data[species_name] = doubletons, doubleton_opportunities, substitution_rates
    
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

pylab.figure(1,figsize=(3.42,5))
fig = pylab.gcf()

outer_grid  = gridspec.GridSpec(2,1, height_ratios=[1,1],hspace=0.1)

sharing_axis = plt.Subplot(fig, outer_grid[0])
fig.add_subplot(sharing_axis)

sharing_axis.set_ylabel('Fraction host pairs $\geq p$')
#sharing_axis.set_xlabel('Marker SNV sharing, $p$') 
sharing_axis.set_xlim([-0.025,1.025])
sharing_axis.set_ylim([0,1])

sharing_axis.spines['top'].set_visible(False)
sharing_axis.spines['right'].set_visible(False)
sharing_axis.get_xaxis().tick_bottom()
sharing_axis.get_yaxis().tick_left()

within_sharing_axis = plt.Subplot(fig, outer_grid[1])
fig.add_subplot(within_sharing_axis)

within_sharing_axis.set_ylabel('Fraction timepoints $\geq p$')
within_sharing_axis.set_xlabel('Marker SNV sharing, $p$') 
within_sharing_axis.set_xlim([-0.025,1.025])
within_sharing_axis.set_ylim([0,1])

within_sharing_axis.spines['top'].set_visible(False)
within_sharing_axis.spines['right'].set_visible(False)
within_sharing_axis.get_xaxis().tick_bottom()
within_sharing_axis.get_yaxis().tick_left()


ds = numpy.logspace(-5,-1,50)

cumulative_doubletons = numpy.zeros_like(ds)
cumulative_doubleton_opportunities = numpy.zeros_like(ds)

low_doubletons = []
low_doubleton_opportunities = []
all_doubletons = []
all_doubleton_opportunities = []

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
    
       

low_doubletons = numpy.array(low_doubletons)
low_doubleton_opportunities = numpy.array(low_doubleton_opportunities)
all_doubletons = numpy.array(all_doubletons)
all_doubleton_opportunities = numpy.array(all_doubleton_opportunities)

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

bootstrapped_low_ps = []
bootstrapped_all_ps = []
bootstrapped_fake_low_ps = []
bootstrapped_fake_all_ps = []
real_all_ps = all_doubletons*1.0/all_doubleton_opportunities 
real_low_ps = low_doubletons*1.0/low_doubleton_opportunities 

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

#xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(bootstrapped_low_ps, min_x=0,max_x=2)
#sharing_axis.step(xs,ns*1.0/ns[0],'r-',label='Low $d_S$ (matched)',zorder=3)

#xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(bootstrapped_all_ps, min_x=0,max_x=2)
#sharing_axis.step(xs,ns*1.0/ns[0],'k-',label='All (matched)',zorder=2)


#xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(bootstrapped_fake_low_ps, min_x=0,max_x=1)
#sharing_axis.step(xs,ns*1.0/ns[0],'r-',label='Low $d_S$ (pooled)',zorder=1,alpha=0.5)


xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(real_all_ps[all_doubleton_opportunities>min_opportunities], min_x=0,max_x=2)
sharing_axis.step(xs,ns*1.0/ns[0],'k-',label='Between hosts (all)',zorder=1) #,alpha=0.5)

xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(real_low_ps[low_doubleton_opportunities>min_opportunities], min_x=0,max_x=2)
sharing_axis.step(xs,ns*1.0/ns[0],'r-',label='Between hosts\n(closely related)',zorder=1) #,alpha=0.5)


#xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(replacement_ps, min_x=0,max_x=2)
#sharing_axis.step(xs,ns*1.0/ns[0],'b-',label='Within-host (rep)',zorder=3,alpha=0.5)


#xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(within_ps, min_x=0,max_x=2)
#sharing_axis.step(xs,ns*1.0/ns[0],'b-',label='Within-host (mod)',zorder=3)



#xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(bootstrapped_fake_all_ps, min_x=0,max_x=1)
#sharing_axis.step(xs,ns*1.0/ns[0],'k-',label='All (pooled)',zorder=0,alpha=0.5)

sharing_axis.legend(loc='upper center',frameon=False)

#######################
#
# Within-host version
#
#######################

sys.stderr.write("Within hosts...\n")

cohorts = ["hmp", "twins"] #, "young_twins"]
countries = ["United States", "United Kingdom"] #, "Western Europe"]
country_cohort_map = {country: cohort for country,cohort in zip(countries,cohorts)}

modification_difference_thresholds = {"hmp": config.modification_difference_threshold, "twins": config.twin_modification_difference_threshold, "young_twins": config.twin_modification_difference_threshold}

replacement_difference_thresholds = {"hmp": config.replacement_difference_threshold, "twins": config.twin_replacement_difference_threshold, "young_twins": config.twin_replacement_difference_threshold}

modification_shared_snps = {cohort: [] for cohort in cohorts}
modification_opportunities = {cohort: [] for cohort in cohorts}
replacement_shared_snps = {cohort: [] for cohort in cohorts}
replacement_opportunities = {cohort: [] for cohort in cohorts}
    
    
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
       
    if len(haploid_samples) < config.within_host_min_haploid_sample_size:
        continue


    same_sample_idxs, same_subject_idxs, diff_subject_idxs = sample_utils.calculate_ordered_subject_pairs(sample_order_map, all_samples)

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
            
            elif (sample_i not in haploid_samples) and (sample_j not in haploid_samples):
                # pair that is non-QP at both timepoints
                qp_counts[cohort][2] += 1
            
            else:
                # pair that is QP at one timepoint and non-QP at another
                qp_counts[cohort][3] += 1
    
    sys.stderr.write("Done!\n")
   
    combined_sample_set = set()
    for cohort in cohorts:
        combined_sample_set.update(qp_sample_sets[cohort])
    combined_samples = list(sorted(combined_sample_set))
    combined_sample_idx_map = {combined_samples[i] : i for i in xrange(0,len(combined_samples))}    
    
    qp_sample_lists = {cohort: list(sorted(qp_sample_sets[cohort])) for cohort in cohorts}
    
    sample_size = len(qp_sample_sets['hmp'])
        
    if sample_size < config.within_host_min_sample_size:
        continue
            
    import calculate_private_snvs
    private_snv_map = calculate_private_snvs.load_private_snv_map(species_name)
    
    import calculate_snp_prevalences
    snv_freq_map = calculate_snp_prevalences.parse_population_freqs(species_name,polarize_by_consensus=True)
    snv_freq_values = snv_freq_map.values()
    
    sys.stderr.write("Loading pre-computed temporal changes...\n")
    temporal_change_map = calculate_temporal_changes.load_temporal_change_map(species_name)
    sys.stderr.write("Done!\n")

    ### Now loop over different cohorts
    for cohort in cohorts:
        
        modification_difference_threshold = modification_difference_thresholds[cohort]
        replacement_difference_threshold = replacement_difference_thresholds[cohort]
        
        desired_samples = qp_sample_lists[cohort]
        
        same_sample_idxs, same_subject_idxs, diff_subject_idxs = sample_utils.calculate_ordered_subject_pairs(sample_order_map, desired_samples)
        
        for sample_pair_idx in xrange(0,len(same_subject_idxs[0])):
        
            
            sample_i = desired_samples[same_subject_idxs[0][sample_pair_idx]] 
            sample_j = desired_samples[same_subject_idxs[1][sample_pair_idx]]
       
            L, perr, mutations, reversions = calculate_temporal_changes.calculate_mutations_reversions_from_temporal_change_map(temporal_change_map, sample_i, sample_j)
        
            nerr = L*perr
        
            num_mutations = len(mutations)
            num_reversions = len(reversions)
            num_snp_changes = num_mutations+num_reversions
        
            gene_L, gene_perr, gains, losses = calculate_temporal_changes.calculate_gains_losses_from_temporal_change_map(temporal_change_map, sample_i, sample_j)
        
            gene_nerr = gene_L*gene_perr
            num_gains = len(gains)
            num_losses = len(losses)
            num_gene_changes = num_gains+num_losses
        
            if (perr<-0.5) or (gene_perr < -0.5):
                continue
        
            if (nerr > max([0.5, 0.1*num_snp_changes])) or (gene_nerr > max([0.5, 0.1*num_gene_changes])):
                continue # Only take things with low-ish FPR
        
            if num_snp_changes==0:
                continue
            
            private_L, private_perr, private_reversions = calculate_temporal_changes.calculate_private_reversions_from_temporal_change_map(temporal_change_map, sample_i, sample_j)
            
            private_nerr = private_perr*private_L
        
            if L==0 or private_L==0:
                continue
            
            if private_nerr > 0.5:
                continue
        
            # print raw numbers
            #print private_L, len(private_reversions) 
            
            if num_snp_changes <= modification_difference_threshold:
                # modification!
                modification_shared_snps[cohort].append(private_L-len(private_reversions))
                modification_opportunities[cohort].append(private_L)
                
            elif num_snp_changes >= replacement_difference_threshold:
                # replacement!
                replacement_shared_snps[cohort].append(private_L-len(private_reversions))
                replacement_opportunities[cohort].append(private_L)
            
            else:
                # weird middle ground
                pass

modification_ps = {}
replacement_ps = {}
for cohort in cohorts:
    modification_shared_snps[cohort] = numpy.array(modification_shared_snps[cohort])
    modification_opportunities[cohort] = numpy.array(modification_opportunities[cohort])
    modification_ps[cohort] = modification_shared_snps[cohort]*1.0/(modification_opportunities[cohort]+(modification_opportunities[cohort]==0))
    
    
    replacement_shared_snps[cohort] = numpy.array(replacement_shared_snps[cohort])
    replacement_opportunities[cohort] = numpy.array(replacement_opportunities[cohort])
    replacement_ps[cohort] = replacement_shared_snps[cohort]*1.0/(replacement_opportunities[cohort]+(replacement_opportunities[cohort]==0))

xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(replacement_ps['hmp'][replacement_opportunities['hmp']>min_opportunities], min_x=0,max_x=2)
within_sharing_axis.step(xs,ns*1.0/ns[0],'-', color='#08519c',linewidth=1, label='Within-host\nreplacements',zorder=2,path_effects=[pe.Stroke(linewidth=5, foreground='#fee0d2'), pe.Normal()])

xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(modification_ps['hmp'][modification_opportunities['hmp']>min_opportunities], min_x=0,max_x=2)
within_sharing_axis.step(xs,ns*1.0/ns[0],'-',color='#08519c',linewidth=1, label='Within-host\nmodifications',zorder=2,path_effects=[pe.Stroke(linewidth=5, foreground='#9ecae1'), pe.Normal()]) 


xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(replacement_ps['twins'][replacement_opportunities['twins']>min_opportunities], min_x=0,max_x=2)
within_sharing_axis.step(xs,ns*1.0/ns[0],'-',label='Twin replacements',zorder=1,color='#8856a7')

xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(modification_ps['twins'][modification_opportunities['twins']>min_opportunities], min_x=0,max_x=2)
within_sharing_axis.step(xs,ns*1.0/ns[0],'-',label='Twin other',zorder=1,alpha=0.5,color='#8856a7')

within_sharing_axis.legend(loc='center',frameon=False)


sys.stderr.write("Saving figure...\t")
fig.savefig('%s/supplemental_low_divergence_marker_sharing.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')
sys.stderr.write("Done!\n")
sys.exit(0)


 