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
from numpy.random import randint, random, choice, multinomial
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

twin_modification_difference_threshold = 1e06
twin_replacement_difference_threshold = 1e06

################################################################################

#####################
#
# Settings for calculation:
#
#####################

min_coverage = config.min_median_coverage
min_sample_size = 3
min_haploid_sample_size = 10

#####
#
# Settings for different cohorts we are looking at 
#
#####
cohorts = ["hmp", "twins", "young_twins"]
countries = ["United States", "United Kingdom", "Western Europe"]
country_cohort_map = {country: cohort for country,cohort in zip(countries,cohorts)}

modification_difference_thresholds = {"hmp": modification_difference_threshold, "twins": twin_modification_difference_threshold, "young_twins": twin_modification_difference_threshold}

replacement_difference_thresholds = {"hmp": replacement_difference_threshold, "twins": twin_replacement_difference_threshold, "young_twins": twin_replacement_difference_threshold}

################################
#
# Set up figures
#
################################


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


################################
#
# Now do calculation
#
################################

hmp_species_qp_counts = {}
twin_species_qp_counts = {}

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
    good_species_list = ["Bacteroides_vulgatus_57955", "Bacteroides_uniformis_57318"]

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

    #all_samples = list(haploid_samples)
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
                #print sample_i, sample_j   
            
            elif (sample_i not in haploid_samples) and (sample_j not in haploid_samples):
                # pair that is non-QP at both timepoints
                qp_counts[cohort][2] += 1
            
            else:
                # pair that is QP at one timepoint and non-QP at another
                qp_counts[cohort][3] += 1
    
    sys.stderr.write("Done!\n")
    
    for cohort in cohorts:    
        print ("%s:" % cohort), qp_counts[cohort][1], "QP pairs,", qp_counts[cohort][2], "non-QP pairs,", qp_counts[cohort][3], "mixed pairs", qp_counts[cohort][0], "species dropouts."
   
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
    
    sys.stderr.write("Proceeding with %d HMP longitudinal comparisons!\n" % (sample_size))
    

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

hmp_haploid_axis.set_xlim([0,250])
twin_haploid_axis.set_xlim([0,140])


### Do stuff for legend
hmp_haploid_axis.barh([-10],[1],linewidth=0,label='QP->QP', color='#08519c')
hmp_haploid_axis.barh([-10],[1],linewidth=0,label='non->non', color='#de2d26')
hmp_haploid_axis.barh([-10],[1],linewidth=0,label='mixed', color='#8856a7')
hmp_haploid_axis.barh([-10],[1],linewidth=0,label='dropout', color='0.7')
hmp_haploid_axis.legend(loc='lower right',frameon=False)


sys.stderr.write("Saving figures...\t")
fig6.savefig('%s/supplemental_temporal_haploid.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight',transparent=True)



    
