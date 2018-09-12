import matplotlib  
matplotlib.use('Agg') 
import sample_utils
import config
import parse_midas_data
import os.path
import pylab
import sys
import numpy

import species_phylogeny_utils
import diversity_utils
import gene_diversity_utils
import calculate_substitution_rates
import clade_utils
import figure_utils

import stats_utils
import matplotlib.colors as colors
import matplotlib.cm as cmx
from math import log10,ceil
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy.random import randint, shuffle

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
low_divergence_threshold = config.between_low_divergence_threshold
min_sample_size = config.between_host_min_sample_size # 46 gives at least 1000 pairs, 33 gives at least 500 (actually 528)
allowed_variant_types = set(['1D','2D','3D','4D'])

num_bootstraps = 10000

clade_Fst = {}
Fst = {}
clade_country_likelihood = {}

divergence_matrices = {}
good_species_list = parse_midas_data.parse_good_species_list()

if debug:
    good_species_list = good_species_list[0:2]

divergence_matrices = {}
sample_names = {}
    
# Load subject and sample metadata
sys.stderr.write("Loading sample metadata...\n")
subject_sample_map = sample_utils.parse_subject_sample_map()
sample_country_map = sample_utils.parse_sample_country_map()
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
    sys.stderr.write("Done!\n")

    snp_substitution_matrix = snp_difference_matrix*1.0/(snp_opportunity_matrix+(snp_opportunity_matrix==0))
    
    divergence_matrices[species_name] = snp_substitution_matrix
    sample_names[species_name] = snp_samples
    
    # Load manually annotated clades
    clade_sets = clade_utils.load_manual_clades(species_name)

    if len(clade_sets)==0:
        continue
    
    clade_idxss = clade_utils.calculate_clade_idxs_from_clade_sets(snp_samples, clade_sets)
        
        
    nonsingleton_clade_idxss = []
    for clade_idxs in clade_idxss:
        if clade_idxs.sum() > 1:
            nonsingleton_clade_idxss.append(clade_idxs)
           
    # Want at least two clades!
    if len(nonsingleton_clade_idxss)>=2:
        all_nonsingleton_clade_idxs = numpy.array([False for x in nonsingleton_clade_idxss[0]])
        for clade_idxs in nonsingleton_clade_idxss:
            all_nonsingleton_clade_idxs = numpy.logical_or(all_nonsingleton_clade_idxs, clade_idxs)
    else:
        all_nonsingleton_clade_idxs = []    
    
    # First do Fst with clades!
    if len(nonsingleton_clade_idxss)>=2:
        # Calculate Fst between manually defined clades
        
        clade_substitution_matrices = []
        clade_ones_matrices = []
        clade_pair_idxss = []
        for clade_idxs in nonsingleton_clade_idxss:
            
            clade_substitution_matrices.append( snp_substitution_matrix[numpy.ix_(numpy.nonzero(clade_idxs)[0], numpy.nonzero(clade_idxs)[0])] )
            
            clade_ones_matrices.append( numpy.ones_like(clade_substitution_matrices[-1]) )
            
            clade_pair_idxss.append(  numpy.triu_indices(clade_substitution_matrices[-1].shape[0], 1) )
         
        all_substitution_matrix =  snp_substitution_matrix[ numpy.ix_(numpy.nonzero(all_nonsingleton_clade_idxs)[0], numpy.nonzero(all_nonsingleton_clade_idxs)[0]) ]
        
        all_ones_matrix = numpy.ones_like(all_substitution_matrix)
        
        all_pair_idxs = numpy.triu_indices(all_substitution_matrix.shape[0], 1)     
    
        within_numerator = sum([clade_substitution_matrices[i][clade_pair_idxss[i]].sum() for i in xrange(0,len(clade_substitution_matrices))])
        
        within_denominator = sum([clade_ones_matrices[i][clade_pair_idxss[i]].sum() for i in xrange(0,len(clade_substitution_matrices))])

        within_rate = within_numerator*1.0/within_denominator
        
        between_rate = (all_substitution_matrix[all_pair_idxs].sum())/(all_ones_matrix[all_pair_idxs].sum())
        
        observed_fst = 1.0 - within_rate/between_rate
        
        clade_Fst[species_name] = (observed_fst, [])

        

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
for species_name in sorted_species_names:
    species_names.append(species_name)
    sample_sizes.append( divergence_matrices[species_name].shape[0] )
        
sys.stderr.write("Postprocessing %d species...\n" % len(species_names))        

####################################################
#
# Set up Figure (3 panels, arranged in 3x1 grid)
#
####################################################

haploid_color = '#08519c'

pylab.figure(1,figsize=(4,1))
fig = pylab.gcf()
# make three panels panels
outer_grid  = gridspec.GridSpec(1,1)

clade_fst_axis = plt.Subplot(fig, outer_grid[0])
fig.add_subplot(clade_fst_axis)

clade_fst_axis.set_ylabel('Fst (clades)')
clade_fst_axis.set_ylim([-0.05,1.05])
clade_fst_axis.set_xlim([-1,len(species_names)])


clade_fst_axis.spines['top'].set_visible(False)
clade_fst_axis.spines['right'].set_visible(False)
clade_fst_axis.get_xaxis().tick_bottom()
clade_fst_axis.get_yaxis().tick_left()

xticks = numpy.arange(0,len(species_names))
#xticklabels = ["%s (%d)" % (species_names[i],sample_sizes[i]) for i in xrange(0,len(sample_sizes))]
xticklabels = ["%s" % figure_utils.get_pretty_species_name(species_names[i]) for i in xrange(0,len(species_names))]

clade_fst_axis.set_xticks(xticks)
clade_fst_axis.set_xticklabels(xticklabels, rotation='vertical',fontsize=4)



# Plot percentiles of divergence distribution
for species_idx in xrange(0,len(species_names)):

    species_name = species_names[species_idx]

    #clade_fst_axis.plot([-1,len(species_names)+1],[0.2,0.2],'k:',linewidth=0.5)
    
    if species_name in clade_Fst:
        observed_fst, bootstrapped_fsts = clade_Fst[species_name]
    
        clade_fst_axis.plot([species_idx], [observed_fst],'r^',markersize=3,markeredgewidth=0)        

sys.stderr.write("Saving figure...\t")
fig.savefig('%s/supplemental_clade_city_correlation.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')
sys.stderr.write("Done!\n")

 