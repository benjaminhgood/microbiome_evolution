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

num_bootstraps = 1000

min_coverage = config.min_median_coverage
alpha = 0.05 # Confidence interval range for rate estimates
low_divergence_threshold = 2e-04
min_change = 0.8
min_sample_size = 33 # 46 gives at least 1000 pairs, 33 gives at least 500 (actually 528)
allowed_variant_types = set(['1D','2D','3D','4D'])

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
sample_subject_map = sample_utils.calculate_sample_subject_map(subject_sample_map)
sample_country_map = sample_utils.parse_sample_country_map()
sys.stderr.write("Done!\n")

all_closest_rates = []
all_pair_rates = []

for species_name in good_species_list:

    sys.stderr.write("Loading haploid samples...\n")
    # Only plot samples above a certain depth threshold that are "haploids"
    haploid_samples = diversity_utils.calculate_haploid_samples(species_name, debug=debug)
    
    snp_samples = []
    for sample_name in haploid_samples:
        if sample_country_map[sample_name]=='United Kingdom':
            snp_samples.append(sample_name)

    if len(snp_samples) < 10:
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
    
    closest_snp_substitution_rates = []
    pair_snp_substitution_rates = []
    for i in xrange(0, snp_opportunity_matrix.shape[0]):
    
        min_substitution_rate = 1e09
    
        for j in xrange(0, snp_opportunity_matrix.shape[0]):
            
            # Don't check yourself
            if i==j:
                continue
        
            
            if snp_opportunity_matrix[i,j]>0.5:
            
                if sample_subject_map[snp_samples[i]] == sample_subject_map[snp_samples[j]]:
                    pair_snp_substitution_rates.append( snp_substitution_matrix[i,j] )
            
                if snp_substitution_matrix[i,j]<min_substitution_rate:
                    min_substitution_rate = snp_substitution_matrix[i,j]
                                    
        closest_snp_substitution_rates.append(min_substitution_rate)
        
    all_closest_rates.extend(closest_snp_substitution_rates) 
    all_pair_rates.extend(pair_snp_substitution_rates)
    
print numpy.sort(all_closest_rates)
print numpy.sort(all_pair_rates)
    
xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(all_closest_rates, min_x=1e-06, max_x=1e09)
pylab.step(xs,ns/ns[0],'-',color='r',linewidth=0.5, alpha=0.5, label='Between-host', where='mid',zorder=2)
    
xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(all_pair_rates, min_x=1e-06, max_x=1e09)
pylab.step(xs,ns/ns[0],'-',color='b',linewidth=0.5, alpha=0.5, label='Twin', where='mid',zorder=2)

pylab.semilogx([1e-07],[1])
pylab.xlim([1e-06,1e-02])
pylab.ylim([0,1])
pylab.savefig('twin_test.pdf',bbox_inches='tight')