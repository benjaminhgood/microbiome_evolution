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

min_coverage = config.min_median_coverage
alpha = 0.5 # Confidence interval range for rate estimates
low_pi_threshold = 1e-03
low_divergence_threshold = 1e-03
min_change = 0.8
min_sample_size = 33 # 46 gives at least 1000 pairs
allowed_variant_types = set(['1D','2D','3D','4D'])

snp_divergences = {}
gene_divergences = {}


good_species_list = parse_midas_data.parse_good_species_list()
if debug:
    good_species_list = good_species_list[0:2]
    

good_species_list = ['Bacteroides_vulgatus_57955', 'Bacteroides_uniformis_57318', 'Alistipes_putredinis_61533']

sys.stderr.write("Loading sample metadata...\n")
subject_sample_map = parse_HMP_data.parse_subject_sample_map()
sample_order_map = parse_HMP_data.parse_sample_order_map()
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
    sys.stderr.write("Calculating snp matrix...\n")
    dummy_samples, snp_difference_matrix, snp_opportunity_matrix = calculate_substitution_rates.calculate_matrices_from_substitution_rate_map(substitution_rate_map, 'core', allowed_samples=snp_samples)
    snp_samples = dummy_samples
    sys.stderr.write("Done!\n")
    
    sys.stderr.write("Calculating gene matrix...\n")
    gene_samples, gene_difference_matrix, gene_opportunity_matrix = calculate_substitution_rates.calculate_matrices_from_substitution_rate_map(substitution_rate_map, 'genes', allowed_samples=snp_samples)
    snp_samples = gene_samples
    sys.stderr.write("Done!\n")
    
    desired_samples = gene_samples
    
    desired_same_sample_idxs, desired_same_subject_idxs, desired_diff_subject_idxs =  parse_midas_data.calculate_subject_pairs( subject_sample_map, desired_samples)

    snp_sample_idx_map = parse_midas_data.calculate_sample_idx_map(desired_samples, snp_samples)
    gene_sample_idx_map = parse_midas_data.calculate_sample_idx_map(desired_samples, gene_samples)
  
    same_subject_snp_idxs = parse_midas_data.apply_sample_index_map_to_indices(snp_sample_idx_map, desired_same_subject_idxs)  
    same_subject_gene_idxs = parse_midas_data.apply_sample_index_map_to_indices(gene_sample_idx_map, desired_same_subject_idxs)  

    diff_subject_snp_idxs = parse_midas_data.apply_sample_index_map_to_indices(snp_sample_idx_map, desired_diff_subject_idxs)  
    diff_subject_gene_idxs = parse_midas_data.apply_sample_index_map_to_indices(gene_sample_idx_map, desired_diff_subject_idxs)  
    
    
    snp_divergences[species_name] = []
    gene_divergences[species_name] = []
    
    for sample_pair_idx in xrange(0,len(diff_subject_snp_idxs[0])):
    
        snp_i = diff_subject_snp_idxs[0][sample_pair_idx]
        snp_j = diff_subject_snp_idxs[1][sample_pair_idx]
        
        # Now do genes
        gene_i = diff_subject_gene_idxs[0][sample_pair_idx]
        gene_j = diff_subject_gene_idxs[1][sample_pair_idx]
        
        if (snp_opportunity_matrix[snp_i,snp_j]>0) and (gene_opportunity_matrix[gene_i, gene_j]>0):
        
            snp_d = snp_difference_matrix[snp_i,snp_j]*1.0/snp_opportunity_matrix[snp_i,snp_j]
            gene_d = gene_difference_matrix[gene_i, gene_j]*1.0/gene_opportunity_matrix[gene_i, gene_j]
            
            if (snp_d < 1e-04) and (gene_d > 1e-02):
                print snp_samples[snp_i], snp_samples[snp_j], gene_samples[gene_i], gene_samples[gene_j], snp_difference_matrix[snp_i, snp_j], snp_opportunity_matrix[snp_i, snp_j], gene_difference_matrix[gene_i, gene_j], gene_opportunity_matrix[gene_i, gene_j]
            
            snp_divergences[species_name].append(snp_d)
            gene_divergences[species_name].append(gene_d)

species_names = []
sample_sizes = []

for species_name in snp_divergences.keys():
    species_names.append(species_name)
    sample_sizes.append( len(snp_divergences[species_name]) )
    
# sort in descending order of sample size
# Sort by num haploids    
sample_sizes, species_names = zip(*sorted(zip(sample_sizes, species_names),reverse=True))
    
sys.stderr.write("Postprocessing %d species...\n" % len(species_names))
        

####################################################
#
# Set up Figure (1 panels, arranged in 1x1 grid)
#
####################################################

haploid_color = '#08519c'

pylab.figure(1,figsize=(3.42,2))
fig = pylab.gcf()
# make three panels panels
outer_grid  = gridspec.GridSpec(1,1)

divergence_axis = plt.Subplot(fig, outer_grid[0])
fig.add_subplot(divergence_axis)

divergence_axis.set_ylabel('Gene divergence')
divergence_axis.set_xlabel('SNP divergence')

divergence_axis.set_xlim([1e-06,1e-01])
divergence_axis.set_ylim([1e-02,1])

divergence_axis.spines['top'].set_visible(False)
divergence_axis.spines['right'].set_visible(False)
divergence_axis.get_xaxis().tick_bottom()
divergence_axis.get_yaxis().tick_left()

# Plot percentiles of divergence distribution
for species_idx in xrange(0,len(species_names)):

    species_name = species_names[species_idx]

    sys.stderr.write("Postprocessing %s...\n" % (species_name))
    
    snp_divergences[species_name] = numpy.array(snp_divergences[species_name])
    gene_divergences[species_name] = numpy.array(gene_divergences[species_name])
    
    if species_name.startswith('Bacteroides_vulgatus'):
        pass
        #divergence_axis.loglog(snp_divergences[species_name], gene_divergences[species_name], 'r.', markersize=2,markeredgewidth=0,zorder=1,label=("%s" % species_name))
    
    divergence_axis.loglog(snp_divergences[species_name], gene_divergences[species_name], '.', markersize=2,alpha=0.5,markeredgewidth=0,zorder=0,label=species_name)

divergence_axis.legend(loc='lower left', frameon=False, numpoints=1,fontsize=4)

sys.stderr.write("Saving figure...\t")
fig.savefig('%s/supplemental_genes_vs_snps.png' % (parse_midas_data.analysis_directory),bbox_inches='tight',dpi=600)
sys.stderr.write("Done!\n")

 