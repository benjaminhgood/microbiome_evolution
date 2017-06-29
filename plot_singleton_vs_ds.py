import matplotlib  
matplotlib.use('Agg') 
import config
import parse_midas_data
import parse_HMP_data


import pylab
import sys
import numpy

import diversity_utils
import gene_diversity_utils

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


mpl.rcParams['font.size'] = 6
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
parser.add_argument("--chunk-size", type=int, help="max number of records to load", default=1000000000)
parser.add_argument('--other-species', type=str, help='Run the script for a different species')

args = parser.parse_args()

debug = args.debug
chunk_size = args.chunk_size
other_species = args.other_species

if other_species:
    species_name = other_species
    other_species_str = "_%s" % species_name
else:
    other_species_str = ""
    
################################################################################

min_coverage = config.min_median_coverage
alpha = 0.5 # Confidence interval range for rate estimates
low_pi_threshold = 1e-03
low_divergence_threshold = 1e-03
min_change = 0.8

# Load subject and sample metadata
sys.stderr.write("Loading sample metadata...\n")
subject_sample_map = parse_HMP_data.parse_subject_sample_map()
sample_country_map = parse_HMP_data.parse_sample_country_map()
sample_phenotype_map = parse_HMP_data.parse_sample_phenotype_map()
sys.stderr.write("Done!\n")
       
# Only plot samples above a certain depth threshold that are "haploids"
snp_samples = diversity_utils.calculate_haploid_samples(species_name, debug=debug)
# Only consider one sample per person
snp_samples = snp_samples[parse_midas_data.calculate_unique_samples(subject_sample_map, sample_list=snp_samples)]
sys.stderr.write("Proceeding with %d haploid samples!\n" % len(snp_samples))

####################################################
#
# Set up Figure (4 panels, arranged in 2x2 grid)
#
####################################################

pylab.figure(1,figsize=(3,3))
fig = pylab.gcf()
outer_grid  = gridspec.GridSpec(2,1,hspace=0.05,height_ratios=[1,1])

d_axis = plt.Subplot(fig, outer_grid[0])
fig.add_subplot(d_axis)

#d_axis.set_xlabel('Synonymous divergence, $d_S$')
d_axis.set_ylabel('Fraction 1D')
d_axis.set_xlim([1e-04,1e-02])

singleton_axis = plt.Subplot(fig, outer_grid[1])
fig.add_subplot(singleton_axis)

singleton_axis.set_xlabel('Synonymous divergence, $d_S$')
singleton_axis.set_ylabel('Fraction 1D')
singleton_axis.set_xlim([1e-04,1e-02])
singleton_axis.set_ylim([0,1])

##############################################################################
#
# Now do calculations
#
##############################################################################

sys.stderr.write("Loading core genes...\n")
core_genes = parse_midas_data.load_core_genes(species_name)
sys.stderr.write("Done! Core genome consists of %d genes\n" % len(core_genes))

# Analyze SNPs, looping over chunk sizes. 
# Clunky, but necessary to limit memory usage on cluster

# Load SNP information for species_name
sys.stderr.write("Loading SNPs for %s...\n" % species_name)
sys.stderr.write("(core genes only...)\n")



pi_matrix_syn = numpy.array([])
avg_pi_matrix_syn = numpy.array([])
snp_difference_matrix = numpy.array([])
snp_opportunity_matrix = numpy.array([])

singletons = [] # (idx, vartype)

final_line_number = 0
while final_line_number >= 0:
    
    sys.stderr.write("Loading chunk starting @ %d...\n" % final_line_number)
    dummy_samples, allele_counts_map, passed_sites_map, final_line_number = parse_midas_data.parse_snps(species_name, debug=debug, allowed_samples=snp_samples,allowed_genes=core_genes,chunk_size=chunk_size,initial_line_number=final_line_number)
    sys.stderr.write("Done! Loaded %d genes\n" % len(allele_counts_map.keys()))
    
    print len(dummy_samples), "dummy samples!"
    
    # Calculate fixation matrix
    sys.stderr.write("Calculating matrix of snp differences...\n")
    chunk_snp_difference_matrix, chunk_snp_opportunity_matrix = diversity_utils.calculate_fixation_matrix(allele_counts_map, passed_sites_map, min_change=min_change, allowed_variant_types=set(['4D']))   # 
    sys.stderr.write("Done!\n")
    
    if snp_difference_matrix.shape[0]==0:
        snp_difference_matrix = numpy.zeros_like(chunk_snp_difference_matrix)*1.0
        snp_opportunity_matrix = numpy.zeros_like(snp_difference_matrix)*1.0
    
    snp_difference_matrix += chunk_snp_difference_matrix
    snp_opportunity_matrix += chunk_snp_opportunity_matrix

    sys.stderr.write("Calculating singletons...\n")
    chunk_singletons = diversity_utils.calculate_singletons(allele_counts_map, passed_sites_map, allowed_genes=core_genes)
    singletons.extend(chunk_singletons)
    sys.stderr.write("Done!\n")

    snp_samples = dummy_samples

snp_substitution_rate = snp_difference_matrix*1.0/(snp_opportunity_matrix+(snp_opportunity_matrix==0))
snp_substitution_rate = numpy.clip(snp_substitution_rate,1e-09,10)

snp_substitution_rate += numpy.identity(len(snp_substitution_rate))*10

sys.stderr.write("Postprocessing vs and ds...\n")
ds = []
vs = []
for idx, variant_type in singletons:
    
    if variant_type=='1D':
        v = 1
    elif variant_type=='4D':
        v = 0
    else:
        continue
        
    d = snp_substitution_rate[idx,:].min()
    
    ds.append(d)
    vs.append(v)

ds = numpy.array(ds)
vs = numpy.array(vs)

sys.stderr.write("Done!\n")   

print len(ds), "total singletons"
print (vs>0.5).sum(), "1D"
print (vs<0.5).sum(), "4D"


# Now plot them. 

xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(ds)
d_axis.step(xs,1-ns*1.0/ns[0],'-')

dstars = numpy.logspace(-4,-2,20)
fraction_nonsynonymous = []

for dstar in dstars:

    less_idxs = (ds<=dstar)
    
    if less_idxs.sum() > 1:
        # some of them!
        
        fraction_nonsynonymous.append(vs[less_idxs].mean())
        
    else:
        fraction_nonsynonymous.append(-1)    

d_axis.set_xticks([])
singleton_axis.semilogx(dstars, fraction_nonsynonymous, 'k.-')


sys.stderr.write("Saving figure...\t")
fig.savefig('%s/singleton_vs_dS%s.pdf' % (parse_midas_data.analysis_directory, other_species_str),bbox_inches='tight')

sys.stderr.write("Done!\n")

 