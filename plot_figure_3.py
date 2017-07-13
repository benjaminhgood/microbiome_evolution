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

args = parser.parse_args()

debug = args.debug
chunk_size = args.chunk_size
memoize = args.memoize

################################################################################

min_coverage = config.min_median_coverage
alpha = 0.5 # Confidence interval range for rate estimates
low_pi_threshold = 1e-03
low_divergence_threshold = 1e-03
min_change = 0.8
min_sample_size = 46
allowed_variant_types = set(['1D','2D','3D','4D'])

intermediate_filename = '%sfigure_3_data.txt' % (parse_midas_data.analysis_directory)


divergence_matrices = {}

if memoize and os.path.isfile(intermediate_filename):


    # load from Disk!  
    file = open(intermediate_filename,"r")
    file.readline() # header
    for line in file:
        items = line.split(";")
        species_name = items[0]
        divergence_matrix = []
        for item in items[1:]:
            divergence_matrix.append([float(subitem) for subitem in  item.strip().split(",")])
            
        divergence_matrix = numpy.array(divergence_matrix)
        divergence_matrices[species_name] = divergence_matrix

    file.close()
    
else:
    
    # Must compete divergence matrix on the fly! 
            
    # Load subject and sample metadata
    sys.stderr.write("Loading sample metadata...\n")
    subject_sample_map = parse_HMP_data.parse_subject_sample_map()
    sys.stderr.write("Done!\n")
    
    good_species_list = parse_midas_data.parse_good_species_list()

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

        sys.stderr.write("Proceeding with %d haploid samples!\n" % len(snp_samples))

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
        final_line_number = 0

        while final_line_number >= 0:
    
            sys.stderr.write("Loading chunk starting @ %d...\n" % final_line_number)
            dummy_samples, allele_counts_map, passed_sites_map, final_line_number = parse_midas_data.parse_snps(species_name, debug=debug, allowed_samples=snp_samples,allowed_genes=core_genes,allowed_variant_types=allowed_variant_types, chunk_size=chunk_size,initial_line_number=final_line_number)
            sys.stderr.write("Done! Loaded %d genes\n" % len(allele_counts_map.keys()))
    
            
            # Calculate fixation matrix
            sys.stderr.write("Calculating matrix of snp differences...\n")
            chunk_snp_difference_matrix, chunk_snp_opportunity_matrix =     diversity_utils.calculate_fixation_matrix(allele_counts_map, passed_sites_map, min_change=min_change)   # 
            sys.stderr.write("Done!\n")
    
            if snp_difference_matrix.shape[0]==0:
                snp_difference_matrix = numpy.zeros_like(chunk_snp_difference_matrix)*1.0
                snp_opportunity_matrix = numpy.zeros_like(snp_difference_matrix)*1.0
    
            snp_difference_matrix += chunk_snp_difference_matrix
            snp_opportunity_matrix += chunk_snp_opportunity_matrix

        snp_substitution_rate = snp_difference_matrix*1.0/(snp_opportunity_matrix+(snp_opportunity_matrix==0))
        snp_substitution_rate[snp_opportunity_matrix==0] = -1
        
        if snp_substitution_rate.shape[0] < min_sample_size:
            continue
            
        divergence_matrices[species_name] = snp_substitution_rate
        
        sys.stderr.write("Done with %s!\n" % species_name) 
    
    sys.stderr.write("Done looping over species!\n")
    
    sys.stderr.write("Writing intermediate file...\n")
    file = open(intermediate_filename,"w")
    # header!
    file.write("Species; d11, d12; d21, d22\n")
    for species_name in divergence_matrices.keys():
        
    
        matrix_strs = []
        for i in xrange(0,divergence_matrices[species_name].shape[0]):
            matrix_strs.append(", ".join([str(d) for d in  divergence_matrices[species_name][:,i]])) 

        output_strs = [species_name]+matrix_strs
        output_str = "; ".join(output_strs)
        file.write("%s\n" % output_str)
    file.close()
    
    
# 
species_names = []
sample_sizes = []

for species_name in divergence_matrices.keys():
    species_names.append(species_name)
    sample_sizes.append( divergence_matrices[species_name].shape[0] )
    
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

pylab.figure(1,figsize=(7,1))
fig = pylab.gcf()
# make three panels panels
outer_grid  = gridspec.GridSpec(1,1)

divergence_axis = plt.Subplot(fig, outer_grid[0])
fig.add_subplot(divergence_axis)

divergence_axis.set_ylabel('Divergence, $d$')
divergence_axis.set_ylim([1e-06,1e-01])
divergence_axis.set_xlim([-1,len(species_names)])

xticks = numpy.arange(0,len(species_names))
#xticklabels = ["%s (%d)" % (species_names[i],sample_sizes[i]) for i in xrange(0,len(sample_sizes))]
xticklabels = ["%s" % (species_names[i]) for i in xrange(0,len(sample_sizes))]

divergence_axis.set_xticks(xticks)
divergence_axis.set_xticklabels(xticklabels, rotation='vertical',fontsize=4)

divergence_axis.spines['top'].set_visible(False)
divergence_axis.spines['right'].set_visible(False)
divergence_axis.get_xaxis().tick_bottom()
divergence_axis.get_yaxis().tick_left()


# Plot percentiles of divergence distribution
for species_idx in xrange(0,len(species_names)):

    species_name = species_names[species_idx]

    sys.stderr.write("Postprocessing %s (%d samples)...\n" % (species_name, divergence_matrices[species_name].shape[0]))
    divergence_matrix = divergence_matrices[species_name]
    divergences = []
    for i in xrange(0, divergence_matrix.shape[0]):
        for j in xrange(i+1, divergence_matrix.shape[0]):
            
            if divergence_matrix[i,j] >= 0:
                
                divergences.append(divergence_matrix[i,j])
    
    divergences = numpy.array(divergences)
    divergences = numpy.clip(divergences,1e-06,1)
    divergences.sort() # ascending by default
    
    log_divergences = numpy.log(divergences)
    
    kernel = gaussian_kde(log_divergences)
    
    n = len(divergences)
    
    
    
    #percentiles = numpy.array([0.001,0.01,0.1,0.25,0.5,0.75,0.9,0.99,0.999])
    percentiles = numpy.array([0.001,0.01,0.1])
    
    quantiles = numpy.array([divergences[long(n*p)] for p in percentiles])
    quantiles = numpy.clip(quantiles,1e-06,1)
    
    theory_log_divergences = numpy.linspace(log_divergences.min(), log_divergences.max()+1,100)
    theory_divergences = numpy.exp(theory_log_divergences)
    theory_pdf = kernel(theory_log_divergences)
    theory_pdf = theory_pdf / theory_pdf.max() * 0.45
    
    
    divergence_axis.fill_betweenx(theory_divergences, species_idx-theory_pdf, species_idx+theory_pdf,linewidth=0,facecolor='#de2d26')
    divergence_axis.semilogy(numpy.ones_like(quantiles)*species_idx, quantiles,'k_',markersize=3)

    

sys.stderr.write("Saving figure...\t")
fig.savefig('%s/figure_3.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')
sys.stderr.write("Done!\n")

 