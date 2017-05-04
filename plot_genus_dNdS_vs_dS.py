import matplotlib  
matplotlib.use('Agg') 
import parse_midas_data
import pylab
import sys
import numpy
import diversity_utils
import stats_utils
import matplotlib.colors as colors
import matplotlib.cm as cmx
from math import log10,ceil
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy.random import randint

mpl.rcParams['font.size'] = 7
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
parser.add_argument("genus_name", help="name of genus to process")
parser.add_argument("--debug", help="Loads only a subset of SNPs for speed", action="store_true")
parser.add_argument("--include-china", help="Includes Chinese subjects from Qin et al (Nature, 2012)", action="store_true")
parser.add_argument("--chunk-size", type=int, help="max number of records to load", default=1000000000)
args = parser.parse_args()

genus_name = args.genus_name
debug = args.debug
chunk_size = args.chunk_size
include_china = args.include_china
################################################################################

##
# Set up figure
# 
#

# Set up figure
fig = plt.figure(figsize=(7, 2.5))

# Set up grids to hold figure panels
outer_grid = gridspec.GridSpec(1, 2, width_ratios=[3,2], wspace=0.1)

main_axis = plt.Subplot(fig, outer_grid[0])
fig.add_subplot(main_axis)

main_axis.set_xlabel('Synonymous divergence, $d_S$')
main_axis.set_ylabel('$dN/dS$ for samples $\leq d_S$')
main_axis.set_xlim([1e-05,1e-01])
main_axis.set_ylim([1e-02,1e01])
 
main_axis.spines['top'].set_visible(False)
main_axis.spines['right'].set_visible(False)
main_axis.get_xaxis().tick_bottom()
main_axis.get_yaxis().tick_left()

#main_axis.semilogx([1e-07,1e-01],[0,0],'k-',linewidth=0.25)
main_axis.semilogx([1e-07,1e-01],[1,1],'-',linewidth=0.25,color='0.7')
#main_axis.semilogx([1e-07,1e-01],[0.5,0.5],'k-',linewidth=0.25)

legend_axis = plt.Subplot(fig, outer_grid[1])
fig.add_subplot(legend_axis)

legend_axis.set_ylim([0,1])
legend_axis.set_xlim([0,1])

legend_axis.spines['top'].set_visible(False)
legend_axis.spines['right'].set_visible(False)
legend_axis.spines['left'].set_visible(False)
legend_axis.spines['bottom'].set_visible(False)

legend_axis.set_xticks([])
legend_axis.set_yticks([])  


good_species_names = parse_midas_data.parse_good_species_list()
species_names = []
for species_name in good_species_names:
    if species_name.startswith(genus_name):
        species_names.append(species_name)

species_names = species_names[0:10]
if debug:
    species_names = species_names[0:2]
print species_names

# Minimum frequency change to count as a fixed difference
# TODO: change this to an argument
min_change = 0.8
# Minimum median coverage of sample to look at
min_coverage = 20

if include_china:
    allowed_countries = set([])
else:
    allowed_countries = set(['United States'])

# Load subject and sample metadata
sys.stderr.write("Loading HMP metadata...\n")
subject_sample_map = parse_midas_data.parse_subject_sample_map()
sample_country_map = parse_midas_data.parse_sample_country_map()
sys.stderr.write("Done!\n")

dss = numpy.logspace(-5,-1,50)

genus_fixation_non = numpy.zeros_like(dss)
genus_opportunities_non = numpy.zeros_like(dss)
               
genus_fixation_syn = numpy.zeros_like(dss)
genus_opportunities_syn = numpy.zeros_like(dss)
        
def calculate_dNdS(n_non, non_sites, n_syn, syn_sites):
    return (n_non*(syn_sites+1.0)/(non_sites+1.0)+1.0)/(n_syn+1.0)

for species_name in species_names:
    
    species_name_items = species_name.split("_")
    
    species_label = "_".join([species_name_items[0]]+species_name_items[1:])
    
    
    
    # Load core gene set
    sys.stderr.write("Loading core genes...\n")
    core_genes = parse_midas_data.load_core_genes(species_name)
    sys.stderr.write("Done! Core genome consists of %d genes\n" % len(core_genes))
 
    
    # Load genomic coverage distributions
    sample_coverage_histograms, samples = parse_midas_data.parse_coverage_distribution(species_name)
    median_coverages = numpy.array([stats_utils.calculate_nonzero_median_from_histogram(sample_coverage_histogram) for sample_coverage_histogram in sample_coverage_histograms])
    sample_coverage_map = {samples[i]: median_coverages[i] for i in xrange(0,len(samples))}

    # Load pi information for species_name
    sys.stderr.write("Loading within-sample diversity for %s...\n" % species_name)
    samples, total_pis, total_pi_opportunities = parse_midas_data.parse_within_sample_pi(species_name, allowed_genes=core_genes, debug=debug)
    sys.stderr.write("Done!\n")
    pis = total_pis/total_pi_opportunities

    median_coverages = numpy.array([sample_coverage_map[samples[i]] for i in xrange(0,len(samples))])

    # Only plot samples above a certain depth threshold that are "haploids"
    desired_samples = samples[(median_coverages>=min_coverage)*(pis<=1e-03)]
    # Restrict to single timepoint single timepoints per person
    unique_subject_idxs = parse_midas_data.calculate_unique_samples(subject_sample_map, desired_samples)
    allowed_country_idxs = parse_midas_data.calculate_country_samples(sample_country_map, desired_samples, allowed_countries=set(['United States']))
    desired_samples = desired_samples[unique_subject_idxs*allowed_country_idxs]
    
    if len(desired_samples) < 2:
        sys.stderr.write("Too few haploid samples for %s.\n" % species_name)
        continue
    else:
        sys.stderr.write("Analyzing %d haploid samples...\n" % len(desired_samples))

    # Load SNP information for species_name
    sys.stderr.write("Loading %s...\n" % species_name)
    
    fixation_matrix_syn = numpy.array([])
    fixation_opportunities_syn = numpy.array([])
    fixation_matrix_non = numpy.array([])
    fixation_opportunities_non = numpy.array([])
    
    
    final_line_number = 0
    while final_line_number >= 0:
    
        sys.stderr.write("Loading chunk starting @ %d...\n" % final_line_number)
        dummy_samples, allele_counts_map, passed_sites_map, final_line_number = parse_midas_data.parse_snps(species_name, debug=debug, allowed_variant_types=set(['1D','4D']), allowed_genes=core_genes, allowed_samples=desired_samples,chunk_size=chunk_size,initial_line_number=final_line_number)
        sys.stderr.write("Done! Loaded %d genes\n" % len(allele_counts_map.keys()))
    
        # Calculate fixation matrix
        sys.stderr.write("Calculating matrix of synonymous differences...\n")
        chunk_snp_difference_matrix, chunk_snp_opportunity_matrix = diversity_utils.calculate_fixation_matrix(allele_counts_map, passed_sites_map, allowed_variant_types=set(['4D']), allowed_genes=core_genes, min_change=min_change)    
        sys.stderr.write("Done!\n")
    
        if fixation_matrix_syn.shape[0]==0:
            fixation_matrix_syn = numpy.zeros_like(chunk_snp_difference_matrix)*1.0
            fixation_opportunities_syn = numpy.zeros_like(fixation_matrix_syn)*1.0
    
        fixation_matrix_syn += chunk_snp_difference_matrix
        fixation_opportunities_syn += chunk_snp_opportunity_matrix

        # Calculate fixation matrix
        sys.stderr.write("Calculating matrix of nonsynonymous differences...\n")
        chunk_snp_difference_matrix, chunk_snp_opportunity_matrix = diversity_utils.calculate_fixation_matrix(allele_counts_map, passed_sites_map, allowed_variant_types=set(['1D']), allowed_genes=core_genes, min_change=min_change)    
        sys.stderr.write("Done!\n")
    
        if fixation_matrix_non.shape[0]==0:
            fixation_matrix_non = numpy.zeros_like(chunk_snp_difference_matrix)*1.0
            fixation_opportunities_non = numpy.zeros_like(fixation_matrix_syn)*1.0
    
        fixation_matrix_non += chunk_snp_difference_matrix
        fixation_opportunities_non += chunk_snp_opportunity_matrix
  
    
    sys.stderr.write("Done!\n")

    # Calculate fraction nonsynonymous  
    dN = fixation_matrix_non/fixation_opportunities_non
    dS = fixation_matrix_syn/fixation_opportunities_syn
    dNplusdS = (dN+dS)
    fraction_nonsynonymous = dN/(dNplusdS+(dNplusdS==0))

    # Calculate synonymous divergence
    dsyn = fixation_matrix_syn/fixation_opportunities_syn
    
    # Calculate which pairs of idxs belong to the same sample, which to the same subject
    # and which to different subjects
    same_sample_idxs, same_subject_idxs, diff_subject_idxs = parse_midas_data.calculate_subject_pairs(subject_sample_map, desired_samples)

    flattened_fixation_matrix_non = fixation_matrix_non[diff_subject_idxs]
    flattened_fixation_opportunities_non = fixation_opportunities_non[diff_subject_idxs]
    flattened_fixation_matrix_syn = fixation_matrix_syn[diff_subject_idxs]
    flattened_fixation_opportunities_syn = fixation_opportunities_syn[diff_subject_idxs]
    
    flattened_dsyn = dsyn[diff_subject_idxs]

    # Calculate fraction nonsynonymous for divergences <= ds
    fraction_nonsynonymouss = []
    total_dNplusdSs = []
    
    total_fixation_non = numpy.zeros_like(dss)
    total_opportunities_non = numpy.zeros_like(dss)
    total_fixation_syn = numpy.zeros_like(dss)
    total_opportunities_syn = numpy.zeros_like(dss)
    
    
    for i in xrange(0,len(dss)):
        ds = dss[i]
        
        desired_idxs = (flattened_dsyn<=ds)

        total_fixation_non[i] = flattened_fixation_matrix_non[desired_idxs].sum()
        total_opportunities_non[i] = flattened_fixation_opportunities_non[desired_idxs].sum()
        total_fixation_syn[i] = flattened_fixation_matrix_syn[desired_idxs].sum()
        total_opportunities_syn[i] = flattened_fixation_opportunities_syn[desired_idxs].sum()
        
    total_fixation_non = numpy.array(total_fixation_non)
    total_opportunities_non = numpy.array(total_opportunities_non)
    total_fixation_syn = numpy.array(total_fixation_syn)
    total_opportunities_syn = numpy.array(total_opportunities_syn)
    total_fixation = total_fixation_non+total_fixation_syn
    
    dNdS = calculate_dNdS(total_fixation_non, total_opportunities_non, total_fixation_syn, total_opportunities_syn)
    
    genus_fixation_non += total_fixation_non
    genus_opportunities_non += total_opportunities_non
    genus_fixation_syn += total_fixation_syn
    genus_opportunities_syn += total_opportunities_syn
    
    line, = main_axis.semilogx(dss[total_fixation>0], dNdS[total_fixation>0],'-', alpha=0.5)
    colorVal = pylab.getp(line,'color')
    legend_axis.plot([-2,-1],[-2,-1],'-',alpha=0.5, label=species_label)

genus_fixation = genus_fixation_non+genus_fixation_syn

genus_dNdS = calculate_dNdS(genus_fixation_non, genus_opportunities_non, genus_fixation_syn, genus_opportunities_syn)
line, = main_axis.loglog(dss[genus_fixation>0], genus_dNdS[genus_fixation>0],'k-', linewidth=1)
colorVal = pylab.getp(line,'color')
legend_axis.plot([-2,-1],[-2,-1],'k-', label='All',linewidth=1)
   
    
legend_axis.legend(loc='upper center',frameon=False,fontsize=5,numpoints=1,ncol=1,handlelength=1)   

pylab.savefig('%s/%s_genus_dNdS_vs_dS.pdf' % (parse_midas_data.analysis_directory,genus_name),bbox_inches='tight')
pylab.savefig('%s/%s_genus_dNdS_vs_dS.png' % (parse_midas_data.analysis_directory,genus_name),bbox_inches='tight')

    
