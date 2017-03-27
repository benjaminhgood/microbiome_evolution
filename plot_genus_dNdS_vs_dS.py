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

if len(sys.argv)>1:
    if len(sys.argv)>2:
        debug=True # debug does nothing in this script
        genus_name=sys.argv[2]
    else:
        debug=False
        genus_name=sys.argv[1]
else:
    sys.stderr.write("Usage: python plot_pNpS_vs_pi.py [debug] species_name")
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

main_axis.set_xlabel('Sequence divergence')
main_axis.set_ylabel('Fraction nonsynonymous')
main_axis.set_xlim([1e-07,1e-01])
main_axis.set_ylim([-0.1,1.1])
 
main_axis.spines['top'].set_visible(False)
main_axis.spines['right'].set_visible(False)
main_axis.get_xaxis().tick_bottom()
main_axis.get_yaxis().tick_left()

main_axis.semilogx([1e-07,1e-01],[0,0],'k-',linewidth=0.25)
main_axis.semilogx([1e-07,1e-01],[1,1],'k-',linewidth=0.25)
main_axis.semilogx([1e-07,1e-01],[0.5,0.5],'k-',linewidth=0.25)

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
print species_names

# Minimum frequency change to count as a fixed difference
# TODO: change this to an argument
min_change = 0.8
# Minimum median coverage of sample to look at
min_coverage = 20

# Load subject and sample metadata
sys.stderr.write("Loading HMP metadata...\n")
subject_sample_map = parse_midas_data.parse_subject_sample_map()
sys.stderr.write("Done!\n")

for species_name in species_names:
    
    species_name_items = species_name.split("_")
    
    species_label = "_".join([species_name_items[0]]+species_name_items[1:])
    
    # Load genomic coverage distributions
    sample_coverage_histograms, samples = parse_midas_data.parse_coverage_distribution(species_name)
    median_coverages = numpy.array([stats_utils.calculate_nonzero_median_from_histogram(sample_coverage_histogram) for sample_coverage_histogram in sample_coverage_histograms])
    sample_coverage_map = {samples[i]: median_coverages[i] for i in xrange(0,len(samples))}

    # Load pi information for species_name
    sys.stderr.write("Loading within-sample diversity for %s...\n" % species_name)
    samples, total_pis, total_pi_opportunities = parse_midas_data.parse_within_sample_pi(species_name, debug=debug)
    sys.stderr.write("Done!\n")
    pis = total_pis/total_pi_opportunities

    median_coverages = numpy.array([sample_coverage_map[samples[i]] for i in xrange(0,len(samples))])

    # Only plot samples above a certain depth threshold that are "haploids"
    desired_samples = samples[(median_coverages>=min_coverage)*(pis<=1e-03)]

    if len(desired_samples) < 2:
        sys.stderr.write("Too few haploid samples for %s.\n" % species_name)
        continue
    else:
        sys.stderr.write("Analyzing %d haploid samples...\n" % len(desired_samples))

    # Load SNP information for species_name
    sys.stderr.write("Loading %s...\n" % species_name)
    dummy_samples, allele_counts_map, passed_sites_map = parse_midas_data.parse_snps(species_name, debug=debug, allowed_samples=desired_samples)
    sys.stderr.write("Done!\n")
    
    # Calculate fixation matrices
    sys.stderr.write("Calculating 4D fixation matrix...\n")
    fixation_matrix_syn, fixation_opportunities_syn = diversity_utils.calculate_fixation_matrix(allele_counts_map, passed_sites_map, allowed_variant_types=set(['4D']), min_change=min_change)
    sys.stderr.write("Calculating 1D fixation matrix...\n")
    fixation_matrix_non, fixation_opportunities_non = diversity_utils.calculate_fixation_matrix(allele_counts_map, passed_sites_map, allowed_variant_types=set(['1D']), min_change=min_change)
    sys.stderr.write("Calculating total fixation matrix...\n")
    fixation_matrix_all, fixation_opportunities_all = diversity_utils.calculate_fixation_matrix(allele_counts_map, passed_sites_map, min_change=min_change)

    sys.stderr.write("Done!\n")

    # Calculate fraction nonsynonymous  
    dN = fixation_matrix_non/fixation_opportunities_non
    dS = fixation_matrix_syn/fixation_opportunities_syn
    dNplusdS = (dN+dS)
    fraction_nonsynonymous = dN/(dNplusdS+(dNplusdS==0))

    # Calculate total divergence
    dtot = fixation_matrix_all/fixation_opportunities_all

    # Calculate which pairs of idxs belong to the same sample, which to the same subject
    # and which to different subjects
    same_sample_idxs, same_subject_idxs, diff_subject_idxs = parse_midas_data.calculate_subject_pairs(subject_sample_map, desired_samples)


    line, = main_axis.semilogx(dtot[diff_subject_idxs], fraction_nonsynonymous[diff_subject_idxs],'.',markersize=2,alpha=0.5,markeredgewidth=0.0)
    colorVal = pylab.getp(line,'color')
    main_axis.semilogx(dtot[same_subject_idxs], fraction_nonsynonymous[same_subject_idxs],'s',alpha=0.5,markersize=2,color=colorVal,markeredgewidth=0.0)
    legend_axis.plot([-2,-1],[-2,-1],'.',markersize=3,alpha=0.5,markeredgewidth=0.0, label=species_label)
    
legend_axis.legend(loc='upper center',frameon=False,fontsize=5,numpoints=1,ncol=1,handlelength=1)   

pylab.savefig('%s/%s_genus_dNdS_vs_dS.pdf' % (parse_midas_data.analysis_directory,genus_name),bbox_inches='tight')
pylab.savefig('%s/%s_genus_dNdS_vs_dS.png' % (parse_midas_data.analysis_directory,genus_name),bbox_inches='tight')

    
