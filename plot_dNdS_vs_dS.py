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
        species_name=sys.argv[2]
    else:
        debug=False
        species_name=sys.argv[1]
else:
    sys.stderr.write("Usage: python plot_pNpS_vs_pi.py [debug] species_name")
################################################################################

# Minimum frequency change to count as a fixed difference
# TODO: change this to an argument
min_change = 0.8
# Minimum median coverage of sample to look at
min_coverage = 20

# Load subject and sample metadata
sys.stderr.write("Loading HMP metadata...\n")
subject_sample_map = parse_midas_data.parse_subject_sample_map()
sys.stderr.write("Done!\n")
    
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

pylab.figure(figsize=(3.42,2))
pylab.xlabel('Sequence divergence')
pylab.ylabel('Fraction nonsynonymous')
pylab.ylim([-0.1,1.1])
pylab.xlim([1e-07,1e-01])
pylab.title(species_name,fontsize=7)

pylab.gca().spines['top'].set_visible(False)
pylab.gca().spines['right'].set_visible(False)
pylab.gca().get_xaxis().tick_bottom()
pylab.gca().get_yaxis().tick_left()


pylab.semilogx([1e-07,1e-01],[0,0],'k-',linewidth=0.25)
pylab.semilogx([1e-07,1e-01],[1,1],'k-',linewidth=0.25)
pylab.semilogx([1e-07,1e-01],[0.5,0.5],'k-',linewidth=0.25)

pylab.semilogx(dtot[diff_subject_idxs], fraction_nonsynonymous[diff_subject_idxs],'r.',markersize=3,alpha=0.5)

pylab.semilogx(dtot[same_subject_idxs], fraction_nonsynonymous[same_subject_idxs],'gs',alpha=0.5,markersize=3)

pylab.savefig('%s/%s_dNdS_vs_dS.pdf' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight')
pylab.savefig('%s/%s_dNdS_vs_dS.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight')

    
