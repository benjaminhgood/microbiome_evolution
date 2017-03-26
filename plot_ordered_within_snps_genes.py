import matplotlib  
matplotlib.use('Agg') 
import parse_midas_data
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

mpl.rcParams['font.size'] = 6
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
    sys.stderr.write("Usage: python command.py [debug] species_name")
################################################################################

min_change = 0.8
min_coverage = 20
alpha = 0.5 # Confidence interval range for rate estimates


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
samples, total_pis, total_pi_opportunities = parse_midas_data.parse_within_sample_pi(species_name, debug)
sys.stderr.write("Done!\n")
pis = total_pis/total_pi_opportunities
clipped_pis = (total_pis+1)/(total_pi_opportunities+1)

median_coverages = numpy.array([sample_coverage_map[samples[i]] for i in xrange(0,len(samples))])

snp_samples = samples[(median_coverages>=min_coverage)]

# Load gene coverage information for species_name
sys.stderr.write("Loading pangenome data for %s...\n" % species_name)
gene_samples, gene_names, gene_presence_matrix, gene_depth_matrix, marker_coverages, gene_reads_matrix = parse_midas_data.parse_pangenome_data(species_name,allowed_samples=snp_samples)
sys.stderr.write("Done!\n")
   
upper_gene_numbers = gene_diversity_utils.calculate_gene_numbers(gene_depth_matrix, marker_coverages, min_copynum=0.1)
lower_gene_numbers = gene_diversity_utils.calculate_gene_numbers(gene_depth_matrix, marker_coverages, min_copynum=0.5)  

desired_samples = gene_samples[marker_coverages>=min_coverage]

# Calculate which pairs of idxs belong to the same sample, which to the same subject
# and which to different subjects
snp_sample_idx_map = parse_midas_data.calculate_sample_idx_map(desired_samples, snp_samples)
gene_sample_idx_map = parse_midas_data.calculate_sample_idx_map(desired_samples, gene_samples)
        
same_sample_pi_plowers = []
same_sample_pi_puppers = []
same_sample_lower_gene_numbers = []
same_sample_upper_gene_numbers = []

for desired_sample_idx in xrange(0,len(desired_samples)):

    snp_idx = snp_sample_idx_map[desired_sample_idx]
    gene_idx = gene_sample_idx_map[desired_sample_idx]
    
    plower,pupper = stats_utils.calculate_poisson_rate_interval(total_pis[snp_idx], total_pi_opportunities[snp_idx], alpha)
    
    same_sample_pi_plowers.append(plower)
    same_sample_pi_puppers.append(pupper)
    same_sample_lower_gene_numbers.append(lower_gene_numbers[gene_idx])
    same_sample_upper_gene_numbers.append(upper_gene_numbers[gene_idx])

# clip lower bounds 
same_sample_pi_plowers = numpy.clip(same_sample_pi_plowers,1e-09,1e09)
same_sample_lower_gene_numbers = numpy.clip(same_sample_lower_gene_numbers,1e-09,1e09)

# Sort lists by ascending lower bound on SNP changes
same_sample_pi_plowers, same_sample_lower_gene_numbers, same_sample_pi_puppers, same_sample_upper_gene_numbers = (numpy.array(x) for x in zip(*sorted(zip(same_sample_pi_plowers, same_sample_lower_gene_numbers, same_sample_pi_puppers, same_sample_upper_gene_numbers))))

# Done calculating... now plot figure!

# Set up figure
fig = plt.figure(figsize=(5, 5))

# Set up grids to hold figure panels
outer_grid = gridspec.GridSpec(1, 2, width_ratios=[1,1], wspace=0.1)

###################
#
# pi Panel
#
###################

pi_axis = plt.Subplot(fig, outer_grid[0])
fig.add_subplot(pi_axis)

pi_axis.set_ylabel('Samples')
pi_axis.set_xlabel('Within sample diversity')
pi_axis.set_xlim([1e-07,9e-02])

###################
#
# Gene number 
#
###################

gene_axis = plt.Subplot(fig, outer_grid[1])
fig.add_subplot(gene_axis)

gene_axis.set_xlabel('Num genes')
#gene_axis.set_xlim([1e-07,9e-02])

y = 0
for pi_plower, pi_pupper, gene_lower, gene_upper in zip(same_sample_pi_plowers, same_sample_pi_puppers, same_sample_lower_gene_numbers, same_sample_upper_gene_numbers):

    y-=1
    pi_axis.semilogx([pi_plower,pi_pupper], [y,y],'b.-',linewidth=0.25,markersize=1.5)
    gene_axis.plot([gene_lower, gene_upper], [y,y], 'b.-',linewidth=0.25, markersize=1.5)
        

pi_axis.set_ylim([y-1,0])    
gene_axis.set_ylim([y-1,0])
pi_axis.set_yticks([])
gene_axis.set_yticks([])

fig.savefig('%s/%s_ordered_within_snps_genes.pdf' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight')
fig.savefig('%s/%s_ordered_within_snps_genes.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight',dpi=300)
