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

# Calculate which pairs of idxs belong to the same sample, which to the same subject
# and which to different subjects
same_sample_idxs, same_subject_idxs, diff_subject_idxs = parse_midas_data.calculate_subject_pairs(subject_sample_map, samples)

# Calculate the smaller and larger of the two pi estimates so we can look at correlation over time
lower_pis = numpy.fmin(clipped_pis[same_subject_idxs[0]],clipped_pis[same_subject_idxs[1]])
upper_pis = numpy.fmax(clipped_pis[same_subject_idxs[0]],clipped_pis[same_subject_idxs[1]])

# Only plot samples above a certain depth threshold that are "haploids"
desired_samples = samples[(median_coverages>=min_coverage)*(pis<=1e-03)]

# Load SNP information for species_name
sys.stderr.write("Loading %s...\n" % species_name)
dummy_samples, allele_counts_map, passed_sites_map = parse_midas_data.parse_snps(species_name, debug=debug, allowed_samples=desired_samples)
sys.stderr.write("Done!\n")
 
# Calculate fixation matrix
sys.stderr.write("Calculating matrix of snp differences...\n")
snp_difference_matrix, snp_opportunity_matrix = diversity_utils.calculate_fixation_matrix(allele_counts_map, passed_sites_map, min_change=min_change)    
sys.stderr.write("Done!\n")
       
# Calculate which pairs of idxs belong to the same sample, which to the same subject
# and which to different subjects
same_sample_snp_idxs, same_subject_snp_idxs, diff_subject_snp_idxs = parse_midas_data.calculate_subject_pairs(subject_sample_map, desired_samples)


same_sample_pi_plowers = []
same_sample_pi_puppers = []
for i in xrange(0,len(pis)):

    if median_coverages[i]<min_coverage:
        continue 
  
    plower,pupper = stats_utils.calculate_poisson_rate_interval(total_pis[i], total_pi_opportunities[i], alpha)
    
    same_sample_pi_plowers.append(plower)
    same_sample_pi_puppers.append(pupper)


# clip lower bounds 
same_sample_pi_plowers = numpy.clip(same_sample_pi_plowers,1e-09,1e09)
# Sort both lists by ascending lower bound on SNP changes
same_sample_pi_plowers, same_sample_pi_puppers = (numpy.array(x) for x in zip(*sorted(zip(same_sample_pi_plowers, same_sample_pi_puppers))))

same_subject_snp_plowers = []
same_subject_snp_puppers = []
for sample_pair_idx in xrange(0,len(same_subject_snp_idxs[0])):
    
    i = same_subject_snp_idxs[0][sample_pair_idx]
    j = same_subject_snp_idxs[1][sample_pair_idx]
    
    plower,pupper = stats_utils.calculate_poisson_rate_interval(snp_difference_matrix[i,j], snp_opportunity_matrix[i,j])
    
    same_subject_snp_plowers.append(plower)
    same_subject_snp_puppers.append(pupper)
    
    
# clip lower bounds 
same_subject_snp_plowers = numpy.clip(same_subject_snp_plowers,1e-09,1e09)
# Sort both lists by ascending lower bound on SNP changes, then gene changes
same_subject_snp_plowers, same_subject_snp_puppers = (numpy.array(x) for x in zip(*sorted(zip(same_subject_snp_plowers, same_subject_snp_puppers))))


diff_subject_snp_plowers = []
diff_subject_snp_puppers = []
for sample_pair_idx in xrange(0,len(diff_subject_snp_idxs[0])):
    
    i = diff_subject_snp_idxs[0][sample_pair_idx]
    j = diff_subject_snp_idxs[1][sample_pair_idx]
    
    plower,pupper = stats_utils.calculate_poisson_rate_interval(snp_difference_matrix[i,j], snp_opportunity_matrix[i,j])
    
    diff_subject_snp_plowers.append(plower)
    diff_subject_snp_puppers.append(pupper)
    
    
# clip lower bounds 
diff_subject_snp_plowers = numpy.clip(diff_subject_snp_plowers,1e-09,1e09)
# Sort both lists by ascending lower bound on SNP changes, then gene changes
diff_subject_snp_plowers, diff_subject_snp_puppers = (numpy.array(x) for x in zip(*sorted(zip(diff_subject_snp_plowers, diff_subject_snp_puppers))))

# Done calculating... now plot figure!

# Set up figure
fig = plt.figure(figsize=(5, 5))

# Set up grids to hold figure panels
outer_grid = gridspec.GridSpec(1, 2, width_ratios=[1,1], wspace=0.1)

###################
#
# SNP Panel
#
###################

snp_axis = plt.Subplot(fig, outer_grid[0])
fig.add_subplot(snp_axis)
fig.suptitle(species_name)

snp_axis.set_ylabel('Sample pairs')
snp_axis.set_xlabel('Substitution rate')
snp_axis.set_xlim([1e-07,9e-02])

snp_axis.semilogx([1e-09,1e-09],[1,1],'g-',label='Within host')
snp_axis.semilogx([1e-09,1e-09],[1,1],'r-',label='Between host')

###################
#
# Pi 
#
###################

pi_axis = plt.Subplot(fig, outer_grid[1])
fig.add_subplot(pi_axis)

pi_axis.set_xlabel('Within sample diversity')
pi_axis.set_xlim([1e-07,9e-02])

y = 0
for snp_plower, snp_pupper in zip(same_subject_snp_plowers, same_subject_snp_puppers):

    y-=1
    
    snp_axis.semilogx([snp_plower,snp_pupper],[y,y],'g.-',linewidth=0.25,markersize=1.5)
        
y-=1
snp_axis.semilogx([1e-09,1e09],[y,y,],'-',linewidth=0.25,color='0.7')


if len(diff_subject_snp_plowers)<=300:
    for snp_plower, snp_pupper in zip(diff_subject_snp_plowers, diff_subject_snp_puppers)[0:100]:

        y-=1
    
        snp_axis.semilogx([snp_plower,snp_pupper],[y,y],'r-',linewidth=0.35)
        
# If more than 300, do three sets of 100
else:

    for snp_plower, snp_pupper in zip(diff_subject_snp_plowers, diff_subject_snp_puppers)[0:100]:

        y-=1
    
        snp_axis.semilogx([snp_plower,snp_pupper],[y,y],'r-',linewidth=0.35)
        

    y-=1
    snp_axis.semilogx([1e-09,1e09],[y,y,],'-',linewidth=0.25,color='0.7')
    
    idxs = randint(0,len(diff_subject_snp_plowers),100)
    idxs.sort()

    for idx in idxs:
    
        snp_plower = diff_subject_snp_plowers[idx]
        snp_pupper = diff_subject_snp_puppers[idx]
        
        y-=1
    
        snp_axis.semilogx([snp_plower,snp_pupper],[y,y],'r-',linewidth=0.35)
        
    # Now do last hundred
    y-=1
    snp_axis.semilogx([1e-09,1e09],[y,y,],'-',linewidth=0.25, color='0.7')
    
    for snp_plower, snp_pupper in zip(diff_subject_snp_plowers, diff_subject_snp_puppers)[-100:]:

        y-=1
    
        snp_axis.semilogx([snp_plower,snp_pupper],[y,y],'r-',linewidth=0.35)
        

snp_axis.set_ylim([y-1,0])
snp_axis.set_yticks([])


# Now redo everything for pi
y = 0
for pi_plower, pi_pupper in zip(same_sample_pi_plowers, same_sample_pi_puppers):

    y-=1
    
    pi_axis.semilogx([pi_plower,pi_pupper],[y,y],'b.-',linewidth=0.25,markersize=1.5)
        


pi_axis.set_ylim([y-1,0])    

pi_axis.set_yticks([])

snp_axis.legend(loc='lower left',frameon=False)

fig.savefig('%s/%s_ordered_pi_within_between.pdf' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight')
fig.savefig('%s/%s_ordered_pi_within_between.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight',dpi=300)

pylab.figure(2,figsize=(3.42,2.5))
pylab.title(species_name)
pylab.loglog([clipped_pis.min()/2,clipped_pis.max()*2],[clipped_pis.min()/2,clipped_pis.max()*2],'k-')
pylab.loglog([clipped_pis.min()/2,clipped_pis.max()*2],[1e-03, 1e-03],'k-')
pylab.loglog(lower_pis,upper_pis,'g.',alpha=0.5)
pylab.xlabel('min(pi1,pi2)')
pylab.ylabel('max(pi1,pi2)')
pylab.xlim([clipped_pis.min()/2,clipped_pis.max()*2])
pylab.ylim([clipped_pis.min()/2,clipped_pis.max()*2])
pylab.savefig('%s/%s_temporal_within_sample_diversity.pdf' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight')
pylab.savefig('%s/%s_temporal_within_sample_diversity.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight',dpi=300)
    
