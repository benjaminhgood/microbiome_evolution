import matplotlib  
matplotlib.use('Agg') 
import parse_midas_data
import parse_HMP_data
import pylab
import sys
import numpy
from numpy.random import normal
#from calculate_pi_matrix import calculate_self_pis
import diversity_utils
import stats_utils
import matplotlib.colors as colors
import matplotlib.cm as cmx
from math import log10,ceil
import matplotlib as mpl

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy.random import randint


import config
import sfs_utils


fontsize = 6
mpl.rcParams['font.size'] = fontsize
mpl.rcParams['lines.linewidth'] = 1.0
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'


min_coverage = config.min_median_coverage

species_name = "Bacteroides_vulgatus_57955"
sample_1 = '700023337'
sample_2 = '700116148'
sample_3 = '700023267'


#could also do this one...
#species_name = Bacteroides_uniformis_57318

# Load subject and sample metadata
sys.stderr.write("Loading sample metadata...\n")
subject_sample_map = parse_HMP_data.parse_subject_sample_map()
sample_country_map = parse_HMP_data.parse_sample_country_map()
sample_order_map = parse_HMP_data.parse_sample_order_map()
sys.stderr.write("Done!\n")
 

# Load SNP information for species_name
sys.stderr.write("Loading SFSs for %s...\t" % species_name)
samples, sfs_map = parse_midas_data.parse_within_sample_sfs(species_name, allowed_variant_types=set(['4D'])) 
sys.stderr.write("Done!\n")

# Load genomic coverage distributions
sample_coverage_histograms, samples = parse_midas_data.parse_coverage_distribution(species_name)
median_coverages = numpy.array([stats_utils.calculate_nonzero_median_from_histogram(sample_coverage_histogram) for sample_coverage_histogram in sample_coverage_histograms])
sample_coverage_map = {samples[i]: median_coverages[i] for i in xrange(0,len(samples))}
samples = numpy.array(samples)

median_coverages = numpy.array([sample_coverage_map[samples[i]] for i in xrange(0,len(samples))])

# Only plot samples above a certain depth threshold
desired_samples = samples[(median_coverages>=min_coverage)]
desired_median_coverages = numpy.array([sample_coverage_map[sample] for sample in desired_samples])

####################################################
#
# Set up Fig. 1 (3 panels, arranged horizontally)
#
####################################################
# This figure spreads them all out

pylab.figure(1,figsize=(5,2))
fig = pylab.gcf()
# make three panels panels
outer_grid  = gridspec.GridSpec(1,3, width_ratios=[1.5, 1, 1.3], wspace=0.25)

##############################################################################
#
# Panel (a). Rank ordered within-host polymorhpism rate for focal species
#
##############################################################################

polymorphism_axis = plt.Subplot(fig, outer_grid[0])
fig.add_subplot(polymorphism_axis)

polymorphism_axis.set_xlabel("Ranked samples (n=%d)" % len(desired_samples))
polymorphism_axis.set_ylabel("Within-sample polymorphism")

polymorphism_axis.set_ylim([1e-06,1e-01])
polymorphism_axis.set_xticks([])

polymorphism_axis.set_title(species_name,fontsize=fontsize)

##############################################################################
#
# Panel (b). Three example (folded) SFSs for focal species
#
##############################################################################

sfs_grid = gridspec.GridSpecFromSubplotSpec(3, 1, height_ratios=[1,1,1],
                subplot_spec=outer_grid[1], hspace=0.25)
                
sfs_axis_1 = plt.Subplot(fig, sfs_grid[0])
fig.add_subplot(sfs_axis_1)

sfs_axis_1.set_title('Sample 1 (D=%d)' % sample_coverage_map[sample_1],fontsize=5,y=0.9)
sfs_axis_1.set_xticks([10*i for i in xrange(0,11)])
sfs_axis_1.set_xticklabels([])
sfs_axis_1.set_xlim([0,50])
sfs_axis_1.set_yticks([])

sfs_axis_2 = plt.Subplot(fig, sfs_grid[1])
fig.add_subplot(sfs_axis_2)

sfs_axis_2.set_title('Sample 2 (D=%d)' % sample_coverage_map[sample_2],fontsize=5,y=0.9)
sfs_axis_2.set_xticks([10*i for i in xrange(0,11)])
sfs_axis_2.set_xticklabels([])
sfs_axis_2.set_xlim([0,50])
sfs_axis_2.set_yticks([])
sfs_axis_2.set_ylabel('Fraction of 4D sites')

sfs_axis_3 = plt.Subplot(fig, sfs_grid[2])
fig.add_subplot(sfs_axis_3)

sfs_axis_3.set_title('Sample 3 (D=%d)' % sample_coverage_map[sample_3],fontsize=5,y=0.9)
sfs_axis_3.set_xlabel('Minor allele freq (%)')

sfs_axis_3.set_xticks([10*i for i in xrange(0,11)])
sfs_axis_3.set_xlim([0,50])
sfs_axis_3.set_yticks([])

##############################################################################
#
# Panel (c). Number of "haploids" for most prevalent species
#
##############################################################################

haploid_axis = plt.Subplot(fig, outer_grid[2])
fig.add_subplot(haploid_axis)

haploid_axis.set_xlabel("Number of samples")

####################################################
#
# Set up Suppplemental Fig (temporal haploid classification)
#
####################################################
# This figure spreads them all out

pylab.figure(2,figsize=(2,3))
fig2 = pylab.gcf()
# make three panels panels
outer_grid2  = gridspec.GridSpec(1,1)

temporal_haploid_axis = plt.Subplot(fig2, outer_grid2[0])
fig2.add_subplot(temporal_haploid_axis)
temporal_haploid_axis.set_xlabel('Timepoint pairs')

###################################
#
# Calculate within polymorphism rates
#
###################################

within_rates = []
within_rate_lowers = []
within_rate_uppers = []

for sample in desired_samples:
    within_sites, between_sites, total_sites = sfs_utils.calculate_polymorphism_rates_from_sfs_map(sfs_map[sample])

    
    within_rate = within_sites*1.0/total_sites
    within_rate_lower, within_rate_upper = stats_utils.calculate_poisson_rate_interval(within_sites, total_sites)
    
    #print within_rate, within_rate_lower, within_rate_upper
    within_rates.append(within_rate)
    within_rate_lowers.append(within_rate_lower)
    within_rate_uppers.append(within_rate_upper)
    
within_rates, within_rate_lowers, within_rate_uppers = (numpy.array(x) for x in zip(*sorted(zip(within_rates, within_rate_lowers, within_rate_uppers),reverse=True)))

within_rate_lowers = numpy.clip(within_rate_lowers, 1e-09,1)
    
for rank_idx in xrange(0,len(within_rates)):
    
    polymorphism_axis.semilogy([rank_idx,rank_idx], [within_rate_lowers[rank_idx],within_rate_uppers[rank_idx]],'b-',linewidth=0.25)


###################################
#
# Plot example SFSs
#
###################################

# Sample 1
fs,pfs = sfs_utils.calculate_binned_sfs_from_sfs_map(sfs_map[sample_1])
df = fs[1]-fs[0]

within_sites, between_sites, total_sites = sfs_utils.calculate_polymorphism_rates_from_sfs_map(sfs_map[sample_1])
between_line = between_sites*1.0/total_sites/((fs>0.2)*(fs<0.5)).sum()
#pmax = numpy.max([pfs[(fs>0.1)*(fs<0.95)].max(), between_line])
pmax = between_line

print between_sites*1.0/total_sites

sfs_axis_1.fill_between([0,20],[0,0],[1,1],color='0.8')

#sfs_axis_1.fill_between([20,100],[0,0],[1,1],color='0.8')
sfs_axis_1.bar((fs-df/2)*100,pfs,width=df,edgecolor='b',color='b')
line, = sfs_axis_1.plot([20,100], [between_line,between_line], 'k-',linewidth=0.35)
line.set_dashes((1.5,1))
sfs_axis_1.set_ylim([0,pmax*3])

# Sample 2
fs,pfs = sfs_utils.calculate_binned_sfs_from_sfs_map(sfs_map[sample_2])
df = fs[1]-fs[0]
within_sites, between_sites, total_sites = sfs_utils.calculate_polymorphism_rates_from_sfs_map(sfs_map[sample_2])
between_line = between_sites*1.0/total_sites/((fs>0.2)*(fs<0.5)).sum()

print between_sites*1.0/total_sites


#pmax = numpy.max([pfs[(fs>0.1)*(fs<0.95)].max(), between_line])
pmax = between_line
sfs_axis_2.fill_between([0,20],[0,0],[1,1],color='0.8')

sfs_axis_2.bar((fs-df/2)*100,pfs,width=df,edgecolor='b',color='b')
line, = sfs_axis_2.plot([20,100], [between_line,between_line], 'k-',linewidth=0.35)
line.set_dashes((1.5,1))

sfs_axis_2.set_ylim([0,pmax*3])

# Sample 3
fs,pfs = sfs_utils.calculate_binned_sfs_from_sfs_map(sfs_map[sample_3])
df = fs[1]-fs[0]
within_sites, between_sites, total_sites = sfs_utils.calculate_polymorphism_rates_from_sfs_map(sfs_map[sample_3])
between_line = between_sites*1.0/total_sites/((fs>0.2)*(fs<0.5)).sum()
#pmax = numpy.max([pfs[(fs>0.1)*(fs<0.95)].max(), between_line])
pmax = between_line

print between_sites*1.0/total_sites
sfs_axis_3.fill_between([0,20],[0,0],[1,1],color='0.8')

sfs_axis_3.bar((fs-df/2)*100,pfs,width=df,edgecolor='b',color='b')
line, = sfs_axis_3.plot([20,100], [between_line,between_line], 'k-',linewidth=0.35)
line.set_dashes((1.5,1))

sfs_axis_3.set_ylim([0,pmax*3])


###################################
#
# Calculate number of haploids across species
#
###################################

good_species_list = parse_midas_data.parse_good_species_list()

species_names = []
num_samples = []
num_haploid_samples = []
ploidy_changes = []

for species_name in good_species_list:
    
    # Load genomic coverage distributions
    sample_coverage_histograms, samples = parse_midas_data.parse_coverage_distribution(species_name)
    median_coverages = numpy.array([stats_utils.calculate_nonzero_median_from_histogram(sample_coverage_histogram) for sample_coverage_histogram in sample_coverage_histograms])
    sample_coverage_map = {samples[i]: median_coverages[i] for i in xrange(0,len(samples))}
    samples = numpy.array(samples)

    median_coverages = numpy.array([sample_coverage_map[samples[i]] for i in xrange(0,len(samples))])

    # Only plot samples above a certain depth threshold
    desired_samples = samples[(median_coverages>=min_coverage)]
    desired_median_coverages = numpy.array([sample_coverage_map[sample] for sample in desired_samples])
    
    if len(desired_samples)<10:
        continue
    
    # Load SNP information for species_name
    sys.stderr.write("Loading SFSs for %s...\t" % species_name)
    samples, sfs_map = parse_midas_data.parse_within_sample_sfs(species_name,     allowed_variant_types=set(['4D'])) 
    sys.stderr.write("Done!\n")

    
    n_haploids = 0
    
    sample_ploidy_map = {}
    
    for sample in desired_samples:
        within_sites, between_sites, total_sites = sfs_utils.calculate_polymorphism_rates_from_sfs_map(sfs_map[sample])
    
        if within_sites <= config.threshold_within_between_fraction*between_sites:
            
            sample_ploidy_map[sample] = 'haploid'
            n_haploids += 1    
        else:
            sample_ploidy_map[sample] = 'polyploid'
            
    
    ploidy_change_map = {'haploid->haploid':0, 'haploid->polyploid':0,'polyploid->haploid':0, 'polyploid->polyploid':0}     

    # Calculate which pairs of idxs belong to the same sample, which to the same subject
    # and which to different subjects
    desired_same_sample_idxs, desired_same_subject_idxs, desired_diff_subject_idxs = parse_midas_data.calculate_ordered_subject_pairs(sample_order_map, desired_samples)

    for sample_pair_idx in xrange(0,len(desired_same_subject_idxs[0])):
    
        initial_sample = desired_samples[desired_same_subject_idxs[0][sample_pair_idx]]
        final_sample = desired_samples[desired_same_subject_idxs[1][sample_pair_idx]]
        
        ploidy_change_str = ("%s->%s" % (sample_ploidy_map[initial_sample], sample_ploidy_map[final_sample])) 
        
        ploidy_change_map[ploidy_change_str] += 1
        
        
        
        
    species_names.append(species_name)
    num_samples.append(len(desired_samples))
    num_haploid_samples.append(n_haploids)
    ploidy_changes.append(ploidy_change_map)
 
# Sort by num haploids    
num_haploid_samples, num_samples, species_names = (numpy.array(x) for x in zip(*sorted(zip(num_haploid_samples, num_samples, species_names),reverse=True)))

# Sort by num samples    
#num_samples, num_haploid_samples, species_names = (numpy.array(x) for x in zip(*sorted(zip(num_samples, num_haploid_samples, species_names),reverse=True)))
    
haploid_haploid_samples = []  
haploid_polyploid_samples = []
polyploid_polyploid_samples = []

    
for species_idx in xrange(0,len(num_haploid_samples)):
        
    print species_names[species_idx], num_haploid_samples[species_idx], num_samples[species_idx]
    print ploidy_changes[species_idx]
    
    haploid_haploid_samples.append( ploidy_changes[species_idx]['haploid->haploid'] )
    haploid_polyploid_samples.append( ploidy_changes[species_idx]['haploid->polyploid'] + ploidy_changes[species_idx]['polyploid->haploid'] )
    polyploid_polyploid_samples.append( ploidy_changes[species_idx]['polyploid->polyploid'] )
    
    if species_idx==19:
        print "Top 20 ^"

haploid_haploid_samples = numpy.array(haploid_haploid_samples)
haploid_polyploid_samples = numpy.array(haploid_polyploid_samples)
polyploid_polyploid_samples = numpy.array(polyploid_polyploid_samples)

num_haploid_samples = num_haploid_samples[:25]
num_samples = num_samples[:25]
species_names = species_names[:25]
haploid_haploid_samples = haploid_haploid_samples[:25]
haploid_polyploid_samples = haploid_polyploid_samples[:25]
polyploid_polyploid_samples = polyploid_polyploid_samples[:25]


ys = 0-numpy.arange(0,len(num_haploid_samples))
width=0.7

haploid_axis.barh(ys, num_samples,color='0.7',linewidth=0)
haploid_axis.barh(ys, num_haploid_samples,color='b',linewidth=0)
haploid_axis.set_xlim([0,300])

haploid_axis.yaxis.tick_right()
haploid_axis.xaxis.tick_bottom()

haploid_axis.set_yticks(ys+0.5)
haploid_axis.set_yticklabels(species_names,fontsize=4)
haploid_axis.set_ylim([-1*len(num_haploid_samples)+1,1])

temporal_haploid_axis.barh(ys, haploid_haploid_samples+haploid_polyploid_samples+polyploid_polyploid_samples,color='r',linewidth=0,label='non/non')
temporal_haploid_axis.barh(ys, haploid_haploid_samples+haploid_polyploid_samples,color='#8856a7',linewidth=0,label='CPS/non')
temporal_haploid_axis.barh(ys, haploid_haploid_samples,color='b',linewidth=0,label='CPS/CPS')

temporal_haploid_axis.yaxis.tick_right()
temporal_haploid_axis.xaxis.tick_bottom()

temporal_haploid_axis.set_yticks(ys+0.5)
temporal_haploid_axis.set_yticklabels(species_names,fontsize=4)
temporal_haploid_axis.set_ylim([-1*len(num_haploid_samples)+1,1])

temporal_haploid_axis.legend(loc='lower right',frameon=False)

sys.stderr.write("Saving figure...\t")
fig.savefig('%s/figure_1.pdf' % parse_midas_data.analysis_directory, bbox_inches='tight')
sys.stderr.write("Done!\n")

sys.stderr.write("Saving figure...\t")
fig2.savefig('%s/supplemental_temporal_haploid.pdf' % parse_midas_data.analysis_directory, bbox_inches='tight')
sys.stderr.write("Done!\n")
