import matplotlib  
matplotlib.use('Agg') 
import parse_midas_data
import sample_utils
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
from numpy.random import randint,binomial


import config
import sfs_utils
import stats_utils
import figure_utils

fontsize = 6
mpl.rcParams['font.size'] = fontsize
mpl.rcParams['lines.linewidth'] = 1.0
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'


min_coverage = config.min_median_coverage

species_name = "Bacteroides_vulgatus_57955"
sample_1 = '700023337' # complicated polyploid
sample_2 = '700101638' #'700171066' #'700096380'  # simple polyploid
sample_3 = '700116148'  # "complicated" haploid
sample_4 = '700023267'  # simple haploid

sample_labels = ['A','B','C','D']

# Make label for samples in left
left_samples = [sample_1, sample_2, sample_3, sample_4]
    

haploid_color = '#08519c'
light_haploid_color = '#6699CC'
diploid_color = '#de2d26' #'#fb6a4a' #
transition_color = '#756bb1'

#could also do this one...
#species_name = Bacteroides_uniformis_57318

# Load subject and sample metadata
sys.stderr.write("Loading sample metadata...\n")
subject_sample_map = sample_utils.parse_subject_sample_map()
sample_order_map = sample_utils.parse_sample_order_map()
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

pylab.figure(1,figsize=(5,2.7))
fig = pylab.gcf()
# make three panels
outer_grid  = gridspec.GridSpec(1,2, width_ratios=[2.7, 1.1], wspace=0.05)

##############################################################################
#
# Panel (a). Rank ordered within-host polymorhpism rate for focal species
#
##############################################################################

focal_grid = gridspec.GridSpecFromSubplotSpec(1, 2, width_ratios=[1,1.6],
                subplot_spec=outer_grid[0], wspace=0.5)

polymorphism_grid = gridspec.GridSpecFromSubplotSpec(2, 1, height_ratios=[0.01,1],
                subplot_spec=focal_grid[1], hspace=0.05)





#depth_axis = plt.Subplot(fig, polymorphism_grid[0])
#fig.add_subplot(depth_axis)

#depth_axis.set_ylabel("Depth")

#depth_axis.set_ylim([1e01,3e03])
#depth_axis.set_xticks([])

#depth_axis.set_title(species_name,fontsize=fontsize)


polymorphism_axis = plt.Subplot(fig, polymorphism_grid[1])
fig.add_subplot(polymorphism_axis)

polymorphism_axis.set_xlabel("Ranked samples (n=%d)" % len(desired_samples))
polymorphism_axis.set_ylabel("Within-sample polymorphism ($0.2 \leq f \leq 0.8$)")

polymorphism_axis.set_ylim([1e-06,2e-01])
polymorphism_axis.set_xticks([])

polymorphism_axis.set_title(figure_utils.get_pretty_species_name(species_name),fontsize=fontsize)

#polymorphism_axis.set_title(species_name,fontsize=fontsize)

##############################################################################
#
# Panel (b). Three example (folded) SFSs for focal species
#
##############################################################################

sfs_grid = gridspec.GridSpecFromSubplotSpec(4, 1, height_ratios=[1,1,1,1],
                subplot_spec=focal_grid[0], hspace=0.25)
                
sfs_axis_1 = plt.Subplot(fig, sfs_grid[0])
fig.add_subplot(sfs_axis_1)

sfs_axis_1.set_title('Sample %s ($\\overline{D}=%d$)' % (sample_labels[0], sample_coverage_map[sample_1]),fontsize=5,y=0.9)
sfs_axis_1.set_xticks([10*i for i in xrange(0,11)])
sfs_axis_1.set_xticklabels([])
sfs_axis_1.set_xlim([50,100])
sfs_axis_1.set_yticks([])
sfs_axis_1.xaxis.tick_bottom()


sfs_axis_2 = plt.Subplot(fig, sfs_grid[1])
fig.add_subplot(sfs_axis_2)

sfs_axis_2.set_title('Sample %s ($\\overline{D}=%d$)' % (sample_labels[1], sample_coverage_map[sample_2]),fontsize=5,y=0.9)
sfs_axis_2.set_xticks([10*i for i in xrange(0,11)])
sfs_axis_2.set_xticklabels([])
sfs_axis_2.set_xlim([50,100])
sfs_axis_2.set_yticks([])
sfs_axis_2.xaxis.tick_bottom()

#sfs_axis_2.set_ylabel('Fraction of 4D sites')

sfs_axis_3 = plt.Subplot(fig, sfs_grid[2])
fig.add_subplot(sfs_axis_3)

sfs_axis_3.set_title('Sample %s ($\\overline{D}=%d$)' % (sample_labels[2], sample_coverage_map[sample_3]),fontsize=5,y=0.9)
sfs_axis_3.set_xticks([10*i for i in xrange(0,11)])
sfs_axis_3.set_xticklabels([])
sfs_axis_3.set_xlim([50,100])
sfs_axis_3.set_yticks([])
sfs_axis_3.set_ylabel('                  Fraction of synonymous sites')
sfs_axis_3.xaxis.tick_bottom()


sfs_axis_4 = plt.Subplot(fig, sfs_grid[3])
fig.add_subplot(sfs_axis_4)

sfs_axis_4.set_title('Sample %s ($\\overline{D}=%d$)' % (sample_labels[3], sample_coverage_map[sample_4]),fontsize=5,y=0.9)
sfs_axis_4.set_xlabel('Major allele freq (%)')

sfs_axis_4.set_xticks([10*i for i in xrange(0,11)])
sfs_axis_4.set_xlim([50,100])
sfs_axis_4.set_yticks([])
sfs_axis_4.xaxis.tick_bottom()

##############################################################################
#
# Panel (c). Number of "haploids" for most prevalent species
#
##############################################################################

haploid_axis = plt.Subplot(fig, outer_grid[1])
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

####################################################
#
# Set up Suppplemental Fig (gene copynum distribution for 4 example samples)
#
####################################################
# This figure spreads them all out

pylab.figure(3,figsize=(7,1.25))
copynum_fig = pylab.gcf()
# make three panels panels
copynum_grid  = gridspec.GridSpec(1,4,width_ratios=[1,1,1,1],wspace=0.05)

copynum_axis_1 = plt.Subplot(copynum_fig, copynum_grid[0])
copynum_fig.add_subplot(copynum_axis_1)

copynum_axis_1.set_title('Sample A (D=%d)' % sample_coverage_map[sample_1],fontsize=5)
copynum_axis_1.set_ylabel('Fraction of genes')
copynum_axis_1.set_xlabel('Estimated gene copynum')

#copynum_axis_1.set_xlim([0,4])
copynum_axis_1.set_yticks([])

copynum_axis_2 = plt.Subplot(copynum_fig, copynum_grid[1])
copynum_fig.add_subplot(copynum_axis_2)

copynum_axis_2.set_title('Sample B (D=%d)' % sample_coverage_map[sample_2],fontsize=5)
copynum_axis_2.set_xlabel('Estimated gene copynum')
#copynum_axis_2.set_xlim([0,4])
copynum_axis_2.set_yticks([])

copynum_axis_3 = plt.Subplot(copynum_fig, copynum_grid[2])
copynum_fig.add_subplot(copynum_axis_3)

copynum_axis_3.set_title('Sample C (D=%d)' % sample_coverage_map[sample_3],fontsize=5)
copynum_axis_3.set_xlabel('Estimated gene copynum')

#copynum_axis_3.set_xlim([0,4])
copynum_axis_3.set_yticks([])

copynum_axis_4 = plt.Subplot(copynum_fig, copynum_grid[3])
copynum_fig.add_subplot(copynum_axis_4)

copynum_axis_4.set_title('Sample D (D=%d)' % sample_coverage_map[sample_4],fontsize=5)
copynum_axis_4.set_xlabel('Estimated gene copynum')

#copynum_axis_4.set_xlim([0,4])
copynum_axis_4.set_yticks([])


pylab.figure(4,figsize=(3.42,3.42))
avg_distance_fig = pylab.gcf()
# make three panels
outer_grid  = gridspec.GridSpec(2, 1, height_ratios=[0.3,1], hspace=0.1)

##############################################################################
#
# Panel (a). Rank ordered within-host polymorhpism rate for focal species
#
##############################################################################

between_axis = plt.Subplot(avg_distance_fig, outer_grid[0])
avg_distance_fig.add_subplot(between_axis)

between_axis.set_ylabel("Avg genetic\ndistance")

between_axis.set_xticks([])

between_axis.set_title(figure_utils.get_pretty_species_name(species_name),fontsize=fontsize)

supplemental_polymorphism_axis = plt.Subplot(avg_distance_fig, outer_grid[1])
avg_distance_fig.add_subplot(supplemental_polymorphism_axis)

supplemental_polymorphism_axis.set_xlabel("Ranked samples (n=%d)" % len(desired_samples))
supplemental_polymorphism_axis.set_ylabel("Within-sample polymorphism")

supplemental_polymorphism_axis.set_ylim([1e-06,2e-01])
supplemental_polymorphism_axis.set_xticks([])


####################################################
#
# Set up Suppplemental Fig (diversity haploid correlation)
#
####################################################

pylab.figure(5,figsize=(3.42,2.5))
correlation_fig = pylab.gcf()
# make three panels panels
outer_grid  = gridspec.GridSpec(1,1)

correlation_axis = plt.Subplot(correlation_fig, outer_grid[0])
correlation_fig.add_subplot(correlation_axis)
correlation_axis.set_xlabel('Fraction non-QP samples')
correlation_axis.set_ylabel('Avg within-host polymorphism')

####################################################
#
# Set up Suppplemental Fig (haploid_distribution)
#
####################################################

pylab.figure(6,figsize=(5,1.7))
haploid_distribution_fig = pylab.gcf()
# make three panels panels
outer_grid  = gridspec.GridSpec(1,2,wspace=0.3,width_ratios=[1,1])

haploid_distribution_axis = plt.Subplot(haploid_distribution_fig, outer_grid[1])
haploid_distribution_fig.add_subplot(haploid_distribution_axis)
haploid_distribution_axis.set_xlabel('# species per sample')
haploid_distribution_axis.set_ylabel('# QP species per sample')

haploid_cdf_axis = plt.Subplot(haploid_distribution_fig, outer_grid[0])
haploid_distribution_fig.add_subplot(haploid_cdf_axis)
haploid_cdf_axis.set_xlabel('Fraction QP species per sample, $p$')
haploid_cdf_axis.set_ylabel('Fraction samples $\geq p$')

###################################
#
# Calculate within polymorphism rates
#
###################################

sample_names = []

between_rates = []

within_rates = []
within_rate_lowers = []
within_rate_uppers = []

median_depths = []
depth_lowers = []
depth_uppers = []

for sample in desired_samples:
    within_sites, between_sites, total_sites = sfs_utils.calculate_polymorphism_rates_from_sfs_map(sfs_map[sample])

    within_rate = within_sites*1.0/total_sites
    between_rate = between_sites*1.0/total_sites    
    between_rates.append(between_rate)
    within_rates.append(within_rate)
    
    # Calculate 95% confidence intervals
    within_rate_lower, within_rate_upper = stats_utils.calculate_poisson_rate_interval(within_sites, total_sites,alpha=0.05)
    within_rate_lowers.append(within_rate_lower)
    within_rate_uppers.append(within_rate_upper)
   
    depths, counts = sfs_utils.calculate_depth_distribution_from_sfs_map(sfs_map[sample])
    dlower, dupper = stats_utils.calculate_IQR_from_distribution(depths, counts)
    dmedian = stats_utils.calculate_median_from_distribution(depths,counts)
    
    depth_lowers.append(dlower)
    depth_uppers.append(dupper)
    median_depths.append(dmedian)
    sample_names.append(sample)
    
# Sort them all in descending order of within-host diversity    
within_rates, within_rate_lowers, within_rate_uppers, between_rates, median_depths, depth_lowers, depth_uppers, sample_names = (numpy.array(x) for x in zip(*sorted(zip(within_rates, within_rate_lowers, within_rate_uppers, between_rates, median_depths,depth_lowers, depth_uppers, sample_names),reverse=True)))

within_rate_lowers = numpy.clip(within_rate_lowers, 1e-09,1)
    
for rank_idx in xrange(0,len(within_rates)):
    
    polymorphism_axis.semilogy([rank_idx,rank_idx], [within_rate_lowers[rank_idx],within_rate_uppers[rank_idx]],'-',color='k',linewidth=0.25)
    supplemental_polymorphism_axis.semilogy([rank_idx,rank_idx], [within_rate_lowers[rank_idx],within_rate_uppers[rank_idx]],'-',color='k',linewidth=0.25)
    
    between_axis.semilogy([rank_idx], [between_rates[rank_idx]],'.',color='k',markersize=2.5,alpha=0.5,markeredgewidth=0)
    #depth_axis.semilogy([rank_idx], [median_depths[rank_idx]],'.',color=haploid_color,markersize=2.5,alpha=0.5,markeredgewidth=0)
    
    if sample_names[rank_idx] in left_samples:
         
        sample_idx = left_samples.index(sample_names[rank_idx])
        
        print "Printing number for", sample_names[rank_idx], sample_idx+1

        polymorphism_axis.text(rank_idx, within_rate_uppers[rank_idx]*1.3,'%s' % sample_labels[sample_idx],horizontalalignment='left', verticalalignment='bottom',fontsize=6)
        
###################################
#
# Plot example SFSs
#
###################################

# Sample 1
fs,pfs = sfs_utils.calculate_binned_sfs_from_sfs_map(sfs_map[sample_1],folding='major')
df = fs[1]-fs[0]

within_sites, between_sites, total_sites = sfs_utils.calculate_polymorphism_rates_from_sfs_map(sfs_map[sample_1])

    

between_line = between_sites*1.0/total_sites/((fs>0.2)*(fs<0.5)).sum()
#pmax = numpy.max([pfs[(fs>0.1)*(fs<0.95)].max(), between_line])
pmax = between_line

within_rate = within_sites*1.0/total_sites
print "Sample 1: within =", within_rate, "avg-distance =", between_sites*1.0/total_sites

sfs_axis_1.fill_between([80,100],[0,0],[1,1],color='0.8')

#sfs_axis_1.fill_between([20,100],[0,0],[1,1],color='0.8')
sfs_axis_1.bar((fs-df/2)*100,pfs,width=df, edgecolor=light_haploid_color, color=light_haploid_color)
#line, = sfs_axis_1.plot([20,80], [between_line,between_line], 'k-',linewidth=0.35)
#line.set_dashes((1.5,1))
sfs_axis_1.set_ylim([0,pmax*3])

# Sample 2
fs,pfs = sfs_utils.calculate_binned_sfs_from_sfs_map(sfs_map[sample_2],folding='major')
df = fs[1]-fs[0]
within_sites, between_sites, total_sites = sfs_utils.calculate_polymorphism_rates_from_sfs_map(sfs_map[sample_2])
between_line = between_sites*1.0/total_sites/((fs>0.2)*(fs<0.5)).sum()

within_rate = within_sites*1.0/total_sites
print "Sample 2: within =", within_rate, "avg-distance =", between_sites*1.0/total_sites


#pmax = numpy.max([pfs[(fs>0.1)*(fs<0.95)].max(), between_line])
pmax = between_line
sfs_axis_2.fill_between([80,100],[0,0],[1,1],color='0.8')

sfs_axis_2.bar((fs-df/2)*100,pfs,width=df,edgecolor=light_haploid_color, color=light_haploid_color)
#line, = sfs_axis_2.plot([20,80], [between_line,between_line], 'k-',linewidth=0.35)
#line.set_dashes((1.5,1))

sfs_axis_2.set_ylim([0,pmax*3])

# Sample 3
fs,pfs = sfs_utils.calculate_binned_sfs_from_sfs_map(sfs_map[sample_3],folding='major')
df = fs[1]-fs[0]
within_sites, between_sites, total_sites = sfs_utils.calculate_polymorphism_rates_from_sfs_map(sfs_map[sample_3])
between_line = between_sites*1.0/total_sites/((fs>0.2)*(fs<0.5)).sum()
#pmax = numpy.max([pfs[(fs>0.1)*(fs<0.95)].max(), between_line])
pmax = between_line

within_rate = within_sites*1.0/total_sites
print "Sample 3: within =", within_rate, "avg-distance =", between_sites*1.0/total_sites
sfs_axis_3.fill_between([80,100],[0,0],[1,1],color='0.8')

sfs_axis_3.bar((fs-df/2)*100,pfs,width=df,edgecolor=haploid_color,color=haploid_color)
#line, = sfs_axis_3.plot([20,80], [between_line,between_line], 'k-',linewidth=0.35)
#line.set_dashes((1.5,1))

sfs_axis_3.set_ylim([0,pmax*3])

# Sample 4
fs,pfs = sfs_utils.calculate_binned_sfs_from_sfs_map(sfs_map[sample_4],folding='major')
df = fs[1]-fs[0]
within_sites, between_sites, total_sites = sfs_utils.calculate_polymorphism_rates_from_sfs_map(sfs_map[sample_4])
between_line = between_sites*1.0/total_sites/((fs>0.2)*(fs<0.5)).sum()
#pmax = numpy.max([pfs[(fs>0.1)*(fs<0.95)].max(), between_line])
pmax = between_line

within_rate = within_sites*1.0/total_sites
print "Sample 4: within =", within_rate, "avg-distance =", between_sites*1.0/total_sites
sfs_axis_4.fill_between([80,100],[0,0],[1,1],color='0.8')

sfs_axis_4.bar((fs-df/2)*100,pfs,width=df,edgecolor=haploid_color,color=haploid_color)
#line, = sfs_axis_4.plot([20,80], [between_line,between_line], 'k-',linewidth=0.35)
#line.set_dashes((1.5,1))

sfs_axis_4.set_ylim([0,pmax*3])

##########################
#
# Calculate&plot copynum distributions for example samples
#
##########################

desired_samples = numpy.array([sample_1,sample_2,sample_3, sample_4])
axes = [copynum_axis_1, copynum_axis_2, copynum_axis_3, copynum_axis_4]

# Load gene coverage information for species_name
sys.stderr.write("Loading pangenome data for %s...\n" % species_name)
gene_samples, gene_names, gene_presence_matrix, gene_depth_matrix, marker_coverages, gene_reads_matrix = parse_midas_data.parse_pangenome_data(species_name,allowed_samples=desired_samples)
sys.stderr.write("Done! Loaded gene info for %d samples\n" % len(gene_samples))

idxs = numpy.array([numpy.nonzero(gene_samples==sample)[0][0] for sample in desired_samples])

# Calculate copynums
gene_copynum_matrix = gene_depth_matrix*1.0/(marker_coverages+(marker_coverages==0))

for idx,axis in zip(idxs, axes):

    gene_copynums = numpy.array(gene_copynum_matrix[:,idx],copy=True)
    #gene_copynums = numpy.clip(gene_copynums, 1e-02,1e02)
    gene_copynums.sort()
    
    gene_copynums = gene_copynums[gene_copynums>1e-02]
    
    
    bins = numpy.logspace(-2,1,100)
    
    #xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(gene_copynums)
    
    #print xs[0], ns[0]
    
    #axis.loglog([0.05,0.05],[1.0/ns[0], 1],'k-',linewidth=0.25)
    #axis.loglog([0.5,0.5],[1.0/ns[0], 1],'k-',linewidth=0.25)
    #axis.loglog([2,2],[1.0/ns[0], 1],'k-',linewidth=0.25)
    #axis.set_ylim([1.0/ns[0],1])
    #axis.step(xs,ns*1.0/ns[0],color='b')
    
    heights,bins,patches = axis.hist(gene_copynums,bins=bins,zorder=1,color=haploid_color,edgecolor=haploid_color)
    
    
    sfs_axis_4.fill_between([0,20],[0,0],[1,1],color='0.8')
    
    ymax = heights.max()*1.1
    axis.fill_between([0.01, 0.05],[0,0],[ymax,ymax],color='0.8',zorder=0)
    axis.fill_between([0.6, 1.2],[0,0],[ymax,ymax],color='0.8',zorder=0)
    axis.semilogx([2],[-1],'k.')
    axis.set_ylim([0,ymax])
    #axis.set_xlim([1e-02,4e00])   
    axis.set_xlim([5e-02,5])
sys.stderr.write("Saving figure...\t")
copynum_fig.savefig('%s/supplemental_copynum_distributions.pdf' % parse_midas_data.analysis_directory, bbox_inches='tight')
sys.stderr.write("Done!\n")


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
avg_within_rates = []
fraction_haploids = []
fraction_polyploids = []

sample_haploid_species = {}
sample_highcoverage_species = {}

for species_name in good_species_list:
    
    desired_samples = diversity_utils.calculate_highcoverage_samples(species_name)
    haploid_samples = diversity_utils.calculate_haploid_samples(species_name)
    n_haploids = len(haploid_samples)
    
    for sample in desired_samples:
        if sample not in sample_highcoverage_species:
            sample_highcoverage_species[sample] = []
            sample_haploid_species[sample] = []
            
        sample_highcoverage_species[sample].append(species_name)
    
    sample_ploidy_map = {sample: 'polyploid' for sample in desired_samples}
    for sample in haploid_samples:
        sample_ploidy_map[sample] = 'haploid'
        sample_haploid_species[sample].append(species_name)
        
    if len(desired_samples)<10:
        continue
    
    # Load SNP information for species_name
    sys.stderr.write("Loading SFSs for %s...\t" % species_name)
    samples, sfs_map = parse_midas_data.parse_within_sample_sfs(species_name,     allowed_variant_types=set(['4D'])) 
    sys.stderr.write("Done!\n")

    
    avg_within_sites = 0.0
    avg_total_sites = 0.0
    for sample in desired_samples:
        within_sites, between_sites, total_sites = sfs_utils.calculate_polymorphism_rates_from_sfs_map(sfs_map[sample])
        avg_within_sites += within_sites
        avg_total_sites += total_sites
        
    avg_within_sites /= len(desired_samples)
    avg_total_sites /= len(desired_samples)
        
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
    avg_within_rates.append( avg_within_sites/avg_total_sites )
    fraction_haploids.append( len(haploid_samples)*1.0/len(desired_samples) )
    fraction_polyploids.append( (len(desired_samples)-len(haploid_samples)+1)*1.0/(len(desired_samples)+1))
 
# Sort by num haploids    
num_haploid_samples, num_samples, ploidy_changes, species_names, avg_within_rates, fraction_haploids, fraction_polyploids = zip(*sorted(zip(num_haploid_samples, num_samples, ploidy_changes, species_names, avg_within_rates, fraction_haploids, fraction_polyploids),reverse=True))

num_haploid_samples = numpy.array(num_haploid_samples)
num_samples = numpy.array(num_samples)
species_names = numpy.array(species_names)
avg_within_rates = numpy.array(avg_within_rates)
fraction_polyploids = numpy.array(fraction_polyploids)
fraction_haploids = numpy.array(fraction_haploids)
total_haploids = num_haploid_samples.sum()

##############
#
# Plot correlation between fraction haploids and avg within polymorphism
#
#############
correlation_axis.loglog(fraction_polyploids[num_samples>=10], avg_within_rates[num_samples>=10], 'o', color=haploid_color,alpha=0.5)

#correlation_axis.loglog([1e-02,1],[1e-04,1e-02],'-',color=haploid_color)  


# Sort by num samples    
num_samples, num_haploid_samples, species_names, ploidy_changes = (numpy.array(x) for x in zip(*sorted(zip(num_samples, num_haploid_samples, species_names, ploidy_changes),reverse=True)))
    
haploid_haploid_samples = []  
haploid_polyploid_samples = []
polyploid_haploid_samples = []
polyploid_polyploid_samples = []

#############
#
# Plot distribution of percentage haploids per sample
#
#############
sample_highcoverage_counts = []
sample_haploid_counts = []
for sample in sample_highcoverage_species.keys():
    sample_highcoverage_counts.append( len(sample_highcoverage_species[sample]))
    sample_haploid_counts.append( len(sample_haploid_species[sample]))
    
sample_highcoverage_counts = numpy.array(sample_highcoverage_counts)
sample_haploid_counts = numpy.array(sample_haploid_counts)

sample_haploid_fractions = sample_haploid_counts*1.0/(sample_highcoverage_counts+(sample_highcoverage_counts==0))

haploid_distribution_axis.plot(sample_highcoverage_counts+normal(0,0.1,size=sample_highcoverage_counts.shape), sample_haploid_counts+normal(0,0.1,size=sample_highcoverage_counts.shape),'.',color=haploid_color,alpha=0.5,markersize=2)    

pavg = sample_haploid_counts.sum()*1.0/sample_highcoverage_counts.sum()

num_bootstraps=100
bootstrapped_haploid_countss = []
for bootstrap_idx in xrange(0,num_bootstraps):
    
    bootstrapped_haploid_countss.append( binomial(sample_highcoverage_counts, pavg) )
    
pooled_bootstrapped_haploid_fractions = []
for bootstrap_idx in xrange(0,num_bootstraps):
    pooled_bootstrapped_haploid_fractions.extend( bootstrapped_haploid_countss[bootstrap_idx][sample_highcoverage_counts>=1]*1.0/sample_highcoverage_counts[sample_highcoverage_counts>=1] )
pooled_bootstrapped_haploid_fractions = numpy.array( pooled_bootstrapped_haploid_fractions )

xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(pooled_bootstrapped_haploid_fractions )
haploid_cdf_axis.step(xs,ns*1.0/ns[0],'-',color='0.7',label='Null')

xs, ns = stats_utils.calculate_unnormalized_survival_from_vector(sample_haploid_fractions[sample_highcoverage_counts>=1])
haploid_cdf_axis.step(xs,ns*1.0/ns[0],'-',color=haploid_color,label='Obs')
haploid_cdf_axis.set_xlim([0,1])

haploid_cdf_axis.legend(loc='upper right',frameon=False,numpoints=1)


#########
#
# Haploid distribution
#
for species_idx in xrange(0,len(num_haploid_samples)):
        
    print species_names[species_idx], num_haploid_samples[species_idx], num_samples[species_idx]
    print ploidy_changes[species_idx]
    
    haploid_haploid_samples.append( ploidy_changes[species_idx]['haploid->haploid'] )
    haploid_polyploid_samples.append( ploidy_changes[species_idx]['haploid->polyploid'] )
    polyploid_haploid_samples.append( ploidy_changes[species_idx]['polyploid->haploid'] )
    polyploid_polyploid_samples.append( ploidy_changes[species_idx]['polyploid->polyploid'] )
    
    if species_idx==19:
        print "Top 20 ^"

haploid_haploid_samples = numpy.array(haploid_haploid_samples)
haploid_polyploid_samples = numpy.array(haploid_polyploid_samples)
polyploid_haploid_samples = numpy.array(polyploid_haploid_samples)
polyploid_polyploid_samples = numpy.array(polyploid_polyploid_samples)

num_haploid_samples = num_haploid_samples[:35]
num_samples = num_samples[:35]
haploid_species_names = species_names[:35]

# get better haploid species names
pretty_haploid_species_names = []
for species_name in haploid_species_names:

    base_name = figure_utils.get_pretty_species_name(species_name)
    pretty_name = base_name
    if pretty_name in pretty_haploid_species_names:
        idx = 1
        while pretty_name in pretty_haploid_species_names:
            idx += 1
            pretty_name = base_name + (" %d" % (idx))
    
    pretty_haploid_species_names.append(pretty_name)
    
haploid_species_names = pretty_haploid_species_names

total_temporal_samples = haploid_haploid_samples+haploid_polyploid_samples+polyploid_haploid_samples+polyploid_polyploid_samples

temporal_species_names = species_names[total_temporal_samples>4]
haploid_haploid_samples = haploid_haploid_samples[total_temporal_samples>4]
haploid_polyploid_samples = haploid_polyploid_samples[total_temporal_samples>4]
polyploid_haploid_samples = polyploid_haploid_samples[total_temporal_samples>4]
polyploid_polyploid_samples = polyploid_polyploid_samples[total_temporal_samples>4]
total_temporal_samples = total_temporal_samples[total_temporal_samples>4]

# sort by total temporal temporal samples
# Sort by num samples    
total_temporal_samples, temporal_species_names, haploid_haploid_samples, haploid_polyploid_samples, polyploid_haploid_samples, polyploid_polyploid_samples = (numpy.array(x) for x in zip(*sorted(zip(total_temporal_samples, temporal_species_names, haploid_haploid_samples, haploid_polyploid_samples, polyploid_haploid_samples, polyploid_polyploid_samples),reverse=True)))
 
ys = 0-numpy.arange(0,len(num_haploid_samples))
width=0.7

haploid_axis.barh(ys, num_haploid_samples,color=haploid_color,linewidth=0,label='QP',zorder=1)
haploid_axis.barh(ys, num_samples,color=light_haploid_color,linewidth=0,label='non-QP',zorder=0)
haploid_axis.set_xlim([0,800])
#haploid_axis.set_xticks([0,100,200,300,400,500,600])
haploid_axis.set_xticks([0,200,400,600,800])

haploid_axis.yaxis.tick_right()
haploid_axis.xaxis.tick_bottom()

haploid_axis.set_yticks(ys+0.5)
haploid_axis.set_yticklabels(haploid_species_names,fontsize=4)
haploid_axis.set_ylim([-1*len(num_haploid_samples)+0.5,1.5])

haploid_axis.tick_params(axis='y', direction='out',length=3,pad=1)

haploid_axis.legend(loc='lower right',frameon=False)

ys = 0-numpy.arange(0,len(temporal_species_names))

haploid_haploid_color = haploid_color
haploid_haploid_alpha=1.0

haploid_polyploid_color = haploid_color
haploid_polyploid_alpha = 0.5

polyploid_haploid_color = diploid_color
polyploid_haploid_alpha = 0.5

polyploid_polyploid_color = diploid_color
polyploid_polyploid_alpha = 1.0

# Plot bars
temporal_haploid_axis.barh(ys, haploid_haploid_samples,linewidth=0,label='QP->QP',zorder=4, color=haploid_haploid_color, alpha=haploid_haploid_alpha)

temporal_haploid_axis.barh(ys, haploid_polyploid_samples, left=haploid_haploid_samples, linewidth=0,label='QP->non',zorder=3, color=haploid_polyploid_color, alpha=haploid_polyploid_alpha)

temporal_haploid_axis.barh(ys, polyploid_haploid_samples, left=haploid_haploid_samples+haploid_polyploid_samples, linewidth=0,label='non->QP',zorder=2, color=polyploid_haploid_color, alpha=polyploid_haploid_alpha)

temporal_haploid_axis.barh(ys, polyploid_polyploid_samples, left=haploid_haploid_samples+haploid_polyploid_samples+polyploid_haploid_samples, linewidth=0,label='non->non',zorder=1, color=polyploid_polyploid_color, alpha=polyploid_polyploid_alpha)

temporal_haploid_axis.yaxis.tick_right()
temporal_haploid_axis.xaxis.tick_bottom()

temporal_haploid_axis.set_yticks(ys+0.5)
temporal_haploid_axis.set_yticklabels(temporal_species_names,fontsize=4)
temporal_haploid_axis.set_ylim([-1*len(temporal_species_names)+1,1])

temporal_haploid_axis.tick_params(axis='y', direction='out',length=3,pad=1)


temporal_haploid_axis.legend(loc='lower right',frameon=False)

sys.stderr.write("%d haploid samples across species\n" % total_haploids)


####
#
# Save figures 
#
####
sys.stderr.write("Saving figure 1...\t")
fig.savefig('%s/figure_1.pdf' % parse_midas_data.analysis_directory, bbox_inches='tight')
sys.stderr.write("Done!\n")

sys.stderr.write("Saving supplemental figures...\t")
#fig2.savefig('%s/supplemental_temporal_haploid.pdf' % parse_midas_data.analysis_directory, bbox_inches='tight')
#sys.stderr.write("Done!\n")

avg_distance_fig.savefig('%s/supplemental_avg_distance.pdf' % parse_midas_data.analysis_directory, bbox_inches='tight')

correlation_fig.savefig('%s/supplemental_haploid_correlation.pdf' % parse_midas_data.analysis_directory, bbox_inches='tight')

haploid_distribution_fig.savefig('%s/supplemental_haploid_distribution_fig.pdf' % parse_midas_data.analysis_directory, bbox_inches='tight')


sys.stderr.write("Done!\n")


