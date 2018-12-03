# LD. Decay of LD within genes across species?

import matplotlib  
matplotlib.use('Agg') 
import sample_utils
import config
import parse_midas_data
import os.path
import pylab
import sys
import numpy
from math import exp

import species_phylogeny_utils
import diversity_utils
import gene_diversity_utils
import calculate_temporal_changes
import calculate_substitution_rates
import calculate_linkage_disequilibria
import calculate_snv_distances
import stats_utils
import sfs_utils
import figure_utils

from scipy.optimize import least_squares, newton, brentq
      
import matplotlib.colors as colors
import matplotlib.cm as cmx
from math import log10,ceil
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy.random import randint, multinomial
import matplotlib.colors as mcolors

from math import log

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
parser.add_argument("--memoize", help="Loads stuff from disk", action="store_true")
parser.add_argument("--chunk-size", type=int, help="max number of records to load", default=1000000000)

args = parser.parse_args()

debug = args.debug
chunk_size = args.chunk_size
memoize = args.memoize

################################################################################

min_sample_size = config.between_host_min_sample_size # 46 gives at least 1000 
low_divergence_threshold = config.between_low_divergence_threshold
allowed_variant_types = set(['4D'])

#focal_speciess = ['Bacteroides_vulgatus_57955', 'Roseburia_inulinivorans_61943']
#focal_speciess = ['Bacteroides_vulgatus_57955', 'Faecalibacterium_prausnitzii_62201']
focal_speciess = ['Bacteroides_vulgatus_57955', 'Akkermansia_muciniphila_55290']

focal_colors = ['b','g']

#supplemental_focal_species = ['Bacteroides_fragilis_54507', 'Alistipes_putredinis_61533', 'Eubacterium_rectale_56927']
supplemental_focal_species = ['Bacteroides_fragilis_54507', 'Parabacteroides_distasonis_56985', 'Alistipes_shahii_62199'] # 'Ruminococcus_bromii_62047']  
good_species_list = parse_midas_data.parse_good_species_list()
if debug:
    good_species_list = good_species_list[0:2]

sys.stderr.write("Loading sample metadata...\n")
subject_sample_map = sample_utils.parse_subject_sample_map()
sample_continent_map = sample_utils.parse_sample_continent_map()
sys.stderr.write("Done!\n")

####################################################
#
# Set up Figure (2 panels, arranged in 2x1 grid)
#
####################################################

pylab.figure(1,figsize=(5,6))
fig = pylab.gcf()
# make three panels panels

outer_grid = gridspec.GridSpec(3,1,height_ratios=[0.9,0.6,1],hspace=0.6)

upper_grid = gridspec.GridSpecFromSubplotSpec(1, 3, width_ratios=[0.15,1,0.15], wspace=0.3, subplot_spec=outer_grid[0])

middle_grid = gridspec.GridSpecFromSubplotSpec(1, 5, width_ratios=[0.1,1,0.5,1,0.1], wspace=0.1, subplot_spec=outer_grid[1])
                
species_grid = gridspec.GridSpecFromSubplotSpec(2, 1, height_ratios=[1.0/1.7,0.7/1.7],
                subplot_spec=outer_grid[2], hspace=0)
                
# Inconsistency axis
inconsistency_axis = plt.Subplot(fig, upper_grid[1])
fig.add_subplot(inconsistency_axis)

inconsistency_axis.spines['top'].set_visible(False)
inconsistency_axis.spines['right'].set_visible(False)
inconsistency_axis.spines['bottom'].set_zorder(22)

inconsistency_axis.get_xaxis().tick_bottom()
inconsistency_axis.get_yaxis().tick_left()

inconsistency_axis.set_xlabel('Maximum divergence age of SNV, $d_B^*$')
inconsistency_axis.set_ylabel('Phylogenetic inconsistency between\nSNVs & core-genome divergence')
inconsistency_axis.set_xlim([2e-05,2e-02])
inconsistency_axis.set_ylim([0,1.05])


passed_species = []
sample_sizes = []
for species_name in good_species_list:

    sys.stderr.write("Loading haploid samples...\n")
    # Only plot samples above a certain depth threshold that are "haploids"
    snp_samples = diversity_utils.calculate_haploid_samples(species_name, debug=debug)
    
    if len(snp_samples) < min_sample_size:
        sys.stderr.write("Not enough haploid samples!\n")
        continue
        
    sys.stderr.write("Calculating unique samples...\n")
    # Only consider one sample per person
    snp_samples = snp_samples[sample_utils.calculate_unique_samples(subject_sample_map, sample_list=snp_samples)]


    if len(snp_samples) < min_sample_size:
        sys.stderr.write("Not enough unique samples!\n")
        continue

    # Load inconsistency data
    sys.stderr.write("(core genes only...)\n")
    snv_distance_map = calculate_snv_distances.load_snv_distance_map(species_name)
    
    ds = numpy.logspace(log10(3e-05),log10(3e-02),50) # 15 points are plotted
    total_snps = numpy.zeros_like(ds)
    inconsistent_snps = numpy.zeros_like(ds)
      
    for location_tuple in snv_distance_map:
        var_type, derived_allele_counts, ancestral_allele_counts, between_d, within_d1, within_d2 = snv_distance_map[location_tuple]
        
        if var_type in allowed_variant_types:
            
            # In this calculation, ds are the list of threshold dB*'s 
            
            within_d = min([within_d1, within_d2])
            good_idxs = (between_d<=ds)
            inconsistent_idxs = good_idxs*(within_d>=2*ds)
            
            total_snps[good_idxs] += 1
            inconsistent_snps[inconsistent_idxs] += 1
    
    fraction_inconsistent = inconsistent_snps*1.0/(total_snps+(total_snps==0))
    
    if species_name in focal_speciess:
        focal_species_idx = focal_speciess.index(species_name)
        color = focal_colors[focal_species_idx]
        linewidth=1
        zorder=2
        alpha=1
    else:
        color = 'r'
        alpha =0.3
        linewidth=0.5
        zorder = 1

    inconsistency_axis.semilogx(ds[total_snps>0], fraction_inconsistent[total_snps>0],'-',color=color,linewidth=linewidth,zorder=zorder,alpha=alpha)

    # Load precomputed LD
    ld_map = calculate_linkage_disequilibria.load_ld_map(species_name)
     
    if len(ld_map)>0:
     
        distances, rsquared_numerators, rsquared_denominators, ns, intergene_distances, intergene_rsquared_numerators, intergene_rsquared_denominators, intergene_ns, control_rsquared_numerators, control_rsquared_denominators, control_n, pi = ld_map[('largest_clade','4D')]
    
        if True:
            passed_species.append(species_name)
            sample_sizes.append(len(snp_samples))
        else:
            sys.stderr.write("%s intergene LD too high: %g (%g)\n" % (species_name, control_rsquared, rsquareds[0])) 

sample_sizes = numpy.array(sample_sizes)

passed_species = species_phylogeny_utils.sort_phylogenetically(passed_species, first_entry=focal_speciess[0], second_sorting_attribute=(-1*sample_sizes))

#passed_species = species_phylogeny_utils.sort_phylogenetically(passed_species)
num_passed_species = len(passed_species)


inconsistency_axis.plot([1],[-1],'-',color=focal_colors[0],linewidth=1, alpha=1,label=figure_utils.get_pretty_species_name(focal_speciess[0], include_number=False))
inconsistency_axis.plot([1],[-1],'-',color=focal_colors[1],linewidth=1, alpha=1,label=figure_utils.get_pretty_species_name(focal_speciess[1], include_number=False))
inconsistency_axis.plot([1],[-1],'r-',linewidth=0.5, alpha=0.3,label='Other species')
inconsistency_axis.legend(loc='lower left',frameon=False,fontsize=5,numpoints=1,handlelength=1)


####################################################
#
# Set up Figure (1 panels, arranged in nx1 grid)
#
####################################################

pylab.figure(3,figsize=(3.42,2))
scaled_fig = pylab.gcf()
# make three panels panels
scaled_outer_grid  = gridspec.GridSpec(1,1)

scaled_axis = plt.Subplot(scaled_fig, scaled_outer_grid[0])
scaled_fig.add_subplot(scaled_axis)
scaled_axis.set_xlim([1e-02,1e02])
scaled_axis.set_ylim([2e-02,2])

focal_example_axes = []
for focal_species_idx in xrange(0,len(focal_speciess)):

    focal_species = focal_speciess[focal_species_idx]

    ####
    #
    # Continue with the main fig
    #
    ####
    focal_example_axis = plt.Subplot(fig, middle_grid[1+2*focal_species_idx])
    fig.add_subplot(focal_example_axis)
    focal_example_axis.set_ylabel('Linkage disequilibrium, $\sigma^2_d$')
    focal_example_axis.set_xlabel('Distance between SNVs, $\ell$')

    focal_example_axis.spines['top'].set_visible(False)
    focal_example_axis.spines['right'].set_visible(False)
    focal_example_axis.spines['bottom'].set_zorder(22)

    focal_example_axis.get_xaxis().tick_bottom()
    focal_example_axis.get_yaxis().tick_left()

    focal_example_axis.set_xlim([2,1e04])
    focal_example_axis.set_ylim([1e-02,1])

    focal_example_axis.text(6e03,4e-03,'Genome-\nwide',     horizontalalignment='center',fontsize='5')
    
    focal_example_axes.append(focal_example_axis)

species_axis = plt.Subplot(fig, species_grid[0])
fig.add_subplot(species_axis)
species_axis.set_xlim([-1.5,num_passed_species-0.5])

species_axis.spines['top'].set_visible(False)
species_axis.spines['right'].set_visible(False)

species_axis.set_ylabel('Linkage Disequilibrium, $\sigma^2_d$')

species_axis.get_xaxis().tick_bottom()
species_axis.get_yaxis().tick_left()

xticks = numpy.arange(0,num_passed_species)
#xticklabels = ["%s (%d)" % (species_names[i],sample_sizes[i]) for i in xrange(0,len(sample_sizes))]
#xticklabels = ["%s" % (passed_species[i]) for i in xrange(0, num_passed_species)]
xticklabels = ["%s" % (figure_utils.get_pretty_species_name(passed_species[i])) for i in xrange(0, num_passed_species)]

species_axis.set_xticks(xticks)
species_axis.set_xticklabels(xticklabels, rotation='vertical',fontsize=4)
 
#species_axis.set_ylim([2e-02,2])

####################################################
#
# Set up Figure (2 panels, arranged in 2x1 grid)
#
####################################################

pylab.figure(4,figsize=(4,1))
rbymu_fig = pylab.gcf()
# make three panels panels

rbymu_outer_grid  = gridspec.GridSpec(1,1)

rbymu_axis = plt.Subplot(rbymu_fig, rbymu_outer_grid[0])
rbymu_fig.add_subplot(rbymu_axis)
rbymu_axis.set_xlim([-1.5,num_passed_species-0.5])
rbymu_axis.set_ylabel('Effective $r/\mu$')
rbymu_axis.spines['top'].set_visible(False)
rbymu_axis.spines['right'].set_visible(False)
rbymu_axis.get_xaxis().tick_bottom()
rbymu_axis.get_yaxis().tick_left()

xticks = numpy.arange(0,num_passed_species)
#xticklabels = ["%s (%d)" % (species_names[i],sample_sizes[i]) for i in xrange(0,len(sample_sizes))]
xticklabels = ["%s" % (figure_utils.get_pretty_species_name(passed_species[i])) for i in xrange(0, num_passed_species)]

rbymu_axis.set_xticks(xticks)
rbymu_axis.set_xticklabels(xticklabels, rotation='vertical',fontsize=4)


####################################################
#
# Set up Figure (n panels, arranged in nx1 grid)
#
####################################################

pylab.figure(2,figsize=(3.42,2*num_passed_species))
fit_fig = pylab.gcf()
# make three panels panels
outer_grid  = gridspec.GridSpec(num_passed_species,1,height_ratios=([1]*num_passed_species),hspace=0.5)

fit_axes = []
for species_idx in xrange(0,num_passed_species):
    species_name = passed_species[species_idx]
    fit_axis = plt.Subplot(fit_fig, outer_grid[species_idx])
    fit_fig.add_subplot(fit_axis)

    fit_axis.set_ylabel(species_name)
    #fit_axis.set_ylim([-35,35])
    fit_axis.set_xlim([1e-02,1e02])
    #fit_axis.loglog([1e,1e04],[1,1],'k:')
    fit_axis.set_xticklabels([])
    fit_axes.append(fit_axis)
 
#####
#
# Supplemental figure (3 example species LD decay)
#
#####
####################################################
#
# Set up Figure (2 panels, arranged in 2x1 grid)
#
####################################################

pylab.figure(5,figsize=(7,1.7))
fig5 = pylab.gcf()
# make three panels panels

outer_grid  = gridspec.GridSpec(1,3,width_ratios=[1,1,1],wspace=0.1)

example_axes = []
for example_idx in xrange(0,3):                

    example_axis = plt.Subplot(fig5, outer_grid[example_idx])
    fig5.add_subplot(example_axis)

    if example_idx==0:
        example_axis.set_ylabel('Linkage disequilibrium, $\sigma^2_d$')
    
    
    example_axis.spines['top'].set_visible(False)
    example_axis.spines['right'].set_visible(False)
    example_axis.spines['bottom'].set_zorder(22)

    example_axis.get_xaxis().tick_bottom()
    example_axis.get_yaxis().tick_left()

    example_axis.set_xlim([2,1e04])
    example_axis.set_ylim([1e-02,1])
    
    if example_idx==1:
        example_axis.set_xlabel('Distance between SNVs, $\ell$')
        example_axis.text(6e03,5.3e-03,'Genome-\nwide', horizontalalignment='center',fontsize='5')
        
        
    #example_axis.set_title(supplemental_focal_species[example_idx],fontsize=5)
    
    example_axes.append(example_axis)

 
def neutral_rsquared(NRs):
    return (10.0+2*NRs)/(22.0+26*NRs+4*NRs*NRs)
    
def normalized_neutral_rsquared(NRs):
    return neutral_rsquared(NRs)/neutral_rsquared(0)
    
def calculate_effective_NR(rsquared_ratio):
    return brentq(lambda x: normalized_neutral_rsquared(x)-rsquared_ratio, 0, 1e09)
     
 
 
    
for species_idx in xrange(0,num_passed_species):
    species_name = passed_species[species_idx]
    # Load precomputed LD
    ld_map = calculate_linkage_disequilibria.load_ld_map(species_name)
     
    if len(ld_map)==0:
        continue
    
    all_distances, all_rsquared_numerators, all_rsquared_denominators, all_ns, all_intergene_distances, all_intergene_rsquared_numerators, all_intergene_rsquared_denominators, all_intergene_ns, all_control_rsquared_numerator, all_control_rsquared_denominator, all_control_n, all_pi = ld_map[('all','4D')]
    all_control_rsquared = all_control_rsquared_numerator/all_control_rsquared_denominator
             
    distances, rsquared_numerators, rsquared_denominators, ns, intergene_distances, intergene_rsquared_numerators, intergene_rsquared_denominators, intergene_ns, control_rsquared_numerator, control_rsquared_denominator, control_n, pi = ld_map[('largest_clade','4D')]
    control_rsquared = control_rsquared_numerator/control_rsquared_denominator
    
    # smooth this stuff:
    smoothed_distances = distances
    window_width = 10**(0.1)
    
    dmins = smoothed_distances/(window_width**0.5)
    dmaxs = smoothed_distances*(window_width**0.5)
    
    smoothed_rsquared_numerators = []
    smoothed_rsquared_denominators = []
    smoothed_counts = []
    
    all_smoothed_rsquared_numerators = []
    all_smoothed_rsquared_denominators = []
    all_smoothed_counts = []
    
    for dmin,dmax in zip(dmins,dmaxs):
        binned_numerators = rsquared_numerators[(distances>=dmin)*(distances<=dmax)]
        binned_denominators = rsquared_denominators[(distances>=dmin)*(distances<=dmax)]
        binned_counts = ns[(distances>=dmin)*(distances<=dmax)]
        smoothed_rsquared_numerators.append( (binned_numerators*binned_counts).sum()/binned_counts.sum() )
        smoothed_rsquared_denominators.append( (binned_denominators*binned_counts).sum()/binned_counts.sum() )
        smoothed_counts.append( binned_counts.sum() )
        
        binned_numerators = all_rsquared_numerators[(distances>=dmin)*(distances<=dmax)]
        binned_denominators = all_rsquared_denominators[(distances>=dmin)*(distances<=dmax)]
        binned_counts = all_ns[(distances>=dmin)*(distances<=dmax)]
        all_smoothed_rsquared_numerators.append( (binned_numerators*binned_counts).sum()/binned_counts.sum() )
        all_smoothed_rsquared_denominators.append( (binned_denominators*binned_counts).sum()/binned_counts.sum() )
        all_smoothed_counts.append( binned_counts.sum() )
        
        
    smoothed_rsquared_numerators = numpy.array( smoothed_rsquared_numerators )
    smoothed_rsquared_denominators = numpy.array( smoothed_rsquared_denominators )
    smoothed_counts = numpy.array( smoothed_counts )
    
    all_smoothed_rsquared_numerators = numpy.array( all_smoothed_rsquared_numerators )
    all_smoothed_rsquared_denominators = numpy.array( all_smoothed_rsquared_denominators )
    all_smoothed_counts = numpy.array( all_smoothed_counts )
    
    early_distances = distances[distances<101]
    early_rsquareds = rsquared_numerators[distances<101]*1.0/rsquared_denominators[distances<101]
    early_ns = ns[distances<101]
    
    early_distances = early_distances[early_ns>0.5]    
    early_rsquareds = early_rsquareds[early_ns>0.5]
    early_ns = early_ns[early_ns>0.5]
    
    distances = smoothed_distances
    rsquareds = smoothed_rsquared_numerators/(smoothed_rsquared_denominators)
    ns = smoothed_counts
    
    distances = distances[ns>0]
    rsquareds = rsquareds[ns>0]
    ns = ns[ns>0]
    
    all_distances = smoothed_distances
    #all_distances = dmins
    all_rsquareds = all_smoothed_rsquared_numerators/(all_smoothed_rsquared_denominators)
    all_ns = all_smoothed_counts
    
    all_distances = all_distances[all_ns>0]
    all_rsquareds = all_rsquareds[all_ns>0]
    all_ns = all_ns[all_ns>0]
    
    if (species_name in focal_speciess) or (species_name in supplemental_focal_species):
        
        if (species_name in focal_speciess):
            focal_species_idx = focal_speciess.index(species_name)
            example_axis=focal_example_axes[focal_species_idx]
            example_idx = -1
            color=focal_colors[focal_species_idx]
        else:
            example_idx = supplemental_focal_species.index(species_name)
            example_axis = example_axes[example_idx]
            color=focal_colors[0]
            
        num_bootstraps = 10
        
        bootstrapped_sigmasquareds = [] # will eventually be a matrix where first index is window_idx and second index is bootstrap index (sorted from lowest to highest)
        # Estimate bootstrap intervals for focal species only
        for dmin,dmax in zip(dmins,dmaxs):
            binned_numerators = rsquared_numerators[(distances>=dmin)*(distances<=dmax)]
            binned_denominators = rsquared_denominators[(distances>=dmin)*(distances<=dmax)]
            binned_counts = ns[(distances>=dmin)*(distances<=dmax)]
            
            total_pairs = binned_counts.sum()
            
            upper_rsquareds = []
            lower_rsquareds = []
            
            if total_pairs>0:
            
                if True: #len(binned_counts)>1:
                    #print total_pairs
                    #print binned_counts
                    ps = binned_counts*1.0/total_pairs
            
                    window_bootstrapped_countss = multinomial(total_pairs,ps,size=num_bootstraps)
            
                    #print window_bootstrapped_countss.shape
                    window_bootstrapped_numerators = (window_bootstrapped_countss*binned_numerators[None,:]).sum(axis=1)*1.0/total_pairs
                    window_bootstrapped_denominators = (window_bootstrapped_countss*binned_denominators[None,:]).sum(axis=1)*1.0/total_pairs
            
                    window_bootstrapped_sigmasquareds = window_bootstrapped_numerators/window_bootstrapped_denominators
            
                    #print window_bootstrapped_sigmasquareds.shape
                    window_bootstrapped_sigmasquareds.sort()
            
                    bootstrapped_sigmasquareds.append(window_bootstrapped_sigmasquareds)
                    
                    #print total_pairs
                    
                else:
                    bootstrapped_sigmasquareds.append([binned_numerators/binned_denominators]*num_bootstraps)
                    
            else:
                
                bootstrapped_sigmasquareds.append([-1]*num_bootstraps)
                
        
        upper_rsquareds = numpy.array([bootstrapped_sigmasquareds[window_idx][long(num_bootstraps*0.95)] for window_idx in xrange(0,len(bootstrapped_sigmasquareds))])
        lower_rsquareds = numpy.array([bootstrapped_sigmasquareds[window_idx][long(num_bootstraps*0.05)] for window_idx in xrange(0,len(bootstrapped_sigmasquareds))])
        
        #print upper_rsquareds-lower_rsquareds
        
        good_distances = (upper_rsquareds>=-0.5)*(lower_rsquareds>=-0.5)
        
        theory_ls = numpy.logspace(0,log10(distances[-1]),100)
        theory_NRs = theory_ls/200.0
        theory_rsquareds = (10+2*theory_NRs)/(22+26*theory_NRs+4*theory_NRs*theory_NRs)
        
        example_axis.loglog(all_distances, all_rsquareds,'-',color='0.7',label='All samples')
        example_axis.fill_between(numpy.array([3.3e03,1e04]),numpy.array([1e-02,1e-02]), numpy.array([1,1]),color='w',zorder=20)

        example_axis.loglog([all_distances[-1],6e03], [all_rsquareds[-1], all_control_rsquared],':',color='0.7',zorder=21)
        example_axis.loglog([6e03], [all_control_rsquared],'o',color='0.7',markersize=3,markeredgewidth=0,zorder=21)
        example_axis.fill_between(distances[good_distances],lower_rsquareds[good_distances], upper_rsquareds[good_distances], linewidth=0, color=color,alpha=0.5)
        example_axis.loglog(distances, rsquareds,'-',color=color,label='Largest clade')
        example_axis.loglog(early_distances, early_rsquareds,'o',color=color,markersize=2,markeredgewidth=0,alpha=0.5)
        
        example_axis.loglog([distances[-1],6e03], [rsquareds[-1], control_rsquared],':',color=color,zorder=21)
        example_axis.loglog([6e03], [control_rsquared],'o',color=color,markersize=3,markeredgewidth=0,zorder=21)
        example_axis.set_title( figure_utils.get_pretty_species_name(species_name),fontsize=6,y=0.95)
        
        example_axis.loglog(theory_ls, theory_rsquareds/theory_rsquareds[0]*3e-01,'k-',linewidth=0.3,zorder=0,label='Neutral')
        
        
        leg = example_axis.legend(loc='lower left',frameon=False, ) #title=figure_utils.get_pretty_species_name(species_name,include_number=False))
        leg._legend_box.align = "left"
        
        #example_axis.set_title(figure_utils.get_pretty_species_name(species_name, include_number=True),fontsize=5)
        
        line, = example_axis.loglog([distances[-1],distances[-1]],[1e-02,1],'k:')
        line.set_dashes((0.5,1))
    
        
        if example_idx>0:
            print "Setting ytick labels"
            example_axis.set_yticklabels([])
    
        example_axis.xaxis.get_major_ticks()[-2].label1.set_visible(False)
        example_axis.xaxis.get_major_ticks()[-2].tick1line.set_visible(False)
                
        for tick_idx in xrange(1,7):
        
            example_axis.xaxis.get_minor_ticks()[-tick_idx].tick1line.set_visible( False)
            example_axis.xaxis.get_minor_ticks()[-tick_idx].tick2line.set_visible( False)
    
   
    #fit_axes[species_idx].loglog(old_distances*old_pi, old_rsquareds,'.-',color='0.7')    
    fit_axes[species_idx].loglog(distances*pi, rsquareds,'k-')
    #fit_axes[species_idx].loglog(theory_distances*pi, theory_rsquareds,'r-')
    fit_axes[species_idx].loglog([1e-02,1e02],[control_rsquared, control_rsquared],'k:')
    
    
    idx_9 = numpy.fabs(distances-9).argmin()
    idx_90 = numpy.fabs(distances-100).argmin()
    idx_900 = numpy.fabs(distances-2000).argmin()
    
    
    scaled_axis.loglog(distances*pi, rsquareds/(rsquareds[idx_9]),'.-',markersize=2,alpha=0.5)
    
    #sys.stderr.write("%s, l* = %g, intergene LD = %g\n" % (species_name, exp(res.x[2]), control_rsquared))
    
    
    
    #epsilons = [0.5,0.25] #,2*control_rsquared/rsquareds[0]]
    #rbymus = []
    #for eps in epsilons:
    
    #    def g(u):
            #print u, model(res.x,u), (rsquareds[0]*eps)
    #        return (model(res.x, u)-(K*eps))
    
    #    ustar = brentq(lambda u: g(u), -5, log(1e09))  
        #ustar = newton(lambda u: g(u), log(100)))
    #    lstar = exp(ustar)
        
    #    rbymu = 2*(1-eps)/eps/pi/lstar
    #    print rbymu
    #    rbymus.append(rbymu)
    
    if species_name in focal_speciess:
        focal_species_idx = focal_speciess.index(species_name)
        color=focal_colors[focal_species_idx]
    else:
        color='r'
    
    species_axis.semilogy([species_idx,species_idx], [rsquareds[idx_900], rsquareds[idx_9]],'-', color=color)
    species_axis.semilogy([species_idx], [rsquareds[idx_9]],'_',markersize=3,color=color) #,markeredgewidth=0)
    species_axis.semilogy([species_idx],[rsquareds[idx_90]],'_',markersize=3, color=color) #,markeredgewidth=0)
    species_axis.semilogy([species_idx], [rsquareds[idx_900]],'_',markersize=3,color=color) #,markeredgewidth=0)
    line, = species_axis.semilogy([species_idx,species_idx], [control_rsquared, rsquareds[idx_900]],':',color=color)
    line.set_dashes((0.5,0.75))
    species_axis.semilogy([species_idx], [control_rsquared],'o',markersize=2,markeredgewidth=0,color=color)


    # now do rbymu estimate
    if rsquareds[-1] < rsquareds[idx_9]/2:
    
        NRstar = calculate_effective_NR(0.5)
        lstar = distances[(rsquareds/rsquareds[idx_9]<=0.5)][0]
        # Old version
        # lstar = distances[numpy.fabs(rsquareds/rsquareds[idx_9]-0.5).argmin()]
        rbymu_2 = NRstar/lstar/pi*2
        
        rbymu_axis.semilogy([species_idx], [rbymu_2],'_',markersize=3,color=color) 
        
        if rsquareds.min() < rsquareds[idx_9]/4:
            critical_fraction = 0.25
        else:
            critical_fraction = rsquareds.min()/rsquareds[idx_9]
            
        if True:
            NRstar = calculate_effective_NR(critical_fraction)
            # get first point where LD/LD(0)<0.5
            lstar = distances[(rsquareds/rsquareds[idx_9]<=critical_fraction)][0]
            # Old version
            # lstar = distances[numpy.fabs(rsquareds/rsquareds[idx_9]-0.25).argmin()]
            rbymu_4 = NRstar/lstar/pi*2
        
            rbymu_axis.semilogy([species_idx,species_idx], [rbymu_4, rbymu_2],'-',color=color)
            rbymu_axis.semilogy([species_idx], [rbymu_4],'_',markersize=3,color=color) 
        

sys.stderr.write("Saving figure...\t")
fit_fig.savefig('%s/supplemental_ld_fits.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')
scaled_fig.savefig('%s/supplemental_ld_shapes.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')
fig.savefig('%s/figure_4.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')
fig5.savefig('%s/supplemental_ld_decay.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')

rbymu_fig.savefig('%s/supplemental_rbymu.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')

sys.stderr.write("Done!\n")

     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
