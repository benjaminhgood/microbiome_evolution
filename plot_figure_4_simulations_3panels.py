import matplotlib  
matplotlib.use('Agg') 
import config
import parse_midas_data
import sample_utils
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

  
# manually include the species for which I have data:
good_species_list=['Bacteroides_fragilis_54507', 'Bacteroides_vulgatus_57955', 'Parabacteroides_distasonis_56985']

passed_species = []

for species_name in good_species_list:

    # Load precomputed LD
    ld_map = calculate_linkage_disequilibria.load_ld_map(species_name)
     
    if len(ld_map)>0:
     
        distances, rsquared_numerators, rsquared_denominators, ns, intergene_distances, intergene_rsquared_numerators, intergene_rsquared_denominators, intergene_ns, control_rsquared_numerators, control_rsquared_denominators, control_n, pi = ld_map[('largest_clade','4D')]
    
        if True:
            passed_species.append(species_name)
        else:
            sys.stderr.write("%s intergene LD too high: %g (%g)\n" % (species_name, control_rsquared, rsquareds[0])) 

#passed_species = species_phylogeny_utils.sort_phylogenetically(passed_species)
num_passed_species = len(passed_species)

####################################################
#
# Set up Figure (3 panels, arranged in 3x1 grid)
#
####################################################

pylab.figure(3,figsize=(7,1.25))

ld_fig = pylab.gcf()
# make three panels panels
ld_grid  = gridspec.GridSpec(1,3,width_ratios=[1,1,1],wspace=0.25)

ld_axis_1 = plt.Subplot(ld_fig, ld_grid[0])
ld_fig.add_subplot(ld_axis_1)

ld_axis_1.set_title(figure_utils.get_pretty_species_name('Bacteroides_fragilis_54507'), fontsize=5)
ld_axis_1.set_ylabel('Linkage disequilibrium, $\sigma^2_d$')
ld_axis_1.set_xlabel('Distance between SNPs, $\ell$')
ld_axis_1.spines['top'].set_visible(False)
ld_axis_1.spines['right'].set_visible(False)
ld_axis_1.spines['bottom'].set_zorder(22)
ld_axis_1.get_xaxis().tick_bottom()
ld_axis_1.get_yaxis().tick_left()  
ld_axis_1.set_xlim([2,1e04])
ld_axis_1.set_ylim([1e-02,1])
ld_axis_1.text(6e03,3.5e-03,'Genome-\nwide', horizontalalignment='center',fontsize='5')



ld_axis_2 = plt.Subplot(ld_fig, ld_grid[1])
ld_fig.add_subplot(ld_axis_2)

ld_axis_2.set_title(figure_utils.get_pretty_species_name('Bacteroides_vulgatus_57955'), fontsize=5)
ld_axis_2.set_xlabel('Distance between SNPs, $\ell$')
ld_axis_2.spines['top'].set_visible(False)
ld_axis_2.spines['right'].set_visible(False)
ld_axis_2.spines['bottom'].set_zorder(22)
ld_axis_2.get_xaxis().tick_bottom()
#ld_axis_2.get_yaxis().tick_left()  
ld_axis_2.set_xlim([2,1e04])
ld_axis_2.set_ylim([1e-02,1])
ld_axis_2.text(6e03,3.5e-03,'Genome-\nwide', horizontalalignment='center',fontsize='5')

ld_axis_3 = plt.Subplot(ld_fig, ld_grid[2])
ld_fig.add_subplot(ld_axis_3)

ld_axis_3.set_title(figure_utils.get_pretty_species_name('Parabacteroides_distasonis_56985'), fontsize=5)
ld_axis_3.set_xlabel('Distance between SNPs, $\ell$')
ld_axis_3.spines['top'].set_visible(False)
ld_axis_3.spines['right'].set_visible(False)
ld_axis_3.spines['bottom'].set_zorder(22)
ld_axis_3.get_xaxis().tick_bottom()
#ld_axis_3.get_yaxis().tick_left()  
ld_axis_3.set_xlim([2,1e04])
ld_axis_3.set_ylim([1e-02,1])
ld_axis_3.text(6e03,3.5e-03,'Genome-\nwide', horizontalalignment='center',fontsize='5')




##################
species_idxs=xrange(0,num_passed_species)
axes = [ld_axis_1, ld_axis_2, ld_axis_3]

for species_idx, axis in zip(species_idxs,axes):
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
    
    #if species_name.startswith('Bacteroides_fragilis'):
    # put a dummy if because I don't feel like re-indenting rest of code:
    if 1==1:

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
        
        #axis.loglog(all_distances, all_rsquareds,'-',color='0.7',label='All samples')
        #axis.fill_between(numpy.array([3.3e03,1e04]),numpy.array([1e-02,1e-02]), numpy.array([1,1]),color='w',zorder=20)

        #axis.loglog([all_distances[-1],6e03], [all_rsquareds[-1], all_control_rsquared],':',color='0.7',zorder=21)
        #axis.loglog([6e03], [all_control_rsquared],'o',color='0.7',markersize=3,markeredgewidth=0,zorder=21)
        
        axis.fill_between(distances[good_distances],lower_rsquareds[good_distances], upper_rsquareds[good_distances], linewidth=0, color='b',alpha=0.5)
        axis.loglog(distances, rsquareds,'b-',label='All isolates')
        axis.loglog(early_distances, early_rsquareds,'bo',markersize=2,markeredgewidth=0,alpha=0.5)
        
        axis.loglog([distances[-1],6e03], [rsquareds[-1], control_rsquared],'b:',zorder=21)
        axis.loglog([6e03], [control_rsquared],'bo',markersize=3,markeredgewidth=0,zorder=21)
        
        axis.loglog(theory_ls, theory_rsquareds/theory_rsquareds[0]*3e-01,'k-',linewidth=0.3,zorder=0,label='Neutral')
        
        
        leg = axis.legend(loc='lower left',frameon=False)
        leg._legend_box.align = "left"
        
        line, = axis.loglog([distances[-1],distances[-1]],[1e-02,1],'k:')
        line.set_dashes((0.5,1))
    
        axis.xaxis.get_major_ticks()[-2].label1.set_visible(False)
        axis.xaxis.get_major_ticks()[-2].tick1line.set_visible(False)
        
        for tick_idx in xrange(1,7):
        
            axis.xaxis.get_minor_ticks()[-tick_idx].tick1line.set_visible(False)
            axis.xaxis.get_minor_ticks()[-tick_idx].tick2line.set_visible(False)
    
        
    

sys.stderr.write("Saving figure...\t")

ld_fig.savefig('%s/figure_4.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')

sys.stderr.write("Done!\n")


