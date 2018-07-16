# LD. Decay of LD within genes across species?

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

  
good_species_list = parse_midas_data.parse_good_species_list()
# manually include the species for which I have data:
good_species_list=['Bacteroides_fragilis_54507', 'Bacteroides_vulgatus_57955', 'Parabacteroides_distasonis_56985']
if debug:
    good_species_list = good_species_list[0:2]

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


####################################################
#
# Set up Figure (2 panels, arranged in 2x1 grid)
#
####################################################

pylab.figure(1,figsize=(7,1.7))
fig = pylab.gcf()
# make three panels panels

outer_grid  = gridspec.GridSpec(1,2,width_ratios=[3,4],wspace=0.2)

species_grid = gridspec.GridSpecFromSubplotSpec(2, 1, height_ratios=[1.0/1.7,0.7/1.7],
                subplot_spec=outer_grid[1], hspace=0)
                

example_axis = plt.Subplot(fig, outer_grid[0])
fig.add_subplot(example_axis)
example_axis.set_ylabel('Linkage disequilibrium, $\sigma^2_d$')
example_axis.set_xlabel('Distance between SNPs, $\ell$')

example_axis.spines['top'].set_visible(False)
example_axis.spines['right'].set_visible(False)
example_axis.spines['bottom'].set_zorder(22)

example_axis.get_xaxis().tick_bottom()
example_axis.get_yaxis().tick_left()

example_axis.set_xlim([2,1e04])
example_axis.set_ylim([1e-02,1])

example_axis.text(6e03,5.3e-03,'Genome-\nwide', horizontalalignment='center',fontsize='5')

species_axis = plt.Subplot(fig, species_grid[0])
fig.add_subplot(species_axis)
species_axis.set_xlim([-1.5,num_passed_species-0.5])
#species_axis.set_ylabel('LD, $\sigma^2_d$')

species_axis.spines['top'].set_visible(False)
species_axis.spines['right'].set_visible(False)


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
    
    if species_name.startswith('Bacteroides_fragilis'):
        
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
        
        
        
        example_axis.fill_between(distances[good_distances],lower_rsquareds[good_distances], upper_rsquareds[good_distances], linewidth=0, color='b',alpha=0.5)
        example_axis.loglog(distances, rsquareds,'b-',label='Largest clade')
        example_axis.loglog(early_distances, early_rsquareds,'bo',markersize=2,markeredgewidth=0,alpha=0.5)
        
        example_axis.loglog([distances[-1],6e03], [rsquareds[-1], control_rsquared],'b:',zorder=21)
        example_axis.loglog([6e03], [control_rsquared],'bo',markersize=3,markeredgewidth=0,zorder=21)
        #example_axis.set_title(figure_utils.get_pretty_species_name(species_name),fontsize=6,y=0.90)
        
        example_axis.loglog(theory_ls, theory_rsquareds/theory_rsquareds[0]*3e-01,'k-',linewidth=0.3,zorder=0,label='Neutral')
        
        
        leg = example_axis.legend(loc='lower left',frameon=False,title=figure_utils.get_pretty_species_name(species_name))
        leg._legend_box.align = "left"
        
        line, = example_axis.loglog([distances[-1],distances[-1]],[1e-02,1],'k:')
        line.set_dashes((0.5,1))
    
        example_axis.xaxis.get_major_ticks()[-2].label1.set_visible(False)
        example_axis.xaxis.get_major_ticks()[-2].tick1line.set_visible(False)
        
        for tick_idx in xrange(1,7):
        
            example_axis.xaxis.get_minor_ticks()[-tick_idx].tick1line.set_visible(False)
            example_axis.xaxis.get_minor_ticks()[-tick_idx].tick2line.set_visible(False)
    
    #
    # Don't need this anymore, leaving for reference
    #
    # prepare for fit to Richard's curve:
    #
    # Y(t) = K/(1+exp(-B*(t-M)))^(1/v)
    #
    #ts = numpy.log(distances)
    #
    #u = ts # independent variable
    #y = rsquareds # dependent variable
    #w = numpy.sqrt(ns*1.0) # (weights)
    #u = u[distances>5.5]
    #y = y[distances>5.5]
    #w = w[distances>5.5]
    
    
    #median_w = numpy.median(w)
    #median_w = (control_n*1.0)**0.5
    #u = numpy.hstack([u, numpy.array(log(5e04))])
    #y = numpy.hstack([y, numpy.array(control_rsquared)])
    #w = numpy.hstack([w, numpy.array([median_w])])
    
    
    #def model(x,u):
    #    return x[0]/numpy.power(1+numpy.exp(-x[1]*(u-x[2])), 1.0/x[3])+control_rsquared
    # residual function
    # x = model params (K,B,M,v)
    #def fun(x,u,y,w):
        
    #    return (model(x,u)-y)*w
    
    # jacobian (matrix of partial derivatives of residual function)
    #def jac(x,u,y,w):
        
    #    J = numpy.zeros((len(u),len(x)))*1.0
        
        # dY/dK
    #    J[:,0] = w/numpy.power(1+numpy.exp(-x[1]*(u-x[2])), 1.0/x[3])
        # dY/dB
    #    J[:,1] =  w*x[0]*(u-x[2])*numpy.exp(-x[1]*(u-x[2]))/x[3]/numpy.power(1+numpy.exp(-x[1]*(u-x[2])), 1.0/x[3]+1)
        # dY/dM
        #J[:,2] =  -1*w*x[0]*x[1]*numpy.exp(-x[1]*(u-x[2]))/x[3]/numpy.power(1+numpy.exp(-x[1]*(u-x[2])), 1.0/x[3]+1)
        # dY/dv
        #J[:,3] = w*x[0]*numpy.log(1+numpy.exp(-x[1]*(u-x[2])))/x[3]/x[3]/numpy.power(1+numpy.exp(-x[1]*(u-x[2])), 1.0/x[3])
        #return J
    
    #x_lower_bounds = numpy.array([1e-09,-numpy.inf,1e-09,1e-09])
    #x_upper_bounds = numpy.array([1-control_rsquared, -1e-09, log(1e05), numpy.inf])
        
    #x0 = numpy.array([0.5,-1.0,log(100),0.5])
    
    
    #res = least_squares(fun, x0, jac=jac, args=(u,y,w), bounds=(x_lower_bounds, x_upper_bounds), verbose=1)
    
    
    #predicted_rsquareds = model(res.x, u)
    #K = res.x[0]
    #print K
    
    #K = predicted_rsquareds[0]
    #K = rsquareds[0]
    
    #theory_distances = numpy.logspace(0,4,100)
    #theory_u = numpy.log(theory_distances)
    
    #theory_rsquareds = model(res.x, theory_u)
    
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
    
    species_axis.semilogy([species_idx,species_idx], [rsquareds[idx_900], rsquareds[idx_9]],'r-')
    species_axis.semilogy([species_idx], [rsquareds[idx_9]],'r_',markersize=3) #,markeredgewidth=0)
    species_axis.semilogy([species_idx],[rsquareds[idx_90]],'r_',markersize=3) #,markeredgewidth=0)
    species_axis.semilogy([species_idx], [rsquareds[idx_900]],'r_',markersize=3) #,markeredgewidth=0)
    line, = species_axis.semilogy([species_idx,species_idx], [control_rsquared, rsquareds[idx_900]],'r:')
    line.set_dashes((0.5,0.75))
    species_axis.semilogy([species_idx], [control_rsquared],'ro',markersize=2,markeredgewidth=0)


    # now do rbymu estimate
    if rsquareds[-1] < rsquareds[idx_9]/2:
    
        NRstar = calculate_effective_NR(0.5)
        lstar = distances[numpy.fabs(rsquareds/rsquareds[idx_9]-0.5).argmin()]
        rbymu_2 = NRstar/lstar/pi*2
        
        rbymu_axis.semilogy([species_idx], [rbymu_2],'r_',markersize=3) 
        
        if rsquareds[-1] < rsquareds[idx_9]/4:
            
            
            NRstar = calculate_effective_NR(0.25)
            lstar = distances[numpy.fabs(rsquareds/rsquareds[idx_9]-0.25).argmin()]
            rbymu_4 = NRstar/lstar/pi*2
        
            rbymu_axis.semilogy([species_idx,species_idx], [rbymu_4, rbymu_2],'r-')
            rbymu_axis.semilogy([species_idx], [rbymu_4],'r_',markersize=3) 
        

sys.stderr.write("Saving figure...\t")
fit_fig.savefig('%s/supplemental_ld_fits.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')
scaled_fig.savefig('%s/supplemental_ld_shapes.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')
fig.savefig('%s/figure_4.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')
rbymu_fig.savefig('%s/supplemental_rbymu.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')

sys.stderr.write("Done!\n")

     
     
     
     
     
     
     
     
     
     
     
     
     
     
     
