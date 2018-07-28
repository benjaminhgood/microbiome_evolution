import matplotlib  
matplotlib.use('Agg') 
import sample_utils
import config
import parse_midas_data


import pylab
import sys
import numpy

import diversity_utils
import gene_diversity_utils
import calculate_substitution_rates

import clade_utils
import figure_utils

import stats_utils
import matplotlib.colors as colors
import matplotlib.cm as cmx
from math import log10,ceil
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy.random import randint

from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster.hierarchy import cophenet
from scipy.cluster.hierarchy import fcluster


mpl.rcParams['font.size'] = 6
mpl.rcParams['lines.linewidth'] = 0.5
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'

examples = ["Bacteroides_vulgatus_57955",  "Bacteroides_stercoris_56735"]

################################################################################
#
# Standard header to read in argument information
#
################################################################################
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--debug", help="Loads only a subset of SNPs for speed", action="store_true")
parser.add_argument("--chunk-size", type=int, help="max number of records to load", default=1000000000)

args = parser.parse_args()

debug = args.debug
chunk_size = args.chunk_size
    
################################################################################

min_coverage = config.min_median_coverage
alpha = 0.05 # FDR for confidence interval range for rate estimates
low_divergence_threshold = 2e-04
allowed_variant_types = set(['1D','2D','3D','4D'])
num_random_pairs = 500

####################################################
#
# Set up Figure (4 panels, arranged in 2x2 grid)
#
####################################################

pylab.figure(1,figsize=(6,2.2))
divergence_fig = pylab.gcf()
pylab.figure(2,figsize=(6,2.2))
dendrogram_fig = pylab.gcf()

divergence_grid = gridspec.GridSpec(2, 2, width_ratios=[4,1], wspace=0.025, height_ratios=[1,1], hspace=0.5)
                 
dendrogram_grid = gridspec.GridSpec(2, 2, width_ratios=[8,1], wspace=0.025, height_ratios=[1,1], hspace=0.5)


snp_axes = []
zoomed_snp_axes = []
###################
#
# SNP change panel
#
###################
for example_idx in xrange(0,len(examples)):

    snp_axis = plt.Subplot(divergence_fig, divergence_grid[example_idx,0])
    divergence_fig.add_subplot(snp_axis)

    snp_axis.set_ylabel('Divergence, $d$')
    snp_axis.set_xlabel('%s QP host pairs (random, $n=%d$)' % (figure_utils.get_abbreviated_species_name(examples[example_idx]), num_random_pairs))
    snp_axis.set_ylim([1e-06,1e-01])
    snp_axis.set_xlim([-512,5])
    snp_axis.set_xticks([])

    snp_axis.spines['top'].set_visible(False)
    snp_axis.spines['right'].set_visible(False)
    snp_axis.get_xaxis().tick_bottom()
    snp_axis.get_yaxis().tick_left()

    snp_axes.append(snp_axis)

    zoomed_snp_axis = plt.Subplot(divergence_fig, divergence_grid[example_idx, 1])
    divergence_fig.add_subplot(zoomed_snp_axis)

    zoomed_snp_axis.semilogy([1,1],[1e-08,1e-08])
    zoomed_snp_axis.set_ylim([1e-06,1e-01])
    zoomed_snp_axis.set_xlim([-110,1])

    zoomed_snp_axis.set_xlabel('100 closest pairs')
    zoomed_snp_axis.set_xticks([])
    zoomed_snp_axis.set_yticklabels([])

    zoomed_snp_axis.spines['top'].set_visible(False)
    zoomed_snp_axis.spines['right'].set_visible(False)
    zoomed_snp_axis.get_xaxis().tick_bottom()
    zoomed_snp_axis.get_yaxis().tick_left()

    zoomed_snp_axes.append(zoomed_snp_axis)
    
##############################################################################
#
# Dendrogram panel
#
##############################################################################

dendrogram_axes = []
###################
#
# SNP change panel
#
###################
for example_idx in xrange(0,len(examples)):

    dendrogram_axis = plt.Subplot(dendrogram_fig, dendrogram_grid[example_idx,0])

    dendrogram_fig.add_subplot(dendrogram_axis)

    dendrogram_axis.set_ylim([1e-06,1e-01])
    dendrogram_axis.set_ylabel('Divergence, $d$')
    dendrogram_axis.set_xlabel('%s QP hosts' % figure_utils.get_abbreviated_species_name(examples[example_idx]))

    dendrogram_axis.set_xticks([])

    dendrogram_axis.spines['top'].set_visible(False)
    dendrogram_axis.spines['right'].set_visible(False)
    dendrogram_axis.spines['bottom'].set_visible(False)

    dendrogram_axis.get_xaxis().tick_bottom()
    dendrogram_axis.get_yaxis().tick_left()
    
    dendrogram_axes.append(dendrogram_axis)
    
##############################################################################
#
# Phylogenetic inconsistency panel
#
##############################################################################
inconsistency_axes = []
for example_idx in xrange(0,len(examples)):

    inconsistency_axis = plt.Subplot(dendrogram_fig, dendrogram_grid[example_idx,1])
    dendrogram_fig.add_subplot(inconsistency_axis)

    inconsistency_axis.set_xlabel('% inconsistent')
    inconsistency_axis.set_yticklabels([])

    inconsistency_axis.spines['top'].set_visible(False)
    inconsistency_axis.spines['right'].set_visible(False)
    inconsistency_axis.get_xaxis().tick_bottom()
    inconsistency_axis.get_yaxis().tick_left()

    inconsistency_axes.append(inconsistency_axis)

##############################################################################
#
# Now do calculations
#
##############################################################################



# Load subject and sample metadata
sys.stderr.write("Loading sample metadata...\n")
subject_sample_map = sample_utils.parse_subject_sample_map()
sample_country_map = sample_utils.parse_sample_country_map()
sample_continent_map = sample_utils.parse_sample_continent_map()

sys.stderr.write("Done!\n")

#continent_color_map = {'North America': '#deebf7', 'Asia': '#de2d26', 'Europe' : 'g'}
continent_color_map = {'North America': '#8da0cb', 'Europe': '#66c2a5', 'Asia' : '#fc8d62'}

       
for example_idx in xrange(0,len(examples)):

    species_name = examples[example_idx]
    snp_axis = snp_axes[example_idx]
    zoomed_snp_axis = zoomed_snp_axes[example_idx]
    dendrogram_axis = dendrogram_axes[example_idx]
    inconsistency_axis = inconsistency_axes[example_idx]
    
    # Only plot samples above a certain depth threshold that are "haploids"
    snp_samples = diversity_utils.calculate_haploid_samples(species_name, debug=debug)
    # Only consider one sample per person
    snp_samples =     snp_samples[sample_utils.calculate_unique_samples(subject_sample_map, sample_list=snp_samples)]
    sys.stderr.write("Proceeding with %d haploid samples!\n" % len(snp_samples))


    sys.stderr.write("Loading pre-computed substitution rates for %s...\n" % species_name)
    substitution_rate_map = calculate_substitution_rates.load_substitution_rate_map(species_name)
    sys.stderr.write("Calculating matrix...\n")
    dummy_samples, snp_difference_matrix, snp_opportunity_matrix = calculate_substitution_rates.calculate_matrices_from_substitution_rate_map(substitution_rate_map, 'core', allowed_samples=snp_samples)
    snp_samples = dummy_samples
    sys.stderr.write("Done!\n")

    snp_substitution_rate = snp_difference_matrix*1.0/(snp_opportunity_matrix+(snp_opportunity_matrix==0))
    snp_substitution_rate = numpy.clip(snp_substitution_rate,1e-09,10)

    sys.stderr.write("Done!\n")   


    sys.stderr.write("Calculating UPGMA dendrogram...\n")
    # calculate compressed distance matrix suitable for agglomerative clustering
    Y = []
    for i in xrange(0,snp_substitution_rate.shape[0]):
        for j in xrange(i+1,snp_substitution_rate.shape[1]):
            Y.append(snp_substitution_rate[i,j]) 
    Y = numpy.array(Y) 
    Z = linkage(Y, method='average')           
    c, coph_dists = cophenet(Z, Y)
    ddata = dendrogram(Z, no_plot=True)
    sys.stderr.write("Done! cophenetic correlation: %g\n" % c)

    # Calculate which pairs of idxs belong to the same sample, which to the same subject
    # and which to different subjects
    same_sample_idxs, same_subject_idxs, diff_subject_idxs =     sample_utils.calculate_subject_pairs(subject_sample_map, snp_samples)

    diff_subject_snp_plowers = []
    diff_subject_snp_puppers = []
    diff_subject_gene_plowers = []
    diff_subject_gene_puppers = []
    between_host_gene_idx_map = {}
    low_divergence_between_host_gene_idx_map = {}
    for sample_pair_idx in xrange(0,len(diff_subject_idxs[0])):
    
        snp_i = diff_subject_idxs[0][sample_pair_idx]
        snp_j = diff_subject_idxs[1][sample_pair_idx]
    
        plower,pupper = stats_utils.calculate_poisson_rate_interval(snp_difference_matrix[snp_i,snp_j], snp_opportunity_matrix[snp_i, snp_j],alpha=0.05)
    
        diff_subject_snp_plowers.append(plower)
        diff_subject_snp_puppers.append(pupper)
            
    # clip lower bounds 
    diff_subject_snp_plowers = numpy.clip(diff_subject_snp_plowers,1e-07,1e09)
    # Sort all four lists by ascending lower bound on SNP changes, then gene changes
    diff_subject_snp_plowers, diff_subject_snp_puppers = (numpy.array(x) for x in zip(*sorted(zip(diff_subject_snp_plowers, diff_subject_snp_puppers))))


    ##############################################################################
    #
    # Now plot figures
    #
    ##############################################################################


    idxs = randint(0,len(diff_subject_snp_plowers), 500)
    idxs.sort()

    y=0
    for idx in idxs:
    
        snp_plower = diff_subject_snp_plowers[idx]
        snp_pupper = diff_subject_snp_puppers[idx]
    
        y-=1
    
        snp_axis.semilogy([y,y],[snp_plower,snp_pupper],'r-',linewidth=0.35)

    y=0
    idxs = numpy.arange(0,100)
    for idx in idxs:
    
        snp_plower = diff_subject_snp_plowers[idx]
        snp_pupper = diff_subject_snp_puppers[idx]
    
        y-=1
    
        zoomed_snp_axis.semilogy([y,y],[snp_plower,snp_pupper],'r-',linewidth=0.35)

        
    zoomed_snp_axis.set_yticklabels([])

    #################################################
    #
    # Plot dendrogram figure
    #
    #######

    # calculate second minimum y value
    ys = []
    xs = []
    for i, d in zip(ddata['icoord'], ddata['dcoord']):
        ys.extend(d)
        xs.extend(i)
    
    xs = list(set(xs))
    xs.sort()
    xs = numpy.array(xs)

    dx = xs[-1]-xs[0]
    xmin = xs[0]-dx*0.025
    xmax = xs[-1]+dx*0.025

    ys = list(set(ys))
    ys.sort()
    ys = numpy.array(ys)

    if ys[0]<1e-09:
        y_penultimin = ys[1]/2
    else:
        y_penultimin = ys[0]/2

    y_penultimax = ys[-1]

    ymin = 1e-06
    #ymin=2e-10
    ymax=1e-01

    yplotmin = 1e-06
    yplotmax = 1e-01


    leaf_xs = []

    for icoord, dcoord in zip(ddata['icoord'], ddata['dcoord']):
        for idx in xrange(0,len(icoord)-1):
            x0 = icoord[idx]
            y0 = dcoord[idx]
            if y0<1e-10:
                y0 = ymin
            x1 = icoord[idx+1]
            y1 = dcoord[idx+1]
            if y1<1e-10:
                y1 = ymin
        
            if (y0==ymin):
                leaf_xs.append(x0)
        
            if (y1==ymin):
                leaf_xs.append(x1)
        
            if False: #(y0<low_divergence_threshold) and (y1<low_divergence_threshold):
                linewidth=0.75
                color='0.4'
            else:
                linewidth=0.3
                color='0.6'
        
            dendrogram_axis.semilogy( [x0,x1],[y0,y1],'-',color=color,linewidth=linewidth)
        
            if (y0==y_penultimax) and (y1==y_penultimax):
                # it's the cross bar that bridges the two most-diverged clades
                # so plot a root branch to the top of the plot
                xavg = (x0+x1)*0.5
            
                dendrogram_axis.semilogy([xavg,xavg],[y_penultimax, ymax],'-',color=color,linewidth=linewidth)

    leaf_xs = list(sorted(set(leaf_xs)))


    for i in xrange(0,len(ddata['ivl'])):
    
        idx = long(ddata['ivl'][i])
        x = leaf_xs[i]
        y = yplotmin
    
        sample = snp_samples[idx]
    
        color = continent_color_map[sample_continent_map[sample]]
        
        dendrogram_axis.plot([x],[y],'o',color=color,markeredgewidth=0,markersize=1.3)
    
    for continent in continent_color_map:
        color = continent_color_map[continent]
        dendrogram_axis.plot( [0],[1e-09],'o',color=color,markeredgewidth=0,markersize=3, label=continent)
        
    if example_idx==0:
        dendrogram_axis.legend(loc='upper right',ncol=3,frameon=False,fontsize=5,numpoints=1)
    
    dendrogram_axis.set_xticks([])
    dendrogram_axis.set_xlim([xmin,xmax])
    
    #line, = dendrogram_axis.plot([xmin,xmax], [low_divergence_threshold, low_divergence_threshold],'-',color='k',linewidth=0.25)
    #line.set_dashes((1,1))


    #dendrogram_axis.plot([xmin,xmax],[7.3e-03,7.3e-03],'k-',linewidth=0.25)

    dendrogram_axis.set_ylim([yplotmin/1.4,yplotmax])

    snp_samples = numpy.array(snp_samples)

    # Calculate range of distances to compute phylogenetic inconsistency for
    
    min_d = snp_substitution_rate[(snp_substitution_rate>=(low_divergence_threshold/10))].min() +1e-07 # a little bit more than the min
    max_d = snp_substitution_rate.max()-1e-07 # a little bit less than the max
    ds = numpy.logspace(log10(min_d),log10(max_d),15) # 15 points are plotted
    print ds

    clade_setss = []
    sys.stderr.write("Assessing phylogenetic consistency...\n")
    sys.stderr.write("Clustering samples at different divergence thresholds...\n")
    for d in ds:
        clade_idx_sets = clade_utils.cluster_samples_within_clades(snp_substitution_rate, d=d)
        clade_sets = [snp_samples[numpy.array(list(clade_idxs))] for clade_idxs in clade_idx_sets]
    
        clade_setss.append(clade_sets)

    dstar = low_divergence_threshold

    manual_clade_sets = clade_utils.load_manual_clades(species_name)
    manual_clade_idxss = clade_utils.calculate_clade_idxs_from_clade_sets(snp_samples, clade_sets)

    clade_idx_sets = clade_utils.cluster_samples_within_clades(snp_substitution_rate, clade_idxss=manual_clade_idxss, d=dstar)
    fine_clade_sets = [snp_samples[numpy.array(list(clade_idxs))] for clade_idxs in clade_idx_sets]

    num_bootstraps = 1
    bootstrapped_fine_clade_idxsss = []

    sys.stderr.write("Loading core genes...\n")
    core_genes = parse_midas_data.load_core_genes(species_name)
    sys.stderr.write("Done! Core genome consists of %d genes\n" % len(core_genes))


    # Load SNP information for species_name
    sys.stderr.write("Loading SNPs for %s...\n" % species_name)


    total_singleton_sites = [0 for d in ds]
    total_polymorphic_sites = [0 for d in ds]
    total_inconsistent_sites = [0 for d in ds]
    total_null_inconsistent_sites = [0 for d in ds]

    total_dstar_singleton_sites = 0
    total_dstar_polymorphic_sites = 0
    total_bootstrapped_singleton_sites = 0
    total_bootstrapped_polymorphic_sites = 0
 
    singleton_variant_types = [{'1D':0, '2D':0, '3D':0, '4D':0} for d in ds]
    polymorphic_variant_types = [{'1D':0, '2D':0, '3D':0, '4D':0} for d in ds]
    inconsistent_variant_types = [{'1D':0, '2D':0, '3D':0, '4D':0} for d in ds]
 
    polymorphic_freqs = [[] for d in ds]
    inconsistent_freqs = [[] for d in ds]
 

    final_line_number = 0
    while final_line_number >= 0:
    
        sys.stderr.write("Loading chunk starting @ %d...\n" % final_line_number)
        dummy_samples, allele_counts_map, passed_sites_map, final_line_number = parse_midas_data.parse_snps(species_name, debug=debug, allowed_variant_types=allowed_variant_types, allowed_samples=snp_samples,chunk_size=chunk_size,initial_line_number=final_line_number, allowed_genes=core_genes)
        snp_samples = dummy_samples
    
        sys.stderr.write("Done! Loaded %d genes\n" % len(allele_counts_map.keys()))
    
        sys.stderr.write("Calculating phylogenetic consistency...\n")
        for i in xrange(0,len(ds)):
        
            clade_sets = clade_setss[i]
            cluster_idxss = clade_utils.calculate_clade_idxs_from_clade_sets(snp_samples, clade_sets)
        
            #print cluster_idxss
        
            chunk_singleton_freqs, chunk_polymorphic_freqs, chunk_inconsistent_freqs, chunk_null_inconsistent_freqs, chunk_singleton_variant_types, chunk_polymorphic_variant_types, chunk_inconsistent_variant_types, chunk_null_variant_types = clade_utils.calculate_phylogenetic_consistency(allele_counts_map, passed_sites_map, cluster_idxss, allowed_genes=core_genes)    
        
            total_singleton_sites[i] += len(chunk_singleton_freqs)
            total_polymorphic_sites[i] += len(chunk_polymorphic_freqs)+len(chunk_singleton_freqs) 
            total_inconsistent_sites[i] += len(chunk_inconsistent_freqs)
            total_null_inconsistent_sites[i] += len(chunk_null_inconsistent_freqs)
    
            polymorphic_freqs[i].extend(chunk_polymorphic_freqs)
            inconsistent_freqs[i].extend(chunk_inconsistent_freqs)
            
            print "Singleton:", chunk_singleton_variant_types
            print "Polymorphic:", chunk_polymorphic_variant_types
            print "Inconsistent:", chunk_inconsistent_variant_types
                
            for variant_type in polymorphic_variant_types[i].keys():
                singleton_variant_types[i][variant_type] += chunk_singleton_variant_types[variant_type]
                polymorphic_variant_types[i][variant_type] += chunk_polymorphic_variant_types[variant_type]
                inconsistent_variant_types[i][variant_type] += chunk_inconsistent_variant_types[variant_type]
    
        sys.stderr.write("Calculating singleton fractions...\n")
        if len(bootstrapped_fine_clade_idxsss)==0:
        
            manual_clade_idxss = clade_utils.calculate_clade_idxs_from_clade_sets(snp_samples, clade_sets)

            fine_clade_idxss = clade_utils.calculate_clade_idxs_from_clade_sets(snp_samples, fine_clade_sets)
        
            for bootstrap_idx in xrange(0,num_bootstraps):
            
                permuted_idxs = clade_utils.permute_idxs_within_clades(manual_clade_idxss)
        
                bootstrapped_fine_clade_idxss = []
                for clade_idxs in fine_clade_idxss:
                    bootstrapped_fine_clade_idxss.append( numpy.array(clade_idxs[permuted_idxs],copy=True) )
                bootstrapped_fine_clade_idxsss.append( bootstrapped_fine_clade_idxss )
            
        chunk_singleton_freqs, chunk_polymorphic_freqs, chunk_inconsistent_freqs, chunk_null_inconsistent_freqs, chunk_singleton_variant_types, chunk_polymorphic_variant_types, chunk_inconsistent_variant_types, chunk_null_variant_types = clade_utils.calculate_phylogenetic_consistency(allele_counts_map, passed_sites_map, fine_clade_idxss, allowed_genes=core_genes)    
        
        total_dstar_singleton_sites += len(chunk_singleton_freqs)
        total_dstar_polymorphic_sites += len(chunk_polymorphic_freqs)+len(chunk_singleton_freqs)    
        
        for bootstrap_idx in xrange(0,num_bootstraps):

            chunk_singleton_freqs, chunk_polymorphic_freqs, chunk_inconsistent_freqs, chunk_null_inconsistent_freqs, chunk_singleton_variant_types, chunk_polymorphic_variant_types, chunk_inconsistent_variant_types, chunk_null_variant_types = clade_utils.calculate_phylogenetic_consistency(allele_counts_map, passed_sites_map, bootstrapped_fine_clade_idxsss[bootstrap_idx], allowed_genes=core_genes)    
        
            total_bootstrapped_singleton_sites += len(chunk_singleton_freqs)
            total_bootstrapped_polymorphic_sites += len(chunk_polymorphic_freqs)+len(chunk_singleton_freqs)       
    
        sys.stderr.write("Done!\n")
        
    total_polymorphic_sites = numpy.array(total_polymorphic_sites)
    total_inconsistent_sites = numpy.array(total_inconsistent_sites)
    total_null_inconsistent_sites = numpy.array(total_null_inconsistent_sites)
  
    fraction_inconsistent = total_inconsistent_sites*1.0/(total_polymorphic_sites+(total_polymorphic_sites==0))
    null_fraction_inconsistent = total_null_inconsistent_sites*1.0/(total_polymorphic_sites+(total_polymorphic_sites==0))

    good_inconsistency_idxs = (total_polymorphic_sites>0)
    
    sys.stderr.write("Done!\n")

    sys.stderr.write("Observed singletons at %g: %d/%d (%g), Expected: %d/%d (%g)\n" % (dstar, total_dstar_singleton_sites, total_dstar_polymorphic_sites, total_dstar_singleton_sites*1.0/(total_dstar_polymorphic_sites+(total_dstar_polymorphic_sites==0)), total_bootstrapped_singleton_sites, total_bootstrapped_polymorphic_sites, total_bootstrapped_singleton_sites*1.0/(total_bootstrapped_polymorphic_sites + (total_bootstrapped_polymorphic_sites==0))))

    inconsistency_axis.semilogy([2,2],[yplotmin/1.4, yplotmax],'-',color='0.7', linewidth=0.25)
    inconsistency_axis.set_ylim([yplotmin/1.4,yplotmax])
    inconsistency_axis.set_xlim([0,1])
    inconsistency_axis.set_yticklabels([])
    inconsistency_axis.set_xticks([0,0.5,1])
    inconsistency_axis.plot(fraction_inconsistent[good_inconsistency_idxs], ds[good_inconsistency_idxs], 'r.-',markersize=2,zorder=1,label='Observed')
    inconsistency_axis.plot(null_fraction_inconsistent[good_inconsistency_idxs], ds[good_inconsistency_idxs], '.-',color='0.7',markersize=2,zorder=0,label='Unlinked')

    if example_idx==0:
        inconsistency_axis.legend(loc='lower right',frameon=False,numpoints=1,fontsize=4)

sys.stderr.write("Saving figure...\t")
divergence_fig.savefig('%s/supplemental_divergence_distribution.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')
dendrogram_fig.savefig('%s/supplemental_trees.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')
sys.stderr.write("Done!\n")

 
