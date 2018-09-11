import matplotlib  
matplotlib.use('Agg') 
import config
import parse_midas_data
import parse_HMP_data


import pylab
import sys
import numpy
import os.path

import diversity_utils
import gene_diversity_utils
import calculate_substitution_rates

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


mpl.rcParams['font.size'] = 4
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
parser.add_argument("--chunk-size", type=int, help="max number of records to load", default=1000000000)

args = parser.parse_args()

debug = args.debug
chunk_size = args.chunk_size

    
################################################################################

min_coverage = config.min_median_coverage
alpha = 0.5 # Confidence interval range for rate estimates
low_pi_threshold = 1e-03
low_divergence_threshold = 1e-03
min_change = 0.8
allowed_variant_types = set(['1D','2D','3D','4D'])
#allowed_variant_types = set(['4D'])


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
alpha = 0.5 # Confidence interval range for rate estimates
low_pi_threshold = 1e-03
low_divergence_threshold = 1e-03
min_change = 0.8
min_sample_size = 33 # 46 gives at least 1000 pairs
allowed_variant_types = set(['1D','2D','3D','4D'])

divergence_matrices = {}
sample_names = {}
good_species_list = parse_midas_data.parse_good_species_list()
if debug:
    good_species_list = good_species_list[0:2]
    
# Load subject and sample metadata
sys.stderr.write("Loading sample metadata...\n")
subject_sample_map = parse_HMP_data.parse_subject_sample_map()
sample_country_map = parse_HMP_data.parse_sample_country_map()
#sample_phenotype_map = parse_HMP_data.parse_sample_phenotype_map()
sys.stderr.write("Done!\n")
 
for species_name in good_species_list:

    sys.stderr.write("Loading haploid samples...\n")
    # Only plot samples above a certain depth threshold that are "haploids"
    snp_samples = diversity_utils.calculate_haploid_samples(species_name, debug=debug)
    
    if len(snp_samples) < min_sample_size:
        sys.stderr.write("Not enough haploid samples!\n")
        continue
        
    sys.stderr.write("Calculating unique samples...\n")
    # Only consider one sample per person
    snp_samples = snp_samples[parse_midas_data.calculate_unique_samples(subject_sample_map, sample_list=snp_samples)]


    if len(snp_samples) < min_sample_size:
        sys.stderr.write("Not enough unique samples!\n")
        continue
        
    # Load divergence matrices 
    sys.stderr.write("Loading pre-computed substitution rates for %s...\n" % species_name)
    substitution_rate_map = calculate_substitution_rates.load_substitution_rate_map(species_name)
    sys.stderr.write("Calculating matrix...\n")
    dummy_samples, snp_difference_matrix, snp_opportunity_matrix = calculate_substitution_rates.calculate_matrices_from_substitution_rate_map(substitution_rate_map, 'core', allowed_samples=snp_samples)
    snp_samples = dummy_samples
    sys.stderr.write("Done!\n")

    snp_substitution_matrix = snp_difference_matrix*1.0/(snp_opportunity_matrix+(snp_opportunity_matrix==0))
    
    divergence_matrices[species_name] = snp_substitution_matrix
    sample_names[species_name] = snp_samples
        

species_names = []
sample_sizes = []

for species_name in divergence_matrices.keys():
    species_names.append(species_name)
    sample_sizes.append( divergence_matrices[species_name].shape[0] )

    
# sort in descending order of sample size
# Sort by num haploids    
sample_sizes, species_names = zip(*sorted(zip(sample_sizes, species_names),reverse=True))
    
sys.stderr.write("Postprocessing %d species...\n" % len(species_names))

####################################################
#
# Set up Figure (4 panels, arranged in 2x2 grid)
#
####################################################
#pylab.figure(1,figsize=(24,1.5*20))
pylab.figure(1,figsize=(40,1.5*len(species_names)))
fig = pylab.gcf()
# make three panels panels
#outer_grid  = gridspec.GridSpec(20,1, height_ratios=[1 for s in species_names], hspace=0.25)
outer_grid  = gridspec.GridSpec(len(species_names),1, height_ratios=[1 for s in species_names], hspace=0.25)


##############################################################################
#
# Dendrogram panels
#
##############################################################################

dendrogram_axes = []
#for species_idx in xrange(20,len(species_names)):
for species_idx in xrange(0,len(species_names)):
    
    dendrogram_axis = plt.Subplot(fig, outer_grid[species_idx])

    fig.add_subplot(dendrogram_axis)
    dendrogram_axes.append(dendrogram_axis)

    dendrogram_axis.set_ylim([1e-06,1e-01])
    dendrogram_axis.set_ylabel(species_names[species_idx])
    
    dendrogram_axis.set_xticks([])

    dendrogram_axis.spines['top'].set_visible(False)
    dendrogram_axis.spines['right'].set_visible(False)
    
    dendrogram_axis.get_xaxis().tick_bottom()
    dendrogram_axis.get_yaxis().tick_left()
    dendrogram_axis


##############################################################################
#
# Now do calculations
#
##############################################################################

# print out the numbers and IDs of each sample for each species
outFile=open(os.path.expanduser('~/ben_nandita_hmp_scripts/clade_definitions.txt'),'w')

#for species_idx in xrange(20,len(species_names)):
for species_idx in xrange(0,len(species_names)):
    species_name = species_names[species_idx]
    
    snp_substitution_rate = divergence_matrices[species_name]
    snp_substitution_rate = numpy.clip(snp_substitution_rate,1e-09,10)
    snp_samples = sample_names[species_name]
    
    dendrogram_axis = dendrogram_axes[species_idx-20]


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
        
            if (y0<2e-04) and (y1<2e-04):
                linewidth=0.75
                color='0.4'
            else:
                linewidth=0.3
                color='0.6'
        
            #print x0, '->', x1, '\t',y0, '->', y1       
            dendrogram_axis.semilogy([x0,x1],[y0,y1],'-',color=color,linewidth=linewidth)
        
            if (y0==y_penultimax) and (y1==y_penultimax):
                # it's the cross bar that bridges the two most-diverged clades
                # so plot a root branch to the top of the plot
                xavg = (x0+x1)*0.5
            
                dendrogram_axis.semilogy([xavg,xavg],[y_penultimax, ymax],'-',color=color,linewidth=linewidth)

    leaf_xs = list(sorted(set(leaf_xs)))

    xticks = []
    xticklabels = []
    samples = []

    print species_name
    outFile.write(species_name +'\n')

    for i in xrange(0,len(ddata['ivl'])):
    
        idx = long(ddata['ivl'][i])
        x = leaf_xs[i]
        y = yplotmin
    
        sample = snp_samples[idx]
        xticks.append(x)
        xticklabels.append(str(i))
        samples.append(sample)
    
        print i, sample
        outFile.write(str(i) + ' ' + sample +'\n')

        if sample_country_map[sample]=='United States':
            color = '#deebf7'
            #if sample_phenotype_map[sample]==0:
            #    color = '#9ecae1'
            #elif sample_phenotype_map[sample]==1:
            #    color = '#3182bd'
            #else:
            #    color = '#deebf7'      
        elif sample_country_map[sample]=='United Kingdom':
            color = '#31a354'
        else:
            color = '#de2d26'
        
        dendrogram_axis.plot([x],[y],'o',color=color,markeredgewidth=0,markersize=2)
    
    dendrogram_axis.plot([0],[1e-09],'o',color='#9ecae1',markeredgewidth=0,markersize=2, label='USA (city 1)')
    dendrogram_axis.plot([0],[1e-09],'o',color='#3182bd',markeredgewidth=0,markersize=2, label='USA (city 2)')
    dendrogram_axis.plot([0],[1e-09],'o',color='#de2d26',markeredgewidth=0,markersize=2, label='China')

    #dendrogram_axis.legend(loc='upper right',ncol=3,frameon=False,fontsize=5,numpoints=1)
    
    new_xmax = xs[0]+(xs[-1]-xs[0])*150.0/len(xticks)
    
    dendrogram_axis.set_xticks([])
    dendrogram_axis.set_xlim([xmin,new_xmax])

    

    #dendrogram_axis.plot([xmin,xmax],[7.3e-03,7.3e-03],'k-',linewidth=0.25)

    dendrogram_axis.set_ylim([yplotmin/1.4,yplotmax])
    dendrogram_axis.set_xticks(xticks)
    dendrogram_axis.set_xticklabels(xticklabels)

sys.stderr.write("Saving figure...\t")
fig.savefig('%s/manual_clade_identification.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')
sys.stderr.write("Done!\n")

 
