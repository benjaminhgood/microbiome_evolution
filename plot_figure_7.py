import matplotlib  
matplotlib.use('Agg') 
import config
import parse_midas_data
import parse_HMP_data
import os.path
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
import matplotlib.colors as mcolors

from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster.hierarchy import cophenet
from scipy.cluster.hierarchy import fcluster


mpl.rcParams['font.size'] = 6
mpl.rcParams['lines.linewidth'] = 0.5
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'

species_name = "Bacteroides_vulgatus_57955"

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

min_coverage = config.min_median_coverage
alpha = 0.5 # Confidence interval range for rate estimates
low_pi_threshold = 1e-03
low_divergence_threshold = 1e-03
min_change = 0.8
min_sample_size = 5
allowed_variant_types = set(['1D','2D','3D','4D'])

clade_divergence_threshold = 1e-02
modification_divergence_threshold = 1e-03
include_high_copynum = False
#include_high_copynum = True


intermediate_filename = '%sfigure_7_data.txt' % (parse_midas_data.analysis_directory)


temporal_change_map = {}

if memoize and os.path.isfile(intermediate_filename):


    # load from Disk!  
    file = open(intermediate_filename,"r")
    file.readline() # header
    for line in file:
        items = line.split(";")
        species_name = items[0]
        temporal_changes = []
        for item in items[1:]:
        
            subitems = item.split(",")
            initial_sample = subitems[0].strip()
            final_sample = subitems[1].strip()
            snp_changes = long(subitems[2])
            gene_changes = long(subitems[3])
        
            temporal_changes.append((initial_sample, final_sample, snp_changes, gene_changes))
            
        temporal_change_map[species_name] = temporal_changes
        
    file.close()
    
else:
    
    # Must compete divergence matrix on the fly! 
            
    # Load subject and sample metadata
    sys.stderr.write("Loading sample metadata...\n")
    subject_sample_map = parse_HMP_data.parse_subject_sample_map()
    sample_order_map = parse_HMP_data.parse_sample_order_map()
    sys.stderr.write("Done!\n")
    
    good_species_list = parse_midas_data.parse_good_species_list()

    for species_name in good_species_list:

        # Only plot samples above a certain depth threshold that are "haploids"
        haploid_samples = diversity_utils.calculate_haploid_samples(species_name, debug=debug)

        if len(haploid_samples) < min_sample_size:
            continue

        same_sample_idxs, same_subject_idxs, diff_subject_idxs = parse_midas_data.calculate_ordered_subject_pairs(sample_order_map, haploid_samples)

        snp_samples = set()
        sample_size = 0        
        for sample_pair_idx in xrange(0,len(same_subject_idxs[0])):
   
            i = same_subject_idxs[0][sample_pair_idx]
            j = same_subject_idxs[1][sample_pair_idx]

            snp_samples.add(haploid_samples[i])
            snp_samples.add(haploid_samples[j])
            
            sample_size += 1
            
        snp_samples = list(snp_samples)

        if sample_size < min_sample_size:
            continue

        sys.stderr.write("Proceeding with %d longitudinal comparisons in %d samples!\n" % (sample_size, len(snp_samples)))

        # Load SNP information for species_name
        sys.stderr.write("Loading SNPs for %s...\n" % species_name)
        sys.stderr.write("(not just core genes...)\n")

        snp_difference_matrix = numpy.array([])
        snp_difference_matrix_mutation = numpy.array([])
        snp_difference_matrix_reversion = numpy.array([])
        snp_opportunity_matrix = numpy.array([])

        final_line_number = 0

        while final_line_number >= 0:
    
            sys.stderr.write("Loading chunk starting @ %d...\n" % final_line_number)
            dummy_samples, allele_counts_map, passed_sites_map, final_line_number = parse_midas_data.parse_snps(species_name, debug=debug, allowed_samples=snp_samples,chunk_size=chunk_size,initial_line_number=final_line_number)
            sys.stderr.write("Done! Loaded %d genes\n" % len(allele_counts_map.keys()))
    
        
            # Calculate fixation matrix
            sys.stderr.write("Calculating matrix of snp differences...\n")
            chunk_snp_difference_matrix_mutation, chunk_snp_difference_matrix_reversion, chunk_snp_opportunity_matrix =     diversity_utils.calculate_fixation_matrix_mutation_reversion(allele_counts_map, passed_sites_map, min_change=min_change)   # 
            sys.stderr.write("Done!\n")
    
            if snp_difference_matrix.shape[0]==0:
                snp_difference_matrix = numpy.zeros_like(chunk_snp_difference_matrix_mutation)*1.0
                snp_difference_matrix_mutation = numpy.zeros_like(snp_difference_matrix)*1.0
                snp_difference_matrix_reversion = numpy.zeros_like(snp_difference_matrix)*1.0
                snp_opportunity_matrix = numpy.zeros_like(snp_difference_matrix)*1.0
    
            snp_difference_matrix += chunk_snp_difference_matrix_mutation+chunk_snp_difference_matrix_reversion
            snp_difference_matrix_mutation += chunk_snp_difference_matrix_mutation
            snp_difference_matrix_reversion += chunk_snp_difference_matrix_reversion
    
    
            snp_opportunity_matrix += chunk_snp_opportunity_matrix

            snp_samples = dummy_samples

        snp_substitution_rate = snp_difference_matrix*1.0/(snp_opportunity_matrix+(snp_opportunity_matrix==0))
        
        # Load gene coverage information for species_name
        sys.stderr.write("Loading pangenome data for %s...\n" % species_name)
        gene_samples, gene_names, gene_presence_matrix, gene_depth_matrix,         marker_coverages, gene_reads_matrix = parse_midas_data.parse_pangenome_data(species_name, allowed_samples=snp_samples)
        sys.stderr.write("Done!\n")

        sys.stderr.write("Loaded gene info for %d samples\n" % len(gene_samples))
      
        # Calculate matrix of number of genes that differ
        sys.stderr.write("Calculating matrix of gene differences...\n")
        gene_gain_matrix, gene_loss_matrix, gene_opportunity_matrix = gene_diversity_utils.calculate_coverage_based_gene_hamming_matrix_gain_loss(gene_depth_matrix, marker_coverages, min_log2_fold_change=4, include_high_copynum=False)

        gene_difference_matrix = gene_gain_matrix + gene_loss_matrix

        # Now need to make the gene samples and snp samples match up
        desired_samples = gene_samples


        num_haploids = len(desired_samples)
     
        # Calculate which pairs of idxs belong to the same sample, which to the same subject
        # and which to different subjects
        desired_same_sample_idxs, desired_same_subject_idxs, desired_diff_subject_idxs = parse_midas_data.calculate_ordered_subject_pairs(sample_order_map, desired_samples)

        sys.stderr.write("%d temporal samples\n" % len(desired_same_subject_idxs[0]))

        snp_sample_idx_map = parse_midas_data.calculate_sample_idx_map(desired_samples, snp_samples)
        gene_sample_idx_map = parse_midas_data.calculate_sample_idx_map(desired_samples, gene_samples)
  

        same_sample_snp_idxs  = parse_midas_data.apply_sample_index_map_to_indices(snp_sample_idx_map,  desired_same_sample_idxs)  
        same_sample_gene_idxs = parse_midas_data.apply_sample_index_map_to_indices(gene_sample_idx_map, desired_same_sample_idxs)  


        same_subject_snp_idxs = parse_midas_data.apply_sample_index_map_to_indices(snp_sample_idx_map, desired_same_subject_idxs)  
        same_subject_gene_idxs = parse_midas_data.apply_sample_index_map_to_indices(gene_sample_idx_map, desired_same_subject_idxs)  

        diff_subject_snp_idxs = parse_midas_data.apply_sample_index_map_to_indices(snp_sample_idx_map, desired_diff_subject_idxs)  
        diff_subject_gene_idxs = parse_midas_data.apply_sample_index_map_to_indices(gene_sample_idx_map, desired_diff_subject_idxs)  

        temporal_changes = []
        
        for sample_pair_idx in xrange(0,len(same_subject_snp_idxs[0])):
#    
            snp_i = same_subject_snp_idxs[0][sample_pair_idx]
            snp_j = same_subject_snp_idxs[1][sample_pair_idx]
    
            initial_sample = snp_samples[snp_i] 
            final_sample = snp_samples[snp_j]
            
            snp_changes = snp_difference_matrix[snp_i,snp_j]
    
            gene_i = same_subject_gene_idxs[0][sample_pair_idx]
            gene_j = same_subject_gene_idxs[1][sample_pair_idx]
    
            if marker_coverages[gene_i]<min_coverage or marker_coverages[gene_j]<min_coverage:
        # can't look at gene changes!
                gene_changes = -1
            else:
                gene_changes = gene_difference_matrix[gene_i, gene_j]
        
            temporal_changes.append((initial_sample, final_sample, snp_changes, gene_changes))        

        temporal_change_map[species_name] = temporal_changes
        
        sys.stderr.write("Done with %s!\n" % species_name) 
    
    sys.stderr.write("Done looping over species!\n")
    
    sys.stderr.write("Writing intermediate file...\n")
    file = open(intermediate_filename,"w")
    # header!
    file.write("Species; initial_sample, final_sample, snp_changes, gene_changes; ...\n")
    for species_name in temporal_change_map.keys():
        
        temporal_change_strs = []
        for i in xrange(0,len(temporal_change_map[species_name])):
        
            temporal_change_strs.append( ("%s, %s, %d, %d" % temporal_change_map[species_name][i]) )
            
        output_strs = [species_name]+temporal_change_strs
        output_str = "; ".join(output_strs)
        file.write("%s\n" % output_str)
    file.close()
    
    

species_names = []
sample_sizes = []

for species_name in temporal_change_map.keys():
    species_names.append(species_name)
    sample_sizes.append( len(temporal_change_map[species_name]) )
    
# sort in descending order of sample size
# Sort by num haploids    
sample_sizes, species_names = zip(*sorted(zip(sample_sizes, species_names),reverse=True))
     
sys.stderr.write("Postprocessing %d species!\n" % len(species_names))


cmap_str = 'YlGnBu'
vmin = -2
vmax = 4
cmap = pylab.get_cmap(cmap_str) 

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=-1):
    if n == -1:
        n = cmap.N
    new_cmap = mcolors.LinearSegmentedColormap.from_list(
         'trunc({name},{a:.2f},{b:.2f})'.format(name=cmap.name, a=minval, b=maxval),
         cmap(numpy.linspace(minval, maxval, n)))
    return new_cmap

cmap = truncate_colormap(cmap, 0.25, 1.0)
cNorm  = colors.Normalize(vmin=0, vmax=vmax)
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cmap)



####################################################
#
# Set up Figure (1 panels, arranged in 1x1 grid)
#
####################################################

pylab.figure(1,figsize=(7,1.5))
fig = pylab.gcf()
# make three panels panels
outer_grid  = gridspec.GridSpec(1,2,width_ratios=[50,1],wspace=0.05)

change_axis = plt.Subplot(fig, outer_grid[0])
fig.add_subplot(change_axis)

change_axis.set_ylabel('Number of samples')
change_axis.set_ylim([-35,35])
change_axis.set_xlim([-1,len(species_names)])
change_axis.plot([-1,len(species_names)],[0,0],'k-')

change_axis.set_yticks([-30,-20,-10,0,10,20,30])
change_axis.set_yticklabels(['30','20','10','0','10','20','30'])

xticks = numpy.arange(0,len(species_names))
#xticklabels = ["%s (%d)" % (species_names[i],sample_sizes[i]) for i in xrange(0,len(sample_sizes))]
xticklabels = ["%s" % (species_names[i]) for i in xrange(0,len(sample_sizes))]

change_axis.set_xticks(xticks)
change_axis.set_xticklabels(xticklabels, rotation='vertical',fontsize=4)

#change_axis.spines['top'].set_visible(False)
#change_axis.spines['right'].set_visible(False)
change_axis.get_xaxis().tick_bottom()
change_axis.get_yaxis().tick_left()

cax = plt.Subplot(fig, outer_grid[1])
fig.add_subplot(cax)
 

##############################################################################
#
# Now do calculations
#
##############################################################################

 
# Plot percentiles of divergence distribution
for species_idx in xrange(0,len(species_names)):

    species_name = species_names[species_idx]
    temporal_changes = temporal_change_map[species_name]
     
    snp_changes = []
    gene_changes = []
    for sample_1, sample_2, num_snps, num_genes in temporal_changes:
         
         snp_changes.append(num_snps)
         
         if num_genes>=0:
             gene_changes.append(num_genes)
             
    snp_changes = numpy.array(snp_changes)
    snp_changes.sort()
    gene_changes = numpy.array(gene_changes)
    gene_changes.sort()
    
    for idx in xrange(0,len(snp_changes)):
        
        if snp_changes[idx]<0.5:
            colorVal='0.7'
        else:
        
            colorVal = scalarMap.to_rgba(log10(snp_changes[idx]))
        
        change_axis.fill_between([species_idx-0.3,species_idx+0.3],[idx,idx],[idx+1.05,idx+1.05],color=colorVal,linewidth=0)
    
    for idx in xrange(0,len(gene_changes)):
        
        if gene_changes[idx]<0.5:
            colorVal='0.7'
        else:
            colorVal = scalarMap.to_rgba(log10(gene_changes[idx]))
        
        change_axis.fill_between([species_idx-0.3,species_idx+0.3],[-idx-1.05,-idx-1.05],[-idx,-idx],color=colorVal,linewidth=0)
            

m = change_axis.scatter([200],[1],c=[0], vmin=0, vmax=vmax, cmap=cmap, marker='^')


cbar = fig.colorbar(m,cax=cax,orientation='vertical', ticks=[0,1,2,3,4])
cbar.set_ticklabels(['$1$','$10$','$10^{2}$','$10^{3}$','$10^{4}$'])
cbar.set_label('Number of changes',rotation=270,labelpad=10)
cl = pylab.getp(cbar.ax, 'ymajorticklabels')
pylab.setp(cl, fontsize=9) 
#fig.text(0.945,0.05,'$\\pi/\\pi_0$',fontsize=12)

cbar.ax.tick_params(labelsize=5)
change_axis.text(20,25,'SNPs',fontsize=5)
change_axis.text(20,-20,'Genes',fontsize=5)


sys.stderr.write("Saving figure...\t")
fig.savefig('%s/figure_7.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')
sys.stderr.write("Done!\n")

 