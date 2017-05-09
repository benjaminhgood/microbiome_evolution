import matplotlib  
matplotlib.use('Agg') 
import parse_midas_data
import parse_patric
import sys
import numpy
from numpy.random import normal
import diversity_utils
import stats_utils
import pylab
import os


# plotting tools
import matplotlib.colors as colors
from matplotlib.colors import LogNorm
import matplotlib.cm as cmx
from math import log10,ceil
import matplotlib as mpl
mpl.rcParams['font.size'] = 8
mpl.rcParams['lines.linewidth'] = 1.0
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'



################################################################################
#
# Standard header to read in argument information
#
################################################################################
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("species_name", help="name of species to process")
parser.add_argument("--debug", help="Loads only a subset of SNPs for speed", action="store_true")
parser.add_argument("--chunk-size", type=int, help="max number of records to load", default=1000000000)
args = parser.parse_args()

species_name = args.species_name
debug = args.debug
chunk_size = args.chunk_size

################################################################################

# get a list of samples
sample_coverage_histograms, samples = parse_midas_data.parse_coverage_distribution(species_name)


# get a list of marker genes
marker_genes = parse_midas_data.load_marker_genes(species_name)

# load the centroid gene map (maps from the centroid to the ref genome)
centroid_gene_map=parse_midas_data.load_centroid_gene_map(species_name)

# invert the centroid gene map:
centroid_gene_map_inverted={}
for key in centroid_gene_map:
    value=centroid_gene_map[key]
    centroid_gene_map_inverted[value]=key

# load the complete centroid gene map
complete_centroid_gene_map=parse_midas_data.load_complete_centroid_gene_map(species_name)


 
# iterate through the marker genes, find all the CNV values for each gene in each centroid
marker_genes_dictionary={}
ref_genome_dictionary={}
for gene in marker_genes:
    #map refernece genome marker gene to the centroid gene
    centroid_gene=centroid_gene_map_inverted[gene]
    # find all the genes in the 95% centroid
    all_genes_in_centroid=complete_centroid_gene_map[centroid_gene]
    # load gene coverage for each gene in the 95% cluster. Each of these genes is its own 99% cluster
    marker_gene_reads=parse_midas_data.parse_99_percent_genes(species_name, samples, all_genes_in_centroid)
    # store all the data for each marker gene in a dictionary for later use
    marker_genes_dictionary[gene]=marker_gene_reads
    # also initialize a dicationary to concatenate all the data for each ref genome across marker genes
    for ref in all_genes_in_centroid:
        ref_genome='.'.join(ref.split('.')[0:2])
        if ref_genome not in ref_genome_dictionary:
            ref_genome_dictionary[ref_genome]=numpy.asarray([])

# fill in the ref_genome_dictionary with info accross all marker genes
for ref_genomes in ref_genome_dictionary.keys():
    for gene in marker_genes:        
        if ref_genomes in marker_genes_dictionary[gene]:
            ref_genome_dictionary[ref_genomes]=numpy.concatenate((ref_genome_dictionary[ref_genomes],marker_genes_dictionary[gene][ref_genomes]))
        else:
            random_ref_genome=marker_genes_dictionary[gene].keys()[0]
            ref_genome_dictionary[ref_genomes]=numpy.concatenate((ref_genome_dictionary[ref_genomes],numpy.zeros_like(marker_genes_dictionary[gene][random_ref_genome])))

# check if there is only 1 ref genome. If so, add a dummy vector of zeros so that I can plot this.
if len(ref_genome_dictionary.keys())==1:
    ref_genomes=ref_genome_dictionary.keys()[0]
    dummy=numpy.zeros_like(ref_genome_dictionary[ref_genomes])
    ref_genome_dictionary['dummy']=dummy

# concatenate all the arrays for each ref genome into a single matrix for plotting later
ref_genome_numpy_array=numpy.asarray([])
for ref_genome in ref_genome_dictionary.keys():
    if len(ref_genome_numpy_array) > 0:
        ref_genome_numpy_array = numpy.row_stack((ref_genome_numpy_array,ref_genome_dictionary[ref_genome]))
    else:
        ref_genome_numpy_array=ref_genome_dictionary[ref_genome]

#plot the distribution of reads for each marker gene vs ref genome
pylab.figure(figsize=(14, 5))   
pylab.gcf().subplots_adjust(bottom=0.5)
plot_no=1
ax=pylab.subplot(1, 1, plot_no) 
pylab.xlabel('marker gene x sample')
pylab.ylabel('reference genome')
pylab.title(species_name)
im=ax.imshow(ref_genome_numpy_array,interpolation='none', origin='low', cmap='jet', extent=[0,len(marker_genes),0,len(all_genes_in_centroid)], aspect='auto', norm=LogNorm(vmin=1e-1, vmax=1e5))
pylab.colorbar(im)
ax.set(yticks=numpy.arange(len(ref_genome_dictionary.keys()))+0.5,yticklabels=ref_genome_dictionary.keys(), xticks=numpy.arange(len(marker_genes)))
ax.set_xticklabels(marker_genes, rotation=90)
pylab.savefig('%s/%s_marker_gene_reads.png' % (parse_midas_data.analysis_directory,species_name),Bbox='tight')
