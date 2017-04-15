import matplotlib  
matplotlib.use('Agg') 
import parse_midas_data
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

mpl.rcParams['font.size'] = 9
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
parser.add_argument("species_name", help="name of species to process")
parser.add_argument("--debug", help="Loads only a subset of SNPs for speed", action="store_true")
parser.add_argument("--chunk-size", type=int, help="max number of records to load", default=1000000000)
args = parser.parse_args()

species_name = args.species_name
debug = args.debug
chunk_size = args.chunk_size
################################################################################

min_change = 0.8
min_coverage = 20
alpha = 0.5 # Confidence interval range for rate estimates


# Load subject and sample metadata
sys.stderr.write("Loading HMP metadata...\n")
subject_sample_map = parse_midas_data.parse_subject_sample_map()
sys.stderr.write("Done!\n")

# Load reference genes
reference_genes = parse_midas_data.load_reference_genes(species_name)
print len(reference_genes), "reference genes"

# Load metaphlan2 genes
metaphlan2_genes = parse_midas_data.load_metaphlan2_genes(species_name)
print len(metaphlan2_genes), "metaphlan2 genes"

# Load MIDAS marker genes
marker_genes = parse_midas_data.load_marker_genes(species_name)
print len(marker_genes), "marker genes"

# Load gene coverage information for species_name
sys.stderr.write("Loading pangenome data for %s...\n" % species_name)
gene_samples, gene_names, gene_presence_matrix, gene_depth_matrix, marker_coverages, gene_reads_matrix = parse_midas_data.parse_pangenome_data(species_name)
sys.stderr.write("Done!\n")

pangenome_genes = set(gene_names)

for marker_gene in marker_genes:
    print marker_gene, marker_gene in pangenome_genes
    
reference_gene_idxs = numpy.array([gene_name in reference_genes for gene_name in gene_names])
metaphlan2_gene_idxs = numpy.array([gene_name in metaphlan2_genes for gene_name in gene_names])
marker_gene_idxs = numpy.array([gene_name in marker_genes for gene_name in gene_names])

print marker_genes

print marker_gene_idxs.sum()

sample_idxs = (parse_midas_data.calculate_unique_samples(subject_sample_map, gene_samples))*(marker_coverages>=min_coverage)

prevalences = gene_diversity_utils.calculate_fractional_gene_prevalences(gene_depth_matrix[:,sample_idxs], marker_coverages[sample_idxs],min_copynum=0.3)

reference_prevalences = prevalences[reference_gene_idxs]
metaphlan2_prevalences = prevalences[metaphlan2_gene_idxs]
marker_prevalences = prevalences[marker_gene_idxs]

print marker_prevalences

pangenome_xs, pangenome_survivals = stats_utils.calculate_unnormalized_survival_from_vector(prevalences, min_x=0, max_x=1)

reference_xs, reference_survivals = stats_utils.calculate_unnormalized_survival_from_vector(reference_prevalences, min_x=0, max_x=1)

metaphlan2_xs, metaphlan2_survivals = stats_utils.calculate_unnormalized_survival_from_vector(metaphlan2_prevalences, min_x=0, max_x=1)

marker_xs, marker_survivals = stats_utils.calculate_unnormalized_survival_from_vector(marker_prevalences, min_x=0, max_x=1)

pylab.figure(1,figsize=(3.42,4))
pylab.title(species_name)

#pylab.step(pangenome_xs, pangenome_survivals/pangenome_survivals[0],label='Pan-genome')
#pylab.step(reference_xs, reference_survivals/reference_survivals[0],label='Reference')
#pylab.step(metaphlan2_xs, metaphlan2_survivals/metaphlan2_survivals[0],label='Metaphlan2')
#pylab.step(marker_xs, marker_survivals/marker_survivals[0],label='MIDAS Marker')
#pylab.ylim([1e-02,1])

pylab.step(pangenome_xs, pangenome_survivals, label='Pangenome')
pylab.step(reference_xs, reference_survivals,label='Reference')
pylab.step(metaphlan2_xs, metaphlan2_survivals,label='Metaphlan2')
pylab.step(marker_xs, marker_survivals,label='MIDAS Marker')
pylab.ylim([10,1e04])

pylab.semilogy([1],[1])

pylab.legend(frameon=False,loc='lower left', fontsize=7)
pylab.xlabel('Prevalence fraction, $p$ (from genes/)')
pylab.ylabel('Fraction of genes >= $p$')

pylab.savefig('%s/%s_prevalence_distribution.pdf' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight',transparent=True)
pylab.savefig('%s/%s_prevalence_distribution.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight',dpi=300)
