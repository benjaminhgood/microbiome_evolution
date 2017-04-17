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

mpl.rcParams['font.size'] = 10
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

desired_samples = set(['700037453', '700037539'])


# Load subject and sample metadata
sys.stderr.write("Loading HMP metadata...\n")
subject_sample_map = parse_midas_data.parse_subject_sample_map()
sys.stderr.write("Done!\n")

# Load marker gene coverages (from species/ folder)
species_coverage_matrix, sample_list, species_list = parse_midas_data.parse_marker_gene_coverages(species_name)
species_marker_coverages = species_coverage_matrix[species_list.index(species_name),:]

# Load gene coverages on reference genome
gene_coverages, samples = parse_midas_data.parse_gene_coverages(species_name)

# Load metaphlan2 genes
metaphlan2_genes = parse_midas_data.load_metaphlan2_genes(species_name)
print len(metaphlan2_genes), "metaphlan2 genes"

metaphlan2_gene_coverages = []
for gene in metaphlan2_genes:
    if gene in gene_coverages:
        metaphlan2_gene_coverages.append( gene_coverages[gene] )
metaphlan2_gene_coverages = numpy.array(metaphlan2_gene_coverages)

median_metaphlan2_coverages = numpy.median(metaphlan2_gene_coverages,axis=0)
#mean_metaphlan2_coverages = metaphlan2_gene_coverages.mean(axis=0)
mean_metaphlan2_coverages = (metaphlan2_gene_coverages*(metaphlan2_gene_coverages>=1)).sum(axis=0)/((metaphlan2_gene_coverages>=1).sum(axis=0))


# Load gene coverage information for species_name
sys.stderr.write("Loading pangenome data for %s...\n" % species_name)
gene_samples, gene_names, gene_presence_matrix, gene_depth_matrix, marker_coverages, gene_reads_matrix = parse_midas_data.parse_pangenome_data(species_name)
sys.stderr.write("Done!\n")

midas_marker_coverage_map = {}
for i in xrange(0,len(gene_samples)):
    midas_marker_coverage_map[gene_samples[i]] = marker_coverages[i]
    
pylab.figure(1,figsize=(22,2))
 
marker_gene_coverages = parse_midas_data.parse_marker_gene_coverage_distribution(species_name)

max_coverages = []
median_coverages = []

marker_genes = sorted(marker_gene_coverages[marker_gene_coverages.keys()[0]].keys())

for sample in desired_samples:
    current_loc = 0
    pooled_coverages = []
    
    line, = pylab.plot([-1],[1],'-')
    color = pylab.getp(line,'color')
    
    print sample
    
    for gene_idx in xrange(0,len(marker_genes)):
        gene_name = marker_genes[gene_idx]
        current_loc += 1
        pylab.plot([current_loc,current_loc],[1e-09,1e09],'k-')
        current_loc += 1
        
        if gene_idx%2==0:
            pylab.text(current_loc+10, 900, gene_name, fontsize=6)
        else:
            pylab.text(current_loc+10, 800, gene_name, fontsize=6)
        
        locs,depths = marker_gene_coverages[sample][gene_name]
        
        if len(locs)==0:
            print gene_name, "has no sites!"
            continue
        
        pylab.plot(locs+current_loc, depths, '-', color=color)
        
        current_loc += locs[-1]
        pooled_coverages.extend(depths)
        
        print gene_name, numpy.median(depths)
        
    
    pooled_coverages = numpy.array(pooled_coverages)
    median_coverage = numpy.median(pooled_coverages[pooled_coverages>0])
    
    pylab.plot([-1],[1],'-',color=color,label=("%s (Median(>0)=%g, genes/=%g)" % (sample, median_coverage, midas_marker_coverage_map[sample])))
    
    
    print "Pooled:", median_coverage
    
    median_coverages.append( median_coverage  )
    max_coverages.append( pooled_coverages.max() )
    
pylab.ylim([0,1000])
pylab.xlim([0.5,current_loc+0.5])
pylab.legend(loc='center right',frameon=False,fontsize=6)

pylab.xlabel('Position of site (marker genes are concatenated)') 
pylab.ylabel('Coverage') 

pylab.savefig('%s/%s_marker_coverages.pdf' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight',transparent=True)
pylab.savefig('%s/%s_marker_coverages.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight',dpi=300)


sys.exit(0)

samples = numpy.array(marker_gene_coverages.keys())

pooled_marker_gene_coverages = []
for sample in samples:
    coverages = []
    for gene_name in marker_gene_coverages[sample].keys():
        locs,depths = marker_gene_coverages[sample][gene_name]
        coverages.extend(depths)
    
    pooled_marker_gene_coverages.append(coverages)
    
pooled_marker_gene_coverages = numpy.array(pooled_marker_gene_coverages)


sys.exit(0)


print marker_gene_coverages.shape


median_marker_coverages = numpy.median(marker_gene_coverages,axis=0)
#mean_marker_coverages = marker_gene_coverages.mean(axis=0)
mean_marker_coverages = (marker_gene_coverages*(marker_gene_coverages>=1)).sum(axis=0)/((marker_gene_coverages>=1).sum(axis=0))
max_marker_coverages = marker_gene_coverages.max(axis=0)
num_zero = (marker_gene_coverages<1).sum(axis=0)

high_coverage_idxs = numpy.nonzero(median_marker_coverages>=100)[0]

print marker_genes
print pooled_marker_coverages[high_coverage_idxs[0]], ":", marker_gene_coverages[:,high_coverage_idxs[0]]
print pooled_marker_coverages[high_coverage_idxs[1]], ":", marker_gene_coverages[:,high_coverage_idxs[1]]

print marker_gene_coverages[1,:]

print num_zero

# Load genomic coverage distributions
sample_coverage_histograms, samples = parse_midas_data.parse_coverage_distribution(species_name)
median_coverages = numpy.array([stats_utils.calculate_nonzero_median_from_histogram(sample_coverage_histogram) for sample_coverage_histogram in sample_coverage_histograms])
sample_coverage_map = {samples[i]: median_coverages[i] for i in xrange(0,len(samples))}


pylab.figure(1)
for i in xrange(0,len(sample_coverage_histograms)):
    sample = samples[i]
    
    if sample not in desired_samples:
        continue
    
    sample_coverage_histogram = sample_coverage_histograms[i]
    x0 = median_coverages[i]
    
    if x0<20:
        continue
    
    #x0 = marker_coverages[i]
    xs, CDFs = stats_utils.calculate_unnormalized_CDF_from_histogram(sample_coverage_histogram)
    pylab.plot(xs, CDFs[-1]-CDFs, '-')

pylab.semilogx([1],[1])
pylab.xlabel('Coverage, D')
pylab.ylabel('Fraction sites with coverage >= D')
#pylab.xlim([1e-01,1e01])
pylab.savefig('%s/%s_genomic_coverage_distribution.pdf' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight',transparent=True)

pylab.figure(2)

median_coverages.sort()

median_coverage_xs, median_coverage_survivals = stats_utils.calculate_unnormalized_survival_from_vector(median_coverages, min_x=0.1, max_x=10000)


# Load gene coverage information for species_name
sys.stderr.write("Loading pangenome data for %s...\n" % species_name)
gene_samples, gene_names, gene_presence_matrix, gene_depth_matrix, marker_coverages, gene_reads_matrix = parse_midas_data.parse_pangenome_data(species_name)
sys.stderr.write("Done!\n")

marker_coverages = numpy.clip(marker_coverages, 2e-01,1e04)

marker_coverages.sort()

marker_coverage_xs, marker_coverage_survivals = stats_utils.calculate_unnormalized_survival_from_vector(marker_coverages, min_x=0.1, max_x=10000)

median_marker_coverage_xs, median_marker_coverage_survivals = stats_utils.calculate_unnormalized_survival_from_vector(median_marker_coverages, min_x=0.1, max_x=10000)
mean_marker_coverage_xs, mean_marker_coverage_survivals = stats_utils.calculate_unnormalized_survival_from_vector(mean_marker_coverages, min_x=0.1, max_x=10000)
pooled_marker_coverage_xs, pooled_marker_coverage_survivals = stats_utils.calculate_unnormalized_survival_from_vector(pooled_marker_coverages, min_x=0.1, max_x=10000)

max_marker_coverage_xs, max_marker_coverage_survivals = stats_utils.calculate_unnormalized_survival_from_vector(max_marker_coverages, min_x=0.1, max_x=10000)
median_metaphlan2_coverage_xs, median_metaphlan2_coverage_survivals = stats_utils.calculate_unnormalized_survival_from_vector(median_metaphlan2_coverages, min_x=0.1, max_x=10000)
mean_metaphlan2_coverage_xs, mean_metaphlan2_coverage_survivals = stats_utils.calculate_unnormalized_survival_from_vector(mean_metaphlan2_coverages, min_x=0.1, max_x=10000)

species_marker_coverage_xs, species_marker_coverage_survivals = stats_utils.calculate_unnormalized_survival_from_vector(species_marker_coverages, min_x=0.1, max_x=10000)


print len(marker_coverages), len(median_coverages)

pylab.title(species_name)
pylab.step(median_coverage_xs, median_coverage_survivals,label='Median(>0) reference')
pylab.step(median_metaphlan2_coverage_xs, median_metaphlan2_coverage_survivals,label='Median metaphlan2')
pylab.step(mean_metaphlan2_coverage_xs, mean_metaphlan2_coverage_survivals,label='Mean(>=1) metaphlan2')
pylab.step(pooled_marker_coverage_xs, pooled_marker_coverage_survivals,label='Median(>0) marker (snps/)')
pylab.step(mean_marker_coverage_xs, mean_marker_coverage_survivals,label='Mean (>=1) marker (snps/)')
pylab.step(marker_coverage_xs, marker_coverage_survivals,label='Marker (genes/)')
pylab.step(species_marker_coverage_xs, species_marker_coverage_survivals,label='Marker (species/)')
pylab.legend(frameon=False,fontsize=7)
pylab.loglog([1],[1e-01])
pylab.ylim([1,1e03])
pylab.legend(frameon=False,fontsize=7)
pylab.xlabel('Coverage, D')
pylab.ylabel('Num samples >= D')


pylab.savefig('%s/%s_median_coverages.pdf' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight',transparent=True)
pylab.savefig('%s/%s_median_coverages.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight',dpi=300)


print " ".join(["%0.1f" % m for m in median_coverages[-10:]])
for sample in sample_coverage_map.keys():
    if sample_coverage_map[sample] > 400:
        print sample, sample_coverage_map[sample]
        
print "---"
print '700037453', sample_coverage_map['700037453']
print '700037539', sample_coverage_map['700037539']
