import matplotlib  
matplotlib.use('Agg') 
import config
import parse_midas_data
import parse_HMP_data
import os.path
import pylab
import sys
import numpy

import species_phylogeny_utils
import diversity_utils
import gene_diversity_utils
import calculate_substitution_rates
import clade_utils

import stats_utils
import matplotlib.colors as colors
import matplotlib.cm as cmx
from math import log10,ceil
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy.random import randint, shuffle

from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster.hierarchy import cophenet
from scipy.cluster.hierarchy import fcluster

from scipy.stats import gaussian_kde

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

num_bootstraps = 10000

clade_Fst = {}
Fst = {}
clade_country_likelihood = {}

divergence_matrices = {}
good_species_list = parse_midas_data.parse_good_species_list()

if debug:
    good_species_list = good_species_list[0:2]

# Load subject and sample metadata
sys.stderr.write("Loading sample metadata...\n")
subject_sample_map = parse_HMP_data.parse_subject_sample_map()
sample_country_map = parse_HMP_data.parse_sample_country_map()
sample_phenotype_map = parse_HMP_data.parse_sample_phenotype_map()
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
    
    # Load manually annotated clades
    clade_sets = clade_utils.load_manual_clades(species_name)

    if len(clade_sets)==0:
        continue
    
    clade_idxss = clade_utils.calculate_clade_idxs_from_clade_sets(snp_samples, clade_sets)
        
        
    nonsingleton_clade_idxss = []
    for clade_idxs in clade_idxss:
        if clade_idxs.sum() > 1:
            nonsingleton_clade_idxss.append(clade_idxs)
           
    # Want at least two clades!
    if len(nonsingleton_clade_idxss)>=2:
        all_nonsingleton_clade_idxs = numpy.array([False for x in nonsingleton_clade_idxss[0]])
        for clade_idxs in nonsingleton_clade_idxss:
            all_nonsingleton_clade_idxs = numpy.logical_or(all_nonsingleton_clade_idxs, clade_idxs)
    else:
        all_nonsingleton_clade_idxs = []    
    
    # First do Fst with clades!
    if len(nonsingleton_clade_idxss)>=2:
        # Calculate Fst between manually defined clades
        
        clade_substitution_matrices = []
        clade_ones_matrices = []
        clade_pair_idxss = []
        for clade_idxs in nonsingleton_clade_idxss:
            
            clade_substitution_matrices.append( snp_substitution_matrix[numpy.ix_(numpy.nonzero(clade_idxs)[0], numpy.nonzero(clade_idxs)[0])] )
            
            clade_ones_matrices.append( numpy.ones_like(clade_substitution_matrices[-1]) )
            
            clade_pair_idxss.append(  numpy.triu_indices(clade_substitution_matrices[-1].shape[0], 1) )
         
        all_substitution_matrix =  snp_substitution_matrix[ numpy.ix_(numpy.nonzero(all_nonsingleton_clade_idxs)[0], numpy.nonzero(all_nonsingleton_clade_idxs)[0]) ]
        
        all_ones_matrix = numpy.ones_like(all_substitution_matrix)
        
        all_pair_idxs = numpy.triu_indices(all_substitution_matrix.shape[0], 1)     
    
        within_numerator = sum([clade_substitution_matrices[i][clade_pair_idxss[i]].sum() for i in xrange(0,len(clade_substitution_matrices))])
        
        within_denominator = sum([clade_ones_matrices[i][clade_pair_idxss[i]].sum() for i in xrange(0,len(clade_substitution_matrices))])

        within_rate = within_numerator*1.0/within_denominator
        
        between_rate = (all_substitution_matrix[all_pair_idxs].sum())/(all_ones_matrix[all_pair_idxs].sum())
        
        observed_fst = 1.0 - within_rate/between_rate
        
        clade_Fst[species_name] = (observed_fst, [])
    
    
    # Load idxs corresponding to US and China
    us_idxs = numpy.array([sample_country_map[sample_name]=='United States' for sample_name in snp_samples])
    china_idxs = numpy.logical_not(us_idxs)
    
    # Make sure there are enough of both countries to do comparison
    if (us_idxs.sum() >= 2) and (china_idxs.sum()>=2):
        
        # Calculate Fst between US and China
        
        us_substitution_matrix = snp_substitution_matrix[numpy.ix_(numpy.nonzero(us_idxs)[0], numpy.nonzero(us_idxs)[0])]
        china_substitution_matrix = snp_substitution_matrix[numpy.ix_(numpy.nonzero(china_idxs)[0], numpy.nonzero(china_idxs)[0])]
    
        us_ones_matrix = numpy.ones_like(us_substitution_matrix)
        china_ones_matrix = numpy.ones_like(china_substitution_matrix)
        all_ones_matrix = numpy.ones_like(snp_substitution_matrix)
    
        us_pair_idxs = numpy.triu_indices(us_substitution_matrix.shape[0], 1)
        china_pair_idxs = numpy.triu_indices(china_substitution_matrix.shape[0], 1)
        all_pair_idxs = numpy.triu_indices(snp_substitution_matrix.shape[0], 1)
        
        observed_fst = 1.0 - (us_substitution_matrix[us_pair_idxs].sum()+china_substitution_matrix[china_pair_idxs].sum())/(us_ones_matrix[us_pair_idxs].sum()+china_ones_matrix[china_pair_idxs].sum())*(all_ones_matrix[all_pair_idxs].sum())/(snp_substitution_matrix[all_pair_idxs].sum())
        
        bootstrapped_fsts = []
        for bootstrap_idx in xrange(0,num_bootstraps):
            
            bootstrapped_us_idxs = numpy.array(us_idxs, copy=True)
            shuffle(bootstrapped_us_idxs)
            bootstrapped_china_idxs = numpy.logical_not(bootstrapped_us_idxs)
            
            us_substitution_matrix = snp_substitution_matrix[numpy.ix_(numpy.nonzero(bootstrapped_us_idxs)[0], numpy.nonzero(bootstrapped_us_idxs)[0])]
            china_substitution_matrix = snp_substitution_matrix[numpy.ix_(numpy.nonzero(bootstrapped_china_idxs)[0], numpy.nonzero(bootstrapped_china_idxs)[0])]
    
            us_ones_matrix = numpy.ones_like(us_substitution_matrix)
            china_ones_matrix = numpy.ones_like(china_substitution_matrix)
            all_ones_matrix = numpy.ones_like(snp_substitution_matrix)
    
            us_pair_idxs = numpy.triu_indices(us_substitution_matrix.shape[0], 1)
            china_pair_idxs = numpy.triu_indices(china_substitution_matrix.shape[0], 1)
            all_pair_idxs = numpy.triu_indices(snp_substitution_matrix.shape[0], 1)
            
            bootstrapped_fst = 1.0 - (us_substitution_matrix[us_pair_idxs].sum()+china_substitution_matrix[china_pair_idxs].sum())/(us_ones_matrix[us_pair_idxs].sum()+china_ones_matrix[china_pair_idxs].sum())*(all_ones_matrix[all_pair_idxs].sum())/(snp_substitution_matrix[all_pair_idxs].sum())
            
            bootstrapped_fsts.append(bootstrapped_fst)
        
        bootstrapped_fsts = numpy.array(bootstrapped_fsts)
        
        Fst[species_name] = (observed_fst, bootstrapped_fsts)
    
        if len(nonsingleton_clade_idxss)<2:
            continue
            
        # Ok, let's get calculating...
        def calculate_LRT(clade_idxss, phenotype_idxs):
        
            ns = numpy.array([clade_idxs.sum() for clade_idxs in clade_idxss])
            n1s = numpy.array([phenotype_idxs[clade_idxs].sum() for clade_idxs in clade_idxss])
        
            n_tot = ns.sum()
            n1_tot = n1s.sum()
        
            pavg = n1_tot*1.0/n_tot
            ps = n1s*1.0/ns
            
            # prevent taking a log of zero
            ps = numpy.clip(ps,1e-05,1-1e-05)
            
            logit_changes = numpy.log(ps/(1-ps)*(1-pavg)/pavg)
            
            ps[logit_changes<1] = pavg
            
            return (n1s*numpy.log(ps/pavg)+(ns-n1s)*numpy.log((1-ps)/(1-pavg))).sum()
        
        observed_LRT = calculate_LRT(nonsingleton_clade_idxss, us_idxs)
        bootstrapped_LRTs = []
        for bootstrap_idx in xrange(0, num_bootstraps):
            
            
            bootstrapped_phenotype_idxs = numpy.array(us_idxs,copy=True)
            desired_phenotypes = list(bootstrapped_phenotype_idxs[all_nonsingleton_clade_idxs])
            shuffle(desired_phenotypes)
            bootstrapped_phenotype_idxs[all_nonsingleton_clade_idxs] = numpy.array(desired_phenotypes)
            #desired_phenotypes = bootstrapped_phenotype_idxs[all_nonsingleton_clade_
            #shuffle(bootstrapped_phenotype_idxs[all_nonsingleton_clade_idxs])
            bootstrapped_LRTs.append( calculate_LRT(nonsingleton_clade_idxss, bootstrapped_phenotype_idxs) )
        
        bootstrapped_LRTs = numpy.array(bootstrapped_LRTs)
        
        print observed_LRT, bootstrapped_LRTs.mean(), (bootstrapped_LRTs>=observed_LRT).mean()
        
        clade_country_likelihood[species_name] = (observed_LRT, bootstrapped_LRTs)
            
            
species_names = []

for species_name in species_phylogeny_utils.sort_phylogenetically(Fst.keys()):
    species_names.append(species_name)
    
# sort in descending order of sample size
# Sort by num haploids    
#sample_sizes, species_names = zip(*sorted(zip(sample_sizes, species_names),reverse=True))
    
sys.stderr.write("Postprocessing %d species...\n" % len(species_names))
        

####################################################
#
# Set up Figure (3 panels, arranged in 3x1 grid)
#
####################################################

haploid_color = '#08519c'

pylab.figure(1,figsize=(4,3))
fig = pylab.gcf()
# make three panels panels
outer_grid  = gridspec.GridSpec(3,1,hspace=0.2,height_ratios=[1,1,1])

clade_fst_axis = plt.Subplot(fig, outer_grid[0])
fig.add_subplot(clade_fst_axis)

clade_fst_axis.set_ylabel('Fst (clades)')
clade_fst_axis.set_ylim([-0.05,1.05])
clade_fst_axis.set_xlim([-1,len(species_names)])

xticks = numpy.arange(0,len(species_names))
#xticklabels = ["%s (%d)" % (species_names[i],sample_sizes[i]) for i in xrange(0,len(sample_sizes))]
xticklabels = ["%s" % (species_names[i]) for i in xrange(0,len(species_names))]

clade_fst_axis.set_xticks(xticks)
clade_fst_axis.set_xticklabels([])

clade_fst_axis.spines['top'].set_visible(False)
clade_fst_axis.spines['right'].set_visible(False)
clade_fst_axis.get_xaxis().tick_bottom()
clade_fst_axis.get_yaxis().tick_left()



fst_axis = plt.Subplot(fig, outer_grid[1])
fig.add_subplot(fst_axis)

fst_axis.set_ylabel('Fst (U.S./China)')
fst_axis.set_ylim([-0.05,0.2])
fst_axis.set_xlim([-1,len(species_names)])

xticks = numpy.arange(0,len(species_names))
#xticklabels = ["%s (%d)" % (species_names[i],sample_sizes[i]) for i in xrange(0,len(sample_sizes))]
xticklabels = ["%s" % (species_names[i]) for i in xrange(0,len(species_names))]

fst_axis.set_xticks(xticks)
fst_axis.set_xticklabels([])

fst_axis.spines['top'].set_visible(False)
fst_axis.spines['right'].set_visible(False)
fst_axis.get_xaxis().tick_bottom()
fst_axis.get_yaxis().tick_left()

lrt_axis = plt.Subplot(fig, outer_grid[2])
fig.add_subplot(lrt_axis)

lrt_axis.set_ylabel('LRT statistic, $\\Delta \\ell$')
lrt_axis.set_ylim([-1,20])
lrt_axis.set_xlim([-1,len(species_names)])

xticks = numpy.arange(0,len(species_names))
#xticklabels = ["%s (%d)" % (species_names[i],sample_sizes[i]) for i in xrange(0,len(sample_sizes))]
xticklabels = ["%s" % (species_names[i]) for i in xrange(0,len(species_names))]

lrt_axis.set_xticks(xticks)
lrt_axis.set_xticklabels(xticklabels, rotation='vertical',fontsize=4)

lrt_axis.spines['top'].set_visible(False)
lrt_axis.spines['right'].set_visible(False)
lrt_axis.get_xaxis().tick_bottom()
lrt_axis.get_yaxis().tick_left()


# Plot percentiles of divergence distribution
for species_idx in xrange(0,len(species_names)):

    species_name = species_names[species_idx]

    clade_fst_axis.plot([-1,len(species_names)+1],[0.2,0.2],'k:',linewidth=0.5)
    
    if species_name in clade_Fst:
        observed_fst, bootstrapped_fsts = clade_Fst[species_name]
    
        clade_fst_axis.plot([species_idx], [observed_fst],'r^',markersize=3,markeredgewidth=0) 
    
    observed_fst, bootstrapped_fsts = Fst[species_name]
    bootstrapped_fsts.sort()
    
    is_significant = ((bootstrapped_fsts>=observed_fst).mean()<0.05)
    
    theory_fsts = numpy.linspace(bootstrapped_fsts.min(), bootstrapped_fsts.max(),100)
    kernel = gaussian_kde(bootstrapped_fsts)
    theory_pdf = kernel(theory_fsts)
    theory_pdf = theory_pdf / theory_pdf.max() * 0.45
    
    fst_axis.fill_betweenx(theory_fsts, species_idx-theory_pdf, species_idx+theory_pdf,linewidth=0,facecolor='0.7')
    if is_significant:
        fst_axis.plot([species_idx],[observed_fst],'rs',markersize=3,markeredgewidth=0)
    else:
        fst_axis.plot([species_idx],[observed_fst],'ro',markersize=3,markeredgewidth=0,alpha=0.5)     
    
    if species_name in clade_country_likelihood:
    
        observed_LRT, bootstrapped_LRTs = clade_country_likelihood[species_name]
        bootstrapped_LRTs.sort()
    
        LRT_scale = 1 #numpy.sqrt(numpy.square(bootstrapped_LRTs).mean())
    
        is_significant = ((bootstrapped_LRTs>=observed_LRT).mean()<0.05)
    
        kernel = gaussian_kde(bootstrapped_LRTs)
    
        theory_LRTs = numpy.linspace(bootstrapped_LRTs.min(), bootstrapped_LRTs.max(),100)
        theory_pdf = kernel(theory_LRTs)
        theory_pdf = theory_pdf / theory_pdf.max() * 0.45
    
        lrt_axis.fill_betweenx(theory_LRTs/LRT_scale, species_idx-theory_pdf, species_idx+theory_pdf,linewidth=0,facecolor='0.7')
        if is_significant:
            lrt_axis.plot([species_idx],[observed_LRT/LRT_scale],'rs',markersize=3,markeredgewidth=0)
        else:
            lrt_axis.plot([species_idx],[observed_LRT/LRT_scale],'ro',markersize=3,markeredgewidth=0,alpha=0.5)
            
        
    

sys.stderr.write("Saving figure...\t")
fig.savefig('%s/supplemental_clade_city_correlation.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')
sys.stderr.write("Done!\n")

 