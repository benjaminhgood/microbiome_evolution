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
import calculate_temporal_changes
import calculate_substitution_rates
import stats_utils
import sfs_utils
    
    
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
min_sample_size = 5

modification_difference_threshold = 100

# Must compete divergence matrix on the fly! 
            
# Load subject and sample metadata
sys.stderr.write("Loading sample metadata...\n")
subject_sample_map = parse_HMP_data.parse_subject_sample_map()
sample_order_map = parse_HMP_data.parse_sample_order_map()
sys.stderr.write("Done!\n")
    
good_species_list = parse_midas_data.parse_good_species_list()
if debug:
    good_species_list = good_species_list[0:2]

num_passed_species = 0

num_temporal_change_map = {}

total_snp_modification_map = {}
total_null_snp_modification_map = {}

total_gene_modification_map = {}
total_null_gene_modification_map = {}

total_syn_snps = 0
total_non_snps = 0
total_syn_opportunities = 0
total_non_opportunities = 0

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



    allowed_sample_set = set(snp_samples)

    if sample_size < min_sample_size:
        continue
    
    
    sys.stderr.write("Proceeding with %d longitudinal comparisons in %d samples!\n" % (sample_size, len(snp_samples)))
    
    sys.stderr.write("Loading SFSs for %s...\t" % species_name)
    dummy_samples, sfs_map = parse_midas_data.parse_within_sample_sfs(species_name, allowed_variant_types=set(['1D','2D','3D','4D'])) 
    sys.stderr.write("Done!\n")

    sys.stderr.write("Loading pre-computed substitution rates for %s...\n" % species_name)
    substitution_rate_map = calculate_substitution_rates.load_substitution_rate_map(species_name)
    sys.stderr.write("Calculating matrix...\n")
    dummy_samples, snp_difference_matrix, snp_opportunity_matrix =    calculate_substitution_rates.calculate_matrices_from_substitution_rate_map(substitution_rate_map, 'all', allowed_samples=snp_samples)
    
    dummy_samples, non_difference_matrix, non_opportunity_matrix =    calculate_substitution_rates.calculate_matrices_from_substitution_rate_map(substitution_rate_map, '1D', allowed_samples=snp_samples)
    dummy_samples, syn_difference_matrix, syn_opportunity_matrix =    calculate_substitution_rates.calculate_matrices_from_substitution_rate_map(substitution_rate_map, '4D', allowed_samples=snp_samples)
    snp_samples = dummy_samples
    
    
    snp_substitution_rate =     snp_difference_matrix*1.0/(snp_opportunity_matrix+(snp_opportunity_matrix==0))
    sys.stderr.write("Done!\n")

    sys.stderr.write("Loading pre-computed temporal changes for %s...\n" % species_name)
    temporal_change_map = calculate_temporal_changes.load_temporal_change_map(species_name)
    sys.stderr.write("Done!\n")

    

    total_snp_modification_map[species_name] = 0
    total_null_snp_modification_map[species_name] = 0
    total_gene_modification_map[species_name] = 0
    total_null_gene_modification_map[species_name] = 0

    temporal_changes = []
    
    same_sample_idxs, same_subject_idxs, diff_subject_idxs = parse_midas_data.calculate_ordered_subject_pairs(sample_order_map, snp_samples)
    
    for sample_pair_idx in xrange(0,len(same_subject_idxs[0])):
#    
        i = same_subject_idxs[0][sample_pair_idx]
        j = same_subject_idxs[1][sample_pair_idx]
    
        sample_i = snp_samples[i] 
        sample_j = snp_samples[j]
        
        if not ((sample_i in allowed_sample_set) and (sample_j in allowed_sample_set)):
            continue
        
        perr, mutations, reversions = calculate_temporal_changes.calculate_mutations_reversions_from_temporal_change_map(temporal_change_map, sample_i, sample_j)
        
        if perr>=0.5:
            
            # Calculate a more fine grained value!
        
            dfs = numpy.array([0.6,0.7,0.8,0.9])
            perrs = diversity_utils.calculate_fixation_error_rate(sfs_map, sample_i, sample_j,dfs=dfs) * snp_opportunity_matrix[i, j]
    
            if (perrs<0.5).any():
                # take most permissive one!
                perr_idx = numpy.nonzero(perrs<0.5)[0][0]
                df = dfs[perr_idx]
                perr = perrs[perr_idx]
            
                # recalculate stuff!    
                perr, mutations, reversions = calculate_temporal_changes.calculate_mutations_reversions_from_temporal_change_map(temporal_change_map, sample_i, sample_j,lower_threshold=(1-df)/2.0, upper_threshold=(1+df)/2.0)
                
            else:
                df = 2
                perr = 1
                mutations = None
                reversions = None
    
        if mutations==None or perr>=0.5:
            num_mutations = 0
            num_reversions = 0
            num_snp_changes = -1
        else:
            num_mutations = len(mutations)
            num_reversions = len(reversions)
            num_snp_changes = num_mutations+num_reversions
    
        
        gene_perr, gains, losses = calculate_temporal_changes.calculate_gains_losses_from_temporal_change_map(temporal_change_map, sample_i, sample_j)
    
        
        if (gains==None) or (gene_perr<-0.5) or (gene_perr>0.5):
            num_gains = 0
            num_losses = 0
            num_gene_changes = -1
        else:
            num_gains = len(gains)
            num_losses = len(losses)
            num_gene_changes = num_gains+num_losses
    
        
    
        if (num_snp_changes<modification_difference_threshold) and (num_snp_changes>-0.5):
            total_snp_modification_map[species_name] += num_snp_changes
            total_null_snp_modification_map[species_name] += perr
            
            total_syn_opportunities += syn_opportunity_matrix[i,j]
            total_non_opportunities += non_opportunity_matrix[i,j]
            
            for snp_change in mutations:
                if snp_change[3]=='4D':
                    total_syn_snps += 1
                elif snp_change[3]=='1D':
                    total_non_snps += 1
                else:
                    pass
        
            for snp_change in reversions:
                if snp_change[3]=='4D':
                    total_syn_snps += 1
                elif snp_change[3]=='1D':
                    total_non_snps += 1
                else:
                    pass
            
            if num_gene_changes > -0.5:
            
                
                total_gene_modification_map[species_name] += num_gene_changes
                total_null_gene_modification_map[species_name] += gene_perr   
        
    
        temporal_changes.append((sample_i, sample_j, num_snp_changes, num_gene_changes))        

    num_passed_species+=1
    
    if len(temporal_changes) > 0:
        num_temporal_change_map[species_name] = temporal_changes
        
    sys.stderr.write("Done with %s!\n" % species_name) 
    
sys.stderr.write("Done looping over species!\n")    

print total_syn_snps*1.0/total_syn_opportunities, total_syn_snps, total_syn_opportunities
print total_non_snps*1.0/total_non_opportunities, total_non_snps, total_non_opportunities

species_names = []
sample_sizes = []


for species_name in num_temporal_change_map.keys():
    species_names.append(species_name)
    sample_sizes.append( len(num_temporal_change_map[species_name]) )
    
# sort in descending order of sample size
# Sort by num haploids    
sample_sizes, species_names = zip(*sorted(zip(sample_sizes, species_names),reverse=True))
     
sys.stderr.write("Postprocessing %d species!\n" % len(species_names))


cmap_str = 'YlGnBu'
vmin = -2
vmax = 3
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
    temporal_changes = num_temporal_change_map[species_name]
    
    
    if total_snp_modification_map[species_name] > 0:
        snp_is_significant = ((total_null_snp_modification_map[species_name]*1.0/total_snp_modification_map[species_name]) < 0.1)
    else:
        snp_is_significant = 0
    
    if total_gene_modification_map[species_name] > 0:
        gene_is_significant = ((total_null_gene_modification_map[species_name]*1.0/total_gene_modification_map[species_name]) < 0.1)
    else:
        gene_is_significant = 0
    
    print species_name, total_snp_modification_map[species_name], total_null_snp_modification_map[species_name], snp_is_significant, total_gene_modification_map[species_name], total_null_gene_modification_map[species_name], gene_is_significant
    
     
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
        
        change_axis.fill_between([species_idx-0.3,species_idx+0.3], [idx,idx],[idx+1.05,idx+1.05],color=colorVal,linewidth=0)
    
        
        if snp_is_significant:
            
            change_axis.text(species_idx, len(snp_changes),'*',fontsize=4)
    
    for idx in xrange(0,len(gene_changes)):
        
        if gene_changes[idx]<0.5:
            colorVal='0.7'
        else:
            colorVal = scalarMap.to_rgba(log10(gene_changes[idx]))
        
        change_axis.fill_between([species_idx-0.3,species_idx+0.3], [-idx-1.05,-idx-1.05],[-idx,-idx],color=colorVal,linewidth=0)
        
        if gene_is_significant:
            
            change_axis.text(species_idx, -len(gene_changes)-3,'*',fontsize=4)
            

m = change_axis.scatter([200],[1],c=[0], vmin=0, vmax=vmax, cmap=cmap, marker='^')


cbar = fig.colorbar(m,cax=cax,orientation='vertical', ticks=[0,1,2,3])
cbar.set_ticklabels(['$1$','$10$','$10^{2}$','$10^{3}$'])
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

 