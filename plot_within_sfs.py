import matplotlib  
matplotlib.use('Agg') 
import parse_midas_data
import pylab
import sys
import numpy
from numpy.random import normal
from calculate_pi_matrix import calculate_self_pis
import diversity_utils
import stats_utils
import matplotlib.colors as colors
import matplotlib.cm as cmx
from math import log10,ceil
import matplotlib as mpl
mpl.rcParams['font.size'] = 8
mpl.rcParams['lines.linewidth'] = 1.0
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'

vmin = log10(1e-04)
vmax = log10(1e-02)
cmap='jet'
jet = cm = pylab.get_cmap(cmap) 
cNorm  = colors.Normalize(vmin=vmin, vmax=vmax)
pi_scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)

species_name=sys.argv[1]



# Load subject and sample metadata
sys.stderr.write("Loading HMP metadata...\n")
subject_sample_map = parse_midas_data.parse_subject_sample_map()
sys.stderr.write("Done!\n")    

pylab.figure()
pylab.xlabel('Minor allele frequency')
pylab.ylabel('SFS')
pylab.xlim([0,0.5])
pylab.ylim([3e-06,3e-02])
pylab.title(species_name)
    

# Load genomic coverage distributions
sample_coverage_histograms, samples = parse_midas_data.parse_coverage_distribution(species_name, combination_type="sample")
median_coverages = numpy.array([stats_utils.calculate_median_from_histogram(sample_coverage_histogram) for sample_coverage_histogram in sample_coverage_histograms])
    
sample_coverage_map = {samples[i]: median_coverages[i] for i in xrange(0,len(samples))}

# Load SNP information for species_name
sys.stderr.write("Loading %s...\n" % species_name)
samples, allele_counts_map, passed_sites_map = parse_midas_data.parse_snps(species_name, combination_type="sample", debug=False)
sys.stderr.write("Done!\n")
    
median_coverages = numpy.array([sample_coverage_map[samples[i]] for i in xrange(0,len(samples))])
    
    
# Calculate full matrix of synonymous pairwise differences
sys.stderr.write("Calculate synonymous pi matrix...\n")
pi_matrix_syn, avg_pi_matrix_syn = diversity_utils.calculate_pi_matrix(allele_counts_map, passed_sites_map, variant_type='4D')

sys.stderr.write("Calculate within person SFS")
sample_freqs, passed_sites = diversity_utils.calculate_sample_freqs( allele_counts_map, passed_sites_map, variant_type='4D')
    
sfss = []
pi_withins = []
bins = numpy.linspace(0.04,0.51,11)
bin_locations = bins[1:]-(bins[1]-bins[0])/2
    
for j in xrange(0,len(samples)):
    
    pi_within = pi_matrix_syn[j,j]
        
    counts,dummy = numpy.histogram(sample_freqs[j],bins)
        
    if counts.sum()<0.5:
        sfs = numpy.zeros_like(bin_locations)
    else:
        sfs = counts*1.0/(passed_sites[j])
        
    if (median_coverages[j]>=100):
        
        sfss.append(sfs)
        pi_withins.append(pi_within)
        
sys.stderr.write("Done!\n")
    
current_max = 0
    
for j in xrange(0,len(sfss)):
    
    if sfss[j].sum()==0:
        continue
    
    colorVal = pi_scalarMap.to_rgba(numpy.log10(pi_withins[j]))
    normalized_sfs = sfss[j]
    pylab.semilogy(bin_locations+normal(0,1)*(bin_locations[1]-bin_locations[0])*0.1, normalized_sfs,'.-',alpha=0.5, color=colorVal)
    
    
    
m = pylab.scatter([-1],[1],c=[-1], vmin=vmin, vmax=vmax, cmap=cmap,    marker='^')
    
fig = pylab.gcf()
cax = fig.add_axes([0.95, 0.1, 0.02, 0.80])
cbar = fig.colorbar(m,cax=cax,orientation='vertical',ticks=[-4,-3,-2])
    
cbar.set_ticklabels(['$10^{-4}$','$10^{-3}$','$10^{-2}$'])
#cl = pylab.getp(cbar.ax, 'ymajorticklabels')
#pylab.setp(cl, fontsize=9) 
fig.text(0.947,0.035,'$\\pi_s$',fontsize=12)
    
pylab.savefig('%s/%s_within_person_sfs.pdf' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight')
    
