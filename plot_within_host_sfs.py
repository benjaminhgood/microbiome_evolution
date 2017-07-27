import matplotlib  
matplotlib.use('Agg') 
import parse_midas_data
import pylab
import sys
import numpy
from numpy.random import normal
#from calculate_pi_matrix import calculate_self_pis
import diversity_utils
import stats_utils
import matplotlib.colors as colors
import matplotlib.cm as cmx
from math import log10,ceil
import matplotlib as mpl

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy.random import randint


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


min_coverage = 20

# Load subject and sample metadata
sys.stderr.write("Loading HMP metadata...\n")
subject_sample_map = parse_midas_data.parse_subject_sample_map()
sys.stderr.write("Done!\n")
       
# Load genomic coverage distributions
sample_coverage_histograms, samples = parse_midas_data.parse_coverage_distribution(species_name)
median_coverages = numpy.array([stats_utils.calculate_nonzero_median_from_histogram(sample_coverage_histogram) for sample_coverage_histogram in sample_coverage_histograms])
sample_coverage_map = {samples[i]: median_coverages[i] for i in xrange(0,len(samples))}
samples = numpy.array(samples)

median_coverages = numpy.array([sample_coverage_map[samples[i]] for i in xrange(0,len(samples))])

# Only plot samples above a certain depth threshold
desired_samples = samples[(median_coverages>=min_coverage)]
desired_median_coverages = numpy.array([sample_coverage_map[sample] for sample in desired_samples])

# Load SNP information for species_name
sys.stderr.write("Loading SFSs for %s...\t" % species_name)
samples, sfs_map = parse_midas_data.parse_within_sample_sfs(species_name, allowed_variant_types=set(['4D'])) 
sys.stderr.write("Done!\n")

# Set up binned SFS
sys.stderr.write("Calculating binned SFSs...\t")
binss = [numpy.arange(0,long(sample_coverage_map[sample])+1)*1.0/long(sample_coverage_map[sample]) for sample in desired_samples]
for i in xrange(0,len(binss)):
    binss[i][0] = -0.01
    binss[i][-1] = 1.01
        
bin_locationss = [binss[i][1:]-(binss[i][2]-binss[i][1])/2 for i in xrange(0,len(desired_samples))]        
countss = [numpy.zeros_like(bin_locationss[i]) for i in xrange(0,len(desired_samples))]

avg_pis = []
# Populate binned SFSs
for i in xrange(0,len(desired_samples)):

    freqs = []
    freq_counts = []
    
    avg_pis.append( diversity_utils.calculate_pi_from_sfs_map(sfs_map[desired_samples[i]]) )
    
    for key in sfs_map[desired_samples[i]].keys():
        D,A = key
        n = sfs_map[desired_samples[i]][key][0]
        
        is_polymorphic_truong = diversity_utils.is_polymorphic_truong(A,D)
        
        #if A==1:
        #    A=0
        #if A==(D-1):
        #    A=D
        
        f = A*1.0/D
        f = min([f,1-f])
        
        if is_polymorphic_truong:
            freqs.append(f)
            freq_counts.append(n)
        
        
    freqs = numpy.array(freqs)
    freq_counts = numpy.array(freq_counts)
    
    bin_idxs = numpy.digitize(freqs,bins=binss[i])
    
    for bin_idx, n in zip(bin_idxs, freq_counts):
        countss[i][bin_idx-1] += n


avg_pis = numpy.array(avg_pis)
           
sys.stderr.write("Done!\n")   

sys.stderr.write("Calculating smoothed SFSs...\n")
fss = []
pfss = []
between_polymorphism_rates = []
within_polymorphism_rates = []
ratios = []
raw_between_polymorphism_rates = []
raw_within_polymorphism_rates = []
bayes_within_polymorphism_rates = []


for i in xrange(0,len(desired_samples)):
    sys.stderr.write("%d\n" % i)
    fs, pfs, p_intermediate, p_poly = diversity_utils.calculate_smoothed_sfs(sfs_map[desired_samples[i]])
    fss.append(fs)
    pfss.append(pfs)
    
    bayes_within_polymorphism_rates.append( p_intermediate )
    
    between_polymorphism_rates.append( pfs[fs>0.8].sum() )
    within_polymorphism_rates.append( pfs[(fs<=0.8)*(fs>=0.2)].sum() )
    
    raw_within_polymorphism_rate, raw_between_polymorphism_rate = diversity_utils.calculate_polymorphism_rates_from_sfs_map(sfs_map[desired_samples[i]])
    raw_within_polymorphism_rates.append(raw_within_polymorphism_rate)
    raw_between_polymorphism_rates.append(raw_between_polymorphism_rate)
    
    ratios.append( raw_within_polymorphism_rate/raw_between_polymorphism_rate )
    #ratios.append( within_polymorphism_rates[-1]/between_polymorphism_rates[-1] )
    
sys.stderr.write("Done!\n")


ratios = numpy.array(ratios)
ratios = numpy.clip(ratios,1e-3,1e3)
between_polymorphism_rates = numpy.array(between_polymorphism_rates)
between_polymorphism_rates = numpy.clip(between_polymorphism_rates,3e-07,2)

within_polymorphism_rates = numpy.array(within_polymorphism_rates)
within_polymorphism_rates = numpy.clip(within_polymorphism_rates,3e-07,2)

raw_within_polymorphism_rates = numpy.array(raw_within_polymorphism_rates)
raw_within_polymorphism_rates = numpy.clip(raw_within_polymorphism_rates,3e-07,2)

raw_between_polymorphism_rates = numpy.array(raw_between_polymorphism_rates)
raw_between_polymorphism_rates = numpy.clip(raw_between_polymorphism_rates,3e-07,2)


# Now sort based on avg pi
sorted_idxs = numpy.arange(0,len(avg_pis))

sorted_avg_pis, sorted_idxs = (numpy.array(x) for x in zip(*sorted(zip(avg_pis,sorted_idxs),reverse=True)))

# This figure has everything on the same one
pylab.figure(1)
pylab.xlabel('Minor allele frequency')
pylab.ylabel('SFS')
pylab.xlim([0,0.5])
pylab.ylim([3e-06,3e-02])
pylab.title(species_name)
   
# This figure spreads them all out
figure_width=3.42
figure_height=2*len(avg_pis)
pylab.figure(2,figsize=(figure_width,figure_height))
fig = pylab.gcf()
# make three panels panels
outer_grid  = gridspec.GridSpec(len(avg_pis), 1, height_ratios=([1]*len(avg_pis)), hspace=0.1)


# This figure spreads them all out
pylab.figure(3,figsize=(figure_width,figure_height))
continuous_fig = pylab.gcf()
# make three panels panels
continuous_outer_grid  = gridspec.GridSpec(len(avg_pis), 1, height_ratios=([1]*len(avg_pis)), hspace=0.1)


sys.stderr.write("Plotting SFSs...\t")

for i in xrange(0,len(avg_pis)):

    pi = sorted_avg_pis[i]
    D = desired_median_coverages[sorted_idxs[i]]
    counts = numpy.array(countss[sorted_idxs[i]])
    bin_locations = bin_locationss[sorted_idxs[i]]
    
    sfs = counts/counts.sum()
    
    # First do figure 1
    colorVal = pi_scalarMap.to_rgba(numpy.log10(pi))
    
    pylab.figure(1)
    pylab.semilogy(bin_locations+normal(0,1)*(bin_locations[1]-bin_locations[0])*0.1, sfs,'.-',alpha=0.5, color=colorVal)
    
    
    # Then do figure 2
    pylab.figure(2)
    
    axis = plt.Subplot(fig, outer_grid[i])
    fig.add_subplot(axis)
    #axis.set_xlim([0,0.5])
    #axis.set_xticks(numpy.arange(0,11)*0.05)
    
    axis.set_xticks(numpy.arange(0,21)*0.05)
    axis.set_xlim([0.05,0.95])
    
    
    axis.set_xticklabels([])
    #axis.set_ylim([0,counts.max()+1])
    
    bin_idxs = (bin_locations>0.05)*(bin_locations<0.95)
    
    axis.plot(bin_locations[bin_idxs], counts[bin_idxs], label=('%d: pi=%g, D=%g' % (i, pi, D)))
    axis.legend(frameon=False,fontsize=7)
    
    # Then do figure 3
    pylab.figure(3)
    
    fs = fss[sorted_idxs[i]]
    pfs = pfss[sorted_idxs[i]]
    within_polymorphism_rate = within_polymorphism_rates[sorted_idxs[i]]
    between_polymorphism_rate = between_polymorphism_rates[sorted_idxs[i]]
    bayes_within_polymorphism_rate = bayes_within_polymorphism_rates[sorted_idxs[i]]
    raw_within_polymorphism_rate = raw_within_polymorphism_rates[sorted_idxs[i]]
    
    other_ratio = ratios[sorted_idxs[i]]
    
    between_polymorphism_line = numpy.ones_like(fs)*between_polymorphism_rate
    between_polymorphism_line /= ((fs>=0.1)*(fs<=0.9)).sum()
    
    within_polymorphism_line = numpy.ones_like(fs)*within_polymorphism_rate
    within_polymorphism_line /= ((fs>=0.1)*(fs<=0.9)).sum()
    
    
    axis = plt.Subplot(continuous_fig, continuous_outer_grid[i])
    continuous_fig.add_subplot(axis)
    #axis.set_xlim([0,0.5])
    #axis.set_xticks(numpy.arange(0,11)*0.05)
    
    axis.set_xticks(numpy.arange(0,21)*0.05)
    axis.set_xlim([0.05,0.95])
    
    if i<(len(desired_samples)-1):
        axis.set_xticklabels([])
    
    f_idxs = (fs>0.05)*(fs<0.95)
    
    #axis.set_ylim([0,pfs[f_idxs].max()])
    
    ratio = within_polymorphism_rate/between_polymorphism_rate
    
    print i, desired_samples[sorted_idxs[i]]
    print pi, ratio, other_ratio, raw_within_polymorphism_rate, bayes_within_polymorphism_rate
    
        
    axis.plot(fs[f_idxs], pfs[f_idxs], 'b-', label=('%d: pi=%g, r=%g, D=%g' % (i, pi, ratio, D)))
    axis.plot(fs[f_idxs], between_polymorphism_line[f_idxs],'r:')
    axis.plot(fs[f_idxs], within_polymorphism_line[f_idxs],'b:')
    
    axis.legend(frameon=False,fontsize=7)
    
     
sys.stderr.write("Done!\n")
   
pylab.figure(1)   
m = pylab.scatter([-1],[1],c=[-1], vmin=vmin, vmax=vmax, cmap=cmap,    marker='^')
    
fig = pylab.gcf()
cax = fig.add_axes([0.95, 0.1, 0.02, 0.80])
cbar = fig.colorbar(m,cax=cax,orientation='vertical',ticks=[-4,-3,-2])
    
cbar.set_ticklabels(['$10^{-4}$','$10^{-3}$','$10^{-2}$'])
#cl = pylab.getp(cbar.ax, 'ymajorticklabels')
#pylab.setp(cl, fontsize=9) 
fig.text(0.947,0.035,'$\\pi_s$',fontsize=12)

pylab.figure(4)
pylab.loglog(avg_pis, ratios, 'k.')
pylab.savefig('%s/%s_ratio_vs_pi.pdf' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight')

pylab.figure(5,figsize=(4,3))
pylab.loglog([1e-07,1e-01],[1e-07,1e-01],'k:')
pylab.fill_between([1e-07,1e-01],[1e-08,1e-08],[1e-08,1e-02],color='0.8')
pylab.loglog(raw_between_polymorphism_rates, raw_within_polymorphism_rates, 'k.',alpha=0.5)
pylab.ylim([1e-07,1e-01])
pylab.xlim([1e-07,1e-01]) 
pylab.ylabel('Intermediate freq polymorphism rate')
pylab.xlabel('High freq polymorphism rate')
pylab.savefig('%s/%s_within_vs_between.pdf' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight')

pylab.figure(6,figsize=(4,3))
pylab.loglog([1e-07,1e-01],[1e-07,1e-01],'k:')
pylab.loglog(raw_within_polymorphism_rates, bayes_within_polymorphism_rates, 'k.',alpha=0.5)
pylab.ylim([1e-07,1e-01])
pylab.xlim([1e-07,1e-01]) 
pylab.xlabel('Raw')
pylab.ylabel('Bayes')
pylab.savefig('%s/%s_within_raw_vs_em.pdf' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight')

sys.stderr.write("Saving figure...\t")
    
pylab.savefig('%s/%s_within_person_sfs.pdf' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight')
pylab.savefig('%s/%s_within_person_sfs.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight')

pylab.figure(2)
pylab.savefig('%s/%s_within_person_sfs_separate.pdf' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight')
pylab.savefig('%s/%s_within_person_sfs_separate.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight')

pylab.figure(3)
pylab.savefig('%s/%s_continuous_within_person_sfs_separate.pdf' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight')
pylab.savefig('%s/%s_continuous_within_person_sfs_separate.png' % (parse_midas_data.analysis_directory,species_name),bbox_inches='tight')
  

  
sys.stderr.write("Done!\n")