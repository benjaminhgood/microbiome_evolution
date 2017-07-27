import matplotlib  
matplotlib.use('Agg') 
import parse_midas_data
import pylab
import sys
import numpy
from numpy.random import normal
import diversity_utils
import gene_diversity_utils
import stats_utils
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

#############################################################################

# Minimum median coverage of sample to look at
min_coverage = 20

sys.stderr.write("Loading core genes...\n")
core_genes = parse_midas_data.load_core_genes(species_name)
sys.stderr.write("Done! %d core genes\n" % len(core_genes))

#################
# Load metadata #
#################

# Load subject and sample metadata
sys.stderr.write("Loading HMP metadata...\n")
subject_sample_map = parse_midas_data.parse_subject_sample_map()
sys.stderr.write("Done!\n")
    
# Load time metadata
subject_sample_time_map_all_samples = parse_midas_data.parse_subject_sample_time_map()

######################
# Load coverage data #
######################

# Load genomic coverage distributions
sample_coverage_histograms, samples = parse_midas_data.parse_coverage_distribution(species_name)
median_coverages = numpy.array([stats_utils.calculate_median_from_histogram(sample_coverage_histogram) for sample_coverage_histogram in sample_coverage_histograms])
sample_coverage_map = {samples[i]: median_coverages[i] for i in xrange(0,len(samples))}
    
# prune time meta data so that the highest coverage sample is retained for those subjects with >1 sample per time pt
subject_sample_time_map = parse_midas_data.prune_subject_sample_time_map(subject_sample_time_map_all_samples,sample_coverage_map)


###############################################################
# Compute Pi within patients to figure out which are haploid  #
###############################################################

# Load pi information for species_name
sys.stderr.write("Loading within-sample diversity for %s...\n" % species_name)
samples, total_pis, total_pi_opportunities = parse_midas_data.parse_within_sample_pi(species_name, allowed_variant_types=set(['4D']), allowed_genes=core_genes, debug=debug)
sys.stderr.write("Done!\n")
pis = total_pis/total_pi_opportunities

######################
# compute median cov #
######################

median_coverages = numpy.array([sample_coverage_map[samples[i]] for i in xrange(0,len(samples))])

##########################################################
# load SNP info
##########################################################

# note that this loads info for all samples. Later the desired samples are selected out. 
sys.stderr.write("Loading %s...\n" % species_name)
samples, allele_counts_map, passed_sites_map, final_line_number = parse_midas_data.parse_snps(species_name, debug)
sys.stderr.write("Done!\n")

###############################################################
# Indexes for SNP samples that have high coverage and low pis #
###############################################################

# Only plot samples above a certain depth threshold that are "haploids"
high_cov_samples_low_pis = samples[(median_coverages>=min_coverage)*(pis<=1e-03)]
desired_samples_low_pis=(median_coverages>=min_coverage)*(pis<=1e-03)

# get the time info for the high_cov_samples_low_pis -- this is for all visno combos and all samples irrespective of coverage and pis
time_pair_idxs, visno1, visno2, day = parse_midas_data.calculate_all_time_pairs(subject_sample_time_map, high_cov_samples_low_pis)


###################################################################
# create subplots to plot the 2D and 1D SFSs side by side
###################################################################                    
    
fig_annotation='low_pis'

# compute the 1D sfs
sys.stderr.write("Calculate within person SFS for low pis...\n")
sample_freqs, passed_sites = diversity_utils.calculate_sample_freqs(allele_counts_map, passed_sites_map, variant_type='4D',fold=False)

sfss= []
bins = numpy.linspace(0.04,0.95,21)
bin_locations = bins[1:]-(bins[1]-bins[0])/2

for j in xrange(0,len(samples)):    
    counts,dummy = numpy.histogram(sample_freqs[j],bins)
    if counts.sum()<0.5:
        sfs = numpy.zeros_like(bin_locations)
    else:
        sfs = counts*1.0/(passed_sites[j])
    sfss.append(sfs)
        
# select the desired samples:
sfss=numpy.asarray(sfss)
sfss=sfss[desired_samples_low_pis]

sys.stderr.write("Done!\n")

# compute 2D SFS, plot, along iwth 2D SFS
# these sfs are polarized based on teh consensus allele. 
# note that desired samples includes the same as what is in high_cov_samples
    
sample_freqs_2D, passed_sites_2D, joint_passed_sites_2D = diversity_utils.calculate_sample_freqs_2D(allele_counts_map, passed_sites_map, desired_samples_low_pis, variant_type='4D', fold=False)

xbins = numpy.linspace(0,1,21) 
ybins = numpy.linspace(0,1,21)

pylab.figure(figsize=(6, len(time_pair_idxs[0])*3))   
color=['b','r']
plot_no=1
for j in range(0, len(time_pair_idxs[0])): 
    
    #plot the 2D SFS
    pylab.subplot(len(time_pair_idxs[0]), 2, plot_no) 
    idx1=time_pair_idxs[0][j]
    idx2=time_pair_idxs[1][j]
    sample_name1=high_cov_samples_low_pis[idx1]
    sample_name2=high_cov_samples_low_pis[idx2]
    freqs_idx1=numpy.array(sample_freqs_2D[idx1])
    freqs_idx2=numpy.array(sample_freqs_2D[idx2])
    joint_passed_idx1=numpy.array(joint_passed_sites_2D[idx1])
    joint_passed_idx2=numpy.array(joint_passed_sites_2D[idx2])
    joint_passed_idx1_idx2=(joint_passed_idx1)*(joint_passed_idx2)
    joint_passed_idx1_idx2=numpy.where(joint_passed_idx1_idx2==True)
    counts, xbins, ybins = numpy.histogram2d(freqs_idx1[joint_passed_idx1_idx2], freqs_idx2[joint_passed_idx1_idx2], bins=(xbins, ybins))
    sfs_2D=counts*1.0/(passed_sites_2D[idx1,idx2])
    pylab.xlabel('time point 1')
    pylab.ylabel('time point 2')
    pylab.xlim([0,1])
    pylab.ylim([0,1])
    pylab.title(sample_name1+', '+sample_name2)
    im=pylab.imshow(sfs_2D.T,interpolation='nearest', origin='low',extent=[xbins[0], xbins[-1], ybins[0], ybins[-1]], norm=LogNorm(vmin=1e-5, vmax=1e-2), cmap='jet')
    pylab.colorbar(im)
    plot_no+=1

    # plot the 1D SFS
    pylab.subplot(len(time_pair_idxs[0]), 2, plot_no)
    pylab.xlabel('Minor allele frequency')
    pylab.ylabel('SFS')
    pylab.xlim([0,1])
    pylab.ylim([3e-06,3e-02])
    colNo=0
    for idx in [idx1,idx2]:
        if sfss[idx].sum()!=0: 
            normalized_sfs = sfss[idx]
            pylab.semilogy(bin_locations+normal(0,1)*(bin_locations[1]-bin_locations[0])*0.1, normalized_sfs,'.-',alpha=0.5, color=color[colNo])
        colNo +=1
    pylab.legend(['first time pt', 'second time pt'],'upper right',prop={'size':6})
    plot_no+=1

pylab.savefig('%s/%s_within_person_2D_sfs_time_pair_polarized_%s.png' % (parse_midas_data.analysis_directory,species_name, fig_annotation),bbox_inches='tight')

#######################
# repeat for high pis #
#######################

# Only plot samples above a certain depth threshold that are "haploids"
high_cov_samples_high_pis = samples[(median_coverages>=min_coverage)*(pis>1e-03)]
desired_samples_high_pis=(median_coverages>=min_coverage)*(pis>1e-03)

# get the time info for the high_cov_samples_low_pis -- this is for all visno combos and all samples irrespective of coverage and pis
time_pair_idxs, visno1, visno2, day = parse_midas_data.calculate_all_time_pairs(subject_sample_time_map, high_cov_samples_high_pis)


###################################################################
# create subplots to plot the 2D and 1D SFSs side by side
###################################################################                    


fig_annotation='high_pis'

# compute the 1D sfs
sys.stderr.write("Calculate within person SFS for high pis...\n")
sample_freqs, passed_sites = diversity_utils.calculate_sample_freqs(allele_counts_map, passed_sites_map, variant_type='4D',fold=False)

sfss= []
bins = numpy.linspace(0.04,0.95,21)
bin_locations = bins[1:]-(bins[1]-bins[0])/2

for j in xrange(0,len(samples)):    
    counts,dummy = numpy.histogram(sample_freqs[j],bins)
    if counts.sum()<0.5:
        sfs = numpy.zeros_like(bin_locations)
    else:
        sfs = counts*1.0/(passed_sites[j])
    sfss.append(sfs)
        
# select the desired samples:
sfss=numpy.asarray(sfss)
sfss=sfss[desired_samples_high_pis]

sys.stderr.write("Done!\n")

# compute 2D SFS, plot, along iwth 2D SFS
# these sfs are polarized based on teh consensus allele. 
# note that desired samples includes the same as what is in high_cov_samples
    
sample_freqs_2D, passed_sites_2D, joint_passed_sites_2D = diversity_utils.calculate_sample_freqs_2D(allele_counts_map, passed_sites_map, desired_samples_high_pis, variant_type='4D', fold=False)

xbins = numpy.linspace(0,1,21) 
ybins = numpy.linspace(0,1,21)

pylab.figure(figsize=(6, len(time_pair_idxs[0])*3))   
color=['b','r']
plot_no=1
for j in range(0, len(time_pair_idxs[0])): 
    
    #plot the 2D SFS
    pylab.subplot(len(time_pair_idxs[0]), 2, plot_no) 
    idx1=time_pair_idxs[0][j]
    idx2=time_pair_idxs[1][j]
    sample_name1=high_cov_samples_high_pis[idx1]
    sample_name2=high_cov_samples_high_pis[idx2]
    freqs_idx1=numpy.array(sample_freqs_2D[idx1])
    freqs_idx2=numpy.array(sample_freqs_2D[idx2])
    joint_passed_idx1=numpy.array(joint_passed_sites_2D[idx1])
    joint_passed_idx2=numpy.array(joint_passed_sites_2D[idx2])
    joint_passed_idx1_idx2=(joint_passed_idx1)*(joint_passed_idx2)
    joint_passed_idx1_idx2=numpy.where(joint_passed_idx1_idx2==True)
    counts, xbins, ybins = numpy.histogram2d(freqs_idx1[joint_passed_idx1_idx2], freqs_idx2[joint_passed_idx1_idx2], bins=(xbins, ybins))
    sfs_2D=counts*1.0/(passed_sites_2D[idx1,idx2])
    pylab.xlabel('time point 1')
    pylab.ylabel('time point 2')
    pylab.xlim([0,1])
    pylab.ylim([0,1])
    pylab.title(sample_name1+', '+sample_name2)
    im=pylab.imshow(sfs_2D.T,interpolation='nearest', origin='low',extent=[xbins[0], xbins[-1], ybins[0], ybins[-1]], norm=LogNorm(vmin=1e-5, vmax=1e-2), cmap='jet')
    pylab.colorbar(im)
    plot_no+=1

    # plot the 1D SFS
    pylab.subplot(len(time_pair_idxs[0]), 2, plot_no)
    pylab.xlabel('Minor allele frequency')
    pylab.ylabel('SFS')
    pylab.xlim([0,1])
    pylab.ylim([3e-06,3e-02])
    colNo=0
    for idx in [idx1,idx2]:
        if sfss[idx].sum()!=0: 
            normalized_sfs = sfss[idx]
            pylab.semilogy(bin_locations+normal(0,1)*(bin_locations[1]-bin_locations[0])*0.1, normalized_sfs,'.-',alpha=0.5, color=color[colNo])
        colNo +=1
    pylab.legend(['first time pt', 'second time pt'],'upper right',prop={'size':6})
    plot_no+=1

pylab.savefig('%s/%s_within_person_2D_sfs_time_pair_polarized_%s.png' % (parse_midas_data.analysis_directory,species_name, fig_annotation),bbox_inches='tight')

#######################
# repeat for any pis #
#######################

# Only plot samples above a certain depth threshold that are "haploids"
high_cov_samples_any_pis = samples[(median_coverages>=min_coverage)]
desired_samples_any_pis=(median_coverages>=min_coverage)

# get the time info for the high_cov_samples_low_pis -- this is for all visno combos and all samples irrespective of coverage and pis
time_pair_idxs, visno1, visno2, day = parse_midas_data.calculate_all_time_pairs(subject_sample_time_map, high_cov_samples_any_pis)



###################################################################
# create subplots to plot the 2D and 1D SFSs side by side
###################################################################                    


fig_annotation='any_pis'

# compute the 1D sfs
sys.stderr.write("Calculate within person SFS for any pis...\n")
sample_freqs, passed_sites = diversity_utils.calculate_sample_freqs(allele_counts_map, passed_sites_map, variant_type='4D',fold=False)

sfss= []
bins = numpy.linspace(0.04,0.95,21)
bin_locations = bins[1:]-(bins[1]-bins[0])/2

for j in xrange(0,len(samples)):    
    counts,dummy = numpy.histogram(sample_freqs[j],bins)
    if counts.sum()<0.5:
        sfs = numpy.zeros_like(bin_locations)
    else:
        sfs = counts*1.0/(passed_sites[j])
    sfss.append(sfs)
        
# select the desired samples:
sfss=numpy.asarray(sfss)
sfss=sfss[desired_samples_any_pis]

sys.stderr.write("Done!\n")

# compute 2D SFS, plot, along iwth 2D SFS
# these sfs are polarized based on teh consensus allele. 
# note that desired samples includes the same as what is in high_cov_samples
    
sample_freqs_2D, passed_sites_2D, joint_passed_sites_2D = diversity_utils.calculate_sample_freqs_2D(allele_counts_map, passed_sites_map, desired_samples_any_pis, variant_type='4D', fold=False)

xbins = numpy.linspace(0,1,21) 
ybins = numpy.linspace(0,1,21)

pylab.figure(figsize=(6, len(time_pair_idxs[0])*3))   
color=['b','r']
plot_no=1
for j in range(0, len(time_pair_idxs[0])): 
    
    #plot the 2D SFS
    pylab.subplot(len(time_pair_idxs[0]), 2, plot_no) 
    idx1=time_pair_idxs[0][j]
    idx2=time_pair_idxs[1][j]
    sample_name1=high_cov_samples_any_pis[idx1]
    sample_name2=high_cov_samples_any_pis[idx2]
    freqs_idx1=numpy.array(sample_freqs_2D[idx1])
    freqs_idx2=numpy.array(sample_freqs_2D[idx2])
    joint_passed_idx1=numpy.array(joint_passed_sites_2D[idx1])
    joint_passed_idx2=numpy.array(joint_passed_sites_2D[idx2])
    joint_passed_idx1_idx2=(joint_passed_idx1)*(joint_passed_idx2)
    joint_passed_idx1_idx2=numpy.where(joint_passed_idx1_idx2==True)
    counts, xbins, ybins = numpy.histogram2d(freqs_idx1[joint_passed_idx1_idx2], freqs_idx2[joint_passed_idx1_idx2], bins=(xbins, ybins))
    sfs_2D=counts*1.0/(passed_sites_2D[idx1,idx2])
    pylab.xlabel('time point 1')
    pylab.ylabel('time point 2')
    pylab.xlim([0,1])
    pylab.ylim([0,1])
    pylab.title(sample_name1+', '+sample_name2)
    im=pylab.imshow(sfs_2D.T,interpolation='nearest', origin='low',extent=[xbins[0], xbins[-1], ybins[0], ybins[-1]], norm=LogNorm(vmin=1e-5, vmax=1e-2), cmap='jet')
    pylab.colorbar(im)
    plot_no+=1

    # plot the 1D SFS
    pylab.subplot(len(time_pair_idxs[0]), 2, plot_no)
    pylab.xlabel('Minor allele frequency')
    pylab.ylabel('SFS')
    pylab.xlim([0,1])
    pylab.ylim([3e-06,3e-02])
    colNo=0
    for idx in [idx1,idx2]:
        if sfss[idx].sum()!=0: 
            normalized_sfs = sfss[idx]
            pylab.semilogy(bin_locations+normal(0,1)*(bin_locations[1]-bin_locations[0])*0.1, normalized_sfs,'.-',alpha=0.5, color=color[colNo])
        colNo +=1
    pylab.legend(['first time pt', 'second time pt'],'upper right',prop={'size':6})
    plot_no+=1

pylab.savefig('%s/%s_within_person_2D_sfs_time_pair_polarized_%s.png' % (parse_midas_data.analysis_directory,species_name, fig_annotation),bbox_inches='tight')


##############################


'''


#######
# old


#folded SFS
sys.stderr.write("Calculate within person SFS...\n")
sample_freqs, passed_sites = diversity_utils.calculate_sample_freqs( allele_counts_map, passed_sites_map, variant_type='4D')

sfss = []
#pi_withins = []
bins = numpy.linspace(0.04,0.51,11)
bin_locations = bins[1:]-(bins[1]-bins[0])/2
    
for j in xrange(0,len(samples)):
    
    #pi_within = pis[j]
        
    counts,dummy = numpy.histogram(sample_freqs[j],bins)
        
    if counts.sum()<0.5:
        sfs = numpy.zeros_like(bin_locations)
    else:
        sfs = counts*1.0/(passed_sites[j])
        
    sfss.append(sfs)
    #pi_withins.append(pi_within)
        

sys.stderr.write("Done!\n")

##########
# Plot:  #
##########

# Folded SFSs

# iterate through time pairs

color=['b','r']
for j in range(0, len(time_pair_idxs[0])):
    pylab.figure()
    pylab.xlabel('Minor allele frequency')
    pylab.ylabel('SFS')
    pylab.xlim([0,0.5])
    pylab.ylim([3e-06,3e-02])
    pylab.title(species_name)
    idx1=time_pair_idxs[0][j]
    idx2=time_pair_idxs[1][j]
    sample_name1=samples[idx1]
    sample_name2=samples[idx2]
    colNo=0
    for idx in [idx1,idx2]:
        if sfss[idx].sum()!=0: 
            normalized_sfs = sfss[idx]
            pylab.semilogy(bin_locations+normal(0,1)*(bin_locations[1]-bin_locations[0])*0.1, normalized_sfs,'.-',alpha=0.5, color=color[colNo])
        colNo +=1
    
    pylab.legend(['first time pt', 'second time pt'],'upper right',prop={'size':6})
    
    pylab.savefig('%s/%s_within_person_sfs_time_pair_folded_%s_%s_low_pis.png' % (parse_midas_data.analysis_directory,species_name, sample_name1, sample_name2),bbox_inches='tight')


#################
#polarized SFS  #
#################

sys.stderr.write("Calculate within person SFS...\n")
sample_freqs, passed_sites = diversity_utils.calculate_sample_freqs( allele_counts_map, passed_sites_map, variant_type='4D',fold=False)

sfss= []
bins = numpy.linspace(0.04,0.95,21)
bin_locations = bins[1:]-(bins[1]-bins[0])/2
    
for j in xrange(0,len(samples)):
    
    counts,dummy = numpy.histogram(sample_freqs[j],bins)
        
    if counts.sum()<0.5:
        sfs = numpy.zeros_like(bin_locations)
    else:
        sfs = counts*1.0/(passed_sites[j])
        
    sfss.append(sfs)
        

sys.stderr.write("Done!\n")



# plot Polarized SFSs

# iterate through time pairs

color=['b','r']
for j in range(0, len(time_pair_idxs[0])):
    pylab.figure()
    pylab.xlabel('Minor allele frequency')
    pylab.ylabel('SFS')
    pylab.xlim([0,1])
    pylab.ylim([3e-06,3e-02])
    pylab.title(species_name)
    idx1=time_pair_idxs[0][j]
    idx2=time_pair_idxs[1][j]
    sample_name1=samples[idx1]
    sample_name2=samples[idx2]
    colNo=0
    for idx in [idx1,idx2]:
        if sfss[idx].sum()!=0: 
            normalized_sfs = sfss[idx]
            pylab.semilogy(bin_locations+normal(0,1)*(bin_locations[1]-bin_locations[0])*0.1, normalized_sfs,'.-',alpha=0.5, color=color[colNo])
        colNo +=1
    
    pylab.legend(['first time pt', 'second time pt'],'upper right',prop={'size':6})
    
    pylab.savefig('%s/%s_within_person_sfs_time_pair_polarized_%s_%s_low_pis.png' % (parse_midas_data.analysis_directory,species_name, sample_name1, sample_name2),bbox_inches='tight')





###################################################################
# Plot 2d SFS for time pairs (polarized)
###################################################################
# get the time info for the snp_samples -- this is for all visno combos
# this time condition on coverage being high and pis being low. 

time_pair_idxs, visno1, visno2, day = parse_midas_data.calculate_all_time_pairs(subject_sample_time_map, high_cov_samples)


# these sfs are polarized based on teh consensus allele. 
# note that desired samples includes the same as what is in high_cov_samples
    
sample_freqs_2D, passed_sites_2D, joint_passed_sites_2D = diversity_utils.calculate_sample_freqs_2D(allele_counts_map, passed_sites_map, desired_samples, variant_type='4D', fold=False)


xbins = numpy.linspace(0,1,21) 
ybins = numpy.linspace(0,1,21)

for j in range(0, len(time_pair_idxs[0])):
    idx1=time_pair_idxs[0][j]
    idx2=time_pair_idxs[1][j]
    sample_name1=high_cov_samples[idx1]
    sample_name2=high_cov_samples[idx2]
    freqs_idx1=numpy.array(sample_freqs_2D[idx1])
    freqs_idx2=numpy.array(sample_freqs_2D[idx2])
    joint_passed_idx1=numpy.array(joint_passed_sites_2D[idx1])
    joint_passed_idx2=numpy.array(joint_passed_sites_2D[idx2])
    joint_passed_idx1_idx2=(joint_passed_idx1)*(joint_passed_idx2)
    joint_passed_idx1_idx2=numpy.where(joint_passed_idx1_idx2==True)
    counts, xbins, ybins = numpy.histogram2d(freqs_idx1[joint_passed_idx1_idx2], freqs_idx2[joint_passed_idx1_idx2], bins=(xbins, ybins))
    sfs_2D=counts*1.0/(passed_sites_2D[idx1,idx2])
    pylab.figure()
    pylab.xlabel('time point 1')
    pylab.ylabel('time point 2')
    pylab.xlim([0,1])
    pylab.ylim([0,1])
    pylab.title(species_name)
    im=pylab.imshow(sfs_2D.T,interpolation='nearest', origin='low',extent=[xbins[0], xbins[-1], ybins[0], ybins[-1]], norm=LogNorm(vmin=1e-5, vmax=1), cmap='jet')
    pylab.colorbar(im)
    pylab.savefig('%s/%s_within_person_2D_sfs_time_pair_polarized_%s_%s_low_pis.png' % (parse_midas_data.analysis_directory,species_name, sample_name1, sample_name2),bbox_inches='tight')



###################################################################
# Plot 2d SFS for time pairs (folded)
###################################################################

    
sample_freqs_2D, passed_sites_2D, joint_passed_sites_2D = diversity_utils.calculate_sample_freqs_2D(allele_counts_map, passed_sites_map, desired_samples, variant_type='4D')


xbins = numpy.linspace(0,1.1,21) 
ybins = numpy.linspace(0,1.1,21)

for j in range(0, len(time_pair_idxs[0])):
    idx1=time_pair_idxs[0][j]
    idx2=time_pair_idxs[1][j]
    sample_name1=high_cov_samples[idx1]
    sample_name2=high_cov_samples[idx2]
    freqs_idx1=numpy.array(sample_freqs_2D[idx1])
    freqs_idx2=numpy.array(sample_freqs_2D[idx2])
    joint_passed_idx1=numpy.array(joint_passed_sites_2D[idx1])
    joint_passed_idx2=numpy.array(joint_passed_sites_2D[idx2])
    joint_passed_idx1_idx2=(joint_passed_idx1)*(joint_passed_idx2)
    joint_passed_idx1_idx2=numpy.where(joint_passed_idx1_idx2==True)
    counts, xbins, ybins = numpy.histogram2d(freqs_idx1[joint_passed_idx1_idx2], freqs_idx2[joint_passed_idx1_idx2], bins=(xbins, ybins))
    sfs_2D=counts*1.0/(passed_sites_2D[idx1,idx2])
    pylab.figure()
    pylab.xlabel('time point 1')
    pylab.ylabel('time point 2')
    pylab.xlim([0,1])
    pylab.ylim([0,1])
    pylab.title(species_name)
    im=pylab.imshow(sfs_2D.T,interpolation='nearest', origin='low',extent=[xbins[0], xbins[-1], ybins[0], ybins[-1]], norm=LogNorm(vmin=1e-5, vmax=1e-2), cmap='jet')
    pylab.colorbar(im)
    pylab.savefig('%s/%s_within_person_2D_sfs_time_pair_folded_%s_%s_low_pis.png' % (parse_midas_data.analysis_directory,species_name, sample_name1, sample_name2),bbox_inches='tight')



'''
