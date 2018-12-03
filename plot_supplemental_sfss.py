import matplotlib  
matplotlib.use('Agg') 
import parse_midas_data
import sample_utils
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
from numpy.random import choice

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from numpy.random import randint,binomial


import config
import sfs_utils
import stats_utils
import figure_utils

fontsize = 6
mpl.rcParams['font.size'] = fontsize
mpl.rcParams['lines.linewidth'] = 1.0
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'


min_coverage = config.min_median_coverage    

haploid_color = '#08519c'
light_haploid_color = '#6699CC'
diploid_color = '#de2d26' #'#fb6a4a' #
transition_color = '#756bb1'

# First twelve species from Fig. 1F
good_species_1 = ['Bacteroides_vulgatus_57955',
'Bacteroides_uniformis_57318',
'Alistipes_putredinis_61533',
'Bacteroides_ovatus_58035',
'Eubacterium_rectale_56927',
'Bacteroides_stercoris_56735',
'Parabacteroides_merdae_56972',
'Bacteroides_xylanisolvens_57185',
'Ruminococcus_bromii_62047',
'Bacteroides_thetaiotaomicron_56941',
'Bacteroides_caccae_53434',
'Alistipes_onderdonkii_55464']

# Next twelve species from Fig. 1F
# Note: since we only have room to show 24 species in two full-page figs
#       and since this supplemental fig focuses on non-QP samples,
#       I removed species that had low numbers of non-QP samples in Fig. 1F
good_species_2 = ['Alistipes_shahii_62199',
'Parabacteroides_distasonis_56985',
#'Barnesiella_intestinihominis_62208',
'Oscillibacter_sp_60799',
'Faecalibacterium_prausnitzii_62201',
#'Akkermansia_muciniphila_55290',
#'Bacteroides_massiliensis_44749',
#'Bacteroides cellulosilyticus_58046',
#'Ruminococcus_bicirculans_59300',
'Prevotella_copri_61740',
'Faecalibacterium_prausnitzii_57453',
#'Dialister_invisus_61905',
#'Bacteroides_fragilis_54507',
#'Alistipes_finegoldii_56071',
'Eubacterium_eligens_61678',
'Faecalibacterium_prausnitzii_61481',
'Faecalibacterium_cf_62236',
#'Bacteroides_plebeius_61623',
'Roseburia_intestinalis_56239',
'Escherichia_coli_58110',
#'Phascolarctobacterium_sp_59817',
'Blautia_wexlerae_56130']
#'Bacteroidales_bacterium_58650']

good_species_lists = [good_species_1, good_species_2]


#could also do this one...
#species_name = Bacteroides_uniformis_57318

# Load subject and sample metadata
sys.stderr.write("Loading sample metadata...\n")
subject_sample_map = sample_utils.parse_subject_sample_map()
sample_order_map = sample_utils.parse_sample_order_map()
sys.stderr.write("Done!\n")

pylab.figure(3,figsize=(7,9))
polymorphism_fig = pylab.gcf() 
 
polymorphism_grid = gridspec.GridSpec(6, 4, height_ratios=[1]*6, width_ratios=[1]*4, hspace=0.25, wspace=0.25)

 
for list_idx in [0,1]:
    figure_idx = list_idx+1
    good_species_list = good_species_lists[list_idx]

    pylab.figure(figure_idx,figsize=(7,9))
    fig = pylab.gcf()    
    species_sfs_grid = gridspec.GridSpec(12, 1, height_ratios=[1]*12, hspace=0.25)

    for species_idx in xrange(0,len(good_species_list)):
        species_name = good_species_list[species_idx]

        species_name_items = species_name.split("_")
        species_name_label_multiline = "\n".join(species_name_items[:-1])
        species_name_label_singleline = " ".join(species_name_items[:-1])
        # Load SNP information for species_name
        sys.stderr.write("Loading SFSs for %s...\t" % species_name)
        samples, sfs_map = parse_midas_data.parse_within_sample_sfs(species_name,     allowed_variant_types=set(['4D'])) 
        sys.stderr.write("Done!\n")

        # Load genomic coverage distributions
        sample_coverage_histograms, samples = parse_midas_data.parse_coverage_distribution(species_name)
        median_coverages = numpy.array([stats_utils.calculate_nonzero_median_from_histogram(sample_coverage_histogram) for sample_coverage_histogram in sample_coverage_histograms])
        sample_coverage_map = {samples[i]: median_coverages[i] for i in xrange(0,len(samples))}
        samples = numpy.array(samples)

        highcoverage_samples = set(diversity_utils.calculate_highcoverage_samples(species_name))
        qp_samples = set(diversity_utils.calculate_haploid_samples(species_name))
        
        non_qp_samples = list(highcoverage_samples - qp_samples)
    
        num_samples = 6
        replace=False
        if len(non_qp_samples)<num_samples:
            replace=True    
        target_samples = choice(non_qp_samples, num_samples,replace)
        
    
        sample_sfs_grid = gridspec.GridSpecFromSubplotSpec(1, num_samples, width_ratios=[1]*6, subplot_spec=species_sfs_grid[species_idx], wspace=0.1)
    
        sample_axes = []
        for sample_idx in xrange(0,num_samples):    
                
            sfs_axis = plt.Subplot(fig, sample_sfs_grid[sample_idx])
            fig.add_subplot(sfs_axis)
            sfs_axis.set_xticks([0.10*i for i in xrange(6,10)])
            sfs_axis.set_xlim([0.50,1.00])
            sfs_axis.set_yticks([])
            sfs_axis.xaxis.tick_bottom()

            if species_idx==11:
                sfs_axis.set_xlabel('Major allele freq')
            else:
                sfs_axis.set_xticklabels([])
            
            sample_axes.append(sfs_axis)

        sample_axes[0].set_ylabel(species_name_label_multiline, fontsize=5)
    
        for sample_idx in xrange(0,num_samples):
            sample = target_samples[sample_idx]
            sfs_axis = sample_axes[sample_idx]
            
            sfs_axis.set_title('$\\overline{D}=%d$' % (sample_coverage_map[sample]),fontsize=5,y=0.9)

            # Sample 1
            fs,pfs = sfs_utils.calculate_binned_sfs_from_sfs_map(sfs_map[sample],folding='major')
            df = fs[1]-fs[0]
    
            within_sites, between_sites, total_sites = sfs_utils.calculate_polymorphism_rates_from_sfs_map(sfs_map[sample])

            between_line = between_sites*1.0/total_sites/((fs>0.2)*(fs<0.5)).sum()
            #pmax = numpy.max([pfs[(fs>0.1)*(fs<0.95)].max(), between_line])
            pmax = between_line
            pmax = numpy.max([pfs[(fs>0.1)*(fs<0.9)].max(), between_line])
            
            within_rate = within_sites*1.0/total_sites
            #print "Sample 1: within =", within_rate, "avg-distance =", between_sites*1.0/total_sites

            sfs_axis.fill_between([0.80,1.00],[0,0],[1,1],color='0.8')

            sfs_axis.bar((fs-df/2),pfs,width=df, edgecolor=light_haploid_color, color=light_haploid_color)
            sfs_axis.set_ylim([0,pmax*3])

        sample_names = []

        between_rates = []

        within_rates = []
        within_rate_lowers = []
        within_rate_uppers = []

        median_depths = []
        depth_lowers = []
        depth_uppers = []
        # Now do polymorphism axis
        for sample in highcoverage_samples:
            within_sites, between_sites, total_sites = sfs_utils.calculate_polymorphism_rates_from_sfs_map(sfs_map[sample])

            within_rate = within_sites*1.0/total_sites
            between_rate = between_sites*1.0/total_sites    
            between_rates.append(between_rate)
            within_rates.append(within_rate)
    
            # Calculate 95% confidence intervals
            within_rate_lower, within_rate_upper = stats_utils.calculate_poisson_rate_interval(within_sites, total_sites,alpha=0.05)
            within_rate_lowers.append(within_rate_lower)
            within_rate_uppers.append(within_rate_upper)
   
            depths, counts = sfs_utils.calculate_depth_distribution_from_sfs_map(sfs_map[sample])
            dlower, dupper = stats_utils.calculate_IQR_from_distribution(depths, counts)
            dmedian = stats_utils.calculate_median_from_distribution(depths,counts)
    
            depth_lowers.append(dlower)
            depth_uppers.append(dupper)
            median_depths.append(dmedian)
            sample_names.append(sample)
    
        # Sort them all in descending order of within-host diversity    
        within_rates, within_rate_lowers, within_rate_uppers, between_rates, median_depths, depth_lowers, depth_uppers, sample_names = (numpy.array(x) for x in zip(*sorted(zip(within_rates, within_rate_lowers, within_rate_uppers, between_rates, median_depths,depth_lowers, depth_uppers, sample_names),reverse=True)))

        within_rate_lowers = numpy.clip(within_rate_lowers, 1e-09,1)
    
        row_idx = ((list_idx*12)+species_idx)/4
        col_idx = species_idx % 4
        polymorphism_axis = plt.Subplot(polymorphism_fig, polymorphism_grid[row_idx, col_idx])
        polymorphism_fig.add_subplot(polymorphism_axis)
        
        print species_name, row_idx, col_idx    
        
        rank_idxs = numpy.arange(0,len(within_rates))
        
        #polymorphism_axis.fill_between(rank_idxs, within_rate_lowers, within_rate_uppers,color='0.7')
        #polymorphism_axis.semilogy(rank_idxs, within_rate_uppers,'-',color='k',linewidth=0.25) 
        #polymorphism_axis.semilogy(rank_idxs, within_rate_lowers,'-',color='k',linewidth=0.25) 
           
        for rank_idx in xrange(0,len(within_rates)):
    
            polymorphism_axis.semilogy([rank_idx,rank_idx], [within_rate_lowers[rank_idx],within_rate_uppers[rank_idx]],'-',color='k',linewidth=0.25)
            
        polymorphism_axis.set_xticks([])
        polymorphism_axis.set_ylim([1e-06,1e-01])
        
        if col_idx>0.5:
            polymorphism_axis.set_yticklabels([])
        else:
            polymorphism_axis.set_ylabel('Polymorphism')
        
        polymorphism_axis.set_title(species_name_label_singleline, fontsize=5,y=1)

        if row_idx==5:
            polymorphism_axis.set_xlabel('Ranked samples')

            
    sys.stderr.write("Saving figure...\t")
    fig.savefig('%s/supplemental_sfs_%d.pdf' % (parse_midas_data.analysis_directory, figure_idx), bbox_inches='tight')
    sys.stderr.write("Done!\n")
    polymorphism_fig.savefig('%s/supplemental_polymorphism.pdf' % (parse_midas_data.analysis_directory), bbox_inches='tight')
    sys.stderr.write("Done!\n")



