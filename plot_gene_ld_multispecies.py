import matplotlib  
matplotlib.use('Agg') 
import os
import parse_midas_data
import pylab
import sys
import numpy
import diversity_utils
from numpy.random import choice

#species_names=['Prevotella_copri_61740','Butyrivibrio_crossotus_61674','Bacteroides_uniformis_57318','Roseburia_intestinalis_56239','Eubacterium_eligens_61678']
species_names=['Alistipes_putredinis_61533','Bacteroides_uniformis_57318']

files={}

for species_name in species_names:
    files[species_name]=numpy.load(os.path.expanduser('~/tmp_intermediate_files/LD_%s.npz' % species_name))


colors=['#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e']
#colors=['r','b','k','g','t']


pylab.figure()
pylab.xlabel('Distance between SNPs')
pylab.ylabel("Ohta and Kimura's $\\sigma^2_d$")
pylab.xlim([1,1e04])
pylab.ylim([1e-02,1])

color_idx=0

handles=[] # this stores the plots so that I can make a legend below

for species_name in species_names:
    
    #plot the decay in LD in a gene
    handles.append(color_idx)
    handles[color_idx], = (pylab.loglog(files[species_name]['distance_bin_locations'][files[species_name]['binned_rsquared_denominators']>0], files[species_name]['binned_rsquareds'][files[species_name]['binned_rsquared_denominators']>0],'k.-', color=colors[color_idx], label=species_names))
    
    # plot LD across genes (control)
    pylab.loglog(files[species_name]['distance_bin_locations'],numpy.ones_like(files[species_name]['distance_bin_locations'])*files[species_name]['control_rsquareds'],'k:', color=colors[color_idx])
    color_idx+=1

# plot the legend
labels=species_names
pylab.legend(handles,labels,'lower left',prop={'size':6})

pylab.savefig('%s/intragene_ld_multispecies.pdf' % (parse_midas_data.analysis_directory), bbox_inches='tight')
pylab.savefig('%s/intragene_ld_multispecies.png' % (parse_midas_data.analysis_directory), bbox_inches='tight', dpi=300)
