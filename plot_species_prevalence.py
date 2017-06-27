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

mpl.rcParams['font.size'] = 7
mpl.rcParams['lines.linewidth'] = 0.5
mpl.rcParams['legend.frameon']  = False
mpl.rcParams['legend.fontsize']  = 'small'



min_marker_coverage = 20
min_prevalence=20

bacteroides_color = '#084594'
alistipes_color = '#B10026'
rest_color = '0.7'

good_species_list = []
prevalences = []
permissive_prevalences = []
    
species_coverage_matrix, samples, species = parse_midas_data.parse_global_marker_gene_coverages()
for i in xrange(0,len(species)):
        
    species_coverages = species_coverage_matrix[i,:]
    prevalence = (species_coverages>=min_marker_coverage).sum()
    if prevalence >= min_prevalence:
        good_species_list.append(species[i])
        prevalences.append(prevalence)
        permissive_prevalences.append( (species_coverages>=10).sum() )
    
pylab.figure(figsize=(5,2))
axis = pylab.gca()

axis.spines['top'].set_visible(False)
axis.spines['right'].set_visible(False)
axis.get_xaxis().tick_bottom()
axis.get_yaxis().tick_left()



prevalences, permissive_prevalences, good_species_list = zip(*sorted(zip(prevalences, permissive_prevalences, good_species_list),reverse=True))
locs = numpy.arange(0,len(prevalences))

pretty_species_names = []
for i in xrange(0,len(prevalences)):

    loc = locs[i]
    prevalence = prevalences[i]
    permissive_prevalence = permissive_prevalences[i]
    species_name = good_species_list[i]
    if species_name.startswith('Bacteroides'):
        color = bacteroides_color
    elif species_name.startswith('Alistipes'):
        color = alistipes_color
    else:
        color = rest_color

    pylab.bar([loc], [permissive_prevalence], width=0.8, linewidth=0,facecolor=color,alpha=0.5)
    pylab.bar([loc], [prevalence], width=0.8, linewidth=0,facecolor=color)
    

    pretty_species_names.append( "%s %s (%s)" % tuple(species_name.split("_")))

axis.set_xticks(locs+0.4)
axis.set_xticklabels(pretty_species_names, rotation='vertical',fontsize=3)

axis.set_xlim([-0.4,len(good_species_list)+0.4])
axis.set_ylabel('Number of samples >20x')

axis.tick_params(axis='y', direction='out',length=3,pad=1)
axis.tick_params(axis='x', direction='out',length=3,pad=1)


pylab.savefig('%s/species_prevalences.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')

    