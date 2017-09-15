import matplotlib  
matplotlib.use('Agg') 
import config
import parse_midas_data
import parse_HMP_data

import pylab
import sys
import numpy
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


# Load marker gene coverages
species_coverage_matrix, samples, species = parse_midas_data.parse_global_marker_gene_coverages()
total_coverages = species_coverage_matrix.sum(axis=0)

relative_frequencies = species_coverage_matrix/(total_coverages[None,:]+(total_coverages[None,:]==0))

pylab.figure(1,figsize=(3.43,2.5))
pylab.xlabel('Top $n$ species')
pylab.ylabel('1 - cumulative abundance')

pylab.gca().spines['top'].set_visible(False)
pylab.gca().spines['right'].set_visible(False)
pylab.gca().get_xaxis().tick_bottom()
pylab.gca().get_yaxis().tick_left()


total_relative_frequencies = relative_frequencies.sum(axis=1)/(total_coverages>0).sum()
total_relative_frequencies = total_relative_frequencies[total_relative_frequencies>=1e-04]
total_relative_frequencies.sort()
total_relative_frequencies = numpy.flipud(total_relative_frequencies)

cumulative_total_frequencies = total_relative_frequencies.cumsum()

# construct rank abundance curves
for j in xrange(0,len(samples)):

    sample_coverages = species_coverage_matrix[:,j]
    sample_freqs = relative_frequencies[:,j]
    
    sample_freqs = sample_freqs[sample_coverages>=1]
    
    # sort in ascending order
    sample_freqs.sort()
    
    # sort in descending order
    sample_freqs = numpy.flipud(sample_freqs)
    cumulative_sample_freqs = numpy.cumsum(sample_freqs)
    
    pylab.loglog(range(1,len(sample_freqs)+1), 1-cumulative_sample_freqs,'-',linewidth=0.5,alpha=0.5)
    
pylab.loglog([1e03,1e03],[1,1],'b-',linewidth=0.5,alpha=0.5,label='Individual hosts')
    
pylab.loglog(range(1,len(total_relative_frequencies)+1), 1-cumulative_total_frequencies,'k-',linewidth=3,label='Pooled')

pylab.xlim([5e-01,2e02])
pylab.ylim([1e-03,1])
pylab.legend(loc='upper right',frameon=False)
pylab.savefig('%scumulative_abundances.pdf' % (parse_midas_data.analysis_directory),bbox_inches='tight')

