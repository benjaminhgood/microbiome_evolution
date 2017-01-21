# This module contains functions and helper methods to filter out sequencing errors from our metagenomic datasets

import numpy
from scipy.stats import binom

def calculate_alt_thresholds(Ds, min_alt=2, max_error_probability=1e-02, pvalue=0.05):

    # returns num alts such that p-value < 0.05
    # in simple binomial model with max error probability
    
    return numpy.fmax(binom.isf(pvalue,Ds, max_error_probability), min_alt)
    
if __name__=='__main__':

    Dmin = 10
    Dmax = 200

    Ds = numpy.arange(Dmin,Dmax+1)
    min_alts = calculate_alt_thresholds(Ds)
    
    import pylab
    pylab.plot(Ds, min_alts, 'k.-')
    pylab.ylim([0,5])
    pylab.show()


    
    
    
    