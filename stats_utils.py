import numpy
from math import log
from scipy.stats import gamma

####
#
# Calculates median from histogram
# histogram is map of value: counts
#
####
def calculate_median_from_histogram(histogram):

    xs, CDF = calculate_CDF_from_histogram(histogram)
    median_idx = numpy.nonzero(CDF>=0.5)[0][0]
    return xs[median_idx]

####
#
# Calculates median from histogram
# histogram is map of value: counts
#
####
def calculate_nonzero_median_from_histogram(histogram):

    

    xs, CDF = calculate_CDF_from_histogram(histogram)
    
    if len(xs)<2:
        return xs[0]

    #if CDF[-1]==0:
    #    return xs[0]
    
    
    if xs[0]<0.5:
    
        if CDF[0]>0.8:
            return xs[0]
    
        CDF -= CDF[0]
        CDF /= CDF[-1]
        
    median_idx = numpy.nonzero(CDF>=0.5)[0][0]
    return xs[median_idx]

####
#
# Calculates median from histogram
# histogram is map of value: counts
#
####
def calculate_thresholded_median_from_histogram(histogram,xmin=0):


    xs, CDF = calculate_CDF_from_histogram(histogram)

    
    
    # Get last index below xmin
    idx = numpy.nonzero(xs>xmin+0.5)[0][0]-1
    
    CDF -= CDF[idx]
    CDF /= CDF[-1]
    CDF = numpy.clip(CDF,0,1e09)
    
    median_idx = numpy.nonzero(CDF>=0.5)[0][0]
    return xs[median_idx]

 

####
#
# Calculates CDF from histogram
# histogram is map of value: counts
#
####
def calculate_unnormalized_CDF_from_histogram(histogram):

    xs = sorted(histogram.keys())
    ns = numpy.array([histogram[x] for x in xs])*1.0
    xs = numpy.array(xs)*1.0
    
    CDF = ns.cumsum()
    return xs, CDF
    
####
#
# Calculates CDF from histogram
# histogram is map of value: counts
#
####
def calculate_CDF_from_histogram(histogram):

    xs = sorted(histogram.keys())
    ns = numpy.array([histogram[x] for x in xs])*1.0
    xs = numpy.array(xs)*1.0
    
    CDF = ns.cumsum()/ns.sum()
    return xs, CDF


def calculate_total_from_histogram(histogram):
    xs = sorted(histogram.keys())
    ns = numpy.array([histogram[x] for x in xs])*1.0
    return ns.sum()
  
####
#
# Calculates unnormalized survival functions (i.e. # of observations >= x)
# from numpy vector of observations
#
####
def calculate_unnormalized_survival_from_vector(xs, min_x=None, max_x=None, min_p=1e-10, eps=1e-10):

    xs = numpy.array(xs)

    if min_x==None:
        min_x = xs.min()-1
    
    if max_x==None:
        max_x = xs.max()+1
        
    unique_xs = set(xs)
    unique_xs.add(min_x)
    unique_xs.add(max_x)
    
    xvalues = []
    num_observations = []
    
    for x in sorted(unique_xs):
        xvalues.append(x)
        num_observations.append( (xs>=x).sum() )
    
    # So that we can plot CDF, SF on log scale
    num_observations[0] -= min_p
    num_observations[1] -= min_p
    num_observations[-1] += min_p    
    
    return numpy.array(xvalues), numpy.array(num_observations)
 
def calculate_IQR_from_distribution(xs, ns):
    
    weights = ns*1.0/ns.sum()
    CDF = numpy.cumsum(weights)
    upper_idx = numpy.nonzero(CDF>=0.75)[0][0]
    lower_idx = numpy.nonzero(CDF>=0.25)[0][0]
    
    return xs[lower_idx], xs[upper_idx]

def calculate_median_from_distribution(xs, ns):
    weights = ns*1.0/ns.sum()
    CDF = numpy.cumsum(weights)
    mid_idx = numpy.nonzero(CDF>=0.5)[0][0]
    return xs[mid_idx]
    
####
#
# Calculates IQR (75-25 percentile) from histogram
# histogram is map of value: counts
#
####
def calculate_IQR_from_histogram(histogram):

    xs, CDF = calculate_CDF_from_histogram(histogram)
    
    upper_idx = numpy.nonzero(CDF>=0.75)[0][0]
    lower_idx = numpy.nonzero(CDF>=0.25)[0][0]
    
    return xs[upper_idx]-xs[lower_idx]
    

    
####
#
# Calculates "confidence intervals" on rate from Poisson distribution 
# based on n>=0 counts at L>=0 sites.
#
####
def calculate_poisson_rate_interval(n,L,alpha=0.5): # by default use a 50% confidence interval
    
    if n<0.5:
        # No counts. Have some info on upper bound, but none on lower bound.
        plower = 0
        pupper = log(2/alpha)/L
    
    else:
        # Posterior distribution is Gamma with shape n-1 and scale 1/L 
        # Get confidence intervals from tail probabilities
        plower = gamma.ppf(alpha/2, n)/L
        pupper = gamma.ppf(1-alpha/2,n)/L
        
    return plower,pupper
        