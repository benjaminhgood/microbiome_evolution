import numpy

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
    
    if xs[0]<0.5:
        CDF -= CDF[0]
        CDF /= CDF[-1]
        
    median_idx = numpy.nonzero(CDF>=0.5)[0][0]
    return xs[median_idx]
 
    
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

####
#
# Calculates unnormalized survival functions (i.e. # of observations >= x)
# from numpy vector of observations
#
####
def calculate_unnormalized_survival_from_vector(observations):

    unique_observations = set(observations)
    unique_observations.add(observations.min()-1)
    unique_observations.add(observations.max()+1)
    
    xvalues = []
    num_observations = []
    
    for x in sorted(unique_observations):
        xvalues.append(x)
        num_observations.append( (observations>=x).sum() )
    
    num_observations[-1] = 1e-12    
    return numpy.array(xvalues), numpy.array(num_observations)
    
    
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
 