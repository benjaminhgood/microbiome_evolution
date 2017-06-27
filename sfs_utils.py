##############################################################################
#
# Utility methods for working with (primarily within-sample) SFSs
#
##############################################################################
import numpy
import config

def calculate_binned_sfs_from_sfs_map(sfs_map):

    alts = []
    depths = []
    counts = []
    for key in sfs_map.keys():
        D,A = key
        n = sfs_map[key][0]
        
        alts.append(A)
        depths.append(D)
        counts.append(n)
    
    alts = numpy.array(alts)
    depths = numpy.array(depths)
    counts = numpy.array(counts)
    weights = counts*1.0/counts.sum()
    
    freqs = alts*1.0/depths
    minor_freqs = numpy.fmin(freqs,1-freqs)
    
    # calculate median depth (or rough approximation)
    sorted_depths, sorted_weights = (numpy.array(x) for x in zip(*sorted(zip(depths, weights))))
    CDF = numpy.cumsum(sorted_weights)
    Dbar = sorted_depths[CDF>0.5][0]
    
    # use this to set up bins
    bins = (numpy.arange(0,Dbar+2)-0.5)/Dbar
    fs = numpy.arange(0,Dbar+1)*1.0/Dbar
    pfs = numpy.zeros_like(fs)
    
    bin_idxs = numpy.digitize(minor_freqs, bins=bins)
    
    for bin_idx, weight in zip(bin_idxs, weights):
        pfs[bin_idx-1] += weight
    
    # should already be normalized, but just to make sure...
    pfs /= pfs.sum()
    
    return fs,pfs

def calculate_polymorphism_rates_from_sfs_map(sfs_map,lower_threshold=0.2,upper_threshold=0.8):
    
    total_sites = 0
    within_sites = 0
    between_sites = 0
    for key in sfs_map.keys():
        D,A = key
        n = sfs_map[key][0]
        reverse_n = sfs_map[key][1]
        
        f = A*1.0/D
        
        total_sites += n
        
        if ((f>lower_threshold) and (f<upper_threshold)):
            # an intermediate frequency site
            within_sites += n
        else:    
            if f>0.5:
                between_sites += (n-reverse_n)
            else:
                between_sites += reverse_n
        
        
    return within_sites, between_sites, total_sites
    