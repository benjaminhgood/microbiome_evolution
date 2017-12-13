##############################################################################
#
# Utility methods for working with (primarily within-sample) SFSs
#
##############################################################################
import numpy
import config
from math import log10

def calculate_binned_sfs_from_sfs_map(sfs_map, bins=[], folding='minor'):

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
    
    if len(bins)==0:
        # use this to set up bins
        bins = (numpy.arange(0,Dbar+2)-0.5)/Dbar
        fs = numpy.arange(0,Dbar+1)*1.0/Dbar
        
    else:
        bins = numpy.array(bins)
        fs = bins[1:]
    
    pfs = numpy.zeros_like(fs)
    bin_idxs = numpy.digitize(minor_freqs, bins=bins)
    
    for bin_idx, weight in zip(bin_idxs, weights):
        pfs[bin_idx-1] += weight
    
    # should already be normalized, but just to make sure...
    pfs /= pfs.sum()
    
    if folding=='major':
        pfs = pfs[::-1]
        fs = (1.0-fs)[::-1]
    
    return fs,pfs
    
def calculate_binned_depth_distribution_from_sfs_map(sfs_map,bins=[],num_bins=30):

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
    
    # calculate median depth (or rough approximation)
    sorted_depths, sorted_weights = (numpy.array(x) for x in zip(*sorted(zip(depths, weights))))
    CDF = numpy.cumsum(sorted_weights)
    Dbar = sorted_depths[CDF>0.5][0]
    
    if len(bins)==0:
        # use this to set up bins
        bins = numpy.logspace(log10(Dbar/8),log10(Dbar*8),num_bins)
        Ds = bins[0:-1]
           
    else:
        bins = numpy.array(bins, copy=True)
        Ds = bins[0:-1]
    
    bins[0] = 0
    bins[-1] = 1e09
    
    pDs = numpy.zeros_like(Ds)
    bin_idxs = numpy.digitize(depths, bins=bins)
    
    for bin_idx, weight in zip(bin_idxs, weights):
        pDs[bin_idx-1] += weight
    
    # should already be normalized, but just to make sure...
    pDs /= pDs.sum()
    
    return bins, Ds, pDs

def calculate_depth_distribution_from_sfs_map(sfs_map):
    
    depth_map = {}
    for key in sfs_map.keys():
        D,A = key
        count = sfs_map[key][0]
        
        if D not in depth_map:
            depth_map[D] = 0
        
        depth_map[D] += count
        
    depths = numpy.array(sorted(depth_map.keys()))
    counts = numpy.array([depth_map[d] for d in depths])
    
    return depths, counts

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
                between_sites += (n-reverse_n) #NRG what does this do?
            else:
                between_sites += reverse_n
        
        
    return within_sites, between_sites, total_sites
    
    
