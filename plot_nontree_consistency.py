import config
import parse_midas_data
import parse_HMP_data
import os.path
import pylab
import sys
import numpy

import diversity_utils
import gene_diversity_utils
import calculate_snv_distances
import stats_utils
from math import log10,ceil,exp
from numpy.random import randint, normal
    
#allowed_variant_types = set(['1D','2D','3D','4D'])
allowed_variant_types = set(['4D'])


if __name__=='__main__':


    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--debug", help="Loads only a subset of SNPs for speed", action="store_true")
    parser.add_argument("--chunk-size", type=int, help="max number of records to load", default=1000000000)
    parser.add_argument("--species", help="Name of specific species to run code on", default="Bacteroides_vulgatus_57955")
    args = parser.parse_args()

    debug = args.debug
    chunk_size = args.chunk_size
    species_name=args.species


    # Load SNP information for species_name
    sys.stderr.write("Loading SNPs for %s...\n" % species_name)
    sys.stderr.write("(core genes only...)\n")
    snv_distance_map = calculate_snv_distances.load_snv_distance_map(species_name)
    
    ds = numpy.logspace(-4,-1.5,15) # 15 points are plotted
    total_snps = numpy.zeros_like(ds)
    inconsistent_snps = numpy.zeros_like(ds)
      
    for location_tuple in snv_distance_map:
        var_type, derived_allele_counts, ancestral_allele_counts, between_d, within_d1, within_d2 = snv_distance_map[location_tuple]
        
        if var_type in allowed_variant_types:
            
            within_d = min([within_d1, within_d2])
            good_idxs = (ds>=between_d)
            inconsistent_idxs = good_idxs*(within_d>=2*ds)
            
            total_snps[good_idxs] += 1
            inconsistent_snps[inconsistent_idxs] += 1
    
    fraction_inconsistent = inconsistent_snps*1.0/(total_snps+(total_snps==0))
    
    print ds
    print total_snps
    print inconsistent_snps
    
    pylab.figure()
    pylab.xlabel('Divergence, $d$')
    pylab.ylabel('Fraction inconsistent')
    pylab.xlim([1e-04,1e-01])
    pylab.ylim([0,1.05])
    
    pylab.semilogx(ds[total_snps>0], fraction_inconsistent[total_snps>0],'b.-',markersize=3)
    pylab.savefig('non_tree_consistency.pdf',bbox_inches='tight')