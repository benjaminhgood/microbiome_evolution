import parse_midas_data
import pylab
import sys
import numpy
import bz2
import calculate_snp_prevalences



################################################################################
#
# Standard header to read in argument information
#
################################################################################
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("species_name", help="name of species to process")
parser.add_argument("--debug", help="Loads only a subset of SNPs for speed", action="store_true")
parser.add_argument("--chunk-size", type=int, help="max number of records to load", default=1000000000)
args = parser.parse_args()

species_name = args.species_name
debug = args.debug
chunk_size = args.chunk_size
################################################################################


# Should we do this? 
sys.stderr.write("Loading core genes...\n")
core_genes = parse_midas_data.load_core_genes(species_name)
sys.stderr.write("Done! %d core genes\n" % len(core_genes))
allowed_genes = core_genes

sys.stderr.write("Loading population freqs...\n")
population_freqs = calculate_snp_prevalences.parse_population_freqs(species_name)
sys.stderr.write("Done! %d SNVs\n" % len(population_freqs))

allowed_variant_type_list = ['1D','2D','3D','4D']
allowed_variant_types = set(allowed_variant_type_list)  

# Open post-processed MIDAS output
snp_file =  bz2.BZ2File("%ssnps/%s/annotated_snps.txt.bz2" % (parse_midas_data.data_directory, species_name),"r")
    
line = snp_file.readline() # header
items = line.split()[1:]
samples = numpy.array([item.strip() for item in items])

# We shouldn't be doing this for raw data 
#samples = parse_midas_data.parse_merged_sample_names(items)
    
site_map = [{} for sample in samples]
for sample_idx in xrange(0,len(samples)):
    site_map[sample_idx] = {variant_type:{} for variant_type in allowed_variant_types}

sys.stderr.write("Calculating within-person SFSs...\n")        
num_sites_processed = 0
for line in snp_file:
    #
    items = line.split()
    # Load information about site
    info_items = items[0].split("|")
    chromosome = info_items[0]
    location = long(info_items[1])
    gene_name = info_items[2]
    variant_type = info_items[3]
    
    if len(info_items) > 5: # for backwards compatability
            polarization = info_items[4]
            pvalue = float(info_items[5])
    else: 
        polarization="?"
        pvalue = float(info_items[4])
        
    #    
    if variant_type not in allowed_variant_types:
        continue
    #    
    if len(allowed_genes)>0 and (gene_name not in allowed_genes):
        continue
    #    
    # Load alt and depth counts
    alts = []
    depths = []
    for item in items[1:]:
        subitems = item.split(",")
        alts.append(long(subitems[0]))
        depths.append(long(subitems[1]))
    alts = numpy.array(alts)
    depths = numpy.array(depths)
    refs = depths-alts
    #print alts
    #print depths
    #
    # population_freq returns the fraction of people for which the alt is the major allele.
    # This is a very important quantity being computed! It is later used for identifying CPS samples.
    if (chromosome, location) in population_freqs:
        population_freq = population_freqs[(chromosome, location)]
    else:
        population_freq = 0
    
    # polarize SFS according to population freq
    if population_freq>0.5:
        alts,refs = refs,alts
        population_freq = 1-population_freq
        
    #    
    for i in xrange(0,len(alts)):
        site = (depths[i],alts[i])
        #
        if site not in site_map[i][variant_type]:
            site_map[i][variant_type][site] = [0,0.0]
        #        
        site_map[i][variant_type][site][0] += 1
        site_map[i][variant_type][site][1] += population_freq # weight of polarization reversals
        #
        #
    num_sites_processed+=1
    #print num_sites_processed
    if num_sites_processed%50000==0:
        sys.stderr.write("%dk sites processed...\n" % (num_sites_processed/1000))   
        if debug:
            break
    
snp_file.close()
sys.stderr.write("Done!\n")
# Write to disk!
sys.stderr.write("Writing output...\n")
# First write (filtered) genome-wide coverage distribution
output_file = bz2.BZ2File("%ssnps/%s/within_sample_sfs.txt.bz2" % (parse_midas_data.data_directory, species_name),"w")
output_file.write("\t".join(["SampleID", "variant_type", "D,A,count,reverse_count", "..."]))
for sample_idx in xrange(0,len(samples)):
    sample = samples[sample_idx]
    for variant_type in allowed_variant_type_list:
        output_file.write("\n")
        output_file.write("\t".join([sample, variant_type]+["%d,%d,%d,%g" % (site[0],site[1],site_map[sample_idx][variant_type][site][0],site_map[sample_idx][variant_type][site][1]) for site in sorted(site_map[sample_idx][variant_type].keys())]))
output_file.close()
sys.stderr.write("Done!\n")