import sys
import bz2
import numpy
import stats_utils

import parse_midas_data
if len(sys.argv) > 1:
    species_name=sys.argv[1]
else:
    species_name=parse_midas_data.debug_species_name

sys.stderr.write("Calculating marker gene coverage distribution for %s...\n" % species_name)

#####
#
# Creates a coverage_distribution.txt.bz2 file in the species snp directory
# 
# In this file, rows are samples, columns are D,count pairs
# samples are guaranteed to be same order as snps_depth.txt.bz2 file
#
# Also creates a gene_coverage.txt.bz2 file in the species snp directory
# In this file, rows are genes, columns are samples, 
# entries are avg coverage of that gene for that sample 
#####

allowed_variant_types = set(["1D","2D","3D","4D"]) # use all types of sites to include most information

# Load list of marker genes
marker_genes = parse_midas_data.load_marker_genes(species_name)

# read depth information straight from snps_depth.txt. No preprocessing
depth_file = bz2.BZ2File("%ssnps/%s/snps_depth.txt.bz2" % (parse_midas_data.data_directory, species_name),"r")

# To get information about which gene it is in
info_file = bz2.BZ2File("%ssnps/%s/snps_info.txt.bz2" % (parse_midas_data.data_directory, species_name),"r")
    
depth_line = depth_file.readline() # header
info_line = info_file.readline()
    
samples = depth_line.split()[1:]

marker_coverages = {gene_name: [] for gene_name in marker_genes}
marker_gene_starts = {gene_name: -1 for gene_name in marker_genes}
marker_locations = {gene_name: [] for gene_name in marker_genes}
        
num_sites_processed = 0
total_marker_sites = 0
    
while True:
    
    num_sites_processed+=1
    
    if num_sites_processed%100000==0:
        sys.stderr.write("Processed %dk sites!\n" % (num_sites_processed/1000))
    
            
    # load next lines
    depth_line = depth_file.readline()
    info_line = info_file.readline()
        
    # quit if file has ended
    if depth_line=="":
        break
        
    # parse site info
    info_items = info_line.split('\t')
    variant_type = info_items[5]
    
    if not variant_type in allowed_variant_types:
        continue
    
    gene_name = info_items[6].strip()
    
    if gene_name not in marker_genes:
        continue
        
    items = depth_line.split('\t')   
    
    subitems = items[0].split("|")
    location = long(subitems[1])
    
    depths = numpy.array([long(item) for item in items[1:]])
    
    if marker_gene_starts[gene_name]<0:
        marker_gene_starts[gene_name]=location
        
    marker_locations[gene_name].append(location-marker_gene_starts[gene_name])
    marker_coverages[gene_name].append(depths)
            
depth_file.close()
    
for gene_name in marker_coverages.keys():
    marker_coverages[gene_name] = numpy.array(marker_coverages[gene_name]).T
    
# Now write output!
        
# First write genome-wide coverage distribution
output_file = bz2.BZ2File("%ssnps/%s/marker_coverage_distribution.txt.bz2" % (parse_midas_data.data_directory, species_name),"w")
output_file.write("SampleID,GeneID\tD,n(D) ...")
for sample_idx in xrange(0,len(samples)):
    sample=samples[sample_idx]
    for gene_name in marker_genes:
        output_file.write("\n")
        output_file.write("\t".join(["%s,%s" % (sample,gene_name)]+["%d,%d" % (marker_locations[gene_name][i], marker_coverages[gene_name][sample_idx,i]) for i in xrange(0,len(marker_locations[gene_name]))]))
    
output_file.close()
 
# Done!
    
