import sys
species=sys.argv[1]

#####
#
# Creates a coverage_distribution.txt.bz2 in the species snp directory
#
# In this file rows are samples, columns are D,count pairs
# samples are guaranteed to be same order as snps_depth.txt.bz2 file
#
#####

from parse_midas_data import *

#default_directory_prefix =  os.path.expanduser("~/Documents/files_too_big_for_Dropbox/midas_output_121816/")
combination_type = "sample"
snp_prefix = get_snp_prefix_from_combination_type(combination_type)

allowed_variant_types = set(["1D","2D","3D","4D"]) # use all types of sites to include most information
    


sys.stderr.write("Processing %s...\n" % species)
 
depth_file = bz2.BZ2File("%ssnps/%s/%ssnps_depth.txt.bz2" % (default_directory_prefix, species, snp_prefix))
info_file = bz2.BZ2File("%ssnps/%s/snps_info.txt.bz2" % (default_directory_prefix, species),"r")
    
depth_line = depth_file.readline() # header
info_line = info_file.readline()
    
samples = depth_line.split()[1:]
    
sample_depth_histograms = {sample: {} for sample in samples}
gene_total_depths = {}
gene_total_sites = {}

num_sites_processed = 0
    
while True:
            
    # load next lines
    depth_line = depth_file.readline()
    info_line = info_file.readline()
        
    # quit if file has ended
    if depth_line=="":
        break
        
    # parse site info
    info_items = info_line.split()
    variant_type = info_items[5]
        
    # make sure it is either a 1D or 4D site
    if not variant_type in allowed_variant_types:
        continue
    
    gene_name = info_items[6]
        
    items = depth_line.split()   
    depths = numpy.array([long(item) for item in items[1:]])
        
    # Add to genome-wide depth distribution
    for sample,D in zip(samples,depths):
        if D not in sample_depth_histograms[sample]:
            sample_depth_histograms[sample][D]=0
        sample_depth_histograms[sample][D] += 1
    
    # Add to gene-specific avg depth
    if gene_name not in gene_total_depths:
        gene_total_depths[gene_name] = numpy.zeros(len(samples))*1.0
        gene_total_sites[gene_name] = 0.0
        
    gene_total_depths[gene_name]+=depths
    gene_total_sites[gene_name]+=1
    
    num_sites_processed+=1
        
    if num_sites_processed%100000==0:
        sys.stderr.write("Processed %dk sites!\n" % (num_sites_processed/1000))
            
depth_file.close()
    
# Now write output!
    
# First write genome-wide coverage distribution
output_file = bz2.BZ2File("%ssnps/%s/%scoverage_distribution.txt.bz2" % (default_directory_prefix, species, snp_prefix),"w")
output_file.write("SampleID\tD,n(D) ...")
for sample in samples:
    output_file.write("\n")
    output_file.write("\t".join([sample]+["%d,%d" % (D,sample_depth_histograms[sample][D]) for D in sorted(sample_depth_histograms[sample].keys())]))
output_file.close()

# Then write gene-specific coverages
output_file = bz2.BZ2File("%ssnps/%s/%sgene_coverage.txt.bz2" % (default_directory_prefix, species, snp_prefix),"w")
output_file.write("\t".join(["Gene"]+samples)) # Header line
for gene_name in sorted(gene_total_depths.keys()):
    avg_depths = gene_total_depths[gene_name]/(gene_total_sites[gene_name]+(gene_total_sites[gene_name]==0))
    output_file.write("\n")
    output_file.write("\t".join([gene_name]+["%0.1f" % D for D in avg_depths]))
output_file.close()
    
    # Done!
    
