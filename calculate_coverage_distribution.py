import sys
import bz2
import numpy

import parse_midas_data
if len(sys.argv) > 1:
    species_name=sys.argv[1]
else:
    species_name=parse_midas_data.debug_species_name

sys.stderr.write("Calculating coverage distribution for %s...\n" % species_name)

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

# These are the default MIDAS parameters. Used to ensure consistency
prevalence_threshold = 0.95 #only sites that pass this prevalence threshold in terms of having at least the prevalence_min_coverage can pass. These sites are stored in the sample_depth_histograms
prevalence_min_coverage = 3

allowed_variant_types = set(["1D","2D","3D","4D"]) # use all types of sites to include most information

depth_file = bz2.BZ2File("%ssnps/%s/snps_depth.txt.bz2" % (parse_midas_data.data_directory, species_name),"r")
info_file = bz2.BZ2File("%ssnps/%s/snps_info.txt.bz2" % (parse_midas_data.data_directory, species_name),"r")
    
depth_line = depth_file.readline() # header
info_line = info_file.readline()
    
samples = depth_line.split()[1:]
    
sample_depth_histograms = {sample: {} for sample in samples} # stores only those sites that pass the prevalence threshold above.
full_sample_depth_histograms = {sample: {} for sample in samples} # stores all sites. 


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
    info_items = info_line.split('\t')
    variant_type = info_items[5]
        
    # make sure it is either a 1D or 4D site
    if not variant_type in allowed_variant_types:
        continue
    
    gene_name = info_items[6]
        
    items = depth_line.split()   
    depths = numpy.array([long(item) for item in items[1:]])
    
    # Manual prevalence filter
    if (depths>=prevalence_min_coverage).sum()*1.0/len(depths)>=prevalence_threshold: 
        # Add to genome-wide depth distribution
        for sample,D in zip(samples,depths):
            if D not in sample_depth_histograms[sample]:
                sample_depth_histograms[sample][D]=0
            sample_depth_histograms[sample][D] += 1
    
    for sample,D in zip(samples,depths):
        if D not in full_sample_depth_histograms[sample]:
            full_sample_depth_histograms[sample][D]=0
        full_sample_depth_histograms[sample][D] += 1
        
    
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
        
# First write (filtered) genome-wide coverage distribution. This is filtered by prevalence
output_file = bz2.BZ2File("%ssnps/%s/coverage_distribution.txt.bz2" % (parse_midas_data.data_directory, species_name),"w")
output_file.write("SampleID\tD,n(D) ...")
for sample in samples:
    output_file.write("\n")
    output_file.write("\t".join([sample]+["%d,%d" % (D,sample_depth_histograms[sample][D]) for D in sorted(sample_depth_histograms[sample].keys())]))
output_file.close()

# Write unfiltered genome-wide coverage distribution
output_file = bz2.BZ2File("%ssnps/%s/full_coverage_distribution.txt.bz2" % (parse_midas_data.data_directory, species_name),"w")
output_file.write("SampleID\tD,n(D) ...")
for sample in samples:
    output_file.write("\n")
    output_file.write("\t".join([sample]+["%d,%d" % (D, full_sample_depth_histograms[sample][D]) for D in sorted(full_sample_depth_histograms[sample].keys())]))
output_file.close()

# Then write gene-specific coverages
output_file = bz2.BZ2File("%ssnps/%s/gene_coverage.txt.bz2" % (parse_midas_data.data_directory, species_name),"w")
output_file.write("\t".join(["Gene"]+samples)) # Header line
for gene_name in sorted(gene_total_depths.keys()):
    avg_depths = gene_total_depths[gene_name]/(gene_total_sites[gene_name]+(gene_total_sites[gene_name]==0))
    output_file.write("\n")
    output_file.write("\t".join([gene_name]+["%0.1f" % D for D in avg_depths]))
output_file.close()
    
# Done!
    
