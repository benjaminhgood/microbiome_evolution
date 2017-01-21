#####
#
# Creates a marker_coverage.txt.bz2 in the species snp directory
#
#####

from parse_midas_data import *

#default_directory_prefix =  os.path.expanduser("~/Documents/files_too_big_for_Dropbox/midas_output_121816/")

#species_names = ["Bacteroides_uniformis_57318"]
#species_names = ["Roseburia_intestinalis_56239", "Eubacterium_eligens_61678"]
species_names = ["Prevotella_copri_61740"]
#species_names=["Butyrivibrio_crossotus_61674"]

for species in species_names:

    sys.stderr.write("Processing %s...\n" % species)
 
    depth_file = bz2.BZ2File("%ssnps/%s/snps_depth.txt.bz2" % (default_directory_prefix, species),"r")
    
    # get list of samples to use
    depth_line = depth_file.readline()
    depth_items = depth_line.split()
    output_samples = depth_items[1:]
    depth_file.close()
    
    coverage_file = bz2.BZ2File("%sspecies/coverage.txt.bz2" % (default_directory_prefix),"r")

    output_coverage_file = bz2.BZ2File("%ssnps/%s/marker_coverage.txt.bz2" % (default_directory_prefix, species),"w")
    
    # get header line
    line = coverage_file.readline()
    
    # get list of samples
    items = line.split()
    samples = items[1:]

    output_idxs = numpy.array([samples.index(sample) for sample in output_samples])
    
    
    # write header lines for output file
    output_coverage_file.write("\t".join([items[0]]+output_samples))
    
    total_depths_defined = False
    while True:
            
        # load next lines
        line = coverage_file.readline()
        
        # quit if file has ended
        if line=="":
            break
        
        items = line.split()
        
        current_species = items[0]
        depths = numpy.array([float(item) for item in items[1:]])[output_idxs]
        
        if current_species==species:
            # write output
            output_coverage_file.write("\n")
            output_coverage_file.write("\t".join([current_species]+["%g" % d for d in depths]))
        
        if not total_depths_defined:
            total_depths = numpy.zeros_like(depths)
            total_depths_defined=True
                
        total_depths += depths
        
    
    output_coverage_file.write("\n")
    output_coverage_file.write("\t".join(["Total"] + ["%g" % d for d in total_depths]))
    
    coverage_file.close()
    output_coverage_file.close()
    
    # Done!
    # no return value
