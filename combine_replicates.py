import sys
import parse_midas_data
species=sys.argv[1]  

# load HMP sample data
subject_sample_map = parse_midas_data.parse_subject_sample_map()

combination_type = "sample" # at the moment, can either be sample or subject
if combination_type == "sample":
    grouping_replicate_map = parse_midas_data.flatten_samples(subject_sample_map)
elif combination_type == "subject":
    grouping_replicate_map = parse_midas_data.flatten_subjects(subject_sample_map)
else:
    sys.stderr.write("Improper combination type!\n")
    sys.exit(1)
    
sys.stderr.write("Combining marker gene coverages in a global file...\n")
parse_midas_data.combine_marker_gene_replicates(grouping_replicate_map, combination_type)
# (problem is, this includes samples that are later excluded from
#  the per-species MIDAS output because their coverage was too low)
  
# Then process things on a species-by-species basis...      

sys.stderr.write("Processing %s...\n" % species)
    
sys.stderr.write("Processsing marker gene coverages...")
parse_midas_data.combine_species_marker_gene_replicates(species, grouping_replicate_map, combination_type)
sys.stderr.write("\tDone!\n")
    
sys.stderr.write("Processing SNPs...\n")
parse_midas_data.combine_replicates(species, grouping_replicate_map, combination_type)
sys.stderr.write("Done!\n")
    
sys.stderr.write("Finished!\n")
