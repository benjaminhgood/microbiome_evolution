###############################################################################
#
# Set up default source and output directories
#
###############################################################################
import os.path 

data_directory = os.path.expanduser("~/ben_nandita_hmp_data/")
analysis_directory = os.path.expanduser("~/ben_nandita_hmp_analysis/")
scripts_directory = os.path.expanduser("~/ben_nandita_hmp_scripts/")
patric_directory = os.path.expanduser("~/patric_db/")
midas_directory = os.path.expanduser("~/midas_db/")

# We use this one to debug because it was the first one we looked at
debug_species_name = 'Bacteroides_uniformis_57318'

good_species_min_coverage = 10
good_species_min_prevalence = 10

min_median_coverage = 20

consensus_lower_threshold = 0.2
consensus_upper_threshold = 0.8
fixation_min_change = 0.8

threshold_within_between_fraction = 0.1
threshold_pi = 1e-03

modification_difference_threshold = 20