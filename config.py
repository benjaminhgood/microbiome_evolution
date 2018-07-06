###############################################################################
#
# Set up default source and output directories
#
###############################################################################
import os.path 
from math import log10

data_directory = os.path.expanduser("~/ben_nandita_hmp_data_060518/")
#data_directory = os.path.expanduser("~/ben_nandita_hmp_data_051318/")
#data_directory = os.path.expanduser("~/ben_nandita_hmp_data/")
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
fixation_min_change = consensus_upper_threshold-consensus_lower_threshold
fixation_log10_depth_ratio_threshold = log10(3)

threshold_within_between_fraction = 0.1
threshold_pi = 1e-03

modification_difference_threshold = 20
#modification_difference_threshold = 100

gainloss_max_absent_copynum = 0.05
gainloss_min_normal_copynum = 0.5
gainloss_max_normal_copynum = 2

core_genome_min_copynum = 0.3
core_genome_max_copynum = 3 # BG: should we use a maximum for "core genome"? I'm going to go w/ yes for now
core_genome_min_prevalence = 0.9
shared_genome_min_copynum = 3

# Default parameters for pipe snps
# (Initial filtering for snps, done during postprocessing)
pipe_snps_min_samples=4
pipe_snps_min_nonzero_median_coverage=5
pipe_snps_lower_depth_factor=0.3
pipe_snps_upper_depth_factor=3

parse_snps_min_freq = 0.05

# Comment this out
from parse_HMP_data import *
# and uncomment this
#from parse_simulated_data import *
# for isolate data