import sys
import numpy
import bz2
import gzip
import config
import os.path


intermediate_filename_template = config.data_directory+"snp_prevalences/%s.txt.gz"
    
# Loading file
def parse_snp_prevalences(desired_species_name):
    
    intermediate_filename = intermediate_filename_template % desired_species_name
        
    snp_prevalences = {}
    
    if not os.path.isfile(intermediate_filename):
        return snp_prevalences
    
    file = gzip.GzipFile(intermediate_filename,"r")
    file.readline()
    for line in file:
        items = line.split(",")
        contig = items[0]
        location = long(items[1])
        population_freq = float(items[2])
        snp_freq = float(items[3])
        
        snp_prevalences[(contig,location)] = snp_freq
                            
    file.close()
    
    return snp_prevalences
    
    
# Loading file
def parse_population_freqs(desired_species_name, polarize_by_consensus=False):
    
    intermediate_filename = intermediate_filename_template % desired_species_name
      
    population_freqs = {}
    
    if not os.path.isfile(intermediate_filename):
        return population_freqs

    file = gzip.GzipFile(intermediate_filename,"r")
    file.readline()
    for line in file:
        items = line.split(",")
        contig = items[0]
        location = long(items[1])
        population_freq = float(items[2])
        snp_freq = float(items[3])
        
        if polarize_by_consensus:
            if population_freq > 0.5:
                population_freq = 1-population_freq
        
        if population_freq==0:
            pass
        else:
            population_freqs[(contig,location)] = population_freq
                            
    file.close()
    
    return population_freqs

if __name__=='__main__': 
    
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

    # Holds panel wide prevalence for each species
    os.system('mkdir -p %ssnp_prevalences' % config.data_directory)

    # Open post-processed MIDAS output
    snp_file =  bz2.BZ2File("%ssnps/%s/annotated_snps.txt.bz2" % (config.data_directory, species_name),"r")
    
    line = snp_file.readline() # header
    items = line.split()[1:]
    samples = numpy.array([item.strip() for item in items])
    
    record_strs = ["Chromosome, Location, AltFreq, SNPFreq"]    
    
    sys.stderr.write("Calculating SNP prevalences...\n")        
    num_sites_processed = 0
    for line in snp_file:
        
        num_sites_processed+=1
        #print num_sites_processed
        if num_sites_processed%50000==0:
            sys.stderr.write("%dk sites processed...\n" % (num_sites_processed/1000))   
            if debug:
                break
        
        items = line.split()
        # Load information about site
        info_items = items[0].split("|")
        chromosome = info_items[0]
        location = long(info_items[1])
        gene_name = info_items[2]
        
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
    
        # First calculate fraction of cohort where alternate (non-reference) allele is the major allele
        population_prevalence = ((alts>=refs)*(depths>0)).sum()
        population_freq = population_prevalence*1.0/(depths>0).sum()
    
        if population_freq>0.5:
            # alternate allele is in the majority
            # re-polarize for now
            alts,refs = refs,alts
            
        # Next calculate fraction of cohort where population minor allele is present at >=10% within-host frequency 
        alt_threshold = numpy.ceil(depths*0.1)+0.5 #at least one read above 10%.
        
        snp_prevalence = ((alts>=alt_threshold)*(depths>0)).sum()
        snp_freq = snp_prevalence*1.0/(depths>0).sum()
        
        if (population_prevalence==0) and (snp_prevalence==0):
            continue
        
        # otherwise record data    
        record_str = "%s, %d, %g, %g" % (chromosome, location, population_freq, snp_freq) 
        
        record_strs.append(record_str)
        
    
    snp_file.close()
    sys.stderr.write("Done!\n")

    # Write to disk!
    intermediate_filename = intermediate_filename_template % species_name
    output_file = gzip.GzipFile(intermediate_filename,"w")
    output_file.write("\n".join(record_strs))
    output_file.close()



