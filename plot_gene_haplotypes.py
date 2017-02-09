import matplotlib  
matplotlib.use('Agg') 
import parse_midas_data
import pylab
import sys
import numpy
import diversity_utils
from numpy.random import choice
import os
species_name=sys.argv[1]



# Load subject and sample metadata
sys.stderr.write("Loading HMP metadata...\n")
subject_sample_map = parse_midas_data.parse_subject_sample_map()
sys.stderr.write("Done!\n")
    
# Load SNP information for species_name
sys.stderr.write("Loading %s...\n" % species_name)
samples, allele_counts_map, passed_sites_map = parse_midas_data.parse_snps(species_name, combination_type="sample", debug=False)
sys.stderr.write("Done!\n")
    

# identify which samples have low piS -- use these to construct haplotypes
pi_matrix_syn, avg_pi_matrix_syn = diversity_utils.calculate_pi_matrix(allele_counts_map, passed_sites_map, variant_type='4D')
low_diversity_samples = (numpy.diag(avg_pi_matrix_syn)<1e-03)# the diagonal of the above has the within-patient pi value. 
unique_samples = parse_midas_data.calculate_unique_samples(subject_sample_map, samples)
desired_samples = unique_samples*low_diversity_samples

for gene_name in allele_counts_map.keys():
    locations_4D = numpy.array([location for chromosome, location in allele_counts_map[gene_name]['4D']['locations']])*1.0
    locations_1D = numpy.array([location for chromosome, location in allele_counts_map[gene_name]['1D']['locations']])*1.0

    #check if there are at least 100 bps in this gene
    gene_length=len(locations_4D) + len(locations_1D)
    if gene_length >=100:

        # create a dictionary keeping track of whether a location is 1D or 4D (key=location, value=[index,1D/4D]).
        location_dictionary={}
        for loc in range(0, len(locations_4D)):
            location_dictionary[locations_4D[loc]]=[loc,'4D']
        for loc in range(0, len(locations_1D)):
            location_dictionary[locations_1D[loc]]=[loc,'1D']

        allele_counts_4D = allele_counts_map[gene_name]['4D']['alleles']
        allele_counts_1D = allele_counts_map[gene_name]['1D']['alleles']

        allele_counts_4D = allele_counts_4D[:,desired_samples,:]
        allele_counts_1D = allele_counts_1D[:,desired_samples,:]

        # Create two files for a gene haplotype:
        # (1) A file with the haplotypes themselves
        # (1) annotation file with numbers 0,1,2,3,4 to indicate:
    
        #0: p=0, site is unmutated
        #1: p=1, site is a fixed syn difference (4D)
        #2: p=1, site is a fixed nonsyn difference (1D)
        #3: p<1, site is a polymorphic syn within patient (4D)
        #4: p>1, site is a polymorphic nonsyn within patient (1D)

        # generate_haplotype produces the following files:
        # (1) tmp_consensus.txt
        # (2) tmp_anno.txt
        diversity_utils.generate_haplotype(allele_counts_4D, allele_counts_1D, location_dictionary)

        # Cluster the haplotypes by identity
        num_samples=sum(desired_samples)
        os.system('python ~/projectBenNandita/H12_H2H1_MIDAS.py tmp_consensus.txt ' + str(num_samples) + ' -o tmp_cluster.txt -g ' + gene_name)  

        # Plot the haplotypes with R
        os.system('mkdir -p ~/ben_nandita_hmp_analysis/hap_plots/' + species_name)
        os.system('Rscript ~/projectBenNandita/visualizeHaplotypesMicrobiome3.R tmp_cluster.txt '+ str(num_samples) + ' ' +  str(gene_length) + ' tmp_consensus.txt tmp_anno.txt ~/ben_nandita_hmp_analysis/hap_plots/' + species_name +'/hap_plot_' + gene_name +'.pdf ' + gene_name)
