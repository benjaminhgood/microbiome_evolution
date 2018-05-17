import numpy
import sys
import config
import gzip

default_core_gene_filename = (config.data_directory+"core_genes.txt.gz")
default_stringent_core_gene_filename = (config.data_directory+"core_genes_stringent.txt.gz")


def parse_core_genes(desired_species_name="", core_gene_filename=default_core_gene_filename):
    
    core_genes = set()
    core_gene_file = gzip.GzipFile(core_gene_filename,"r")
    for line in core_gene_file:
        
        items = line.split(":")
        if len(items)<2:
            continue
            
        species_name = items[0].strip()
        gene_names = [subitem.strip() for subitem in items[1].split(",")]
        
        if (species_name==desired_species_name) or (desired_species_name==""):
            core_genes.update(gene_names)
            
    core_gene_file.close()    
    return core_genes

def parse_core_gene_map(core_gene_filename=default_core_gene_filename):
    
    core_genes = {}
    
    core_gene_file = gzip.GzipFile(core_gene_filename,"r")
    for line in core_gene_file:
        
        items = line.split(":")
        if len(items)<2:
            continue
            
        species_name = items[0].strip()
        gene_names = [subitem.strip() for subitem in items[1].split(",")]
        core_genes[species_name] = gene_names
    
    core_gene_file.close()    
    return core_genes

# Actually calculate the core genes
if __name__=='__main__':
    
    import parse_midas_data
    
    pangenome_species = parse_midas_data.parse_good_species_list()
 
    cmin = config.core_genome_min_copynum
    cmax = config.core_genome_max_copynum  
    min_good_fraction = config.core_genome_min_prevalence
    min_coverage = 5 # (for assessing core genome, use a lower coverage value)
    
    output_filename = default_core_gene_filename
    output_file = gzip.GzipFile(output_filename,"w")

    stringent_output_filename = default_stringent_core_gene_filename
    stringent_output_file = gzip.GzipFile(stringent_output_filename,"w")


    for species_name in pangenome_species:
      
        # Load gene coverage information for species_name
        sys.stderr.write("Loading pangenome data for %s...\n" % species_name)
        gene_samples, gene_names, gene_presence_matrix, gene_depth_matrix, marker_coverages, gene_reads_matrix = parse_midas_data.parse_pangenome_data(species_name)
        sys.stderr.write("Done!\n")  
        
        if len(marker_coverages)==0:
            continue
            
        high_coverage_idxs = (marker_coverages>=min_coverage)

        if high_coverage_idxs.sum() < 0.5:
            continue

        gene_names = numpy.array(gene_names)
        gene_samples = gene_samples[high_coverage_idxs]
        marker_coverages = marker_coverages[high_coverage_idxs]
        gene_depth_matrix = gene_depth_matrix[:,high_coverage_idxs] 
        gene_copynum_matrix = gene_depth_matrix*1.0/(marker_coverages+(marker_coverages==0))

        # Load reference genes
        sys.stderr.write("Loading genes on reference genome..\n")
        reference_genes = parse_midas_data.load_reference_genes(species_name)
        sys.stderr.write("Done!\n")

        reference_gene_idxs = numpy.array([gene_name in reference_genes for gene_name in gene_names])
        
        # calculating good genes
        good_idxs = (((gene_copynum_matrix>=cmin)*(gene_copynum_matrix<=cmax)).sum(axis=1)*1.0/len(marker_coverages) >= min_good_fraction)
        good_gene_names = gene_names[good_idxs*reference_gene_idxs]
        sys.stderr.write("%d core genes out of %d\n" % (len(good_gene_names), len(gene_names)))
        output_file.write("%s: %s\n" % (species_name, ", ".join([gene_name for gene_name in good_gene_names])))
        
        # calculating good genes w/ stringent definition (100%)
        bad_idxs = (gene_copynum_matrix<config.gainloss_max_absent_copynum).sum(axis=1) > 0.5
        good_idxs = numpy.logical_not(bad_idxs)
        good_gene_names = gene_names[good_idxs*reference_gene_idxs]
        sys.stderr.write("%d stringent core genes out of %d\n" % (len(good_gene_names), len(gene_names)))
        stringent_output_file.write("%s: %s\n" % (species_name, ", ".join([gene_name for gene_name in good_gene_names])))
    
    output_file.close()
    stringent_output_file.close()
    
    
    