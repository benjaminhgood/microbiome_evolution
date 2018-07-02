import numpy
import sys
import config
import gzip
import os.path
import midas_db_utils

default_external_shared_gene_filename = (config.data_directory+"external_shared_genes.txt.gz")
default_external_core_gene_filename = (config.data_directory+"external_core_genes.txt.gz")
default_external_stringent_core_gene_filename = (config.data_directory+"core_genes_stringent.txt.gz")

default_shared_gene_filename = (config.data_directory+"shared_genes.txt.gz")
default_core_gene_filename = (config.data_directory+"core_genes.txt.gz")
default_stringent_core_gene_filename = (config.data_directory+"core_genes_stringent.txt.gz")

def parse_core_genes(desired_species_name="", core_gene_filename=default_core_gene_filename, external_core_gene_filename=default_external_core_gene_filename, external_filtering=True):
    
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
    
    external_core_genes = set()
    if os.path.isfile(external_core_gene_filename):
        external_core_gene_file = gzip.GzipFile(external_core_gene_filename,"r")
        
        for line in external_core_gene_file:
        
            items = line.split(":")
            if len(items)<2:
                continue
            
            species_name = items[0].strip()
            gene_names = [subitem.strip() for subitem in items[1].split(",")]
        
            if (species_name==desired_species_name) or (desired_species_name==""):
                external_core_genes.update(gene_names)
            
        external_core_gene_file.close() 
    
    if external_filtering and len(external_core_genes)>0:
        # some externally provided core genes
        core_genes = (core_genes & external_core_genes)
       
    return core_genes
    
    
def parse_shared_genes(desired_species_name="", shared_gene_filename=default_shared_gene_filename, external_shared_gene_filename=default_external_shared_gene_filename, external_filtering=True):
    
    shared_genes = set()
    shared_gene_file = gzip.GzipFile(shared_gene_filename,"r")
    for line in shared_gene_file:
        
        items = line.split(":")
        if len(items)<2:
            continue
            
        species_name = items[0].strip()
        gene_names_str = items[1].strip()
        if gene_names_str.startswith('N/A'): # Wasn't enough pangenome data to detect shared genes
            gene_names = []
        else:
            gene_names = [subitem.strip() for subitem in gene_names_str.split(",")]
        
        if (species_name==desired_species_name) or (desired_species_name==""):
            shared_genes.update(gene_names)
            
    shared_gene_file.close() 
    
    external_shared_genes = set()
    if os.path.isfile(external_shared_gene_filename):
        external_shared_gene_file = gzip.GzipFile(external_shared_gene_filename,"r")
        
        for line in external_shared_gene_file:
        
            items = line.split(":")
            if len(items)<2:
                continue
            
            species_name = items[0].strip()
            gene_names_str = items[1].strip()
            if gene_names_str.startswith('N/A'): # Wasn't enough pangenome data to detect shared genes
                gene_names = []
            else:
                gene_names = [subitem.strip() for subitem in gene_names_str.split(",")]
        
            if (species_name==desired_species_name) or (desired_species_name==""):
                external_shared_genes.update(gene_names)
            
        external_shared_gene_file.close() 
    
    if external_filtering and len(external_shared_genes)>0:
        # some externally provided core genes
        shared_genes = (shared_genes | external_shared_genes)
       
    return shared_genes

def parse_non_shared_reference_genes(desired_species_name="", shared_gene_filename=default_shared_gene_filename, external_shared_gene_filename=default_external_shared_gene_filename, external_filtering=True):
    import parse_midas_data
    shared_genes = parse_shared_genes(desired_species_name, shared_gene_filename, external_shared_gene_filename, external_filtering)
    reference_genes = parse_midas_data.load_reference_genes(desired_species_name)
    return set(reference_genes)-shared_genes
    
def get_good_pangenome_samples(marker_coverages, gene_copynum_matrix):

    cmin = config.core_genome_min_copynum
    cmax = config.core_genome_max_copynum  

    # Load reference genes
    num_reference_genes = len(parse_midas_data.load_reference_genes(species_name))
    
    num_present_genes = (gene_copynum_matrix>cmin).sum(axis=0)
    num_high_genes = (gene_copynum_matrix>cmax).sum(axis=0)
    
    good_sample_idxs = (num_present_genes>0.3*num_reference_genes)*(num_high_genes<0.3*num_present_genes)
    return good_sample_idxs
    
# Actually calculate the core genes
if __name__=='__main__':
    
    import parse_midas_data
    
    pangenome_species = parse_midas_data.parse_good_species_list()
 
    cmin = config.core_genome_min_copynum
    cmax = config.core_genome_max_copynum  
    shared_cmin = config.shared_genome_min_copynum
    
    min_good_fraction = config.core_genome_min_prevalence
    min_coverage = 5 # (for assessing core genome, we'll use a lower coverage value than when we look at real changes)
    
    output_filename = default_core_gene_filename
    output_file = gzip.GzipFile(output_filename,"w")

    stringent_output_filename = default_stringent_core_gene_filename
    stringent_output_file = gzip.GzipFile(stringent_output_filename,"w")
    
    shared_output_file = gzip.GzipFile(default_shared_gene_filename,"w")
    
    for species_name in pangenome_species:
      
        # Load reference genes
        sys.stderr.write("Loading genes on reference genome..\n")
        reference_genes = midas_db_utils.load_reference_genes(species_name)
        sys.stderr.write("Done!\n")

        # Load reference genes
        sys.stderr.write("Loading shared genes from midas db..\n")
        midas_shared_genes = midas_db_utils.parse_midas_shared_genes(species_name)             
        sys.stderr.write("Done!\n")


        bad_pangenome_data = False
                        
        # Load gene coverage information for species_name
        sys.stderr.write("Loading pangenome data for %s...\n" % species_name)
        gene_samples, gene_names, gene_presence_matrix, gene_depth_matrix, marker_coverages, gene_reads_matrix = parse_midas_data.parse_pangenome_data(species_name)
        sys.stderr.write("Done!\n")  
        
        if len(marker_coverages)==0:
            bad_pangenome_data = True
        else:        
            
            high_coverage_idxs = (marker_coverages>=min_coverage)

            if high_coverage_idxs.sum() < 0.5:
                bad_pangenome_data = True

        if bad_pangenome_data:
            # Just use reference genes
            shared_gene_names = sorted(midas_shared_genes)
            core_gene_names = sorted(reference_genes - midas_shared_genes)
            stringent_gene_names = sorted(reference_genes - midas_shared_genes)
      
        else:    
        
            gene_names = numpy.array(gene_names)
            gene_samples = gene_samples[high_coverage_idxs]
            marker_coverages = marker_coverages[high_coverage_idxs]
            gene_depth_matrix = gene_depth_matrix[:,high_coverage_idxs] 
            gene_copynum_matrix = gene_depth_matrix*1.0/(marker_coverages+(marker_coverages==0))

            good_sample_idxs = get_good_pangenome_samples(marker_coverages, gene_copynum_matrix)
            bad_sample_idxs = numpy.logical_not(good_sample_idxs)
                    
            sys.stderr.write("%d bad samples!\n" % bad_sample_idxs.sum())
            
            gene_samples = gene_samples[good_sample_idxs]
            marker_coverages = marker_coverages[good_sample_idxs]
            gene_copynum_matrix = gene_copynum_matrix[:,good_sample_idxs]
            
            reference_gene_idxs = numpy.array([gene_name in reference_genes for gene_name in gene_names])
            
            midas_shared_idxs = numpy.array([gene_name in midas_shared_genes for gene_name in gene_names])
        
            # These are genes that have coverage >=3x normal in some sample. This are candidates for being linked to another species.
            # (they could also be multi-copy genes, but we can't look at much on these genes anyway, so might as well toss them out)
            shared_idxs = ((gene_copynum_matrix>shared_cmin).sum(axis=1)>0.5)
            
            sys.stderr.write("%d putative shared genes\n" % (shared_idxs.sum()))
            sys.stderr.write("%d shared genes in db\n" % (midas_shared_idxs.sum()))
            
            # Now union with those we identified from midas db
            shared_idxs = numpy.logical_or(shared_idxs, midas_shared_idxs)
            non_shared_idxs = numpy.logical_not(shared_idxs)
            shared_gene_names = gene_names[shared_idxs]
            sys.stderr.write("%d shared genes out of %d\n" % (len(shared_gene_names), len(gene_names)))
            sys.stderr.write("(%d out of %d on reference genome)\n" % ((shared_idxs*reference_gene_idxs).sum(), reference_gene_idxs.sum())) 
            # calculating good genes
            good_idxs = (((gene_copynum_matrix>=cmin)*(gene_copynum_matrix<=cmax)).sum(axis=1)*1.0/len(marker_coverages) >= min_good_fraction)
            core_gene_names = gene_names[good_idxs*reference_gene_idxs*non_shared_idxs]
            sys.stderr.write("%d core genes out of %d\n" % (len(core_gene_names), len(gene_names)))
            sys.stderr.write("%d w/ relaxed criteria\n" % (good_idxs*reference_gene_idxs).sum())
            
            # calculating good genes w/ stringent definition (100%)
            bad_idxs = (gene_copynum_matrix<config.gainloss_max_absent_copynum).sum(axis=1) > 0.5
            good_idxs = numpy.logical_not(bad_idxs)
            stringent_gene_names = gene_names[good_idxs*reference_gene_idxs*non_shared_idxs]
            #sys.stderr.write("%d stringent core genes out of %d\n" % (len(stringent_gene_names), len(gene_names)))
        
        # Write output to file!
        shared_output_file.write("%s: %s\n" % (species_name, ", ".join([gene_name for gene_name in shared_gene_names])))
        output_file.write("%s: %s\n" % (species_name, ", ".join([gene_name for gene_name in core_gene_names])))
        stringent_output_file.write("%s: %s\n" % (species_name, ", ".join([gene_name for gene_name in stringent_gene_names])))
    
    shared_output_file.close()
    output_file.close()
    stringent_output_file.close()
    
    
    