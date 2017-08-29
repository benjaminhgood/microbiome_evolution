import numpy
import sys
import bz2
import gzip
import os.path 
import stats_utils
from math import floor, ceil
import gene_diversity_utils
import parse_midas_data

#########################################################################################
#
# Read in the Kegg info for a given speceis
# What is returned are all pathways for a given species rather than just the pathways for the genes in gene_names.
#
#########################################################################################

def load_kegg_annotations(gene_names):
    
    # dictionary to store the kegg ids (gene_id -> [[kegg_id, description]])
    kegg_ids={}
    

    genomes_visited=[] #check if I have already loaded the genome for this gene
    for i in range(0, len(gene_names)):
        genome_id='.'.join(gene_names[i].split('.')[0:2])
        if genome_id not in genomes_visited:
            genomes_visited.append(genome_id)
            file= bz2.BZ2File("%skegg/%s.kegg.txt.bz2" % (parse_midas_data.patric_directory, genome_id),"r")
            file.readline() #header  
            for line in file:
                if line.split('\t')[0]!='':
                    if line.strip() != "":
                        items = line.split("\t")
                        gene_name=items[0].strip().split('|')[1]
                        kegg_ids[gene_name]=[]
                        kegg_pathway_tmp=items[1].strip().split(';')
                        if len(kegg_pathway_tmp)>0 and kegg_pathway_tmp[0] !='':
                            for i in range(0, len(kegg_pathway_tmp)):
                                kegg_ids[gene_name].append(kegg_pathway_tmp[i].split('|'))
                        elif kegg_pathway_tmp[0] =='':
                            kegg_ids[gene_name].append(['',''])
    return kegg_ids

########################################
def load_spgenes_annotations(gene_names):
    
    # dictionary to store the special gene ids (gene_id -> [property,product])
    spgenes_ids={}
    genomes_visited=[]
    for gene_name in gene_names: 
        genome_id='.'.join(gene_name.split('.')[0:2])
        if genome_id not in genomes_visited:
            genomes_visited.append(genome_id)
            file= gzip.open("%spatric_spgene/%s.PATRIC.spgene.tab.gz" % (parse_midas_data.patric_directory, genome_id),"r")
            file.readline() #header  
            for line in file:
                if line.strip() != "":
                    items = line.split("\t")
                    gene_name=items[2].strip().split('|')[1]
                    product=items[6]
                    property=items[7]
                    spgenes_ids[gene_name]=[[property,product]]
    
    return spgenes_ids

def load_antibiotic_resistance_genes(species_name):
    
    # get pangenome genome for species_name
    
    pangenome_genes = parse_midas_data.load_pangenome_genes(species_name)
    spgenes_ids = load_spgenes_annotations(pangenome_genes)
    
    antibiotic_resistance_genes = set([])
    
    for gene_name in spgenes_ids.keys():
    
        if spgenes_ids[gene_name][0][0] == 'Antibiotic Resistance':
            antibiotic_resistance_genes.add(gene_name)
            
    return antibiotic_resistance_genes

def load_virulence_factors(species_name):
    
    # get pangenome genome for species_name
    
    pangenome_genes = parse_midas_data.load_pangenome_genes(species_name)
    spgenes_ids = load_spgenes_annotations(pangenome_genes)
    
    virulence_genes = set([])
    
    for gene_name in spgenes_ids.keys():
    
        if spgenes_ids[gene_name][0][0] == 'Virulence Factor':
            virulence_genes.add(gene_name)
            
    return virulence_genes



#################################################################
#
# Load individual gene names from patric
# This returns a dictionary with all gene names for the genomes included in the gene_names object.
# In the main code, I will pull out the actual gene names. 
#
#################################################################

def load_patric_gene_descriptions(gene_ids):
    
    # dictionary to store all gene names (gene_id -> )
    gene_descriptions={}

    genomes_visited=[] #check if I have already loaded the genome for this gene
    for i in range(0, len(gene_ids)):
        genome_id='.'.join(gene_ids[i].split('.')[0:2])
        if genome_id not in genomes_visited:
            genomes_visited.append(genome_id)
            file=gzip.open('/pollard/shattuck0/snayfach/databases/PATRIC/genomes/%s/%s.PATRIC.features.tab.gz' % (genome_id, genome_id), 'r') # update once I pull out gene names for Ben?
            file.readline() #header  
            for line in file:
                items = line.strip().split("\t")
                if items[0] !='':
                    if items[5] !='' and len(items)>14: # sometimes entries are blank
                        gene_id =  items[5].split('|')[1] # id of gene
                        gene_description = items[14] # what the gene does
                        gene_descriptions[gene_id] = gene_description # load into the dictionary

    return gene_descriptions
    


######################################################################################
# 
# Create a new genome.features.gz file for running MIDAS on a different representative genome for SNP calling. 
#
####################################################################################
def new_genome_features_file(genome_id, outFN):
    
    pollard_patric_dir='/pollard/shattuck0/snayfach/databases/PATRIC/genomes'

    #outFile=gzip.open('/pollard/home/ngarud/BenNanditaProject/MIDAS_ref_genome_test/genome_features_files/%s_features.gz' % genome_id,"w")

    outFile=gzip.open(outFN, "w")    
    outFile.write("gene_id\tscaffold_id\tstart\tend\tstrand\tgene_type\tfunctions\n")
    
    for genome_part in ['cds','rna']:
        file= gzip.open("%s/%s/%s.PATRIC.%s.tab.gz" % (pollard_patric_dir,genome_id, genome_id, genome_part),"r")
        file.readline() #header

        for line in file:
            items = line.split("\t")
            gene_id=items[5].strip().split('|')[1]
            scaffold_id=items[2].strip()
            start=items[9].strip()
            end=items[10].strip()
            strand=items[11].strip()
            gene_type=items[4].strip()
            if len(items)>15:
                functions=items[15].strip()
            else:
                functions=''
            outFile.write(gene_id +'\t' +scaffold_id +'\t'+start +'\t' + end +'\t' +strand +'\t' +gene_type +'\t' +functions +'\n')

######################################################################################
# 
# Read in the genome_metadata file. This tells you where genomes in PATRIC came from 
# for example, if we want genomes that are part of the HMP reference panel, this is the file to look into
#
####################################################################################
def get_HMP_reference_genomes():
    
    genome_metadata = open("/pollard/shattuck0/snayfach/databases/PATRIC/metadata/genome_metadata")
    HMP_genomes={}
    for line in genome_metadata:
        items=line.strip().split('\t')
        if len(items) == 65: # sometimes a line is shorter, this creates problems
            genome_id=items[0]
            species=items[1]
            annotation=items[64]
            contigs=items[27]
            genome_length=items[29]
            body_part = items[35]
            host = items[45]
            if 'Reference genome for the Human Microbiome Project' in annotation:
                HMP_genomes[genome_id] = [int(contigs), int(genome_length), body_part, host]
        
    return HMP_genomes

#######################    

if __name__=='__main__':

    pass
