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
#
#########################################################################################

def load_kegg_annotations(gene_names):
    
    # dictionary to store the kegg ids (gene_id -> [[kegg_id, description]])
    kegg_ids={}
    

    genomes_visited=[]
    for i in range(0, len(gene_names)):
        genome_id='.'.join(gene_names[i].split('.')[0:2])
        if genome_id not in genomes_visited:
            genomes_visited.append(genome_id)
            file= bz2.BZ2File("%skegg/%s.kegg.txt.bz2" % (parse_midas_data.patric_directory, genome_id),"r")
            file.readline() #header  
            file.readline() #blank line
            for line in file:
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

#######################    

if __name__=='__main__':

    pass