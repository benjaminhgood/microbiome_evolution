import numpy
import sys
import bz2
import gzip
import os.path 
import stats_utils
from math import floor, ceil
import gene_diversity_utils
import parse_midas_data
import config

#########################################################################################
#
# Read in the Kegg info for a given speceis
# What is returned are all pathways for a given species.
#
#########################################################################################

def load_kegg_annotations(genome_ids):
    
    # dictionary to store the kegg ids (gene_id -> [[kegg_id, description]])
    kegg_ids={}
    

    genomes_visited=[] #check if I have already loaded the genome for this gene
    for genome_id in genome_ids:
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
# This returns a dictionary with all gene names for the genomes included in the genome_ids object.
# In the main code, I will pull out the actual gene names. 
#
#################################################################

def load_patric_gene_descriptions(genome_ids,  allowed_genes=[]):
    
    allowed_genes = set(allowed_genes)
    # dictionary to store all gene names (gene_id -> )
    gene_descriptions={}

    for genome_id in genome_ids:
        file=gzip.open('%s/features/%s.PATRIC.features.tab.gz' % (config.patric_directory, genome_id), 'r') 
        file.readline() #header  
        for line in file:
            items = line.strip().split("\t")
            if items[0] !='':
                if items[5] !='' and len(items)>14: # sometimes entries are blank
                    gene_id =  items[5].split('|')[1] # id of gene
                    if gene_id in allowed_genes:
                        gene_description = items[14] # what the gene does
                        gene_descriptions[gene_id] = gene_description # load into the dictionary

    return gene_descriptions
    
###########################################################
# 
# Categorize PATRIC gene descriptions by  regular expression
#
###########################################################
import operator
def cluster_patric_gene_descriptions(gene_descriptions):
    
    
    #iterate through and alphabetically categorize genes based on their string identity. If there are at most 2 string mismatches with the previous string, then clump it with that string.

    gene_categories={}  # key=gene_name, value=number of genes in this category
    gene_category_map={} # key=gene_id, value=category 
    # I'm making this map so that we can later look up which category a patric id belongs in. 

    prev_gene='' # I will do a regexp with this as I iterate and keep updating
    gene_categories[prev_gene]=0 #initialize
    #iterate through alphabetically (faster)
    for item in sorted(gene_descriptions.items(), key=operator.itemgetter(1)):
        gene_id=item[0]
        gene=item[1]
        hamming_distance= hamming(gene, prev_gene) 
        if hamming_distance<=2:
            gene_categories[prev_gene]+=1
            gene_category_map[gene_id]=prev_gene
        else:
            # sometimes the alpha sort doesn't take care of corner cases. Iterate through the whole list again to find an existing match if possible. 
            found_category=False
            for existing_gene in gene_categories.keys():
                hamming_distance= hamming(gene, existing_gene)
                if hamming_distance <=2 and found_category==False:
                    gene_categories[existing_gene] +=1
                    gene_category_map[gene_id]=existing_gene
                    found_category=True
            # if no match is found, create a new category
            if found_category==False:
                gene_categories[gene]=1
                gene_category_map[gene_id]=gene
                prev_gene=gene

    return gene_categories, gene_category_map

#########
# hamming distance between two strings:

import itertools

def hamming(str1, str2):
  diff=sum(itertools.imap(str.__ne__, str1, str2))
  diff += abs(len(str1) - len(str2)) # above doesn't take into account difference in string length
  return diff
  
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
            # NRG: added 'accn|' to match the headers in the fasta file (09/06/17)
            outFile.write(gene_id +'\t' +'accn|'+scaffold_id +'\t'+start +'\t' + end +'\t' +strand +'\t' +gene_type +'\t' +functions +'\n')

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
            # This annotation is one of a few! Missed a few genomes (NRG, 09/06/17)
            if 'Reference genome for the Human Microbiome Project' in annotation:
                HMP_genomes[genome_id] = [int(contigs), int(genome_length), body_part, host]
        
    return HMP_genomes

#######################    

if __name__=='__main__':

    pass
