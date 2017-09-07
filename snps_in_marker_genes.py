import parse_midas_data
import diversity_utils

# read in the centroid fasta file:
species_name='Bacteroides_vulgatus_57955' 
centroid_fastas = parse_midas_data.load_centroid_fasta(species_name)

# marker gene of interest: B000096

B_vul_gene='435590.9.peg.1355' # this is the reference genome gene and is B. vulgatus
#centroid_genes=['483217.6.peg.1552','997875.3.peg.2812'] # both are B. dorei
B_dorei_gene='483217.6.peg.1552' # compare iwith one B. dorei genome because the two dorei genomes are identical at positions that differ from B. vul.

B_vul_fasta=centroid_fastas[B_vul_gene]
B_dorei_fasta=centroid_fastas[B_dorei_gene]


alignment=diversity_utils.find_snps_in_gene_pair(B_vul_fasta,B_dorei_fasta)
print alignment
