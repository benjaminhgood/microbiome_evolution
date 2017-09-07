import parse_midas_data


# read in the centroid fasta file:
species_name='Bacteroides_vulgatus_57955' 
centroid_fastas = parse_midas_data.load_centroid_fasta(species_name)

# marker gene of interest: B000096

ref_gene='435590.9.peg.1355' # this is the reference genome gene and is B. vulgatus
#centroid_genes=['483217.6.peg.1552','997875.3.peg.2812'] # both are B. dorei
centroid_genes=['483217.6.peg.1552']# iterating through one B. dorei genome because the two dorei genomes are identical at positions that differ from B. vul.

# fasta sequence of B. vulgatus gene:
ref_fasta=centroid_fastas[ref_gene]

alignment={}
# key=bp 
# value=[B. vul, B. dorei]

for bp in range(0, len(ref_fasta)):
    ref_bp=ref_fasta[bp]
    for gene in centroid_genes:
        fasta=centroid_fastas[gene]
        if fasta[bp] != ref_bp:
            alignment[bp]=[ref_bp, fasta[bp]]

print alignment
